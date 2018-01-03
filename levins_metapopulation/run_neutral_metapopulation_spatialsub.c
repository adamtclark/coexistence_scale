#include <R.h>
#include <Rmath.h>
#include <stdio.h>
#include "error.c"
double make_expnential_event(double lambda);
void minpos(double eventtimes_c[], double eventtimes_m [], int gridsize, double event_info[]);

//scheduling functions
void EventInit();
void EventStartTime(double t0);
void EventSchedule(int n, double te);
void EventCancel(int n);
int cancel1(int n, int i);
int EventNext();
void EventRenumber(int n, int m);
int sort(int list[], int p, int n, int (*c)(int,int));
int isort(int n);
int imerge(int p, int q);

//Initialize scheduling functions
#define PINIT  if(run1==0) EventInit();
#define indiv 2000000   //Maximum population size.
#define INDIV indiv       //Maximum population size.

#define PEMPTY -1      //Marker for bins containing no linkages.
#define TN (INDIV+0)   //Maximum number of time bins.
#define PN (INDIV+3)   //Maximum number of time bin forward indexes.
#define TW    20       //Time width of all bins combined (for optimization).
#define TINF 1E10      //Close enough to infinity for the scheduler.
double t;              //Current time, last dispatched event.

static int run1 = 0;   //Flag to detect if the routine is being reused.

static double T[PN];      //Time for each scheduled event.
static int P[PN];      //Forward indexes within bins, ending with zero.
static int Q[TN];      //First index for the bin, with zero for empty bins.

static int Qn  = TN;   //Number of elements in 'Q'.
static double Qw  = TW;   //Interval of time represented for each cycle in 'Q'.
static int Qi  = 0;    //Index of the immediate time bin.
static int Qo  = 1;    //Flag set if the immediate bin is in order.
static int Qe  = 0;    //Number of events in all bins.

static double Qt0 = 0;    //Earliest time representable this cycle in 'Q'.
static double Qt1 = TW;   //Earliest time beyond this cycle in 'Q'.

static int *Pl, pcurr, pprev, count;
static int (*order)(int,int);

void run_neutral_metapopulation_spatialsub(double *ptmax, int *pgridsize, int *pnsp, int xylim[], int destroyed[], //grid
		double c_sptraits[], double m_sptraits[], int abundances[], int colsites[], int *pncolsites, //traits
		double eventtimes_c[], double eventtimes_m[], //events
		int speciesid[], //species
		double output[], int *pnsteps, //outputs
		int *ptalktime, //tell us what time it is?
		int abundances_sub[], int sites_sub[], int *pnsites_sub, double output_sub[]) { //spatial subset

	int i=0;
	double counter = 0;
	int tmpsp=0;

	double tnow = 0; //Current time
	double tmax = *ptmax; //Maximum time for simulation
	int gridsize = *pgridsize; //Total size of trid (length x width)
	int nsp = *pnsp; //Number of species
	int nsteps = *pnsteps; //How many time steps to record?
	double dt = tmax/nsteps; //dt between recording time steps
	double trecord = dt; //Time at which to record system state
	int istep = 0; //How many times have we saved system state?
	int ncolsites = *pncolsites;
	int talktime = *ptalktime;
	int nsites_sub = *pnsites_sub;

	int cadj = 0; //plus what to go from R to C?

	int sppos = 0; //position in species list
	int eventtype = 0; //type of event that is happeneing

	GetRNGstate();  // Load seed for random number generation, R base

	//Output time and species abundances
	output_sub[0+cadj] = 0; //First time step
	for(i=0; i<nsp; i++) {
	 output_sub[(nsteps+1)*(i+1)+cadj] = abundances_sub[i+cadj]; //Each species		
	}

	output[0+cadj] = 0; //First time step
	for(i=0; i<nsp; i++) {
	 output[(nsteps+1)*(i+1)+cadj] = abundances[i+cadj]; //Each species		
	}

	istep=istep+1;

	//Loop until we reach tmax
	int totabund=0; //Stop if we run out of individuals

	for(i=0; i<nsp; i++) {
	 totabund=totabund+abundances[i+cadj];
	}

	//set up events
	EventStartTime(tnow);
	for(i=0; i<gridsize; i++) {
		if(eventtimes_c[i]!=0 || eventtimes_m[i]!=0) {
			EventSchedule(i+1, eventtimes_c[i]);
			EventSchedule(i+gridsize+1, eventtimes_m[i]);
		}
	}


	//start main loop
	while((tnow<tmax) && (totabund>0)) {

	  //find minimum next position
	  double event_info[3] = {0, 0, 0};

	  //Find next event
	  int nl = EventNext();

	  if(nl==0) {
	  	break; //out of events if all are done
	  }

	  nl=nl-1;

	  if(nl>=gridsize && nl <= (2*gridsize-1)) {
	  	event_info[0]=1; //death event
	  	event_info[2]=(nl-gridsize); //event position
	  	event_info[1]=t; //event time
	  } else if(nl<gridsize) {
	  	event_info[0]=2; //birth event
	  	event_info[2]=nl; //event position
	  	event_info[1]=t; //event time
	  } else {
	  	fprintf(stderr, "ERROR - event too large, = %d \n", nl);
	  	break; //error - nl > nind
	  }
	  
	 tnow = event_info[1+cadj]; //time now
	 sppos = (int) event_info[2+cadj]; //position of species with event
	 eventtype = (int) event_info[0+cadj]; //type of event

	 if(eventtype==1) { //death event
	   
	   //remove individual for abundane list
	   abundances[speciesid[sppos+cadj]+cadj]=abundances[speciesid[sppos+cadj]+cadj]-1; //remove one from species abundance
	   totabund=totabund-1;
	                                          
	   //remove events attached to species
	   EventCancel(sppos+1); //cancel birth

	   eventtimes_c[sppos+cadj]=0;
	   eventtimes_m[sppos+cadj]=0;
	   speciesid[sppos+cadj]=nsp;
	   
	 } else if(eventtype==2) { //birth event
	   //generate propagule
	   //Remove birth event that just happened
	   eventtimes_c[sppos+cadj]=0;
	   
	   //generate new birth event for parent
	   eventtimes_c[sppos+cadj]=tnow+make_expnential_event(c_sptraits[speciesid[sppos+cadj]+cadj]);

	   EventSchedule(sppos+1, eventtimes_c[sppos+cadj]);
	   
	   //make new birth location
	   int currentx = floor(sppos/xylim[0+cadj]); // current x location of parent
	   int currenty = sppos-xylim[0+cadj]*currentx; // current y location of parent
	   
	   int newlocation=floor(runif(0,1)*ncolsites);
	   int newx = colsites[newlocation+cadj]+currentx;
	   int newy = colsites[newlocation+ncolsites+cadj]+currenty;
	   
	   //wrap taurus
	   while(newx < 0) {
	     newx = xylim[0+cadj]-newx;
	   }
	   while(newx >= xylim[0+cadj]) {
	     newx = newx - xylim[0+cadj];
	   }
	   while(newy < 0) {
	     newy = xylim[1+cadj]-newy;
	   }
	   while(newy >= xylim[1+cadj]) {
	     newy = newy - xylim[1+cadj];
	   }
	   //translate xy back into list position
	   int newpos = xylim[0+cadj]*newx+newy;
	   //note x is number of columns, y is number of rows
	   
	   ////check for habitat destruction, etc.
	   int trigger=0; //will switch to 1 if location cannot be colonized
	   if(destroyed[newpos+cadj]==1) {
	     trigger=1;
	   } else if(speciesid[newpos+cadj]!=nsp) { ////check for superior competitor (or same species)
	     trigger=1;
	   }

	   ////if successful, generate birth and death events for offspring
	   if(trigger==0) {
	     //add species
	     eventtimes_c[newpos+cadj]=tnow+make_expnential_event(c_sptraits[speciesid[sppos+cadj]+cadj]);
	     eventtimes_m[newpos+cadj]=tnow+make_expnential_event(m_sptraits[speciesid[sppos+cadj]+cadj]);

	     EventSchedule(newpos+1, eventtimes_c[newpos+cadj]);
	     EventSchedule(newpos+gridsize+1, eventtimes_m[newpos+cadj]);

	     speciesid[newpos+cadj]=speciesid[sppos+cadj];
	     
	     abundances[speciesid[newpos+cadj]+cadj]=abundances[speciesid[newpos+cadj]+cadj]+1;
	     totabund=totabund+1;
	   }
	 }
	   
	 //save state
	 if((tnow > trecord)){
	   trecord=tnow+dt;
	   
	   output[istep+cadj] = tnow; //First time step
	   output_sub[istep+cadj] = tnow; //First time step
	  
	   for(i=0; i<nsp; i++) {
	     output[(nsteps+1)*(i+1)+istep+cadj] = abundances[i+cadj]; //Each species		
	   }

	   //record spatial subset
	   for(i=0; i<nsp; i++) {
	     abundances_sub[i+cadj] = 0;
	   }

	   for(i=0; i<nsites_sub; i++) {
	   	 tmpsp = speciesid[sites_sub[i]];
	   	 if(tmpsp<nsp) {
	   	 	abundances_sub[tmpsp] = abundances_sub[tmpsp]+1;
	   	 }
	   }

	   for(i=0; i<nsp; i++) {
	   	 output_sub[(nsteps+1)*(i+1)+istep+cadj] = abundances_sub[i+cadj]; //Each species		
	   }

	   istep=istep+1;

	   if(talktime==1) {
	   	fprintf(stderr, "time = %f \n", trecord);
	   }
	 }

	 counter=counter+1;
	 if(counter/100 == floor(counter/100)) {
	   //Allow for user interrupt
	     R_CheckUserInterrupt();
	 }
	}
}

double make_expnential_event(double lambda) {
  //Inverse exponential distribution for waiting time with rate lambda
  
  double x = runif(0,1);
  double tm = log(-x+1)/(-lambda);
  return(tm);
}

//Scheduling functions
void EventInit() {
  int i;

  for(i=0; i<PN; i++) P[i] = PEMPTY;
  if(run1==0) { run1 = 1; return; }

  for(i=0; i<PN; i++) T[i] = 0;
  for(i=0; i<TN; i++) Q[i] = 0;

  Qn  = TN; Qw  = TW; Qi  = 0;
  Qo  = 1;  Qe  = 0;
  Qt0 = 0;  Qt1 = TW;

  t = 0;
}

void EventStartTime(double t0) {
  PINIT;                                     //Initialize if necessary.

  if(Qe) Error(742.);                        //Make sure the bins are empty.

  Qt0 = t0-(Qw/Qn)/2;                        //Set the time boundaries, leaving
  Qt1 = Qt0+Qw;                              //room for rounding errors.

  t = t0;                                    //Set the global time.
}

void EventSchedule(int n, double te) {
  int i; double tr;

  PINIT;                                     //Initialize if necessary.

  if(n<1||n>=PN)   Error1(734.1, "n=",n);    //Check the index and make sure an
  if(P[n]!=PEMPTY) Error1(735.1, "n=",n);    //event is not already scheduled
  if(te>=TINF)     Error1(737.0, "n=",n);    //is not in the infinite future,

  if(te<t) Error2(738.1, "t=",t, ">",te);    //and is not in the past.

  T[n] = te;                                 //Record the time of the new event.

  tr = (te-Qt0)/Qw; tr -= (int)tr;           //Convert the time to a bin number
  i  = tr*Qn; if(i==Qi) Qo = 0;              //and mark for sorting if needed.

  P[n] = Q[i]; Q[i] = n;                     //Add the event to the list for
  Qe += 1;                                   //that bin and increment the number
} 

void EventCancel(int n) {
  int i;
  double tr;

  PINIT;                                     //Initialize if necessary.

  if(n<1||n>=PN)   Error1(734.2, "n=",n);    //Check the index and make sure an
  if(P[n]==PEMPTY) Error1(736.2, "n=",n);    //event is scheduled.

  tr = (T[n]-Qt0)/Qw; tr -= (int)tr;         //Convert the time to a bin number,
  i  = tr*Qn;                                //modulo the duration of the cycle.

  if(cancel1(n, i)) return;                  //Remove it from its normal bin.

  i = (i-1+Qn) % Qn;                         //If it is not there, check the
  if(cancel1(n, i)) return;                  //bin below (from rounding error).

  i = (i+2+Qn) % Qn;                         //Finally, check the bin above.
  if(cancel1(n, i)) return;

  Error1(818., "n=",n);                      //If the specified event was not in
}                                            //the list, signal an error.

int cancel1(int n, int i) {
  int j, jp;

  for(jp=0,j=Q[i]; j>0; jp=j,j=P[j])         //Scan the list of pending events
    if(j==n)                                 //in this bin and remove the
    { if(jp>0) P[jp] = P[j];                 //specified event. (The average
      else     Q[i]  = P[j];                 //number events in a non-empty bin
      P[j] = PEMPTY;                         //is 1.5)
      Qe -= 1; if(Qe<0)
        Error2(819., "n=",n, " bin=",i);
      return 1; }
  return 0;
}

//sorting routines
int o1(int i, int j) {
  double w = T[i]-T[j]; return w<0? -1: w>0? 1: 0;
}

int EventNext() {
  int j;

  PINIT;                                     //Initialize if necessary.

  while(Qe>0)
  { for(; Qi<Qn; Qo=0,Qi++)                  //Advance to the next non-empty
    { j = Q[Qi]; if(j==0) continue;          //bin.

      if(Qo==0)                              //Sort the bin if it may be
        j = Q[Qi] = sort(P,j,0,o1), Qo = 1;  //necessary (usually sorts 1 or 2).

      if(T[j]<Qt1)                           //If the event belongs to this
      { if(P[j]==PEMPTY) Error(820.1);       //pass, remove it from the list,
        Q[Qi] = P[j];                        //decrement the number of events,
        P[j] = PEMPTY; Qe -= 1;              //advance the global time, and
        t = T[j]; return j; } }              //return the event's index.

    Qi = 0; Qt0 += Qw; Qt1 = Qt0+Qw; }       //Circle back to the first bin.

  return 0;                                  //Signal completion of all events.
}

void EventRenumber(int n, int m) {
  if(n<1||n>=PN) Error1(734.3, "n=",n);      //Check the indexes and make sure
  if(m<1||m>=PN) Error1(734.4, "n=",m);      //they are in range.

  if(n!=m)
  { T[n] = T[m];                             //Transfer the time.
    EventCancel(m);                          //Cancel the old number.
    EventSchedule(n, T[n]); }                //Reschedule as the new number.
}

int sort(int list[], int p, int n, int (*c)(int,int)) {
  int i;

  Pl = list;                                  //Record the location of the list
  order = c;                                 //and of the ordering function.

  if(n==0) for(i=p; i; i=Pl[i],n++);          //Count the number of elements and
  if(n==0||p==0) return 0;                   //return empty lists immediately.

  if(n==1) return p;                         //Do the same for single elements.

  if(n==2)                                   //If the list contains only two
  { if(order(p, Pl[p])<=0) return p;          //elements, sort it by inspection
    i = Pl[p]; Pl[i] = p; Pl[p] = 0;            //(One and two element lists are
    return i; }                              //the most common in applications
                                             //like hashing.)

  pcurr = p;                                 //Otherwise sort the list with the
  return isort(n);                           //full algorithm.
}

int isort(int n) {
  int wp1, wp2, count1;

// SORT A SINGLE ELEMENT

  if(n<=1)                                   //If a single element has been
  { if(pcurr==0) return 0;                   //requested, initialize variables
    wp1 = pcurr; count = 0;                  //and check for error in count.

    do                                       //Scan forward in the list to
    { pprev = pcurr; count += 1;             //find the longest string that
      if ((pcurr=Pl[pcurr])==0) return wp1; } //is already in order.
    while(order(pprev, pcurr)<=0);

    Pl[pprev] = 0;                            //Break this string from the
    return wp1; }                            //rest of the list and return.

// SORT MORE THAN ONE ELEMENT

  wp1 = isort(n/2);                          //Sort the first part of the list
  if(n<=count) return wp1;                   //and return if the required
  count1 = count;                            //amount was fortuitously sorted.

  wp2 = isort(n-count);                      //If it was not, then sort what
  count += count1;                           //remains to be sorted on this
  return imerge(wp1, wp2);                   //call and merge the two sublists.
}

int imerge(int p, int q) {
 int pbeg;

 if(p==0) return q;                          //Handle cases where one or more
 if(q==0) return p;                          //of the lists is null.

 pbeg = p;                                   //Save an index to the beginning
 if(order(p,q)<=0) goto scanp;               //of the list and enter the
 pbeg = q;                                   //proper scanning routine.

 scans:                                      //Scan for a secondary element
   do { pprev = q;                           //greater than or equal to the
        if ((q=Pl[q])==0) goto ends; }        //current primary element and mend
   while(order(p,q)>0);                      //the secondary list.
   Pl[pprev] = p;

 scanp:                                      //Scan for a primary element
   do { pprev = p;                           //greater than the current
        if ((p=Pl[p])==0) goto endp; }        //secondary element, mend the
   while(order(p,q)<=0);                     //primary list, and repeat.
   Pl[pprev] = q; goto scans;

 endp: p = q;                                //Attach any remaining elements,
 ends: Pl[pprev] = p;                         //which need not be scanned and
   return pbeg;                              //return the merged list.
}