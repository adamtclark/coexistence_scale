// SIS - Spatial Invasion Simulator
//
// WS Harpole, v3 rewritten in C 7/1/2007 from R prototype
//	
//	The Spatial Invasion Simulator (SIS) is a spatially explicit, two-species lattice model. To explore alternative restoration and invasion scenarios, the two species are designated as Exotic and Native.
//	Individuals of each species occupy discrete cells on a square lattice. At each time step individuals have some probability of mortality (100% if they are annuals). Unoccupied cells are colonized by either species, depending on their dispersal mode, their abundance, and the positive soil feedback state generated by the previous cell occupant.
//	Dispersal can be either local (nearest 8 neighbors) or global (all individuals in the lattice can send propagules to all cells). Boundaries are absorbing and the lattice is updated synchronously.
//	Positive soil feedback decreases the probability of the other species’ establishment and survival. The strength of the feedback increases with time of occupancy of the species up to the maximum value.
//	The Model assumes the following:
//	- The environment is homogenous (there is no variation in the geography, soil type, etc.) and that species are identical in their environmental requirements.
//	- Net feedback is positive (or zero, i.e., positive feedback effects are greater than or equal to negative feedback effects).
//	- There is no long-term seedbanking – all propagules come from the previous time step.
//
//	Alternative management scenario options include:
//	- Solarizing the soil (setting the initial soil feedback state).
//	- Control of the exotic species (increasing exotic mortality).
//	- Augmenting native seed production (decreasing exotic/native seed ratio).
//	- Including an external source of exotic seeds (one edge of the lattice provides a persistent source of seeds).
//
// 
// updated to run in R programming language with 13 parameter arguments (plus input/output vectors), described in comments below

#include <R.h>
#include <Rmath.h>
#include <stdio.h>

void sis_R (int *psp1dis, int *psp2dis, int *psp1fb, int *psp2fb, int *psp1m, int *psp2m, int *pscenario,
	int *pinitabund, double *pseed, int *pedge, int *pdim, int *ptmax, double *ppr_nocol, int outmat[], int outmap0[], int outmap[])  {
	//sp1dis, sp2dis;		species dispersal: 0=local, 1=global
	//sp1fb, sp2fb;			species feedback strength: 0=none, 100=overwhelmingly strong
	//sp1m, sp2m;			species percent mortality: 0=immortal, <100=perennial, 100=annual
	//scenario;				initial soil feedback state: 0=neutral, 1=exotic (spp1), 2=native (spp2)
	//initabund; 			initial percent natives
	//seed; 				ratio of exotic to native seed production: 1= default (values: 1000, 100, 10, 1, 0.1, 0.01, 0.001)
	//edge;					0=absorbing boundaries, 1=one edge all exotics (spp1)
	//dim;					dimension of matrix (not counting the edges)
	//tmax;					time steps
	//ppr_nocol;			probability reduction in colonization success (allows for failed colonization events)
	//outmat				matrix for storing abundances through time
	//outmap0				initial positions for all individuals
	//outmap 				final positions for all individuals

	// de-pointerize input from R...
	int sp1dis = *psp1dis;
	int sp2dis = *psp2dis;
	int sp1fb = abs(*psp1fb);
	int sp2fb = abs(*psp2fb);
	int sp1m = *psp1m;
	int sp2m = *psp2m;
	int scenario = *pscenario;
	int initabund = *pinitabund;
	double seed = *pseed;
	int edge = *pedge;
	int dim = *pdim;
	int tmax = *ptmax;
	double pr_nocol= *ppr_nocol;
	
	GetRNGstate();  // Load seed for random number generation, R base

	// import user defined parameters
	// Note: current error handling only checks for out-of-bounds values - doesn't change flow control
	//		species 1=exotic, species 2=native
	//			if( argc<14) fprintf( outlog, "ERROR: too few parameters %d\n", argc );
	
	int trigger = 0; // triggers exit if parameter inputs are outside of bounds

	if( !( sp1dis == 0 || sp1dis == 1) ) {trigger=1; 			fprintf(stderr, "ERROR: sp1dis must be 0 or 1");}
	if( !( sp2dis == 0 || sp2dis == 1) ) {trigger=1; 			fprintf(stderr, "ERROR: sp2dis must be 0 or 1");}
	if( ( sp1fb < 0 || sp1fb > 100) ) {trigger=1;	  			fprintf(stderr, "ERROR: sp1fb must be >= -100 and <=100");}
	if( ( sp2fb < 0 || sp2fb > 100) ) {trigger=1;  				fprintf(stderr, "ERROR: sp2fb must be >= -100 and <=100");}
	if( ( sp1m < 0 || sp1m > 100) ) {trigger=1;  				fprintf(stderr, "ERROR: sp1m must be >= 0 and <=100");}
	if( ( sp2m < 0 || sp2m > 100) ) {trigger=1; 				fprintf(stderr, "ERROR: sp2m must be >= 0 and <=100");}
	if(!(scenario==0||scenario==1||scenario==2)) {trigger=1; 	fprintf(stderr, "ERROR: scenario must be 0, 1, or 2");}
	if( ( initabund <= 0 || initabund > 100) ) {trigger=1;  	fprintf(stderr, "ERROR: initabund must be >0 and <=100");}
	if( ( seed < 0.001 || seed > 1000) ) {trigger=1;  			fprintf(stderr, "ERROR: seed must be between 0.001 and 1000");}
	if( !( edge == 0 || edge == 1) ) {trigger=1;  				fprintf(stderr, "ERROR: edge must be 0 or 1]");}
	if( ( dim <= 0 || dim > 1000) ) {trigger=1;  				fprintf(stderr, "ERROR: dim must be > 0 and <=1000]");}
	if( ( tmax <= 0 || tmax > 10000) ) {trigger=1;  			fprintf(stderr, "ERROR: tmax must be > 0 and <=10000");}
	if( ( pr_nocol < 0 || pr_nocol > 1) ) {trigger=1;  			fprintf(stderr, "ERROR: pr_nocol must be >= 0 and <=1");}

	if(trigger == 0) { // only run function if no error
		// declarations of internal variables
		int i, j, k, t;				// matrix and time indexes
		double s1, s2;				// species 1 and 2 survival probabilities
		double c1, c2;				// species 1 and 2 colonization probabilities
		double p1, p2, p3;				// species 1 and 2 establishment probabilities
		int xt0[*pdim+2][*pdim+2];	// time t generation, a dim x dim lattice with 1-cell edge boundaries
		int xt1[*pdim+2][*pdim+2];	// time t+1 generation, a dim x dim lattice with 1-cell edge boundaries
		int near[3];				// counts of 8 nearest neighbors+cell[i][j]: empty, spp1, and spp2
		int all[(*ptmax)+1][3];			// counts of all cells: empty, spp1, and spp2
		int state[*pdim+2][*pdim+2];// cell state, indicating effect of plant-soil feedback
		double tmpp;				// stores random number
		int stepsize = 10;			// size of step to take in soil feedbacks

		// get direction of feedbacks
		int signsp1fb; // sign of feedback
		int signsp2fb; // sign of feedback

		if((*psp1fb)<0) {
			signsp1fb = -1;
		} else {
			signsp1fb = 1;
		}
		if((*psp2fb)<0) {
			signsp2fb = -1;
		} else {
			signsp2fb = 1;
		}


		// initialize species matrix
		for (i =1; i< dim+1 ; ++i) {
			for (j=1; j< dim+1 ; ++j) {
			xt1[i][j] = 1 ; // set to species 1
			if ( (runif(0,1)) < initabund/100. )
				xt1[i][j] = 2; // set to species 2 with P=initial abundance/100
				switch (scenario) {
					case 0:  // initially neutral soils
						state[i][j] = 0; 
						break;
					case 1:  // initially exotic soils (species 1 values are NEGATIVE)
						state[i][j] = -sp1fb; 
						break;
					case 2:  // initially native soils
						state[i][j] = sp2fb; 
						break;
				}
			}
		}
			
		// absorbing edges
		for (i =0; i< dim+2 ; ++i) {
			xt1[0][i]	  =0;
			if( edge==1 && i<dim+1 ) xt1[0][i]=1; // set one edge to all exotics (spp1)
			xt1[dim+1][i] =0;
			xt1[i][0]	  =0;
			xt1[i][dim+1] =0;
		}
		
		// save initial matrix
		for (i =0; i<dim+2 ; ++i) {
			for (j=0; j<dim+2 ; ++j){
				outmap0[i*(dim+2)+j] = xt1[i][j];
			}
		}


		// initialize global abundance vs time array
		for (t =0; t<=tmax; ++t){
			all[t][0]=0;
			all[t][1]=0;
			all[t][2]=0;
		}


		// run model for tmax time steps	
		for( t=0; t< tmax; ++t){
			
			// determine global species abundance and set previous time step = to current
			for (i =0; i< dim+2 ; ++i) {
				for (j=0; j< dim+2 ; ++j) {
					// set previous t0 matrix to current t1 matrix 
					xt0[i][j] = xt1[i][j] ;
					// global abundance			
					if (xt0[i][j]==0){ all[t][0]++ ; }
					if (xt0[i][j]==1){ all[t][1]++ ; }
					if (xt0[i][j]==2){ all[t][2]++ ; }
				}
			}

			// update current matrix (mortality then colonization)
			for (i =1; i<dim+1 ; ++i) {
				for (j=1; j<dim+1 ; ++j){
				
					// count 8 nearest neighbors plus current cell using previous time t0 matrix (empty, spp1 or spp2)
					// BUT only calculate if dispersal is local
					if((sp1dis== 0 && xt0[i][j]==1)||(sp2dis== 0 && xt0[i][j]==2)) { 
						near[0]=0;
						near[1]=0;
						near[2]=0;
						for (k =0; k<9; ++k){
							switch ( xt0[i+(k - 3*(k/3) -1)][j+(k/3-1)] ) { 
								case 0: near[0]++ ; 
								break;
								case 1: near[1]++ ; 
								break;
								case 2: near[2]++ ; 
								break;
							}
						}
					}
					
					// feedback state effect on survival probabilities
					//s1= s2= 1. ;
					//if( state[i][j]>0 ) // positive state value decreases s1, given positive feedback
						//s1= (1-abs(state[i][j])/100.) ;
					//if ( state[i][j]<0 ) // negative state value decreases s2, given positive feedback
						//s2= (1-abs(state[i][j])/100.) ;

					// code follows paper:
					// given POSITIVE feedback, s1 is maximized when state = -100; s2 is maximized when state = 100;
					// given NEGATIVE feedback, s1 is maximized when state = 100; s2 is maximized when state = -100;
					// NOTE: species 1 always makes soils more negative; species 2 always makes soils more positive
					s1= (1-signsp1fb*(state[i][j])/100.)/2 ;  
					s2= (1+signsp2fb*(state[i][j])/100.)/2 ;  
						
					// Species survival
					switch( xt1[i][j] ) {
						case 0:      // cell is empty, move feedback back towards zero
							if( state[i][j] > 0 ) { // decrement cell state in favor of species 1, up to max state 0
								state[i][j] = state[i][j]-stepsize  ; // species 1 has NEGATIVE values
								if(state[i][j] < 0) {
									state[i][j] = 0;
								}
							} else if( state[i][j] < 0 ) { // increment cell state in favor of species 1, up to max state 0
								state[i][j] = state[i][j]+stepsize  ;
								if(state[i][j] > 0) {
									state[i][j] = 0;
								}
							}
						break;

						case 1:      // species 1 occupies cell, m% mortality + additional feedback effect
							if((runif(0,1)) < (sp1m/100. + (1-s1)) )
							{
								xt1[i][j]= 0;
							}
							if( state[i][j] > -sp1fb ) // decrement cell state in favor of species 1, up to (max state -10)
								state[i][j] = state[i][j]-stepsize  ; // species 1 has NEGATIVE values
						break; 
						
						case 2: 	// species 2
							if((runif(0,1)) < (sp2m/100. + (1-s2)) )
								xt1[i][j]= 0;
							if( state[i][j] < sp2fb ) // increment cell state in favor of species 2, up to (max state +10)
								state[i][j] = state[i][j]+stepsize  ; // species 2 has POSITIVE values
						break;
					}

					
					if( xt1[i][j]==0){ 	// empty cell: colonization probability proportional to local or global neighborhood
									//		weighted by ratio of exotic to native seed production
							if(sp1dis== 0)  //local
								c1= seed * near[1]/9. ;
							else  //global
								c1=  seed * (1.* all[t][1]) / (all[t][1] + all[t][2]) ; 
							
							if(sp2dis==0)  //local
								c2= near[2]/9. ;
							else //global
								c2= ( 1.* all[t][2]) / (all[t][1] + all[t][2]) ; 
								
							// probabilities of establishment: Psurvival x Pcolonization
							p1= c1*s1/(c1*s1+c2*s2) ;
							p2= 1-p1 ;
							
							// add probability of no colonization
							p1=(1-pr_nocol)*p1;
							p2=(1-pr_nocol)*p2;
							p3= 1-p1-p2;

							tmpp=runif(0,1);

							if ( tmpp < p1) {
								xt1[i][j] = 1; // set to species 1
							} else if( (tmpp>=p1 & tmpp < (p1+p2)) ) {
								xt1[i][j] = 2; // set to species 2
							} else {
								xt1[i][j] = 0; // set to empty
							}
					}
				}
			}

			R_CheckUserInterrupt();
		}

		// get last time step
		for (i =0; i< dim+2 ; ++i) {
			for (j=0; j< dim+2 ; ++j) {
				// set previous t0 matrix to current t1 matrix 
				xt0[i][j] = xt1[i][j] ;
				// global abundance			
				if (xt0[i][j]==0){ all[t][0]++ ; }
				if (xt0[i][j]==1){ all[t][1]++ ; }
				if (xt0[i][j]==2){ all[t][2]++ ; }
			}
		}

		// save abundance vs. time
		for(t=0; t<=tmax; ++t){
			outmat[t] = all[t][0];
			outmat[(tmax+1)+t] = all[t][1];
			outmat[(tmax+1)*2+t] = all[t][2];
		}

		// save final matrix 
		for (i =0; i<dim+2 ; ++i) {
			for (j=0; j<dim+2 ; ++j){
				//fprintf(stderr, "i = %d ", i);
				outmap[i*(dim+2)+j] = (xt1[i][j]);
			}
		}
	}


	PutRNGstate();  // Unload seed for random number generation, R base
}

