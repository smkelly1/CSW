#include <math.h>
#include "csw.h"

void timestep_uv(void)
{ 
	int i, j, n;
	double gamma;
	
	#ifdef HIGH_PASS
		double alpha;
	#endif
	
	gamma=DT/12;
		
	for(n=0; n<NM; n++){
		for(j=0; j<NY; j++){			
			for(i=0; i<NX; i++){					
				
				// Update U
				if (H[j+1][i]>H_MIN && H[j+1][i+1]>H_MIN && -65<(lat[j+1])*180/M_PI) {	
	
					// Save the current data
					U1[n][j+1][i+1]=U[n][j+1][i+1];
	
					// Integrate to find the new data
					U[n][j+1][i+1]=U[n][j+1][i+1]+gamma*(23*Fu[n][j][i]-16*Fu_1[n][j][i]+5*Fu_2[n][j][i]);
					
					// Old forcing becomes very old forcing 
					Fu_2[n][j][i]=Fu_1[n][j][i];
	
					// Current forcing becomes old forcing 
					Fu_1[n][j][i]=Fu[n][j][i];		
					
					// Remove the exponential running average
					#ifdef HIGH_PASS
						alpha=fabs(DT/(NUM_PERIODS*2*M_PI/f[j]));
						U_low1[n][j][i]=alpha*U[n][j][i]+(1-alpha)*U_low1[n][j][i];
						U_low2[n][j][i]=alpha*(U[n][j][i]-U_low1[n][j][i])+(1-alpha)*U_low2[n][j][i];
						U_low3[n][j][i]=alpha*(U[n][j][i]-U_low1[n][j][i]-U_low2[n][j][i])+(1-alpha)*U_low3[n][j][i];
						U[n][j][i]=U[n][j][i]-U_low1[n][j][i]-U_low2[n][j][i]-U_low3[n][j][i];	
					#endif		
				
				}
				
				// Update V
				if (H[j][i+1]>H_MIN && H[j+1][i+1]>H_MIN && -65<(lat[j+1])*180/M_PI) {

					// Save the current data
					V1[n][j+1][i+1]=V[n][j+1][i+1];
					
					// Integrate to find the new data
					V[n][j+1][i+1]=V[n][j+1][i+1]+gamma*(23*Fv[n][j][i]-16*Fv_1[n][j][i]+5*Fv_2[n][j][i]);				
					
					// Old forcing becomes very old forcing 
					Fv_2[n][j][i]=Fv_1[n][j][i];
					
					// Current forcing becomes old forcing 
					Fv_1[n][j][i]=Fv[n][j][i];			
										
					// Remove the exponential running average
					#ifdef HIGH_PASS
						alpha=fabs(DT/(NUM_PERIODS*2*M_PI/f[j]));
						V_low1[n][j][i]=alpha*V[n][j][i]+(1-alpha)*V_low1[n][j][i];
						V_low2[n][j][i]=alpha*(V[n][j][i]-V_low1[n][j][i])+(1-alpha)*V_low2[n][j][i];
						V_low3[n][j][i]=alpha*(V[n][j][i]-V_low1[n][j][i]-V_low2[n][j][i])+(1-alpha)*V_low3[n][j][i];
						V[n][j][i]=V[n][j][i]-V_low1[n][j][i]-V_low2[n][j][i]-V_low3[n][j][i];	
					#endif
								
				}
				
			}
		}		
	} // end n-loop
		
}
