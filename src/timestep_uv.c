#include <math.h>
#include "csw.h"

void timestep_uv(void)
{ 
	int i, j, n;
	double gamma;
	
	gamma=DT/24;
		
	for(n=0; n<NM; n++){
		for(j=0; j<NY; j++){			
			for(i=0; i<NX; i++){					
				
				// Update U
				if (H[j+1][i]>H_MIN && H[j+1][i+1]>H_MIN && -70<(lat[j+1])*180/M_PI) {	
	
					// Save the current data
					U1[n][j+1][i+1]=U[n][j+1][i+1];
	
					// Integrate to find the new data
					U[n][j+1][i+1]=U[n][j+1][i+1]+gamma*(55*Fu[n][j][i]-59*Fu_1[n][j][i]+37*Fu_2[n][j][i]-9*Fu_3[n][j][i]);

					// Old forcing becomes very old forcing 
					Fu_3[n][j][i]=Fu_2[n][j][i];
					Fu_2[n][j][i]=Fu_1[n][j][i];
					Fu_1[n][j][i]=Fu[n][j][i];	
				}
				
				// Update V
				if (H[j][i+1]>H_MIN && H[j+1][i+1]>H_MIN && -70<(lat[j+1])*180/M_PI) {

					// Save the current data
					V1[n][j+1][i+1]=V[n][j+1][i+1];
					
					// Integrate to find the new data
					V[n][j+1][i+1]=V[n][j+1][i+1]+gamma*(55*Fv[n][j][i]-59*Fv_1[n][j][i]+37*Fv_2[n][j][i]-9*Fv_3[n][j][i]);
					
					// Old forcing becomes very old forcing 
					Fv_3[n][j][i]=Fv_2[n][j][i];
					Fv_2[n][j][i]=Fv_1[n][j][i];
					Fv_1[n][j][i]=Fv[n][j][i];
				}
				
			}
		}		
	} // end n-loop
			
}
