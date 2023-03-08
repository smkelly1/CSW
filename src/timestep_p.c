#include <math.h>
#include "csw.h"

void timestep_p(void)
{ 
	int i, j, n;
	double gamma;
	
	#ifdef FLAG_GROWTH
		double alpha;
	#endif

	gamma=DT/24;

	for(n=0; n<NM; n++){
		for(j=0; j<NY; j++){
			for(i=0; i<NX; i++){

				if (H[j+1][i+1]>H_MIN && -70<(lat[j+1])*180/M_PI) {
					
					// Save the current data
					p1[n][j+1][i+1]=p[n][j+1][i+1];
					
					// Integrate to find the new data			
					p[n][j+1][i+1]=p[n][j+1][i+1]+gamma*(55*Fp[n][j][i]-59*Fp_1[n][j][i]+37*Fp_2[n][j][i]-9*Fp_3[n][j][i]);

					// Old forcing becomes very old forcing 
					Fp_3[n][j][i]=Fp_2[n][j][i];
					Fp_2[n][j][i]=Fp_1[n][j][i];				
					Fp_1[n][j][i]=Fp[n][j][i];					
					
					// Identify any exponential growth (don't remove, but use this field to flag bad data later)
					#ifdef FLAG_GROWTH
						alpha=fabs(DT/(NUM_PERIODS*2*M_PI/f[j]));
						p_low[n][j][i]=alpha*p[n][j][i]+(1-alpha)*p_low[n][j][i];
					#endif	
												
				}
				
			}
		}
	} 
}
