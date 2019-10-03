#include <math.h>
#include "csw.h"

void timestep_p(void)
{ 
	int i, j, n;
	double gamma;
	
	#ifdef R_MASK
		double c_crit, mask;
	#endif

	gamma=DT/12;

	for(n=0; n<NM; n++){
		for(j=0; j<NY; j++){

			// Slowest resolvable wave on the grid (Adcroft et al. 1999) 
			#ifdef R_MASK
				c_crit=fabs(f[j+1])*A*cos(lat[j+1])*DX;  
			#endif

			for(i=0; i<NX; i++){

				if (H[j+1][i+1]>H_MIN) {

					// Create a sponge with an R_MASK decay time scale where there is insufficient wave resolution
					#ifdef R_MASK
						mask=1-fmax(0,(1-c[n][j+1][i+1]/c_crit)*DT*R_MASK);  
					#endif
					
					// Save the current data
					p1[n][j+1][i+1]=p[n][j+1][i+1];
					
					// Integrate to find the new data
					#ifdef R_MASK
						p[n][j+1][i+1]=(p1[n][j+1][i+1]+gamma*(23*Fp[n][j][i]-16*Fp_1[n][j][i]+5*Fp_2[n][j][i]))*mask;
					#else
						p[n][j+1][i+1]=p1[n][j+1][i+1]+gamma*(23*Fp[n][j][i]-16*Fp_1[n][j][i]+5*Fp_2[n][j][i]);
					#endif

					// Old forcing becomes very old forcing 
					Fp_2[n][j][i]=Fp_1[n][j][i];

					// Current forcing becomes old forcing 
					Fp_1[n][j][i]=Fp[n][j][i];

					// Average to find the mid-point data
					p1[n][j+1][i+1]=(p[n][j+1][i+1]+p1[n][j+1][i+1])/2;

				}
				
			}
		}
	} 
}
