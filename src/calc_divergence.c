#include <math.h>
#include "csw.h"

void calc_divergence(void)
{
	int i, j, n;
	double cos01, cos1, cos12;
	double gamma;
	
	#ifdef MODECOUPLE
		int m;
	#endif

	
	for(j=0; j<NY; j++){

		// Only compute cos once for each latitude
		cos01=cos((lat[j]+lat[j+1])/2);
		cos1=cos(lat[j+1]);
		cos12=cos((lat[j+1]+lat[j+2])/2);
	
		gamma=1/(A*cos1*DX);
		
		for(n=0; n<NM; n++){
			for(i=0; i<NX; i++){
						
				if (H[j+1][i+1]>H_MIN) {
					
					// Volume divergence
					#ifdef SPHERE
						Fp[n][j][i]=-c[n][j+1][i+1]*c[n][j+1][i+1]
							*((U1[n][j+1][i+2]-U1[n][j+1][i+1])
							+(cos12*V1[n][j+2][i+1]-cos01*V1[n][j+1][i+1]))*gamma/H[j+1][i+1];	
					#else	
						Fp[n][j][i]=-c[n][j+1][i+1]*c[n][j+1][i+1]/H[j+1][i+1]
							*((U1[n][j+1][i+2]-U1[n][j+1][i+1])
							+(V1[n][j+2][i+1]-V1[n][j+1][i+1]))/DX;
					#endif
									
					// Topographic coupling
					#ifdef MODECOUPLE

						// Note: This code assumes T is a depth integral (not depth average), hence the division by H^2.
						for(m=0; m<NM; m++){
							Fp[n][j][i]=Fp[n][j][i]+
								(T.x[m][n][j+1][i+1]*(U1[m][j+1][i+1]+U1[m][j+1][i+2])
								+T.y[m][n][j+1][i+1]*(V1[m][j+1][i+1]+V1[m][j+2][i+1]))
								*c[n][j+1][i+1]*c[n][j+1][i+1]/(2*H[j+1][i+1]*H[j+1][i+1]);
						}							
										
					#endif
										
				} // end-if: land mask
				
			}
		}
	}	
}
