#include <math.h>
#include "csw.h"

void calc_divergence(void)
{
	int i, j, m, n;
	double cos01, cos1, cos12;
	double gamma;

	#if defined(KAPPA) || defined(FILE_KAPPA)
		double kappa0, invDX2;
		invDX2=1/(DX*DX);
		
		#ifdef SPHERE
			double dpdy[2], dpdx2, dpdy2;
		#endif
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

					// Diffusion
					#if defined(KAPPA) || defined(FILE_KAPPA)
						// File values over-rides constant value
						#ifdef FILE_KAPPA
							kappa0=kappa[n][j+1][i+1];
						#else
							kappa0=KAPPA;
						#endif

						#ifdef SPHERE
							dpdx2=1/(A*A*cos1*cos1)*(p[n][j+1][i+2]-2*p[n][j+1][i+1]+p[n][j+1][i]);	

							dpdy[0]=cos01*(p[n][j+1][i+1]-p[n][j][i+1]);
							dpdy[1]=cos12*(p[n][j+2][i+1]-p[n][j+1][i+1]);
							dpdy2=1/(A*A*cos1)*(dpdy[1]-dpdy[0]);

							Fp[n][j][i]=Fp[n][j][i]+kappa0*(dpdx2+dpdy2)*invDX2;
						#else
							Fp[n][j][i]=Fp[n][j][i]+kappa0*((p[n][j+1][i+2]-2*p[n][j+1][i+1]+p[n][j+1][i])
								+(p[n][j+2][i+1]-2*p[n][j+1][i+1]+p[n][j][i+1]))*invDX2;
						#endif
					#endif

				} // end-if: land mask

			}
		}
	}
}
