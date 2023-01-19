#include <math.h>
#include "csw.h"

void calc_divergence(void)
{
	int i, j, m, n;
	double cos01, cos1, cos12;
	double gamma;

	#ifdef MODECOUPLE
		double cm2, cn2, phi_m, phi_n;
	#endif

	#if defined(KAPPA) || defined(FILE_KAPPA)
		double kappa0, invDX2;
		invDX2=1/(DX*DX);
		
		#ifdef SPHERE
			double dpdy[2], dpdx2, dpdy2;
		#endif
	#endif

	#if defined(R) || defined(FILE_R)
		double r0;
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
						if (H[j+1][i+1]>H_MIN_COUPLE) {	
							
							#ifdef NO_ANTARCTIC
								if(lat[j+1]>-60*M_PI/180) { // Recall: lat is in radians
							#endif	
									// Note: See Zaron et al. (2020; JPO)
									cn2=c[n][j+1][i+1]*c[n][j+1][i+1];
									phi_n=phi_bott[n][j+1][i+1];
										
									for(m=0; m<NM; m++){								
										if (m==n) {
											Fp[n][j][i]=Fp[n][j][i]+
												((U1[n][j+1][i+1]+U1[n][j+1][i+2])/2*dHdx[j][i]
												+(V1[n][j+1][i+1]+V1[n][j+2][i+1])/2*dHdy[j][i])
												*0.5*(1-phi_n*phi_n)
												*cn2/(H[j+1][i+1]*H[j+1][i+1]);
										}
										else {
											cm2=c[m][j+1][i+1]*c[m][j+1][i+1];
											phi_m=phi_bott[m][j+1][i+1];

											Fp[n][j][i]=Fp[n][j][i]+
												((U1[m][j+1][i+1]+U1[m][j+1][i+2])/2*dHdx[j][i]
												+(V1[m][j+1][i+1]+V1[m][j+2][i+1])/2*dHdy[j][i])
												*cm2/(cn2-cm2)*phi_m*phi_n
												*cn2/(H[j+1][i+1]*H[j+1][i+1]);
										}
									}
							#ifdef NO_ANTARCTIC
								}
							#endif
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

						#ifdef GAMMA
							kappa0=GAMMA*kappa0;
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
					
					// Linear drag
					#if defined(R) || defined(FILE_R)
						// File values override constant value
						#ifdef FILE_R
							r0=r[n][j+1][i+1];
						#else
							r0=R;
						#endif

						Fp[n][j][i]=Fp[n][j][i]-r0*p[n][j+1][i+1];
					#endif

				} // end-if: land mask

			}
		}
	}
}
