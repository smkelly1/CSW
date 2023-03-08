#include <math.h>
#include "csw.h"

void calc_divergence(void)
{
	int i, j, m, n;
	double cos01, cos1, cos12;
	double gamma;
	double cm2, cn2, phi_m, phi_n;
	double tmp2;
	
	#if defined(KAPPA) || defined(FILE_KAPPA)
		double kappa0; 
		double invDX, invDX2;
		invDX=1/DX;
		invDX2=1/(DX*DX);
		
		double dpdy[2], dpdx2, dpdy2;
	#endif

	#ifdef CORIOLIS
		double wave_res;
	#endif
	
	#if defined(R) || defined(FILE_R)
		double r0;
	#endif

	////////////////////////////////////////////////////////////////////
	// Average velocity in time
	for(n=0; n<NM; n++){
		for(j=0; j<NY+2; j++){
			for(i=0; i<NX+2; i++){
				U1[n][j][i]=(U[n][j][i]+U1[n][j][i])/2;
				V1[n][j][i]=(V[n][j][i]+V1[n][j][i])/2;
			}
		}
	}

	////////////////////////////////////////////////////////////////////
	// Start pressure forcing 
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
					Fp[n][j][i]=-c[n][j+1][i+1]*c[n][j+1][i+1]
						*((U1[n][j+1][i+2]-U1[n][j+1][i+1])
						+(cos12*V1[n][j+2][i+1]-cos01*V1[n][j+1][i+1]))*gamma/H[j+1][i+1];
					
					// Topographic coupling 
					#ifdef MODECOUPLE
						if (H[j+1][i+1]>H_MIN_COUPLE) {
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
									
									tmp2=((U1[m][j+1][i+1]+U1[m][j+1][i+2])/2*dHdx[j][i]
										 +(V1[m][j+1][i+1]+V1[m][j+2][i+1])/2*dHdy[j][i])
										 *cm2/(cn2-cm2)*phi_m*phi_n
										 *cn2/(H[j+1][i+1]*H[j+1][i+1]);
									
									if (isfinite(tmp2)) {
										Fp[n][j][i]=Fp[n][j][i]+tmp2;
									}
								}
							}
						}
					#endif           

					// Compute wave resolution for damping
					#ifdef CORIOLIS
						wave_res=fabs(2*c[n][j+1][i+1]/(f[j+1]*DX*A*cos1));													
					#endif

					// Diffusion
					#if defined(KAPPA) || defined(FILE_KAPPA)
						// File values over-rides constant value
						#ifdef FILE_KAPPA
							kappa0=kappa[n][j+1][i+1];
							#ifdef GAMMA
								kappa0=GAMMA*kappa0;
							#endif
						#else
							kappa0=KAPPA;
						#endif
						
						// Increase diffusion and add linear-damping in under-resolved regions 
						//#ifdef CORIOLIS
							//if (wave_res<2) {
								//kappa0=(KAPPA_MAX-kappa0)*(2-wave_res)/2+kappa0; 	
							//}							
						//#endif
						
						dpdx2=1/(A*A*cos1*cos1)*(p[n][j+1][i+2]-2*p[n][j+1][i+1]+p[n][j+1][i]);	

						dpdy[0]=cos01*(p[n][j+1][i+1]-p[n][j][i+1]);
						dpdy[1]=cos12*(p[n][j+2][i+1]-p[n][j+1][i+1]);
						dpdy2=1/(A*A*cos1)*(dpdy[1]-dpdy[0]);

						Fp[n][j][i]=Fp[n][j][i]+kappa0*(dpdx2+dpdy2)*invDX2;						
					#endif
					
					// Linear drag
					#if defined(R) || defined(FILE_R)
					
						// File values override constant value
						#ifdef FILE_R
							r0=r[n][j+1][i+1];
						#else
							r0=fmin(R/(c[n][j+1][i+1]*c[n][j+1][i+1]/2),R_MAX);
						#endif
						
						#ifdef CORIOLIS
							r0=fmin(r0,fabs(f[j+1])); // Nothing larger than f
							if (wave_res<2) {
								r0=fmax(fabs(f[j+i])*(2-wave_res)/2,r0); // Increase damping for low resolution	
							}							
						#endif
						
						Fp[n][j][i]=Fp[n][j][i]-r0*p[n][j+1][i+1];
					#endif

				} // end-if: land mask

			}
		}
	}
}
