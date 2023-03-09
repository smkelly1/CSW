#include <math.h>
#include "csw.h"

void calc_divergence(void)
{
	int i, j, m, n;
	double cos01, cos1, cos12;
	double gamma;
	double cm2, cn2, phi_m, phi_n;
	double tmp2;
	
	#ifdef TIDE_FORCING
		double complex phase;
	#endif 
	
	////////////////////////////////////////////////////////////////////
	// Extrapolate U and V to t+dt/2 using 3rd Order Adams Bashforth 
	for(n=0; n<NM; n++){
		for(j=0; j<NY+2; j++){
			for(i=0; i<NX+2; i++){
				// Staggered AB3			
				//U1[n][j][i]=(U[n][j][i]+U1[n][j][i])/2;
				//UE[n][j][i]=(23*U1[n][j][i]-16*U2[n][j][i]+5*U3[n][j][i])/12;

				//V1[n][j][i]=(V[n][j][i]+V1[n][j][i])/2;
				//VE[n][j][i]=(23*V1[n][j][i]-16*V2[n][j][i]+5*V3[n][j][i])/12;
												
				// ROMS
				UE[n][j][i]=(1.5+BETA)*U[n][j][i]-(0.5+2.0*BETA)*U1[n][j][i]+BETA*U2[n][j][i];
				VE[n][j][i]=(1.5+BETA)*V[n][j][i]-(0.5+2.0*BETA)*V1[n][j][i]+BETA*V2[n][j][i];
			}
		}
	}

	////////////////////////////////////////////////////////////////////
	// Start pressure forcing (divergence & topographic coupling)
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
						*((UE[n][j+1][i+2]-UE[n][j+1][i+1])
						+(cos12*VE[n][j+2][i+1]-cos01*VE[n][j+1][i+1]))*gamma/H[j+1][i+1];
					
					// Compute these only once
					#if defined(MODECOUPLE) || defined(TIDAL_FORCING)
						cn2=c[n][j+1][i+1]*c[n][j+1][i+1];
						phi_n=phi_bott[n][j+1][i+1];
					#endif
					
					// Topographic coupling, see Zaron et al. (2020; JPO)
					#ifdef MODECOUPLE
						if (H[j+1][i+1]>H_MIN_COUPLE) {															
							for(m=0; m<NM; m++){								
								if (m==n) {
									Fp[n][j][i]=Fp[n][j][i]+
										((UE[n][j+1][i+1]+UE[n][j+1][i+2])/2*dHdx[j][i]
										+(VE[n][j+1][i+1]+VE[n][j+2][i+1])/2*dHdy[j][i])
										*0.5*(1-phi_n*phi_n)
										*cn2/(H[j+1][i+1]*H[j+1][i+1]);
								}
								else {
									cm2=c[m][j+1][i+1]*c[m][j+1][i+1];
									phi_m=phi_bott[m][j+1][i+1];
									
									tmp2=((UE[m][j+1][i+1]+UE[m][j+1][i+2])/2*dHdx[j][i]
										 +(VE[m][j+1][i+1]+VE[m][j+2][i+1])/2*dHdy[j][i])
										 *cm2/(cn2-cm2)*phi_m*phi_n
										 *cn2/(H[j+1][i+1]*H[j+1][i+1]);
									
									if (isfinite(tmp2)) {
										Fp[n][j][i]=Fp[n][j][i]+tmp2;
									}
								}
							}
						}
					#endif 
					
					// Barotropic tidal forcing
					#ifdef TIDAL_FORCING
						if(H[j+1][i+1]>H_MIN_FORCE) {							
							for(m=0; m<NC; m++){  // Cycle through tidal frequencies 
								phase=cexp(I*ITGF.omega[m]*(t+DT/2)); // Compute forcing at t+dt/2 				
								Fp[n][j][i]=Fp[n][j][i]-cn2*phi_n*creal(ITGF.F[m][j][i]*phase);
							}
						}
					#endif				

				} // end-if: land mask

			}
		}
	}
	
	////////////////////////////////////////////////////////////////////
	// Time step pressure 
	for(n=0; n<NM; n++){
		for(j=0; j<NY+2; j++){
			for(i=0; i<NX+2; i++){
				p3[n][j][i]=p2[n][j][i];
				p2[n][j][i]=p1[n][j][i];
				p1[n][j][i]=p[n][j][i];		
			}
		}
	}
	
	for(n=0; n<NM; n++){
		for(j=0; j<NY; j++){
			for(i=0; i<NX; i++){
				if (H[j+1][i+1]>H_MIN && -70<(lat[j+1]*180/M_PI)) {									
					// Integrate to find the new data			
					p[n][j+1][i+1]=p[n][j+1][i+1]+DT*Fp[n][j][i];
					//p[n][j+1][i+1]=p[n][j+1][i+1]+gamma*(55*Fp[n][j][i]-59*Fp_1[n][j][i]+37*Fp_2[n][j][i]-9*Fp_3[n][j][i]);
					//p[n][j+1][i+1]=p[n][j+1][i+1]+gamma*(23*Fp[n][j][i]-16*Fp_1[n][j][i]+5*Fp_2[n][j][i]);

					// Old forcing becomes very old forcing 
					//Fp_3[n][j][i]=Fp_2[n][j][i];
					//Fp_2[n][j][i]=Fp_1[n][j][i];				
					//Fp_1[n][j][i]=Fp[n][j][i];	
				}
			}
		}
	}
	
}
