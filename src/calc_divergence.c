#include <math.h>
#include "csw.h"

void calc_divergence(void)
{
	int i, j, m, n;
	double cos01, cos1, cos12;
	double cm2, cn2, phi_m, phi_n, tmp2;
	
	#if (defined(R) || defined(R_MAX) || defined(FILE_R))
		double r0;
	#endif
	
	#if (defined(KAPPA) || defined(KAPPA_MAX) || defined(FILE_KAPPA))
		double kappaX[2], kappaY[2], Fdx[2], Fdy[2];
	#endif
	
	////////////////////////////////////////////////////////////////////
	// Extrapolate U and V to t+dt/2 using 3rd Order Adams Bashforth 
	for(n=0; n<NM; n++){
		for(j=0; j<NY+2; j++){
			for(i=0; i<NX+2; i++){										
				UE[n][j][i]=(1.5+BETA)*U[n][j][i]-(0.5+2.0*BETA)*U1[n][j][i]+BETA*U2[n][j][i];
				VE[n][j][i]=(1.5+BETA)*V[n][j][i]-(0.5+2.0*BETA)*V1[n][j][i]+BETA*V2[n][j][i];
				#if defined(R) || defined(R_MAX) || defined(KAPPA)
					pE[n][j][i]=(1.5+BETA)*p[n][j][i]-(0.5+2.0*BETA)*p1[n][j][i]+BETA*p2[n][j][i]; 
				#endif	
				
				// Shift old pressures to make room for new p calculation
				p3[n][j][i]=p2[n][j][i];
				p2[n][j][i]=p1[n][j][i];
				p1[n][j][i]=p[n][j][i];						
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

		for(i=0; i<NX; i++){
			if (H[j+1][i+1]>H_MIN) {				
				for(n=0; n<NM; n++){

					// Volume divergence
					Fp[n][j][i]=-((UE[n][j+1][i+2]-UE[n][j+1][i+1])
						+(cos12*VE[n][j+2][i+1]-cos01*VE[n][j+1][i+1]))/(A*cos1*DX);
					
					// Need this for several terms
					cn2=c[n][j+1][i+1]*c[n][j+1][i+1];
					phi_n=phi_bott[n][j+1][i+1];
					
					// Topographic coupling, see Zaron et al. (2020; JPO)
					#ifdef MODECOUPLE
						if (H[j+1][i+1]>H_MIN_COUPLE) {															
							for(m=0; m<NM; m++){								
								if (m==n) {
									Fp[n][j][i]=Fp[n][j][i]+
										((UE[n][j+1][i+1]+UE[n][j+1][i+2])/2*dHdx[j][i]
										+(VE[n][j+1][i+1]+VE[n][j+2][i+1])/2*dHdy[j][i])
										*0.5*(1-phi_n*phi_n)/H[j+1][i+1];
								}
								else {
									cm2=c[m][j+1][i+1]*c[m][j+1][i+1];
									phi_m=phi_bott[m][j+1][i+1];
									
									tmp2=((UE[m][j+1][i+1]+UE[m][j+1][i+2])/2*dHdx[j][i]
										 +(VE[m][j+1][i+1]+VE[m][j+2][i+1])/2*dHdy[j][i])
										 *cm2/(cn2-cm2)*phi_m*phi_n/H[j+1][i+1];
									
									if (isfinite(tmp2)) {
										Fp[n][j][i]=Fp[n][j][i]+tmp2;
									}
								}
							}
						}
					#endif 
					
					// Barotropic tidal forcing (at t+DT/2)
					#ifdef TIDE_FORCING
						for(m=0; m<NC; m++){  // Cycle through tidal frequencies 
							Fp[n][j][i]=Fp[n][j][i]-phi_n*creal(ITGF.F[m][j][i]*cexp(I*ITGF.omega[m]*(t+DT/2)));
						}
					#endif				

					// Multiply the pressure forcing by c_n^2/H
					Fp[n][j][i]=Fp[n][j][i]*cn2/H[j+1][i+1];


					////////////////////////////////////////////////////
					// Frictional forces (do not divide by c_n^2/H)
					Fp_eps[n][j][i]=0;

					// Diffusion
					#if (defined(KAPPA) || defined(DAMP_GROWTH) || defined(FILE_KAPPA))				
											
						#if defined(FILE_KAPPA) || defined(DAMP_GROWTH)
							kappaX[0]=(kappa[j+1][i]+kappa[j+1][i+1])/2;
							kappaX[1]=(kappa[j+1][i+1]+kappa[j+1][i+2])/2;
							kappaY[0]=(kappa[j][i+1]+kappa[j+1][i+1])/2;
							kappaY[1]=(kappa[j+1][i+1]+kappa[j+2][i+1])/2;
						#else
							kappaX[0]=fmin(KAPPA,KAPPA_MAX);
							kappaX[1]=fmin(KAPPA,KAPPA_MAX);
							kappaY[0]=fmin(KAPPA,KAPPA_MAX);
							kappaY[1]=fmin(KAPPA,KAPPA_MAX);
						#endif
					
						// Compute diffusive fluxes
						Fdx[0]=-kappaX[0]*(pE[n][j+1][i+1]-pE[n][j+1][i])/(A*cos1*DX);
						Fdx[1]=-kappaX[1]*(pE[n][j+1][i+2]-pE[n][j+1][i+1])/(A*cos1*DX);

						Fdy[0]=-kappaY[0]*(pE[n][j+1][i+1]-pE[n][j][i+1])/(A*DX);
						Fdy[1]=-kappaY[1]*(pE[n][j+2][i+1]-pE[n][j+1][i+1])/(A*DX);
						
						// Compute diffusive flux convergence
						Fp_eps[n][j][i]=Fp_eps[n][j][i]-((Fdx[1]-Fdx[0])+(cos12*Fdy[1]-cos01*Fdy[0]))/(A*cos1*DX);											
					#endif

					// Linear drag
					#if (defined(R) || defined(R_MAX) || defined(FILE_R))
									
						#ifdef FILE_R
							r0=fmin(r[n][j+1][i+1],R_MAX);
						#endif
						
						#ifdef R
							r0=R;
						#endif	
						
						#ifdef RC2
							r0=fmin(R/cn2,R_MAX);
						#endif
						
						#ifdef DAMP_GROWTH							
							// set r0 to increase linearly from background value to R_MAX over the SSH range THRESHOLD to MAX 
							if ((etaA[j+1][i+1]>AMP_THRESHOLD) || (fabs(etaG[j+1][i+1])>GROWTH_THRESHOLD)) {
								r0=r0+fmax(
									(etaA[j+1][i+1]-AMP_THRESHOLD)*(R_MAX-r0)/(AMP_MAX-AMP_THRESHOLD),
									(fabs(etaG[j+1][i+1])-GROWTH_THRESHOLD)*(R_MAX-r0)/(GROWTH_MAX-GROWTH_THRESHOLD)); 
							} 			
														
							// If r0 is bigger than R_MAX, use R_MAX and don't count dissipation
							if (r0>R_MAX){ 
								Fp[n][j][i]=Fp[n][j][i]-R_MAX*pE[n][j+1][i+1];
								r0=0;
							}
						#endif					
					
						Fp_eps[n][j][i]=Fp_eps[n][j][i]-r0*pE[n][j+1][i+1];						
					#endif
					
					// Add the frictional forces to the other forces
					Fp[n][j][i]=Fp[n][j][i]+Fp_eps[n][j][i];
				
				} // end n loop 
			} // end-if: land mask
		} // end i loop
	} // end j loop
	
	
	////////////////////////////////////////////////////////////////////
	// Time step pressure 
	for(j=0; j<NY; j++){
		for(i=0; i<NX; i++){
			for(n=0; n<NM; n++){				
				p[n][j+1][i+1]=p[n][j+1][i+1]+DT*Fp[n][j][i];
			} // end-for: n		
		} // end-for: i
	} // end-for: j
	
}
