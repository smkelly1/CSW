#include <math.h>
#include "csw.h"

void calc_divergence(void)
{
	int i, j, m, n;
	double cos01, cos1, cos12;
	double cm2, cn2, phi_m, phi_n, tmp2;

	#ifdef DAMP_GROWTH
		double alpha;
	#endif 	
	
	#if (defined(R) || defined(R_MAX) || defined(FILE_R))
		double r0;
	#endif
	
	#if (defined(KAPPA) || defined(KAPPA_MAX) || defined(FILE_KAPPA))
		double kappa0, dpdy[2], dpdx2, dpdy2;
	#endif
	
	////////////////////////////////////////////////////////////////////
	// Extrapolate U and V to t+dt/2 using 3rd Order Adams Bashforth 
	for(n=0; n<NM; n++){
		for(j=0; j<NY+2; j++){
			for(i=0; i<NX+2; i++){
				#ifdef AB4
					// Staggered AB4
					#if (!defined(ENERGY) && !defined(FLUX) && !defined(WORK) && !defined(WRITE_SSH) && !defined(WRITE_TRANSPORT))	// Otherwise this was already done in calc_diagnostics.c	
						UE[n][j][i]=(U[n][j][i]+UE[n][j][i])/2;
						VE[n][j][i]=(V[n][j][i]+VE[n][j][i])/2;
						#if defined(R) || defined(R_MAX) || defined(KAPPA)
							pE[n][j][i]=p[n][j][i]; 
						#endif
					#endif
				#else							
					// ROMS
					UE[n][j][i]=(1.5+BETA)*U[n][j][i]-(0.5+2.0*BETA)*U1[n][j][i]+BETA*U2[n][j][i];
					VE[n][j][i]=(1.5+BETA)*V[n][j][i]-(0.5+2.0*BETA)*V1[n][j][i]+BETA*V2[n][j][i];
					#if defined(R) || defined(R_MAX) || defined(KAPPA)
						pE[n][j][i]=(1.5+BETA)*p[n][j][i]-(0.5+2.0*BETA)*p1[n][j][i]+BETA*p2[n][j][i]; 
					#endif	
					
					// Shift old pressures to make room for new p calculation
					p3[n][j][i]=p2[n][j][i];
					p2[n][j][i]=p1[n][j][i];
					p1[n][j][i]=p[n][j][i];						
				#endif
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
					#ifdef TIDAL_FORCING
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
					#if (defined(KAPPA) || defined(KAPPA_MAX) || defined(FILE_KAPPA))
						dpdx2=1/(A*A*cos1*cos1)*(pE[n][j+1][i+2]-2*pE[n][j+1][i+1]+pE[n][j+1][i]);	

						dpdy[0]=cos01*(pE[n][j+1][i+1]-pE[n][j][i+1]);
						dpdy[1]=cos12*(pE[n][j+2][i+1]-pE[n][j+1][i+1]);
						dpdy2=1/(A*A*cos1)*(dpdy[1]-dpdy[0]);
						
						kappa0=0;
						
						#ifdef FILE_KAPPA
							kappa0=fmin(GAMMA*kappa[n][j+1][i+1],KAPPA_MAX);
						#endif
						
						#ifdef KAPPA
							kappa0=fmin(KAPPA,KAPPA_MAX);
						#endif

						#ifdef DAMP_GROWTH 	// If killing growth, use KAPPA_MAX and don't count drag in F_eps
							if (0<flag_growth[j+1][i+1]) {
								Fp[n][j][i]=Fp[n][j][i]+KAPPA_MAX*(dpdx2+dpdy2)/(DX*DX);
								kappa0=0;
							}
						#endif

						Fp_eps[n][j][i]=Fp_eps[n][j][i]+kappa0*(dpdx2+dpdy2)/(DX*DX);						
					#endif

					// Linear drag
					#if (defined(R) || defined(R_MAX) || defined(FILE_R))
						r0=0;
						
						#ifdef FILE_R
							r0=fmin(r[n][j+1][i+1],R_MAX);
						#endif
						
						#ifdef R
							r0=fmin(R/cn2,R_MAX);
						#endif
						
						#ifdef DAMP_GROWTH 	// If killing growth, use R_MAX and don't count drag in F_eps
							if (0<flag_growth[j+1][i+1]) {
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
			//if (H[j+1][i+1]>H_MIN && -70<(lat[j+1]*180/M_PI)) {	
				
				for(n=0; n<NM; n++){
					#ifdef AB4
						pE[n][j+1][i+1]=p[n][j+1][i+1];					
						p[n][j+1][i+1]=p[n][j+1][i+1]+DT*(55*Fp[n][j][i]-59*Fp1[n][j][i]+37*Fp2[n][j][i]-9*Fp3[n][j][i])/24;
						Fp3[n][j][i]=Fp2[n][j][i];
						Fp2[n][j][i]=Fp1[n][j][i];				
						Fp1[n][j][i]=Fp[n][j][i];	
					#else								
						p[n][j+1][i+1]=p[n][j+1][i+1]+DT*Fp[n][j][i];						
					#endif				
				}
				
				#ifdef DAMP_GROWTH				
					// Sum the modes
					eta[j][i]=0;
					for(n=0; n<NM; n++){
						eta[j][i]=eta[j][i]+p[n][j+1][i+1]*phi_surf[n][j][i]/9.81;
					}
								
					// Low pass filter (triple running average)
					alpha=fabs(DT/(NUM_PERIODS*2*M_PI/f[j]));
					eta1[j][i]=alpha*eta[j][i]+(1-alpha)*eta1[j][i];
					eta2[j][i]=alpha*(eta[j][i]-eta1[j][i])+(1-alpha)*eta2[j][i];
					eta3[j][i]=alpha*(eta[j][i]-eta1[j][i]-eta2[j][i])+(1-alpha)*eta3[j][i];
					
					// Flag the baddies 
					if (0==flag_growth[j+1][i+1] && GROWTH_THRESHOLD<fabs(eta1[j][i]+eta2[j][i]+eta3[j][i])) {
						flag_growth[j+0][i+0]=1;
						flag_growth[j+1][i+0]=1;
						flag_growth[j+2][i+0]=1;
						flag_growth[j+0][i+1]=1;
						flag_growth[j+1][i+1]=1;
						flag_growth[j+2][i+1]=1;
						flag_growth[j+0][i+2]=1;
						flag_growth[j+1][i+2]=1;
						flag_growth[j+2][i+2]=1;
						if (0<i) {
							flag_growth[j+0][i-1]=1;
							flag_growth[j+1][i-1]=1;
							flag_growth[j+2][i-1]=1;
						}
						if (i<(NX-1)) {
							flag_growth[j+0][i+3]=1;
							flag_growth[j+1][i+3]=1;
							flag_growth[j+2][i+3]=1;
						}
						if (0<j) {
							flag_growth[j-1][i+0]=1;
							flag_growth[j-1][i+1]=1;
							flag_growth[j-1][i+2]=1;
						}
						if (j<(NY-1)) {
							flag_growth[j+3][i+0]=1;
							flag_growth[j+3][i+1]=1;
							flag_growth[j+3][i+2]=1;
						}
					}
					// Let the bad boys back out to play after rehab? Not so fast, most growth is in problem areas.
					//if (0<flag_growth[j+1][i+1] && fabs(eta1[j][i]+eta2[j][i]+eta3[j][i])<0.01*GROWTH_THRESHOLD) {
					//	flag_growth[j+1][i+1]=0;
					//}					
				#endif
				
			//} // end-if: land mask
		} // end-for: i
	} // end-for: j
	
}
