#include <math.h>
#include "csw.h"

void calc_forces(void)
{
	int i, j, m, n;
	double cos1, cn2, cm2, phi_n, phi_m, tmp2;

	#if (defined(NU) || defined(NU_MAX) || defined(FILE_NU))
		double cos0, cos01, cos12;
		double nuX[2], nuY[2], Fdx[2], Fdy[2];
	#endif

	#if (defined(CORIOLIS) || defined(CD))
		double f01;
		double wC, wSW, wNW, wSE, wNE;
	#endif

	#if (defined(R) || defined(R_MAX) || defined(FILE_R))
		double r0;
	#endif

	#ifdef CD
		double Hu, Hv, absU, absV;
	#endif

	#ifdef DAMP_GROWTH
		double alpha;
	#endif 

	////////////////////////////////////////////////////////////////////
	// Compute provisional pressure 
	for(n=0; n<NM; n++){
		for(j=0; j<NY+2; j++){
			for(i=0; i<NX+2; i++){				
				pE[n][j][i]=(0.5+GAMMA+2.0*EPSILON)*p[n][j][i]+(0.5-2.0*GAMMA-3.0*EPSILON)*p1[n][j][i]+GAMMA*p2[n][j][i]+EPSILON*p3[n][j][i];
				
				// Shift old velocities to make room for new velocity calculations 
				U3[n][j][i]=U2[n][j][i];
				U2[n][j][i]=U1[n][j][i];
				U1[n][j][i]=U[n][j][i];	
				
				V3[n][j][i]=V2[n][j][i];
				V2[n][j][i]=V1[n][j][i];
				V1[n][j][i]=V[n][j][i];				
			}
		}
	}

	////////////////////////////////////////////////////////////////////
	// Identify problem areas
	#ifdef DAMP_GROWTH				
		for(j=0; j<NY+2; j++){
			for(i=0; i<NX+2; i++){		
				// Get total SSH
				eta[j][i]=0;
				for(n=0; n<NM; n++){
					eta[j][i]=eta[j][i]+p[n][j][i]*phi_surf[n][j][i]/9.81;
				}
							
				// Run the low pass filter three times to get rid of all secular Growth
				alpha=DT/(NUM_PERIODS*2*M_PI/fmax(fabs(f[j]),2.5325E-05)); // Use the local inertial frequency, or f(10), which ever is higher. 
				eta1[j][i]=alpha*eta[j][i]+(1-alpha)*eta1[j][i];
				eta2[j][i]=alpha*(eta[j][i]-eta1[j][i])+(1-alpha)*eta2[j][i];
				eta3[j][i]=alpha*(eta[j][i]-eta1[j][i]-eta2[j][i])+(1-alpha)*eta3[j][i];
				
				// Low-pass the growth so that diffusivity adjustments are smooth 
				etaG[j][i]=alpha*(eta1[j][i]+eta2[j][i]+eta3[j][i])+(1-alpha)*etaG[j][i];
				
				// Find the amplitude of the waves (pi/2 relates the average to the amplitude). 
				etaA[j][i]=alpha*fabs(0.5*M_PI*(eta[j][i]-eta1[j][i]-eta2[j][i]-eta3[j][i]))+(1-alpha)*etaA[j][i]; 
				
				// Get a metric for how much extra damping should happen
				flag_growth[j][i]=fmax((etaA[j][i]-AMP_THRESHOLD)/(AMP_MAX-AMP_THRESHOLD),
										(fabs(etaG[j][i])-GROWTH_THRESHOLD)/(GROWTH_MAX-GROWTH_THRESHOLD)); 	
				flag_growth[j][i]=fmax(0,flag_growth[j][i]); // Negative values are below the thresholds
				flag_growth[j][i]=fmin(flag_growth[j][i],1); // Values greater than one exceed max growth
				
				// set nu and kappa to increase linearly from background values to MAX over the SSH range THRESHOLD to MAX
				nu[j][i]=NU+flag_growth[j][i]*(NU_MAX-NU);
				kappa[j][i]=KAPPA+flag_growth[j][i]*(KAPPA_MAX-KAPPA);
			
			} // end-for: i
		} // end-for: j
	#endif

	////////////////////////////////////////////////////////////////////
	// Forces in U-direction
	for(j=0; j<NY; j++){

		cos1=cos(lat[j+1]);
		
		#if (defined(NU) || defined(NU_MAX) || defined(FILE_NU))
			cos01=cos((lat[j]+lat[j+1])/2);
			cos12=cos((lat[j+1]+lat[j+2])/2);
		#endif

		for(i=0; i<NX; i++){				
			if (H[j+1][i]>H_MIN && H[j+1][i+1]>H_MIN) {
				for(n=0; n<NM; n++){

					// Pressure gradient force 
					Fu[n][j][i]=-(H[j+1][i]+H[j+1][i+1])/2*(pE[n][j+1][i+1]-pE[n][j+1][i])/(A*cos1*DX);
						
					// Coriolis
					#ifdef CORIOLIS
						if (0<fabs(f[j+1])){
							wC=sqrt(fabs(f[j+1]))/fmax((H[j+1][i]+H[j+1][i+1])/2,1);
							wSW=sqrt(fabs((f[j]+f[j+1])/2))/fmax((H[j][i]+H[j+1][i])/2,1);
							wNW=sqrt(fabs((f[j+1]+f[j+2])/2))/fmax((H[j+1][i]+H[j+2][i])/2,1);
							wSE=sqrt(fabs((f[j]+f[j+1])/2))/fmax((H[j][i+1]+H[j+1][i+1])/2,1);
							wNE=sqrt(fabs((f[j+1]+f[j+2])/2))/fmax((H[j+1][i+1]+H[j+2][i+1])/2,1);
							Fu[n][j][i]=Fu[n][j][i]+f[j+1]/(4*wC)*(wSW*VE[n][j+1][i]+wNW*VE[n][j+2][i]+wSE*VE[n][j+1][i+1]+wNE*VE[n][j+2][i+1]);
							//Fu[n][j][i]=Fu[n][j][i]+f[j+1]*(VE[n][j+1][i]+VE[n][j+2][i]+VE[n][j+1][i+1]+VE[n][j+2][i+1])/4;
						}
					#endif

					// Wind forcing
					#ifdef WIND_FORCING
						Fu[n][j][i]=Fu[n][j][i]+(tau_x[j+1][i]+tau_x[j+1][i+1])*(phi_surf[n][j+1][i]+phi_surf[n][j+1][i+1])/(4*RHO);
					#endif

					// Topographic coupling, see Zaron et al. (2020; JPO)
					#ifdef MODECOUPLE
						if (H[j+1][i]>H_MIN_COUPLE && H[j+1][i+1]>H_MIN_COUPLE) {
							cn2=(c[n][j+1][i]+c[n][j+1][i+1])*(c[n][j+1][i]+c[n][j+1][i+1])/4;
							phi_n=(phi_bott[n][j+1][i]+phi_bott[n][j+1][i+1])/2;
				
							for(m=0; m<NM; m++){
								if (m==n) {
									Fu[n][j][i]=Fu[n][j][i]-
										(1-phi_n*phi_n)/2*dHdx_u[j][i]
										*(pE[n][j+1][i]+pE[n][j+1][i+1])/2;
								}
								else {
									cm2=(c[m][j+1][i]+c[m][j+1][i+1])*(c[m][j+1][i]+c[m][j+1][i+1])/4;
									phi_m=(phi_bott[m][j+1][i]+phi_bott[m][j+1][i+1])/2;

									tmp2=-cn2/(cm2-cn2)*phi_m*phi_n*dHdx_u[j][i]
										*(pE[m][j+1][i]+pE[m][j+1][i+1])/2;
									
									if (isfinite(tmp2)) {
										Fu[n][j][i]=Fu[n][j][i]+tmp2;
									}
								}
							}
						}
					#endif 

					////////////////////////////////////////////////////
					// Frictional forces
					Fu_eps[n][j][i]=0;

					// Horizontal diffusion
					#if (defined(NU) || defined(DAMP_GROWTH) || defined(FILE_NU))
											
						#if defined(FILE_NU) || defined(DAMP_GROWTH)
							nuX[0]=nu[j+1][i];
							nuX[1]=nu[j+1][i+1];
							nuY[0]=(nu[j][i]+nu[j+1][i]+nu[j][i+1]+nu[j+1][i+1])/4;
							nuY[1]=(nu[j+1][i]+nu[j+2][i]+nu[j+1][i+1]+nu[j+2][i+1])/4;
						#else
							nuX[0]=fmin(NU,NU_MAX);
							nuX[1]=fmin(NU,NU_MAX);
							nuY[0]=fmin(NU,NU_MAX);
							nuY[1]=fmin(NU,NU_MAX);
						#endif						
						
						// Compute diffusive fluxes
						Fdx[0]=-nuX[0]*(UE[n][j+1][i+1]-UE[n][j+1][i])/(A*cos1*DX);
						Fdx[1]=-nuX[1]*(UE[n][j+1][i+2]-UE[n][j+1][i+1])/(A*cos1*DX);

						Fdy[0]=-nuY[0]*(UE[n][j+1][i+1]-UE[n][j][i+1])/(A*DX);
						Fdy[1]=-nuY[1]*(UE[n][j+2][i+1]-UE[n][j+1][i+1])/(A*DX);
						
						// Compute diffusive flux convergence
						Fu_eps[n][j][i]=Fu_eps[n][j][i]-((Fdx[1]-Fdx[0])+(cos12*Fdy[1]-cos01*Fdy[0]))/(A*cos1*DX);											
					#endif

					// Linear drag
					#if (defined(R) || defined(R_MAX) || defined(FILE_R))
						r0=0;

						#ifdef FILE_R
							r0=fmin((r[n][j+1][i]+r[n][j+1][i+1])/2,R_MAX);
						#endif						

						#ifdef R
							r0=R;
						#endif					

						#ifdef RC2
							r0=fmin(R/((c[n][j+1][i]+c[n][j+1][i+1])/2*(c[n][j+1][i]+c[n][j+1][i+1])/2),R_MAX);
						#endif
					
						#ifdef DAMP_GROWTH 	
							if (((etaA[j+1][i]+etaA[j+1][i+1])/2>AMP_THRESHOLD) || ((etaG[j+1][i]+etaG[j+1][i+1])/2>GROWTH_THRESHOLD)) {
								r0=r0+fmax(
									((etaA[j+1][i]+etaA[j+1][i+1])/2-AMP_THRESHOLD)*(R_MAX-r0)/(AMP_MAX-AMP_THRESHOLD),
									(fabs(etaG[j+1][i]+etaG[j+1][i+1])/2-GROWTH_THRESHOLD)*(R_MAX-r0)/(GROWTH_MAX-GROWTH_THRESHOLD)); 
							}											
							
							// Use NU_MAX and don't count dissipation
							if (r0>R_MAX){
								Fu[n][j][i]=Fu[n][j][i]-R_MAX*UE[n][j+1][i+1];					
								r0=0;
							}							
						#endif						
												
						Fu_eps[n][j][i]=Fu_eps[n][j][i]-r0*UE[n][j+1][i+1];					
					#endif

					// Quadratic bottom drag
					#ifdef CD
						Hu=(H[j+1][i]+H[j+1][i+1])/2;
						wSW=1/fmax((H[j][i]+H[j+1][i])/2,1);
						wNW=1/fmax((H[j+1][i]+H[j+2][i])/2,1);
						wSE=1/fmax((H[j][i+1]+H[j+1][i+1])/2,1);
						wNE=1/fmax((H[j+1][i+1]+H[j+2][i+1])/2,1);					
						absU=cabs(UE[n][j+1][i+1]+I*(Hu/4)*(wSW*VE[n][j+1][i]+wNW*VE[n][j+2][i]+wSE*VE[n][j+1][i+1]+wNE*VE[n][j+2][i+1]));
						Fu_eps[n][j][i]=Fu_eps[n][j][i]-CD*absU*UE[n][j+1][i+1]/(Hu*Hu);
					#endif

					// Add the frictional forces to the other forces
					Fu[n][j][i]=Fu[n][j][i]+Fu_eps[n][j][i];

				} // end m loop				
			} // end-if: land mask
		} // end i loop	
	} // end j loop	


	////////////////////////////////////////////////////////////////////
	// Forces in V-direction
	for(j=0; j<NY; j++){

		#if (defined(NU) || defined(FILE_NU))
			cos01=cos((lat[j]+lat[j+1])/2);
			cos0=cos(lat[j]);      
			cos1=cos(lat[j+1]);
		#endif

		#ifdef CORIOLIS
			f01=(f[j]+f[j+1])/2;
		#endif

		for(i=0; i<NX; i++){
			if (H[j][i+1]>H_MIN && H[j+1][i+1]>H_MIN) {
				for(n=0; n<NM; n++){

					// Pressure gradient force	
					Fv[n][j][i]=-(H[j][i+1]+H[j+1][i+1])/2*(pE[n][j+1][i+1]-pE[n][j][i+1])/(A*DX);
													
					// Coriolis
					#ifdef CORIOLIS
						if (0<fabs(f01)) {
							wC=sqrt(fabs(f01))/fmax((H[j][i+1]+H[j+1][i+1])/2,1);
							wSW=sqrt(fabs(f[j]))/fmax((H[j][i]+H[j][i+1])/2,1);	
							wSE=sqrt(fabs(f[j]))/fmax((H[j][i+1]+H[j][i+2])/2,1);
							wNW=sqrt(fabs(f[j+1]))/fmax((H[j+1][i]+H[j+1][i+1])/2,1);
							wNE=sqrt(fabs(f[j+1]))/fmax((H[j+1][i+1]+H[j+1][i+2])/2,1);
							Fv[n][j][i]=Fv[n][j][i]-f01/(4*wC)*(wSW*UE[n][j][i+1]+wSE*UE[n][j][i+2]+wNW*UE[n][j+1][i+1]+wNE*UE[n][j+1][i+2]);
							//Fv[n][j][i]=Fv[n][j][i]-f01*(UE[n][j][i+1]+UE[n][j][i+2]+UE[n][j+1][i+1]+UE[n][j+1][i+2])/4;
						}
					#endif

					// Wind forcing
					#ifdef WIND_FORCING
						Fv[n][j][i]=Fv[n][j][i]+(tau_y[j][i+1]+tau_y[j+1][i+1])*(phi_surf[n][j][i+1]+phi_surf[n][j+1][i+1])/(4*RHO);
					#endif
				
					// Topographic coupling, see Zaron et al. (2020; JPO) 
					#ifdef MODECOUPLE
						if (H[j][i+1]>H_MIN_COUPLE && H[j+1][i+1]>H_MIN_COUPLE) {
							cn2=(c[n][j][i+1]+c[n][j+1][i+1])*(c[n][j][i+1]+c[n][j+1][i+1])/4;
							phi_n=(phi_bott[n][j][i+1]+phi_bott[n][j+1][i+1])/2;
				
							for(m=0; m<NM; m++){
								if (m==n) {
									Fv[n][j][i]=Fv[n][j][i]-
										(1-phi_n*phi_n)/2*dHdy_v[j][i]
										*(pE[n][j][i+1]+pE[n][j+1][i+1])/2;									
								}
								else {
									cm2=(c[m][j][i+1]+c[m][j+1][i+1])*(c[m][j][i+1]+c[m][j+1][i+1])/4;
									phi_m=(phi_bott[m][j][i+1]+phi_bott[m][j+1][i+1])/2;

									tmp2=-cn2/(cm2-cn2)*phi_m*phi_n*dHdy_v[j][i]
										*(pE[m][j][i+1]+pE[m][j+1][i+1])/2;

									if (isfinite(tmp2)) {
										Fv[n][j][i]=Fv[n][j][i]+tmp2;
									}
								}
							}													
						}
					#endif
					
					////////////////////////////////////////////////////
					// Frictional forces
					Fv_eps[n][j][i]=0;

					// Horizontal diffusion
					#if (defined(NU) || defined(DAMP_GROWTH) || defined(FILE_NU))					
						
						#if defined(FILE_NU) || defined(DAMP_GROWTH)
							nuX[0]=(nu[j][i]+nu[j+1][i]+nu[j][i+1]+nu[j+1][i+1])/4;
							nuX[1]=(nu[j+1][i]+nu[j+2][i]+nu[j+1][i+1]+nu[j+2][i+1])/4;
							nuY[0]=nu[j][i+1];
							nuY[1]=nu[j+1][i+1];
						#else
							nuX[0]=fmin(NU,NU_MAX);
							nuX[1]=fmin(NU,NU_MAX);
							nuY[0]=fmin(NU,NU_MAX);
							nuY[1]=fmin(NU,NU_MAX);
						#endif						
						
						// Compute diffusive fluxes
						Fdx[0]=-nuX[0]*(VE[n][j+1][i+1]-VE[n][j+1][i])/(A*cos01*DX);
						Fdx[1]=-nuX[1]*(VE[n][j+1][i+2]-VE[n][j+1][i+1])/(A*cos01*DX);

						Fdy[0]=-nuY[0]*(VE[n][j+1][i+1]-VE[n][j][i+1])/(A*DX);
						Fdy[1]=-nuY[1]*(VE[n][j+2][i+1]-VE[n][j+1][i+1])/(A*DX);
						
						// Compute diffusive flux convergence
						Fv_eps[n][j][i]=Fv_eps[n][j][i]-((Fdx[1]-Fdx[0])+(cos1*Fdy[1]-cos0*Fdy[0]))/(A*cos01*DX);
					#endif

					// Linear drag
					#if (defined(R) || defined(R_MAX) || defined(FILE_R))
						r0=0;
						
						#ifdef FILE_R
							r0=fmin((r[n][j][i+1]+r[n][j+1][i+1])/2,R_MAX);
						#endif
						
						#ifdef R
							r0=R;
						#endif	
						
						#ifdef RC2
							r0=fmin(R/((c[n][j][i+1]+c[n][j+1][i+1])/2*(c[n][j][i+1]+c[n][j+1][i+1])/2),R_MAX);
						#endif
						
						#ifdef DAMP_GROWTH 	// If killing growth, use R_MAX and don't count drag in F_eps							
							if (((etaA[j][i+1]+etaA[j+1][i+1])/2>AMP_THRESHOLD) || ((etaG[j][i+1]+etaG[j+1][i+1])/2>GROWTH_THRESHOLD)) {
								r0=r0+fmax(
									((etaA[j][i+1]+etaA[j+1][i+1])/2-AMP_THRESHOLD)*(R_MAX-r0)/(AMP_MAX-AMP_THRESHOLD),
									(fabs(etaG[j][i+1]+etaG[j+1][i+1])/2-GROWTH_THRESHOLD)*(R_MAX-r0)/(GROWTH_MAX-GROWTH_THRESHOLD)); 
							}											
							
							// Use NU_MAX and don't count dissipation
							if (r0>R_MAX){
								Fv[n][j][i]=Fv[n][j][i]-R_MAX*VE[n][j+1][i+1];						
								r0=0;
							}
						#endif					
					
						Fv_eps[n][j][i]=Fv_eps[n][j][i]-r0*VE[n][j+1][i+1];						
					#endif

					// Quadratic bottom drag
					#ifdef CD
						Hv=(H[j][i+1]+H[j+1][i+1])/2;
						wSW=1/fmax((H[j][i]+H[j][i+1])/2,1);	
						wSE=1/fmax((H[j][i+1]+H[j][i+2])/2,1);
						wNW=1/fmax((H[j+1][i]+H[j+1][i+1])/2,1);
						wNE=1/fmax((H[j+1][i+1]+H[j+1][i+2])/2,1);
						absV=cabs((Hv/4)*(wSW*UE[n][j][i+1]+wSE*UE[n][j][i+2]+wNW*UE[n][j+1][i+1]+wNE*UE[n][j+1][i+2])+I*VE[n][j+1][i+1]);
						Fv_eps[n][j][i]=Fv_eps[n][j][i]-CD*absV*VE[n][j+1][i+1]/(Hv*Hv);
					#endif

					// Add the frictional forces to the other forces
					Fv[n][j][i]=Fv[n][j][i]+Fv_eps[n][j][i];

				} // end m-loop				
			} // end-if: land mask
		} // end i-loop	
	} // end j-loop	

	////////////////////////////////////////////////////////////////////
	// Time step U and V	
	for(j=0; j<NY; j++){			
		for(i=0; i<NX; i++){					
			for(n=0; n<NM; n++){					
				U[n][j+1][i+1]=U[n][j+1][i+1]+DT*Fu[n][j][i];					
				V[n][j+1][i+1]=V[n][j+1][i+1]+DT*Fv[n][j][i];
			} // end n-loop			
		} // end i-loop	
	} // end j-loop

}
