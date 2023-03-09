#include <math.h>
#include "csw.h"

void calc_forces(int Na)
{
	int i, j, m, n;
	double cos0, cos01, cos1, cos12;
	double tmp2;
	double gamma;
	
	double invDX, invDX2;
	invDX=1/DX;
	invDX2=1/(DX*DX);

	double cn2, cm2, phi_n, phi_m;

	#ifdef CORIOLIS
		double wave_res;
	#endif

	#if defined(NU) || defined(FILE_NU)
		double nu0;
		double dUdy[2], dUdx2, dUdy2;
		double dVdy[2], dVdx2, dVdy2;
	#endif

	// Define variables that might be needed
	#if defined(CORIOLIS) || defined(CD) 
		double f01;
		double wC, wSW, wNW, wSE, wNE;
	#endif

	#if defined(R) || defined(FILE_R)
		double r0;
	#endif

	#ifdef CD
		double Hu, Hv, absU, absV;
	#endif

	#ifdef WORK
		double dFdx, dFdy;
		double gamma;

		#if defined(KAPPA) || defined(FILE_KAPPA)
			double kappa0;
			double dpdy[2], dpdx2, dpdy2;
		#endif		
	#endif 

	#ifdef WRITE_TRANSPORT
		float tmp_ave;
	#endif 

	////////////////////////////////////////////////////////////////////
	// Compute provisional pressure 
	for(n=0; n<NM; n++){
		for(j=0; j<NY+2; j++){
			for(i=0; i<NX+2; i++){
				
				// Staggered AB3
				//p1[n][j][i]=(p[n][j][i]+p1[n][j][i])/2;	
				//pE[n][j][i]=(23*p1[n][j][i]-16*p2[n][j][i]+5*p3[n][j][i])/12;
				//UE[n][j][i]=U[n][j][i];
				//VE[n][j][i]=V[n][j][i];
				
				// ROMS
				pE[n][j][i]=(0.5+GAMMA+2.0*EPSILON)*p[n][j][i]+(0.5-2.0*GAMMA-3.0*EPSILON)*p1[n][j][i]+GAMMA*p2[n][j][i]+EPSILON*p3[n][j][i];	
			}
		}
	}

	////////////////////////////////////////////////////////////////////
	// Forces in U-direction
	for(j=0; j<NY; j++){

		cos01=cos((lat[j]+lat[j+1])/2);
		cos1=cos(lat[j+1]);
		cos12=cos((lat[j+1]+lat[j+2])/2);

		for(n=0; n<NM; n++){
			for(i=0; i<NX; i++){
				
				if (H[j+1][i]>H_MIN && H[j+1][i+1]>H_MIN) {
					// Pressure gradient force (computed as flux divergence + topographic coupling term so that dHdx_u can be limited if necessary)
					Fu[n][j][i]=-(H[j+1][i]+H[j+1][i+1])/2*(pE[n][j+1][i+1]-pE[n][j+1][i])*invDX/(A*cos1);
					//Fu[n][j][i]=-(H[j+1][i+1]*pE[n][j+1][i+1]-H[j+1][i]*pE[n][j+1][i])*invDX/(A*cos1)
					//			+dHdx_u[j][i]*(pE[n][j+1][i]+pE[n][j+1][i+1])/2;
					
					// Coriolis
					#ifdef CORIOLIS
						if (0<fabs(f[j+1])){
							wC=sqrt(fabs(f[j+1]))/fmax((H[j+1][i]+H[j+1][i+1])/2,1);
							wSW=sqrt(fabs((f[j]+f[j+1])/2))/fmax((H[j][i]+H[j+1][i])/2,1);
							wNW=sqrt(fabs((f[j+1]+f[j+2])/2))/fmax((H[j+1][i]+H[j+2][i])/2,1);
							wSE=sqrt(fabs((f[j]+f[j+1])/2))/fmax((H[j][i+1]+H[j+1][i+1])/2,1);
							wNE=sqrt(fabs((f[j+1]+f[j+2])/2))/fmax((H[j+1][i+1]+H[j+2][i+1])/2,1);
							Fu[n][j][i]=Fu[n][j][i]+f[j+1]/(4*wC)*(wSW*VE[n][j+1][i]+wNW*VE[n][j+2][i]+wSE*VE[n][j+1][i+1]+wNE*VE[n][j+2][i+1]);
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

					// Compute wave resolution
					#ifdef CORIOLIS
						wave_res=fabs((c[n][j+1][i]+c[n][j+1][i+1])/(f[j+1]*DX*A*cos1));												
					#endif

					// Horizontal diffusion
					#if defined(NU) || defined(FILE_NU)
						#ifdef FILE_NU
							nu0=EFFICIENCY*(nu[n][j+1][i]+nu[n][j+1][i+1])/2;
						#else
							nu0=NU;
						#endif									
						dUdx2=1/(A*A*cos1*cos1)*(UE[n][j+1][i+2]-2*UE[n][j+1][i+1]+UE[n][j+1][i]);	
						dUdy[0]=cos01*(UE[n][j+1][i+1]-UE[n][j][i+1]);
						dUdy[1]=cos12*(UE[n][j+2][i+1]-UE[n][j+1][i+1]);
						dUdy2=1/(A*A*cos1)*(dUdy[1]-dUdy[0]);
						Fu_eps[n][j][i]=Fu_eps[n][j][i]+nu0*(dUdx2+dUdy2)*invDX2;						
					#endif

					// Linear drag
					#if defined(R) || defined(FILE_R)
						#ifdef FILE_R
							r0=(r[n][j+1][i]+r[n][j+1][i+1])/2;
						#else
							r0=fmin(R/((c[n][j+1][i]+c[n][j+1][i+1])/2*(c[n][j+1][i]+c[n][j+1][i+1])/2),R_MAX);
						#endif
						
						#ifdef CORIOLIS
							r0=fmin(r0,fabs(f[j+1])); // Nothing larger than f
							if (wave_res<2) {
								r0=fmax(fabs(f[j+i])*(2-wave_res)/2,r0); // Increase damping for poor resolution	
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

				} // end-if: land mask
				
			}
		}
	}


	////////////////////////////////////////////////////////////////////
	// Forces in V-direction
	for(j=0; j<NY; j++){

		#if (defined(NU) || defined(FILE_NU))
			cos0=cos(lat[j]);      
		#endif

		cos01=cos((lat[j]+lat[j+1])/2);
		cos1=cos(lat[j+1]);

		#ifdef CORIOLIS
			f01=(f[j]+f[j+1])/2;
		#endif

		for(n=0; n<NM; n++){
			for(i=0; i<NX; i++){

				if (H[j][i+1]>H_MIN && H[j+1][i+1]>H_MIN) {
					// Pressure gradient force	(computed as flux divergence + topographic coupling term so that dHdx_u can be limited if necessary)			
					Fv[n][j][i]=-(H[j][i+1]+H[j+1][i+1])/2*(pE[n][j+1][i+1]-pE[n][j][i+1])*invDX/A;
					//Fv[n][j][i]=-(H[j+1][i+1]*pE[n][j+1][i+1]-H[j][i+1]*pE[n][j][i+1])*invDX/A
					//			+dHdy_v[j][i]*(pE[n][j][i+1]+pE[n][j+1][i+1])/2;
								
					// Coriolis
					#ifdef CORIOLIS
						if (0<fabs(f01)) {
							wC=sqrt(fabs(f01))/fmax((H[j][i+1]+H[j+1][i+1])/2,1);
							wSW=sqrt(fabs(f[j]))/fmax((H[j][i]+H[j][i+1])/2,1);	
							wSE=sqrt(fabs(f[j]))/fmax((H[j][i+1]+H[j][i+2])/2,1);
							wNW=sqrt(fabs(f[j+1]))/fmax((H[j+1][i]+H[j+1][i+1])/2,1);
							wNE=sqrt(fabs(f[j+1]))/fmax((H[j+1][i+1]+H[j+1][i+2])/2,1);
							Fv[n][j][i]=Fv[n][j][i]-f01/(4*wC)*(wSW*UE[n][j][i+1]+wSE*UE[n][j][i+2]+wNW*UE[n][j+1][i+1]+wNE*UE[n][j+1][i+2]);
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

					// Wave resolution
					#ifdef CORIOLIS
						wave_res=fabs((c[n][j][i+1]+c[n][j+1][i+1])/(f01*DX*A*cos1));													
					#endif

					// Horizontal diffusion
					#if defined(NU) || defined(FILE_NU)
						#ifdef FILE_NU
							nu0=EFFICIENCY*(nu[n][j][i+1]+nu[n][j+1][i+1])/2;
						#else
							nu0=NU;
						#endif												
						dVdx2=1/(A*A*cos01*cos01)*(VE[n][j+1][i+2]-2*VE[n][j+1][i+1]+VE[n][j+1][i]);
						dVdy[0]=cos0*(VE[n][j+1][i+1]-VE[n][j][i+1]);
						dVdy[1]=cos1*(VE[n][j+2][i+1]-VE[n][j+1][i+1]);
						dVdy2=1/(A*A*cos01)*(dVdy[1]-dVdy[0]);
						Fv_eps[n][j][i]=Fv_eps[n][j][i]+nu0*(dVdx2+dVdy2)*invDX2;						
					#endif

					// Linear drag
					#if defined(R) || defined(FILE_R)
						#ifdef FILE_R
							r0=(r[n][j][i+1]+r[n][j+1][i+1])/2;
						#else
							r0=fmin(R/((c[n][j][i+1]+c[n][j+1][i+1])/2*(c[n][j][i+1]+c[n][j+1][i+1])/2),R_MAX);
						#endif
						
						#ifdef CORIOLIS
							r0=fmin(r0,fabs(f01)); // Nothing larger than f
							if (wave_res<2) {
								r0=fmax(fabs(f01)*(2-wave_res)/2,r0); // Increase damping for poor resolution	
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

				} // end-if: land mask
				
			}
		}
	}


	////////////////////////////////////////////////////////////////////
	// Compute Energy Diagnostics
	#if defined(ENERGY) || defined(FLUX) || defined(WORK) 

		for(j=0; j<NY; j++){

			cos01=cos((lat[j]+lat[j+1])/2);
			cos1=cos(lat[j+1]);
			cos12=cos((lat[j+1]+lat[j+2])/2);

			#ifdef WORK
				gamma=invDX/(A*cos(lat[j+1]));
			#endif

			for(n=0; n<NMW; n++){
				for(i=0; i<NX; i++){

					if (H[j+1][i+1]>H_MIN) {

						#ifdef ENERGY
							KE[n][j][i]=KE[n][j][i]+(float)(0.125*RHO*((U[n][j+1][i+1]+U[n][j+1][i+2])*(U[n][j+1][i+1]+U[n][j+1][i+2])
								+(V[n][j+1][i+1]+V[n][j+2][i+1])*(V[n][j+1][i+1]+V[n][j+2][i+1]))/H[j+1][i+1]); // the extra factor of 1/4 comes from averaging U^2
							PE[n][j][i]=PE[n][j][i]+(float)(0.5*RHO*H[j+1][i+1]*p1[n][j+1][i+1]*p1[n][j+1][i+1]/(c[n][j+1][i+1]*c[n][j+1][i+1]));
						#endif

						#ifdef FLUX
							up[n][j][i]=up[n][j][i]+(float)(0.5*RHO*(U[n][j+1][i+1]+U[n][j+1][i+2])*p1[n][j+1][i+1]);
							vp[n][j][i]=vp[n][j][i]+(float)(0.5*RHO*(V[n][j+1][i+1]+V[n][j+2][i+1])*p1[n][j+1][i+1]);
						#endif

						#ifdef WORK
							// Energy-flux divergence
							dFdx=(U[n][j+1][i+2]*(p1[n][j+1][i+1]+p1[n][j+1][i+2])-U[n][j+1][i+1]*(p1[n][j+1][i]+p1[n][j+1][i+1]));
							dFdy=(cos12*V[n][j+2][i+1]*(p1[n][j+1][i+1]+p1[n][j+2][i+1])-cos01*V[n][j+1][i+1]*(p1[n][j][i+1]+p1[n][j+1][i+1]));
							divF[n][j][i]=divF[n][j][i]+(float)(0.5*RHO*(dFdx+dFdy)*gamma); // the factor of 1/2 comes from averaging p

							// Wind forcing
							#ifdef WIND_FORCING
								W[n][j][i]=W[n][j][i]+(float)(0.5*(
									tau_x[j+1][i+1]*(U[n][j+1][i+1]+U[n][j+1][i+2])
									+tau_y[j+1][i+1]*(V[n][j+1][i+1]+V[n][j+2][i+1]))
									*phi_surf[n][j+1][i+1]/H[j+1][i+1]);	
									// The factor of 1/2 comes from averaging U
							#endif

							// Compute dissipation
							D[n][j][i]=D[n][j][i]-(float)(RHO*0.25*((Fu_eps[n][j][i]+Fu_eps[n][j][i+1])*(U[n][j+1][i+1]+U[n][j+1][i+2])
								+(Fv_eps[n][j][i]+Fv_eps[n][j+1][i])*(V[n][j+1][i+1]+V[n][j+2][i+1]))/H[j+1][i+1]);	

							// Compute diffusion
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
								
								dpdx2=1/(A*A*cos1*cos1)*(p1[n][j+1][i+2]-2*p1[n][j+1][i+1]+p1[n][j+1][i]);	

								dpdy[0]=cos01*(p1[n][j+1][i+1]-p1[n][j][i+1]);
								dpdy[1]=cos12*(p1[n][j+2][i+1]-p1[n][j+1][i+1]);
								dpdy2=1/(A*A*cos1)*(dpdy[1]-dpdy[0]);

								D[n][j][i]=D[n][j][i]-(float)(RHO*H[j+1][i+1]*kappa0*(dpdx2+dpdy2)*invDX2*p1[n][j+1][i+1]/(c[n][j+1][i+1]*c[n][j+1][i+1]));							
							#endif
							
							// Linear drag on pressure
							#if defined(R) || defined(FILE_R)
								// File values override constant value
								#ifdef FILE_R
									r0=r[n][j+1][i+1];
								#else
									r0=R;
								#endif

								D[n][j][i]=D[n][j][i]+(float)(RHO*H[j+1][i+1]*r0*p1[n][j+1][i+1]*p1[n][j+1][i+1]/(c[n][j+1][i+1]*c[n][j+1][i+1]));
							#endif

							// Only compute scattering if there is mode coupling.
							#ifdef MODECOUPLE					
								if (H[j+1][i+1]>H_MIN_COUPLE) {
									cn2=c[n][j+1][i+1]*c[n][j+1][i+1];
									phi_n=phi_bott[n][j+1][i+1];
										
									for(m=0; m<NM; m++){									
										if (m!=n){
											cm2=c[m][j+1][i+1]*c[m][j+1][i+1];
											phi_m=phi_bott[m][j+1][i+1];
							
											Cn[n][j][i]=Cn[n][j][i]
												+(float)(0.5*RHO*(
												cm2/(cn2-cm2)*phi_m*phi_n
												*p1[n][j+1][i+1]
												*((U[m][j+1][i+1]+U[m][j+1][i+2])*dHdx[j][i]+(V[m][j+1][i+1]+V[m][j+2][i+1])*dHdy[j][i])
												-cn2/(cm2-cn2)*phi_n*phi_m
												*p1[m][j+1][i+1]
												*((U[n][j+1][i+1]+U[n][j+1][i+2])*dHdx[j][i]+(V[n][j+1][i+1]+V[n][j+2][i+1])*dHdy[j][i]))
												/H[j+1][i+1]); 
												//the factor of 1/2 comes from averaging U and V
										}										
									}									
								}							
							#endif 

						#endif // end WORK if
					} // land mask
				}
			}
		}

	#endif // end diagnostics if

	////////////////////////////////////////////////////////////////////
	// Estimate SSH amplitude and phase
	#ifdef WRITE_SSH
		for(n=0; n<NMW; n++){
			for(j=0; j<NY; j++){
				for(i=0; i<NX; i++){
					if (p1[n][j+1][i+1]>SSH_amp[n][j][i]) {
						SSH_amp[n][j][i]=p1[n][j+1][i+1];
						SSH_phase[n][j][i]=Na;
					}
				}
			}
		}
	#endif
	
	////////////////////////////////////////////////////////////////////
	// Estimate transport amplitude and phase
	#ifdef WRITE_TRANSPORT
		for(n=0; n<NMW; n++){
			for(j=0; j<NY; j++){
				for(i=0; i<NX; i++){

					tmp_ave=(float)((U[n][j+1][i+1]+U[n][j+1][i+2])/2);
					if (tmp_ave>U_amp[n][j][i]) {
						U_amp[n][j][i]=tmp_ave;
						U_phase[n][j][i]=Na;
					}

					tmp_ave=(float)((V[n][j+1][i+1]+V[n][j+2][i+1])/2);
					if (tmp_ave>V_amp[n][j][i]) {
						V_amp[n][j][i]=tmp_ave;
						V_phase[n][j][i]=Na;
					}
				}
			}
		}
	#endif
	
	
	////////////////////////////////////////////////////////////////////
	// Time step pressure (first shift every point in time, then update the interior points with good forcing	
	for(n=0; n<NM; n++){
		for(j=0; j<NY+2; j++){
			for(i=0; i<NX+2; i++){
				U3[n][j][i]=U2[n][j][i];
				U2[n][j][i]=U1[n][j][i];
				U1[n][j][i]=U[n][j][i];	
				
				V3[n][j][i]=V2[n][j][i];		
				V2[n][j][i]=V1[n][j][i];
				V1[n][j][i]=V[n][j][i];
			}
		}
	}
	
	for(n=0; n<NM; n++){
		for(j=0; j<NY; j++){			
			for(i=0; i<NX; i++){					
				
				// Update U
				if (H[j+1][i]>H_MIN && H[j+1][i+1]>H_MIN && -70<(lat[j+1]*180/M_PI)) {		
					// Integrate to find the new data
					U[n][j+1][i+1]=U[n][j+1][i+1]+DT*Fu[n][j][i];
					//U[n][j+1][i+1]=U1[n][j+1][i+1]+gamma*(23*Fu[n][j][i]-16*Fu_1[n][j][i]+5*Fu_2[n][j][i]);
					//U[n][j+1][i+1]=U[n][j+1][i+1]+gamma*(55*Fu[n][j][i]-59*Fu_1[n][j][i]+37*Fu_2[n][j][i]-9*Fu_3[n][j][i]);

					// Old forcing becomes very old forcing 
					//Fu_3[n][j][i]=Fu_2[n][j][i];
					//Fu_2[n][j][i]=Fu_1[n][j][i];
					//Fu_1[n][j][i]=Fu[n][j][i];	
				}
				
				// Update V
				if (H[j][i+1]>H_MIN && H[j+1][i+1]>H_MIN && -70<(lat[j+1]*180/M_PI)) {
					
					// Integrate to find the new data
					V[n][j+1][i+1]=V[n][j+1][i+1]+DT*Fv[n][j][i];
					//V[n][j+1][i+1]=V1[n][j+1][i+1]+gamma*(23*Fv[n][j][i]-16*Fv_1[n][j][i]+5*Fv_2[n][j][i]);
					//V[n][j+1][i+1]=V[n][j+1][i+1]+gamma*(55*Fv[n][j][i]-59*Fv_1[n][j][i]+37*Fv_2[n][j][i]-9*Fv_3[n][j][i]);
					
					// Old forcing becomes very old forcing 
					//Fv_3[n][j][i]=Fv_2[n][j][i];
					//Fv_2[n][j][i]=Fv_1[n][j][i];
					//Fv_1[n][j][i]=Fv[n][j][i];
				}
				
			}
		}		
	} // end n-loop
	

}
