#include <math.h>
#include "csw.h"

void calc_diagnostics(int Na)
{
	int i, j, m, n;
	double cos01, cos1, cos12, cn2, cm2, phi_n, phi_m;

	#ifdef WORK
		double dFdx, dFdy;
	#endif

	#ifdef WRITE_TRANSPORT
		float tmpf;
	#endif 

	////////////////////////////////////////////////////////////////////
	// Use t-dt/2 variables 
	for(n=0; n<NM; n++){
		for(j=0; j<NY+2; j++){
			for(i=0; i<NX+2; i++){
				#ifdef AB4
					UE[n][j][i]=(U[n][j][i]+UE[n][j][i])/2;	
					VE[n][j][i]=(V[n][j][i]+VE[n][j][i])/2;	
					pE[n][j][i]=p[n][j][i];
				#else
					UE[n][j][i]=(U[n][j][i]+U1[n][j][i])/2;	
					VE[n][j][i]=(V[n][j][i]+V1[n][j][i])/2;		
					pE[n][j][i]=(p[n][j][i]+p1[n][j][i])/2;		
				#endif
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

			for(i=0; i<NX; i++){
				if (H[j+1][i+1]>H_MIN) {
				#ifdef DAMP_GROWTH
					if (flag_growth[j+1][i+1]<1) {
				#endif
					for(n=0; n<NMW; n++){

						#ifdef ENERGY
							KE[n][j][i]=KE[n][j][i]+(float)(0.125*RHO*((UE[n][j+1][i+1]+UE[n][j+1][i+2])*(UE[n][j+1][i+1]+UE[n][j+1][i+2])
								+(VE[n][j+1][i+1]+VE[n][j+2][i+1])*(VE[n][j+1][i+1]+VE[n][j+2][i+1]))/H[j+1][i+1]); // the extra factor of 1/4 comes from averaging U^2
							PE[n][j][i]=PE[n][j][i]+(float)(0.5*RHO*H[j+1][i+1]*pE[n][j+1][i+1]*pE[n][j+1][i+1]/(c[n][j+1][i+1]*c[n][j+1][i+1]));
						#endif

						#ifdef FLUX
							up[n][j][i]=up[n][j][i]+(float)(0.5*RHO*(UE[n][j+1][i+1]+UE[n][j+1][i+2])*pE[n][j+1][i+1]);
							vp[n][j][i]=vp[n][j][i]+(float)(0.5*RHO*(VE[n][j+1][i+1]+VE[n][j+2][i+1])*pE[n][j+1][i+1]);
						#endif

						#ifdef WORK
							// Energy-flux divergence
							dFdx=(UE[n][j+1][i+2]*(pE[n][j+1][i+1]+pE[n][j+1][i+2])-UE[n][j+1][i+1]*(pE[n][j+1][i]+pE[n][j+1][i+1]));
							dFdy=(cos12*VE[n][j+2][i+1]*(pE[n][j+1][i+1]+pE[n][j+2][i+1])-cos01*VE[n][j+1][i+1]*(pE[n][j][i+1]+pE[n][j+1][i+1]));
							divF[n][j][i]=divF[n][j][i]+(float)(0.5*RHO*(dFdx+dFdy)/(A*cos1*DX)); // the factor of 1/2 comes from averaging p

							// Wind forcing
							#ifdef WIND_FORCING
								W[n][j][i]=W[n][j][i]+(float)(0.5*(
									tau_x[j+1][i+1]*(UE[n][j+1][i+1]+UE[n][j+1][i+2])
									+tau_y[j+1][i+1]*(VE[n][j+1][i+1]+VE[n][j+2][i+1]))
									*phi_surf[n][j+1][i+1]/H[j+1][i+1]);	
									// The factor of 1/2 comes from averaging U
							#endif

							// Compute dissipation
							D[n][j][i]=D[n][j][i]-(float)(RHO*0.25*((Fu_eps[n][j][i]+Fu_eps[n][j][i+1])*(UE[n][j+1][i+1]+UE[n][j+1][i+2])
								+(Fv_eps[n][j][i]+Fv_eps[n][j+1][i])*(VE[n][j+1][i+1]+VE[n][j+2][i+1]))/H[j+1][i+1])
								-(float)(RHO*H[j+1][i+1]*Fp_eps[n][j][i]*pE[n][j+1][i+1]/(c[n][j+1][i+1]*c[n][j+1][i+1]));	

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
												*pE[n][j+1][i+1]
												*((UE[m][j+1][i+1]+UE[m][j+1][i+2])*dHdx[j][i]+(VE[m][j+1][i+1]+VE[m][j+2][i+1])*dHdy[j][i])
												-cn2/(cm2-cn2)*phi_n*phi_m
												*pE[m][j+1][i+1]
												*((UE[n][j+1][i+1]+UE[n][j+1][i+2])*dHdx[j][i]+(VE[n][j+1][i+1]+VE[n][j+2][i+1])*dHdy[j][i]))
												/H[j+1][i+1]); 
												//the factor of 1/2 comes from averaging U and V
										}										
									}									
								}							
							#endif // end MODECOUPLE 
						#endif // end WORK 
						
					} // end m loop	
				#ifdef DAMP_GROWTH
					} // growth mask
				#endif			
				} // end-if: land mask
			} // end i loop	
		} // end j loop	

	#endif // end diagnostics 

	////////////////////////////////////////////////////////////////////
	// Estimate SSH amplitude and phase
	#ifdef WRITE_SSH
		for(n=0; n<NMW; n++){
			for(j=0; j<NY; j++){
				for(i=0; i<NX; i++){
					if (p[n][j+1][i+1]>SSH_amp[n][j][i]) {
						SSH_amp[n][j][i]=p[n][j+1][i+1];
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

					tmpf=(float)((U[n][j+1][i+1]+U[n][j+1][i+2])/2);
					if (tmpf>U_amp[n][j][i]) {
						U_amp[n][j][i]=tmpf;
						U_phase[n][j][i]=Na;
					}

					tmpf=(float)((V[n][j+1][i+1]+V[n][j+2][i+1])/2);
					if (tmpf>V_amp[n][j][i]) {
						V_amp[n][j][i]=tmpf;
						V_phase[n][j][i]=Na;
					}
				}
			}
		}
	#endif

}
