#include <math.h>
#include "csw.h"

void calc_forces(int Na)
{
	int i, j, n;
	double Hu, Hv;
	double cos01, cos1, cos12;
	
	double invDX;
	invDX=1/DX;
	
	// Define variables that might be needed
	#ifdef MODECOUPLE
		int m;
	#endif
	
	#ifdef CORIOLIS
		double f01;
	#endif
		
	#if defined(AX) && defined(SPHERE)
		double dUdy[2], dUdx2, dUdy2;
		double dVdy[2], dVdx2, dVdy2;
	#endif
	
	#ifdef CD
		double absU, absV;
	#endif
	
	#ifdef WORK
		double dFdx, dFdy;
		double gamma;
	#endif 
	
	#ifdef AX
		double invDX2;
		invDX2=1/(DX*DX);
		
		#ifdef SPHERE
			double cos0;
		#endif 		
	#endif
	
	////////////////////////////////////////////////////////////////////
	// Forces in U-direction


	for(j=0; j<NY; j++){
	
		cos01=cos((lat[j]+lat[j+1])/2);
		cos1=cos(lat[j+1]);
		cos12=cos((lat[j+1]+lat[j+2])/2);
			
		for(n=0; n<NM; n++){			
			for(i=0; i<NX+1; i++){
				
				if (H[j+1][i]>H_MIN && H[j+1][i+1]>H_MIN) {
				
					// Velocity node depth
					Hu=(H[j+1][i]+H[j+1][i+1])/2;
					
					// Pressure gradient force
					#ifdef SPHERE
						Fu[n][j][i]=-Hu*(p1[n][j+1][i+1]-p1[n][j+1][i])*invDX/(A*cos1);
					#else
						Fu[n][j][i]=-Hu*(p1[n][j+1][i+1]-p1[n][j+1][i])*invDX;
					#endif
					
					// Coriolis
					#ifdef CORIOLIS						
						Fu[n][j][i]=Fu[n][j][i]+f[j+1]*(V[n][j+1][i]+V[n][j+2][i]+V[n][j+1][i+1]+V[n][j+2][i+1])/4;					
					#endif
												
					// Topographic coupling
					#ifdef MODECOUPLE

						// Note: This code assumes T is a depth integral (not depth average), hence there is no multiplication by Hu.
						for(m=0; m<NM; m++){								
							Fu[n][j][i]=Fu[n][j][i]-(T.x[n][m][j+1][i]+T.x[n][m][j+1][i+1])*(p1[m][j+1][i]+p1[m][j+1][i+1])/4;
						}
								
					#endif
					
					////////////////////////////////////////////////////
					// Frictional forces
					Fu_eps[n][j][i]=0;
					
					// Horizontal diffusion
					#ifdef AX
						#ifdef SPHERE
							dUdx2=1/(A*A*cos1*cos1)*(U[n][j+1][i+2]-2*U[n][j+1][i+1]+U[n][j+1][i]);	
														
							dUdy[0]=cos01*(U[n][j+1][i+1]-U[n][j][i+1]);
							dUdy[1]=cos12*(U[n][j+2][i+1]-U[n][j+1][i+1]);
							dUdy2=1/(A*A*cos1)*(dUdy[1]-dUdy[0]);
						
							Fu_eps[n][j][i]=Fu_eps[n][j][i]+AX*(dUdx2+dUdy2)*invDX2;
						#else
							Fu_eps[n][j][i]=Fu_eps[n][j][i]+AX*((U[n][j+1][i+2]-2*U[n][j+1][i+1]+U[n][j+1][i])
								+(U[n][j+2][i+1]-2*U[n][j+1][i+1]+U[n][j][i+1]))*invDX2;
						#endif
					#endif
					
					// Linear drag
					#ifdef R
						Fu_eps[n][j][i]=Fu_eps[n][j][i]-R*U[n][j+1][i+1];
					#endif
					
					// Quadratic bottom drag
					#ifdef CD
						absU=cabs(U[n][j+1][i+1]+I*(V[n][j+1][i]+V[n][j+2][i]+V[n][j+1][i+1]+V[n][j+2][i+1])/4);
						Fu_eps[n][j][i]=Fu_eps[n][j][i]-CD*absU*U[n][j+1][i+1]/(Hu*Hu);
					#endif									
							
					// Add the frictional forces to the other forces					
					Fu[n][j][i]=Fu[n][j][i]+Fu_eps[n][j][i];
					
				} // end-if: land mask				
				
			}
		}
	}


	////////////////////////////////////////////////////////////////////
	// Forces in V-direction
	
	
	for(j=0; j<NY+1; j++){
			      
		#if defined(AX) && defined(SPHERE)	      
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
	
					// Velocity node depth
					Hv=(H[j][i+1]+H[j+1][i+1])/2;
					
					// Pressure gradient force
					#ifdef SPHERE
						Fv[n][j][i]=-Hv*(p1[n][j+1][i+1]-p1[n][j][i+1])*invDX/A;
					#else
						Fv[n][j][i]=-Hv*(p1[n][j+1][i+1]-p1[n][j][i+1])*invDX;
					#endif
				
					// Coriolis
					#ifdef CORIOLIS
						Fv[n][j][i]=Fv[n][j][i]-f01*(U[n][j][i+1]+U[n][j][i+2]+U[n][j+1][i+1]+U[n][j+1][i+2])/4;
					#endif
															
					// Topographic coupling
					#ifdef MODECOUPLE
						
						// Note: This code assumes T is a depth integral (not depth average)
						for(m=0; m<NM; m++){
							Fv[n][j][i]=Fv[n][j][i]-(T.y[n][m][j][i+1]+T.y[n][m][j+1][i+1])*(p1[m][j][i+1]+p1[m][j+1][i+1])/4;
						}
												
					#endif
					
					
					////////////////////////////////////////////////////
					// Frictional forces
					Fv_eps[n][j][i]=0;

					// Horizontal diffusion
					#ifdef AX
						#ifdef SPHERE							
							dVdx2=1/(A*A*cos01*cos01)*(V[n][j+1][i+2]-2*V[n][j+1][i+1]+V[n][j+1][i]);
										
							dVdy[0]=cos0*(V[n][j+1][i+1]-V[n][j][i+1]);
							dVdy[1]=cos1*(V[n][j+2][i+1]-V[n][j+1][i+1]);
							dVdy2=1/(A*A*cos01)*(dVdy[1]-dVdy[0]);
												
							Fv_eps[n][j][i]=Fv_eps[n][j][i]+AX*(dVdx2+dVdy2)*invDX2;	
						#else	
							Fv_eps[n][j][i]=Fv_eps[n][j][i]+AX*((V[n][j+1][i+2]-2*V[n][j+1][i+1]+V[n][j+1][i])
								+(V[n][j+2][i+1]-2*V[n][j+1][i+1]+V[n][j][i+1]))*invDX2;
						#endif
					#endif
					
					// Linear drag
					#ifdef R
						Fv_eps[n][j][i]=Fv_eps[n][j][i]-R*V[n][j+1][i+1];
					#endif
					
					// Quadratic bottom drag
					#ifdef CD
						absV=cabs((U[n][j][i+1]+U[n][j][i+2]+U[n][j+1][i+1]+U[n][j+1][i+2])/4+I*V[n][j+1][i+1]);
						Fv_eps[n][j][i]=Fv_eps[n][j][i]-CD*absV*V[n][j+1][i+1]/(Hv*Hv);
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
			cos12=cos((lat[j+1]+lat[j+2])/2);
			
			#ifdef WORK
				gamma=invDX/(A*cos(lat[j+1]));
			#endif
		
			for(n=0; n<NMW; n++){	
				for(i=0; i<NX; i++){					
					
					#ifdef ENERGY
						KE[n][j][i]=KE[n][j][i]+(float)(0.125*RHO*((U[n][j+1][i+1]+U[n][j+1][i+2])*(U[n][j+1][i+1]+U[n][j+1][i+2])
							+(V[n][j+1][i+1]+V[n][j+2][i+1])*(V[n][j+1][i+1]+V[n][j+2][i+1]))/H[j+1][i+1]); // the extra factor of 1/4 comes from averaging U^2
						PE[n][j][i]=PE[n][j][i]+(float)(0.5*RHO*H[j+1][i+1]*(p1[n][j+1][i+1]*p1[n][j+1][i+1]/(c[n][j+1][i+1]*c[n][j+1][i+1])));
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
						
						// Compute dissipation								
						D[n][j][i]=D[n][j][i]
							-(float)(RHO*0.25*((Fu_eps[n][j][i]+Fu_eps[n][j][i+1])*(U[n][j+1][i+1]+U[n][j+1][i+2])
								+(Fv_eps[n][j][i]+Fv_eps[n][j+1][i])*(V[n][j+1][i+1]+V[n][j+2][i+1]))/H[j+1][i+1]);	
										
						// Only compute scattering if there is mode coupling.			
						#ifdef MODECOUPLE	
							
							for(m=0; m<NM; m++){				
								
								Cn[n][j][i]=Cn[n][j][i]
									+(float)(0.5*RHO*(
									 (T.x[m][n][j+1][i+1]*(U[m][j+1][i+1]+U[m][j+1][i+2])+T.y[m][n][j+1][i+1]*(V[m][j+1][i+1]+V[m][j+2][i+1]))*p1[n][j+1][i+1]
									-(T.x[n][m][j+1][i+1]*(U[n][j+1][i+1]+U[n][j+1][i+2])+T.y[n][m][j+1][i+1]*(V[n][j+1][i+1]+V[n][j+2][i+1]))*p1[m][j+1][i+1])/H[j+1][i+1]);
									 // the factor of 1/2 comes from averaging U and V
							}
							
						#endif // end modecouple if				
					#endif
				}
			}
		}		
		
	#endif // end diagnostics if
	
	////////////////////////////////////////////////////////////////////
	// Estimate SSH amplitude and phase
	#ifdef SSH
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
	
}
