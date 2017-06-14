#include <math.h>
#include <csw.h>


void calc_ITGF(double t)
{
	int i, j, k, n;
	double complex phaseFp;
	
	#ifdef WORK
		double complex phaseFt;
	#endif 
	
	for(k=0; k<NC; k++){ // Cycle through frequencies 

		phaseFp=cexp(I*ITGF.omega[k]*(t+DT/2)); // Compute forcing at time=t+dt/2 (i.e., the time of p)
		
		#ifdef WORK			
			phaseFt=cexp(I*ITGF.omega[k]*t); // Compute generation at time=t (i.e., the time of p1)
		#endif
	
		for(n=0; n<NM; n++){
			for(j=0; j<NY; j++){
				for(i=0; i<NX; i++){					
					
					if (H[j+1][i+1]>H_MIN_FORCE) {					
					
						// Add ITGF to pressure forcing at the current time=t-dt (we're trying to update the variables to time=t) 
						Fp[n][j][i]=Fp[n][j][i]-c[n][j+1][i+1]*c[n][j+1][i+1]*phi_bott[n][j+1][i+1]*creal(ITGF.F[k][j][i]*phaseFp);
					
						#ifdef WORK			
							// Compute internal-tide generation between the last two time steps
							C0[n][j][i]=C0[n][j][i]-(float)(RHO*H[j+1][i+1]*creal(ITGF.F[k][j][i]*phaseFt)*p1[n][j+1][i+1]*phi_bott[n][j+1][i+1]);							
						#endif					
					
					}
												
				}
			}
		}
	}
		
}
