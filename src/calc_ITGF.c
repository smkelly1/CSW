#include <math.h>
#include <csw.h>


void calc_ITGF(double t)
{
	int i, j, k, n;
	double complex phaseFp;
	
	#ifdef DIAGNOSTICS
		double complex phaseFt;
		double Ft;
	#endif 
	
	for(k=0; k<NC; k++){ // Cycle through frequencies 

		phaseFp=cexp(I*ITGF.omega[k]*(t+DT/2)); // Compute forcing at time=t+dt/2 (i.e., the time of p)
		
		#ifdef DIAGNOSTICS			
			phaseFt=cexp(I*ITGF.omega[k]*t); // Compute generation at time=t (i.e., the time of p1)
		#endif
	
		for(n=0; n<NM; n++){
			for(j=0; j<NY; j++){
				for(i=0; i<NX; i++){					
					
					#ifdef NO_FORCING_MASK
						if (mask.p[n][j][i]>0.999) {
					#endif
					
						// Add ITGF to pressure forcing at the current time=t-dt (we're trying to update the variables to time=t) 
						Fp[n][j][i]=Fp[n][j][i]-c[n][j+1][i+1]*c[n][j+1][i+1]*creal((ITGF.pr[k][n][j][i]+I*ITGF.pi[k][n][j][i])*phaseFp);
					
						#ifdef DIAGNOSTICS			
							// Compute internal-tide generation between the last two time steps
							Ft=-creal((ITGF.pr[k][n][j][i]+I*ITGF.pi[k][n][j][i])*phaseFt);		
							C0[n][j][i]=C0[n][j][i]+(float)(RHO*H[j+1][i+1]*Ft*p1[n][j+1][i+1]);							
						#endif
					
					# ifdef NO_FORCING_MASK	
					}
					#endif
												
				}
			}
		}
	}
		
}
