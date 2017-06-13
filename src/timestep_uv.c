#include <math.h>
#include <csw.h>

void timestep_uv(void)
{ 
	int i, j, n;
	double gamma;
	
	#ifdef R_MASK
		double c_crit, mask;
	#endif
	
	gamma=DT/12;
		
	for(n=0; n<NM; n++){
		
		// Update U
		for(j=0; j<NY; j++){
			for(i=0; i<NX+1; i++){
				
				// Slowest resolvable wave on the grid (Adcroft et al. 1999) 
				#ifdef R_MASK
					c_crit=fabs(f[j+1])*A*cos(lat[j+1])*DX;  
				#endif
				
				if (H[j+1][i]>0 && H[j+1][i+1]>0) {
	
					// Create a sponge with an R_MASK decay time scale where there is insufficient wave resolution
					#ifdef R_MASK
						mask=1-fmax(0,(1-((c[n][j+1][i]+c[n][j+1][i+1])/2)/c_crit)*DT*R_MASK);  
					#endif
	
					// Save the current data
					U1[n][j+1][i+1]=U[n][j+1][i+1];
	
					// Integrate to find the new data
					#ifdef R_MASK
						U[n][j+1][i+1]=(U1[n][j+1][i+1]+gamma*(23*Fu[n][j][i]-16*Fu_1[n][j][i]+5*Fu_2[n][j][i]))*mask;
					#else
						U[n][j+1][i+1]=U1[n][j+1][i+1]+gamma*(23*Fu[n][j][i]-16*Fu_1[n][j][i]+5*Fu_2[n][j][i]);
					#endif 
					
					// Old forcing becomes very old forcing 
					Fu_2[n][j][i]=Fu_1[n][j][i];
	
					// Current forcing becomes old forcing 
					Fu_1[n][j][i]=Fu[n][j][i];
	
					// Average to find the mid-point data
					U1[n][j+1][i+1]=(U[n][j+1][i+1]+U1[n][j+1][i+1])/2;
				
				}
			}
		}
		
		// Update V
		for(j=0; j<NY+1; j++){
			for(i=0; i<NX; i++){
				
				// Slowest resolvable wave on the grid (Adcroft et al. 1999) 
				#ifdef R_MASK
					c_crit=fabs((f[j]+f[j+1])/2)*A*cos(lat[j+1])*DX;  
				#endif
				
				if (H[j][i+1]>0 && H[j+1][i+1]>0) {

					// Create a sponge with an R_MASK decay time scale where there is insufficient wave resolution
					#ifdef R_MASK
						mask=1-fmax(0,(1-((c[n][j][i+1]+c[n][j+1][i+1])/2)/c_crit)*DT*R_MASK);  
					#endif

					// Save the current data
					V1[n][j+1][i+1]=V[n][j+1][i+1];
					
					// Integrate to find the new data
					#ifdef R_MASK
						V[n][j+1][i+1]=(V1[n][j+1][i+1]+gamma*(23*Fv[n][j][i]-16*Fv_1[n][j][i]+5*Fv_2[n][j][i]))*mask;
					#else						
						V[n][j+1][i+1]=V1[n][j+1][i+1]+gamma*(23*Fv[n][j][i]-16*Fv_1[n][j][i]+5*Fv_2[n][j][i]);				
					#endif			
					
					// Old forcing becomes very old forcing 
					Fv_2[n][j][i]=Fv_1[n][j][i];
					
					// Current forcing becomes old forcing 
					Fv_1[n][j][i]=Fv[n][j][i];
					
					// Average to find the mid-point data
					V1[n][j+1][i+1]=(V[n][j+1][i+1]+V1[n][j+1][i+1])/2;
				
				}
			}
		}
		
	} // end n-loop
		
}
