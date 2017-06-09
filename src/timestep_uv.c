#include <csw.h>

void timestep_uv(void)
{ 
	int i, j, n;
	double gamma;
	
	gamma=DT/12;
		
	for(n=0; n<NM; n++){
		
		// Update U
		for(j=0; j<NY; j++){
			for(i=0; i<NX+1; i++){
				
				if (mask.U[n][j][i]>1E-10) { // not sure if this is faster
	
					// Save the current data
					U1[n][j+1][i+1]=U[n][j+1][i+1];
	
					// Integrate to find the new data
					U[n][j+1][i+1]=(U1[n][j+1][i+1]+gamma*(23*Fu[n][j][i]-16*Fu_1[n][j][i]+5*Fu_2[n][j][i]))*mask.U[n][j][i];
					
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
				
				if (mask.V[n][j][i]>1E-10) {

					// Save the current data
					V1[n][j+1][i+1]=V[n][j+1][i+1];
					
					// Integrate to find the new data
					V[n][j+1][i+1]=(V1[n][j+1][i+1]+gamma*(23*Fv[n][j][i]-16*Fv_1[n][j][i]+5*Fv_2[n][j][i]))*mask.V[n][j][i];
									
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
