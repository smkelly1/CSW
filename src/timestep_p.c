#include <csw.h>

void timestep_p(void)
{ 
	int i, j, n;
	double gamma;
	
	gamma=DT/12;

	for(n=0; n<NM; n++){
		for(j=0; j<NY; j++){
			for(i=0; i<NX; i++){
				
				//if (mask.p[n][j][i]>1E-10) {

					// Save the current data
					p1[n][j+1][i+1]=p[n][j+1][i+1];
					
					// Integrate to find the new data
					p[n][j+1][i+1]=(p1[n][j+1][i+1]+gamma*(23*Fp[n][j][i]-16*Fp_1[n][j][i]+5*Fp_2[n][j][i]))*mask.p[n][j][i];
									
					// Old forcing becomes very old forcing 
					Fp_2[n][j][i]=Fp_1[n][j][i];
					
					// Cuttent forcing becomes old forcing 
					Fp_1[n][j][i]=Fp[n][j][i];
					
					// Average to find the mid-point data
					p1[n][j+1][i+1]=(p[n][j+1][i+1]+p1[n][j+1][i+1])/2;
				
				//}
				
			}
		}
	} 
}
