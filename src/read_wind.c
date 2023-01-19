#include <math.h>
#include <netcdf.h>
#include "csw.h"

void read_wind(int sF, int rank)
{
	int i, j, ncid, varid, status;
	
	float slab1[NYW][1];
	float slab2[NYW][2];
	float slab3[NYW][NXW];

	int iW, jW;
	double dlon, dlat, w1, w2;

	size_t start_1D[1], count_1D[1];
	size_t start_3D[3], count_3D[3];

	// Open the file
	if((status = nc_open(FILE_WIND, NC_NOWRITE, &ncid)))
		ERR(status);
	
	// Load grid
	if((status = nc_inq_varid(ncid, "lon", &varid)))
		ERR(status);

	start_1D[0]=NXW-1;
	count_1D[0]=1;
	if((status = nc_get_vara_double(ncid, varid, start_1D, count_1D, &WIND.lon[0])))
		ERR(status);
	WIND.lon[0]=WIND.lon[0]-360;
	
	start_1D[0]=0;
	count_1D[0]=NXW;
	if((status = nc_get_vara_double(ncid, varid, start_1D, count_1D, &WIND.lon[1])))
		ERR(status);
		
	start_1D[0]=0;
	count_1D[0]=2;
	if((status = nc_get_vara_double(ncid, varid, start_1D, count_1D, &WIND.lon[NXW+1])))
		ERR(status);
	WIND.lon[NXW+1]=WIND.lon[NXW+1]+360;
	WIND.lon[NXW+2]=WIND.lon[NXW+2]+360;

	if((status = nc_inq_varid(ncid, "lat", &varid)))
		ERR(status);
		
	start_1D[0]=0;
	count_1D[0]=NYW;
	if((status = nc_get_vara_double(ncid, varid, start_1D, count_1D, &WIND.lat[0])))
		ERR(status);
			
	// Convert degrees to radians
	for(i=0; i<NXW+3; i++) {
		WIND.lon[i]=WIND.lon[i]*M_PI/180;
	}
	
	for(j=0; j<NYW; j++) {
		WIND.lat[j]=WIND.lat[j]*M_PI/180;
	}
					
	// Load East wind	
	if((status = nc_inq_varid(ncid, "TAUX", &varid)))
		ERR(status);
	
	start_3D[0]=sF;
	start_3D[1]=0;
	start_3D[2]=NXW-1;
	count_3D[0]=1;
	count_3D[1]=NYW;		
	count_3D[2]=1;		
	if((status = nc_get_vara_float(ncid, varid, start_3D, count_3D, &slab1[0][0])))
		ERR(status);

	start_3D[0]=sF;
	start_3D[1]=0;
	start_3D[2]=0;
	count_3D[0]=1;
	count_3D[1]=NYW;		
	count_3D[2]=2;			
	if((status = nc_get_vara_float(ncid, varid, start_3D, count_3D, &slab2[0][0])))
		ERR(status);
		
	start_3D[0]=sF;
	start_3D[1]=0;
	start_3D[2]=0;
	count_3D[0]=1;
	count_3D[1]=NYW;		
	count_3D[2]=NXW;
	if((status = nc_get_vara_float(ncid, varid, start_3D, count_3D, &slab3[0][0])))
		ERR(status);	
	
	// Glue the slabs
	for(j=0; j<NYW; j++){
		WIND.tau_x[j][0]=slab1[j][0];		
		WIND.tau_x[j][NXW+1]=slab2[j][0];		
		WIND.tau_x[j][NXW+2]=slab2[j][1];
		for(i=0; i<NXW; i++){
			WIND.tau_x[j][i+1]=slab3[j][i];		
		}		
	}
	
	// Load North wind	
	if((status = nc_inq_varid(ncid, "TAUY", &varid)))
		ERR(status);
	
	start_3D[0]=sF;
	start_3D[1]=0;
	start_3D[2]=NXW-1;
	count_3D[0]=1;
	count_3D[1]=NYW;		
	count_3D[2]=1;		
	if((status = nc_get_vara_float(ncid, varid, start_3D, count_3D, &slab1[0][0])))
		ERR(status);

	start_3D[0]=sF;
	start_3D[1]=0;
	start_3D[2]=0;
	count_3D[0]=1;
	count_3D[1]=NYW;		
	count_3D[2]=2;			
	if((status = nc_get_vara_float(ncid, varid, start_3D, count_3D, &slab2[0][0])))
		ERR(status);
		
	start_3D[0]=sF;
	start_3D[1]=0;
	start_3D[2]=0;
	count_3D[0]=1;
	count_3D[1]=NYW;		
	count_3D[2]=NXW;
	if((status = nc_get_vara_float(ncid, varid, start_3D, count_3D, &slab3[0][0])))
		ERR(status);
		
	// Glue the slabs
	for(j=0; j<NYW; j++){
		WIND.tau_y[j][0]=slab1[j][0];		
		WIND.tau_y[j][NXW+1]=slab2[j][0];		
		WIND.tau_y[j][NXW+2]=slab2[j][1];
		for(i=0; i<NXW; i++){
			WIND.tau_y[j][i+1]=slab3[j][i];		
		}		
	}
		
	// Close the file
	if ((status = nc_close(ncid)))
		ERR(status);
	
	// Zero out extreme values
	for(j=0; j<NYW; j++){
		for(i=0; i<NXW+3; i++){
			if(WIND.tau_x[j][i]<-10 || 10<WIND.tau_x[j][i]){WIND.tau_x[j][i]=0;}		
			if(WIND.tau_y[j][i]<-10 || 10<WIND.tau_y[j][i]){WIND.tau_y[j][i]=0;}
			//#ifdef NO_ANTARCTIC
			//	if(WIND.lat[j]<-60*M_PI/180){
			//		WIND.tau_x[j][i]=0;
			//		WIND.tau_y[j][i]=0;
			//	}
			//#endif						
		}		
	}
			
	// Bilinear interpolation (see Numerical Recipes in C p. 123)
	dlon=WIND.lon[1]-WIND.lon[0];
	dlat=WIND.lat[1]-WIND.lat[0];
			
	jW=0; 
	for(j=0; j<NY+2; j++){
		
		// Find the first wind grid point past the model grid point
		while(WIND.lat[jW+1]<lat[j]) {++jW;}
		
		// Latitude weight 
		w2=(lat[j]-WIND.lat[jW])/dlat;  

		iW=0;
		for(i=0; i<NX+2; i++){
			
			// Find the first wind grid point past the model grid point
			if(lon[i]<M_PI){
				// Eastern Hemisphere
				while(WIND.lon[iW+1]<lon[i]) {++iW;}

				w1=(lon[i]-WIND.lon[iW])/dlon; // Longitude weight 
			}
			else{
				// Western Hemisphere (add 360 degrees to WIND.lon)
				if(WIND.lon[iW+1]>0) {iW=0;} // if the wind index is in the Eastern Hemisphere, reset the index to zero (this is only true right after we cross the Grand Meridian)				
				while((WIND.lon[iW+1]+2*M_PI)<lon[i]) {++iW;}
				
				w1=(lon[i]-(WIND.lon[iW]+2*M_PI))/dlon; // Longitude weight 
			}
								
			tau_x[j][i]=(double)((1-w1)*(1-w2)*WIND.tau_x[jW][iW]
				+w1*(1-w2)*WIND.tau_x[jW][iW+1]
				+w1*w2*WIND.tau_x[jW+1][iW+1]
				+(1-w1)*w2*WIND.tau_x[jW+1][iW]);		
				
			tau_y[j][i]=(double)((1-w1)*(1-w2)*WIND.tau_y[jW][iW]
				+w1*(1-w2)*WIND.tau_y[jW][iW+1]
				+w1*w2*WIND.tau_y[jW+1][iW+1]
				+(1-w1)*w2*WIND.tau_y[jW+1][iW]);
						
		}		
	}
	
	// Print progress
	if (rank == 0) {
		printf ("Wind %d \n",sF);
	}
	
}
