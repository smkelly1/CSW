#include <netcdf.h>
#include <math.h>
#include <csw.h>

void write_diagnostics(int sD, int Na, int rank)
{
	int i, j, n; 
	int ncid, status, varid;
	char name[100];

	size_t start[]={sD-1, 0, 0, 0};
	size_t count[]={1, NM, NY, NX};

	////////////////////////////////////////////////////////////////////
	// Average diagnostic fields  
	for(n=0; n<NM; n++){
		for(j=0; j<NY; j++){
			for(i=0; i<NX; i++){
				
				#ifdef ENERGY
					KE[n][j][i]=KE[n][j][i]/Na;
					PE[n][j][i]=PE[n][j][i]/Na;
				#endif
							
				#ifdef FLUX
					up[n][j][i]=up[n][j][i]/Na;						
					vp[n][j][i]=vp[n][j][i]/Na;
				#endif
							
				#ifdef WORK
					C0[n][j][i]=C0[n][j][i]/Na;					
					Cn[n][j][i]=Cn[n][j][i]/Na;
					D[n][j][i]=D[n][j][i]/Na;
					divF[n][j][i]=divF[n][j][i]/Na;					
	
					tmp[n][j][i]=C0[n][j][i]+Cn[n][j][i]-divF[n][j][i]-D[n][j][i]; // residual/error
				#endif
			}
		}
	}


	////////////////////////////////////////////////////////////////////
	// Get file name
  	sprintf(name,FILE_OUT ".%03d.nc",rank);    

	// Open the file
	if ((status = nc_open(name, NC_WRITE, &ncid)))
		ERR(status);

	// Write fields	
	#ifdef ENERGY
		if ((status = nc_inq_varid(ncid, "KE", &varid)))
			ERR(status);
	
		if ((status = nc_put_vara_float(ncid, varid, start, count, &KE[0][0][0])))
			ERR(status);
	    
	   	if ((status = nc_inq_varid(ncid, "PE", &varid)))
			ERR(status);
	
		if ((status = nc_put_vara_float(ncid, varid, start, count, &PE[0][0][0])))
			ERR(status);
	#endif
	
	#ifdef FLUX
	   	if ((status = nc_inq_varid(ncid, "up", &varid)))
			ERR(status);
	
		if ((status = nc_put_vara_float(ncid, varid, start, count, &up[0][0][0])))
			ERR(status);
			
		if ((status = nc_inq_varid(ncid, "vp", &varid)))
			ERR(status);
	
		if ((status = nc_put_vara_float(ncid, varid, start, count, &vp[0][0][0])))
			ERR(status);
	#endif
	
	#ifdef WORK
		if ((status = nc_inq_varid(ncid, "C0", &varid)))
			ERR(status);
	
		if ((status = nc_put_vara_float(ncid, varid, start, count, &C0[0][0][0])))
			ERR(status);
			
	   	if ((status = nc_inq_varid(ncid, "Cn", &varid)))
			ERR(status);
	
		if ((status = nc_put_vara_float(ncid, varid, start, count, &Cn[0][0][0])))
			ERR(status);
			
		if ((status = nc_inq_varid(ncid, "D", &varid)))
			ERR(status);
	
		if ((status = nc_put_vara_float(ncid, varid, start, count, &D[0][0][0])))
			ERR(status);
	
		if ((status = nc_inq_varid(ncid, "divF", &varid)))
			ERR(status);
	
		if ((status = nc_put_vara_float(ncid, varid, start, count, &divF[0][0][0])))
			ERR(status);
	    
	    if ((status = nc_inq_varid(ncid, "error", &varid)))
			ERR(status);
	
		if ((status = nc_put_vara_float(ncid, varid, start, count, &tmp[0][0][0])))
			ERR(status);
    #endif    
        
    // Write time  
    if ((status = nc_inq_varid(ncid, "period", &varid)))
		ERR(status);

	if ((status = nc_put_vara_int(ncid, varid, &start[0], &count[0], &sD)))
		ERR(status);
        
        
	// Close the file
	if ((status = nc_close(ncid)))
		ERR(status);
		
		
	////////////////////////////////////////////////////////////////////	
	// Clear diagnostic fields
	for(n=0; n<NM; n++){
		for(j=0; j<NY; j++){
			for(i=0; i<NX; i++){
				
				#ifdef ENERGY
					KE[n][j][i]=0;
					PE[n][j][i]=0;
				#endif
				
				#ifdef FLUX
					up[n][j][i]=0;
					vp[n][j][i]=0;
				#endif
				
				#ifdef WORK
					C0[n][j][i]=0;
					Cn[n][j][i]=0;
					divF[n][j][i]=0;
					D[n][j][i]=0;
				#endif
				
			}
		}
	}
	
	
}
