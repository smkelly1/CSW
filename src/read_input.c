#include <math.h>
#include <netcdf.h>
#include <csw.h>

void read_input(int rank)
{
	
	int x0, y0, ncid, varid, status;
	
	// Define start points
	y0=(int)(floor(rank/NPX)*NY); // y start
	x0=(int)(floor(rank-floor(rank/NPX)*NPX)*NX); // x start
	
	size_t start_mask[]={0, y0, x0}; 

	// Define counts to read 
	size_t count_mask[]={NM, NY, NX}; 
	
	// Define variables that might be needed
	#ifdef IT_FORCING	
		size_t start_ITGF[]={0, 0, y0, x0};
		size_t count_ITGF[]={NC, NM, NY, NX};
	#endif


	// Open the file
	if ((status = nc_open(FILE_IN, NC_NOWRITE, &ncid)))
		ERR(status);
		
		
	// Masks	
	if ((status = nc_inq_varid(ncid, "mask_p", &varid)))
		ERR(status);

	if ((status = nc_get_vara_double(ncid, varid, start_mask, count_mask, &mask.p[0][0][0])))
		ERR(status);

	if ((status = nc_inq_varid(ncid, "mask_u", &varid)))
		ERR(status);

	count_mask[2]=NX+1;
	if ((status = nc_get_vara_double(ncid, varid, start_mask, count_mask, &mask.U[0][0][0])))
		ERR(status);
	count_mask[2]=NX;
		
	if ((status = nc_inq_varid(ncid, "mask_v", &varid)))
		ERR(status);

	count_mask[1]=NY+1;
	if ((status = nc_get_vara_double(ncid, varid, start_mask, count_mask, &mask.V[0][0][0])))
		ERR(status);
	count_mask[1]=NY;

			
	// Forcing	
	#ifdef IT_FORCING	
			
		if ((status = nc_inq_varid(ncid, "ITGF_omega", &varid)))
			ERR(status);

		if ((status = nc_get_var_double(ncid, varid, &ITGF.omega[0])))
			ERR(status);		
		
		if ((status = nc_inq_varid(ncid, "ITGF_pr", &varid)))
			ERR(status);

		if ((status = nc_get_vara_double(ncid, varid, start_ITGF, count_ITGF, &ITGF.pr[0][0][0][0])))
			ERR(status);	
			
		if ((status = nc_inq_varid(ncid, "ITGF_pi", &varid)))
			ERR(status);

		if ((status = nc_get_vara_double(ncid, varid, start_ITGF, count_ITGF, &ITGF.pi[0][0][0][0])))
			ERR(status);
			
	#endif
		

}
