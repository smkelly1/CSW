#include <netcdf.h>
#include "csw.h"


void init_output(int rank)
{
	int  ncid, varid, status, dimid[4];
	char name[100];

	// Get name of output file
	sprintf(name,FILE_OUT ".%03d.nc",rank);    

	// Create file 
	if ((status = nc_create(name, NC_CLOBBER, &ncid)))
		ERR(status);

	// Create dimensions  
	if ((status = nc_def_dim(ncid, "x", NX, &dimid[3])))
		ERR(status);

	if ((status = nc_def_dim(ncid, "y", NY, &dimid[2])))
		ERR(status);

	if ((status = nc_def_dim(ncid, "mode", NMW, &dimid[1])))
		ERR(status);

	if ((status = nc_def_dim(ncid, "time", NC_UNLIMITED, &dimid[0])))
		ERR(status);

	// Define the variables
	if ((status = nc_def_var(ncid, "yday", NC_DOUBLE, 1, &dimid[0], &varid)))
		ERR(status);

	// Write pressure
	#ifdef WRITE_PRESSURE
		if ((status = nc_def_var(ncid, "p", NC_FLOAT, 4, dimid, &varid)))
			ERR(status);
	#endif

	// Write velocity
	#ifdef WRITE_VELOCITY
		if ((status = nc_def_var(ncid, "u", NC_FLOAT, 4, dimid, &varid)))
		ERR(status);
	
		if ((status = nc_def_var(ncid, "v", NC_FLOAT, 4, dimid, &varid)))
		ERR(status);
	#endif

	// Define the diagnostics 
	#ifdef WRITE_SSH
		if ((status = nc_def_var(ncid, "SSH_amp", NC_FLOAT, 4, dimid, &varid)))
		ERR(status);

		if ((status = nc_def_var(ncid, "SSH_phase", NC_FLOAT, 4, dimid, &varid)))
		ERR(status);
	#endif

	#ifdef WRITE_TRANSPORT
		if ((status = nc_def_var(ncid, "U_amp", NC_FLOAT, 4, dimid, &varid)))
		ERR(status);

		if ((status = nc_def_var(ncid, "U_phase", NC_FLOAT, 4, dimid, &varid)))
		ERR(status);
		
		if ((status = nc_def_var(ncid, "V_amp", NC_FLOAT, 4, dimid, &varid)))
		ERR(status);

		if ((status = nc_def_var(ncid, "V_phase", NC_FLOAT, 4, dimid, &varid)))
		ERR(status);
	#endif

	#ifdef ENERGY
		if ((status = nc_def_var(ncid, "KE", NC_FLOAT, 4, dimid, &varid)))
			ERR(status);
    
		if ((status = nc_def_var(ncid, "PE", NC_FLOAT, 4, dimid, &varid)))
			ERR(status);
	#endif
	
	#ifdef FLUX
		if ((status = nc_def_var(ncid, "up", NC_FLOAT, 4, dimid, &varid)))
			ERR(status);

		if ((status = nc_def_var(ncid, "vp", NC_FLOAT, 4, dimid, &varid)))
			ERR(status);
	#endif

	#ifdef WORK  
		if ((status = nc_def_var(ncid, "C0", NC_FLOAT, 4, dimid, &varid)))
			ERR(status);

		if ((status = nc_def_var(ncid, "Cn", NC_FLOAT, 4, dimid, &varid)))
			ERR(status);

		if ((status = nc_def_var(ncid, "D", NC_FLOAT, 4, dimid, &varid)))
			ERR(status);

		if ((status = nc_def_var(ncid, "divF", NC_FLOAT, 4, dimid, &varid)))
			ERR(status);

		if ((status = nc_def_var(ncid, "error", NC_FLOAT, 4, dimid, &varid)))
			ERR(status);
	#endif

	#if defined(ENERGY) || defined(FLUX) || defined(WORK) || defined(WRITE_SSH) || defined(WRITE_TRANSPORT)
		if ((status = nc_def_var(ncid, "period", NC_INT, 1, &dimid[0], &varid)))
			ERR(status);
	#endif 


	// Close the file
	if ((status = nc_close(ncid)))
		ERR(status);
}
