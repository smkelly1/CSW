#include <math.h>
#include <netcdf.h>
#include "csw.h"

void read_grid(int rank)
{

	int x0, y0, ncid, varid, status;

	// Define start points 
	y0=(int)(floor(rank/NPX)*NY); // y start
	x0=(int)(floor(rank-floor(rank/NPX)*NPX)*NX); // x start

	size_t start_H[]={y0, x0};
	size_t start_c[]={0, y0, x0};

	// Define counts to read 
	size_t count_H[]={NY+2, NX+2};
	size_t count_c[]={NM, NY+2, NX+2};
	
	// Define variables that might be needed
	#ifdef SSH
		size_t start_phi0[]={0, y0, x0};
		size_t count_phi0[]={NMW, NY+2, NX+2};
	#endif

	#ifdef MODECOUPLE
		size_t start_T[]={0, 0, y0, x0};
		size_t count_T[]={NM, NM, NY+2, NX+2};
	#endif
		
	#if defined(SPHERE) || defined(CORIOLIS)
		int j;
		size_t start_lat[]={y0};
		size_t count_lat[]={NY+2};
	#endif
		
		
	// Open the file
	if ((status = nc_open(FILE_GRID, NC_NOWRITE, &ncid)))
		ERR(status);
	
				
	// Read the grid	
	if ((status = nc_inq_varid(ncid, "H", &varid)))
		ERR(status);

	if ((status = nc_get_vara_double(ncid, varid, start_H, count_H, &H[0][0])))
		ERR(status);
		
	if ((status = nc_inq_varid(ncid, "c", &varid)))
		ERR(status);

	if ((status = nc_get_vara_double(ncid, varid, start_c, count_c, &c[0][0][0])))
		ERR(status);
	
	#ifdef IT_FORCING
		if ((status = nc_inq_varid(ncid, "phi_bott", &varid)))
			ERR(status);

		if ((status = nc_get_vara_double(ncid, varid, start_c, count_c, &phi_bott[0][0][0])))
			ERR(status);
	#endif
	
	#ifdef SSH
		if ((status = nc_inq_varid(ncid, "phi_surf", &varid)))
			ERR(status);

		if ((status = nc_get_vara_double(ncid, varid, start_phi0, count_phi0, &phi_surf[0][0][0])))
			ERR(status);
	#endif
	
	// Read the topographic coupling coefficients if necessary
	#ifdef MODECOUPLE
		
		if ((status = nc_inq_varid(ncid, "T_x", &varid)))
			ERR(status);

		if ((status = nc_get_vara_double(ncid, varid, start_T, count_T, &T.x[0][0][0][0])))
			ERR(status);
	
		if ((status = nc_inq_varid(ncid, "T_y", &varid)))
			ERR(status);

		if ((status = nc_get_vara_double(ncid, varid, start_T, count_T, &T.y[0][0][0][0])))
			ERR(status);
			
	#endif
	
	// Read latitude in degrees and convert to radians
	#if defined(SPHERE) || defined(CORIOLIS)
	
		if ((status = nc_inq_varid(ncid, "lat", &varid)))
			ERR(status);

		if ((status = nc_get_vara_double(ncid, varid, start_lat, count_lat, &lat[0])))
			ERR(status);
			
		for(j=0; j<NY+2; j++){
			lat[j]=lat[j]*M_PI/180;
			f[j]=2*(2*M_PI)/(24*3600)*sin(lat[j]);	
		}	
	
	#endif
	
}	
