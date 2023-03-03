#!/bin/bash

# Start date
#model_date=19871001 # Ocean Storms, 100 era 31 days
#model_date=20190526 # NISKINE, 400 era, 31 days
model_date=20190420  # HYCOM, 400 era, 60 days
#model_date=19861201 # Winter 1987 + 1 month spin-up, 100 era, 121 days
#model_date=19870601 # Summer 1987 + 1 month spin-up, 100 era, 122 days
era=400
for day in {1..60}; do # Days to download 
	
	# Download files
	wget --user=XXXX --password=XXXX --content-disposition "https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/MERRA2/M2T1NXFLX.5.12.4/"${model_date:0:4}"/"${model_date:4:2}"/MERRA2_"$era".tavg1_2d_flx_Nx."$model_date".nc4.nc?TAUX[0:1:23][0:1:360][0:1:575],TAUY[0:1:23][0:1:360][0:1:575],lat[0:1:360],lon[0:1:575],time[0:1:23]" --waitretry=1

	# make time the record dimension then remove original file
	ncks -O --mk_rec_dmn time MERRA2_"$era".tavg1_2d_flx_Nx."$model_date".nc4.nc "$model_date".nc
	rm MERRA2_"$era".tavg1_2d_flx_Nx."$model_date".nc4.nc

	# Update the date
	model_date=$(date -d "$model_date +1 days" +%Y%m%d)
done

# combine all the files 
ncrcat -O "${model_date:0:2}"*.nc out.nc
