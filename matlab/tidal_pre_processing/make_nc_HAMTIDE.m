clear; 
ii=complex(0,1);

% Specify file names
file.out='m2UV_HAMTIDE.nc';

% start output file
mode = netcdf.getConstant('CLOBBER'); % Erase existing file
mode = bitor(mode,netcdf.getConstant('64BIT_OFFSET')); % Allow giant files
temp=netcdf.create(file.out,mode); % Create file
netcdf.close(temp);
ncwriteatt(file.out,'/','Source','HAMTIDE11a') % Add metadata 

% Read data
file.in='HAMcurrent11a_m2.nc';
lon=ncread(file.in,'LON'); lon=lon(1:end-1);
lat=ncread(file.in,'LAT');

Nx=length(lon);
Ny=length(lat);

UAMP=ncread(file.in,'UAMP')/100; UAMP=UAMP(1:end-1,:);
UPHA=ncread(file.in,'UPHA'); UPHA=UPHA(1:end-1,:);

VAMP=ncread(file.in,'VAMP')/100; VAMP=VAMP(1:end-1,:);
VPHA=ncread(file.in,'VPHA'); VPHA=VPHA(1:end-1,:);

% Write data
nccreate(file.out,'lon','Dimensions',{'lon',Nx});
ncwriteatt(file.out,'lon', '_CoordinateAxisType','lon');
ncwrite(file.out,'lon',lon);

nccreate(file.out,'lat','Dimensions',{'lat',Ny});
ncwriteatt(file.out,'lat', '_CoordinateAxisType','lat');
ncwrite(file.out,'lat',lat);

nccreate(file.out,'UAMP','Dimensions',{'lon',Nx,'lat',Ny});
ncwriteatt(file.out,'UAMP','long_name','Eastward velocity amplitude m/s');
ncwrite(file.out,'UAMP',UAMP);    
  
nccreate(file.out,'UPHA','Dimensions',{'lon',Nx,'lat',Ny});
ncwriteatt(file.out,'UPHA','long_name','Eastward transport Greenwich phase [deg]');
ncwrite(file.out,'UPHA',UPHA);

nccreate(file.out,'VAMP','Dimensions',{'lon',Nx,'lat',Ny});
ncwriteatt(file.out,'VAMP','long_name','Northward velocity amplitude m/s');
ncwrite(file.out,'VAMP',VAMP);    
  
nccreate(file.out,'VPHA','Dimensions',{'lon',Nx,'lat',Ny});
ncwriteatt(file.out,'VPHA','long_name','Northward transport Greenwich phase [deg]');
ncwrite(file.out,'VPHA',VPHA);
