clear; 
ii=complex(0,1);

% Specify file names
file.out='m2UV_FES.nc';

% start output file
mode = netcdf.getConstant('CLOBBER'); % Erase existing file
mode = bitor(mode,netcdf.getConstant('64BIT_OFFSET')); % Allow giant files
temp=netcdf.create(file.out,mode); % Create file
netcdf.close(temp);
ncwriteatt(file.out,'/','Source','FES2014') % Add metadata 

% Read U  and u-grid
file.in='fes2014a_currents/eastward_velocity/m2.nc';
UAMP=ncread(file.in,'Ua')/100;
UPHA=ncread(file.in,'Ug');
UPHA(UPHA<0)=UPHA(UPHA<0)+360;
[Nx Ny]=size(UAMP);

lon=ncread(file.in,'lon');
lat=ncread(file.in,'lat');

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

clear U UAMP UPHA

% Read V  
file.in='fes2014a_currents/northward_velocity/m2.nc';
VAMP=ncread(file.in,'Va')/100;
VPHA=ncread(file.in,'Vg');
VPHA(VPHA<0)=VPHA(VPHA<0)+360;
[Nx Ny]=size(VAMP);

nccreate(file.out,'VAMP','Dimensions',{'lon',Nx,'lat',Ny});
ncwriteatt(file.out,'VAMP','long_name','Northward velocity amplitude m/s');
ncwrite(file.out,'VAMP',VAMP);    
  
nccreate(file.out,'VPHA','Dimensions',{'lon',Nx,'lat',Ny});
ncwriteatt(file.out,'VPHA','long_name','Nothward transport Greenwich phase [deg]');
ncwrite(file.out,'VPHA',VPHA);  

