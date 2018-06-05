clear; 
ii=complex(0,1);

% Specify file names
file.in='./DATA/uv.m2_tpxo8_atlas_30c.nc';
file.out='m2UV_TPXO.nc';

% start output file
mode = netcdf.getConstant('CLOBBER'); % Erase existing file
mode = bitor(mode,netcdf.getConstant('64BIT_OFFSET')); % Allow giant files
temp=netcdf.create(file.out,mode); % Create file
netcdf.close(temp);
ncwriteatt(file.out,'/','Source','TPXO 8 atlas') % Add metadata 

% Nodal corrections f and u
day=0;
omegaN=125.0445-0.05295377*(day-51544.4993);% mean longitude of ascending lunar node
p=83.3535+0.11140353*(day-51544.4993);% mean longitude of lunar perigee
sinn = sin(omegaN*pi/180);
cosn = cos(omegaN*pi/180);
sin2n = sin(2*omegaN*pi/180);
cos2n = cos(2*omegaN*pi/180);
sin3n = sin(3*omegaN*pi/180);

omega=1.405189e-04;
phase=1.731557546;
u=atan((-.03731*sinn+.00052*sin2n)./(1.-.03731*cosn+.00052*cos2n))/pi/180;                                   
f=sqrt((1.-.03731*cosn+.00052*cos2n).^2+(.03731*sinn-.00052*sin2n).^2); 

% Read and write grid
lon=ncread(file.in,'lon_v');
lat=ncread(file.in,'lat_u');

Ny=length(lat);
Nx=length(lon);

nccreate(file.out,'lon','Dimensions',{'lon',Nx});
ncwriteatt(file.out,'lon', '_CoordinateAxisType','lon');
ncwrite(file.out,'lon',lon);

nccreate(file.out,'lat','Dimensions',{'lat',Ny});
ncwriteatt(file.out,'lat', '_CoordinateAxisType','lat');
ncwrite(file.out,'lat',lat);

% Read and write U
U=double(ncread(file.in,'uRe'))'+ii*double(ncread(file.in,'uIm'))';
U(end+1,:)=U(1,:);
U=(U(1:end-1,:)+U(2:end,:))/2;
U=U.'*f*exp(ii*u/180*pi)/(100^2); % Transpose, phase correction and convert to m^2/s
U(U==0)=NaN;
UAMP=abs(U);
UPHA=-angle(U)/pi*180;UPHA(UPHA<0)=UPHA(UPHA<0)+360;

nccreate(file.out,'UAMP','Dimensions',{'lon',Nx,'lat',Ny});
ncwriteatt(file.out,'UAMP','long_name','Eastward transport amplitude m^2/s');
ncwrite(file.out,'UAMP',UAMP.');    
  
nccreate(file.out,'UPHA','Dimensions',{'lon',Nx,'lat',Ny});
ncwriteatt(file.out,'UPHA','long_name','Eastward transport Greenwich phase [deg]');
ncwrite(file.out,'UPHA',UPHA.');  

clear U UAMP UPHA

% Read and write V
V=double(ncread(file.in,'vRe'))'+ii*double(ncread(file.in,'vIm'))';
V(:,end+1)=V(:,end);
V=(V(:,1:end-1)+V(:,2:end))/2;
V=V.'*f*exp(ii*u/180*pi)/(100^2); % Transpose, phase correction and convert to m^2/s
V(V==0)=NaN;
VAMP=abs(V);
VPHA=-angle(V)/pi*180;VPHA(VPHA<0)=VPHA(VPHA<0)+360;

nccreate(file.out,'VAMP','Dimensions',{'lon',Nx,'lat',Ny});
ncwriteatt(file.out,'VAMP','long_name','Northward transport amplitude m^2/s');
ncwrite(file.out,'VAMP',VAMP.');    
  
nccreate(file.out,'VPHA','Dimensions',{'lon',Nx,'lat',Ny});
ncwriteatt(file.out,'VPHA','long_name','Nothward transport Greenwich phase [deg]');
ncwrite(file.out,'VPHA',VPHA.');  

