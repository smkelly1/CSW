% Save delta_c2 and c0 for Coupled Shallow Water model (CSW)
clear

ii=complex(0,1);
warning off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set input parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
res=25; 
Nm=1;
strat='WOA';
bathy='SS';

% Input files
fid.grid=['../../../18-6_grids/',num2str(res),'th_deg_',strat,'_',bathy,'_grid.nc'];
fid.HYCOM=['../../../HYCOM/2015_c.nc'];

% Output files
fid.delta_c2='../../delta_c2.nc';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load CSW grid
H=ncread(fid.grid,'H');
c=ncread(fid.grid,'c');
c=c(:,:,1:Nm);

lon=ncread(fid.grid,'lon');
lat=ncread(fid.grid,'lat');
[Nx Ny]=size(H);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write delta_c2 for model grid

% Load HYCOM
HYCOM.lon=ncread(fid.HYCOM,'lon');
HYCOM.lat=ncread(fid.HYCOM,'lat');
HYCOM.delta_c2=ncread(fid.HYCOM,'delta_c_var');

% Compute variance and interpolate
delta_c2=interp2(HYCOM.lon,HYCOM.lat',HYCOM.delta_c2',lon,lat')';

% Initiate file
mode = netcdf.getConstant('CLOBBER');
mode = bitor(mode,netcdf.getConstant('NETCDF4'));
temp=netcdf.create(fid.delta_c2,mode);
netcdf.close(temp);

% Write fields
nccreate(fid.delta_c2,'lon','Dimensions',{'x',Nx});
nccreate(fid.delta_c2,'lat','Dimensions',{'y',Ny});
nccreate(fid.delta_c2,'H','Dimensions',{'x',Nx,'y',Ny});
nccreate(fid.delta_c2,'c','Dimensions',{'x',Nx,'y',Ny,'mode',Nm});
nccreate(fid.delta_c2,'delta_c2','Dimensions',{'x',Nx,'y',Ny,'mode',Nm});

ncwrite(fid.delta_c2,'lon',lon);
ncwrite(fid.delta_c2,'lat',lat);
ncwrite(fid.delta_c2,'H',H);
ncwrite(fid.delta_c2,'c',c);
ncwrite(fid.delta_c2,'delta_c2',delta_c2);

