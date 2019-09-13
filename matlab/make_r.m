% Compute drag parameterizations for Coupled Shallow Water model (CSW)
clear

ii=complex(0,1);
warning off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set input parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
res=25;
strat='WOA';
bathy='SS';
Nm=8;
omega=2*pi/(12.42*3600);
r_min=0;
r_max=1/(24*3600); % Min and max values of r

fid.grid=['../../../18-6_grids/',num2str(res),'th_deg_',strat,'_',bathy,'_grid.nc'];
fid.EKE=['../../../AVISO/AVISO.nc'];

folder.r='../../';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load grid
H=ncread(fid.grid,'H');
c=ncread(fid.grid,'c');
f=abs(ncread(fid.grid,'f'));
lon=ncread(fid.grid,'lon');
lat=ncread(fid.grid,'lat');
[Nx Ny]=size(H);

% Load AVISO
AVISO.EKE=ncread(fid.EKE,'EKE');
AVISO.H2=ncread(fid.EKE,'H2');
AVISO.H=ncread(fid.EKE,'H');
AVISO.N_bott=ncread(fid.EKE,'N_bott');
AVISO.lon=ncread(fid.EKE,'lon');
AVISO.lat=ncread(fid.EKE,'lat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EKE drag (diffusion) from Klocker and Abernathey
alpha=real(sqrt(1-f.^2/omega^2));
L=c(:,:,1)./f+40e3; % Eddy length scale
L(L>125e3)=125e3;
k=omega*alpha./c(:,:,1);% Tide wave number
k(k>2*pi/100e3)=0;
EKE=interp2(AVISO.lon,AVISO.lat',AVISO.EKE',lon,lat')';
u_rms=sqrt(2*EKE);
gamma=0.15; % Mixing efficiency (determined by best fit of zonal average) 
rEKE=gamma*u_rms.*L.*k.^2;

bad=(repmat(abs(lat)'<10,[Nx 1]) & u_rms>0.25 | rEKE>r_max); % Need to deal with high EKE at equator
rEKE(bad)=r_max;
rEKE(rEKE<r_min)=r_min;
rEKE(isnan(rEKE))=0;

% Compute horizontal viscosity
Ax=rEKE.*k.^(-2);
Ax(isinf(Ax) | isnan(Ax) | H<100)=0;

% Wave drag from Jayne and St. Laurent
H2=interp2(AVISO.lon,AVISO.lat',AVISO.H2',lon,lat')';
N_bott=interp2(AVISO.lon,AVISO.lat',AVISO.N_bott',lon,lat')';
N_bott(N_bott>0.01)=0.01;
kappa=2*pi*1e-4; % Tuned factor from Jayne and St. Laureng, still works well though
rWAVE=1/2*kappa*N_bott.*H2.*H.^(-1);
rWAVE(H<100)=0;
rWAVE(rWAVE>r_max)=r_max;
rWAVE(isnan(rWAVE) | isinf(rWAVE))=0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write NetCDF files
% EKE drag
fid.r='r_EKE.nc';

currentFolder=pwd;
cd(folder.r);

mode = netcdf.getConstant('CLOBBER');
mode = bitor(mode,netcdf.getConstant('NETCDF4'));
temp=netcdf.create(fid.r,mode);

% define dimensions
netcdf.defDim(temp,'x',Nx);
netcdf.defDim(temp,'y',Ny);
netcdf.defDim(temp,'mode',Nm);
netcdf.endDef(temp);

% Close file and return to data directory
netcdf.close(temp);
cd(currentFolder);

% Create and write fields
nccreate([folder.r,fid.r],'lon','Dimensions',{'x',Nx});
nccreate([folder.r,fid.r],'lat','Dimensions',{'y',Ny});
nccreate([folder.r,fid.r],'r','Dimensions',{'x',Nx,'y',Ny,'mode',Nm});

ncwrite([folder.r,fid.r],'lon',lon);
ncwrite([folder.r,fid.r],'lat',lat);
ncwrite([folder.r,fid.r],'r',repmat(rEKE,[1 1 Nm]));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wave drag only
fid.r='r_wave.nc';

currentFolder=pwd;
cd(folder.r);
	
mode = netcdf.getConstant('CLOBBER');
mode = bitor(mode,netcdf.getConstant('NETCDF4'));
temp=netcdf.create(fid.r,mode);

% define dimensions
netcdf.defDim(temp,'x',Nx);
netcdf.defDim(temp,'y',Ny);
netcdf.defDim(temp,'mode',Nm);
netcdf.endDef(temp);

% Close file and return to data directory
netcdf.close(temp);
cd(currentFolder);

% Create and write fields
nccreate([folder.r,fid.r],'lon','Dimensions',{'x',Nx});
nccreate([folder.r,fid.r],'lat','Dimensions',{'y',Ny});
nccreate([folder.r,fid.r],'r','Dimensions',{'x',Nx,'y',Ny,'mode',Nm});

ncwrite([folder.r,fid.r],'lon',lon);
ncwrite([folder.r,fid.r],'lat',lat);
ncwrite([folder.r,fid.r],'r',repmat(rWAVE,[1 1 Nm]));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Total drag
fid.r='r.nc';

currentFolder=pwd;
cd(folder.r);

mode = netcdf.getConstant('CLOBBER');
mode = bitor(mode,netcdf.getConstant('NETCDF4'));
temp=netcdf.create(fid.r,mode);

% define dimensions
netcdf.defDim(temp,'x',Nx);
netcdf.defDim(temp,'y',Ny);
netcdf.defDim(temp,'mode',Nm);
netcdf.endDef(temp);

% Close file and return to data directory
netcdf.close(temp);
cd(currentFolder);

% Create and write fields
nccreate([folder.r,fid.r],'lon','Dimensions',{'x',Nx});
nccreate([folder.r,fid.r],'lat','Dimensions',{'y',Ny});
nccreate([folder.r,fid.r],'r','Dimensions',{'x',Nx,'y',Ny,'mode',Nm});

ncwrite([folder.r,fid.r],'lon',lon);
ncwrite([folder.r,fid.r],'lat',lat);
ncwrite([folder.r,fid.r],'r',repmat(rWAVE+rEKE,[1 1 Nm]));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Horizontal viscosity
fid.r='Ax.nc';

currentFolder=pwd;
cd(folder.r);

mode = netcdf.getConstant('CLOBBER');
mode = bitor(mode,netcdf.getConstant('NETCDF4'));
temp=netcdf.create(fid.r,mode);

% define dimensions
netcdf.defDim(temp,'x',Nx);
netcdf.defDim(temp,'y',Ny);
netcdf.defDim(temp,'mode',Nm);
netcdf.endDef(temp);

% Close file and return to data directory
netcdf.close(temp);
cd(currentFolder);

% Create and write fields
nccreate([folder.r,fid.r],'lon','Dimensions',{'x',Nx});
nccreate([folder.r,fid.r],'lat','Dimensions',{'y',Ny});
nccreate([folder.r,fid.r],'Ax','Dimensions',{'x',Nx,'y',Ny,'mode',Nm});

ncwrite([folder.r,fid.r],'lon',lon);
ncwrite([folder.r,fid.r],'lat',lat);
ncwrite([folder.r,fid.r],'Ax',repmat(Ax,[1 1 Nm]));
