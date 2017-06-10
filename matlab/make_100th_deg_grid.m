% Make 100th deg grid from 50th degree grid


% If this is run on MSI, add paths that are sometimes lost using qsub
MSI=1;

if MSI
	addpath(genpath('/home/kellys/smkelly/software/matlab_libraries/sam_ware'))
	addpath(genpath('/home/kellys/smkelly/software/matlab_libraries/seawater'))
	addpath(genpath('/home/kellys/smkelly/software/data_products/SS_topo'))
	addpath(genpath('/home/kellys/smkelly/software/data_products/OTPS'))
	addpath(genpath('/home/kellys/smkelly/software/data_products/WOA13'))
	folder='~/simulations/SWOT/global_grids/';
else
	folder='~/simulations/SWOT/global_grids/';
end

% Variables that define the grid
dx=1/100; % Degrees
fid='100th_deg';
fid_old='50th_deg';

% Calculation parameters 
Nm=8; % Want Nm=8 eventually

% Grid parameters (probably leave these alone)
H_min=16; % want H_min=16 eventually
H_max=6000;
latlims=[-80 66];
a0=6371e3; % radius of Earth

% Define output grid
dy=dx;
lon=dx/2:dx:360; 
lat=latlims(1)+dy/2:dy:latlims(2);

% Add buffer points
lon=[lon(1)-dx; lon'; lon(end)+dx];
lat=[lat(1)-dy; lat'; lat(end)+dy];

Nx=length(lon);
Ny=length(lat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start the grid file
disp('Glueing grid');

folder='./';
fid_grid=[fid,'_grid.nc'];
name=[folder,fid_grid];

%Create file
mode = netcdf.getConstant('CLOBBER');
mode = bitor(mode,netcdf.getConstant('NETCDF4'));
ncid=netcdf.create(name,mode);

% define dimensions
netcdf.defDim(ncid,'x',Nx);
netcdf.defDim(ncid,'y',Ny);
netcdf.defDim(ncid,'mode',Nm);
netcdf.endDef(ncid)

% Close file and return to data directory
netcdf.close(ncid);

% Create variables
% Dimensions
nccreate(name,'Nx');
nccreate(name,'Ny');
nccreate(name,'Nm');

% Grid 
nccreate(name,'dy');
nccreate(name,'dx');
nccreate(name,'dz');

nccreate(name,'lon','Dimensions',{'x',Nx});
nccreate(name,'lat','Dimensions',{'y',Ny}); 

% Data Variables
nccreate(name,'H','Dimensions',{'x',Nx,'y',Ny});
nccreate(name,'f','Dimensions',{'x',Nx,'y',Ny});

nccreate(name,'c','Dimensions',{'x',Nx,'y',Ny,'mode',Nm});
nccreate(name,'phi_surf','Dimensions',{'x',Nx,'y',Ny,'mode',Nm});
nccreate(name,'phi_bott','Dimensions',{'x',Nx,'y',Ny,'mode',Nm});

nccreate(name,'T_x','Dimensions',{'x',Nx,'y',Ny,'mode',Nm,'mode',Nm});
nccreate(name,'T_y','Dimensions',{'x',Nx,'y',Ny,'mode',Nm,'mode',Nm});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Write the variables we already know
ncwrite(name,'Nx',Nx);
ncwrite(name,'Ny',Ny);
ncwrite(name,'Nm',Nm);

ncwrite(name,'lon',lon);
ncwrite(name,'lat',lat);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Start loading in 50th degree grid and interpolating
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid_grid=[fid_old,'_grid.nc'];
name_old=[folder,fid_grid];

% get old grid
lon_old=ncread(name_old,'lon');
lat_old=ncread(name_old,'lat');

% Depth
H=ncread(name_old,'H');
H=interp2(lon_old,lat_old',H.',lon,lat').';
H(H<H_min)=0;
H(H>H_max)=H_max;
ncwrite(name,'H',H);

% f
tmp=ncread(name_old,'f');
tmp=interp2(lon_old,lat_old',tmp.',lon,lat').';
ncwrite(name,'f',tmp);

% c, phi_surf, and phi_bott
for n=1:Nm
	tmp=ncread(name_old,'c',[1 1 n],[Inf Inf 1]);
	tmp=interp2(lon_old,lat_old',tmp.',lon,lat').';
	tmp(H==0)=0;
	ncwrite(name,'c',tmp,[1 1 n]);
	
	tmp=ncread(name_old,'phi_surf',[1 1 n],[Inf Inf 1]);
	tmp=interp2(lon_old,lat_old',tmp.',lon,lat').';
	tmp(H==0)=0;
	ncwrite(name,'phi_surf',tmp,[1 1 n]);
	
	
	tmp=ncread(name_old,'phi_bott',[1 1 n],[Inf Inf 1]);
	tmp=interp2(lon_old,lat_old',tmp.',lon,lat').';
	tmp(H==0)=0;
	ncwrite(name,'phi_bott',tmp,[1 1 n]);
end

% The coupling coefficients
for n=1:Nm
	for m=1:Nm
		tmp=ncread(name_old,'T_x',[1 1 m n],[Inf Inf 1 1]);
		tmp=interp2(lon_old,lat_old',tmp.',lon,lat').';
		tmp(H==0)=0;
		ncwrite(name,'T_x',tmp,[1 1 m n]);
		
		tmp=ncread(name_old,'T_y',[1 1 m n],[Inf Inf 1 1]);
		tmp=interp2(lon_old,lat_old',tmp.',lon,lat').';
		tmp(H==0)=0;
		ncwrite(name,'T_y',tmp,[1 1 m n]);
	end
end



