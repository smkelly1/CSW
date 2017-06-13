%% Coupled Shallow Water model (CSW)
clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set input parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MSI=1; 

fid.grid='../../../17-6_global_grids/50th_deg_grid.nc';
fid.tides='../../50th_deg_tides.nc';

Nc=1;               % Number of tidal constituents 
H_min=16;			% Set minimum depth

% Add paths on MSI that are sometimes lost using qsub
if MSI
	addpath(genpath('/home/kellys/smkelly/software/matlab_libraries/sam_ware'))
	addpath(genpath('/home/kellys/smkelly/software/matlab_libraries/seawater'))
	addpath(genpath('/home/kellys/smkelly/software/data_products/OTPS'))
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                    Technical code below here
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Loading grid');

Nx=ncread(fid.grid,'Nx')-2;
Ny=ncread(fid.grid,'Ny')-2;

H=ncread(fid.grid,'H'); H=H(2:end-1,2:end-1);
lon=ncread(fid.grid,'lon'); lon=lon(2:end-1);
lat=ncread(fid.grid,'lat'); lat=lat(2:end-1); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Creating NetCDF input file');

% start netcdf file
mode = netcdf.getConstant('CLOBBER');
mode = bitor(mode,netcdf.getConstant('NETCDF4'));
ncid=netcdf.create(fid.tides,mode);

% define dimensions
netcdf.defDim(ncid,'x',Nx);
netcdf.defDim(ncid,'y',Ny);
netcdf.defDim(ncid,'con',Nc);
netcdf.endDef(ncid)

% Close file and return to data directory
netcdf.close(ncid);

% Create variables and write small variables
nccreate(fid.tides,'Nx');
ncwrite(fid.tides,'Nx',Nx);

nccreate(fid.tides,'Ny');
ncwrite(fid.tides,'Ny',Ny);

nccreate(fid.tides,'Nc');
ncwrite(fid.tides,'Nc',Nc);

nccreate(fid.tides,'lon','Dimensions',{'x',Nx});
ncwrite(fid.tides,'lon',lon);

nccreate(fid.tides,'lat','Dimensions',{'y',Ny});      
ncwrite(fid.tides,'lat',lat);  

% Initiate tides in NetCDF file 
nccreate(fid.tides,'omega','Dimensions',{'con',Nc});
nccreate(fid.tides,'Ur','Dimensions',{'x',Nx,'y',Ny,'con',Nc},'datatype','single');
nccreate(fid.tides,'Ui','Dimensions',{'x',Nx,'y',Ny,'con',Nc},'datatype','single');
nccreate(fid.tides,'Vr','Dimensions',{'x',Nx,'y',Ny,'con',Nc},'datatype','single');
nccreate(fid.tides,'Vi','Dimensions',{'x',Nx,'y',Ny,'con',Nc},'datatype','single');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shallow water mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Masking tides in shallow regions');
  
% Mask out shallow regions
H(H<H_min)=0;

% Mask out the Caspian Sea (seems to create a CFL instability)
bad.lon=45<lon & lon<60;
bad.lat=35<lat & lat<50;
H(bad.lon,bad.lat)=0;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Surface-tide forcing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Loading TPXO surface tides');

ii=complex(0,1);
warning off;
	
% Obtain TPXO surface tide velocities
[U,V,~,omega]=TPXO(lon,lat,[1],datenum(2015,1,1));

if length(omega)~=Nc
    disp('Error, wrong number of constituents')
end

% Wrap the tidal fields at the Grand Meridian
if isnan(U(1,1,1))
	U(1,:,:)=(U(end,:,:)+U(2,:,:))/2;
	V(1,:,:)=(V(end,:,:)+V(2,:,:))/2;
end

% Convert transport to velocity 
for k=1:Nc
	U(:,:,k)=U(:,:,k)./H;
	V(:,:,k)=V(:,:,k)./H;   
end

% Remove bad data
U(isinf(U) | isnan(U))=0+ii*0;
V(isinf(V) | isnan(V))=0+ii*0;

% Limit surface tides to 1 m/s
thresh=1; 
fast=thresh<abs(U);
U(fast)=U(fast).*thresh./abs(U(fast));
fast=thresh<abs(V);
V(fast)=V(fast).*thresh./abs(V(fast));

% Round the tidal frequency and write to file
period=(2*pi/omega)/3600; % in hours
period=round(period*100)/100; % Round to second decimal place (i.e., 12.42 for the M2 tide)
omega=(2*pi)/(period*3600);

% Write to file
ncwrite(fid.tides,'omega',omega);
ncwrite(fid.tides,'Ur',real(U));
ncwrite(fid.tides,'Ui',imag(U));
ncwrite(fid.tides,'Vr',real(V));
ncwrite(fid.tides,'Vi',imag(V));


