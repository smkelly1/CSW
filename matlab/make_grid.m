% Make a netCDF grid file
clear

% Variables that define the grid
res=100; % Degrees, want 1/100 at some point

% Define the base grid
fid_in='/home/smkelly/experiments/NISKINE/CSW/22-12_grid/SS_WOA_grid.nc';

% Define the new grid goes
folder='/home/smkelly/experiments/NISKINE/CSW/22-12_grid/';
filename=[num2str(res),'th_deg_grid.nc'];
fid_out=[folder,filename]; 

% Grid parameters (probably leave these alone)
latlims=[-80 66];
a0=6371e3; % radius of Earth

% Define output grid
dx=1/res;
dy=dx;
lon=dx/2:dx:360;
lat=latlims(1)+dy/2:dy:latlims(2);

% Obtain bathymetry at full resolution
bathy.lon=ncread(fid_in,'lon');
bathy.lat=ncread(fid_in,'lat');

% Smooth to desired resolution (Note: you need to pad the Grand Meridian)
Ns=ceil((1/res)/median(diff(bathy.lon)))*3; % The extra factor of 3 is for stability

% Load, smooth, and interpolate depth
tmp=ncread(fid_in,'H');
tmp=[tmp(end-Ns+1:end,:);tmp;tmp(1:Ns,:)];
tmp=AVE2D(tmp,Ns);
tmp=tmp(Ns+1:end-Ns,:);

% Need to re-arrange grid from [-180, 180] to [0, 360]
right=(bathy.lon<0);
left=(bathy.lon>=0);
tmp=[tmp(left,:); tmp(right,:)];
bathy.lon=[bathy.lon(left); bathy.lon(right)+360];

H=interp2(bathy.lon,bathy.lat',tmp',lon',lat)';
H(H<0)=0; % Remove land
clear tmp

% Add one buffer point on each edge
H=[H(end,:); H; H(1,:)];
H=[H(:,1) H H(:,end)];

dlon=mean(diff(lon));
lon=[lon(1)-dlon lon lon(end)+dlon];

dlat=mean(diff(lat));
lat=[lat(1)-dlat lat lat(end)+dlat];

[Nx,Ny]=size(H);
tmp=ncinfo(fid_in,'c');
Nm=tmp.Dimensions(3).Length;

%% Start NetCDF file
mode=netcdf.getConstant('CLOBBER');
mode=bitor(mode,netcdf.getConstant('NETCDF4'));
ncid=netcdf.create(fid_out,mode);

% define dimensions
netcdf.defDim(ncid,'x',Nx);
netcdf.defDim(ncid,'y',Ny);
netcdf.defDim(ncid,'mode',Nm);
netcdf.endDef(ncid)
netcdf.close(ncid);

% Create variables
% Dimensions
nccreate(fid_out,'Nx');
nccreate(fid_out,'Ny');
nccreate(fid_out,'Nm');
nccreate(fid_out,'Ns');

% Grid 
nccreate(fid_out,'lon','Dimensions',{'x',Nx});
nccreate(fid_out,'lat','Dimensions',{'y',Ny}); 

% Data Variables
nccreate(fid_out,'H','Dimensions',{'x',Nx,'y',Ny});
nccreate(fid_out,'f','Dimensions',{'x',Nx,'y',Ny});

nccreate(fid_out,'c','Dimensions',{'x',Nx,'y',Ny,'mode',Nm});
nccreate(fid_out,'phi_surf','Dimensions',{'x',Nx,'y',Ny,'mode',Nm});
nccreate(fid_out,'phi_bott','Dimensions',{'x',Nx,'y',Ny,'mode',Nm});

% Write the variables we already know
ncwrite(fid_out,'Nx',Nx);
ncwrite(fid_out,'Ny',Ny);
ncwrite(fid_out,'Nm',Nm);
ncwrite(fid_out,'Ns',Ns);

ncwrite(fid_out,'lon',lon);
ncwrite(fid_out,'lat',lat);

ncwrite(fid_out,'H',H);
ncwrite(fid_out,'f',repmat(2*(7.292e-5)*sin(lat/180*pi),[Nx 1]));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Write c
tmp=ncread(fid_in,'c');

% Smooth to desired resolution (Note: you need to pad the Grand Meridian to get this right)
tmp=[tmp(end-Ns+1:end,:,:);tmp;tmp(1:Ns,:,:)];
for i=1:Nm
    tmp(:,:,i)=AVE2D(tmp(:,:,i),Ns);
end
tmp=tmp(Ns+1:end-Ns,:,:);

% Need to re-arrange grid from [-180, 180] to [0, 360]
tmp=[tmp(left,:,:); tmp(right,:,:)];

tmp2=zeros([Nx-2 Ny-2 Nm]);
for i=1:Nm
    tmp2(:,:,i)=interp2(bathy.lon,bathy.lat',tmp(:,:,i)',lon(2:end-1)',lat(2:end-1))';
end
clear tmp

% Add one buffer point on each edge
tmp2=[tmp2(end,:,:); tmp2; tmp2(1,:,:)];
tmp2=[tmp2(:,1,:) tmp2 tmp2(:,end,:)];

% write to NetCDF
ncwrite(fid_out,'c',tmp2);
clear tmp2

%% Write phi_surf
tmp=ncread(fid_in,'phi_surf');

% Smooth to desired resolution (Note: you need to pad the Grand Meridian to get this right)
tmp=[tmp(end-Ns+1:end,:,:);tmp;tmp(1:Ns,:,:)];
for i=1:Nm
    tmp(:,:,i)=AVE2D(tmp(:,:,i),Ns);
end
tmp=tmp(Ns+1:end-Ns,:,:);

% Need to re-arrange grid from [-180, 180] to [0, 360]
tmp=[tmp(left,:,:); tmp(right,:,:)];

tmp2=zeros([Nx-2 Ny-2 Nm]);
for i=1:Nm
    tmp2(:,:,i)=interp2(bathy.lon,bathy.lat',tmp(:,:,i)',lon(2:end-1)',lat(2:end-1))';
end
clear tmp

% Add one buffer point on each edge
tmp2=[tmp2(end,:,:); tmp2; tmp2(1,:,:)];
tmp2=[tmp2(:,1,:) tmp2 tmp2(:,end,:)];

% write to NetCDF
ncwrite(fid_out,'phi_surf',tmp2);
clear tmp2

%% Write phi_bott
tmp=ncread(fid_in,'phi_bott');

% Smooth to desired resolution (Note: you need to pad the Grand Meridian to get this right)
tmp=[tmp(end-Ns+1:end,:,:);tmp;tmp(1:Ns,:,:)];
for i=1:Nm
    tmp(:,:,i)=AVE2D(tmp(:,:,i),Ns);
end
tmp=tmp(Ns+1:end-Ns,:,:);

% Need to re-arrange grid from [-180, 180] to [0, 360]
tmp=[tmp(left,:,:); tmp(right,:,:)];

tmp2=zeros([Nx-2 Ny-2 Nm]);
for i=1:Nm
    tmp2(:,:,i)=interp2(bathy.lon,bathy.lat',tmp(:,:,i)',lon(2:end-1)',lat(2:end-1))';
end
clear tmp

% Add one buffer point on each edge
tmp2=[tmp2(end,:,:); tmp2; tmp2(1,:,:)];
tmp2=[tmp2(:,1,:) tmp2 tmp2(:,end,:)];

% write to NetCDF
ncwrite(fid_out,'phi_bott',tmp2);
clear tmp2

