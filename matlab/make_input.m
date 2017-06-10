%% Coupled Shallow Water model (CSW)
clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set input parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Setting paths and parameters');

MSI=0; 
fid.grid='../../../17-6_global_grids/10th_deg_grid.nc';
fid.in='../../10th_deg_in.nc';

Nm=2;				% Number of modes to include in the input file
dt=12.42*3600/100;	% Approximate time step for sponge layer
H_min=16;			% Turn off forcing and mask for shallow water


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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Loading grid');

Nx=ncread(fid.grid,'Nx')-2;
Ny=ncread(fid.grid,'Ny')-2;

dx=ncread(fid.grid,'dx');
dy=ncread(fid.grid,'dy');

H=ncread(fid.grid,'H'); 
f=ncread(fid.grid,'f');
lon=ncread(fid.grid,'lon');
lat=ncread(fid.grid,'lat');    

c=ncread(fid.grid,'c'); 
c=c(:,:,1:Nm);
phi_bott=ncread(fid.grid,'phi_bott');
phi_bott=phi_bott(:,:,1:Nm);

% Set depth to zero where c_1=0;
H(c(:,:,1)==0)=0;

% Set depth to zero if it's less than the minimum depth
H(H<H_min)=0;

% Mask out the Caspian Sea (seems to create a CFL instability)
bad.lon=45<lon & lon<60;
bad.lat=35<lat & lat<50;
H(bad.lon,bad.lat)=0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Land/boundary masks (Note: we don't use border data for H, f, or c here)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Generating land mask');

% Initiate with 1 (good/ocean) everywhere
mask.u=ones([Nx+1 Ny]);
mask.v=ones([Nx Ny+1]);
mask.p=ones([Nx Ny]);

% Set velocities on land at at boundaries to zero
mask.u(1:end-1,:)=min(H(2:end-1,2:end-1)~=0,mask.u(1:end-1,:));
mask.u(2:end,:)=min(H(2:end-1,2:end-1)~=0,mask.u(2:end,:));

mask.v(:,1:end-1)=min(H(2:end-1,2:end-1)~=0,mask.v(:,1:end-1));
mask.v(:,2:end)=min(H(2:end-1,2:end-1)~=0,mask.v(:,2:end));

% Add in mode dimension
mask.u=repmat(mask.u,[1 1 Nm]);
mask.v=repmat(mask.v,[1 1 Nm]);
mask.p=repmat(mask.p,[1 1 Nm]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add Coriolis "wave resolution" mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Generating wave resolution mask');

a=6371e3; % Radius of the earth
dx_f=repmat(a*cos(lat(2:end-1)'*pi/180)*dx/360*2*pi,[Nx 1]); % grid spacing in m at each latitude

f_crit=abs(f);
c_crit=repmat(f_crit(2:end-1,2:end-1).*dx_f/2,[1 1 Nm]);

maskf=1-c(2:end-1,2:end-1,:)./(2*c_crit); % "High res" Adcroft 1999
maskf(maskf<0)=0;
maskf=1-maskf*dt/(1*3600); % Create a sponge with a 1 hour decay time scale

% Use the Coriolis mask if it provides more damping
mask.p=min(maskf,mask.p);

mask.u(1:end-1,:,:)=min(maskf,mask.u(1:end-1,:,:));
mask.u(2:end,:,:)=min(maskf,mask.u(2:end,:,:));

mask.v(:,1:end-1,:)=min(maskf,mask.v(:,1:end-1,:));
mask.v(:,2:end,:)=min(maskf,mask.v(:,2:end,:));    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Surface-tide forcing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Loading TPXO surface tides');

ii=complex(0,1);
warning off;

% Convert dx and longitude to radians
a=6371e3; % Radius of the earth
lat_R=lat(2:end-1)/180*pi;
dx_R=dx/180*pi;
dy_R=dy/180*pi;

% Compute Reduced slope
dHdx=1./(a*repmat(cos(lat_R)',[Nx 1])).*(H(3:end,2:end-1)-H(1:end-2,2:end-1))/(2*dx_R)./H(2:end-1,2:end-1);
dHdy=1/a*(H(2:end-1,3:end)-H(2:end-1,1:end-2))/(2*dy_R)./H(2:end-1,2:end-1);
	
% Obtain TPXO surface tide velocities
Nc=1;
[U V junk ITGF.omega]=TPXO(lon(2:end-1),lat(2:end-1),[1],datenum(2015,1,1));
if isnan(U(1,1,1))
	U(1,:,:)=(U(end,:,:)+U(2,:,:))/2;
	V(1,:,:)=(V(end,:,:)+V(2,:,:))/2;
end

% Round the tidal frequency
period=(2*pi/ITGF.omega)/3600; % in hours
period=round(period*100)/100; % Round to second decimal place (i.e., 12.42 for the M2 tide)
ITGF.omega=(2*pi)/(period*3600);

% Compute velocities
U=U./repmat(H(2:end-1,2:end-1),[1 1 Nc]);
V=V./repmat(H(2:end-1,2:end-1),[1 1 Nc]);
	
% Remove bad slopes and surface tides
U(isinf(U) | isnan(U))=0;
V(isinf(V) | isnan(V))=0;
dHdx(isnan(dHdx) | isinf(dHdx))=0;
dHdy(isnan(dHdy) | isinf(dHdy))=0;

% Limit surface tides to 1 m/s
thresh=1; 
fast=thresh<abs(U);
U(fast)=U(fast).*thresh./abs(U(fast));
fast=thresh<abs(V);
V(fast)=V(fast).*thresh./abs(V(fast));

% Ignore surface tides in less than 50 m depth
U(repmat(H(2:end-1,2:end-1),[1 1 Nc])<H_min)=0+ii*0;
V(repmat(H(2:end-1,2:end-1),[1 1 Nc])<H_min)=0+ii*0;
	
% Compute ITGF forcing
ITGF.p=U.*repmat(dHdx,[1 1 Nc])+V.*repmat(dHdy,[1 1 Nc]); 

% Multiply by phi_bottom 
ITGF.p=repmat(permute(ITGF.p,[1 2 4 3]),[1 1 Nm 1]).*repmat(phi_bott(2:end-1,2:end-1,:),[1 1 1 Nc]);

% Smooth at grid scale and ignore forcing in less than 50 m depth
for j=1:Nc
    for i=1:Nm
		ITGF.p(:,:,i,j)=AVE2D(ITGF.p(:,:,i,j),3);
    end
end

% Turn off tidal forcing in locations with insufficient wave resolution 
ITGF.p(repmat(mask.p,[1 1 1 Nc])~=1)=0+ii*0;	

% One last check for bad forcing
ITGF.p(isnan(ITGF.p) | isinf(ITGF.p))=0+ii*0;  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Write NetCDF file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Writing NetCDF input file');

% start netcdf file
mode = netcdf.getConstant('CLOBBER');
mode = bitor(mode,netcdf.getConstant('NETCDF4'));
ncid=netcdf.create(fid.in,mode);

% define dimensions
netcdf.defDim(ncid,'x',Nx);
netcdf.defDim(ncid,'xu',Nx+1);
netcdf.defDim(ncid,'xH',Nx+2);
netcdf.defDim(ncid,'y',Ny);
netcdf.defDim(ncid,'yv',Ny+1);
netcdf.defDim(ncid,'yH',Ny+2);
netcdf.defDim(ncid,'mode',Nm);
netcdf.defDim(ncid,'con',Nc);
netcdf.endDef(ncid)

% Close file and return to data directory
netcdf.close(ncid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create variables and write data

% Dimmensions  
nccreate(fid.in,'Nx');
ncwrite(fid.in,'Nx',Nx);

nccreate(fid.in,'Ny');
ncwrite(fid.in,'Ny',Ny);

nccreate(fid.in,'Nm');
ncwrite(fid.in,'Nm',Nm);


% Grid spacing
nccreate(fid.in,'dx');
ncwrite(fid.in,'dx',dx);

nccreate(fid.in,'dy');
ncwrite(fid.in,'dy',dy);


% Masks
nccreate(fid.in,'mask_u','Dimensions',{'xu',Nx+1,'y',Ny,'mode',Nm});
ncwrite(fid.in,'mask_u',mask.u);

nccreate(fid.in,'mask_v','Dimensions',{'x',Nx,'yv',Ny+1,'mode',Nm});
ncwrite(fid.in,'mask_v',mask.v);

nccreate(fid.in,'mask_p','Dimensions',{'x',Nx,'y',Ny,'mode',Nm});
ncwrite(fid.in,'mask_p',mask.p);


% Lat/lon for spherical grids
nccreate(fid.in,'lon','Dimensions',{'xH',Nx+2});
ncwrite(fid.in,'lon',lon);

nccreate(fid.in,'lat','Dimensions',{'yH',Ny+2});      
ncwrite(fid.in,'lat',lat);  


% Internal-tide Forcing
nccreate(fid.in,'ITGF_omega','Dimensions',{'con',Nc});
ncwrite(fid.in,'ITGF_omega',ITGF.omega);

nccreate(fid.in,'ITGF_pr','Dimensions',{'x',Nx,'y',Ny,'mode',Nm,'con',Nc});
ncwrite(fid.in,'ITGF_pr',real(ITGF.p));

nccreate(fid.in,'ITGF_pi','Dimensions',{'x',Nx,'y',Ny,'mode',Nm,'con',Nc});
ncwrite(fid.in,'ITGF_pi',imag(ITGF.p));

