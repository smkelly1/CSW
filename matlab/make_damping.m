% Compute drag parameterizations for Coupled Shallow Water model (CSW)
clear

ii=complex(0,1);
warning off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set input parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
res=25; 
Nm=8; % Max
strat='WOA';
bathy='SS';

r_min=0;
r_max=1/(24*3600); % Min and max values of r

nu_min=0;
nu_max=8000; % The decay-scale in hours is roughly 1/(nu_max*((2*pi/(175e3))^2))/3600 ~ 24 h

kappa_min=0;
kappa_max=8000; 

% Input files
fid.grid=['../../../18-6_grids/',num2str(res),'th_deg_',strat,'_',bathy,'_grid.nc'];
fid.EKE=['../../../AVISO/AVISO.nc'];
fid.delta_c2=['../../../HYCOM/2015_c.nc'];

% Output files
fid.r='../../r.nc';
fid.nu='../../nu.nc';
fid.kappa='../../kappa.nc';

% Flags
flag.r=1;
flag.nu=1;
flag.kappa=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load grid
H=ncread(fid.grid,'H');
c=ncread(fid.grid,'c');
f=abs(ncread(fid.grid,'f'));
lon=ncread(fid.grid,'lon');
lat=ncread(fid.grid,'lat');
[Nx Ny]=size(H);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write spatially variable diffusion
if flag.kappa
	% Load HYCOM
	HYCOM.lon=ncread(fid.delta_c2,'lon');
	HYCOM.lat=ncread(fid.delta_c2,'lat');
	%HYCOM.c=ncread(fid.delta_c2,'c');
	%HYCOM.delta_c2=ncread(fid.delta_c2,'delta_c2');
	HYCOM.delta_c2=ncread(fid.delta_c2,'delta_c_var');

	% Compute variance and interpolate
	delta_c2=interp2(HYCOM.lon,HYCOM.lat',HYCOM.delta_c2',lon,lat')';

	% Compute diffusivity
	L=c(:,:,1)./f+40e3; % Eddy length scale
	L(L>125e3)=125e3;
	kappa=1/2*delta_c2.*L./c(:,:,1);

	% Filter high diffusivities
	kappa(kappa>kappa_max)=kappa_max;
	kappa(kappa<kappa_min)=kappa_min;
	kappa(isinf(kappa) | isnan(kappa) | H<100)=0;

	% Initiate file
	mode = netcdf.getConstant('CLOBBER');
	mode = bitor(mode,netcdf.getConstant('NETCDF4'));
	temp=netcdf.create(fid.kappa,mode);
	netcdf.close(temp);

	% Write fields
	nccreate(fid.kappa,'lon','Dimensions',{'x',Nx});
	nccreate(fid.kappa,'lat','Dimensions',{'y',Ny});
	nccreate(fid.kappa,'kappa','Dimensions',{'x',Nx,'y',Ny,'mode',Nm});

	ncwrite(fid.kappa,'lon',lon);
	ncwrite(fid.kappa,'lat',lat);
	ncwrite(fid.kappa,'kappa',repmat(kappa,[1 1 Nm]));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write spatially variable viscosity
if flag.nu
	% Load AVISO
	AVISO.EKE=ncread(fid.EKE,'EKE');
	AVISO.lon=ncread(fid.EKE,'lon');
	AVISO.lat=ncread(fid.EKE,'lat');

	% EKE drag (diffusion) from Klocker and Abernathey
	L=c(:,:,1)./f+40e3; % Eddy length scale
	L(L>125e3)=125e3;
	EKE=interp2(AVISO.lon,AVISO.lat',AVISO.EKE',lon,lat')';
	u_rms=sqrt(2*EKE);
	gamma=0.15; % Mixing efficiency (determined by best fit of zonal average) 
	nu=gamma*u_rms.*L;

	bad=(repmat(abs(lat)'<10,[Nx 1]) & u_rms>0.25 | nu>nu_max); % Need to deal with high EKE at equator
	nu(bad)=nu_max;
	nu(nu<nu_min)=nu_min;
	nu(isinf(nu) | isnan(nu) | H<100)=0;

	% Initiate file
	mode = netcdf.getConstant('CLOBBER');
	mode = bitor(mode,netcdf.getConstant('NETCDF4'));
	temp=netcdf.create(fid.nu,mode);
	netcdf.close(temp);

	% Write fields
	nccreate(fid.nu,'lon','Dimensions',{'x',Nx});
	nccreate(fid.nu,'lat','Dimensions',{'y',Ny});
	nccreate(fid.nu,'nu','Dimensions',{'x',Nx,'y',Ny,'mode',Nm});

	ncwrite(fid.nu,'lon',lon);
	ncwrite(fid.nu,'lat',lat);
	ncwrite(fid.nu,'nu',repmat(nu,[1 1 Nm]));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wave drag from Jayne and St. Laurent
if flag.r
	% Load AVISO
	AVISO.H2=ncread(fid.EKE,'H2');
	AVISO.N_bott=ncread(fid.EKE,'N_bott');
	AVISO.lon=ncread(fid.EKE,'lon');
	AVISO.lat=ncread(fid.EKE,'lat');

	H2=interp2(AVISO.lon,AVISO.lat',AVISO.H2',lon,lat')';
	N_bott=interp2(AVISO.lon,AVISO.lat',AVISO.N_bott',lon,lat')';
	N_bott(N_bott>0.01)=0.01;
	kappa0=2*pi*1e-4; % Tuned factor from Jayne and St. Laurent, still works well though
	r=1/2*kappa0*N_bott.*H2.*H.^(-1);

	r(H<100)=0;
	r(r>r_max)=r_max;
	r(isnan(r) | isinf(r))=0;

	% Initiate file
	mode = netcdf.getConstant('CLOBBER');
	mode = bitor(mode,netcdf.getConstant('NETCDF4'));
	temp=netcdf.create(fid.r,mode);
	netcdf.close(temp);

	% Write fields
	nccreate(fid.r,'lon','Dimensions',{'x',Nx});
	nccreate(fid.r,'lat','Dimensions',{'y',Ny});
	nccreate(fid.r,'r','Dimensions',{'x',Nx,'y',Ny,'mode',Nm});

	ncwrite(fid.r,'lon',lon);
	ncwrite(fid.r,'lat',lat);
	ncwrite(fid.r,'r',repmat(r,[1 1 Nm]));
end

