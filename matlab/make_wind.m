% Create NetCDF wind file
clear;

fid.MERRA2='./wind_pre_processing/out_HYCOM.nc';
fid.wind='HYCOM_CSW.nc';
folder.wind='../../';
cutoff=0.8;
equator=[10 15];
antarctic=[-70 -65];

% Load MERRA2 wind
lon=ncread(fid.MERRA2,'lon');        
lat=ncread(fid.MERRA2,'lat');        
taux=ncread(fid.MERRA2,'TAUX');
tauy=ncread(fid.MERRA2,'TAUY');
[Nx,Ny,Nt]=size(taux);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fourier filter the wind

% 24-h 1/2 cosine taper
mask=ones([1 1 Nt]);
mask(1,1,1:24)=1/2-1/2*cos((1:24)/24*pi);
mask(1,1,Nt-23:Nt)=1/2+1/2*cos((Nt-23:Nt)/24*pi);
mask=repmat(mask,[Nx Ny 1]);
taux=taux.*mask;
tauy=tauy.*mask;

dt=3600;
omegaN=2*pi/(2*dt);
dw=2*pi/(Nt*dt);
omega=-omegaN:dw:omegaN-dw;
omega=fftshift(omega);

omega=repmat(permute(omega,[3 1 2]),[Nx Ny 1]);
f=repmat(sw_f(lat)',[Nx 1 Nt]);

tmp=fft(taux,[],3);
tmp(abs(omega)<cutoff*abs(f))=0;
taux=real(ifft(tmp,[],3));

tmp=fft(tauy,[],3);
tmp(abs(omega)<cutoff*abs(f))=0;
tauy=real(ifft(tmp,[],3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notch the equator and 60 south
mask=ones([1 Ny]);
mask(abs(lat)<equator(1))=0;

good=equator(1)<=lat & lat<=equator(2);
mask(good)=1/2-1/2*cos((lat(good)-equator(1))/(equator(2)-equator(1))*pi);

good=-equator(2)<=lat & lat<=-equator(1);
mask(good)=1/2-1/2*cos((lat(good)+equator(1))/(equator(2)-equator(1))*pi);

good=antarctic(1)<=lat & lat<=antarctic(2);
mask(good)=1/2-1/2*cos((lat(good)-antarctic(1))/(antarctic(2)-antarctic(1))*pi);

good=lat<=antarctic(1);
mask(good)=0;

mask=repmat(mask,[Nx 1 Nt]);
taux=taux.*mask;
tauy=tauy.*mask;


%% start netcdf file
disp('Writing wind');

name=[folder.wind,fid.wind];

mode = netcdf.getConstant('CLOBBER'); % Careful, this overwrites any existing file with the same name! 
mode = bitor(mode,netcdf.getConstant('NETCDF4'));
temp=netcdf.create(name,mode);
netcdf.close(temp);

% The array sizes are useful below
[Nx,Ny,Nt]=size(taux);

% write the INPUT fields to files. 
nccreate(name,'lon','Dimensions',{'x',Nx});
ncwrite(name,'lon',lon);

nccreate(name,'lat','Dimensions',{'y',Ny});
ncwrite(name,'lat',lat);

nccreate(name,'TAUX','Dimensions',{'x',Nx,'y',Ny,'time',inf});
ncwrite(name,'TAUX',taux);

nccreate(name,'TAUY','Dimensions',{'x',Nx,'y',Ny,'time',inf});
ncwrite(name,'TAUY',tauy);

