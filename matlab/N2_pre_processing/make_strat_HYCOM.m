clear

folder='/home/smkelly/data/data_products/HYCOM/monthly_averages/';

fid=['hycom_GLBS0.50_53X_archMN.1995_01_temp.nc'];
lon=ncread([folder,fid],'lon');
lat=ncread([folder,fid],'lat');
z=ncread([folder,fid],'depth');
year=1994:2015;

Nx=length(lon);
Ny=length(lat);
Nz=length(z);
Nt=length(year)*12;

% Load data
T=NaN([Nx Ny Nz Nt]);
S=NaN([Nx Ny Nz Nt]);
ind=1;
for i=1:length(year)
    for j=1:12
        fid=['hycom_GLBS0.50_53X_archMN.',num2str(year(i)),'_',num2str(j,'%02d'),'_temp.nc'];
        T(:,:,:,ind)=ncread([folder,fid],'water_temp');

        fid=['hycom_GLBS0.50_53X_archMN.',num2str(year(i)),'_',num2str(j,'%02d'),'_sali.nc'];
        S(:,:,:,ind)=ncread([folder,fid],'salinity');

        ind=ind+1;
    end
end

%% Find median temperature and salinity at each depth 
T0=nanmedian(T,4);
S0=nanmedian(S,4);

%%
z_N2=(z(2:end)+z(1:end-1))/2;
Nz=length(z_N2);

N2=NaN([Nx Ny Nz]);
for j=1:Ny
    p=sw_pres(z,lat(j));
    for i=1:Nx
        N2(i,j,:)=sw_bfrq(squeeze(S0(i,j,:)),squeeze(T0(i,j,:)),p,lat(j));
    end
end

%% Create a netCDF file  
fid='strat_HYCOM.nc';
mode = netcdf.getConstant('CLOBBER'); % Erase existing file
mode = bitor(mode,netcdf.getConstant('64BIT_OFFSET')); % Allow giant files
temp=netcdf.create(fid,mode);
netcdf.close(temp);

% Write grid
nccreate(fid,'lon','Dimensions',{'lon',Nx});
ncwriteatt(fid,'lon', '_CoordinateAxisType','lon');
ncwrite(fid,'lon',lon);

nccreate(fid,'lat','Dimensions',{'lat',Ny});
ncwriteatt(fid,'lat', '_CoordinateAxisType','lat');
ncwrite(fid,'lat',lat);

nccreate(fid,'depth','Dimensions',{'depth',Nz});
ncwriteatt(fid,'depth', '_CoordinateAxisType','depth');
ncwrite(fid,'depth',z_N2);

% Write N2
nccreate(fid,'N2','Dimensions',{'lon',Nx,'lat',Ny,'depth',Nz});
ncwriteatt(fid,'N2','units','1/s^2');
ncwriteatt(fid,'N2','long_name','Buoyancy frequency squared');
ncwrite(fid,'N2',N2);








