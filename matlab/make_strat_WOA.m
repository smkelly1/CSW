clear

name='~/software/data_products/WOA13/woa13_decav_t00an01v2.dat';

fid=fopen(name,'r');
t0=fscanf(fid,'%g');
fclose(fid);
t0=reshape(t0,[360 180 102]);
t0(t0<0 | 99<t0)=NaN;

name='~/software/data_products/WOA13/woa13_decav_s00an01v2.dat';

fid=fopen(name,'r');
s0=fscanf(fid,'%g');
fclose(fid);
s0=reshape(s0,[360 180 102]);
s0(s0<0 | 99<s0)=NaN;

%%
z0=[0:5:100 125:25:500 550:50:2000 2100:100:5500];

strat.lat=[-89.5:1:89.5];
strat.lon=[.5:1:359.5];
strat.s=s0;
strat.t=t0;
strat.z0=z0;
for j=1:size(t0,2);
    p0=sw_pres(z0',strat.lat(j));
    for i=1:size(t0,1);
       strat.N2(i,j,:)=sw_bfrq(squeeze(s0(i,j,:)),squeeze(t0(i,j,:)),p0,strat.lat(j));
    end
end
strat.z=(z0(2:end)+z0(1:end-1))/2;

%save strat strat;

%% Create a netCDF file  
[Nx Ny Nz]=size(strat.N2);

fid='strat_WOA.nc';
mode = netcdf.getConstant('CLOBBER'); % Erase existing file
mode = bitor(mode,netcdf.getConstant('64BIT_OFFSET')); % Allow giant files
temp=netcdf.create(fid,mode);
netcdf.close(temp);

% Write grid
nccreate(fid,'lon','Dimensions',{'lon',Nx});
ncwriteatt(fid,'lon', '_CoordinateAxisType','lon');
ncwrite(fid,'lon',strat.lon);

nccreate(fid,'lat','Dimensions',{'lat',Ny});
ncwriteatt(fid,'lat', '_CoordinateAxisType','lat');
ncwrite(fid,'lat',strat.lat);

nccreate(fid,'depth','Dimensions',{'depth',Nz});
ncwriteatt(fid,'depth', '_CoordinateAxisType','depth');
ncwrite(fid,'depth',strat.z);

% Write N2
nccreate(fid,'N2','Dimensions',{'lon',Nx,'lat',Ny,'depth',Nz});
ncwriteatt(fid,'N2','units','1/s^2');
ncwriteatt(fid,'N2','long_name','Buoyancy frequency squared');
ncwrite(fid,'N2',strat.N2);



