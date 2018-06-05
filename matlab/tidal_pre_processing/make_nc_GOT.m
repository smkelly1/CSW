clear; 

% Specify file names
ncfile='m2UV_GOT.nc';
dfile='m2UV.d';

% start netcdf file
mode = netcdf.getConstant('CLOBBER'); % Erase existing file
mode = bitor(mode,netcdf.getConstant('64BIT_OFFSET')); % Allow giant files
temp=netcdf.create(ncfile,mode); % Create file
netcdf.close(temp);
ncwriteatt(ncfile,'/','Source','GOT: M2 transport, R Ray 17-Aug-07') % Add metadata 

% Read ASCII file
fid=fopen(dfile);
var_name={'UAMP','UPHA','VAMP','VPHA'};
for i=1:4
    
    % Process header
    name=fgetl(fid); % Name of variable	
    author=fgetl(fid); 
    dims=str2num(fgetl(fid)); % Matrix dimensions
    latlims=str2num(fgetl(fid)); % Latitude limits
    lonlims=str2num(fgetl(fid)); % Longitude limits
    bad_flag=str2num(fgetl(fid));bad_flag=bad_flag(1); % Bad data identifier
    dat_format=fgetl(fid); % Data format
    
    % Read data
    Ny=dims(1); 
    Nx=dims(2);  
    dat=fscanf(fid,'%f',[Nx Ny]);
    dat(dat==bad_flag)=NaN; % Change bad data to NaN
    
    % Finish reading the last line
    blank=fgetl(fid);
    
    % Create and write grid
    if i==1
        dlat=diff(latlims)/(Ny-1);
        lat=latlims(1):dlat:latlims(2);
        
        dlon=diff(lonlims)/(Nx-1);
        lon=lonlims(1):dlat:lonlims(2);
        lon(end+1)=lon(1)+360;
        lon=(lon(1:end-1)+lon(2:end))/2;
        %lon(end+1)=lon(end)+dlon;
        
        nccreate(ncfile,'lon','Dimensions',{'lon',Nx});
        ncwriteatt(ncfile,'lon', '_CoordinateAxisType','lon');
        ncwrite(ncfile,'lon',lon);
        
        nccreate(ncfile,'lat','Dimensions',{'lat',Ny});
        ncwriteatt(ncfile,'lat', '_CoordinateAxisType','lat');
        ncwrite(ncfile,'lat',lat);
    end
        
    % Write data
    if i<3
        nccreate(ncfile,var_name{i},'Dimensions',{'lon',Nx,'lat',Ny});
        dat(end+1,:)=dat(1,:);
        dat=(dat(1:end-1,:)+dat(2:end,:))/2;
    else
        nccreate(ncfile,var_name{i},'Dimensions',{'lon',Nx,'lat',Ny+1});
        dat=[dat(:,1) dat dat(:,end)];
        dat=(dat(:,1:end-1)+dat(:,2:end))/2;
    end
    ncwriteatt(ncfile,var_name{i},'long_name',name);
    ncwrite(ncfile,var_name{i},dat);    
end

