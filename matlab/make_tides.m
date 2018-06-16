% Coupled Shallow Water model (CSW)
clear

ii=complex(0,1);
warning off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set input parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
res=10;

tide_list={'TPXO','GOT','HAMTIDE','FES'};
bathy_list={'SS','GEBCO'};

for i=1:length(tide_list)
	for j=1:length(bathy_list)
		% Choose the inputs
		tide.name=tide_list{i}; % can be TPXO, GOT, HAMTIDE, or FES
		tide.bathy=bathy_list{j}; % can be SS or GEBCO
		tide.folder='../../../tides/';
		
		% Output files
		fid.bathy=['../../../18-6_grids/',num2str(res),'th_deg_',tide.bathy,'_bathy.nc'];
		fid.tides=['../../',tide.name,'_',tide.bathy,'.nc'];
		
		Nc=1;       % Number of tidal constituents 
		H_min=16;	% Set minimum depth for tides
		thresh=1;   % Maximum tidal velocity 
		
		
		%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%
		%                    Technical code below here
		%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		disp('Loading grid');
		
		Nx=ncread(fid.bathy,'Nx');
		Ny=ncread(fid.bathy,'Ny');
		
		H=ncread(fid.bathy,'H'); 
		lon=ncread(fid.bathy,'lon'); 
		lat=ncread(fid.bathy,'lat'); 
		
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Shallow water mask
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		disp('Masking tides in shallow regions');
		  
		% Mask out shallow regions
		H(H<H_min)=0;
		
		% Mask out the Caspian Sea (seems to create a CFL instability)
		bad.lon=45<lon & lon<60;
		bad.lat=35<lat & lat<50;
		H(bad.lon,bad.lat)=0;
		
		%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Surface-tide forcing
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		disp('Loading tides');
		
		tide.lon=ncread([tide.folder,'m2UV_',tide.name,'.nc'],'lon');
		tide.lat=ncread([tide.folder,'m2UV_',tide.name,'.nc'],'lat');
		
		UAMP=ncread([tide.folder,'m2UV_',tide.name,'.nc'],'UAMP');
		VAMP=ncread([tide.folder,'m2UV_',tide.name,'.nc'],'VAMP');
		if strcmp(tide.name,'FES') | strcmp(tide.name,'HAMTIDE')
		    H0=interp2([lon(end)-360; lon; lon(1)+360],lat',[H(end,:); H; H(1,:)]',tide.lon,tide.lat')';
		    UAMP=H0.*UAMP;
		    VAMP=H0.*VAMP;
		end
		UPHA=ncread([tide.folder,'m2UV_',tide.name,'.nc'],'UPHA');
		VPHA=ncread([tide.folder,'m2UV_',tide.name,'.nc'],'VPHA');
		
		U=UAMP.*exp(ii*UPHA/180*pi);
		U(isnan(U))=0;
		
		clear UAMP UPHA
		U=interp2([tide.lon(end)-360; tide.lon; tide.lon(1)+360],tide.lat',[U(end,:); U; U(1,:)].',lon,lat').';
		
		V=VAMP.*exp(ii*VPHA/180*pi);
		V(isnan(V))=0;
		clear VAMP VPHA
		V=interp2([tide.lon(end)-360; tide.lon; tide.lon(1)+360],tide.lat',[V(end,:); V; V(1,:)].',lon,lat').';
		 
		% Convert transport to velocity 
		for k=1:Nc
			U(:,:,k)=U(:,:,k)./H;
			V(:,:,k)=V(:,:,k)./H;   
		end
		
		% Remove bad data
		U(isinf(U) | isnan(U))=0+ii*0;
		V(isinf(V) | isnan(V))=0+ii*0;
		
		% Limit surface tides to 1 m/s
		fast=thresh<abs(U);
		U(fast)=U(fast).*thresh./abs(U(fast));
		fast=thresh<abs(V);
		V(fast)=V(fast).*thresh./abs(V(fast));
		
		for k=1:Nc
		    tmp=U(:,:,k);
		    tmp(H<H_min)=0;
		    U(:,:,k)=tmp;
		    
		    tmp=V(:,:,k);
		    tmp(H<H_min)=0;
		    V(:,:,k)=tmp;
		end
		
		% Round the tidal frequency and write to file
		omega0=1.405189e-04;
		period=(2*pi/omega0)/3600; % in hours
		period=round(period*100)/100; % Round to second decimal place (i.e., 12.42 for the M2 tide)
		omega=(2*pi)/(period*3600);
		
		% Write to file
		ncwrite(fid.tides,'omega',omega);
		ncwrite(fid.tides,'Ur',real(U));
		ncwrite(fid.tides,'Ui',imag(U));
		ncwrite(fid.tides,'Vr',real(V));
		ncwrite(fid.tides,'Vi',imag(V));
	end
end

