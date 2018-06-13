% Make a netCDF grid file

% Variables that define the grid
res=25; % Degrees, want 1/100 at some point
folder='~/simulations/SWOT/18-6_grids/';

% Automatically run for both SS and GEBCO bathymetry
for i=1:2

	if i==1
		source='SS';
		bathy_path='/home/kellys/smkelly/software/data_products/SS_topo';
	elseif i==2
		source='GEBCO';
		bathy_path='/home/kellys/smkelly/software/data_products/GEBCO';
	end

	% Define where the grid goes
	fid=[num2str(res),'th_deg_',source,'_bathy.nc'];

	% Grid parameters (probably leave these alone)
	latlims=[-80 66];
	a0=6371e3; % radius of Earth

	% Define output grid
	dx=1/res; 
	dy=dx;
	lon0=dx/2:dx:360; 
	lat0=latlims(1)+dy/2:dy:latlims(2);

	% Obtain bathymetry at full resolution
	if strcmp(source,'SS')
		% Load raw data
		addpath(genpath(bathy_path))
		[bathy.lat,bathy.lon,bathy.H]=satbath(1,[-81 81],[0 360]);

		% Pad grand meridian for interpolation 
		bathy.lon=[bathy.lon(:,end-1)-360 bathy.lon];
		bathy.lat=[bathy.lat(:,end) bathy.lat];
		bathy.H=[bathy.H(:,end-1) bathy.H];
		
		% Interpolate onto desired grid
		H=interp2(bathy.lon,bathy.lat,-bathy.H,lon0',lat0)';
		 
	elseif strcmp(source,'GEBCO')
		% Load raw data
		bathy.H=ncread([bathy_path,'/GEBCO_2014_2D.nc'],'elevation');
		bathy.lon=ncread([bathy_path,'/GEBCO_2014_2D.nc'],'lon');
		bathy.lat=ncread([bathy_path,'/GEBCO_2014_2D.nc'],'lat');
		
		% Need to re-arrange grid from [-180, 180] to [0, 360]
		right=(bathy.lon<0);
		left=(bathy.lon>0);
		bathy.H=[bathy.H(left,:); bathy.H(right,:)];
		bathy.lon=[bathy.lon(left); bathy.lon(right)+360];
		
		% Pad grand meridian for interpolation 
		bathy.lon=[bathy.lon(end)-360; bathy.lon; bathy.lon(1)+360];
		bathy.H=[bathy.H(end,:); bathy.H; bathy.H(1,:)];
		
		% Interpolate onto desired grid
		H=interp2(bathy.lon,bathy.lat',-bathy.H',lon0',lat0)'; 
	end	
	clear bathy
	
	% Remove mountains
	H(H<0)=0;
	
	% Smooth at the grid scale
	H2=[H(end,:); H; H(1,:)];
	H2=AVE2D(H2,3);
	H=H2(2:end-1,:);
	
	% Compute gradients
	%[Nx2 Ny2]=size(H);
	%H2=[H(end,:); H; H(1,:)];
	%dHdx=repmat(1./(a0*cos(lat0/180*pi)),[Nx2 1]).*(H2(3:end,:)-H2(1:end-2,:))/(2*dx/180*pi);
	
	%dHdy=zeros([Nx2 Ny2]);
	%dHdy(:,2:end-1)=1/a0*(H(:,3:end)-H(:,1:end-2))/(2*dy/180*pi);
	
	% Write to NetCDF file
	currentFolder=pwd;
	cd(folder);
	
	mode = netcdf.getConstant('CLOBBER');
	mode = bitor(mode,netcdf.getConstant('NETCDF4'));
	temp=netcdf.create(fid,mode);
	
	% define dimensions
	netcdf.defDim(temp,'x',size(H,1));
	netcdf.defDim(temp,'y',size(H,2));
	netcdf.endDef(temp)
	
	% Close file and return to data directory
	netcdf.close(temp);
	cd(currentFolder);
	
	% Create and write fields
	nccreate([folder,fid],'Nx');
	nccreate([folder,fid],'Ny');
	nccreate([folder,fid],'lon','Dimensions',{'x',size(H,1)});
	nccreate([folder,fid],'lat','Dimensions',{'y',size(H,2)});
	nccreate([folder,fid],'H','Dimensions',{'x',size(H,1),'y',size(H,2)});
	%nccreate([folder,fid],'dHdx','Dimensions',{'x',size(H,1),'y',size(H,2)});
	%nccreate([folder,fid],'dHdy','Dimensions',{'x',size(H,1),'y',size(H,2)});
	
	ncwrite([folder,fid],'Nx',size(H,1));
	ncwrite([folder,fid],'Ny',size(H,2));
	ncwrite([folder,fid],'lon',lon0);
	ncwrite([folder,fid],'lat',lat0);
	ncwrite([folder,fid],'H',H);
	%ncwrite([folder,fid],'dHdx',dHdx);
	%ncwrite([folder,fid],'dHdy',dHdy);

end
