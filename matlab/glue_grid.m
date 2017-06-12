% Write the individual grid NetCDF files to one giant NetCDF file.
clear

MSI=0;

% Add paths that are sometimes lost using qsub
if MSI
	addpath(genpath('/home/kellys/smkelly/software/matlab_libraries/sam_ware'))
	addpath(genpath('/home/kellys/smkelly/software/matlab_libraries/seawater'))
	addpath(genpath('/home/kellys/smkelly/software/data_products/SS_topo'))
	addpath(genpath('/home/kellys/smkelly/software/data_products/OTPS'))
	addpath(genpath('/home/kellys/smkelly/software/data_products/WOA13'))
	folder='~/simulations/SWOT/global_grids/';
else
	folder='~/experiments/SWOT/proc/CSW/global_grids/';
end

fid.bathy='./25th_deg_bathy.nc';
fid.grid='./25th_deg_grid.nc';

N_subgrids=160; % Number of individual grid files
Nm=8;  % Number of modes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load bathy file to get Nx and Ny dimensions
H=ncread(fid.bathy,'H');
lon=ncread(fid.bathy,'lon');
lat=ncread(fid.bathy,'lat');

% Add buffer regions
H=[H(end,:); H; H(1,:)];
H=[H(:,1) H H(:,end)];

dlon=mean(diff(lon));
lon=[lon(1)-dlon; lon; lon(end)+dlon];

dlat=mean(diff(lat));
lat=[lat(1)-dlat; lat; lat(end)+dlat];

[Nx Ny]=size(H);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start the grid file
disp('Glueing grid');

%folder='./';
%fid_grid=[fid,'_grid.nc'];
%name=[folder,fid_grid];

%Create file
mode = netcdf.getConstant('CLOBBER');
mode = bitor(mode,netcdf.getConstant('NETCDF4'));
ncid=netcdf.create(fid.grid,mode);

% define dimensions
netcdf.defDim(ncid,'x',Nx);
netcdf.defDim(ncid,'y',Ny);
netcdf.defDim(ncid,'mode',Nm);
netcdf.endDef(ncid)

% Close file and return to data directory
netcdf.close(ncid);

% Create variables
% Dimensions
nccreate(fid.grid,'Nx');
nccreate(fid.grid,'Ny');
nccreate(fid.grid,'Nm');

% Grid 
nccreate(fid.grid,'dy');
nccreate(fid.grid,'dx');
nccreate(fid.grid,'dz');

nccreate(fid.grid,'lon','Dimensions',{'x',Nx});
nccreate(fid.grid,'lat','Dimensions',{'y',Ny}); 

% Data Variables
nccreate(fid.grid,'H','Dimensions',{'x',Nx,'y',Ny});
nccreate(fid.grid,'f','Dimensions',{'x',Nx,'y',Ny});

nccreate(fid.grid,'c','Dimensions',{'x',Nx,'y',Ny,'mode',Nm});
nccreate(fid.grid,'phi_surf','Dimensions',{'x',Nx,'y',Ny,'mode',Nm});
nccreate(fid.grid,'phi_bott','Dimensions',{'x',Nx,'y',Ny,'mode',Nm});

nccreate(fid.grid,'T_x','Dimensions',{'x',Nx,'y',Ny,'mode',Nm,'mode',Nm});
nccreate(fid.grid,'T_y','Dimensions',{'x',Nx,'y',Ny,'mode',Nm,'mode',Nm});
   
% Write the variables we already know
ncwrite(fid.grid,'Nx',Nx);
ncwrite(fid.grid,'Ny',Ny);
ncwrite(fid.grid,'Nm',Nm);

ncwrite(fid.grid,'lon',lon);
ncwrite(fid.grid,'lat',lat);

ncwrite(fid.grid,'H',H);
ncwrite(fid.grid,'f',repmat(sw_f(lat)',[Nx 1]));

clear H

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write data to file
Ny0=floor((length(lat)-2)/N_subgrids);

for id=1:N_subgrids
    
    % Get x-indices for the slab
    ind_y=(id-1)*Ny0+2;
        
    % Define slab name
   	fid.slab=[fid.grid(1:end-3),num2str(id,'%03d'),'.nc'];
	
    if id==1  % Get grid spacings
        dx=ncread(fid.slab,'dx');
        dy=ncread(fid.slab,'dy');
        dz=ncread(fid.slab,'dz');
        
        Nm2=ncread(fid.slab,'Nm');
        if Nm~=Nm2
           disp('Wrong number of modes') 
        end

        ncwrite(fid.grid,'dx',dx);
        ncwrite(fid.grid,'dy',dy);
        ncwrite(fid.grid,'dz',dz);  
    end
    
    % Read and write variables 
    tmp=ncread(fid.slab,'c');
    tmp=tmp(:,:,:);
    ncwrite(fid.grid,'c',tmp,[2 ind_y(1) 1]);
        
    tmp=ncread(fid.slab,'phi_surf');
    tmp=tmp(:,:,:);
    ncwrite(fid.grid,'phi_surf',tmp,[2 ind_y(1) 1]);
    
    tmp=ncread(fid.slab,'phi_bott');
    tmp=tmp(:,:,:);
    ncwrite(fid.grid,'phi_bott',tmp,[2 ind_y(1) 1]);
    
    tmp=ncread(fid.slab,'T_x');
    tmp=tmp(:,:,:,:);
    ncwrite(fid.grid,'T_x',tmp,[2 ind_y(1) 1 1]);
        
    tmp=ncread(fid.slab,'T_y');
    tmp=tmp(:,:,:,:);
    ncwrite(fid.grid,'T_y',tmp,[2 ind_y(1) 1 1]);
    
    PROGRESS_BAR(id,1:N_subgrids)
end

% Now fill in the borders
tmp=ncread(fid.grid,'c',[2 1 1],[1 Ny Nm]);
ncwrite(fid.grid,'c',tmp,[Nx 1 1]);

tmp=ncread(fid.grid,'c',[Nx-1 1 1],[1 Ny Nm]);
ncwrite(fid.grid,'c',tmp,[1 1 1]);

tmp=ncread(fid.grid,'c',[1 Ny-1 1],[Nx 1 Nm]);
ncwrite(fid.grid,'c',tmp,[1 Ny 1]);

tmp=ncread(fid.grid,'c',[1 2 1],[Nx 1 Nm]);
ncwrite(fid.grid,'c',tmp,[1 1 1]);


tmp=ncread(fid.grid,'phi_bott',[2 1 1],[1 Ny Nm]);
ncwrite(fid.grid,'phi_bott',tmp,[Nx 1 1]);

tmp=ncread(fid.grid,'phi_bott',[Nx-1 1 1],[1 Ny Nm]);
ncwrite(fid.grid,'phi_bott',tmp,[1 1 1]);

tmp=ncread(fid.grid,'phi_bott',[1 Ny-1 1],[Nx 1 Nm]);
ncwrite(fid.grid,'phi_bott',tmp,[1 Ny 1]);

tmp=ncread(fid.grid,'phi_bott',[1 2 1],[Nx 1 Nm]);
ncwrite(fid.grid,'phi_bott',tmp,[1 1 1]);


tmp=ncread(fid.grid,'phi_surf',[2 1 1],[1 Ny Nm]);
ncwrite(fid.grid,'phi_surf',tmp,[Nx 1 1]);

tmp=ncread(fid.grid,'phi_surf',[Nx-1 1 1],[1 Ny Nm]);
ncwrite(fid.grid,'phi_surf',tmp,[1 1 1]);

tmp=ncread(fid.grid,'phi_surf',[1 Ny-1 1],[Nx 1 Nm]);
ncwrite(fid.grid,'phi_surf',tmp,[1 Ny 1]);

tmp=ncread(fid.grid,'phi_surf',[1 2 1],[Nx 1 Nm]);
ncwrite(fid.grid,'phi_surf',tmp,[1 1 1]);


tmp=ncread(fid.grid,'T_x',[2 1 1 1],[1 Ny Nm Nm]);
ncwrite(fid.grid,'T_x',tmp,[Nx 1 1 1]);

tmp=ncread(fid.grid,'T_x',[Nx-1 1 1 1],[1 Ny Nm Nm]);
ncwrite(fid.grid,'T_x',tmp,[1 1 1 1]);

tmp=ncread(fid.grid,'T_x',[1 Ny-1 1 1],[Nx 1 Nm Nm]);
ncwrite(fid.grid,'T_x',tmp,[1 Ny 1 1]);

tmp=ncread(fid.grid,'T_x',[1 2 1 1],[Nx 1 Nm Nm]);
ncwrite(fid.grid,'T_x',tmp,[1 1 1 1]);


tmp=ncread(fid.grid,'T_y',[2 1 1 1],[1 Ny Nm Nm]);
ncwrite(fid.grid,'T_x',tmp,[Nx 1 1 1]);

tmp=ncread(fid.grid,'T_y',[Nx-1 1 1 1],[1 Ny Nm Nm]);
ncwrite(fid.grid,'T_y',tmp,[1 1 1 1]);

tmp=ncread(fid.grid,'T_y',[1 Ny-1 1 1],[Nx 1 Nm Nm]);
ncwrite(fid.grid,'T_y',tmp,[1 Ny 1 1]);

tmp=ncread(fid.grid,'T_y',[1 2 1 1],[Nx 1 Nm Nm]);
ncwrite(fid.grid,'T_y',tmp,[1 1 1 1]);


% Lastly set missing values to 0 and/or smooth Tx and Ty at grid scale
if 1

    for n=1:Nm
        for m=1:Nm
			tmp=ncread(fid.grid,'T_x',[1 1 m n],[Nx Ny 1 1]);  
            tmp(tmp>100)=0;
            %tmp=[tmp(end,:); tmp; tmp(1,:)];
            %tmp=AVE2D_v2(tmp,3);
            %tmp=tmp(2:end-1,:); 
			ncwrite(fid.grid,'T_x',tmp,[1 1 m n]);
			
			tmp=ncread(fid.grid,'T_y',[1 1 m n],[Nx Ny 1 1]);  
            tmp(tmp>100)=0;
            %tmp=[tmp(end,:); tmp; tmp(1,:)];
            %tmp=AVE2D_v2(tmp,3);
            %tmp=tmp(2:end-1,:); 
			ncwrite(fid.grid,'T_y',tmp,[1 1 m n]);
        end
    end
     	
end
