% Write the individual grid NetCDF files to one giant NetCDF file.
clear

MSI=1;

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

fid='50th_deg';
N_subgrids=160; % Number of individual grid files
Nm=8;  % Number of modes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load bathy file to get Nx and Ny dimensions
folder='./';
fid_bathy=[fid,'_bathy.nc'];
name=[folder,fid_bathy];

H=ncread(name,'H');
lon=ncread(name,'lon');
lat=ncread(name,'lat');

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

folder='./';
fid_grid=[fid,'_grid.nc'];
name=[folder,fid_grid];

%Create file
mode = netcdf.getConstant('CLOBBER');
mode = bitor(mode,netcdf.getConstant('NETCDF4'));
ncid=netcdf.create(name,mode);

% define dimensions
netcdf.defDim(ncid,'x',Nx);
netcdf.defDim(ncid,'y',Ny);
netcdf.defDim(ncid,'mode',Nm);
netcdf.endDef(ncid)

% Close file and return to data directory
netcdf.close(ncid);

% Create variables
% Dimensions
nccreate(name,'Nx');
nccreate(name,'Ny');
nccreate(name,'Nm');

% Grid 
nccreate(name,'dy');
nccreate(name,'dx');
nccreate(name,'dz');

nccreate(name,'lon','Dimensions',{'x',Nx});
nccreate(name,'lat','Dimensions',{'y',Ny}); 

% Data Variables
nccreate(name,'H','Dimensions',{'x',Nx,'y',Ny});
nccreate(name,'f','Dimensions',{'x',Nx,'y',Ny});

nccreate(name,'c','Dimensions',{'x',Nx,'y',Ny,'mode',Nm});
nccreate(name,'phi_surf','Dimensions',{'x',Nx,'y',Ny,'mode',Nm});
nccreate(name,'phi_bott','Dimensions',{'x',Nx,'y',Ny,'mode',Nm});

nccreate(name,'T_x','Dimensions',{'x',Nx,'y',Ny,'mode',Nm,'mode',Nm});
nccreate(name,'T_y','Dimensions',{'x',Nx,'y',Ny,'mode',Nm,'mode',Nm});
   
% Write the variables we already know
ncwrite(name,'Nx',Nx);
ncwrite(name,'Ny',Ny);
ncwrite(name,'Nm',Nm);

ncwrite(name,'lon',lon);
ncwrite(name,'lat',lat);

ncwrite(name,'H',H);
ncwrite(name,'f',repmat(sw_f(lat)',[Nx 1]));

clear H

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write data to file
Ny0=floor((length(lat)-2)/N_subgrids);

for id=1:N_subgrids
    
    % Get x-indices for the slab
    ind_y=(id-1)*Ny0+2;
        
    % Define slab name
    folder='./slabs/';
    fid_slab=[fid,'_grid',num2str(id,'%03d'),'.nc'];
    name_slab=[folder,fid_slab];

    if id==1  % Get grid spacings
        dx=ncread(name_slab,'dx');
        dy=ncread(name_slab,'dy');
        dz=ncread(name_slab,'dz');
        
        Nm2=ncread(name_slab,'Nm');
        if Nm~=Nm2
           disp('Wrong number of modes') 
        end

        ncwrite(name,'dx',dx);
        ncwrite(name,'dy',dy);
        ncwrite(name,'dz',dz);  
    end
    
    % Read and write variables 
    tmp=ncread(name_slab,'c');
    tmp=tmp(:,:,:);
    ncwrite(name,'c',tmp,[2 ind_y(1) 1]);
        
    tmp=ncread(name_slab,'phi_surf');
    tmp=tmp(:,:,:);
    ncwrite(name,'phi_surf',tmp,[2 ind_y(1) 1]);
    
    tmp=ncread(name_slab,'phi_bott');
    tmp=tmp(:,:,:);
    ncwrite(name,'phi_bott',tmp,[2 ind_y(1) 1]);
    
    tmp=ncread(name_slab,'T_x');
    tmp=tmp(:,:,:,:);
    ncwrite(name,'T_x',tmp,[2 ind_y(1) 1 1]);
        
    tmp=ncread(name_slab,'T_y');
    tmp=tmp(:,:,:,:);
    ncwrite(name,'T_y',tmp,[2 ind_y(1) 1 1]);
    
    PROGRESS_BAR(id,1:N_subgrids)
end

% Now fill in the borders
tmp=ncread(name,'c',[2 1 1],[1 Ny Nm]);
ncwrite(name,'c',tmp,[Nx 1 1]);

tmp=ncread(name,'c',[Nx-1 1 1],[1 Ny Nm]);
ncwrite(name,'c',tmp,[1 1 1]);

tmp=ncread(name,'c',[1 Ny-1 1],[Nx 1 Nm]);
ncwrite(name,'c',tmp,[1 Ny 1]);

tmp=ncread(name,'c',[1 2 1],[Nx 1 Nm]);
ncwrite(name,'c',tmp,[1 1 1]);


tmp=ncread(name,'phi_bott',[2 1 1],[1 Ny Nm]);
ncwrite(name,'phi_bott',tmp,[Nx 1 1]);

tmp=ncread(name,'phi_bott',[Nx-1 1 1],[1 Ny Nm]);
ncwrite(name,'phi_bott',tmp,[1 1 1]);

tmp=ncread(name,'phi_bott',[1 Ny-1 1],[Nx 1 Nm]);
ncwrite(name,'phi_bott',tmp,[1 Ny 1]);

tmp=ncread(name,'phi_bott',[1 2 1],[Nx 1 Nm]);
ncwrite(name,'phi_bott',tmp,[1 1 1]);


tmp=ncread(name,'phi_surf',[2 1 1],[1 Ny Nm]);
ncwrite(name,'phi_surf',tmp,[Nx 1 1]);

tmp=ncread(name,'phi_surf',[Nx-1 1 1],[1 Ny Nm]);
ncwrite(name,'phi_surf',tmp,[1 1 1]);

tmp=ncread(name,'phi_surf',[1 Ny-1 1],[Nx 1 Nm]);
ncwrite(name,'phi_surf',tmp,[1 Ny 1]);

tmp=ncread(name,'phi_surf',[1 2 1],[Nx 1 Nm]);
ncwrite(name,'phi_surf',tmp,[1 1 1]);


tmp=ncread(name,'T_x',[2 1 1 1],[1 Ny Nm Nm]);
ncwrite(name,'T_x',tmp,[Nx 1 1 1]);

tmp=ncread(name,'T_x',[Nx-1 1 1 1],[1 Ny Nm Nm]);
ncwrite(name,'T_x',tmp,[1 1 1 1]);

tmp=ncread(name,'T_x',[1 Ny-1 1 1],[Nx 1 Nm Nm]);
ncwrite(name,'T_x',tmp,[1 Ny 1 1]);

tmp=ncread(name,'T_x',[1 2 1 1],[Nx 1 Nm Nm]);
ncwrite(name,'T_x',tmp,[1 1 1 1]);


tmp=ncread(name,'T_y',[2 1 1 1],[1 Ny Nm Nm]);
ncwrite(name,'T_x',tmp,[Nx 1 1 1]);

tmp=ncread(name,'T_y',[Nx-1 1 1 1],[1 Ny Nm Nm]);
ncwrite(name,'T_y',tmp,[1 1 1 1]);

tmp=ncread(name,'T_y',[1 Ny-1 1 1],[Nx 1 Nm Nm]);
ncwrite(name,'T_y',tmp,[1 Ny 1 1]);

tmp=ncread(name,'T_y',[1 2 1 1],[Nx 1 Nm Nm]);
ncwrite(name,'T_y',tmp,[1 1 1 1]);


% Lastly set missing values to 0 and/or smooth Tx and Ty at grid scale
if 1

    for n=1:Nm
        for m=1:Nm
			tmp=ncread(name,'T_x',[1 1 m n],[Nx Ny 1 1]);  
            tmp(tmp>100)=0;
            %tmp=[tmp(end,:); tmp; tmp(1,:)];
            %tmp=AVE2D_v2(tmp,3);
            %tmp=tmp(2:end-1,:); 
			ncwrite(name,'T_x',tmp,[1 1 m n]);
			
			tmp=ncread(name,'T_y',[1 1 m n],[Nx Ny 1 1]);  
            tmp(tmp>100)=0;
            %tmp=[tmp(end,:); tmp; tmp(1,:)];
            %tmp=AVE2D_v2(tmp,3);
            %tmp=tmp(2:end-1,:); 
			ncwrite(name,'T_y',tmp,[1 1 m n]);
        end
    end
     	
end
