function make_grid(id_node,N_subgrids,N_processors,N_threads)
% Matlab function to compute c_n and T_mn
%
% The script presently smooths WOA N^2 over 5 degrees and fill gaps with "fillmissing.m"
%
% Explanation of input parameters:
%
%id_node       Start ID under this instance, starts at 0
%N_subgrids    Total number of subgrids, typically 36 or 160
%N_processors  Number of processors under this instance, typically 8 on Itasca
%N_threads     Number of threads per processor, typically 1

% First slab to work on
id_start=1+N_processors*id_node;

% Add path to raw bathymetry data
addpath(genpath('/home/kellys/smkelly/software/data_products/SS_topo'))

% Define where the grid goes
folder='~/simulations/SWOT/18-5_global_grids/';

% Variables that define the grid
dx=1/25; % Degrees, want 1/100 at some point
fid='25th_deg';
make_bathy=0;

% Calculation parameters
fid_strat=['./strat_HYCOM.nc'];
Nm=8; % Want Nm=8 eventually
Nm0=128; % Number of structure functions to solve (want Nm0=128)
dz=1; % Want dz=1 eventually

% Grid parameters (probably leave these alone)
H_min=16; % want H_min=16 eventually
H_max=6000;
latlims=[-80 66];
a0=6371e3; % radius of Earth

% Define output grid
dy=dx;
lon0=dx/2:dx:360; 
lat0=latlims(1)+dy/2:dy:latlims(2);
z=(dz/2):dz:H_max;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create bathymetry if necessary 
fid_bathy=[fid,'_bathy.nc'];
name_bathy=[folder,fid_bathy];

if make_bathy
    
    % Obtain bathymetry at full (1/60 degree, i.e., 1 minute in longitude) resolution
    [bathy.lat,bathy.lon,bathy.H]=satbath(1,[-81 81],[0 360]);
    bathy.lon=[bathy.lon(:,end-1)-360 bathy.lon];
    bathy.lat=[bathy.lat(:,end) bathy.lat];
    bathy.H=[bathy.H(:,end-1) bathy.H];
    
    % Interpolate onto desired grid
    H=interp2(bathy.lon,bathy.lat,-bathy.H,lon0',lat0)'; 
    clear bathy

    % Remove mountains
    H(H<0)=0;
    
    % Smooth at the grid scale
    H2=[H(end,:); H; H(1,:)];
    H2=AVE2D(H2,3);
    H=H2(2:end-1,:);
    
    % Compute gradients
    [Nx2 Ny2]=size(H);
    H2=[H(end,:); H; H(1,:)];
    dHdx=repmat(1./(a0*cos(lat0/180*pi)),[Nx2 1]).*(H2(3:end,:)-H2(1:end-2,:))/(2*dx/180*pi);

    dHdy=zeros([Nx2 Ny2]);
    dHdy(:,2:end-1)=1/a0*(H(:,3:end)-H(:,1:end-2))/(2*dy/180*pi);
    
    % Write to NetCDF file
    currentFolder=pwd;
    cd(folder);
    
    mode = netcdf.getConstant('CLOBBER');
    mode = bitor(mode,netcdf.getConstant('NETCDF4'));
    temp=netcdf.create(fid_bathy,mode);
    
    % define dimensions
    netcdf.defDim(temp,'x',size(H,1));
    netcdf.defDim(temp,'y',size(H,2));
    netcdf.endDef(temp)

	% Close file and return to data directory
    netcdf.close(temp);
    cd(currentFolder);
       
   	% Create and write fields
    nccreate(name_bathy,'Nx');
    nccreate(name_bathy,'Ny');
    nccreate(name_bathy,'lon','Dimensions',{'x',size(H,1)});
    nccreate(name_bathy,'lat','Dimensions',{'y',size(H,2)});
    nccreate(name_bathy,'H','Dimensions',{'x',size(H,1),'y',size(H,2)});
    nccreate(name_bathy,'dHdx','Dimensions',{'x',size(H,1),'y',size(H,2)});
    nccreate(name_bathy,'dHdy','Dimensions',{'x',size(H,1),'y',size(H,2)});

    ncwrite(name_bathy,'Nx',size(H,1));
    ncwrite(name_bathy,'Ny',size(H,2));
    ncwrite(name_bathy,'lon',lon0);
    ncwrite(name_bathy,'lat',lat0);
    ncwrite(name_bathy,'H',H);
    ncwrite(name_bathy,'dHdx',dHdx);
    ncwrite(name_bathy,'dHdy',dHdy);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine size of the subgrids
Nx=length(lon0);
Ny0=floor(length(lat0)/N_subgrids);
Nz=length(z);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start parallel enviornment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while isempty(gcp('nocreate'))
	try
		parpool(N_processors); % Start the workers
	end
end

% Get longitude for each processor
lon=[lon0(1)-dx lon0 lon0(end)+dx];


% Start processor loop
parfor id_subindex=1:N_processors
    
    % Set the maximum number of threads
    warning('off','MATLAB:maxNumCompThreads:Deprecated')
    maxNumCompThreads(N_threads);

    % Define sub grid size and indicies
    id=(id_start-1)+id_subindex; 
        
    % Find indicies to compute (split by latitude)
    if id==N_subgrids
        ind_y=(id-1)*Ny0:length(lat0);
    else
        ind_y=(1:Ny0)+(id-1)*Ny0;
    end
    Ny=length(ind_y);
    
    % Load bathymetry
    H=ncread(name_bathy,'H',[1 ind_y(1)],[Nx Ny]);
    
    % Now add boundary points
    if ind_y(1)==1
            H_N=ncread(name_bathy,'H',[1 ind_y(1)],[Nx 1]);
            H_S=ncread(name_bathy,'H',[1 ind_y(end)+1],[Nx 1]);
    elseif ind_y(end)==length(lat0)
            H_N=ncread(name_bathy,'H',[1 ind_y(1)-1],[Nx 1]);
            H_S=ncread(name_bathy,'H',[1 ind_y(end)],[Nx 1]);
    else
            H_N=ncread(name_bathy,'H',[1 ind_y(1)-1],[Nx 1]);
            H_S=ncread(name_bathy,'H',[1 ind_y(end)+1],[Nx 1]);
    end
    H=[H_N H H_S];
    H=[H(end,:); H; H(1,:)];
    
    % Find lat for local slab
    lat=lat0(ind_y);
    lat=[lat(1)-dx lat lat(end)+dx];

    % Eliminate deep and shallow regions
    H(H>H_max)=H_max;
    H(H<H_min)=0;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % start grid NetCDF file
    fid_grid=[fid,'_grid',num2str(id,'%03d'),'.nc'];
    name=[folder,fid_grid];
    
    % figure out where to start
    if  exist(name, 'file') ~= 2
        
        disp('Making new grid file');
        ind_x_start=1;
        
        % Create file
        currentFolder=pwd;
        cd(folder);
        
        mode = netcdf.getConstant('CLOBBER');
        mode = bitor(mode,netcdf.getConstant('NETCDF4'));
        temp=netcdf.create(fid_grid,mode);
                
        % define dimensions
        netcdf.defDim(temp,'x',Nx);
        netcdf.defDim(temp,'y',Ny);
        netcdf.defDim(temp,'mode',Nm);
        netcdf.endDef(temp)

		% Close file and return to data directory
        netcdf.close(temp);        
        cd(currentFolder);
        
        % Write data
        % Dimmensions
        nccreate(name,'Nx');
        nccreate(name,'Ny');
        nccreate(name,'Nm');
        
        % Grid
        nccreate(name,'dy');
        nccreate(name,'dx');
        nccreate(name,'dz');
        
        nccreate(name,'lon','Dimensions',{'x',Nx});
        nccreate(name,'lat','Dimensions',{'y',Ny});
        
        % Variables
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
        
        ncwrite(name,'dx',dx);
        ncwrite(name,'dy',dy);
        ncwrite(name,'dz',dz);
        
        ncwrite(name,'lon',lon(2:end-1));
        ncwrite(name,'lat',lat(2:end-1));
        
        ncwrite(name,'H',H(2:end-1,2:end-1));
        %ncwrite(name,'f',repmat(sw_f(lat(2:end-1)),[Nx 1]));
        ncwrite(name,'f',repmat(2*(7.292e-5)*sin(lat(2:end-1)/180*pi),[Nx 1]));
             
    else
        % Figure out where to start
        tmp=ncread(name,'c');
        ind_x_start=find(mean(tmp(:,:,1),2)<1000,1,'last');
        tmp=[];
        
        if isempty(ind_x_start); ind_x_start=1; end;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Now compute c and Tx
    if ind_x_start<Nx
        disp(['Starting from index ',num2str(ind_x_start),' of ',num2str(Nx)]);
           
        strat.N2=ncread(fid_strat,'N2');
        strat.z=ncread(fid_strat,'depth');
        strat.lon=ncread(fid_strat,'lon')';
        strat.lat=ncread(fid_strat,'lat')';

        
        % Smooth the stratification over three degrees
        z0=strat.z;
        N20=NaN([Nx+2 Ny+2 length(z0)]);
        lon_tmp=[strat.lon(end-1)-360 strat.lon(end)-360 strat.lon strat.lon(1)+360 strat.lon(2)+360];
        for i=1:length(z0)

            % Extend N2 slightly
            N2_tmp=[strat.N2(end-1,:,i); strat.N2(end,:,i); strat.N2(:,:,i); strat.N2(1,:,i); strat.N2(2,:,i)];
            N2_tmp(N2_tmp<0 | N2_tmp>1e-2)=NaN;
                
            % Smooth WOA
            if strcmp(fid_strat,'./strat_WOA.nc')
                N2_tmp=AVE2D_v2(N2_tmp,5);
            end
            
            % Take nearest neighbor (at same latitude) for missing data
            bad=(N2_tmp<0 | isnan(N2_tmp) | N2_tmp>1e-2);
            N2_tmp(bad)=NaN;
            N2_tmp=fillmissing(N2_tmp,'nearest');

            % Interpolate onto grid
            N20(:,:,i)=interp2(lon_tmp',strat.lat,N2_tmp.',lon',lat).';
            
            % Take nearest neighbor (at same logintude) for missing data
            % (typically near Antarctica because HYCOM doesn't go that far south)
            N20(:,:,i)=fillmissing(N20(:,:,i)','nearest')';
        end
        strat=[];% This is really a clear command
        
        % Initiate temporary arrays
        phi_W=NaN([Nz Nm Ny+2]);
        phi_C=NaN([Nz Nm Ny+2]);
        phi_E=NaN([Nz Nm Ny+2]);
        for i=ind_x_start:Nx+2
            
            % Pre-allocate output arrays
            c=zeros([1 Ny+2 Nm]);
            phi_surf=zeros([1 Ny+2 Nm]);
            phi_bott=zeros([1 Ny+2 Nm]);
            Tx=zeros([1 Ny+2 Nm Nm]);
            Ty=zeros([1 Ny+2 Nm Nm]);
            
            for j=1:Ny+2
                
                % Compute c, phi(0), and phi(-H) at current longitude
                if H(i,j)>H_min
                    
                    % Identify the proper number of modes to find
                    Nm0_H=min([Nm0 floor(H(i,j)/(2*dz))+1]);
                    Nm_H=min([Nm floor(H(i,j)/(16*dz))]); % This is the old condition
                    %Nm_H=min([Nm floor(H(i,j)/(8*dz))]);
                    
                    N2=squeeze(N20(i,j,:));
                    good=isfinite(N2);
                    
                    if sum(good)>=2 % Note: strat.z(4)=17.5 for WOA or strat.z(2)=10 m for HYCOM                   
                        % Interpolate coarse WOA N2 onto fine depth grid
                        N2=interp1(z0(good),N2(good),z).';
                        N2(N2<1e-8)=1e-8;
                        ind_z=find(isfinite(N2),1,'first');
                        N2(1:ind_z)=N2(ind_z);
                        ind_z=find(isfinite(N2),1,'last');
                        N2(ind_z:end)=N2(ind_z);
                        
                        % Compute the modes
                        ind_z=round(H(i,j)/dz);
                        [PHI,c(1,j,1:Nm_H)]=MODES_FAST(dz,N2(1:ind_z),Nm_H,Nm0_H);
                        
                        % Write to array
                        phi_surf(1,j,1:Nm_H)=PHI(1,:)';
                        phi_bott(1,j,1:Nm_H)=PHI(end,:)';
                        phi_E(1:ind_z,1:Nm_H,j)=PHI;
                    end
                end
                
                % Compute Tx and Ty at previous longitude and latitude
                if (ind_x_start+2)<=i && 3<=j && H(i-1,j-1)>0
                    
                    % Find correct vertical integration limit
                    ind_z=find(isfinite(phi_E(:,1,j-1)+phi_C(:,1,j-1)+phi_W(:,1,j-1)),1,'last');
                    if ~isempty(ind_z)
                        
                        % Mode at location of interest
                        PHI=phi_C(ind_z,:,j-1);
                        
                        % Gradient in spherical coordinates
                        dPHIdx=1/(a0*cos(lat(j-1)/180*pi))*(phi_E(ind_z,:,j-1)-phi_W(ind_z,:,j-1))/(2*dx/180*pi);
                        
                        % Integral for coupling coefficient
                        Tx(1,j-1,:,:)=dPHIdx'*PHI*dz;
                    end
                    
                    % Find correct vertical integration limit
                    ind_z=find(isfinite(phi_C(:,1,j-2)+phi_C(:,1,j-1)+phi_C(:,1,j)),1,'last');
                    if ~isempty(ind_z)
                        
                        % Mode at location of interest
                        PHI=phi_C(ind_z,:,j-1);
                        
                        % Gradient in spherical coordinates
                        dPHIdy=1/a0*(phi_C(ind_z,:,j)-phi_C(ind_z,:,j-2))/(2*dy/180*pi);
                        
                        % Integral for coupling coefficient
                        Ty(1,j-1,:,:)=dPHIdy'*PHI*dz;
                    end
                end
            end
            
            % Write output to file
            if (ind_x_start+1)<=i && i<=Nx+1 % Must be on the second loop to write c
                ncwrite(name,'c',c(:,2:end-1,:),[i-1 1 1]);
                ncwrite(name,'phi_surf',phi_surf(:,2:end-1,:),[i-1 1 1]);
                ncwrite(name,'phi_bott',phi_bott(:,2:end-1,:),[i-1 1 1]);
            end
            
            if (ind_x_start+2)<=i % Must be on the third loop to write Tx
                Tx(isnan(Tx))=0;
                Ty(isnan(Ty))=0;
                ncwrite(name,'T_x',Tx(:,2:end-1,:,:),[i-2 1 1 1]);
                ncwrite(name,'T_y',Ty(:,2:end-1,:,:),[i-2 1 1 1]);
            end
            
            % Update phi arrays for next longitude
            phi_W=phi_C;
            phi_C=phi_E;
            phi_E=NaN([Nz Nm Ny+2]);
            
        end
        disp('Clean finish');
        
    else
        disp('Slab was already complete');
    end
    
end % end SPMD

% Close parallel pool
delete(gcp('nocreate'));

return;

