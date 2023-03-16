% Make a grid file for global simulations
%
% This script requires access to:
% (1) WOA23 temperature and salinity annual average NetCDF files OR the pre-made strat.mat 
% (2) Smith & Sandwell topography: topo_24.1.nc
% (3) Argo mixed layer climatology: Argo_mixedlayers_monthlyclim_04142022.nc

% The code computes the vertical modes on the native bathymetry grid and saves
% (1) Bottom depth
% (2) Eigenspeeds
% (3) phi_surf (the horizontal-velocity/pressure vertical mode's value at the top)
% (4) phi_bott (the horizontal-velocity/pressure vertical mode's value at the bottom)
%  
% Note 1: phi_surf is necessary for wind forcing. It is actually an integral of the phi and the surface stress profile 
% Note 2: phi_bott is necessary for tidal forcing and topographic coupling see Zaron et al. (2022) JPO
% Note 3: The idea is that all of these values can be smoothed and interpolated onto a grid of any resolution. 
%
% smkelly@d.umn.edu 1/10/2023

clear

% User-defined inputs
season=4; % 0 for annual mean, 1 for NH winter, 2 for NH summer, 3 for June, 4 for October 
Nm=16;
Nm0=128; % 128 seems good
dz=1;
H_max=6000;
z=(dz/2):dz:H_max;

if season==0
    fid.grid='SS_WOA_grid.nc';
elseif season==1
    fid.grid='SS_WOA_winter_grid.nc';
elseif season==2
    fid.grid='SS_WOA_summer_grid.nc';
elseif season==3
    fid.grid='SS_WOA_JUN_grid.nc';
elseif season==4
    fid.grid='SS_WOA_OCT_grid.nc';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load World Ocean Atlas (1/4 degree)
if 1
    if season==0
        load('strat.mat');
    elseif season==1
        load('strat_winter.mat');
    elseif season==2
        load('strat_summer.mat');
    elseif season==3
        load('strat_JUN.mat');
    elseif season==4
        load('strat_OCT.mat');
    end
    
else
    if season==0
        strat.lon=ncread('woa23_decav91C0_t00_04.nc','lon');
        strat.lat=ncread('woa23_decav91C0_t00_04.nc','lat');
        strat.z=ncread('woa23_decav91C0_t00_04.nc','depth');
        strat.T=ncread('woa23_decav91C0_t00_04.nc','t_an');
        strat.S=ncread('woa23_decav91C0_s00_04.nc','s_an');
    elseif season==1
        strat.lon=ncread('woa23_decav91C0_t13_04.nc','lon');
        strat.lat=ncread('woa23_decav91C0_t13_04.nc','lat');
        strat.z=ncread('woa23_decav91C0_t13_04.nc','depth');
        strat.T=ncread('woa23_decav91C0_t13_04.nc','t_an');
        strat.S=ncread('woa23_decav91C0_s13_04.nc','s_an');
    elseif season==2
        strat.lon=ncread('woa23_decav91C0_t15_04.nc','lon');
        strat.lat=ncread('woa23_decav91C0_t15_04.nc','lat');
        strat.z=ncread('woa23_decav91C0_t15_04.nc','depth');
        strat.T=ncread('woa23_decav91C0_t15_04.nc','t_an');
        strat.S=ncread('woa23_decav91C0_s15_04.nc','s_an');
    elseif season==3
        strat.lon=ncread('woa23_decav91C0_t06_04.nc','lon');
        strat.lat=ncread('woa23_decav91C0_t06_04.nc','lat');
        strat.z=ncread('woa23_decav91C0_t06_04.nc','depth');
        strat.T=ncread('woa23_decav91C0_t06_04.nc','t_an');
        strat.S=ncread('woa23_decav91C0_s06_04.nc','s_an');
    elseif season==4
        strat.lon=ncread('woa23_decav91C0_t10_04.nc','lon');
        strat.lat=ncread('woa23_decav91C0_t10_04.nc','lat');
        strat.z=ncread('woa23_decav91C0_t10_04.nc','depth');
        strat.T=ncread('woa23_decav91C0_t10_04.nc','t_an');
        strat.S=ncread('woa23_decav91C0_s10_04.nc','s_an');
    end
    [Nx Ny Nz]=size(strat.T);
    
    disp('Calculating buoyancy frequency')
    strat.z_ave=(strat.z(1:end-1)+strat.z(2:end))/2;
    strat.N2=NaN([Nx Ny Nz-1]);
    for j=1:Ny
        p=sw_pres(strat.z,strat.lat(j));
        for i=1:Nx
            strat.N2(i,j,:)=sw_bfrq(squeeze(strat.S(i,j,:)),squeeze(strat.T(i,j,:)),p,strat.lat(j));
        end
        PROGRESS_BAR(j,1:Ny)
    end
    
    % Now extrapolate and fill gaps using a 2D moving average
    tmp=strat.N2;
    for i=1:Nx
        indx=i-4:i+4;
        indx(indx<1)=indx(indx<1)+Nx;
        indx(Nx<indx)=indx(Nx<indx)-Nx;
            
        for j=1:Ny            
            indy=j-4:j+4;
            indy=indy(0<indy & indy<Ny+1);
            
            for k=1:Nz-1
                if isnan(strat.N2(i,j,k))
                    tmp2=strat.N2(indx,indy,k);
                    tmp(i,j,k)=nanmean(tmp2(:));
                end
            end
        end
        PROGRESS_BAR(i,1:Nx)
    end
    strat.N2=tmp;
    
    % Pad at the dateline
    strat.lon=[strat.lon(end)-360; strat.lon; strat.lon(1)+360];
    strat.N2=[strat.N2(end,:,:); strat.N2; strat.N2(1,:,:)];
    
    % Save for quicker calculations later
    strat.z=strat.z_ave;
    strat=rmfield(strat,'z_ave');
    strat=rmfield(strat,'T');
    strat=rmfield(strat,'S');
    
    if season==0
        save('strat','strat');
    elseif season==1
        save('strat_winter','strat');
    elseif season==2
        save('strat_summer','strat');
    elseif season==3
        save('strat_JUN','strat');
    elseif season==4
        save('strat_OCT','strat');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load mixed layer depth
mix.lon=ncread('Argo_mixedlayers_monthlyclim_04142022.nc','lon');
mix.lat=ncread('Argo_mixedlayers_monthlyclim_04142022.nc','lat');
mix.H=ncread('Argo_mixedlayers_monthlyclim_04142022.nc','mld_da_mean');

% Now extrapolate and fill gaps using a 2D moving average
[~,Nx,Ny]=size(mix.H);
tmp=mix.H;
for i=1:Nx
    indx=i-2:i+2;
    indx(indx<1)=indx(indx<1)+Nx;
    indx(indx>Nx)=indx(indx>Nx)-Nx;
    
    for j=1:Ny
        indy=j-2:j+2;
        indy=indy(0<indy & indy<Ny+1);
        
        for k=1:12
            if isnan(tmp(k,i,j))
                tmp2=mix.H(k,indx,indy);
                tmp(k,i,j)=nanmean(tmp2(:));
            end
        end
    end
end
mix.H=permute(tmp,[2 3 1]);

if season==0
    mix.H0=nanmean(mix.H,3);
    mix.zT=mix.H0;
elseif season==1
    mix.H0=nanmean(mix.H(:,:,1:3),3);
    mix.zT=mix.H0;
elseif season==2
    mix.H0=nanmean(mix.H(:,:,7:9),3);
    mix.zT=mix.H0;
elseif season==3
    mix.H0=mix.H(:,:,6);
     mix.zT=mix.H0;
    mix.H0(100:220,90:180)=10;
    mix.zT(100:220,90:180)=40;
elseif season==4
    mix.H0=mix.H(:,:,10);
    mix.zT=mix.H0;
    mix.H0([1:82 287:360],90:180)=35;
    mix.zT([1:82 287:360],90:180)=70;
end

% Smooth slightly
tmp=[mix.H0(end,:); mix.H0; mix.H0(1,:)];
tmp=AVE2D(tmp,3);
mix.H0=tmp(2:end-1,:);

tmp=[mix.zT(end,:); mix.zT; mix.zT(1,:)];
tmp=AVE2D(tmp,3);
mix.zT=tmp(2:end-1,:);

% Pad at the dateline
mix.lon=[mix.lon(end)-360; mix.lon; mix.lon(1)+360];
mix.H0=[mix.H0(end,:,:); mix.H0; mix.H0(1,:,:)];
mix.zT=[mix.zT(end,:,:); mix.zT; mix.zT(1,:,:)];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load SS topography (1/60 degree)
lon=ncread('topo_24.1.nc','lon');
lat=ncread('topo_24.1.nc','lat');
H=-ncread('topo_24.1.nc','z');
H(H<0)=0;

% Interpolate mixedlayer depth
H_mix=interp2(mix.lon,mix.lat',mix.H0',lon,lat')';
H_mix(isnan(H_mix) | H_mix<10)=10;

zT=interp2(mix.lon,mix.lat',mix.zT',lon,lat')';
zT(isnan(zT) | zT<10)=10;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start a NetCDF file
disp('Writing grid');
[Nx Ny]=size(H);

if 1
mode = netcdf.getConstant('CLOBBER');
mode = bitor(mode,netcdf.getConstant('NETCDF4'));
temp=netcdf.create(fid.grid,mode);
netcdf.close(temp);

nccreate(fid.grid,'lon','Dimensions',{'x',Nx});
ncwrite(fid.grid,'lon',lon);

nccreate(fid.grid,'lat','Dimensions',{'y',Ny});
ncwrite(fid.grid,'lat',lat);

nccreate(fid.grid,'H','Dimensions',{'x',Nx,'y',Ny});
ncwrite(fid.grid,'H',H);

nccreate(fid.grid,'H_mix','Dimensions',{'x',Nx,'y',Ny});
ncwrite(fid.grid,'H_mix',H_mix);

nccreate(fid.grid,'c','Dimensions',{'x',Nx,'y',Ny,'n',Nm});
nccreate(fid.grid,'phi_stress','Dimensions',{'x',Nx,'y',Ny,'n',Nm});
nccreate(fid.grid,'phi_surf','Dimensions',{'x',Nx,'y',Ny,'n',Nm});
nccreate(fid.grid,'phi_bott','Dimensions',{'x',Nx,'y',Ny,'n',Nm});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate modes
Nxs=length(strat.lon);
Nys=length(strat.lat);
[Nx Ny]=size(H);
Nz=size(strat.z);

for i=1:Nx    
    c=zeros([1 Ny Nm]);
    phi_stress=zeros([1 Ny Nm]);
    phi_surf=zeros([1 Ny Nm]);
    phi_bott=zeros([1 Ny Nm]);
        
    % Interpolate stratification
    indx=[-1 0 1]+dsearchn(strat.lon,lon(i));
    N20=NaN([Ny Nz]);
    for k=1:Nz
        tmp=strat.N2(indx,:,k);
        N20(:,k)=interp2(strat.lon(indx),strat.lat,tmp',lon(i),lat);
    end
    
    for j=1:Ny           
        
        if H(i,j)>5       
            good=isfinite(N20(j,:)');
            if sum(good)>=2
                % Interpolate coarse N2 onto fine depth grid
                N2=interp1(strat.z(good),N20(j,good)',z).';
                N2(N2<1e-8)=1e-8;
                ind_z=find(isfinite(N2),1,'first');
                N2(1:ind_z)=N2(ind_z);
                ind_z=find(isfinite(N2),1,'last');
                N2(ind_z:end)=N2(ind_z);

                % Identify the proper number of modes to compute
                Nm0_ij=min([Nm0 floor(H(i,j)/(2*dz))+1]);
                Nm_ij=min([Nm Nm0_ij]);
                
                % Compute the modes
                ind_z=min([round(H(i,j)/dz) H_max/dz]);
                [PHI,c(1,j,1:Nm_ij)]=MODES_FAST(dz,N2(1:ind_z),Nm_ij,Nm0_ij);
                
                % Phi_surf is really an integral of the stress profile
                dZdz=z*0;
                alpha=1-(H_mix(i,j)/zT(i,j))^2;

                keep=z<zT(i,j);
                dZdz(keep)=2/(alpha*zT(i,j))*(1-z(keep)/zT(i,j));

                keep=z<=H_mix(i,j);
                dZdz(keep)=2/(zT(i,j)+H_mix(i,j));                                            
                
                phi_stress(1,j,1:Nm_ij)=sum(dZdz(1:ind_z)'.*PHI)'*dz; % This is the projection of the surface stress divergence onto the mode
                
                % % Legacy code: Just a mixed layer
                %dZdz2=z*0; 
                %dZdz2(z<=H_mix(i,j))=1/H_mix(i,j);
                %phi_stress(1,j,1:Nm_ij)=sum(dZdz(1:ind_z)'.*PHI)'*dz; 

                % The average value through the mixed layer (useful for calculating surface aka mixed-layer velocity)
                keep=z(1:ind_z)<=H_mix(i,j);                
                phi_surf(1,j,1:Nm_ij)=mean(PHI(keep,:))';
                
                % The un-weighted bottom value
                phi_bott(1,j,1:Nm_ij)=PHI(end,:)';
            end
            
        end
    end
    ncwrite(fid.grid,'c',c,[i 1 1]);
    ncwrite(fid.grid,'phi_stress',phi_stress,[i 1 1]);
    ncwrite(fid.grid,'phi_surf',phi_surf,[i 1 1]);
    ncwrite(fid.grid,'phi_bott',phi_bott,[i 1 1]);
   
    PROGRESS_BAR(i,1:Nx);
end


