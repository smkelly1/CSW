% Load and integrate output data
clear

fid.out='25th_deg_out.';
fid.in='../25th_deg_in.nc';
fid.grid='../../17-6_global_grids/25th_deg_grid.nc';

res=1/25;
NM=4;
NPX=4;
NPY=1;

H=ncread(fid.grid,'H');H=H(2:end-1,2:end-1);
lon=ncread(fid.grid,'lon');lon=lon(2:end-1);
lat=ncread(fid.grid,'lat');lat=lat(2:end-1);
mask.u=ncread(fid.in,'mask_u',[1 1 1],[Inf Inf NM]);mask.u=(mask.u(1:end-1,:,:)+mask.u(2:end,:,:))/2;
mask.v=ncread(fid.in,'mask_v',[1 1 1],[Inf Inf NM]);mask.v=(mask.v(:,1:end-1,:)+mask.v(:,2:end,:))/2;
mask.p=ncread(fid.in,'mask_p',[1 1 1],[Inf Inf NM]);
mask.all=mask.u==1 & mask.v==1 & mask.p==1;

[NX0 NY0]=size(H);
NX=NX0/NPX;
NY=NY0/NPY;

a=6371e3;
dx=res/180*pi*a*repmat(cos(lat'/180*pi),[NX*NPX 1]);
dy=res/180*pi*a;

%% 
clear int;
cycle=[9];
for j=1:length(cycle)
    
    ind.t=cycle(j);
    NT=1;
    
    p=zeros([NX*NPX NY*NPY NM]);
    C=zeros([NX*NPX NY*NPY NM]);
    Cn=zeros([NX*NPX NY*NPY NM]);
    D=zeros([NX*NPX NY*NPY NM]);
    divF=zeros([NX*NPX NY*NPY NM]);
    error=zeros([NX*NPX NY*NPY NM]);
    
    for n=0:NPX*NPY-1
        start.x=floor(n-floor(n/NPX)*NPX)*NX;
        start.y=floor(n/NPX)*NY;
        p(start.x+[1:NX],start.y+[1:NY],:,:)=nanmean(ncread([fid.out,num2str(n,'%03d'),'.nc'],'p',[1 1 1 ind.t],[NX NY NM NT]),4);
        C(start.x+[1:NX],start.y+[1:NY],:,:)=nanmean(ncread([fid.out,num2str(n,'%03d'),'.nc'],'C0',[1 1 1 ind.t],[NX NY NM NT]),4);
        Cn(start.x+[1:NX],start.y+[1:NY],:,:)=nanmean(ncread([fid.out,num2str(n,'%03d'),'.nc'],'Cn',[1 1 1 ind.t],[NX NY NM NT]),4);
        D(start.x+[1:NX],start.y+[1:NY],:,:)=nanmean(ncread([fid.out,num2str(n,'%03d'),'.nc'],'D',[1 1 1 ind.t],[NX NY NM NT]),4);
        divF(start.x+[1:NX],start.y+[1:NY],:,:)=nanmean(ncread([fid.out,num2str(n,'%03d'),'.nc'],'divF',[1 1 1 ind.t],[NX NY NM NT]),4);
        error(start.x+[1:NX],start.y+[1:NY],:,:)=nanmean(ncread([fid.out,num2str(n,'%03d'),'.nc'],'error',[1 1 1 ind.t],[NX NY NM NT]),4);
    end
    
    int.C(j)=nansum(nansum(nansum(C(:,:,:),3).*dx.*dy))*1e-9;
    int.Cn(j)=nansum(nansum(nansum(Cn(:,:,1),3).*dx.*dy))*1e-9;
    int.D(j)=nansum(nansum(nansum(D(:,:,:),3).*dx.*dy))*1e-9;
    
    divF(mask.all(:,:,1:NM)~=1)=0;
    int.divF(j)=nansum(nansum(nansum(divF(:,:,:),3).*dx.*dy))*1e-9;

    error(mask.all(:,:,1:NM)~=1)=0;
    int.error(j)=nansum(nansum(nansum(error(:,:,:),3).*dx.*dy))*1e-9;
end

%%
plot(cycle,int.C,'k',cycle,int.D,'b',cycle,int.divF,'m',cycle,int.error,'r--')

%%
tmp=C;
%tmp(tmp<0)=0;
Clims=[-.1 .1];

figure(1);clf;colormap(REDBLUE5);

subplot(2,1,1)
pcolor(lon,lat,AVE2D(divF(:,:,1),1)');shading flat;caxis(Clims);colorbar;hold on;
contour(lon,lat,H',[1 1],'k');

subplot(2,1,2)
pcolor(lon,lat,AVE2D(tmp(:,:,1),1)');shading flat;caxis(Clims);colorbar;hold on;
contour(lon,lat,H',[1 1],'k');

linkaxes;

%%
plims=[-100 100];

figure(1);clf;colormap(REDBLUE5);
pcolor(lon,lat,p(:,:,1)');shading flat;caxis(plims);colorbar;hold on;
contour(lon,lat,H',[1 1],'k');

