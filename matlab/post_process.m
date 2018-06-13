% This is meant to be flexible post_processing script

%fid.out='./run_HAMTIDE_r3/25th_deg_out.';
%fid.grid='../18-5_global_grids/25th_deg_grid.nc';

%res=1/25;
%NM=1;
%NPX=4;
%NPY=4;

%H=ncread(fid.grid,'H');H=H(2:end-1,2:end-3);
%lon=ncread(fid.grid,'lon');lon=lon(2:end-1);
%lat=ncread(fid.grid,'lat');lat=lat(2:end-3);

%[NX0 NY0]=size(H);
%NX=NX0/NPX;
%NY=floor(NY0/NPY);

%a=6371e3;
%dx=res/180*pi*a*repmat(cos(lat'/180*pi),[NX*NPX 1]);
%dy=res/180*pi*a;

%%% 
%cycle=[20 30];
%p=zeros([NX*NPX NY*NPY length(cycle)]);
%C=zeros([NX*NPX NY*NPY length(cycle)]);
%Cn=zeros([NX*NPX NY*NPY length(cycle)]);
%D=zeros([NX*NPX NY*NPY length(cycle)]);
%divF=zeros([NX*NPX NY*NPY length(cycle)]);
%error=zeros([NX*NPX NY*NPY length(cycle)]);

%for j=1:length(cycle)
    
    %ind.t=cycle(j);
    %NT=1;
        
    %for n=0:NPX*NPY-1
        %start.x=floor(n-floor(n/NPX)*NPX)*NX;
        %start.y=floor(n/NPX)*NY;
        %C(start.x+[1:NX],start.y+[1:NY],j)=ncread([fid.out,num2str(n,'%03d'),'.nc'],'C0',[1 1 1 ind.t],[NX NY NM NT]);
        %Cn(start.x+[1:NX],start.y+[1:NY],j)=ncread([fid.out,num2str(n,'%03d'),'.nc'],'Cn',[1 1 1 ind.t],[NX NY NM NT]);
        %D(start.x+[1:NX],start.y+[1:NY],j)=ncread([fid.out,num2str(n,'%03d'),'.nc'],'D',[1 1 1 ind.t],[NX NY NM NT]);
        %divF(start.x+[1:NX],start.y+[1:NY],j)=ncread([fid.out,num2str(n,'%03d'),'.nc'],'divF',[1 1 1 ind.t],[NX NY NM NT]);
        %error(start.x+[1:NX],start.y+[1:NY],j)=ncread([fid.out,num2str(n,'%03d'),'.nc'],'error',[1 1 1 ind.t],[NX NY NM NT]);
    %end
%end

%clear int;

%tmp1=(C(:,:,end)-C(:,:,1))./C(:,:,end);
%tmp1(isinf(tmp1))=0;
%tmp2=(C(:,:,end)-C(:,:,1));

%mask=(H>100 & ~(abs(tmp1)>0.1 & abs(tmp2)>0.01));
%int.C100=squeeze(nansum(nansum(mask.*C(:,:,end).*dx.*dy))*1e-9);
%int.Cn100=squeeze(nansum(nansum(mask.*Cn(:,:,end).*dx.*dy))*1e-9);
%int.divF100=squeeze(nansum(nansum(mask.*divF(:,:,end).*dx.*dy))*1e-9);
%int.D100=squeeze(nansum(nansum(mask.*D(:,:,end).*dx.*dy))*1e-9);
%int.err100=squeeze(nansum(nansum(mask.*error(:,:,end).*dx.*dy))*1e-9);

%mask=(H>1000 & ~(abs(tmp1)>0.1 & abs(tmp2)>0.01));
%int.C=squeeze(nansum(nansum(mask.*C(:,:,end).*dx.*dy))*1e-9);
%int.Cn=squeeze(nansum(nansum(mask.*Cn(:,:,end).*dx.*dy))*1e-9);
%int.divF=squeeze(nansum(nansum(mask.*divF(:,:,end).*dx.*dy))*1e-9);
%int.D=squeeze(nansum(nansum(mask.*D(:,:,end).*dx.*dy))*1e-9);
%int.err=squeeze(nansum(nansum(mask.*error(:,:,end).*dx.*dy))*1e-9)










fid.out='./run_GOT_r30/25th_deg_out.';
%fid.out='./run_HAMTIDE_r3/25th_deg_out.';
%fid.in='../10th_deg_in.nc';
%fid.grid='../../../18-5_global_grids/25th_deg_grid.nc';
fid.grid='../18-5_global_grids/25th_deg_grid.nc';

res=1/25;
NM=1;
NPX=8;
NPY=8;

H=ncread(fid.grid,'H');H=H(2:end-1,2:end-3);
lon=ncread(fid.grid,'lon');lon=lon(2:end-1);
lat=ncread(fid.grid,'lat');lat=lat(2:end-3);

%mask=(H>500);
mask=(H>1000);
C_thresh=Inf;%0.5;
%C_thresh=Inf;%0.5;
%C_thresh=Inf;%0.5;

[NX0 NY0]=size(H);
NX=NX0/NPX;
NY=floor(NY0/NPY);

a=6371e3;
dx=res/180*pi*a*repmat(cos(lat'/180*pi),[NX*NPX 1]);
dy=res/180*pi*a;

%% 
clear int;
cycle=100;%[1:30];
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
        %p(start.x+[1:NX],start.y+[1:NY],:,:)=nanmean(ncread([fid.out,num2str(n,'%03d'),'.nc'],'p',[1 1 1 ind.t],[NX NY NM NT]),4);
        C(start.x+[1:NX],start.y+[1:NY],:,:)=nanmean(ncread([fid.out,num2str(n,'%03d'),'.nc'],'C0',[1 1 1 ind.t],[NX NY NM NT]),4);
        Cn(start.x+[1:NX],start.y+[1:NY],:,:)=nanmean(ncread([fid.out,num2str(n,'%03d'),'.nc'],'Cn',[1 1 1 ind.t],[NX NY NM NT]),4);
        D(start.x+[1:NX],start.y+[1:NY],:,:)=nanmean(ncread([fid.out,num2str(n,'%03d'),'.nc'],'D',[1 1 1 ind.t],[NX NY NM NT]),4);
        divF(start.x+[1:NX],start.y+[1:NY],:,:)=nanmean(ncread([fid.out,num2str(n,'%03d'),'.nc'],'divF',[1 1 1 ind.t],[NX NY NM NT]),4);
        error(start.x+[1:NX],start.y+[1:NY],:,:)=nanmean(ncread([fid.out,num2str(n,'%03d'),'.nc'],'error',[1 1 1 ind.t],[NX NY NM NT]),4);
    end
    
    C(C>C_thresh)=C_thresh;    
    %C(abs(C)>C_thresh)=NaN;    
    int.C(j,:)=squeeze(nansum(nansum(mask.*C.*dx.*dy))*1e-9);
    int.Cn(j,:)=squeeze(nansum(nansum(mask.*Cn.*dx.*dy))*1e-9);
    int.D(j,:)=squeeze(nansum(nansum(mask.*D.*dx.*dy))*1e-9);
    
    % Use a mask to avoid shallow regions
    int.divF(j,:)=squeeze(nansum(nansum(mask.*divF.*dx.*dy))*1e-9);
    int.error(j,:)=squeeze(nansum(nansum(mask.*error.*dx.*dy))*1e-9);
end

%%
%plot(cycle,int.C,'k',cycle,int.D,'r',cycle,int.divF,'b',cycle,int.error,'k--')

plot(cycle,sum(int.C,2),'k',cycle,sum(int.D,2),'r',cycle,sum(int.divF,2),'b',cycle,sum(int.error,2),'k--')

Nd=length(int.C);
cycle=1:Nd;
plot(cycle,sum(int.C,2),'k',cycle,sum(int.D,2),'r',cycle,sum(int.divF,2),'b',cycle,sum(int.error,2),'k--',cycle,sum(int.Cn*1000,2),'r--')

%%
tmp=error;
%tmp(tmp<0)=0;
Clims=[-.1 .1];

figure(1);clf;colormap(REDBLUE5);

subplot(2,1,1)
imagesc(lon,lat,AVE2D(divF(:,:,1),5)');shading flat;caxis(Clims);colorbar;hold on;set(gca,'ydir','normal')
contour(lon,lat,H',[1 1],'k');

subplot(2,1,2)
imagesc(lon,lat,AVE2D(tmp(:,:,1),5)');shading flat;caxis(Clims);colorbar;hold on;set(gca,'ydir','normal')
contour(lon,lat,H',[1 1],'k');

linkaxes;

%%
plims=[-100 100];

figure(1);clf;colormap(REDBLUE5);
pcolor(lon,lat,p(:,:,1)');shading flat;caxis(plims);colorbar;hold on;
contour(lon,lat,H',[1 1],'k');

%%
for i=32;%1:22
ind.t=i;
NT=1;

p=zeros([NX*NPX NY*NPY]);
u=zeros([NX*NPX NY*NPY]);
v=zeros([NX*NPX NY*NPY]);
C0=zeros([NX*NPX NY*NPY]);

for n=0:NPX*NPY-1
    start.x=floor(n-floor(n/NPX)*NPX)*NX;
    start.y=floor(n/NPX)*NY;
    p(start.x+[1:NX],start.y+[1:NY])=ncread([fid.out,num2str(n,'%03d'),'.nc'],'p',[1 1 1 ind.t],[NX NY 1 1]);
    u(start.x+[1:NX],start.y+[1:NY])=ncread([fid.out,num2str(n,'%03d'),'.nc'],'u',[1 1 1 ind.t],[NX NY 1 1]);
    v(start.x+[1:NX],start.y+[1:NY])=ncread([fid.out,num2str(n,'%03d'),'.nc'],'v',[1 1 1 ind.t],[NX NY 1 1]);
    C0(start.x+[1:NX],start.y+[1:NY])=ncread([fid.out,num2str(n,'%03d'),'.nc'],'C0',[1 1 1 ind.t],[NX NY 1 1]);
end

imagesc(p');set(gca,'ydir','normal');caxis([-1 1]*100);colormap(REDBLUE5)
%imagesc(u');set(gca,'ydir','normal');caxis([-1 1]*0.01);colormap(REDBLUE5)
%imagesc(C0');set(gca,'ydir','normal');caxis([-1 1]*0.01);colormap(REDBLUE5)
pause(1);
end


%tmp=C;
%tmp(abs(tmp)>1)=1;
%mask=H>500;
%squeeze(nansum(nansum(mask.*tmp.*dx.*dy))*1e-9)

