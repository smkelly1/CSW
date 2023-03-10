% This is meant to be flexible post_processing script
clear

fid.out='./diag';
fid.grid='../../22-12_grid/10th_deg_OCT_grid.nc';

% Last cycle to analyze
period=23;

% Info about grid
res=10;
NM=8;
NPX=2;
NPY=2;

if res==25
	NX0=9000; % Get from csw.h
    NY0=3650;
elseif res==10
	NX0=3600; % Get from csw.h
    NY0=1460;
 elseif res==5
	NX0=1800; % Get from csw.h
    NY0=730;
end
NX=NX0/NPX;
NY=NY0/NPY;

% Load bathymetry
H=ncread(fid.grid,'H');H=H(2:end-1,1+[1:NY0]);
lon=ncread(fid.grid,'lon');lon=lon(2:end-1);
lat=ncread(fid.grid,'lat');lat=lat(1+[1:NY0]);

a=6371e3;
dx=(1/res)/180*pi*a*repmat(cos(lat'/180*pi),[NX*NPX 1]);
dy=(1/res)/180*pi*a;

% Mask points where low-pass amplitude is more than 1 cm
% eta_low=zeros([NX*NPX NY*NPY]);
% for n=0:NPX*NPY-1
%     start.x=floor(n-floor(n/NPX)*NPX)*NX;
%     start.y=floor(n/NPX)*NY;
%     ssh0(start.x+[1:NX],start.y+[1:NY])=ncread([fid.out,'.',num2str(n,'%03d'),'.nc'],'flag_growth',[1 1 1 period],[NX NY 1 1]);
% end
out.mask=(H>1000);% & ssh0<0.01);
out.mask100=(H>100);% & ssh0<0.01);

% Extract
clear int 

int.pct_bad=(sum(sum(H>1000))-sum(sum(out.mask)))/sum(sum(H>1000))*100;
int100.pct_bad=(sum(sum(H>100))-sum(sum(out.mask100)))/sum(sum(H>100))*100;

out.E=zeros([NX*NPX NY*NPY]);
out.KE=zeros([NX*NPX NY*NPY]);
out.W=zeros([NX*NPX NY*NPY]);
out.C=zeros([NX*NPX NY*NPY]);
out.Cn=zeros([NX*NPX NY*NPY]);
out.D=zeros([NX*NPX NY*NPY]);
out.divF=zeros([NX*NPX NY*NPY]);
out.error=zeros([NX*NPX NY*NPY]);

time=1:period;
for j=1:length(time) 
    for i=1:NM
        for n=0:NPX*NPY-1
            start.x=floor(n-floor(n/NPX)*NPX)*NX;
            start.y=floor(n/NPX)*NY;
            
            out.E(start.x+[1:NX],start.y+[1:NY],:)=ncread([fid.out,'.',num2str(n,'%03d'),'.nc'],'KE',[1 1 i time(j)],[NX NY 1 1])+ncread([fid.out,'.',num2str(n,'%03d'),'.nc'],'PE',[1 1 i time(j)],[NX NY 1 1]);
            out.KE(start.x+[1:NX],start.y+[1:NY],:)=ncread([fid.out,'.',num2str(n,'%03d'),'.nc'],'KE',[1 1 i time(j)],[NX NY 1 1]);
            out.W(start.x+[1:NX],start.y+[1:NY],:)=ncread([fid.out,'.',num2str(n,'%03d'),'.nc'],'W',[1 1 i time(j)],[NX NY 1 1]);
            out.C(start.x+[1:NX],start.y+[1:NY],:)=ncread([fid.out,'.',num2str(n,'%03d'),'.nc'],'C0',[1 1 i time(j)],[NX NY 1 1]);
            out.Cn(start.x+[1:NX],start.y+[1:NY],:)=ncread([fid.out,'.',num2str(n,'%03d'),'.nc'],'Cn',[1 1 i time(j)],[NX NY 1 1]);
            out.D(start.x+[1:NX],start.y+[1:NY],:)=ncread([fid.out,'.',num2str(n,'%03d'),'.nc'],'D',[1 1 i time(j)],[NX NY 1 1]);
            out.divF(start.x+[1:NX],start.y+[1:NY],:)=ncread([fid.out,'.',num2str(n,'%03d'),'.nc'],'divF',[1 1 i time(j)],[NX NY 1 1]);
            out.error(start.x+[1:NX],start.y+[1:NY],:)=ncread([fid.out,'.',num2str(n,'%03d'),'.nc'],'error',[1 1 i time(j)],[NX NY 1 1]);
        end
        
        int.E(i,j)=squeeze(nansum(nansum(out.mask.*out.E.*dx.*dy))*1e-9);
        int.KE(i,j)=squeeze(nansum(nansum(out.mask.*out.KE.*dx.*dy))*1e-9);
        int.W(i,j)=squeeze(nansum(nansum(out.mask.*out.W.*dx.*dy))*1e-9);
        int.C(i,j)=squeeze(nansum(nansum(out.mask.*out.C.*dx.*dy))*1e-9);
        int.Cn(i,j)=squeeze(nansum(nansum(out.mask.*out.Cn.*dx.*dy))*1e-9);
        int.divF(i,j)=squeeze(nansum(nansum(out.mask.*out.divF.*dx.*dy))*1e-9);
        int.D(i,j)=squeeze(nansum(nansum(out.mask.*out.D.*dx.*dy))*1e-9);
        int.error(i,j)=squeeze(nansum(nansum(out.mask.*out.error.*dx.*dy))*1e-9);
        
        int100.E(i,j)=squeeze(nansum(nansum(out.mask100.*out.E.*dx.*dy))*1e-9);
        int100.KE(i,j)=squeeze(nansum(nansum(out.mask100.*out.E.*dx.*dy))*1e-9);
        int100.W(i,j)=squeeze(nansum(nansum(out.mask100.*out.W.*dx.*dy))*1e-9);
        int100.C(i,j)=squeeze(nansum(nansum(out.mask100.*out.C.*dx.*dy))*1e-9);
        int100.Cn(i,j)=squeeze(nansum(nansum(out.mask100.*out.Cn.*dx.*dy))*1e-9);
        int100.divF(i,j)=squeeze(nansum(nansum(out.mask100.*out.divF.*dx.*dy))*1e-9);
        int100.D(i,j)=squeeze(nansum(nansum(out.mask100.*out.D.*dx.*dy))*1e-9);
        int100.error(i,j)=squeeze(nansum(nansum(out.mask100.*out.error.*dx.*dy))*1e-9);
    end
end
int.time=time;
int100.time=time;
out.period=period;

out.H=H;
out.lon=lon;
out.lat=lat;

%plot(int.time,int.C-int.divF,int.time,int.D,int.time,int.error)
%plot(time,int100.C-int100.divF,time,int100.D)
%plot(int.time,int.E)

%save(fid.out,'out','int','int100','-v7.3');

%%
KE=mean(int.KE,2);
W=sum(int.KE,2);

GM=1./([1:NM].^2+3^2)/sum(1./([1:NM].^2+3^2));

figure(1);clf;

subplot(2,1,1)
plot(1:NM,log10(W),'k-+',1:NM,log10(5e9*GM),'r-')

subplot(2,1,2)
plot(1:NM,log10(KE),'k-+',1:NM,log10(1e8*GM),'r-')



%%

W=sum(int.W(:,1:10),2)*(32*3600)*1e-6;
D=sum(int.D(:,1:10),2)*(32*3600)*1e-6;
Cn=sum(int.Cn(:,1:10),2)*(32*3600)*1e-6;
divF=sum(int.divF(:,1:10),2)*(32*3600)*1e-6;
E=int.E(:,10)*1e-6;
residual=W+Cn-D-divF-E;

figure(1);clf;

subplot(3,1,1)
n=1;
bar(1:6,[W(n) D(n) Cn(n) divF(n) E(n) residual(n)],'k')
ylabel('PJ')

subplot(3,1,2)
n=2;
bar(1:6,[W(n) D(n) Cn(n) divF(n) E(n) residual(n)],'k')
ylabel('PJ')

subplot(3,1,3)
n=3:8;
bar(1:6,[sum(W(n)) sum(D(n)) sum(Cn(n)) sum(divF(n)) sum(E(n)) sum(residual(n))],'k')
ylabel('PJ')
set(gca,'xtick',[1:6],'xticklabel',{'Wind Work','Dissipation','Topo. scattering','Coastal dissipation','Remaining energy','Residual'})

%%

bar(1:6,[W D Cn divF E residual],'stacked')
set(gca,'xtick',[1:6],'xticklabel',{'Wind Work','Dissipation','Topo. scattering','Coastal dissipation','Remaining energy','Residual'})


