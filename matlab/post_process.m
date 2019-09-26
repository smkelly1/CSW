% This is meant to be flexible post_processing script

fid.out='./out';
fid.grid='../../18-6_grids/25th_deg_SS_bathy.nc';

% Last cycle to analyze
cycle=50;

% Info about grid
res=1/25;
NM=1;
NPX=4;
NPY=4;
NX0=9000; % Get from csw.h
NY0=3648;
NX=NX0/NPX;
NY=NY0/NPY;

% Load bathymetry
H=ncread(fid.grid,'H');H=H(:,1+[1:NY0]);
lon=ncread(fid.grid,'lon');
lat=ncread(fid.grid,'lat');lat=lat(1+[1:NY0]);

a=6371e3;
dx=res/180*pi*a*repmat(cos(lat'/180*pi),[NX*NPX 1]);
dy=res/180*pi*a;

% Create mask to determine exponential growth
Ai=zeros([NX*NPX NY*NPY]);
out.SSH_AMP=zeros([NX*NPX NY*NPY]);
out.SSH_phase=zeros([NX*NPX NY*NPY]);
out.U_AMP=zeros([NX*NPX NY*NPY]);
out.U_phase=zeros([NX*NPX NY*NPY]);
out.V_AMP=zeros([NX*NPX NY*NPY]);
out.V_phase=zeros([NX*NPX NY*NPY]);
for n=0:NPX*NPY-1
    start.x=floor(n-floor(n/NPX)*NPX)*NX;
    start.y=floor(n/NPX)*NY;
    
    Ai(start.x+[1:NX],start.y+[1:NY])=ncread([fid.out,'.',num2str(n,'%03d'),'.nc'],'SSH_amp',[1 1 1 cycle-10],[NX NY NM 1]);
    
    out.SSH_amp(start.x+[1:NX],start.y+[1:NY])=ncread([fid.out,'.',num2str(n,'%03d'),'.nc'],'SSH_amp',[1 1 1 cycle],[NX NY NM 1]);
    out.SSH_phase(start.x+[1:NX],start.y+[1:NY])=ncread([fid.out,'.',num2str(n,'%03d'),'.nc'],'SSH_phase',[1 1 1 cycle],[NX NY NM 1]);
    
    out.U_amp(start.x+[1:NX],start.y+[1:NY])=ncread([fid.out,'.',num2str(n,'%03d'),'.nc'],'U_amp',[1 1 1 cycle],[NX NY NM 1]);
    out.U_phase(start.x+[1:NX],start.y+[1:NY])=ncread([fid.out,'.',num2str(n,'%03d'),'.nc'],'U_phase',[1 1 1 cycle],[NX NY NM 1]);
    
    out.V_amp(start.x+[1:NX],start.y+[1:NY])=ncread([fid.out,'.',num2str(n,'%03d'),'.nc'],'V_amp',[1 1 1 cycle],[NX NY NM 1]);
    out.V_phase(start.x+[1:NX],start.y+[1:NY])=ncread([fid.out,'.',num2str(n,'%03d'),'.nc'],'V_phase',[1 1 1 cycle],[NX NY NM 1]);
    
    error(start.x+[1:NX],start.y+[1:NY])=ncread([fid.out,'.',num2str(n,'%03d'),'.nc'],'error',[1 1 1 cycle],[NX NY NM 1]);
end

% Mask points where amplitude changes by more than 1 cm and 10%
tmp1=(out.SSH_amp-Ai)./Ai;
tmp2=(out.SSH_amp-Ai);
%out.mask=(H>1000 & ~(abs(tmp1)>0.1 & abs(tmp2)>0.005) & error<0.1);
out.mask=(H>1000 & ~(abs(tmp1)>0.5 & abs(tmp2)*100>2) & error<0.1);
out.mask100=(H>100 & ~(abs(tmp1)>0.1 & abs(tmp2)>0.01) & error<0.1);

% Extract
clear int 

int.pct_bad=(sum(sum(H>1000))-sum(sum(out.mask)))/sum(sum(H>1000))*100;
int100.pct_bad=(sum(sum(H>100))-sum(sum(out.mask100)))/sum(sum(H>100))*100;

out.E=zeros([NX*NPX NY*NPY]);
out.KE=zeros([NX*NPX NY*NPY]);
out.C=zeros([NX*NPX NY*NPY]);
out.Cn=zeros([NX*NPX NY*NPY]);
out.D=zeros([NX*NPX NY*NPY]);
out.divF=zeros([NX*NPX NY*NPY length(cycle)]);
out.error=zeros([NX*NPX NY*NPY length(cycle)]);

time=1:cycle;
for j=1:length(time)        
    for n=0:NPX*NPY-1
        start.x=floor(n-floor(n/NPX)*NPX)*NX;
        start.y=floor(n/NPX)*NY;
        
        out.E(start.x+[1:NX],start.y+[1:NY])=ncread([fid.out,'.',num2str(n,'%03d'),'.nc'],'KE',[1 1 1 time(j)],[NX NY NM 1])+ncread([fid.out,'.',num2str(n,'%03d'),'.nc'],'PE',[1 1 1 time(j)],[NX NY NM 1]);
        out.KE(start.x+[1:NX],start.y+[1:NY])=ncread([fid.out,'.',num2str(n,'%03d'),'.nc'],'KE',[1 1 1 time(j)],[NX NY NM 1]);
        out.C(start.x+[1:NX],start.y+[1:NY])=ncread([fid.out,'.',num2str(n,'%03d'),'.nc'],'C0',[1 1 1 time(j)],[NX NY NM 1]);
        out.Cn(start.x+[1:NX],start.y+[1:NY])=ncread([fid.out,'.',num2str(n,'%03d'),'.nc'],'Cn',[1 1 1 time(j)],[NX NY NM 1]);
        out.D(start.x+[1:NX],start.y+[1:NY])=ncread([fid.out,'.',num2str(n,'%03d'),'.nc'],'D',[1 1 1 time(j)],[NX NY NM 1]);
        out.divF(start.x+[1:NX],start.y+[1:NY])=ncread([fid.out,'.',num2str(n,'%03d'),'.nc'],'divF',[1 1 1 time(j)],[NX NY NM 1]);
        out.error(start.x+[1:NX],start.y+[1:NY])=ncread([fid.out,'.',num2str(n,'%03d'),'.nc'],'error',[1 1 1 time(j)],[NX NY NM 1]);
    end
    
    int.E(j)=squeeze(nansum(nansum(out.mask.*out.E.*dx.*dy))*1e-9);
    int.KE(j)=squeeze(nansum(nansum(out.mask.*out.KE.*dx.*dy))*1e-9);
    int.C(j)=squeeze(nansum(nansum(out.mask.*out.C.*dx.*dy))*1e-9);
    int.Cn(j)=squeeze(nansum(nansum(out.mask.*out.Cn.*dx.*dy))*1e-9);
    int.divF(j)=squeeze(nansum(nansum(out.mask.*out.divF.*dx.*dy))*1e-9);
    int.D(j)=squeeze(nansum(nansum(out.mask.*out.D.*dx.*dy))*1e-9);
    int.error(j)=squeeze(nansum(nansum(out.mask.*out.error.*dx.*dy))*1e-9); 

    int100.E(j)=squeeze(nansum(nansum(out.mask100.*out.E.*dx.*dy))*1e-9);
    int100.KE(j)=squeeze(nansum(nansum(out.mask100.*out.E.*dx.*dy))*1e-9);
    int100.C(j)=squeeze(nansum(nansum(out.mask100.*out.C.*dx.*dy))*1e-9);
    int100.Cn(j)=squeeze(nansum(nansum(out.mask100.*out.Cn.*dx.*dy))*1e-9);
    int100.divF(j)=squeeze(nansum(nansum(out.mask100.*out.divF.*dx.*dy))*1e-9);
    int100.D(j)=squeeze(nansum(nansum(out.mask100.*out.D.*dx.*dy))*1e-9);
    int100.error(j)=squeeze(nansum(nansum(out.mask100.*out.error.*dx.*dy))*1e-9);    
end
int.time=time;
int100.time=time;
out.cycle=cycle;

out.H=H;
out.lon=lon;
out.lat=lat;

%plot(int.time,int.C-int.divF,int.time,int.D,int.time,int.error)
%plot(time,int100.C-int100.divF,time,int100.D)
%plot(int.time,int.E)

save(fid.out,'out','int','int100','-v7.3');
