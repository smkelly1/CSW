% This is meant to be flexible post_processing script

fid.csw='./out';
fid.SSH='./SSH';
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
for n=0:NPX*NPY-1
    start.x=floor(n-floor(n/NPX)*NPX)*NX;
    start.y=floor(n/NPX)*NY;
    
    Ai(start.x+[1:NX],start.y+[1:NY])=ncread([fid.csw,'.',num2str(n,'%03d'),'.nc'],'SSH_amp',[1 1 1 cycle-10],[NX NY NM 1]);
    
    out.SSH_amp(start.x+[1:NX],start.y+[1:NY])=ncread([fid.csw,'.',num2str(n,'%03d'),'.nc'],'SSH_amp',[1 1 1 cycle],[NX NY NM 1]);
    out.SSH_phase(start.x+[1:NX],start.y+[1:NY])=ncread([fid.csw,'.',num2str(n,'%03d'),'.nc'],'SSH_phase',[1 1 1 cycle],[NX NY NM 1]);
end

% Mask points where amplitude changes by more than 1 cm and 10%
tmp1=(out.SSH_amp-Ai)./Ai;
tmp2=(out.SSH_amp-Ai);
out.mask=(H>1000 & ~(abs(tmp1)>0.5 & abs(tmp2)*100>2));
out.mask100=(H>100 & ~(abs(tmp1)>0.1 & abs(tmp2)>0.01));

% Save the rest of the data
out.cycle=cycle;
out.H=H;
out.lon=lon;
out.lat=lat;

save(fid.SSH,'out','-v7.3');
