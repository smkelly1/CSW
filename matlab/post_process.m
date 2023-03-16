% This is meant to be flexible post_processing script
clear

fid.out='./diag';
%fid.grid='../../22-12_grid/10th_deg_OCT_grid.nc';
fid.grid='../../22-12_grid/10th_deg_JUN_grid.nc';

% Last cycle to analyze
Nt=31;

% Info about grid
res=10;
Nm=8;
NPX=2;
NPY=2;
DT=24*3600;

if res==25
	Nx=9000; % Get from csw.h
    Ny=3650;
elseif res==10
	Nx=3600; % Get from csw.h
    Ny=1460;
 elseif res==5
	Nx=1800; % Get from csw.h
    Ny=730;
end
NXp=Nx/NPX;
NYp=Ny/NPY;

% Load bathymetry
H_tot=ncread(fid.grid,'H'); 
H_tot=H_tot(:,1:Ny+2); 
H=H_tot(2:end-1,2:end-1);
lon=ncread(fid.grid,'lon');lon=lon(2:end-1);
lat=ncread(fid.grid,'lat');lat=lat(1+[1:Ny]);

a=6371e3;
dx=(1/res)/180*pi*a*repmat(cos(lat'/180*pi),[Nx 1]);
dy=(1/res)/180*pi*a;

% Load energy terms
clear int 

out.KE=zeros([Nx Ny]);
out.PE=zeros([Nx Ny]);
out.W=zeros([Nx Ny]);
out.C=zeros([Nx Ny]);
out.Cn=zeros([Nx Ny]);
out.D=zeros([Nx Ny]);
out.coast=zeros([Nx Ny]);
out.growth=zeros([Nx Ny]);
out.up=zeros([Nx Ny]);
out.vp=zeros([Nx Ny]);
out.bad=zeros([Nx Ny]);

for k=1:Nt 
    for n=1:Nm
        for p=0:(NPX*NPY-1)
            start.x=floor(p-floor(p/NPX)*NPX)*NXp;
            start.y=floor(p/NPX)*NYp;
            
            out.PE(start.x+[1:NXp],start.y+[1:NYp])=ncread([fid.out,'.',num2str(p,'%03d'),'.nc'],'PE',[1 1 n k],[NXp NYp 1 1]);
            out.KE(start.x+[1:NXp],start.y+[1:NYp])=ncread([fid.out,'.',num2str(p,'%03d'),'.nc'],'KE',[1 1 n k],[NXp NYp 1 1]);
            out.W(start.x+[1:NXp],start.y+[1:NYp])=ncread([fid.out,'.',num2str(p,'%03d'),'.nc'],'W',[1 1 n k],[NXp NYp 1 1]);
            out.C(start.x+[1:NXp],start.y+[1:NYp])=ncread([fid.out,'.',num2str(p,'%03d'),'.nc'],'C0',[1 1 n k],[NXp NYp 1 1]);
            out.Cn(start.x+[1:NXp],start.y+[1:NYp])=ncread([fid.out,'.',num2str(p,'%03d'),'.nc'],'Cn',[1 1 n k],[NXp NYp 1 1]);
            out.D(start.x+[1:NXp],start.y+[1:NYp])=ncread([fid.out,'.',num2str(p,'%03d'),'.nc'],'D',[1 1 n k],[NXp NYp 1 1]);
            out.up(start.x+[1:NXp],start.y+[1:NYp])=ncread([fid.out,'.',num2str(p,'%03d'),'.nc'],'up',[1 1 n k],[NXp NYp 1 1]);
            out.vp(start.x+[1:NXp],start.y+[1:NYp])=ncread([fid.out,'.',num2str(p,'%03d'),'.nc'],'vp',[1 1 n k],[NXp NYp 1 1]);
            out.bad(start.x+[1:NXp],start.y+[1:NYp])=ncread([fid.out,'.',num2str(p,'%03d'),'.nc'],'flag_growth',[1 1 k],[NXp NYp 1]);
        end
        out.E=out.KE+out.PE;
        out.D(out.D<0)=0;
        out.D(:,end)=0;
               
        % Figure out onshore flux at 1000m isobath
        out.coast=out.coast*0;        
        mask=(1000<H);
        mask=[mask(end,:); mask; mask(1,:)];
        mask=[mask(:,1) mask mask(:,end)];
        for i=1:Nx
            for j=1:Ny
                if mask(i+1,j+1)==1
                    if mask(i,j+1)==0
                        out.coast(i,j)=out.coast(i,j)-out.up(i,j)*dy;
                    end
                    if mask(i+2,j+1)==0
                        out.coast(i,j)=out.coast(i,j)+out.up(i,j)*dy;
                    end
                    if mask(i+1,j)==0
                        out.coast(i,j)=out.coast(i,j)-out.vp(i,j)*dx(i,j);
                    end
                    if mask(i+1,j+2)==0
                        out.coast(i,j)=out.coast(i,j)+out.vp(i,j)*dx(i,j);
                    end
                end
            end
        end
        
        % Figure out energy flux out of bad regions
        out.growth=out.growth*0;
        mask2=(out.bad==0 & out.E<5000);
        mask2=[mask2(end,:); mask2; mask2(1,:)];
        mask2=[mask2(:,1) mask2 mask2(:,end)];
        for i=1:Nx
            for j=1:Ny
                if mask(i+1,j+1)==1 && mask2(i+1,j+1)==1
                    if mask(i,j+1)==1 && mask2(i,j+1)==0
                        out.growth(i,j)=out.growth(i,j)-out.up(i,j)*dy;
                    end
                    if mask(i+2,j+1)==1 && mask2(i+2,j+1)==0
                        out.growth(i,j)=out.growth(i,j)+out.up(i,j)*dy;
                    end
                    if mask(i+1,j)==1 && mask2(i+1,j)==0
                        out.growth(i,j)=out.growth(i,j)-out.vp(i,j)*dx(i,j);
                    end
                    if mask(i+1,j+2)==1 && mask2(i+1,j+2)==0
                        out.growth(i,j)=out.growth(i,j)+out.vp(i,j)*dx(i,j);
                    end
                end
            end
        end
        
        out.mask=(1000<H & out.bad==0 & out.E<5000);
        int.E(n,k)=squeeze(nansum(nansum(out.mask.*out.E.*dx.*dy))*1e-9);
        int.KE(n,k)=squeeze(nansum(nansum(out.mask.*out.KE.*dx.*dy))*1e-9);
        int.PE(n,k)=squeeze(nansum(nansum(out.mask.*out.KE.*dx.*dy))*1e-9);
        int.W(n,k)=squeeze(nansum(nansum(out.mask.*out.W.*dx.*dy))*1e-9);
        int.C(n,k)=squeeze(nansum(nansum(out.mask.*out.C.*dx.*dy))*1e-9);
        int.Cn(n,k)=squeeze(nansum(nansum(out.mask.*out.Cn.*dx.*dy))*1e-9);
        int.coast(n,k)=squeeze(nansum(nansum(out.mask.*out.coast))*1e-9);
        int.growth(n,k)=squeeze(nansum(nansum(out.mask.*out.growth))*1e-9);
        int.D(n,k)=squeeze(nansum(nansum(out.mask.*out.D.*dx.*dy))*1e-9);
    end
end
int.yday=(1:Nt)*DT/(24*3600);
out.period=1:Nt;

int.pct_bad=(sum(sum(H>1000))-sum(sum(out.mask)))/sum(sum(H>1000))*100;

out.H=H;
out.lon=lon;
out.lat=lat;

%save(fid.out,'out','int','int100','-v7.3');


%% Plot bar graphs of energy terms by mode

ind=1:31;
W=sum(int.W(:,ind),2)*(24*3600)*1e-6;
D=sum(int.D(:,ind),2)*(24*3600)*1e-6;
Cn=sum(int.Cn(:,ind),2)*(24*3600)*1e-6;
E=(int.E(:,ind(end))-int.E(:,ind(1)))*1e-6;
%E=int.E(:,ind(end))*1e-6;
coast=sum(int.coast(:,ind),2)*(24*3600)*1e-6;
growth=sum(int.growth(:,ind),2)*(24*3600)*1e-6;
%E=E-growth;
error=W+Cn-(D+coast+E);


ylims=[-10 160];

figure(1);clf;

subplot(3,1,1)
n=1:2;
bar(1:6,[sum(W(n)) sum(Cn(n)) sum(D(n))  sum(coast(n)) sum(E(n)) sum(error(n))],'k')
ylabel('PJ')
%ylim(ylims);

subplot(3,1,2)
n=3:4;
bar(1:6,[sum(W(n)) sum(Cn(n)) sum(D(n)) sum(coast(n)) sum(E(n)) sum(error(n))],'k')
ylabel('PJ')
%ylim(ylims);

subplot(3,1,3)
n=5:8;
bar(1:6,[sum(W(n)) sum(Cn(n)) sum(D(n)) sum(coast(n)) sum(E(n)) sum(error(n))],'k')
ylabel('PJ')
%ylim(ylims);
set(gca,'xtick',[1:6],'xticklabel',{'Wind Work','Topo. scattering','Dissipation','Flux to coast','Remaining energy','Residual'})

%%

bar(1:6,[W D Cn coast E error],'stacked')
set(gca,'xtick',[1:6],'xticklabel',{'Wind Work','Dissipation','Topo. scattering','Flux to coast','Remaining energy','Residual'})


