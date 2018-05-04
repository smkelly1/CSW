function [U V eta omega]=TPXO(lon,lat,conID,day,vel_flag)
% function [U V eta omega]=TPXO(lon,lat,conID,day,vel_flag)
%
% Reads the NetCDF TPXOO8 atlas and returns the amplitudes of the tidal 
% transports and elevations at each lon and lat. The amplitudes can be 
% used to predict the tides for a given constituent using the formula
%
% U(x,y,t)=Real[U(x,y)*exp(ii*omega*t)]
%
% where t is seconds from the start of the input day
%
% INPUTS
% lon    [m x 1] or [m x n]	 Vector of longitudes [0 360]
% lat    [n x 1] or [m x n]	 Vector of latitudes [-90 90]
% conID  [k x 1]             Constituents to include (e.g., see list below, [1 3] returns M2 and K1)
% day	 [1 x 1]             The serial day (use 'datenum') for which the phase is set to zero
% vel_flag (0 or 1)          If 0 (default) return transports, if 1 return velocity
%
% OUTPUTS
% U, V   [m x n x k]  Complex amplitude of transport [m^2/s]
% eta    [m x n x k]  Complex amplitude of tidal elevation [m] 
% omega	 [k x 1]      The tidal frequency in rad/s
%
% written by Sam Kelly 5-1-2014 (smkelly@d.umn.edu)
% 
% based on Laurie Padman's TMD code, the Oregon State Tidal prediction 
% software (OTPS), Richard Ray's notes, and David Cartwright's calculations  

ii=complex(0,1);

if nargin<5
    vel_flag=0;
    if nargin<4
        day=0;        
        if nargin<3
            conID=[];
        end
    end
end

names={'m2','s2','k1','o1','n2','p1','k2','q1','m4'};

fid={'uv.m2_tpxo8_atlas_30c.nc',...
    'uv.s2_tpxo8_atlas_30c.nc',...
    'uv.k1_tpxo8_atlas_30c.nc',...
    'uv.o1_tpxo8_atlas_30c.nc',...
    'uv.n2_tpxo8_atlas_30c.nc',...
    'uv.p1_tpxo8_atlas_30c.nc',...
    'uv.k2_tpxo8_atlas_30c.nc',...
    'uv.q1_tpxo8_atlas_30c.nc',...
    'uv.m4_tpxo8_atlas_30c.nc'};

fid2={'hf.m2_tpxo8_atlas_30c.nc',...
    'hf.s2_tpxo8_atlas_30c.nc',...
    'hf.k1_tpxo8_atlas_30c.nc',...
    'hf.o1_tpxo8_atlas_30c.nc',...
    'hf.n2_tpxo8_atlas_30c.nc',...
    'hf.p1_tpxo8_atlas_30c.nc',...
    'hf.k2_tpxo8_atlas_30c.nc',...
    'hf.q1_tpxo8_atlas_30c.nc',...
    'hf.m4_tpxo8_atlas_30c.nc'};

fid3='grid_tpxo8atlas_30.nc';


if isempty(conID); conID=1:length(names); end
Nc=length(conID);

j=1;
for i=conID
    % Interpolate U
    lon0=ncread(fid{i},'lon_u');
    lat0=ncread(fid{i},'lat_u');
    ind.x=max([1 find(lon0<min(min(lon)),1,'last')]):min([find(max(max(lon))<lon0,1,'first') length(lon0)]);
    ind.y=find(lat0<min(min(lat)),1,'last'):find(max(max(lat))<lat0,1,'first');
    
    if isempty(ind.x)
        ind.x=1:length(lon0);
    end
    
    U0=double(ncread(fid{i},'uRe',[ind.y(1) ind.x(1)],[length(ind.y) length(ind.x)]))...
        +ii*double(ncread(fid{i},'uIm',[ind.y(1) ind.x(1)],[length(ind.y) length(ind.x)]));
    
    [loni lati]=meshgrid(lon0(ind.x),lat0(ind.y));
    try
        U(:,:,j)=interp2(loni,lati,U0,lon',lat).'/(100*100);
    catch
        U(:,:,j)=interp2(loni,lati,U0,lon,lat)/(100*100);
    end
    
    % Interpolate V
    lon0=ncread(fid{i},'lon_v');
    lat0=ncread(fid{i},'lat_v');
    ind.x=max([1 find(lon0<min(min(lon)),1,'last')]):min([find(max(max(lon))<lon0,1,'first') length(lon0)]);
    ind.y=find(lat0<min(min(lat)),1,'last'):find(max(max(lat))<lat0,1,'first');
      
    V0=double(ncread(fid{i},'vRe',[ind.y(1) ind.x(1)],[length(ind.y) length(ind.x)]))...
        +ii*double(ncread(fid{i},'vIm',[ind.y(1) ind.x(1)],[length(ind.y) length(ind.x)]));
    
    [loni lati]=meshgrid(lon0(ind.x),lat0(ind.y));
    try
        V(:,:,j)=interp2(loni,lati,V0,lon',lat).'/(100*100);
    catch
        V(:,:,j)=interp2(loni,lati,V0,lon,lat)/(100*100);
    end
    
    % Interpolate eta
    lon0=ncread(fid2{i},'lon_z');
    lat0=ncread(fid2{i},'lat_z');
    ind.x=max([1 find(lon0<min(min(lon)),1,'last')]):min([find(max(max(lon))<lon0,1,'first') length(lon0)]);
    ind.y=find(lat0<min(min(lat)),1,'last'):find(max(max(lat))<lat0,1,'first');
   
    eta0=double(ncread(fid2{i},'hRe',[ind.y(1) ind.x(1)],[length(ind.y) length(ind.x)]))...
        +ii*double(ncread(fid2{i},'hIm',[ind.y(1) ind.x(1)],[length(ind.y) length(ind.x)]));
    
    [loni lati]=meshgrid(lon0(ind.x),lat0(ind.y));
    try
        eta(:,:,j)=interp2(loni,lati,eta0,lon',lat).'/1000;
    catch
        eta(:,:,j)=interp2(loni,lati,eta0,lon,lat)/1000;
    end     
    
    % Adjust amplitude and phase to emphemeris
    if day~=0
        if i==9; k=21; else k=i; end
        [weight(j) omega(j)]=TPXO_PHASE(day,k);
        U(:,:,j)=weight(j)*U(:,:,j);
        V(:,:,j)=weight(j)*V(:,:,j);
        eta(:,:,j)=weight(j)*eta(:,:,j);
    end
    j=j+1;
end

% Convert transports to velocities 
if vel_flag
    % Interpolate depth
    lon0=ncread(fid3,'lon_u');
    lat0=ncread(fid3,'lat_u');
    ind.x=find(lon0<min(min(lon)),1,'last'):find(max(max(lon))<lon0,1,'first');
    ind.y=find(lat0<min(min(lat)),1,'last'):find(max(max(lat))<lat0,1,'first');
    
    H=double(ncread(fid3,'hu',[ind.y(1) ind.x(1)],[length(ind.y) length(ind.x)]));
    
    [loni lati]=meshgrid(lon0(ind.x),lat0(ind.y));
    try
        H=interp2(loni,lati,H,lon',lat).';
    catch
        H=interp2(loni,lati,H,lon,lat);
    end
    U=U./repmat(H,[1 1 size(U,3)]);
    
    lon0=ncread(fid3,'lon_v');
    lat0=ncread(fid3,'lat_v');
    ind.x=find(lon0<min(min(lon)),1,'last'):find(max(max(lon))<lon0,1,'first');
    ind.y=find(lat0<min(min(lat)),1,'last'):find(max(max(lat))<lat0,1,'first');
    
    H=double(ncread(fid3,'hv',[ind.y(1) ind.x(1)],[length(ind.y) length(ind.x)]));
    
    [loni lati]=meshgrid(lon0(ind.x),lat0(ind.y));
    try
        H=interp2(loni,lati,H,lon',lat).';
    catch
        H=interp2(loni,lati,H,lon,lat);
    end
    V=V./repmat(H,[1 1 size(U,3)]);
end


return;

function [weight omega names]=TPXO_PHASE(day,conID)
% function [weight omega names]=TPXO_PHASE(day,conID)
%
% Return the weights and frequencies so that TPXO amplitudes can be used 
% to predict tides using 
%
% TS=Real[weight*Amp(x,y)*exp(ii*omega*t)]
%
% where Amp is extracted from TPXO and t is seconds from the start of 
% the input day
%
% INPUT
% day	 	A number or vector in serial days (see 'help datenum')
%
% OUTPUT
% weight    A complex weight that shifts the magnitude and phase of the TPXO amplitude
% omega		The tidal frequency in rad/s
% names		The names of the constituents
%
% written by Sam Kelly 5-1-2014 (smkelly@d.umn.edu)
% 
% based on Laurie Padman's TMD code, the Oregon State Tidal prediction 
% software (OTPS), Richard Ray's notes, and David Cartwright's 
% calculations

% Constituent names (Everything is in this order)
names={'m2  ';'s2  ';'k1  ';'o1  '; ...
    'n2  ';'p1  ';'k2  ';'q1  '; ...
    '2n2 ';'mu2 ';'nu2 ';'l2  '; ...
    't2  ';'j1  ';'no1 ';'oo1 '; ...
    'rho1';'mf  ';'mm  ';'ssa ';'m4  '};

% Frequencies
omega=[1.405189e-04;1.454441e-04;7.292117e-05;6.759774e-05; ...
    1.378797e-04;7.252295e-05;1.458423e-04;6.495854e-05; ...
    1.352405e-04;1.355937e-04;1.382329e-04;1.431581e-04; ...
    1.452450e-04;7.556036e-05;7.028195e-05;7.824458e-05; ...
    6.531174e-05;0.053234e-04;0.026392e-04;0.003982e-04;2.810377e-04]';

% Astronomical phase (relative to t0 = 1 Jan 0:00 1992]
phase=[ 1.731557546;0.000000000;0.173003674;1.558553872;...
    6.050721243;6.110181633;3.487600001;5.877717569;
    4.086699633;3.463115091;5.427136701;0.553986502;
    0.052841931;2.137025284;2.436575100;1.929046130;
    5.254133027;1.756042456;1.964021610;3.487600001;
    3.463115091]';

% Nodal corrections f and u
omegaN=125.0445-0.05295377*(day-51544.4993);% mean longitude of ascending lunar node
p=83.3535+0.11140353*(day-51544.4993);% mean longitude of lunar perigee
sinn = sin(omegaN*pi/180);
cosn = cos(omegaN*pi/180);
sin2n = sin(2*omegaN*pi/180);
cos2n = cos(2*omegaN*pi/180);
sin3n = sin(3*omegaN*pi/180);

f(1)=sqrt((1.-.03731*cosn+.00052*cos2n).^2+(.03731*sinn-.00052*sin2n).^2); 
f(2)=1;
f(3)=sqrt((1.+.1158*cosn-.0029*cos2n).^2+(.1554*sinn-.0029*sin2n).^2); 
f(4)=sqrt((1.0+0.189*cosn-0.0058*cos2n).^2+(0.189*sinn-0.0058*sin2n).^2);
f(5)=sqrt((1.-.03731*cosn+.00052*cos2n).^2+(.03731*sinn-.00052*sin2n).^2);
f(6)=1;
f(7)=sqrt((1.+.2852*cosn+.0324*cos2n).^2+(.3108*sinn+.0324*sin2n).^2); 
f(8)=sqrt((1.+.188*cosn).^2+(.188*sinn).^2);
f(9)=sqrt((1.-.03731*cosn+.00052*cos2n).^2+(.03731*sinn-.00052*sin2n).^2);
f(10)=sqrt((1.-.03731*cosn+.00052*cos2n).^2+(.03731*sinn-.00052*sin2n).^2);
f(11)=sqrt((1.-.03731*cosn+.00052*cos2n).^2+(.03731*sinn-.00052*sin2n).^2);
temp1=1.-0.25*cos(2*p*pi/180)-0.11*cos((2*p-omegaN)*pi/180)-0.04*cosn;
temp2=0.25*sin(2*p*pi/180)+0.11*sin((2*p-omegaN)*pi/180)+ 0.04*sinn;
f(12)=sqrt(temp1.^2 + temp2.^2);             
f(13)=1;
f(14)=sqrt((1.+.169*cosn).^2+(.227*sinn).^2);
tmp1=1.36*cos(p*pi/180)+.267*cos((p-omegaN)*pi/180);% Ray's
tmp2=0.64*sin(p*pi/180)+.135*sin((p-omegaN)*pi/180);
f(15)=sqrt(tmp1.^2 + tmp2.^2);                
f(16)=sqrt((1.0+0.640*cosn+0.134*cos2n).^2+(0.640*sinn+0.134*sin2n).^2 );
f(17)=sqrt((1.+.188*cosn).^2+(.188*sinn).^2);
f(18)=1.043 + 0.414*cosn;          
f(19)=1-0.130*cosn;           
f(20)=1;
f(21)=(1.-.03731*cosn+.00052*cos2n).^2+(.03731*sinn-.00052*sin2n).^2;

u(1)=atan((-.03731*sinn+.00052*sin2n)./(1.-.03731*cosn+.00052*cos2n))/pi/180;                                   
u(2)=0;
u(3)=atan((-.1554*sinn+.0029*sin2n)./(1.+.1158*cosn-.0029*cos2n))/pi/180;
u(4)=10.8*sinn-1.3*sin2n+0.2*sin3n;
u(5)=atan((-.03731*sinn+.00052*sin2n)./(1.-.03731*cosn+.00052*cos2n))/pi/180;
u(6)=0;
u(7)=atan(-(.3108*sinn+.0324*sin2n)./(1.+.2852*cosn+.0324*cos2n))/pi/180;
u(8)=atan(.189*sinn./(1.+.189*cosn))/pi/180;
u(9)=atan((-.03731*sinn+.00052*sin2n)./(1.-.03731*cosn+.00052*cos2n))/pi/180;
u(10)=atan((-.03731*sinn+.00052*sin2n)./(1.-.03731*cosn+.00052*cos2n))/pi/180;
u(11)=atan((-.03731*sinn+.00052*sin2n)./(1.-.03731*cosn+.00052*cos2n))/pi/180;
u(12) = atan(-temp2./temp1)/pi/180 ;
u(13) = 0;
u(14) = atan(-.227*sinn./(1.+.169*cosn))/pi/180; 
u(15) = atan2(tmp2,tmp1)/pi/180;
u(16) = atan(-(.640*sinn+.134*sin2n)./(1.+.640*cosn+.134*cos2n))/pi/180;  
u(17) = atan(.189*sinn./(1.+.189*cosn))/pi/180; 
u(18) = -23.7*sinn + 2.7*sin2n - 0.4*sin3n;
u(19) = 0;
u(20) = 0;
u(21) = (atan((-.03731*sinn+.00052*sin2n)./ ...
           (1.-.03731*cosn+.00052*cos2n))/pi/180)*2;
       
% Make the weight
ii=complex(0,1);

omega=omega(conID);
weight=f(conID).*exp(ii*(phase(conID)+u(conID)*pi/180+omega*(day-datenum(1992,1,1))*86400));

return;











