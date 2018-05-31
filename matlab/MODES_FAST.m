function [PHI C]=MODES_FAST(dz,N2,Nm,Nm0)
% USAGE: [PHI C]=MODES_FAST(dz,N2,Nm,Nm0)
% Obtain vertical modes for arbitrary stratification with (or without) a free surface
%
% INPUTS:
% dz  [1  x 1]   vertical spacing (must be uniform)
% N2  [Nz x 1]   stratification 
% Nm  [1  x 1]   number of internal modes to extract 
% Nm0 [1  x 1]   number of intenral modes to use as the basis for the spectral method 
%
% OUTPUTS:
% PHI [Nz x Nm]  pressure and velocity structure eigenfunctions 
% C   [Nm x 1 ]  eigenspeeds 
%
% A good reference for spectral and pseudo-spectral methods is:
% J. P. Boyd (2001) Chebyshev and Fourier Spectral Methods: Second Revised Edition, Dover, 688p.
%
% Sam Kelly, May 2018 (smkelly@d.umn.edu)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure out how N2 couples the modes (using analytical solutions to integrals)
N=length(N2); % Number of vertical grid points

% Find cosine series coefficients 
N2s=[N2(N:-1:1); N2(1:end-1)];
xN=real(fft(N2s))/N;
xN=xN(2:2*Nm0+1).*(-1).^(1:2*Nm0);
x0=mean(N2);

% Reconstruct N2 to check that we have the cosine coefficients
%N2r=xN(1)+cos([0.5:N]'*[1:Nm0]*pi/N)*xN(2:end)';
%N2r=cos([0.5:N]'*[0:Nm0]*pi/N)*xN';
%plot(N2,z,'k',N2r,z,'r--');

% Compute the coupling coefficients using analytical integrals
A=zeros([Nm0 Nm0]);
for m=1:Nm0
    for n=1:m-1
        % The only nonzer integrals are when k=m-n and k=m+n (note: -(m+n)>0 and can be igored. Also, if m-n>0 then n-m<0 and can be ignored)
        A(m,n)=(xN(m-n)-xN(m+n))/(2*m*pi*n*pi);
    end
    % Find diagonal terms and (divide by 2, since they are doubled later)
    A(m,m)=(2*x0-xN(2*m))/(4*m*pi*m*pi);
end
A=(dz*N)*(dz*N)*(A+A');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve eigenvalue problem ("x" are the expansion coefficients for the spectral basis functions)
[x,C]=eig(A,'nobalance','vector');

% Sort modes by eigenspeed and retain Nm fastest modes
C(C<0)=0;
C=sqrt(C); % eigenspeed
[C,ind]=sort(C,1,'descend');
C=C(1:Nm);
x=x(:,ind(1:Nm));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruct the p modes using expansion coefficients
x0=zeros([2*N Nm]);
x0(2:Nm0+1,:)=x;
PHI=sqrt(2)*real(ifft(x0,'symmetric'))*N;
PHI=PHI(1:N,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the p modes as positive at the surface
PHI(:,PHI(1,:)<0)=-PHI(:,PHI(1,:)<0);


return