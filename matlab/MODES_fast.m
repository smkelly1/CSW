function [PHI C]=MODES_fast(dz,N2,Nm,Nm0)
% USAGE: [PHI C]=MODES_fast(dz,N2,Nm,Nm0)
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
% Sam Kelly, May 2017 (smkelly@d.umn.edu)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define spectral basis functions (i.e., rigid-lid functions in constant stratification)
N=length(N2);      % Number of vertical grid points
H=dz*N;            % Water depth

alpha0=(1:Nm0)*pi; % Define constant-N modes (alpha0_n = n*pi only for a rigid lid)
z=(1:N)*dz-dz/2;   % The z grid is define in the center of the dz sections
phi_w=sqrt(2)*sin((1-z'/H)*alpha0)*sparse(1:Nm0,1:Nm0,1./alpha0)*H;
phi_p=sqrt(2)*cos((1-z'/H)*alpha0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define A matrix (depth-varying stratification couples the constant N2 modes); Note: Although the typical finite-difference matrix is sparse (banded), this matrix of coupling coefficients is typically dense.
A=phi_w'*sparse(1:N,1:N,N2)*phi_w/N;

% Solve eigenvalue problem ("x" are the expansion coefficients for the spectral basis functions)
[x,C]=eig(A);

% Sort modes by eigenspeed and retain Nm fastest modes
C=diag(C);
C(C<0)=0;
C=sqrt(C); % eigenspeed
[C,ind]=sort(C,1,'descend');
C=C(1:Nm);
x=x(:,ind(1:Nm));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruct the p modes using expansion coefficients
PHI=phi_p*x;

% Normalize the p modes (Doesn't appear to be necessary...)
%A=sqrt(sum(PHI.^2)/N);
%PHI=PHI.*repmat(1./A,[N 1]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the p modes as positive at the surface
PHI(:,PHI(1,:)<0)=-PHI(:,PHI(1,:)<0);


return



