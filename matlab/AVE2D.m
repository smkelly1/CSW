function [out]=AVE2D(in,N,win)
% AVE2D smooths a 2D data array 
%
% AVE2D(in,N,win) returns symmetrically smoothed data using a 2D 
% convolution with the window win  
%
% INPUTS
%     in - data series profile [1xM]
%     N - number of scans/points over which to average, Default 4
% OPTIONAL
%     win - choose the window type from matlab windows:
%          Hanning - default
%          Bartlett, Blackman, Boxcar, Hamming, etc
% 
% Jan-06 (skelly@coas.oregonstate.edu)
if nargin<3
    win='hanning';    
end
eval(['win=',win,'(N);']);

weight=win*win.';
weight=weight/(sum(sum(weight)));

m = ~isnan(in);
n = ~isnan(weight);

in(~m) = 0;
out = conv2(in,weight,'same'); % 2-dim convolution
N = conv2(real(m),real(n),'same'); % normalization term
c = conv2(ones(size(in)),ones(size(weight)),'same'); % correction of normalization

out=out.*c./N;

return
