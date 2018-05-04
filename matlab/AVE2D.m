function [series_s]=AVE2D(series,N,win)
% AVE2D smooths data signal from a SeaBird CTD, unbinned
%
% AVE2D(SERIES,N,[mywindow]) returns a symmetrically smoothed data 
% series signal with the same size as PRESS_SMO.  Smooths using a MYWINDOW 
% with N elements (over N scans or points of CTD data).  Seasoft recommends 
% a boxcar with at least N=4 for a pressure signal before using their Loop 
% Edit.
%
% INPUTS
%     SERIES data series profile [1xM]
%     N number of scans/points over which to average, Default 4
% OPTIONAL
%     MYWINDOW - choose the window type from matlab windows:
%          Hanning - default
%          Bartlett, Blackman, Boxcar, Hamming, etc
% 
% Jan-06 (skelly@coas.oregonstate.edu)
if nargin<3
    win='hanning';
    if nargin<2
        N=4;
    end
end
eval(['win=',win,'(N);']);
%c=win/sum(win); % so sums to 1;

c=win*win.';
c=c/(sum(sum(c)));
%c=ones(N)/N;

flipflag=0;  % To make sure series is a column vector
if size(series,2)~=1 
    flipflag=1;
    series = series.'; 
end 
series_s=conv2(series,c,'same');

if flipflag
    series_s = series_s.';
end

series_s(1:ceil(N/2)) = series(1:ceil(N/2));
leni=length(series_s);
series_s(leni-[0:ceil(N/2)]) = series(leni-[0:ceil(N/2)]);

return