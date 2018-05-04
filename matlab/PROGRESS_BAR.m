function PROGRESS_BAR(ind,all_inds)
% PROGRESS_BAR(ind,all_inds)
% Make a progress bar 
%
% Example call:
% for i=1:50 
%   PROGRESS_BAR(i,1:50)
% end
%
% 
global counters

% Initiate
if ind==all_inds(1)    
    disp(['0%       10%       20%       30%       40%       50%       60%       70%       80%       90%      100%'])
    counters.percent=0; 
    counters.N=1;
end

% Print progress
temp=round(counters.N/length(all_inds)*100);
while counters.percent<=temp
    fprintf('|');
    counters.percent=counters.percent+1;    
end
counters.N=counters.N+1;

% Closeout
if ind==all_inds(end)
    clear global counters
    fprintf('\n'); 
end

return
