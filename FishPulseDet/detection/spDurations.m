function [start, stop] = spDurations(Indices, MergeThreshold,idxMax)
% [start, duration] = spDurations(Indices)
% Given a vector of indices into another vector, determine
% the starting point of each distinct region and how long it is.
%
% example:
% [start, stop, dur] = spDurations([17:20, 50:52, 75:80])
% returns
% start = [17, 50, 75]
% stop = [20, 52, 80]
% duration = [4, 3, 6]
%Indices2 = Indices(Indices>5 & Indices<(idxMax-5));
if isempty(Indices)
    stop = [];
    start = [];
    return
end
% find skips in sequence
diffs = diff(Indices)';

% 1st index is always a start 
% last index is always a stop
% indices whose first difference is greater than one denote a 
% start or stop boundary.
%start = Indices([1; find(diffs > 1) + 1]);
startPositions = [1; find(diffs > MergeThreshold) + 1];
start = Indices(startPositions);
if length(startPositions) > 1
    stopPositions = [startPositions(2:end) - 1; length(Indices)];
    stop = Indices(stopPositions);
else
    stop = min(Indices(startPositions)+1,idxMax);
end



