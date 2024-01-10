function is_tracked = isTracked(cTimelapse,timepoints)
% is_tracked = isTracked(cTimelapse,timepoints)
% 
% checks if cTimelapse has been tracked upto and including hte timepoints
% 'timepoints'.
if nargin<2 
    timepoints = cTimelapse.timepointsToProcess;
end

is_tracked = all(~cellfun('isempty',{cTimelapse.cTimepoint(timepoints).trapLocations}));

end