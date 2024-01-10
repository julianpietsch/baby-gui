function max_cell_label = returnMaxCellLabel(cTimelapse,traps,timepoint_range)
% max_cell_num = getMaxCellNum(cTimelapse,traps,timepoint_range)
%
% intended to replace trapMaxCell, which is ever problematic. returns the
% maximum cell label occurring in the timepoints timepoint_range. defaults
% to cTimelapse.timepointsToProcess.

if nargin<2 || isempty(traps)
    traps = cTimelapse.defaultTrapIndices;
end

if nargin<3 || isempty(timepoint_range)
    
    timepoint_range = cTimelapse.timepointsToProcess;
    
end

max_cell_label = zeros(size(traps));
for ti = 1:length(traps)
    for tp = timepoint_range
        
        max_cell_label(ti) = max([cTimelapse.cTimepoint(tp).trapInfo(traps(ti)).cellLabel max_cell_label(ti)]);
        
    end
end
