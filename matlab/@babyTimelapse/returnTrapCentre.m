function [ trap_centre ] = returnTrapCentre( cTimelapse,trap_indices,timepoint )
% [ trap_centre ] = returnTrapCentre( cTimelapse,trap_indices,timepoint )

% can handle an TrapIndex as an array, in which case returns
% column of form [x's   y's].
% If the image does not contain traps it returns [0 0] for every
% trap requested.
% if trap_index is empty it returns all of them.

if isempty(trap_indices)
    
    trap_indices = cTimelapse.defaultTrapIndices;
    
end

if cTimelapse.trapsPresent
    trap_centre = round([[cTimelapse.cTimepoint(timepoint).trapLocations(trap_indices).xcenter]' [cTimelapse.cTimepoint(timepoint).trapLocations(trap_indices).ycenter]']);
    
else
    trap_centre = repmat([0 0],length(trap_indices),1);
end

end

