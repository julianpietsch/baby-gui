function clearTrapInfo( cTimelapse )
% clearTrapInfo( cTimelapse ) 
%
% fairly needless function to clear away the trapInfo completely from
% the cTimelapse and return it to an unblemished state.

[cTimelapse.cTimepoint(:).trapInfo] = deal(cTimelapse.trapInfoTemplate);
[cTimelapse.cTimepoint(:).trapLocations] = deal([]);
cTimelapse.timepointsProcessed = false(1,max(cTimelapse.timepointsToProcess));

end

