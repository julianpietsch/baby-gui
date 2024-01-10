function TrapsToRemove = identifyExcludedTraps(cDisplay,trapLocations,trapLocationsNotToRemove )
% TrapsToRemove = identifyExcludedTraps(cDisplay,trapLocations )
%
% identify which traps in trapLocation are in the exclude zones of
% cDisplay.
% trapLocations defaults to cDisplay.trapLocations and is a structure array
% with the fields xcenter, ycenter (as in babyTimelapse.trapLocations).
%
% TrapsToRemove is a  list of indices of the trapLocations (i.e. entries in
% the structure array) that occur

if nargin<2 || isempty(trapLocations)
    trapLocations = cDisplay.trapLocations;
end

if nargin<3
    trapLocationsNotToRemove = [];
end
    
% bit ugly, but don't want to remove traps that are in the
% trapLocationsNotToRemove (these are typically traps added by a user on
% another instance of the same GUI.

if ~isempty(trapLocationsNotToRemove) 
    [~,TrapsNotToRemove] = ismember([[trapLocations(:).xcenter]' [trapLocations(:).ycenter]'],...
        [[trapLocationsNotToRemove(:).xcenter]' [trapLocationsNotToRemove(:).ycenter]'],...
        'rows');
    
    %only remove traps marked to remove and not to put back.
    trapLocations = trapLocations(~TrapsNotToRemove);

end

TrapsToRemove = [];
if~isempty(trapLocations) && ~isempty(trapLocations(1).xcenter)
    for trapi = 1:length(trapLocations)
        for zonei = 1:size(cDisplay.ExclusionZones,1)
            
            if trapLocations(trapi).xcenter>=cDisplay.ExclusionZones(zonei,1) && ...
                    trapLocations(trapi).xcenter<=cDisplay.ExclusionZones(zonei,1) + cDisplay.ExclusionZones(zonei,3) && ...
                    trapLocations(trapi).ycenter>=cDisplay.ExclusionZones(zonei,2) && ...
                    trapLocations(trapi).ycenter<=cDisplay.ExclusionZones(zonei,2) + cDisplay.ExclusionZones(zonei,4);
                TrapsToRemove = [TrapsToRemove trapi];
            end
        end
        
    end
end



end

