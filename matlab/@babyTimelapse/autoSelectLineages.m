function autoSelectLineages(cTimelapse,traps)
if nargin<2 || isempty(traps)
    traps = 1:size(cTimelapse.cellMothers,1);
end

% Ensure that cellsToPlot and cellMothers have the same label size
nlbl_toplot = size(cTimelapse.cellsToPlot,2);
nlbl_mother = size(cTimelapse.cellMothers,2);
if nlbl_toplot < nlbl_mother
    cTimelapse.cellsToPlot(:,end+1:nlbl_mother) = 0;
elseif nlbl_mother < nlbl_toplot
    cTimelapse.cellMothers(:,end+1:nlbl_toplot) = 0;
end

% First select all daughter cells
cTimelapse.cellsToPlot(traps,:) = cTimelapse.cellMothers(traps,:) > 0;
% Then select all mothers
for t=traps
    mothers = cTimelapse.cellMothers(t,:);
    mothers = unique(mothers(mothers>0));
    cTimelapse.cellsToPlot(t,mothers) = true;
end
end