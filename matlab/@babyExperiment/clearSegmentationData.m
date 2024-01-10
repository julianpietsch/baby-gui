function clearSegmentationData(cExperiment,poses)
%CLEAR_SEGMENTATION_DATA deletes outlines in all timelapses
%
%   Loops over cTimelapses to reset trapInfo structs, lineage assignments and
%   selected cells (cellsToPlot)
if nargin<2 || isempty(poses)
    poses = 1:numel(cExperiment.dirs);
end

for p=poses
    cTimelapse = cExperiment.loadCurrentTimelapse(p);

    % Clear lineage and selected cell info
    cTimelapse.cellMothers(:) = 0;
    cTimelapse.cellsToPlot(:) = false;

    % Reset trapInfo struct arrays
    tps = cTimelapse.timepointsToProcess;
    reftp = tps(1);
    tmplt = cTimelapse.createTrapInfoTemplate;
    tmplt = tmplt(ones(size(cTimelapse.cTimepoint(reftp).trapLocations)));
    [cTimelapse.cTimepoint(tps).trapInfo] = deal(tmplt);

    cExperiment.saveTimelapseExperiment([],false);
end

cExperiment.posSegmented(:) = false;
cExperiment.posTracked(:) = false;
cExperiment.saveExperiment;

end
