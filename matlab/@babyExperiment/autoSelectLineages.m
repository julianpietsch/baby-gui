function autoSelectLineages(cExperiment,poses)
if nargin<2 || isempty(poses)
    poses = find(cExperiment.posTracked);
end

if size(cExperiment.cellsToPlot,3)>1
    cExperiment.cellsToPlot=cell(1);
    for i=1:length(cExperiment.posTracked)
        cExperiment.cellsToPlot{i}=sparse(zeros(1,1));
    end
end

% Start logging protocol

cExperiment.logger.start_protocol('auto-selecting any cells in a lineage',length(poses));
try
    
    for i=1:length(poses)
        pos=poses(i);
        cTimelapse=cExperiment.returnTimelapse(pos);
        cTimelapse.autoSelectLineages;
        cExperiment.cTimelapse=cTimelapse;
        cExperiment.cellsToPlot{i}=cTimelapse.cellsToPlot;
        if i==length(poses)
            cExperiment.saveTimelapseExperiment(pos);
        else
            cExperiment.saveTimelapse(pos);
        end
    end

% Finish logging protocol
cExperiment.logger.complete_protocol;
catch err
    cExperiment.logger.protocol_error;
    rethrow(err);
end
end