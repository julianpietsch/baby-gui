function selectCellsToPlotAutomatic(cExperiment,positionsToCheck,params)
if nargin<2
    positionsToCheck=find(cExperiment.posTracked);
    %     positionsToTrack=1:length(cExperiment.dirs);
end

if nargin<3 || isempty(params)
    params = cExperiment.cellAutoSelectParams;
    if isempty(params)
        params = cExperiment.selectCellSelectParamsGUI();
    end
    if isfield(params,'anyInLineage') && params.anyInLineage
        cExperiment.autoSelectLineages(positionsToCheck);
        return
    end
end

if size(cExperiment.cellsToPlot,3)>1
    cExperiment.cellsToPlot=cell(1);
    for i=1:length(cExperiment.posTracked)
        cExperiment.cellsToPlot{i}=sparse(zeros(1,1));
    end
end

%for backcompatibility with scripts that don't use these parameters
if ~isfield(params,'maximumNumberOfCells')
    params.maximumNumberOfCells = Inf;
end
if ~isfield(params,'includeDaughters')
    params.includeDaughters = false;
end

%% Run the tracking on the timelapse

% Start logging protocol

cExperiment.logger.add_arg('Fraction of timelapse that cells are present for',params.fraction);
cExperiment.logger.add_arg('Number of frames a cell must be present',params.duration);
cExperiment.logger.add_arg('Cell must appear by frame',params.framesToCheck);
cExperiment.logger.add_arg('Cell must still be present by frame',params.framesToCheckEnd);
cExperiment.logger.add_arg('Maximum number of cells',params.maximumNumberOfCells);
cExperiment.logger.add_arg('Include daughters?',params.includeDaughters);
cExperiment.logger.start_protocol('autoselecting cells',length(positionsToCheck));

try

    for i=1:length(positionsToCheck)
        if params.maximumNumberOfCells>=0
            experimentPos=positionsToCheck(i);
            cTimelapse=cExperiment.returnTimelapse(experimentPos);
            cTimelapse.automaticSelectCells(params);
            params.maximumNumberOfCells = max(params.maximumNumberOfCells - full(sum(cTimelapse.cellsToPlot(:))),0);
            cExperiment.cTimelapse=cTimelapse;
            cExperiment.cellsToPlot{i}=cTimelapse.cellsToPlot;
            %         cExperiment.saveTimelapseExperiment(experimentPos);
            if i==length(positionsToCheck)
                cExperiment.saveTimelapseExperiment(experimentPos);
            else
                cExperiment.saveTimelapse(experimentPos);
            end
        end
    end

% Finish logging protocol
cExperiment.logger.complete_protocol;
catch err
    cExperiment.logger.protocol_error;
    rethrow(err);
end

end
