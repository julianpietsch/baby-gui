function trackCells(cExperiment,positionsToTrack,cellMovementThresh)
%trackCells(cExperiment,positionsToTrack,cellMovementThresh)
%-------------------------------------------------------------

if isempty(cExperiment.saveFolder)
    cExperiment.saveFolder=cExperiment.rootFolder;
end


if nargin<2
    positionsToTrack=find(cExperiment.posSegmented);
%     positionsToTrack=1:length(cExperiment.dirs);
end

if nargin<3
    prompt = {['Max change in position and radius before a cell is classified as a new cell. A larger number is more lenient, and you will be more likely' ...
        ' to not have interruptions in the tracking for a cell. At the same time though, you will be more likely to identify unrelated cells as the same.' ...
        ' This is especially true for daughters.']};
    dlg_title = 'Tracking Threshold';
    num_lines = 1;
    def = {'5'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    cellMovementThresh=str2double(answer{1});
end

%% Start logging protocol
cExperiment.logger.add_arg('Tracking threshold',cellMovementThresh);
cExperiment.logger.start_protocol('tracking cells',length(positionsToTrack));
try

%% Run the tracking on the timelapse
for i=1:length(positionsToTrack)
    experimentPos=positionsToTrack(i);
    cTimelapse=loadCurrentTimelapse(cExperiment,experimentPos);
    cTimelapse.trackCells(cellMovementThresh);
    cExperiment.posTracked(experimentPos)=1;
    cExperiment.cTimelapse=cTimelapse;
    if i==length(positionsToTrack)
        cExperiment.saveTimelapseExperiment(experimentPos);
    else
        cExperiment.saveTimelapse(experimentPos);
    end
end

%% Finish logging protocol
cExperiment.logger.complete_protocol;
catch err
    cExperiment.logger.protocol_error;
    rethrow(err);
end

end
