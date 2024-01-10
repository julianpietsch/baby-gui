function mmRetrackCells(cExperiment,poses,varargin)

if nargin<2 || isempty(poses)
    poses = 1:length(cExperiment.dirs);
end

%% Start logging protocol
cExperiment.logger.start_protocol('MM retracking',length(poses));
try

%% Run segmentation on each timelapse
cExperiment.posTracked(poses) = false;
for pos=poses
    cTimelapse = cExperiment.loadCurrentTimelapse(pos);
    cTimelapse.mmRetrackCells(varargin{:});
    cExperiment.saveTimelapseExperiment([],false);
    cExperiment.posTracked(pos) = true;
end

cExperiment.saveExperiment;

%% Finish logging protocol
cExperiment.logger.complete_protocol;
catch err
    cExperiment.logger.protocol_error;
    rethrow(err);
end

end