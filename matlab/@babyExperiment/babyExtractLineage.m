function babyExtractLineage(cExperiment,poses,assign,mintps,calcEllipse)

if nargin<2 || isempty(poses)
    poses = 1:length(cExperiment.dirs);
end
if nargin<3, assign = true; end
if nargin<4, mintps = []; end
if nargin<5, calcEllipse = []; end

%% Start logging protocol
cExperiment.logger.start_protocol('baby lineage extraction',length(poses));
try

%% Run segmentation on each timelapse
for pos=poses
    cTimelapse = cExperiment.loadCurrentTimelapse(pos);
    if assign
        cTimelapse.assignMothers([],mintps);
    end
    cTimelapse.extractMotherDaughterArea(calcEllipse);
    cExperiment.saveTimelapseExperiment([],false);
end

%% Finish logging protocol
cExperiment.logger.complete_protocol;
catch err
    cExperiment.logger.protocol_error;
    rethrow(err);
end

end