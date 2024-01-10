function babySegment(cExperiment,varargin)

ip = inputParser;
isnat = @(x) all(isnumeric(x)) && all(round(x)==x) && all(x>=0);
ip.addOptional('poses',[],@(x) isempty(x) || (isvector(x) && isnat(x)));
ip.addParameter('copyconfig',true,@(x) islogical(x) && isscalar(x));
ip.KeepUnmatched = true;
ip.parse(varargin{:});

poses = ip.Results.poses;
copyconfig = ip.Results.copyconfig;
unmatched = [fieldnames(ip.Unmatched),struct2cell(ip.Unmatched)]';

if isempty(poses)
    poses = 1:length(cExperiment.dirs);
end

if copyconfig
    assert(~strcmp(cExperiment.babyBrain.status,'asleep'),...
        'A valid BabyBrain configuration needs to be specified');
end

%% Start logging protocol
cExperiment.logger.start_protocol('baby segmentation',length(poses));
try

%% Run segmentation on each timelapse
for pos=poses
    cTimelapse = cExperiment.loadCurrentTimelapse(pos);
    if copyconfig
        cTimelapse.babyBrain = cExperiment.babyBrain;
    end
    cTimelapse.babySegment(unmatched{:});
    cExperiment.saveTimelapseExperiment([],false);
    cExperiment.posSegmented(pos) = true;
    cExperiment.posTracked(pos) = true; % BABY segments and tracks
end

cExperiment.saveExperiment;

%% Finish logging protocol
cExperiment.logger.complete_protocol;
catch err
    cExperiment.logger.protocol_error;
    rethrow(err);
end

end
