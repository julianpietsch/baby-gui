function cTimelapse=returnTimelapse(cExperiment,timelapseNum)
%cTimelapse=returnTimelapse(cExperiment,timelapseNum)
%
% loads a cTimelapse. timelapseNum should be a single number
% indicating which position in cExperiment.dirs to load.
%
% See also, EXPERIMENTTRACKING.LOADCURRENTTIMELAPSE

cExperiment.dataset.ensure_session;
posName = cExperiment.dirs{timelapseNum};
fileName = [posName 'cTimelapse_' cExperiment.rootFolder '.mat'];

% Download/update the specified cTimelapse if necessary:
localFile = cExperiment.server.downloadFile(cExperiment.dataset.dataset,fileName);

% Load the cTimelapse
l1 = load(localFile);
cTimelapse = l1.cTimelapse;

was_disco = isprop(cTimelapse,'recovered_data') && ...
    exist('disco_baby_map','file') == 2;
if was_disco
    % Partially backwards compatible map for old DISCO-style cTimelapses. 
    % This can only be called if the 'baby-gui-migration' repo is on the 
    % path and only saves trap template data from cellVision, so there will
    % be data loss for DISCO experiments.
    new_class_name = disco_baby_map(class(cTimelapse));
    recovered_data = cTimelapse.recovered_data;
    assert(startsWith(new_class_name,'babyTimelapse') ...
        && exist(new_class_name,'class')==8,'bad load_class_name!');
    cTimelapse = eval([new_class_name,'([],true)']);
    cTimelapse.copyprops(recovered_data);
    cTimelapse.trapTemplates = cExperiment.trapTemplates;
end

assert(isequal(posName,cTimelapse.dataset.pos),...
    'saved cTimelapse has a different position name than stored in cExperiment');

%Ensure the timelapse has the channels list for the dataset
cTimelapse.channelNames = cExperiment.channelNames;
cTimelapse.metadata = cExperiment.metadata;
cTimelapse.metadata.posname = cExperiment.dirs{timelapseNum};
cTimelapse.metadata.exptid = cExperiment.id;

% In either case, once the timelapse is successfully loaded, trigger a
% PositionChanged event to notify babyLogging
babyLogging.changePos(cExperiment,timelapseNum,cTimelapse);

% populate these transient properties for when new cells need to be
% detected.
cTimelapse.imcache = cExperiment.imcache;
end
