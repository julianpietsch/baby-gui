function cTimelapse=returnTimelapse(cExperiment,timelapseNum)
%cTimelapse=returnTimelapse(cExperiment,timelapseNum)
%
% loads a cTimelapse. timelapseNum should be a single number
% indicating which position to load. Note, timelapseNum indicates an index
% in cExperiment.dirs to load, so depending on the ordering of the
% directories in dirs cExperiment.loadCurrentTimelapse(2) will not
% necessarily load the cTimlapse related to directory pos2, and will in
% general load pos10 - this is due to alphabetic ordering.
%
% the distinction between this method, and
% BABYEXPERIMENT.LOADCURRENTTIMELAPSE is that this one does not
% populate the cTimelapse and currentPos fields of the babyExperiment
% object, and as such should just be used to return a cTimelapse.
%
% See also, BABYEXPERIMENT.LOADCURRENTTIMELAPSE

l1 = load([cExperiment.saveFolder,filesep,cExperiment.dirs{timelapseNum},'cTimelapse']);
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
