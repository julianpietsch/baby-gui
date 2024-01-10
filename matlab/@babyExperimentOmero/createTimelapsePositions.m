function createTimelapsePositions(cExperiment,searchString,positionsToLoad,pixelSize,image_rotation,timepointsToLoad,traps_present)
% createTimelapsePositions(cExperiment,searchString,positionsToLoad,pixelSize,image_rotation,timepointsToLoad,traps_present)
%
% goes through the folders in the rootDir of cExperiment and creates a
% babyTimelapse object for each one, instantiating the cTimepoint
% structure array using the files containing the string searchString (see
% babyTimelapse.loadTimelapse for details). Any input not provided is
% defined by GUI.
%
% inputs are all those passed to loadTimelapse method of babyTimelapse in
% creating each timelapse.
%
% See also BABYTIMELAPSE.LOADTIMELAPSE

dataset = cExperiment.dataset;

chNames = cExperiment.channelNames;
if isempty(searchString)
    ch = menu('Choose channel used in segmentation (brightfield/DIC images)',chNames);
    searchString=chNames{ch};
end

if ischar(positionsToLoad)&&strcmp(positionsToLoad,'all')
    positionsToLoad=1:numel(dataset.posNames);
end

assert(all(ismember(cExperiment.dirs,dataset.posNames)),...
    'cExperiment.dirs specifies positions that do not exist in dataset');

if nargin<4
    pixelSize=[];
end

if nargin<5
    image_rotation=[];
end

%timepoints to load functionality has been superceded by
%timepointstoProcess, which is done after the object has been created and
%just limits processing.
if nargin<6
    timepointsToLoad=[];
end

if nargin<7
    traps_present = cExperiment.trapsPresent;
end



cExperiment.searchString=searchString;
cExperiment.pixelSize=pixelSize;
cExperiment.image_rotation=image_rotation;
cExperiment.timepointsToLoad=timepointsToLoad;
cExperiment.trapsPresent = traps_present;
% Start adding arguments to experiment creation protocol log:

cExperiment.logger.add_arg('Omero experiment name',cExperiment.rootFolder);
cExperiment.logger.add_arg('Temporary working directory',cExperiment.saveFolder);

if isempty(cExperiment.timepointsToLoad)
    cExperiment.logger.add_arg('Timepoints to load','all');
else
    cExperiment.logger.add_arg('Timepoints to load',cExperiment.timepointsToLoad);
end
if ~isempty(cExperiment.pixelSize)
    
end
% The other arguments are added and the protocol started after the first
% call to loadTimelapse below...

try
       
    % Load timelapses
    for i=1:length(positionsToLoad)
        currentPos=positionsToLoad(i);
        dataset.pos = cExperiment.dirs{currentPos};
        cExperiment.cTimelapse = babyTimelapseOmero(dataset,cExperiment.channelNames);
        
        % Trigger a PositionChanged event to notify babyLogging
        babyLogging.changePos(cExperiment,currentPos,cExperiment.cTimelapse);
        cExperiment.cTimelapse.metadata = cExperiment.metadata;
        cExperiment.cTimelapse.metadata.posname = cExperiment.dirs{currentPos};
        cExperiment.cTimelapse.trapTemplates = cExperiment.trapTemplates;
        cExperiment.cTimelapse.loadTimelapse(cExperiment.searchString,cExperiment.image_rotation,cExperiment.trapsPresent,cExperiment.timepointsToLoad,cExperiment.pixelSize);
        
        % After the first call to loadTimelapse, the arguments should now all
        % be set, so start logging the creation protocol:
        if i==1
            cExperiment.pixelSize = cExperiment.cTimelapse.pixelSize;
            cExperiment.image_rotation = cExperiment.cTimelapse.image_rotation;
            cExperiment.trapsPresent = cExperiment.cTimelapse.trapsPresent;
            cExperiment.timepointsToProcess = cExperiment.cTimelapse.timepointsToProcess;
            
            cExperiment.logger.add_arg('Default segmentation channel',cExperiment.searchString);
            cExperiment.logger.add_arg('Traps present',cExperiment.trapsPresent);
            cExperiment.logger.add_arg('Image rotation',cExperiment.image_rotation);
            cExperiment.logger.add_arg('Pixel size',cExperiment.pixelSize);
            cExperiment.logger.start_protocol('creating new experiment',length(positionsToLoad));
        else
            if ~isequal(cExperiment.pixelSize,cExperiment.cTimelapse.pixelSize)
                warning('pixel size differs between positions');
            end
            if ~isequal(cExperiment.image_rotation,cExperiment.cTimelapse.image_rotation)
                warning('image_rotation differs between positions');
            end
            if ~isequal(cExperiment.trapsPresent,cExperiment.cTimelapse.trapsPresent)
                warning('whether traps are present differs between positions');
            end
            if ~isequal(cExperiment.timepointsToProcess,cExperiment.cTimelapse.timepointsToProcess)
                warning('timepointsToProcess differs between positions');
            end
        end
        
        % The false input tells this function not to save the cExperiment 
        % each time. Will speed it up a bit.
        cExperiment.saveTimelapseExperiment(currentPos,false);
    end
    
    % set a housekeeping variables:
    
    %whether a position has had traps tracked
    cExperiment.posTrapsTracked = false(size(cExperiment.dirs));
    
    %whether a position has been segmented
    cExperiment.posSegmented = false(size(cExperiment.dirs));
    
    %whether a position has had cells tracked
    cExperiment.posTracked = false(size(cExperiment.dirs));
    
    % this is true when it is appropriate to reset old trapInfo
    cExperiment.clearOldTrapInfo = false(size(cExperiment.dirs));
    
    cExperiment.saveExperiment;
    
    % if experiment has no traps, the trap tracking must still be run to
    % initialise the trapsInfo. This causes no end of confusion, so I have
    % done it here automatically.
    if ~cExperiment.trapsPresent
        cExperiment.trackTrapsInTime(positionsToLoad);
    end
    
    % Finish logging protocol
    cExperiment.logger.complete_protocol;
catch err
    cExperiment.logger.protocol_error;
    rethrow(err);
end

end
