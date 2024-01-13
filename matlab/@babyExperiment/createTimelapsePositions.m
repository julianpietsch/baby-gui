function createTimelapsePositions(cExperiment,varargin)
%createTimelapsePositions Create babyTimelapses for each position/point
%
% Creates a babyTimelapse object (or suitable subclass) for each identified
% position/point. Various configuration options are set in the process,
% e.g., pixel size, image rotation/flip and trap templates. Some inputs
% that are not provided may result in a GUI prompt to provide them.
%
% Additional arguments to this function are passed through
% BABYTIMELAPSE.LOADTIMELAPSE and then further on to
% BABYTIMELAPSE.INITIALIZEIMAGEPROPERTIES. See those functions for details.

ip = inputParser;
ip.addOptional('poses',[],@(x) isempty(x) || (isnumeric(x) && isvector(x)));
ip.KeepUnmatched = true;
ip.parse(varargin{:});

loadargs = [fieldnames(ip.Unmatched),struct2cell(ip.Unmatched)]';
poses = ip.Results.poses;

if isempty(poses)
    poses=1:numel(cExperiment.dirs);
end

% Start adding arguments to experiment creation protocol log:
cExperiment.logger.add_arg('Root folder',cExperiment.rootFolder);
cExperiment.logger.add_arg('Save folder',cExperiment.saveFolder);

% The other arguments are added and the protocol started after the first
% call to loadTimelapse below...

try
    % Load timelapses
    for i=1:length(poses)
        currentPos=poses(i);
        cExperiment.cTimelapse = cExperiment.newTimelapse(currentPos);
        cExperiment.cTimelapse.metadata = cExperiment.metadata;
        cExperiment.cTimelapse.metadata.posname = cExperiment.dirs{currentPos};

        % Trigger a PositionChanged event to notify babyLogging
        babyLogging.changePos(cExperiment,currentPos,cExperiment.cTimelapse);

        cExperiment.cTimelapse.loadTimelapse(loadargs{:})

        if i==1
            cExperiment.pixelSize=cExperiment.cTimelapse.pixelSize;
            cExperiment.image_flipud=cExperiment.cTimelapse.image_flipud;
            cExperiment.image_rotation=cExperiment.cTimelapse.image_rotation;
            cExperiment.trapTemplates = cExperiment.cTimelapse.trapTemplates;
            cExperiment.trapTemplateChannel = cExperiment.cTimelapse.trapTemplateChannel;
            cExperiment.trapsPresent = cExperiment.cTimelapse.trapsPresent;
            cExperiment.timepointsToProcess = cExperiment.cTimelapse.timepointsToProcess;

            % After the first call to loadTimelapse, arguments should now
            % be set and copied to all other timelapses
            loadargs = {'PixelSize',cExperiment.pixelSize,...
                'FlipImage',cExperiment.image_flipud,...
                'ImageRotation',cExperiment.image_rotation,...
                'TrapTemplate',cExperiment.trapTemplates,...
                'TrapTemplateChannel',cExperiment.trapTemplateChannel,...
                'TrapsPresent',cExperiment.trapsPresent,...
                'TimepointsToProcess',cExperiment.timepointsToProcess};
            
            extraChannels = ~ismember(cExperiment.cTimelapse.channelNames,cExperiment.channelNames);
            if any(extraChannels)
                cExperiment.channelNames(end+1:end+sum(extraChannels)) = ...
                    cExperiment.cTimelapse.channelNames(extraChannels);
                if numel(cExperiment.cTimelapse.channelNames) == 1
                    % Assume this was the SearchString specified in
                    % babyTimelapse.loadTimelapse and add to args
                    loadargs(end+1:end+2) = {'SearchString',cExperiment.cTimelapse.channelNames{1}};
                end
            end

            % Also start logging the creation protocol with set args:
            cExperiment.logger.add_arg('Pixel size',cExperiment.pixelSize);
            cExperiment.logger.add_arg('Image flip',cExperiment.image_flipud);
            cExperiment.logger.add_arg('Image rotation',cExperiment.image_rotation);
            cExperiment.logger.add_arg('Traps present',cExperiment.trapsPresent);
            if cExperiment.trapsPresent
                cExperiment.logger.add_arg('Trap template channel',...
                    cExperiment.channelNames{cExperiment.trapTemplateChannel});
            end
            cExperiment.logger.start_protocol('creating new experiment',length(poses));
        else
            if ~isequal(cExperiment.pixelSize,cExperiment.cTimelapse.pixelSize)
                warning('pixel size differs between positions');
            end
            if ~isequal(cExperiment.image_rotation,cExperiment.cTimelapse.image_rotation)
                warning('image_rotation differs between positions');
            end
            if ~isequal(cExperiment.image_flipud,cExperiment.cTimelapse.image_flipud)
                warning('image_flipud differs between positions');
            end
            if ~isequal(cExperiment.trapsPresent,cExperiment.cTimelapse.trapsPresent)
                warning('whether traps are present differs between positions');
            end
            if ~isequal(cExperiment.trapTemplates,cExperiment.cTimelapse.trapTemplates)
                warning('trap templates differ between positions');
            end
            if ~isequal(cExperiment.trapTemplateChannel,cExperiment.cTimelapse.trapTemplateChannel)
                warning('trap template channel differs between positions');
            end
            if ~isequal(cExperiment.timepointsToProcess,cExperiment.cTimelapse.timepointsToProcess)
                warning('timepointsToProcess differs between positions');
            end
        end

        cExperiment.saveTimelapseExperiment(currentPos);
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

    % If experiment has no traps, the trap tracking must still be run to
    % initialise the trapsInfo. This causes no end of confusion, so it is
    % done here automatically.
    if ~cExperiment.trapsPresent
        cExperiment.trackTrapsInTime(poses);
    end

    % Finish logging protocol
    cExperiment.logger.complete_protocol;
catch err
    cExperiment.logger.protocol_error;
    rethrow(err);
end

end
