function trapStack = getTrapTimepoints(this,pos,trap,channel,...
    timepoints,stack_class,timepoints_are_unprocessed)
%ImageCache.getTrapTimepoints Retrieve images for a single trap
%
%trapStack = imcache.getTrapTimepoints(pos,trap,channel,
%         timepoints,stack_class,timepoints_are_unprocessed)
%
%Arguments:
%   - pos: single position index (leave empty for most recent)
%   - trap: single trap index
%   - channel: channel name as a char array
%   - timepoints: array of time point indices
%   - stack_class (optional): the image output type, one of
%       'raw' (original units), 'double' (normalised between
%       0 and 1; the default), 'uint8', or 'uint16'
%   - timepoints_are_unprocessed (optional): whether the time
%       points are are indices into the filtered/processed
%       array or not (default is false).
%
%Returns:
%   - trapStack: array of images with time points along dim 3

if isempty(pos), pos = this.current_position; end
assert(isscalar(pos),...
    'this function only returns traps for a single position');
assert(isscalar(trap),'this function only returns timepoints for a single trap');
if nargin<7 || isempty(timepoints_are_unprocessed)
    timepoints_are_unprocessed = false;
end
if nargin<6 || isempty(stack_class)
    stack_class = 'double';
end
assert(ismember(stack_class,{'raw','double','uint8','uint16'}),...
    '"stack_class" must be one of "double", "uint8" or "uint16"');

if ~isKey(this.posMap,pos), this.addPos(pos); end
ipos = this.posMap(pos);
% Ensure the specified position is loaded:
this.loadTimelapse(pos);

if timepoints_are_unprocessed
    assert(all(ismember(timepoints,this.timepointIndices)),...
        'some specified time points are missing in this "ImageCache"');
    missingtps = setdiff(1:max(this.timepointIndices),this.timepointIndices);
    [~,invtps] = sort([this.timepointIndices,missingtps]);
    timepoints = invtps(timepoints);
end

if ~isKey(this.channelMap,channel)
    this.addPosChannel(pos,channel);
end
ichannel = this.channelMap(channel);

% Ensure that the appropriate caches are initialised
% Determine the correct channel number from the timelapse
channelNum = find(strcmp(this.cTimelapse.channelNames,channel),1);
assert(~isempty(channelNum),this.Errors.NotAvailable,...
    'The channel %s is not available at position %u',channel,pos);
this.init_cache(ipos,ichannel);

% Load the timepoints into the cache if necessary
itrap = this.trapMaps{ipos}(trap);
areloaded = this.loadedtps{ipos,ichannel}{itrap}(timepoints);
if any(~areloaded)
    if ~isempty(this.savedir)
        % Attempt first to load this trap from the local cache
        try
            this.loadLocalTrapStack(pos,trap,channel);
            % refresh the array of loaded timepoints
            areloaded = this.loadedtps{ipos,ichannel}{itrap}(timepoints);
        catch err
            if ~ismember(err.identifier,{this.Errors.NotSaved,...
                    this.Errors.BadImage})
                rethrow(err);
            end
        end
    end
    for tp=timepoints(~areloaded)
        this.loadTimepointWithoutChecks(ipos,ichannel,channelNum,tp);
    end
end

% Return the trapStack from the cache
trapStack = this.imcache{ipos,ichannel}{itrap}(:,:,timepoints);
if ~strcmp(stack_class,this.storage_class)
    switch stack_class
        case 'double'
            trapStack = double(trapStack)/double(intmax(this.storage_class));
        case 'raw'
            traprange = this.imrange{ipos,ichannel}{itrap};
            trapStack = this.imraw(trapStack,traprange(1),traprange(2));
        otherwise
            traprange = this.imrange{ipos,ichannel}{itrap};
            trapStack = this.imraw(trapStack,traprange(1),traprange(2));
            trapStack = this.imnorm(trapStack,traprange(1),traprange(2),stack_class);
    end
end
end
