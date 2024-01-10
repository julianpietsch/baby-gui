function trapStack = getTimepointTraps(this,pos,timepoint,...
    channel,traps,stack_class,timepoint_is_unprocessed)
%ImageCache.getTimepointTraps Retrieve trap images at a single
%time point
%
%trapStack = imcache.getTimepointTraps(pos,timepoint,channel,
%         traps,stack_class,timepoint_is_unprocessed)
%
%Arguments:
%   - pos: single position index (leave empty for most recent)
%   - timepoint: single time point index
%   - channel: channel name as a char array
%   - traps: trap indices to return (leave empty for all)
%   - stack_class (optional): the image output type, one of
%       'raw' (original units), 'double' (normalised between
%       0 and 1; the default), 'uint8', or 'uint16'
%   - timepoint_is_unprocessed (optional): whether the time
%       point is an index into the filtered/processed
%       array or not (default is false).
%
%Returns:
%   - trapStack: array of images with time points along dim 3

if isempty(pos), pos = this.current_position; end
assert(isscalar(pos),...
    'this function only returns traps for a single position');
assert(isscalar(timepoint),...
    'this function only returns traps for a single timepoint');
if nargin<7 || isempty(timepoint_is_unprocessed)
    timepoint_is_unprocessed = false;
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

if timepoint_is_unprocessed
    timepoint = find(timepoint==this.timepointIndices,1);
    assert(~isempty(timepoint),...
        'the specified time point is missing in this "ImageCache"');
end

if nargin<5 || isempty(traps)
    traps = this.trapMaps{ipos}.keys;
    traps = [traps{:}];
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

% Load the timepoint into the cache if necessary
itraps = arrayfun(@(trap) this.trapMaps{ipos}(trap), traps);
areloaded = arrayfun(@(itrap) ...
    this.loadedtps{ipos,ichannel}{itrap}(timepoint), itraps);
if any(~areloaded)
    if ~isempty(this.savedir)
        for t=1:length(traps)
            % Attempt first to load this trap from the local cache
            try
                this.loadLocalTrapStack(pos,traps(t),channel);
                % refresh the array of loaded timepoints
                areloaded = this.loadedtps{ipos,ichannel}{itraps(t)}(timepoint);
            catch err
                if ~ismember(err.identifier,{this.Errors.NotSaved,...
                        this.Errors.BadImage})
                    rethrow(err);
                end
            end
        end
    end
    if any(~areloaded)
        this.loadTimepointWithoutChecks(ipos,ichannel,channelNum,timepoint);
    end
end

% Return the trapStack from the cache
output_class = stack_class;
if ~ismember(output_class,{'uint8','uint16'})
    output_class = 'double';
end
trapStack = zeros([this.stackSize(1:2),length(itraps)],output_class);
for t=1:length(itraps)
    trapStack(:,:,t) = this.imcache{ipos,ichannel}{itraps(t)}(:,:,timepoint);
end
if ~strcmp(stack_class,this.storage_class)
    switch stack_class
        case 'double'
            trapStack = double(trapStack)/double(intmax(this.storage_class));
        case 'raw'
            for t=1:length(itraps)
                traprange = this.imrange{ipos,ichannel}{itraps(t)};
                trapStack(:,:,t) = this.imraw(trapStack(:,:,t),traprange(1),traprange(2));
            end
        otherwise
            for t=1:length(itraps)
                traprange = this.imrange{ipos,ichannel}{itraps(t)};
                trapStack(:,:,t) = this.imraw(trapStack(:,:,t),traprange(1),traprange(2));
                trapStack(:,:,t) = this.imnorm(trapStack(:,:,t),traprange(1),traprange(2),stack_class);
            end
    end
end
end
