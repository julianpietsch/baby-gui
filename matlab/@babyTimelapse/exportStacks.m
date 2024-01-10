function exportStacks(cTimelapse,exportdir,channel,timepoints,traps)

assert(isdir(exportdir),'Export directory is invalid');

assert(isstruct(cTimelapse.metadata) && isfield(cTimelapse.metadata,'exptid'),...
    'The cTimelapse metadata must include the experiment ID "exptid"');
exptid = cTimelapse.metadata.exptid;
assert(ischar(exptid) && isvector(exptid),'The cTimelapse "exptid" is invalid');

assert(isfield(cTimelapse.metadata,'posname'),...
    'The cTimelapse metadata must include the position name "posname"');
posname = cTimelapse.metadata.posname;
assert(ischar(posname) && isvector(posname),'The cTimelapse "posname" is invalid');

if nargin<3 || isempty(channel)
    channel = 1;
end

if nargin<4 || isempty(timepoints)
    timepoints = 1:length(cTimelapse.cTimepoint);
end

if nargin<5 || isempty(traps)
    traps = 1:length(cTimelapse.cTimepoint(1).trapInfo);
end

chNames = cTimelapse.channelNames;
nch = length(chNames);

if isnumeric(channel) && channel>=1 && channel<=nch && round(channel)==channel
    channelNum = channel;
    channel = chNames{channelNum};
else
    assert(ischar(channel) && isvector(channel) && ismember(channel,chNames),...
        '"channel" argument must be a channel index or name');
    channelNum = find(strcmp(cTimelapse.channelNames,channel),1);
end

for tp=timepoints
    % Assume 16 bit raw image, but process to get correct registration
    % (i.e., rotation/offset) with the trap definitions
    tpIm = cTimelapse.returnSingleTimepoint(tp,channelNum,'stack');
    trapIms = zeros([cTimelapse.trapImSize,length(traps),size(tpIm,3)],'uint16');
    for z=1:size(tpIm,3)
        trapIms(:,:,:,z) = round(cTimelapse.returnTrapsFromImage(tpIm(:,:,z),tp,traps));
    end
    trapIms = permute(trapIms,[1,2,4,3]);
    for ti=1:length(traps)
        fname = fullfile(exportdir,sprintf('%s_trap%03u_tp%04u_%s.png',...
            posname,traps(ti),tp,channel));
        
        % Fill description structure
        desc = struct();
        desc.experimentID = exptid;
        desc.position = posname;
        desc.trap = traps(ti);
        desc.tp = tp;
        desc.channel = channel;
        desc.pixel_size = cTimelapse.pixelSize;
        
        save_tiled(trapIms(:,:,:,ti),fname,desc);
    end
end
    
end