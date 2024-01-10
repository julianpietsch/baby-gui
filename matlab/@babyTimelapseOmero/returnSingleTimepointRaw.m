function timepointIm=returnSingleTimepointRaw(cTimelapse,timepoint,channel)
% timepointIm = returnSingleTimepointRaw(cTimelapse,timepoint,channel)
%
% returns the raw image stack before corrections (rescaling, shifting,
% background correction,projection) are made. Intended as a background
% function, generally better to use TIMELAPSETRAPS.RETURNSINGLETIMEPOINT
%
% Important to note that in the case of numerous 'stack' channels this will
% return a stack, as oppose to TIMELAPSETRAPS.RETURNSINGLETIMEPOINT which
% generally returns a projected image.
%
% This version of the function is for experiments stored on an OMERO server,
% so images are downloaded from a database. 
%
% See also TIMELAPSETRAPS.RETURNSINGLETIMEPOINT, OMERODATASET


if iscell(cTimelapse.channelNames)
    channelName=cTimelapse.channelNames{channel};
else
    channelName=cTimelapse.channelNames;
    channel = 1;
end

ds = cTimelapse.dataset;
assert(ds.posNum==cTimelapse.posNum,...
    'cTimelapse OmeroDataset has the wrong position set');

zrange = []; % select all Z sections by default
% Always prefer to pick a raw channel name
chInd = find(strcmp(ds.channelNames,channelName));
if isempty(chInd)
    % Otherwise attempt to parse the Z section
    match = regexp(channelName,'^(.*)_(\d{3})$','tokens','once');
    if ~isempty(match)
        chInd = find(strcmp(ds.channelNames,match{1}));
        zrange = str2double(match{2});
        assert(~isnan(zrange) && zrange>0 && zrange<=ds.imageSize(3),...
            'cTimelapse specifies a Z section that does not exist on OMERO');
    end
end
if numel(chInd)>1
    error('Cannot unambiguously identify the cTimelapse channel on OMERO');
end

if ~isempty(chInd) && ~ds.hasZstacks(chInd), zrange = 1; end

if isempty(chInd)
    if isempty(zrange)
        timepointIm = zeros(ds.imageSize(1:3,'uint16'));
    else
        timepointIm = zeros(ds.imageSize(1:2,'uint16'));
    end
else
    timepointIm = ds.getHypercube('Z',zrange,'C',chInd,'T',timepoint);
end

%Images are flipped in both directions on Omero upload - so if the
%data was segmented from a folder before being uploaded/converted
%then it should be flipped to ensure data is consistent with any
%segmentation results
if strcmp(cTimelapse.segmentationSource,'Folder')
    timepointIm=rot90(timepointIm);
    timepointIm=flipud(timepointIm);
end

% A bug in the microscope code meant that some channels were saved
% flipped relative to brightfield. The flipchannels property can be
% used to specify which channels should be flipped to match the
% orientation of the brightfield image.
if ~isempty(cTimelapse.flipchannels) && ...
        length(cTimelapse.flipchannels)>channel && ...
        cTimelapse.flipchannels(channel)
    timepointIm=flipud(timepointIm);
end

if isempty(cTimelapse.rawImSize)
    cTimelapse.rawImSize = [size(timepointIm,1),size(timepointIm,2)];
end

end
