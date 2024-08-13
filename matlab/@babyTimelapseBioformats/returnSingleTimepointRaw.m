function timepointIm = returnSingleTimepointRaw(cTimelapse,timepoint,channel)
% timepointIm = returnSingleTimepointRaw(cTimelapse,timepoint,channel)
%
% returns the raw image stack before corrections (rescaling, shifting,
% background correction,projection) are made. Intended as a background
% function, generally better to use BABYTIMELAPSE.RETURNSINGLETIMEPOINT
%
% Important to note that in the case of numerous 'stack' channels this will
% return a stack, as oppose to BABYTIMELAPSE.RETURNSINGLETIMEPOINT which
% generally returns a projected image (unless the 'stack' 'type' is
% chosen).
%
% INPUTS
% timepoint   :   number indicating the desired timepoint (will access this
%                 element of the cTimepoint array).
% channel     :   number (default 1) indicating which of the channels in
%                 cTimelapse.channelNames to use to identify appropriate
%                 files to load.
% 
% A timepoint and channel are specified, and the channel index is used to
% get a channel string from the cTimelapse.channelNames field. all
% filenames associated with that timepoint (i.e. in cTimepoint.filename
% cell array) that contain the channel string are identified. If there is
% only one this is loaded, whereas if there are more than one they are all
% loaded and put into a stack.
%
% The function also takes a number of other liberties. If the timelapseDir
% is set to 'ignore' it takes the filenames to be absolute (this can be
% done using the method makeFileNamesAbsolute). If this is not the case, it
% constructs the file name from timelapseDir and filename{i}. In doing so
% it takes the liberty of resaving the file name with only the relative
% path (i.e. it throws away anything behind the last / or \).
%
% before loading the existence of the putative file is checked, and if it
% doesn't exist the user is asked to specify a new file location. This is
% done using a special GUI that rechecks whether the file has been found
% every 5 seconds, this is done in case the disk was just temporarily
% disconnected. For cExperiments this is best done using the
% changeRootDirsAll method. This is not done if the timelapseDir is 'ignore'
% and filenames are absolute, since in this case it is very hard to work
% out where the files are (since they could be in different folders).
%
% See also BABYTIMELAPSE.RETURNSINGLETIMEPOINT,
% BABYTIMELAPSE.MAKEFILENAMESABSOLUTE,
% BABYEXPERIMENT.CHANGEROOTDIRALLTIMELAPSES

if nargin<3 || isempty(channel), channel = 1; end

channelName = cTimelapse.channelNames{channel};
reader = cTimelapse.reader;
zrange = []; % select all Z sections by default
% Always prefer to pick a raw channel name
chInd = find(strcmp(reader.channels,channelName));
if isempty(chInd)
    % Otherwise attempt to parse the Z section
    match = regexp(channelName,'^(.*)_(\d{3})$','tokens','once');
    if ~isempty(match)
        chInd = find(strcmp(reader.channels,match{1}));
        zrange = str2double(match{2});
        assert(~isnan(zrange) && zrange>0 && zrange<=reader.imageSize(3),...
            'cTimelapse specifies a Z section that does not exist');
    end
end
if numel(chInd)>1
    error('Cannot unambiguously identify the cTimelapse channel from file');
end

if isempty(zrange)
    timepointIm = reader.getTimepoint(timepoint,'C',chInd);
else
    timepointIm = reader.getTimepoint(timepoint,'C',chInd,'Z',zrange);
end
cTimelapse.rawImSize = [size(timepointIm,1),size(timepointIm,2)];

end
