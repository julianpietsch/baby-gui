function trackTrapsInTime(cExperiment,positionToTrack,first_timepoint,last_timepoint,varargin)
% function trackTrapsInTime(cExperiment,positionToTrack,first_timepoint,last_timepoint)
%
% tracks the traps in time from first_timepoint to last_timepoint, or if
% first_timepoint is 'Timepoints', then last_timepoint should be either
% empty, or a vector specifying the time points to process.
% If traps are present, requires that they have been selected.
% 
% If there are no traps it populates the trapInfo and trapLocations fields
% of cTimelapse.cTimepoint(:) as though there were one large trap of size
% cTimelapse.imSize.
%
% See also, CTRAPSELECTDISPLAY, BABYEXPERIMENT.IDENTIFYTRAPSTIMELAPSES,
% BABYTIMELAPSE.IDENTIFYTRAPLOCATIONSSSINGLETP

if nargin<2 || isempty(positionToTrack)
    positionToTrack = 1:numel(cExperiment.dirs);
end

if nargin<3 || isempty(first_timepoint)
    first_timepoint = min(cExperiment.timepointsToProcess);
end

if nargin<4, last_timepoint = []; end

if isequal(first_timepoint,'Timepoints')
    timepoints = last_timepoint;
else
    if isempty(last_timepoint)
        last_timepoint = max(cExperiment.timepointsToProcess);
    end
    timepoints = first_timepoint:last_timepoint;
end

% Start logging protocol
cExperiment.logger.start_protocol('tracking traps in time',length(positionToTrack));
try

for posi = 1:numel(positionToTrack)
    pos = positionToTrack(posi);
    cExperiment.loadCurrentTimelapse(pos);
    cExperiment.cTimelapse.trackTrapsThroughTime(timepoints,varargin{:});
    cExperiment.saveTimelapseExperiment([],false);
    cExperiment.posTrapsTracked(pos) = true;
end

cExperiment.saveExperiment;

% Finish logging protocol
cExperiment.logger.complete_protocol;

catch err
    cExperiment.logger.protocol_error;
    rethrow(err);
end

end
