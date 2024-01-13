function loadTimelapse(cTimelapse,varargin)
% loadTimelapse Populate cTimepoint and initialise image properties
%
% Determines how many timepoints there are in the timelapse from the
% OMERODATASET object associated with this cTimelapse.
%
% INPUTS
% 
% - cTimelapse: BABYTIMELAPSEOMERO object
% 
% Extra arguments are passed through to
% BABYTIMELAPSE.INITIALIZEIMAGEPROPERTIES.
%
% See also BABYEXPERIMENT.CREATETIMELAPSEPOSITIONS

ip = inputParser;
ip.addParameter('TimepointsToProcess',[],@(x) isempty(x) || ...
    (isvector(x) && isnumeric(x)));
ip.KeepUnmatched = true;
ip.parse(varargin{:});

initargs = [fieldnames(ip.Unmatched),struct2cell(ip.Unmatched)]';

cTimepointTemplate = cTimelapse.cTimepointTemplate;
cTimelapse.cTimepoint = cTimepointTemplate;

% Make sure cTimepoint and timepointsToProcess have the correct lengths
ntimepoints = cTimelapse.dataset.imageSize(5);
cTimelapse.cTimepoint(ntimepoints).filename=[];

timepointsToProcess = ip.Results.TimepointsToProcess;
if isempty(timepointsToProcess)
    timepointsToProcess = 1:ntimepoints;
end
cTimelapse.timepointsToProcess = timepointsToProcess;

cTimelapse.timepointsProcessed = false(1,ntimepoints);

%Load first timepoint of this cTimelapse to fill out the remaining details
image = cTimelapse.returnSingleTimepointRaw(timepointsToProcess(1),1);

cTimelapse.initializeImageProperties(image,initargs{:});

end

