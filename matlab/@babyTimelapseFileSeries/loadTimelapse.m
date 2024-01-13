function loadTimelapse(cTimelapse,varargin)
% loadTimelapse Populate cTimepoint and initialise image properties
%
% Determines how many timepoints there are in the timelapse from the
% IMAGEREADERFILESERIES object associated with this cTimelapse.
%
% INPUTS
% 
% - cTimelapse: BABYTIMELAPSEFILESERIES object
% - PixelSize: Size of a pixel in microns. If left empty or unspecified,
% attempts to fill this property from the pixelSize property of the
% IMAGEREADERFILESERIES object.
% 
% Extra arguments are passed through to
% BABYTIMELAPSE.INITIALIZEIMAGEPROPERTIES.
%
% See also BABYEXPERIMENT.CREATETIMELAPSEPOSITIONS

ip = inputParser;
ip.addParameter('TimepointsToProcess',[],@(x) isempty(x) || ...
    (isvector(x) && isnumeric(x)));
ip.addParameter('PixelSize',[],@(x) isempty(x) || ...
    (isscalar(x) && isnumeric(x)));
ip.KeepUnmatched = true;
ip.parse(varargin{:});

initargs = [fieldnames(ip.Unmatched),struct2cell(ip.Unmatched)]';
pixel_size = ip.Results.PixelSize;
if isempty(pixel_size)
    pixel_size = cTimelapse.reader.pixelSize;
end

cTimepointTemplate = cTimelapse.cTimepointTemplate;
cTimelapse.cTimepoint = cTimepointTemplate;

ntps = cTimelapse.reader.imageSize(5);
for tp=1:ntps
    cTimelapse.cTimepoint(tp) = cTimepointTemplate;
    cTimelapse.cTimepoint(tp).trapLocations=[];
end

timepointsToProcess = ip.Results.TimepointsToProcess;
if isempty(timepointsToProcess)
    timepointsToProcess = 1:ntps;
end
cTimelapse.timepointsToProcess = timepointsToProcess;

cTimelapse.timepointsProcessed = false(1,ntps);

image = cTimelapse.returnSingleTimepointRaw(timepointsToProcess(1),1);

cTimelapse.initializeImageProperties(image,'PixelSize',pixel_size,initargs{:});

end
