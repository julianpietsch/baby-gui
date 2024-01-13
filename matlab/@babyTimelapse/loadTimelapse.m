function loadTimelapse(cTimelapse,varargin)
% loadTimelapse Populate cTimepoint and initialise image properties
%
% Determines how many timepoints there are in the timelapse by identifying
% images who's name contains the searchString.
%
% INPUTS
% 
% cTimelapse            -  object of the babyTimelapse class.
% searchString          -  string.the string that appears in each image
%                          associated with a particular timepoint
% 
% Extra arguments are passed through to
% BABYTIMELAPSE.INITIALIZEIMAGEPROPERTIES.
%
% seaches through the timelapseDir for filenames containing the string
% searchString. Uses the ordered list of these to populate the cTimepoints
% - one cTimepoint for each matching file.
%
% See also BABYEXPERIMENT.CREATETIMELAPSEPOSITIONS

ip = inputParser;
ip.addParameter('SearchString','',@(x) ischar(x) && isrow(x));
ip.addParameter('TimepointsToProcess',[],@(x) isempty(x) || ...
    (isvector(x) && isnumeric(x)));
ip.KeepUnmatched = true;
ip.parse(varargin{:});

initargs = [fieldnames(ip.Unmatched),struct2cell(ip.Unmatched)]';
searchString = ip.Results.SearchString;

if isempty(searchString)
    searchString = inputdlg('Enter the string to search for the brightfield/DIC images','SearchString',1,{'Brightfield_002'});
    searchString = searchString{1};
end

%get names of all files in the timelapseDir folder
cTimelapse.channelNames={searchString};
folder=cTimelapse.timelapseDir;
tempdir=dir(folder);
names={tempdir(:).name};
files=sort(names);

% Read images into timelapse class
% Timelapse is a seletion of images from a file. These images must be
% loaded in the correct order from low to high numbers to ensure that the
% cell tracking performs correctly, and they must be rotated to ensure the
% trap correctly aligns with the images

cTimepointTemplate = cTimelapse.cTimepointTemplate;
cTimelapse.cTimepoint = cTimepointTemplate;

timepoint_index=0;
for n = 1:length(files)
    % check file name is an image and is not a hidden unix file
    if (~isempty(strfind(files{n},'tif'))|| ~isempty(strfind(files{n},'png')) || ~isempty(strfind(files{n},'TIF')))...
            && isempty(regexp(files{n},'^\.','once'))
        if ~isempty(strfind(files{n},searchString))
            
            cTimelapse.cTimepoint(timepoint_index+1) = cTimepointTemplate;
            cTimelapse.cTimepoint(timepoint_index+1).filename{end+1}=files{n};
            
            cTimelapse.cTimepoint(timepoint_index+1).trapLocations=[];
            timepoint_index=timepoint_index+1;
        end
    end
end

timepointsToProcess = ip.Results.TimepointsToProcess;
if isempty(timepointsToProcess)
    timepointsToProcess = 1:timepoint_index;
end
cTimelapse.timepointsToProcess = timepointsToProcess;

cTimelapse.timepointsProcessed = false(1,timepoint_index);

% load the first image
image = cTimelapse.returnSingleTimepointRaw(timepointsToProcess(1),1);

cTimelapse.initializeImageProperties(image,initargs{:});

end

