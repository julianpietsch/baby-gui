function loadTimelapse(cTimelapse,searchString,image_rotation,trapsPresent,timepointsToLoad,pixel_size,image_flipud)
% loadTimelapse(cTimelapse,searchString,image_rotation,trapsPresent,timepointsToLoad,pixel_size)
%
% populates the cTimpoint field, determining how many timepoints there are
% in the timelapse by identifying images who's name contains the searchString.
%
% INPUTS
% 
% cTimelapse            -  object of the babyTimelapse class.
% searchString          -  string.the string that appears in each image
%                          associated with a particular timepoint
% image_rotation        -  counter clockwise rotation of images (in
%                          degrees) to perform when an image is requested.
%                          Generally rotated to align traps with template.
% trapsPresent          -  a boolean that states whether there are traps
%                          present. Used at various stages of the
%                          processing.
% timepointsToLoad      -  scalar. Declares the maximum number of
%                          timepoints to load. Largely superceded by
%                          timepointsToProcess - set later and just limits
%                          the timepoints to segment and extract.
% pixel_size            -  width of a pixel in the image in micrometers.
%                          Default is 0.262 - the value for the swainlab
%                          miscroscopes at 60X.
%
%
% seaches through the timelapseDir for filenames containing the string
% searchString. Uses the ordered list of these to populate the cTimepoints
% - one cTimepoint for each matching file.
%
% expects images to be png,tif or TIF format.
%
% other properties:
%   - rotation
%   - trapsPresent
%   - imSize
%   - rawImSize
% are also populated, by GUI if not provided.
%
% See also BABYEXPERIMENT.CREATETIMELAPSEPOSITIONS


%get names of all files in the timelapseDir folder
cTimelapse.channelNames={searchString};
folder=cTimelapse.timelapseDir;
tempdir=dir(folder);
names={tempdir(:).name};
files=sort(names);
folder=[folder filesep];

% Read images into timelapse class
% Timelapse is a seletion of images from a file. These images must be
% loaded in the correct order from low to high numbers to ensure that the
% cell tracking performs correctly, and they must be rotated to ensure the
% trap correctly aligns with the images

cTimepointTemplate = cTimelapse.cTimepointTemplate;

cTimelapse.cTimepoint = cTimepointTemplate;


timepoint_index=0;
for n = 1:length(files);
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

cTimelapse.timepointsToProcess = 1:timepoint_index;
cTimelapse.timepointsProcessed = false(1,timepoint_index);

if nargin>=6 && ~isempty(timepointsToLoad)
    if max(timepointsToLoad)>length(cTimelapse.cTimepoint)
        timepointsToLoad=timepointsToLoad(timepointsToLoad<=length(cTimelapse.cTimepoint));
    end
    cTimelapse.cTimepoint=cTimelapse.cTimepoint(timepointsToLoad);
end

% load this first image via imread since the following property 
image = cTimelapse.returnSingleTimepointRaw(1,1);

cTimelapse.initializeImageProperties(image,image_rotation,trapsPresent,pixel_size,image_flipud);


end

