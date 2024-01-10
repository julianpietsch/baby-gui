function loadTimelapse(cTimelapse,searchString,image_rotation,trapsPresent,timepointsToLoad,pixel_size)
% loadTimelapse(cTimelapse,searchString,image_rotation,trapsPresent,timepointsToLoad,pixel_size)
%
% populates the cTimpoint field, determining how many timepoints there are
% in the timelapse by identifying images who's name contains the searchString.
%
% INPUTS
% 
% cTimelapse            -  object of the babyTimelapseOmero class
% searchString          -  string.the string that appears in each image
%                          associated with a particular timepoint
% image_rotation        -  counter clockwise rotation of images (in
%                          degrees) to perform when an image is requested.
%                          Generally rotated to align traps with the template.
% trapsPresent          -  a boolean that states whether there are traps
%                          present. Used at various stages of the
%                          processing.
% timepointsToLoad      -  unused in Omero case but kept for compatibility.
% pixel_size            -  width of a pixel in the image in micrometers.
%                          Default is 0.262 - the value for the swainlab
%                          miscroscopes at 60X.
% 
% 
% Populates the cTimepoints structure using information from Omero.
%
% other properties:
%   - rotation
%   - trapsPresent
%   - imSize
%   - rawImSize
% are also populated, by GUI if not provided.
%
% See also BABYEXPERIMENTOMERO.CREATETIMELAPSEPOSITIONS

cTimepointTemplate = cTimelapse.cTimepointTemplate;

cTimelapse.cTimepoint = cTimepointTemplate;

ntimepoints = cTimelapse.dataset.imageSize(5);
% Make sure cTimepoint and timepointsToProcess have the correct lengths
cTimelapse.cTimepoint(ntimepoints).filename=[];
cTimelapse.timepointsToProcess = 1:ntimepoints;

%Load first timepoint of this cTimelapse to fill out the remaining
%details
image = cTimelapse.returnSingleTimepointRaw(1,find(strcmp(cTimelapse.channelNames,searchString),1));

cTimelapse.initializeImageProperties(image,image_rotation,trapsPresent,pixel_size);


end

