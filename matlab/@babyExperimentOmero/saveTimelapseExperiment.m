function saveTimelapseExperiment(cExperiment,currentPos,saveCE)
% saveTimelapseExperiment(cExperiment,currentPos,saveCE)
%
%   Saves cExperiment.cTimelapse as .mat file and uploads to OMERO server.
%   
%   File is first saved locally to the cExperiment.saveFolder and the
%   cTimelapse is given a suffix specified by cExperiment.rootFolder.
%
%   If currentPos is not provided, cExperiment.currentPos (populated when
%   BABYEXPERIMENT.LOADCURRENTTIMELAPSE is called) is used. It will be
%   empty if cExperiment.cTimelape has been replaced by a non-identical 
%   object (see BABYEXPERIMENT.SET.CTIMELAPSE)
% 
%   Third input is boolean - saveCE: logical - if true,
%   save the cExperiment file as well as the timelapse. Defaults to false.
%
%   See also, BABYEXPERIMENT.LOADCURRENTTIMELAPSE

if nargin<2 || isempty(currentPos)
    currentPos = cExperiment.currentPos;
end

if nargin<3
    saveCE=true;
end

%Save code for Omero loaded cTimelapses - upload cExperiment file to
%Omero database. Use the alternative method saveExperiment if you want
%to save only the cExperiment file.

cTimelapse = cExperiment.cTimelapse;

% Name of cTimelapse position should be equivalent to the indexed position
% name from cExperiment.dirs
posName = cExperiment.dirs{currentPos};
if ~strcmp(cTimelapse.dataset.pos,posName)
    error('"currentPos" does not match the current cTimelapse');
end

fileName=[cExperiment.saveFolder filesep posName 'cTimelapse_' cExperiment.rootFolder '.mat'];

cT_description = 'cTimelapse file uploaded by @babyExperiment.saveTimelapseExperiment';

srv = cExperiment.server;
ds = cExperiment.dataset.dataset;

% If fileAnnotation IDs have not already been set for this cTimelapse,
% determine them before saving the object.
if isempty(cTimelapse.fileAnnotation_id)
    fileAnnotations = getDatasetFileAnnotations(srv.session,ds);
    % Make a 'dummy' update call to ensure the files have an associated
    % fileAnnotation:
    cT_fA = srv.updateFile(ds,fileName,'dummy',true,...
        'dsFiles',fileAnnotations,'description',cT_description);
    cTimelapse.fileAnnotation_id = cT_fA.getId().getValue();
else
    % The ID is already known, so retrieve fileAnnotations directly:
    cT_fA = getFileAnnotations(srv.session,cTimelapse.fileAnnotation_id);
end

% Save cTimelapse object, then restore image and OmeroDatabase objects
save(fileName,'cTimelapse');

% Update the file on Omero:
srv.updateFile(ds,fileName,'dsFiles',cT_fA,'description',cT_description);

if saveCE, cExperiment.saveExperiment; end

end
