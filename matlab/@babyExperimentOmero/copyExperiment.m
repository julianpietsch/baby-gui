function copyExperiment(cExperiment,new_suffix,delete_extracted)
% BABYEXPERIMENTOMERO.COPYEXPERIMENT
%
%   cExperiment.copyExperiment(cExperiment,new_suffix,delete_extracted)
%   Updates the suffix in the cExperiment and cTimelapse files and re-saves
%   as new files on Omero.

%First check what suffix names exist for this dataset:
expNames = cExperiment.dataset.cExperiments;

if nargin<2 || isempty(new_suffix)
    %Need to create a new cExperiment with a name different from
    %any of the existing ones, so bring up a GUI:
    new_suffix = cExperiment.dataset.getValidExpNameGUI;
    if ~new_suffix, return; end % user cancelled
end

if nargin<3 || isempty(delete_extracted)
    delete_extracted = false;
end

assert(ischar(new_suffix) & isvector(new_suffix),...
    'The new experiment suffix must be specified as a char vector');
assert(~ismember({new_suffix},expNames),...
    'That experiment suffix is already in use');

old_suffix = cExperiment.rootFolder;
    
% Log the copy command to the old log:
logmsg(cExperiment,'=====================');
logmsg(cExperiment,'Copying experiment from suffix "%s" to "%s"',old_suffix,new_suffix);
logmsg(cExperiment,'---------------------');

% The save methods for both cExperiment and cTimelapse create new files on
% Omero if they do not already exist, so just need to clear references to
% annotation ID and resave all cTimelapses:
for diri=1:length(cExperiment.dirs)
    cExperiment.rootFolder = old_suffix;
    cTimelapse = cExperiment.loadCurrentTimelapse(diri);
    cTimelapse.fileAnnotation_id = [];
    cExperiment.rootFolder = new_suffix;
    if delete_extracted
        cTimelapse.extractedData = [];
    end
    cExperiment.saveTimelapseExperiment([],false);
end
% Loop is finished with rootFolder as new_suffix

% Close the logger file to ensure that it gets updated to the new file name
cExperiment.logger.close_logfile;

% Log the copy command to a new log file in the new location:
logmsg(cExperiment,'=====================');
logmsg(cExperiment,'Experiment copied from suffix "%s" to "%s"',old_suffix,new_suffix);
logmsg(cExperiment,'---------------------');

cExperiment.fileAnnotation_id = [];
cExperiment.logFileAnnotation_id = [];
if delete_extracted
    cExperiment.cellInf = [];
end

% Save the experiment (which also uploads log file):
cExperiment.saveExperiment;

end

