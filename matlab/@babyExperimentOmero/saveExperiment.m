function saveExperiment(cExperiment)
% SAVEEXPERIMENT(cExperiment)
% Uploads changes to the experiment to the Omero database

%cExperiment file

fileName = [cExperiment.saveFolder filesep 'cExperiment_' cExperiment.rootFolder '.mat'];
if iscell(fileName)
    fileName=fileName{:}; % equivalent to fileName{1};
end

server = cExperiment.server;
ds = cExperiment.dataset.dataset;

cE_description = 'cExperiment file uploaded by @babyExperiment.saveTimelapseExperiment';
lF_description = 'cExperiment log file uploaded by @babyExperiment.saveExperiment';

% If fileAnnotation IDs have not already been set for this experiment,
% determine them before saving the object. Start by getting the
% fileAnnotations for this dataset:
if isempty(cExperiment.fileAnnotation_id) || isempty(cExperiment.logFileAnnotation_id)
    fileAnnotations = getDatasetFileAnnotations(server.session,ds);
else
    % Both IDs are already known, so retrieve fileAnnotations directly:
    fileAnnotations = getFileAnnotations(server.session,...
        [cExperiment.fileAnnotation_id,cExperiment.logFileAnnotation_id]);
end
fA_Ids = arrayfun(@(x) x.getId().getValue(),fileAnnotations);

% Now check if IDs are specified, and if not, make a 'dummy' update call to
% ensure the files have an associated fileAnnotation:
if isempty(cExperiment.fileAnnotation_id)
    cE_fA = server.updateFile(ds,fileName,'dummy',true,...
        'dsFiles',fileAnnotations,'description',cE_description);
    cExperiment.fileAnnotation_id = cE_fA.getId().getValue();
else
    cE_fA = fileAnnotations(fA_Ids==cExperiment.fileAnnotation_id);
end
logFileName = fullfile(cExperiment.logger.file_dir,cExperiment.logger.file_name);
if isempty(cExperiment.logFileAnnotation_id)
    if ~isempty(cExperiment.logger), cExperiment.logger.close_logfile; end
    lF_fA = server.updateFile(ds,logFileName,'dummy',true,...
        'dsFiles',fileAnnotations,'description',lF_description);
    cExperiment.logFileAnnotation_id = lF_fA.getId().getValue();
else
    lF_fA = fileAnnotations(fA_Ids==cExperiment.logFileAnnotation_id);
end

% cTimelapse is saved in another file, but keep handle to restore later
cTimelapse = cExperiment.cTimelapse;
cExperiment.cTimelapse = [];

save(fileName,'cExperiment');

%Restore the cExperiment object
cExperiment.cTimelapse = cTimelapse;

% Update the files on Omero:
server.updateFile(ds,fileName,'dsFiles',cE_fA,...
    'description',cE_description);
% Only update the log file if it exists (some people may have
% shouldLog=false):
if exist(logFileName,'file')==2
    server.updateFile(ds,logFileName,'dsFiles',lF_fA,...
        'description',lF_description);
end

end
