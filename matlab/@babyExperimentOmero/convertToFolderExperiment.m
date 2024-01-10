function cExperiment=convertToFolderExperiment(cExperimentOm,saveFolder)
%convertToFolderExperiment Convert Omero cExperiments to a folder version
%   cExperiment = cExperimentOm.convertSegmented(saveFolder)
%   converts an babyExperimentOmero object to an babyExperiment
%   object saved in a local folder.
%       - cExperiment: converted babyExperiment object
%       - saveFolder (optional): the folder in which to save the folder 
%           cExperiment; if omitted or empty, will prompt user for a save
%           location.

if nargin<2 || isempty(saveFolder)
    msg = 'Specify where you want to save the folder cExperiment...';
    fprintf('%s\n',msg);
    saveFolder=uigetdir('',msg);
    if isequal(saveFolder,0)
        error('A save folder must be specified for conversion...');
    end
end

if ~isdir(saveFolder)
    error('A valid save folder must be specified for conversion...');
end

%Create a new bare babyExperiment object
cExperiment = babyExperiment(true);
%Copy properties from the input cExperiment. Only copy ones that are
%not supposed to be different between the two versions:
ignoreFields = {'rootFolder','saveFolder'};
cExperiment.copyprops(cExperimentOm,ignoreFields);
cExperiment.saveFolder = saveFolder;
cExperiment.rootFolder = saveFolder;

% If the cExperiment has a log file, copy it to the save folder
old_log_file = fullfile(cExperimentOm.logger.file_dir,...
    cExperimentOm.logger.file_name);
if exist(old_log_file,'file')==2
    [success,msg] = copyfile(old_log_file,...
        fullfile(cExperiment.logger.file_dir,...
        cExperiment.logger.file_name));
    if ~success
        warning(msg);
    end
end

% In any case, log this conversion to the folder cExperiment's
% log file before attempting remaining conversion and upload:
cExperiment.logger.add_arg('saveFolder',cExperiment.saveFolder);
cExperiment.logger.start_protocol('converting to folder cExperiment',...
    length(cExperimentOm.dirs));

try
    % First save the cExperiment file
    cExperiment.saveExperiment;
    
    for d=1:length(cExperimentOm.dirs)
        % Load the template cTimelapseOm
        cTimelapseOm = cExperimentOm.returnTimelapse(d);
        
        % Create a new bare babyTimelapse object
        cTimelapse=babyTimelapse(true);
        
        % Trigger a PositionChanged event to notify cExperiment.logger
        babyLogging.changePos(cExperiment,d,cTimelapse);
        
        % Copy properties from cTimelapseOm to cTimelapse
        % Only copy ones that are not supposed to be different between the 
        % two versions
        ignoreFields = {'defaultExtractParameters'};
        cTimelapse.copyprops(cTimelapseOm,ignoreFields);
        cTimelapse.timelapseDir = fullfile(saveFolder,cExperimentOm.dirs{d});
        
        % Save the cTimelapse
        cExperiment.cTimelapse=cTimelapse;
        cExperiment.saveTimelapseExperiment(d,false);
    end
    
    % Finish logging protocol
    cExperiment.logger.complete_protocol;
catch err
    cExperiment.logger.protocol_error;
    rethrow(err);
end

end
