function [cExperiment,suffix] = loadcExperiment(this,dataset,suffix)
%LOADCEXPERIMENT Return a cExperiment loaded from Omero
%   cExperiment = obj.loadcExperiment(dataset,suffix) returns the
%   cExperiment with the specified suffix for the specified dataset (either
%   an Omero Dataset, its integer ID, or the unique ID string). If the 
%   fileAnnotation of the cExperiment is known, this can be passed instead 
%   of suffix; WARNING: in this case, no check will be performed to ensure
%   that the fileAnnotation belongs to the specified dataset.

if ~isa(dataset,'omero.model.DatasetI')
    if ischar(dataset)
        exptid = dataset;
        dataset = this.getDsListFromTag(exptid);
    elseif isnumeric(dataset)
        dataset = getDatasets(this.session,dataset);
    else
        error('dataset is of wrong type...');
    end
    if isempty(dataset)
        % experiment is not in the database
        error('The specified dataset "%s" does not exist.', exptid)
    else
        if length(dataset)>1
            warning('More than one dataset was found with that ID: choosing the first found.');
        end
        dataset = dataset(1);
    end
end

if isa(suffix,'omero.model.FileAnnotationI')
    % We do not need to retrieve fileAnnotations:
    fileAnnotation = suffix;
    fileName = char(fileAnnotation.getFile().getName().getValue());
    % Validate the file name and retrieve the suffix
    if isempty(regexp(fileName,'^cExperiment_.*\.mat$','once'))
        error('The specified fileAnnotation does not point to a valid cExperiment');
    else
        suffix = regexprep(fileName,'^cExperiment_(.*)\.mat$','$1');
    end
else
    % Retrieve fileAnnotations for this dataset
    [expNames,fileAnnotations] = this.listcExperiments(dataset);
    if sum(strcmp(expNames,suffix))~=1
        error('The specified experiment does not exist for this dataset or is not unique.');
    end
    fileAnnotation = fileAnnotations(strcmp(expNames,suffix));
end
    
PathName = this.downloadDir(dataset);

% Download the cExperiment:
FilePath = this.downloadFile(dataset,fileAnnotation,PathName);

l1 = load(FilePath);
cExperiment=l1.cExperiment;
% Note: Conversion for cases in which experiments were created from Omero but
% as babyExperiment objects, not babyExperimentOmero objects should be handled
% by the loadobj method of babyExperiment.

% Ensure that the cExperiment save folder is updated for the current
% computer:
cExperiment.saveFolder = PathName;

% Ensure the cExperiment log file is downloaded if it exists
if ~isempty(cExperiment.logFileAnnotation_id)
    lfa = getFileAnnotations(this.session,cExperiment.logFileAnnotation_id);
    this.downloadFile(dataset,lfa,PathName);
end

end
