function [suffixes,fileAnnotations] = listcExperiments(this,dataset)
%LISTCEXPERIMENTS List the cExperiments for the given dataset
%   suffixes = OmeroDatabase.listcExperiments(dataset) returns a cellstr of
%   the cExperiments found for the specified dataset (either an Omero 
%   Dataset object, its integer ID, or the unique ID string).
%
%   [~,fileAnnotations] = OmeroDatabase.listcExperiments(dataset) returns
%   instead the fileAnnotations of the available cExperiments.

if ischar(dataset)
    exptid = dataset;
    dataset = this.getDsListFromTag(exptid);
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

% Get the fileAnnotations for this dataset. Note that the following Omero
% function accepts either numeric IDs or Dataset objects:
fileAnnotations = getDatasetFileAnnotations(this.session,dataset);

% Get the names of the files
faNames = arrayfun(@(x) char(x.getFile().getName().getValue()),...
    fileAnnotations,'Uni',0);

% Filter to cExperiments:
cExps = ~cellfun(@isempty, regexp(faNames,'^cExperiment_.*\.mat$'));
fileAnnotations = fileAnnotations(cExps);
faNames = faNames(cExps);

% Determine the suffixes:
suffixes = regexprep(faNames,'^cExperiment_(.*)\.mat$','$1');

end

