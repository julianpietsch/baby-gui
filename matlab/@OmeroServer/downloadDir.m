function dsDir = downloadDir(this,dataset)
%DOWNLOADDIR Return the download directory for the specified dataset
%   If the directory for the dataset (specified either as an Omero Dataset 
%   object or its integer ID) does not already exist, it will be created. 
%   This function obeys the user Nursery, and will set the download directory
%   to the temporary directory if a local cache is not desired.

if ~isdir(this.cachePath)
    error('The "OmeroTemp" dir could not be found.');
end

if isnumeric(dataset) && floor(dataset)==dataset && dataset>0
    dsId = uint32(dataset); % force to be an unsigned integer
elseif isa(dataset,'omero.model.DatasetI')
    dsId = dataset.getId().getValue();
elseif isa(dataset,'OmeroDataset')
    dsId = dataset.dataset.getId().getValue();
else
    error('"dataset" must be an Omero Dataset object or its integer ID');
end

% Derive the name of the dataset cache directory:
dsDir = sprintf('dataset_%u',dsId);

% Check if the user wants files from this dataset cached locally:
sp = Nursery.active;
rowsel = strcmp(sp.localOmeroCache.dirname,dsDir);
if sum(rowsel)==0
    % Add a new row for this dataset since it's not in the table
    newrow = sp.defaultOmeroCacheRow;
    newrow.dirname = {dsDir};
    if ~isempty(sp.localOmeroCache)
        sp.localOmeroCache(end+1,:) = newrow;
    else
        sp.localOmeroCache = newrow;
    end
    rowsel = size(sp.localOmeroCache,1);
end
if ~sp.localOmeroCache.cExperiment(rowsel)
    % User does not want a local cache for this experiment
    dsDir = 'tmp';
end

if ~isdir(fullfile(this.cachePath,dsDir))
    % The folder does not yet exist, so create it
    mkdir(this.cachePath,dsDir);
end

% Return the full path to the directory
dsDir = fullfile(this.cachePath,dsDir);

end
