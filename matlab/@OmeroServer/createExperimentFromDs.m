function cExperiment = createExperimentFromDs(this,dataset,varargin)
%createExperimentFromDs Create cExperiment for an Omero data set

ip = inputParser;
ip.addOptional('expName',[],@(x) isempty(x) || (isrow(x) && ischar(x)));
ip.addParameter('pixelSize',[],@(x) isempty(x) || (isscalar(x) && isnumeric(x)));
ip.addParameter('trapChannel','',@(x) isempty(x) || (isrow(x) && ischar(x)));
ip.addParameter('useTraps',[],@(x) isempty(x) || (isscalar(x) && islogical(x)));
ip.addParameter('rotation',[],@(x) isempty(x) || (isscalar(x) && isnumeric(x)));
ip.parse(varargin{:});

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

expName = ip.Results.expName;
if isempty(expName)
    if ~isa(dataset,'OmeroDataset')
        dataset = OmeroDataset(dataset,'use_omero_cache',true);
    end
    expName = dataset.getValidExpNameGUI;
end

% Create a new cExperiment from the Omero dataset
cExperiment=babyExperimentOmero(dataset,expName);
cExperiment.segmentationSource='Omero';%Data will be retrieved from the Omero database for segmentation

% createTimelapsePositions given with explicit arguments so that all
% positions are loaded.
cExperiment.createTimelapsePositions(ip.Results.trapChannel,'all',...
    ip.Results.pixelSize,ip.Results.rotation,[],ip.Results.useTraps);

% Upload the cExperiment file to the database
cExperiment.saveExperiment;
end