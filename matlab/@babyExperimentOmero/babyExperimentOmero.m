classdef babyExperimentOmero < babyExperiment
    %BABYEXPERIMENTOMERO Subclass of babyExperiment for use with OMERO
    %
    %   babyExperimentOmero(DATASET) creates a new experiment for the
    %   OmeroDataset DATASET with default experiment ID of '001'.
    %   
    %   babyExperimentOmero(DATASET,EXPNAME) creates a new experiment with
    %   ID given by the char vector EXPNAME.
    %   
    %   babyExperimentOmero(true) creates a 'bare' experiment with no
    %   defined fields (primarily for use by loadobj method).
    %
    % Code originally written by Ivan Clark and edited by Elco Bakker.
    % Subsequently heavily modified by Julian Pietsch.
    %
    % See also BABYEXPERIMENT, BABYTIMELAPSEOMERO
    
    properties
        dataset % OmeroDataset object
        fileAnnotation_id %id number of the file annotation used to save the babyExperimentOmero object in the database
        logFileAnnotation_id %id number of the file annotation used to save the babyExperimentOmero log file in the database
        segmentationSource='' %Flag to determine where the source data was obtained during segmentation, can be 'Omero', 'Folder' or empty (if not segmented). Data segemented from a file folder must be flipped both vertically and horizontally to match the segmentation results
        archivedChannelNames = {} %Channel names from folder cExperiment; empty unless converted from folder cExperiment
    end
    
    properties (Dependent)
        server % OmeroServer object
    end
    
    methods
        function cExperiment = babyExperimentOmero(dataset,expName)
            %babyExperimentOmero Create a new babyExperimentOmero object
             
            % call super class constructor such that it does not initialise
            % anything.
            cExperiment@babyExperiment(true);
            
            % If initialised with single input true, make an empty object
            % (used in load function)
            if nargin==1 && islogical(dataset) && dataset
                % if folder is true, cExperiment returned bare for load
                % function.
                return
            end
            
            if ~isa(dataset,'OmeroDataset')
                dataset = OmeroDataset(dataset,'use_omero_cache',true);
            end
            %Define the Omero properties of cExperiment
            cExperiment.dataset = dataset;
            server = cExperiment.server;
            
            % Default to enabling logging:
            cExperiment.shouldLog=true;
            % NB: constructor no longer needs to instantiate the logger
            % property, see babyExperiment
            
            %Experiment is being initialized from an Omero dataset - folder is an omero.model.DatasetI object
            if nargin<3 || isempty(expName)
                expName = '001';
            end
            if iscell(expName)
                expName=expName{:};
            end
            cExperiment.rootFolder = expName;
            cExperiment.saveFolder = server.downloadDir(dataset);
            
            %Record the user who is creating the cExperiment
            if ispc
                cExperiment.creator=getenv('USERNAME');
            else
                [~, cExperiment.creator] = system('whoami');
            end
            %Initialize records of positions segmented and tracked
            cExperiment.posTrapsTracked=0;
            cExperiment.posSegmented=0;
            cExperiment.posTracked=0;
            
            cExperiment.dirs=dataset.posNames;
            
            cExperiment.cellsToPlot=cell(1);
            
            % Take meta data from dataset object
            cExperiment.metadata = dataset.meta;
            if isfield(cExperiment.metadata,'logfiledata')
                logfiledata = cExperiment.metadata.logfiledata;
                logfields = fieldnames(logfiledata);
                for f=1:numel(logfields)
                    if ~isfield(cExperiment.metadata,logfields{f})
                        cExperiment.metadata.(logfields{f}) = logfiledata.(logfields{f});
                    end
                end
            end
            if isfield(cExperiment.metadata,'acqfiledata')
                cExperiment.metadata.acq = cExperiment.metadata.acqfiledata;
                cExperiment.metadata = rmfield(cExperiment.metadata,'acqfiledata');
            end
            
            nposes = numel(dataset.posNames);
            allBaseChannelNames = cell(nposes,1);
            allBaseChannelHasStacks = cell(nposes,1);
            maxnz = 1;
            for p=1:nposes
                dataset.posNum = p;
                allBaseChannelNames{p} = dataset.channelNames;
                allBaseChannelHasStacks{p} = dataset.hasZstacks;
                maxnz = max(maxnz,dataset.imageSize(3));
            end
            allBaseChannelNames = vertcat(allBaseChannelNames{:});
            allBaseChannelHasStacks = vertcat(allBaseChannelHasStacks{:});
            allBaseChannelHasStacks = unique(allBaseChannelNames(allBaseChannelHasStacks));
            allBaseChannelNames = unique(allBaseChannelNames,'stable');
            channelNames = cell(numel(allBaseChannelNames),1);
            for c=1:numel(channelNames)
                if ismember(allBaseChannelNames{c},allBaseChannelHasStacks)
                    channelNames{c} = cell(maxnz,1);
                    for z=1:maxnz
                         channelNames{c}{z} = sprintf('%s_%03u',allBaseChannelNames{c},z);
                    end
                end
            end
            
            cExperiment.channelNames = vertcat(allBaseChannelNames{:},channelNames{:});
            
        end
        
        function val = get.server(this)
            val = this.dataset.database;
        end
    end
    
    methods (Access={?babyExperiment,?babyExperimentOmero,...
            ?OmeroDatabase,?babyExperimentSamples,...
            ?babyExperimentBioformats,?babyExperimentFileSeries,...
            ?babyExperimentOmero})
        function propNames = copyprops(cExperiment,TemplateExperiment,omit)
            %COPYPROPS Copy all properties from a cExperiment into this one
            %   This function can copy both public and private properties.
            %   Use OMIT to specify a cellstr of properties that will not
            %   be copied. This function gets used in the loadobj method
            %   and also by the convertSegmented method of the 
            %   OmeroDatabase class.
            
            if nargin<3 || isempty(omit), omit = {}; end
            if ~iscellstr(omit)
                error('The "omit" argument must be a cellstr.');
            end
            
            copied = copyprops@babyExperiment(cExperiment,TemplateExperiment,omit);
            
            % Only populate copyable fields occuring in both this object
            % and the template object:
            propNames = intersect(...
                getCopyableProperties(cExperiment,'babyExperimentOmero'),...
                getCopyableProperties(TemplateExperiment,'babyExperimentOmero'));
            % Omit requested properties:
            propNames = setdiff(propNames,omit);
            
            % Omit properties copied by parent:
            copyNames = setdiff(propNames,copied);
            
            % Copy all properties/fields to this cExperiment:
            for f = 1:numel(copyNames)
                cExperiment.(copyNames{f}) = TemplateExperiment.(copyNames{f});
            end
        end
    end

    methods (Static)
        function cExperiment = loadobj(load_structure)
            % Store the class name in the load_structure before passing to
            % parent load routine
            if isstruct(load_structure)
                load_structure.load_class_name = 'babyExperimentOmero';
            end
            cExperiment = loadobj@babyExperiment(load_structure);
        end
    end

    methods (Access=protected)
        function cTimelapse = newTimelapse(cExperiment,pos)
            ds = cExperiment.dataset;
            assert(all(ismember(cExperiment.dirs,ds.posNames)),...
                'cExperiment.dirs specifies positions that do not exist in dataset');
            ds.pos = cExperiment.dirs{pos};
            cTimelapse = babyTimelapseOmero(ds,cExperiment.channelNames);
        end
    end
end