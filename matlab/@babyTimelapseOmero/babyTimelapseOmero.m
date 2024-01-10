classdef babyTimelapseOmero<babyTimelapse
    % BABYTIMELAPSEOMERO A subclass of babyTimelapse to process images
    % coming from an OMERO Database.
    % See also BABYEXPERIMENTOMERO, BABYTIMELAPSE
    
    properties
        dataset % an OmeroDataset with the position correctly set
        posNum % a copy of the posNum
        segmentationSource='';%Flag to determine where the source data was obtained during segmentation, can be 'Omero', 'Folder' or empty (if not segmented). Data segemented from a file folder must be flipped both vertically and horizontally to match the segmentation results
        archivedChannelNames={};%Channel names from folder cTimelapse; empty unless converted from folder cTimelapse
        fileAnnotation_id%id number of the file annotation used to save the babyTimelapseOmero object in the database
        flipchannels=[]%logical array specifying which channels should be flipped in returnSingleTimepointRaw
    end
    properties (Dependent)
        server
    end
    
    methods
        function cTimelapse = babyTimelapseOmero(dataset,channelNames)
            %BABYTIMELAPSEOMERO instantiate cTimelapse from OmeroDataset
            %
            %   The position is inferred from the state of the OmeroDataset
            %   when provided here. For consistency of channel number 
            %   across all positions, the channel names can optionally 
            %   be provided.
            %
            % Most of the actual setting up is done by
            % BABYTIMELAPSEOMERO.LOADTIMELAPSE
            %
            % See also, BABYTIMELAPSEOMERO.LOADTIMELAPSE
            
            % call babyTimelapse constructor as though loading (i.e. to
            % make a bare object).
            cTimelapse@babyTimelapse([],true);
            
            if nargin>=2 && islogical(channelNames) && channelNames
                return
            end
            
            if nargin<2 || isempty(channelNames)
                channelNames = dataset.channelNames;
            end
            
            % Infer the position number from the current posNum of the
            % provided dataset:
            cTimelapse.posNum = dataset.posNum;
            
            % Ensure that we work with a copy of the dataset
            cTimelapse.dataset = OmeroDataset(dataset.dataset,...
                'meta',dataset.meta,'use_omero_cache',true);
            cTimelapse.dataset.posNum = cTimelapse.posNum;
            
            cTimelapse.channelNames = channelNames;
            
            cTimelapse.cellsToPlot=sparse(100,1e3);
        end
        
        function name = getName(cTimelapse)
            % name = getName(cTimelapse)
            % sometimes you want to have an identifiable name for a
            % babyTimelapse object, for figure names and such.
            name = cTimelapse.dataset.name;
        end
        
        function val = get.server(cTimelapse)
            val = cTimelapse.dataset.database;
        end
    end
    
    methods (Access={?babyTimelapse,?babyTimelapseOmero,...
            ?OmeroDatabase,?babyExperiment,?babyExperimentOmero,...
            ?babyTimelapseSamples,?babyTimelapseBioformats,...
            ?babyTimelapseFileSeries,?babyTimelapseOmero})
        function propNames = copyprops(cTimelapse,TemplateTimelapse,omit)
            %COPYPROPS Copy all properties from cTimelapse into this one
            %   This function can copy both public and private properties.
            %   Use OMIT to specify a cellstr of properties that will not
            %   be copied. This function gets used in the loadobj method
            %   and also by the convertSegmented method of the 
            %   OmeroDatabase class.
            
            if nargin<3 || isempty(omit), omit = {}; end
            if ~iscellstr(omit)
                error('The "omit" argument must be a cellstr.');
            end
            
            copied = copyprops@babyTimelapse(cTimelapse,TemplateTimelapse,omit);
            
            % Only populate copyable fields occuring in both this object
            % and the template object:
            propNames = intersect(...
                getCopyableProperties(cTimelapse,'babyTimelapseOmero'),...
                getCopyableProperties(TemplateTimelapse,'babyTimelapseOmero'));
            % Omit properties copied by parent:
            propNames = setdiff(propNames,copied);
            % Omit requested properties:
            propNames = setdiff(propNames,omit);
            
            % Copy all properties/fields to this cTimelapse:
            for f = 1:numel(propNames)
                cTimelapse.(propNames{f}) = TemplateTimelapse.(propNames{f});
            end
        end
    end

    methods (Static)
        function cTimelapse = loadobj(load_structure)
            % Store the class name in the load_structure before passing to
            % parent load routine
            if isstruct(load_structure)
                load_structure.load_class_name = 'babyTimelapseOmero';
            end
            cTimelapse = loadobj@babyTimelapse(load_structure);
        end
    end
end
