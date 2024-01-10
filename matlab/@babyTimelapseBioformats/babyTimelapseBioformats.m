classdef babyTimelapseBioformats<babyTimelapse
    % BABYTIMELAPSEBIOFORMATS A subclass of babyTimelapse to process 
    % any image formats that can be loaded by OME BioFormats.
    % See also BABYEXPERIMENTBIOFORMATS, BABYTIMELAPSE
    
    properties
        posNum % index to the sequence for this cTimelapse in the imageFile
        stackingArgs % arguments to pass to the reader if files are stacked
    end
    
    properties (Dependent)
        reader
    end
    
    properties (Access=private,Transient)
        reader_val
    end
    
    methods
        

        function cTimelapse=babyTimelapseBioformats(imgfile,varargin)
            % cTimelapse=babyTimelapseBioformats(imgfile,varargin)
            % Instantiate cTimelapse from an image file using OME Bioformats.
            %
            % varargin{1} can be a logical that will make the constructor
            % run nothing if it is true. this was done to be able to write
            % nice load functions.
            %
            % Most of the actual setting up is done by
            % BABYTIMELAPSEBIOFORMATS.LOADTIMELAPSE
            %
            % See also, BABYTIMELAPSEBIOFORMATS.LOADTIMELAPSE
            
            % call babyTimelapse constructor as though loading (i.e. to
            % make a bare object).
            cTimelapse@babyTimelapse([],true);
            
            if nargin>=2 && islogical(varargin{1})
                return
            end
            
            cTimelapse.timelapseDir = imgfile;
            cTimelapse.posNum = varargin{1};
            if nargin>2
                cTimelapse.stackingArgs = varargin{2};
            end
            
            if nargin<4 || isempty(varargin{3})
                channelNames = cTimelapse.reader.channels;
            else
                channelNames = varargin{3};
            end
            cTimelapse.channelNames = channelNames;
            cTimelapse.cellsToPlot = sparse(100,1e3);
        end
        
        function name = getName(cTimelapse)
            % name = getName(cTimelapse)
            % sometimes you want to have an identifiable name for a
            % babyTimelapse object, for figure names and such.
            name = char(cTimelapse.reader.posName);
        end
        
        function val = get.reader(cTimelapse)
            if isempty(cTimelapse.reader_val) && isstruct(cTimelapse.stackingArgs)
                fvals = struct2cell(cTimelapse.stackingArgs);
                if ~all(cellfun(@isempty,fvals))
                    stackargs = [fieldnames(cTimelapse.stackingArgs),fvals]';
                    cTimelapse.reader_val = ImageReaderBioformatsStack(...
                        cTimelapse.timelapseDir,stackargs{:});
                end
            end
            if isempty(cTimelapse.reader_val)
                cTimelapse.reader_val = ImageReaderBioformats(cTimelapse.timelapseDir);
            end
            val = cTimelapse.reader_val;
            val.pos = cTimelapse.posNum;
        end
        
        function releaseReader(cTimelapse)
            delete(cTimelapse.reader_val);
            cTimelapse.reader_val = [];
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
                getCopyableProperties(cTimelapse,'babyTimelapseBioformats'),...
                getCopyableProperties(TemplateTimelapse,'babyTimelapseBioformats'));
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
                load_structure.load_class_name = 'babyTimelapseBioformats';
            end
            cTimelapse = loadobj@babyTimelapse(load_structure);
        end

        function cTimelapse_save = saveobj(cTimelapse_in)
            % cTimelapse_save = saveobj(cTimelapse_in)
            % currently just runs BABYTIMELAPSE saveobj method, but could
            % be modified to do different things.
            
            cTimelapse_save = saveobj@babyTimelapse(cTimelapse_in);
            
        end
    end
end
