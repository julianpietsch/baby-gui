classdef babyTimelapseFileSeries<babyTimelapse
    % BABYTIMELAPSEFILESERIES A subclass of babyTimelapse to process 
    % time series built from series of image files.
    % See also BABYEXPERIMENTFILESERIES, BABYTIMELAPSE
    
    properties
        posNum % index to the position for this cTimelapse in the reader
        readerArgs
    end
    
    properties (Dependent)
        reader
    end
    
    properties (Access=private,Transient)
        reader_val
    end
    
    methods
        

        function cTimelapse=babyTimelapseFileSeries(rootfolder,varargin)
            % cTimelapse=babyTimelapseFileSeries(rootfolder,varargin)
            % Instantiate cTimelapse from an image file series.
            %
            % varargin{1} can be a logical that will make the constructor
            % run nothing if it is true. this was done to be able to write
            % nice load functions.
            %
            % Most of the actual setting up is done by
            % BABYTIMELAPSEFILESERIES.LOADTIMELAPSE
            %
            % See also, BABYTIMELAPSEFILESERIES.LOADTIMELAPSE
            
            NoAction = false;
            if nargin>=2 && islogical(varargin{1})
                NoAction = varargin{1};
            end
            
            % call babyTimelapse constructor as though loading (i.e. to
            % make a bare object).
            cTimelapse@babyTimelapse([],true);
            
            if NoAction, return; end

            ip = inputParser;
            ip.addRequired('rootfolder',@(x) ischar(x) && isrow(x));
            ip.addRequired('posnum',@(x) isscalar(x) && isnumeric(x));
            ip.addParameter('channelNames',[],@(x) isempty(x) || (iscellstr(x) || isstring));
            ip.KeepUnmatched = true;
            ip.parse(rootfolder,varargin{:});
            
            cTimelapse.readerArgs = [fieldnames(ip.Unmatched),struct2cell(ip.Unmatched)]';
            
            cTimelapse.timelapseDir = rootfolder;
            cTimelapse.posNum = ip.Results.posnum;
            channelNames = ip.Results.channelNames;
            if isempty(channelNames)
                channelNames = cTimelapse.reader.channels;
            else
                channelNames = cellstr(channelNames);
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
            if isempty(cTimelapse.reader_val)
                cTimelapse.reader_val = ImageReaderFileSeries(...
                    cTimelapse.timelapseDir,cTimelapse.readerArgs{:});
            end
            val = cTimelapse.reader_val;
            val.pos = cTimelapse.posNum;
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
                getCopyableProperties(cTimelapse,'babyTimelapseFileSeries'),...
                getCopyableProperties(TemplateTimelapse,'babyTimelapseFileSeries'));
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
                load_structure.load_class_name = 'babyTimelapseFileSeries';
            end
            cTimelapse = loadobj@babyTimelapse(load_structure);
        end
    end
end
