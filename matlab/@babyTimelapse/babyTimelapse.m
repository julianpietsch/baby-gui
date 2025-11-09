classdef babyTimelapse<MemAware
    % BABYTIMELAPSE This is the base class used for processing a time-lapse of
    % images in the BABY GUI. BABYEXPERIMENT organises multiple BABYTIMELAPSE
    % objects, typically one for each position/point in an experiment. It
    % stores the location of the images for each channel at each timepoint,
    % processing parameters and all the segmentation results (which are
    % then compiled together in BABYEXPERIMENT). 
    % Each timepoint is stored as an entry in the cTimepoint structure
    % array. This contains:
    %
    % filename:      cell array of the names of the files associated with
    %                that timepoint
    %
    % after processing it also stores:
    %
    % trapLocations: the location of the traps. 
    % trapInfo:      a structure that holds all the information about the
    %                location, label and outline of each cell in each trap
    %
    % See also BABYEXPERIMENT
    
    properties
        timelapseDir %Location in which files are stored. set to 'ignore' to use absolute file names
        cTimepoint  %structure array. This contains:
                    %
                    %filename:      cell array of the names of the files associated with that timepoint
                    %
                    %after processing it also stores:
                    %
                    %trapLocations: the location of the traps.
                    %trapInfo:      a structure that holds all the information about the
                    %               location, label and outline of each cell in each trap
                    %trapMaxCell :  the maximum cell label for the cells in that trap
                    %

        pixelSize = 0.263 % the real size of pixels in the image. default of 0.263 is for swainlab microscopes at 60x.
        image_flipud=false % whether to flip the image
        image_rotation % rotation only recommended in multiples of 90 degrees
        rawImSize % size of the images before any rescaling or rotation.
        scaledImSize % size of images after rescaling but before rotation.
        imSize %the size of the images after rotation.
                          
        trapsPresent % a boolean whether traps are present or not in the image
        trapTemplates % templates and associated info for identifying traps; empty if there are no traps
        trapTemplateChannel % channel used in cTrapSelectDisplay to identify traps
        
        cellsToPlot %Array indicating which cells to extract data for. row = trap num, col is cell tracking number
        timepointsProcessed %a logical array of timepoints which have been processed
        timepointsToProcess %list of timepoints that should be processed (i.e. checked for cells and what not)
        extractedData %a structure array of sparse arrays of extracted statistics for the cell. One entry for each channel extracted, one row in the array for each cell with zeros where the cell is absent.
        channelNames % cell array of channel names. Used in selecting the file(s) to load when a particular channel is selected in returnImage method

        cellMothers = sparse(100,1e3)
        offset = [0 0] %a n x 2 offset of each channel compared to DIC. So [0 0; x1 y1; x2 y2]. Positive shifts left/up.
        BackgroundCorrection = {[]} %correction matrix for image channels. If non empty, returnSingleTimepoint will '.multiply' the image by this matrix.
                                    %this is applied BEFORE rescaling, on the grounds that the background correction is generally found from the raw images. 
        BackgroundOffset = {[]} %scalar offset to be used with BackgroundCorrection matrix. If non empty, returnSingleTimepoint will subtract this offset 
                                %before multiplying by the correction matrix, then add it back after applying
                                %the correction matrix. This is to stop the flatfield correction inflating the noise
                                %where readings are low.
                                %Adding back is mostly to keep things compatible with older data. 
        extractionParameters = babyTimelapse.defaultExtractParameters;
        %parameters for the extraction of cell Data, a function handle and
        %a parameter structure which the function makes use of.
        
        metadata = []
    
       
    end
    
    properties(Dependent = true)
        % not real properties,calculated from other data.
        defaultTrapDataTemplate % a sparse array of the right size for holding semgmentation data
        cellInfoTemplate % template for the cellInfo structure.
        trapInfoTemplate % template for the trapInfo structure
        cTimepointTemplate % template for the cTimepoint structure.
        cTrapSize % struct specifying bounding-box half-width and half-height (bb_width & bb_height); derived from size of trap templates
        trapImSize % uses the cTrapSize property to give the size of the image.
        nzstacks % a Map specifying the number of z stacks for each channel
        babyBrain % options for and status of segmentation by the baby
    end
    
    properties (Access=private)
        babyBrain_val
    end
    
    properties(SetAccess = immutable)
        
    end
    
    properties(Constant)
        defaultExtractParameters = struct('extractFunction',@extractCellDataStandardParfor,...
            'functionParameters',struct('type','max','channels','all','nuclearMarkerChannel',NaN,'maxPixOverlap',5,'maxAllowedOverlap',25));
    end
    
    properties (Transient)
        % Transient properties won't be saved
        logger % *optional* handle to an babyLogging object to keep a log
        imcache % *optional* handle to an ImageCache object for caching images in memory and on disk
        temporaryImageStorage=struct('channel',-1,'images',[]); %this is to store the loaded images from a single channel (ie BF) into memory
        %This allows the cell tracking and curating things to happen a
        %whole lot faster and easier. This way you can just modify the
        %returnTimepoint file to check to see if something is loaded.
        % - channel
        % - images
    end
    
   
    
    events
        LogMsg
        TimepointChanged
        TrapChanged
    end
    
    methods
        
        function cTimelapse=babyTimelapse(folder,varargin)
            % cTimelapse=babyTimelapse(folder,varargin)
            % instantiate a babyTimelapse object from a folder containing
            % images. If folder is empty it is requested by uigetdir, and
            % it becomes the timelapseDir.
            % 
            % varargin{1} can be a logical that will make the constructor
            % run nothing if it is true. this was done to be able to write
            % nice load functions.
            %
            % Most of the actual setting up is done by
            % BABYTIMELAPSE.LOADTIMELAPSE
            %
            % See also, BABYTIMELAPSE.LOADTIMELAPSE
            
            if nargin>=2 && islogical(varargin{1})
                NoAction = varargin{1};
            else
                NoAction = false;
            end
            
            if ~NoAction
                if nargin<1 || isempty(folder)
                    folder=uigetdir(pwd,'Select the folder containing the images associated with this timelapse');
                    fprintf('\n    Select the folder containing the images associated with this timelapse\n');
                end
                cTimelapse.timelapseDir=folder;
                cTimelapse.cellsToPlot=sparse(100,1e3);
            end
        end
            
        function name = getName(cTimelapse)
            % name = getName(cTimelapse)
            % sometimes you want to have an identifiable name for a
            % babyTimelapse object, for figure names and such.
            try
                if strcmp(cTimelapse.timelapseDir,'ignore')
                    name = cTimelapse.cTimepoint(1).filename{1};
                else
                    name = [cTimelapse.timelapseDir '/'];
                end
                
                % get section of this path between second to last /|\ and
                % last /|\ in a reasonably robust way.
                locs = regexp(name,'[\\|/]','start');
                name = name(max(locs(max(length(locs)-2,1))+1,1):max(locs(end)-1,1));
                
            catch
                name = [];
            end
        end
        
        
        function cTimelapseOUT = copy(cTimelapseIN)
        %cTimelapseOUT = copy(cTimelapseIN)
        % make a new cTimelapse object with all the same field values. 
            cTimelapseOUT = babyTimelapse([],true);
            
            FieldNames = fields(cTimelapseIN);
            
            for i = 1:numel(FieldNames)
                m = findprop(cTimelapseIN,FieldNames{i});
                if ~ismember(m.SetAccess,{'immutable','none'}) || m.Dependent
                    cTimelapseOUT.(FieldNames{i}) = cTimelapseIN.(FieldNames{i});
                end
            end

            
        end
        
        function trapInfo_struct = createTrapInfoTemplate(cTimelapse,data_template,avg_ncells)
            % trapInfo_struct =
            % createTrapInfoTemplate(cTimelapse,data_template)
            %
            % create strandard empty trapInfo structure for use in
            % intialising trapInfo.
            %
            % data template is optional. should be a sparse array. If not
            % it throws an error. default is spares of size cTrapSize. If
            % this is empty it just uses an empty array.
            
            if nargin<3, avg_ncells = 16; end
            if nargin<2, data_template = []; end
            
            if isempty(data_template)
                data_template = cTimelapse.defaultTrapDataTemplate;
                % By default we expect a low number of cells per trap; use
                % the avg_ncells argument to increase nzmax
                centers_template = sparse([],[],false(0),...
                    size(data_template,1),size(data_template,2),avg_ncells);
            elseif issparse(data_template)
                centers_template = data_template;
            else
                error('data_template should be a sparse array');
            end
            
            trapInfo_struct = cTimelapse.trapInfoTemplate;
            trapInfo_struct.cell = cTimelapse.cellInfoTemplate;
            trapInfo_struct.segCenters = centers_template;
            trapInfo_struct.cell.segmented = data_template;
            
        end
        
        function default_trap_indices = defaultTrapIndices(cTimelapse,tp)
            % default_trap_indices = defaultTrapIndices(cTimelapse,tp=1)
            % return the default trap indices to run anything over.
            if nargin<2
                tp = cTimelapse.timepointsToProcess(1);
            end
            default_trap_indices = 1:length(cTimelapse.cTimepoint(tp).trapInfo);
        end
        
        function data_template = get.defaultTrapDataTemplate(cTimelapse)
            % data_template = defaultTrapDataTemplate(cTimelapse)
            % returns a sparse array of the default size for populating
            % cell and trapInfo structures. Used at various points in the
            % code where these things need to be populated.
            % for trap containing cTimelapses, this is the trapSize.
            % for those without traps, it is the image size.
            data_template_size = cTimelapse.trapImSize;
            if ~isempty(data_template_size)
                data_template = sparse([],[],false(0),...
                    data_template_size(1),data_template_size(2),...
                    min(data_template_size)); % nzmax
                % i.e., expect average radius of cell to be about 1/3 the
                % size of the smallest trap dimension
            else
                data_template = sparse([],[],false(0));
            end
            
        end

        function val = get.cTrapSize(cTimelapse)
            if ~isempty(cTimelapse.trapsPresent) && cTimelapse.trapsPresent
                assert(isstruct(cTimelapse.trapTemplates) && ...
                    isfield(cTimelapse.trapTemplates,'positiveExamples'),...
                    'Cannot identify traps, no trap templates are available');
    
                positiveExamples = cTimelapse.trapTemplates.positiveExamples;
                val = struct();
                val.bb_width = ceil((size(positiveExamples,2)-1)/2);
                val.bb_height = ceil((size(positiveExamples,1)-1)/2);
            else
                val = [];
            end
        end

        function trapImSize = get.trapImSize(cTimelapse)
            % size of the trap if traps present or imSize if not.
            trapImSize = [];
            if ~isempty(cTimelapse.trapsPresent)
                if cTimelapse.trapsPresent &&  ~isempty(cTimelapse.cTrapSize)
                    trapImSize = 2*[cTimelapse.cTrapSize.bb_height cTimelapse.cTrapSize.bb_width] + 1;
                elseif   ~cTimelapse.trapsPresent &&  ~isempty(cTimelapse.imSize)
                    trapImSize = cTimelapse.imSize;
                end
            end
        end
        
        function set.trapImSize(~,~)    
            % do nothing, just to stop errors
        end
        
        function set.defaultTrapDataTemplate(~,~)    
            % do nothing, just to stop errors
        end

        function set.cellInfoTemplate(~,~)
            % do nothing, just to stop errors
        end
        
        function set.trapInfoTemplate(~,~)
            % do nothing, just to stop errors
        end
        function set.cTimepointTemplate(~,~)
            % do nothing, just to stop errors
        end
        
        function cellInfoTemplate = get.cellInfoTemplate(cTimelapse)
             segTemplate = cTimelapse.defaultTrapDataTemplate;
            
             cellInfoTemplate = struct('cellCenter',[],...
                       'cellRadius',[],...
                       'segmented',segTemplate,...
                       'cellRadii',[],...
                       'cellAngle',[]);
        end
        
        function trapInfoTemplate = get.trapInfoTemplate(cTimelapse)
            trapInfoTemplate = struct('segCenters',[],...
                'cell',cTimelapse.cellInfoTemplate, ...
                'cellsPresent',0,'cellLabel',[]);
        end
        
        function cTimepointTemplate =get.cTimepointTemplate(cTimelapse)
            cTimepointTemplate = struct('filename',[],'trapLocations',[],...
                'trapInfo',cTimelapse.trapInfoTemplate,'trapMaxCell',[]);
        end
        
        function val = get.logger(cTimelapse)
            val = cTimelapse.logger;
        end
        
        function val = get.nzstacks(cTimelapse)
            val = containers.Map;
            meta = cTimelapse.metadata;
            if isempty(meta) || ...
                    ~all(isfield(meta,{'acq','posname'})) || ...
                    ~all(isfield(meta.acq,{'channels',...
                    'zsections','positions'}))
                return;
            end
            
            acq = meta.acq;
            chnames = acq.channels.names;
            if size(acq.positions,1) > 0
                posind = strcmp(acq.positions.name,meta.posname);
                if sum(posind) ~= 1 && isprop(cTimelapse,'reader')
                    posind = strcmp(acq.positions.name,cTimelapse.reader.posName);
                end
                active = table2array(acq.positions(posind,acq.channels.names)) ~= 0;
            else
                active = true(size(acq.channels.names));
            end
            chnames = chnames(active);
            
            nz = acq.zsections.sections;
            for c=1:numel(chnames)
                chname = chnames{c};
                if acq.channels.zsect(strcmp(acq.channels.names,chname))
                    val(chname) = nz;
                else
                    val(chname) = 1;
                end
            end
        end
        
        function val = get.babyBrain(cTimelapse)
            if isempty(cTimelapse.babyBrain_val)
                cTimelapse.babyBrain_val = BabyBrain;
                cTimelapse.babyBrain_val.set_channel_nzstacks(...
                    cTimelapse.nzstacks);
            end
            val = cTimelapse.babyBrain_val;
        end
        
        function set.babyBrain(cTimelapse,val)
            assert(isa(val,'BabyBrain'),...
                '"babyBrain" must be specified as a "BabyBrain" object');
            % Copy only the configuration
            cTimelapse.babyBrain.copy_channel_map(val);
            cTimelapse.babyBrain.config = val.config;
        end
        
        function resetBrain(cTimelapse)
            cTimelapse.babyBrain_val = [];
        end
    end
    
    methods (Access={?babyTimelapse,?babyTimelapseOmero,...
            ?OmeroDatabase,?babyExperiment,?babyExperimentOmero,...
            ?babyTimelapseSamples,?babyTimelapseBioformats,...
            ?babyTimelapseFileSeries,?babyTimelapseOmero})
        function propNames = copyprops(cTimelapse,TemplateTimelapse,omit)
            %COPYPROPS Copy all properties from a cTimelapse into this one
            %   This function can copy both public and private properties.
            %   Use OMIT to specify a cellstr of properties that will not
            %   be copied. This function gets used in the loadobj method,
            %   by OmeroDatabase.convertSegmented and by 
            %   babyExperimentOmero.convertToFolderExperiment
            
            if nargin<3 || isempty(omit), omit = {}; end
            if ~iscellstr(omit)
                error('The "omit" argument must be a cellstr.');
            end
            
            % Only populate copyable fields occuring in both this object
            % and the template object:
            propNames = intersect(...
                getCopyableProperties(cTimelapse,'babyTimelapse'),...
                getCopyableProperties(TemplateTimelapse,'babyTimelapse'));
            % Omit requested properties
            propNames = setdiff(propNames,omit);
            
            % Copy all properties/fields to this cTimelapse:
            for f = 1:numel(propNames)
                cTimelapse.(propNames{f}) = TemplateTimelapse.(propNames{f});
            end
        end
    end
    
    methods (Static)
        function cTimelapse = loadobj(load_structure)
            % cTimelapse = loadobj(load_structure)
            % load function to help maintain back compatability and take
            % care of fiddly loading behaviour
            
            if isstruct(load_structure)
                % The 'load_structure' argument could be for this class or 
                % any of its children. We work this out based on the 
                % 'load_class_name' property set by the children.
                if isfield(load_structure,'load_class_name')
                    lcname = load_structure.load_class_name;
                    assert(startsWith(lcname,'babyTimelapse') ...
                        && exist(lcname,'class')==8,'bad load_class_name!');
                    cTimelapse = eval([lcname,'([],true)']);
                    load_structure = rmfield(load_structure,'load_class_name');
                else
                    % Default to instantiating according to this class
                    cTimelapse = babyTimelapse([],true);
                end
            else
                cTimelapse = eval([class(load_structure),'([],true)']);
            end
                
            cTimelapse.copyprops(load_structure);
            
            % Back compatibility checks and what not
            %when a new field is added this load operation should be
            %updated to populate the field appropriately and maintain back
            %compatibility.
            
            if isempty(cTimelapse.timepointsToProcess)
                
                cTimelapse.timepointsToProcess = 1:length(cTimelapse.cTimepoint);
                
            end
            
            if length(cTimelapse.BackgroundCorrection)==1 && isempty(cTimelapse.BackgroundCorrection{1})
                cTimelapse.BackgroundCorrection = {};
                cTimelapse.BackgroundCorrection(1:length(cTimelapse.channelNames)) = {[]};
            end
            
            if size(cTimelapse.offset,1)<length(cTimelapse.channelNames)
                cTimelapse.offset(end+1:length(cTimelapse.channelNames),:) = 0;
            end
            
            if isempty(cTimelapse.scaledImSize)
                cTimelapse.scaledImSize = cTimelapse.imSize;
            end

        end
        
        function cTimelapse_save = saveobj(cTimelapse_in)
            
            cTimelapse_save = cTimelapse_in.copy;
            
        end
        
        cTimelapse = loadFrom(cTimelapse_filepath)
    end
end

