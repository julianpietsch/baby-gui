classdef ImageCache < MemAware
    %ImageCache Handles caching and retrieval of trap images
    %   This class provides an interface to access trap image time series
    %   derived from a babyExperiment object.
    %
    % TODO: include filters to exclude traps/positions that do not
    % contain cells / cells of interest
    % TODO: loop over cTimelapses to fill the trapMaps
    
    
    properties (Transient)
        channels = {} % Names of the channels accessible from this cache
        posNames = {} % Names of positions accessible from this cache
        trapMaps = {}
    end
    
    properties (Dependent)
        savedir
        cExperiment
        progress_bar
        storage_class
        storage_model
        stackSize % standard dimensions of each trap stack for current pos
        stackSizes
        timepointIndices
                
        % The cache structure is stored using the following properties
        posMap
        channelMap
    end
    
    properties (Access=private)
        savedir_val = []
        storage_class_val = 'uint8'
        storage_model_val = ''
    end
    
    properties (Transient,Access=private)
        imcache = {}
        imrange = {}
        loadedtps = {}
        cExperiment_val = []
        progress_bar_val = []
        current_position = []
        cTimelapse = []
        local_cache_id = ''
        stackSizes_val = {}
        timepointIndices_val = {}
        posMap_val = []
        channelMap_val = []
    end
    
    properties (Constant)
        SegmentationImageChannels = {'segoutline','segarea'}
        
        % Error message types that get raised by this class
        Errors = struct(...
            'BadParams','BABY:ImageCache:badParams',...
            'BadImage','BABY:ImageCache:badImage',...
            'BadTimelapse','BABY:ImageCache:badTimelapse',...
            'NoExperiment','BABY:ImageCache:noExperiment',...
            'NotSaved','BABY:ImageCache:notSaved',...
            'NotAvailable','BABY:ImageCache:notAvailable',...
            'WrongType','BABY:ImageCache:wrongType');
    end
    
    events
        InfoUpdated % Event triggered when ImageCache has been modified; notifies parent objects to autosave
    end
    
    methods
        function this = ImageCache(cExperiment,varargin)
            %ImageCache Handles caching and retrieval of trap images
            ip = inputParser;
            ip.addRequired('cExperiment',@(x) isa(x,'babyExperiment') && isscalar(x));
            ip.addParameter('LocalCachePath','',@(x) ischar(x));
            ip.addParameter('storage_class','uint8',@(x) ischar(x) && isvector(x));
            ip.addParameter('storage_model','',@(x) isempty(x) || (ischar(x) && isvector(x)));
            ip.parse(cExperiment,varargin{:});
            
            this.cExperiment = cExperiment;
            
            cachepath = ip.Results.LocalCachePath;
            if ~isempty(cachepath)
                this.openLocalCache(cachepath);
            end
        end
        
        function clearcache(this)
            this.imcache = {};
        end
        
        function sz = cachesize(this)
            multipliers = struct('uint8',8,'uint16',16,'double',64);
            sz = prod(this.stackSize(:)) * multipliers.(this.storage_class)...
                * sum(cellfun(@length,this.imcache(:)));
        end
        
        function val = get.stackSize(this)
            val = [];
            if ~isempty(this.current_position)
                val = this.stackSizes{this.current_position};
            end
        end
        
        function val = get.stackSizes(this)
            if isempty(this.stackSizes_val)
                this.stackSizes_val = cell(1,numel(this.cExperiment.dirs));
            end
            val = this.stackSizes_val;
        end
        function set.stackSizes(this,val)
            assert(numel(val)==numel(this.cExperiment.dirs));
            this.stackSizes_val = val;
        end
        
        function val = get.timepointIndices(this)
            if isempty(this.timepointIndices_val)
                this.timepointIndices_val = cell(1,numel(this.cExperiment.dirs));
            end
            val = this.timepointIndices_val;
        end
        function set.timepointIndices(this,val)
            assert(numel(val)==numel(this.cExperiment.dirs));
            this.timepointIndices_val = val;
        end
        
        function loadTimelapse(this,pos,preloaded_cTimelapse)
            if nargin>2 && ~isempty(preloaded_cTimelapse)
                assert(isa(preloaded_cTimelapse,'babyTimelapse'),this.Errors.BadParams,...
                    'The "preloaded_cTimelapse" arg must be a "babyTimelapse" object');
                assert(pos>=1 && pos<=length(this.cExperiment.dirs),this.Errors.BadParams,...
                    'This cExperiment does not have position %u.',pos);
                this.cTimelapse = preloaded_cTimelapse;
                this.current_position = pos;
            elseif isempty(this.cTimelapse) || ...
                    (~isempty(this.cTimelapse) && this.current_position~=pos)
                this.cTimelapse = this.cExperiment.returnTimelapse(pos);
                this.current_position = pos;
            end
        end
        
        function addPos(this,pos)
            %ImageCache.addPos Helper function to add empty pos
            %   imcache.addPos(pos) loads the cTimelapse with index
            %   "pos" and extracts information to create storage arrays
            
            % Sanity checks:
            assert(pos>=1 && pos<=length(this.cExperiment.dirs),this.Errors.BadParams,...
                'This cExperiment does not have position %u.',pos);
            
            if isKey(this.posMap,pos) && ~isempty(this.stackSizes{pos})
                % Position has already been added
                return
            end
            
            this.loadTimelapse(pos);
            
            if any(this.cTimelapse.timepointsProcessed)
                tInd = find(this.cTimelapse.timepointsProcessed);
            else
                tInd = this.cTimelapse.timepointsToProcess;
            end
            if isempty(this.timepointIndices{pos})
                this.timepointIndices{pos} = tInd;
            end
            
            trapKeys = unique(arrayfun(@(x) length(x.trapLocations),...
                this.cTimelapse.cTimepoint(this.timepointIndices{pos})));
            assert(length(trapKeys)==1,this.Errors.BadTimelapse,...
                'Cannot have timepoints with different numbers of traps');
            trapKeys = 1:trapKeys;
            
            if this.cTimelapse.trapsPresent
                trapHeight = 2*this.cTimelapse.cTrapSize.bb_height+1;
                trapWidth = 2*this.cTimelapse.cTrapSize.bb_width+1;
            else
                % Traps not present
                trapHeight = this.cTimelapse.imSize(1);
                trapWidth = this.cTimelapse.imSize(2);
            end
            nTimepoints = length(this.timepointIndices{pos});
            sSize = [trapHeight,trapWidth,nTimepoints];
            if isempty(this.stackSizes{pos})
                this.stackSizes{pos} = sSize;
            end
            
            this.posNames{end+1} = this.cExperiment.dirs{pos};
            ipos = length(this.posNames);
            this.posMap(pos) = ipos;
            this.channels = unique([this.cTimelapse.channelNames(:);this.channels(:)]);
            this.trapMaps{ipos} = containers.Map(trapKeys,uint32(1:length(trapKeys)));
            
            % Information has been updated, so trigger InfoUpdated event:
            notify(this,'InfoUpdated');
        end
        
        function addPosChannel(this,pos,channel)
            %ImageCache.addPosChannel Helper function to add empty channel
            %   imcache.addPosChannel(pos,channel) loads the cTimelapse 
            %   with index "pos" and then initialises this.imcache with
            %   storage for that channel's images.
            
            assert(~isempty(this.cExperiment),this.Errors.NoExperiment,...
                'This ImageCache must be assigned a cExperiment before images can be loaded.');
            assert(isKey(this.posMap,pos),this.Errors.BadParams,...
                'The position needs to be added before a channel is added.');
            
            % Ensure the specified position is loaded:
            this.loadTimelapse(pos);
            
            if ischar(channel), channel = {channel}; end
            assert(iscellstr(channel), this.Errors.BadParams,...
                '"channel" must be specified as a char or cellstr');
            channel = reshape(channel,1,[]); % Ensure cellstr is vertical
            
            updated = false;
            
            for c = 1:length(channel)
                chname = channel{c};
                
                % Check that the timelapse has the specified channel
                channelNum = find(strcmp(this.cTimelapse.channelNames,chname),1);
                assert(~isempty(channelNum),this.Errors.NotAvailable,...
                    'The channel %s is not available at position %u',chname,pos);

                if ~isKey(this.channelMap,chname)
                    ichannel = this.channelMap.Count+1;
                    this.channelMap(chname) = ichannel;
                    updated = true;
                end
            end
            
            if updated
                % Information has been updated, so trigger InfoUpdated event:
                notify(this,'InfoUpdated');
            end
        end
        
        function loadTimepoint(this,pos,channel,timepoint)
            %ImageCache.loadTimepoint Load a timepoint into memory
            %   imcache.loadTimepoint(pos,channel,timepoint) loads the 
            %   cTimelapse with index "pos" and then reads the image for 
            %   "channel" at "timepoint" into memory.
            
            assert(~isempty(this.cExperiment),this.Errors.NoExperiment,...
                'This ImageCache must be assigned a cExperiment before images can be loaded.');
            assert(isKey(this.posMap,pos),this.Errors.BadParams,...
                'The position needs to be added before images can be loaded.');
            assert(1<=pos && pos<=size(this.imcache,1),this.Errors.BadParams,...
                'The position needs to be added before images can be loaded.');
            assert(isKey(this.channelMap,channel),this.Errors.BadParams,...
                'The channel needs to be added before images can be loaded.');
            
            ipos = this.posMap(pos);
            ichannel = this.channelMap(channel);
            assert(1<=ichannel && ichannel<=size(this.imcache,2),this.Errors.BadParams,...
                'The channel needs to be added before images can be loaded.');
            
            % Ensure the specified position is loaded:
            this.loadTimelapse(pos);
            
            % Ensure that the appropriate caches are initialised
            % Determine the correct channel number from the timelapse
            channelNum = find(strcmp(this.cTimelapse.channelNames,channel),1);
            assert(~isempty(channelNum),this.Errors.NotAvailable,...
                'The channel %s is not available at position %u',channel,pos);
            this.init_cache(ipos,ichannel);
            this.loadTimepointWithoutChecks(ipos,ichannel,channelNum,timepoint);
        end
        
        function loadPosChannel(this,pos,channel,timepoints,progress_bar)
            this.addPos(pos);
            this.addPosChannel(pos,channel);
            
            close_progress = false;
            if nargin<5 || ~isa(progress_bar,'Progress')
                % Initialise a progress bar
                progress_bar = Progress();
                % Centre the dialog box
                progress_bar.frame.setLocationRelativeTo([]);
                % Set the progress bar title
                progress_bar.frame.setTitle('Loading images...');
                close_progress = true;
            end
            
            ipos = this.posMap(pos);
            ichannel = this.channelMap(channel);
            
            % Ensure the specified position is loaded:
            this.loadTimelapse(pos);
            
            % Ensure that the cache is resized if needed
            % Determine the correct channel number from the timelapse
            channelNum = find(strcmp(this.cTimelapse.channelNames,channel),1);
            assert(~isempty(channelNum),this.Errors.NotAvailable,...
                'The channel %s is not available at position %u',channel,pos);
            this.init_cache(ipos,ichannel);
            
            assert(1<=ichannel && ichannel<=size(this.imcache,2),this.Errors.BadParams,...
                'The channel needs to be added before images can be loaded.');
            
            if nargin<4 || isempty(timepoints)
                nTimepoints = length(this.timepointIndices{pos});
                timepoints = 1:nTimepoints;
            end
            
            % Initialise bar to loop over each timepoint:
            progress_bar.push_bar('Reading time point',1,numel(timepoints));
            for t=1:numel(timepoints)
                timepoint = timepoints(t);
                trapsloaded = cellfun(@(trap) trap(timepoint),...
                    this.loadedtps{ipos,ichannel});
                if any(~trapsloaded)
                    this.loadTimepointWithoutChecks(ipos,ichannel,channelNum,timepoint);
                end
                progress_bar.set_val(t);
            end
            
            % Clean up progress bar:
            progress_bar.pop_bar; % finished all timepoints
            if close_progress
                progress_bar.frame.dispose;
            end
        end
        
        function createLocalCache(this,cachepath)
            if nargin<2 || isempty(cachepath)
                msg = 'Choose folder in which to store local cache...';
                fprintf([msg,'\n']);
                cachepath = uigetdir(fullfile(this.cExperiment.saveFolder,'..'),msg);
                if isequal(cachepath,0), return; end % user cancelled
            end
            
            assert(isdir(cachepath),this.Errors.BadParams,...
                'The directory %s does not exist',cachepath);
            
            cachefilename = fullfile(cachepath,'ImageCache.txt');
            if exist(cachefilename,'file')==2
                answer = questdlg(['A local cache already exists in that ',...
                    'directory. Do you want to overwrite it?'],...
                    'Overwrite existing cache?','Yes','No','No');
                if ~strcmp(answer,'Yes'), return; end % user cancelled
            end
            
            fh = fopen(cachefilename,'w');
            fprintf(fh,'%s\n',this.cExperiment.id);
            fclose(fh);
            this.local_cache_id = this.cExperiment.id;
            
            [basedir,reldir] = fileparts(cachepath);
            this.savedir_val = BABYutil.LinkedFile(reldir,basedir,...
                'Filter','/');
        end
        
        function openLocalCache(this,cachepath)
            if nargin<2 || isempty(cachepath)
                msg = 'Choose folder in which local cache is stored...';
                fprintf([msg,'\n']);
                cachepath = uigetdir(fullfile(this.cExperiment.saveFolder,'..'),msg);
                if isequal(cachepath,0), return; end % user cancelled
            end
            
            assert(isdir(cachepath),this.Errors.BadParams,...
                'The directory "%s" does not exist',cachepath);
            cachefilename = fullfile(cachepath,'ImageCache.txt');
            
            assert(isequal(exist(cachefilename,'file'),2),...
                this.Errors.BadParams,...
                '"%s" does not contain an ImageCache.txt file',cachepath);
            
            fh = fopen(cachefilename,'r');
            cacheid = fgetl(fh);
            fclose(fh);
            this.local_cache_id = cacheid;
            
            assert(strcmp(cacheid,this.cExperiment.id),this.Errors.BadParams,...
                'The local cache ID "%s" does not match that of the cExperiment "%s"',...
                cacheid,this.cExperiment.id);
            
            [basedir,reldir] = fileparts(cachepath);
            this.savedir_val = BABYutil.LinkedFile(reldir,basedir,...
                'Filter','/');
        end
        
        function saveLocalTrapStack(this,pos,trap,channel)
            assert(~isempty(this.savedir),this.Errors.NotSaved,...
                'The local cache has not yet been assigned');
            
            % Determine the correct pos index
            if ~isKey(this.posMap,pos), this.addPos(pos); end
            ipos = this.posMap(pos);
            
            % Ensure the specified position is loaded:
            this.loadTimelapse(pos);
            
            % Collate segmentation images if requested, otherwise just
            % gather the cached microscope images:
            imagemap = [];
            if ismember(channel,{'segoutline','segarea'})
                trapStack = zeros(this.stackSize,'uint8');
                for t=1:length(this.timepointIndices{pos})
                    tp = this.timepointIndices{pos}(t);
                    trapInfo = this.cTimelapse.cTimepoint(tp).trapInfo(trap);
                    if trapInfo.cellsPresent
                        tpim = trapStack(:,:,t);
                        for c=1:length(trapInfo.cellLabel)
                            segim = trapInfo.cell(c).segmented;
                            if strcmp(channel,'segarea')
                                segim = imfill(full(segim),'holes');
                            end
                            tpim(segim) = trapInfo.cellLabel(c);
                        end
                        trapStack(:,:,t) = tpim;
                    end
                end
                stackrange = [min(trapStack(:)), max(trapStack(:))];
                imagemap = [0,0,0;jet(max(1,double(stackrange(2))))];
                save_class = 'uint8';
            else
                % Retrieve the entire time series in one stack:
                trapStack = this.getTrapTimepoints(pos,trap,channel,...
                    1:length(this.timepointIndices{pos}),this.storage_class);
                
                % The previous call ensures that the following have been loaded:
                ichannel = this.channelMap(channel);
                itrap = this.trapMaps{ipos}(trap);
                
                stackrange = this.imrange{ipos,ichannel}{itrap};
                save_class = this.storage_class;
            end
            
            n_stacks = size(trapStack,3);
            n_per_row = ceil(sqrt(n_stacks));
            n_rows = floor(n_stacks/n_per_row);
            last_row = mod(n_stacks,n_per_row);
            imwidth = size(trapStack,1);
            imheight = size(trapStack,2);
            
            if last_row>0
                n_rows = n_rows+1;
            end
            tiled_image = zeros(n_per_row*imwidth,...
                n_rows*imheight,save_class);
            
            for i=1:n_rows
                if i==n_rows && last_row>0
                    ncols = last_row;
                else
                    ncols = n_per_row;
                end
                for j=1:ncols
                    xslice = (j-1)*imwidth+(1:imwidth);
                    yslice = (i-1)*imheight+(1:imheight);
                    tiled_image(xslice,yslice) = trapStack(:,:,(i-1)*n_per_row+j);
                end
            end
            
            % Create a sub-dir for this position if it doesn't exist:
            if ~isdir(fullfile(this.savedir,this.posNames{ipos}))
                mkdir(this.savedir,this.posNames{ipos});
            end
            filename = sprintf('trap_%03u_%s.png',trap,channel);
            filepath = fullfile(this.savedir,this.posNames{ipos},filename);
            description = sprintf(...
                ['n_tiles: %u; tile_array: %ux%u; image_dim: %ux%u; ',...
                'range: [%.1f,%.1f]; exptID: %s; posID: %u; ',...
                'trapID: %u; channel: %s'],...
                n_stacks, n_per_row, n_rows, imwidth, imheight,...
                stackrange(1), stackrange(2), this.cExperiment.id, ...
                pos, trap, channel);
            if isempty(imagemap)
                imwrite(tiled_image,filepath,'png','Description',description);
            else
                imwrite(tiled_image,imagemap,filepath,'png',...
                    'Transparency',[0,ones(1,size(imagemap,1)-1)],...
                    'Description',description);
            end
        end
        
        function loadLocalTrapStack(this,pos,trap,channel)
            assert(~isempty(this.savedir),this.Errors.NotSaved,...
                'The local cache has not yet been assigned');
            
            if ~isKey(this.posMap,pos), this.addPos(pos); end
            ipos = this.posMap(pos);
            % Ensure the specified position is loaded:
            this.loadTimelapse(pos);
            
            if ~isKey(this.channelMap,channel)
                this.addPosChannel(pos,channel);
            end
            ichannel = this.channelMap(channel);
            
            % Ensure that the appropriate caches are initialised
            this.init_cache(ipos,ichannel);
            
            % Load the timepoints for this trap into the cache
            itrap = this.trapMaps{ipos}(trap);
            
            % Create a sub-dir for this position if it doesn't exist:
            assert(isdir(fullfile(this.savedir,this.posNames{ipos})),...
                this.Errors.NotSaved,...
                'The position "%s" has not been saved to the cache',...
                this.posNames{ipos});
            
            filename = sprintf('trap_%03u_%s.png',trap,channel);
            filepath = fullfile(this.savedir,this.posNames{ipos},filename);
            if ~isequal(exist(filepath,'file'),2) && isa(this.cExperiment,'babyExperimentOmero')
                % Try the non-Omero channel format
                channel = regexprep(channel,'_([1-9])$','_00$1');
                filename = sprintf('trap_%03u_%s.png',trap,channel);
                filepath = fullfile(this.savedir,this.posNames{ipos},filename);
            end
            assert(isequal(exist(filepath,'file'),2),this.Errors.NotSaved,...
                'The channel "%s" for trap "%u" in position "%s" has not been saved to the cache',...
                channel,trap,this.posNames{ipos});
            
            iminf = imfinfo(filepath,'png');
            description = iminf.Description;
            parsed_desc = regexp(description,...
                ['n_tiles: (?<n_tiles>\d+); ',...
                'tile_array: (?<n_per_row>\d+)x(?<n_rows>\d+); ',...
                'image_dim: (?<imwidth>\d+)x(?<imheight>\d+); ',...
                'range: \[(?<stackmin>[-.0-9]+),(?<stackmax>[-.0-9]+)\]; ',...
                'exptID: (?<exptID>[^;]+); ',...
                'posID: (?<posID>\d+); trapID: (?<trapID>\d+); ',...
                'channel: (?<channel>[^;]+)$'],'names');
            assert(length(parsed_desc)==1,this.Errors.BadImage,...
                'Encountered badly formatted image: pos %u, trap %u, channel: %s',...
                pos, trap, channel);
            
            assert(strcmp(parsed_desc.exptID,this.cExperiment.id),this.Errors.BadImage,...
                'cExperiment ID does not match image: pos %u, trap %u, channel: %s',...
                pos, trap, channel);
            assert(sscanf(parsed_desc.posID,'%u')==pos,this.Errors.BadImage,...
                'Position number does not match image''s (pos %u, trap %u, channel: %s)',...
                pos, trap, channel);
            assert(sscanf(parsed_desc.trapID,'%u')==trap,this.Errors.BadImage,...
                'Trap number does not match image''s (pos %u, trap %u, channel: %s)',...
                pos, trap, channel);
            assert(strcmp(parsed_desc.channel,channel),this.Errors.BadImage,...
                'Channel name does not match image''s (pos %u, trap %u, channel: %s)',...
                pos, trap, channel);
            
            n_stacks = sscanf(parsed_desc.n_tiles,'%u');
            n_per_row = sscanf(parsed_desc.n_per_row,'%u');
            n_rows = sscanf(parsed_desc.n_rows,'%u');
            imwidth = sscanf(parsed_desc.imwidth,'%u');
            imheight = sscanf(parsed_desc.imheight,'%u');
            stackrange = [sscanf(parsed_desc.stackmin,'%f'),...
                sscanf(parsed_desc.stackmax,'%f')];
            
            tiled_image = imread(filepath,'png');
            last_row = mod(n_stacks,n_per_row);

            % Recreate trap stack
            trapStack = zeros(imwidth,imheight,n_stacks,'uint8');
            for i=1:n_rows
                if i==n_rows && last_row>0
                    ncols = last_row;
                else
                    ncols = n_per_row;
                end
                for j=1:ncols
                    xslice = (j-1)*imwidth+(1:imwidth);
                    yslice = (i-1)*imheight+(1:imheight);
                    trapStack(:,:,(i-1)*n_per_row+j) = tiled_image(xslice,yslice);
                end
            end
            
            assert(all(size(trapStack)==size(this.imcache{ipos,ichannel}{itrap})),...
                this.Errors.BadImage,['Trap stack dimensions do not match cExperiment''s ',...
                '(pos %u, trap %u, channel: %s)'], pos, trap, channel);
            
            % Overwrite the existing cache
            this.imcache{ipos,ichannel}{itrap} = trapStack;
            this.imrange{ipos,ichannel}{itrap} = stackrange;
            this.loadedtps{ipos,ichannel}{itrap}(:) = true;
        end
        
        function saveLocalPosChannel(this,pos,channel,progress_bar)
            if ismember(channel,{'segoutline','segarea'})
                % Determine the correct pos index
                if ~isKey(this.posMap,pos), this.addPos(pos); end
                ipos = this.posMap(pos);
                % Ensure the specified position is loaded:
                this.loadTimelapse(pos);
            else 
                if ~isKey(this.posMap,pos) || ~isKey(this.channelMap,channel)
                    this.loadPosChannel(pos,channel,[],progress_bar);
                end
                ipos = this.posMap(pos);
                ichannel = this.channelMap(channel);

                if ichannel<1 || ichannel>size(this.loadedtps,2) || ...
                        isempty(this.loadedtps{ipos,ichannel})
                    this.loadPosChannel(pos,channel,[],progress_bar);
                end

                % Load the timepoints into the cache if necessary
                areloaded = cellfun(@all,this.loadedtps{ipos,ichannel});
                if ~all(areloaded)
                    this.loadPosChannel(pos,channel,[],progress_bar);
                end
            end

            close_progress = false;
            if nargin<4 || ~isa(progress_bar,'Progress')
                % Initialise a progress bar
                progress_bar = Progress();
                % Centre the dialog box
                progress_bar.frame.setLocationRelativeTo([]);
                % Set the progress bar title
                progress_bar.frame.setTitle('Saving images...');
                close_progress = true;
            end
            
            trapKeys = this.trapMaps{ipos}.keys;
            
            % Save images for all traps
            % Start a progress bar for saving traps
            progress_bar.push_bar('Writing trap',1,length(trapKeys));
            % For each trap, save a cached image
            for trapindex = 1:length(trapKeys)
                trap = trapKeys{trapindex};
                this.saveLocalTrapStack(pos,trap,channel)
                progress_bar.set_val(trapindex);
            end
            
            % Clean up progress bar:
            progress_bar.pop_bar; % finished all traps
            if close_progress
                progress_bar.frame.dispose;
            end
        end
        
        function saveLocal(this,poses,channels,progress_bar)
            close_progress = false;
            
            if nargin<2 || isempty(poses)
                poses = 1:length(this.cExperiment.dirs);
            end
            
            if nargin<3 || isempty(channels)
                channels = this.cExperiment.channelNames;
            end
            
            if nargin<4 || ~isa(progress_bar,'Progress')
                % Initialise a progress bar
                progress_bar = Progress();
                % Centre the dialog box
                progress_bar.frame.setLocationRelativeTo([]);
                % Set the progress bar title
                progress_bar.frame.setTitle('Saving image cache...');
                close_progress = true;
            end
            
            % Start a progress bar for each position
            progress_bar.push_bar('Position',1,length(poses));
            for p=1:length(poses)
                progress_bar.set_val(p);
                pos = poses(p);
                progress_bar.push_bar('Channel',1,length(channels));
                for c=1:length(channels)
                    progress_bar.set_val(c);
                    channel = channels{c};
                    this.saveLocalPosChannel(pos,channel,progress_bar);
                end
                this.clearcache;
                progress_bar.pop_bar; % finished all channels
            end
            
            % Clean up progress bar:
            progress_bar.pop_bar; % finished all positions
            if close_progress
                progress_bar.frame.dispose;
            end
        end
        
        function val = get.savedir(this)
            val = '';
            if ~isempty(this.savedir_val) && ...
                    isa(this.savedir_val, 'BABYutil.LinkedFile') && ...
                    ~strcmp(this.local_cache_id,'USER_CANCELLED')
                try
                    val = this.savedir_val.filename;
                catch err
                    if strcmp(err.identifier,BABYutil.LinkedFile.errUserCancel)
                        this.local_cache_id = 'USER_CANCELLED';
                    else
                        rethrow(err)
                    end
                end
                if isempty(this.local_cache_id)
                    this.openLocalCache(val);
                end
            end
        end
        
        function val = get.cExperiment(this)
            val = this.cExperiment_val;
            assert(isa(val,'babyExperiment'),this.Errors.NoExperiment,...
                'This ImageCache must be assigned a valid cExperiment.');
        end
        
        function set.cExperiment(this,val)
            if ~isa(val,'babyExperiment') || ~isscalar(val)
                error(this.Errors.WrongType,...
                    'The cExperiment property must be set to a single "babyExperiment" object');
            end
            this.cExperiment_val = val;
        end
        
        function val = get.storage_class(this)
            val = this.storage_class_val;
        end
        
        function set.storage_class(this,val)
            if ~ismember(val,{'uint8','uint16','double'})
                error('"storage_class" must be one of "uint8", "uint16" or "double"');
            end
            this.storage_class_val = val;
        end
        
        function val = get.storage_model(this)
            if isempty(this.storage_model_val)
                this.storage_model_val = 'bytrap';
            end
            val = this.storage_model_val;
        end
        
        function set.storage_model(this,val)
            if ~isempty(this.storage_model_val)
                error('"storage_model" has already been set for this cache and cannot be changed');
            end
            if ~ismember(val,{'bytrap','bytimepoint'})
                error('"storage_model" must be one of "bytrap" or "bytimepoint"');
            end
            this.storage_class_val = val;
        end
        
        function val = get.posMap(this)
            if isempty(this.posMap_val)
                this.posMap_val = containers.Map('KeyType','double','ValueType','uint32');
            end
            val = this.posMap_val;
        end
        function set.posMap(this,val)
            assert(isa(val,'containers.Map') ...
                && isequal(val.KeyType,'double') ...
                && isequal(val.ValueType,'uint32'),...
                'posMap must be a Map with KeyType "double" and ValueType "uint32"');
            this.posMap_val = val;
        end
        
        function val = get.channelMap(this)
            if isempty(this.channelMap_val)
                this.channelMap_val = containers.Map('KeyType','char','ValueType','uint32');
            end
            val = this.channelMap_val;
        end
        function set.channelMap(this,val)
            assert(isa(val,'containers.Map') ...
                && isequal(val.KeyType,'char') ...
                && isequal(val.ValueType,'uint32'),...
                'channelMap must be a Map with KeyType "char" and ValueType "uint32"');
            this.channelMap_val = val;
        end
    end
    
    methods (Access=private)
        function init_cache(this,ipos,ichannel)
            % This function assumes that addPos and addPosChannel have both
            % been run to obtain an ipos and ichannel
            
            % Initialise the image cache if it hasn't been:
            if isempty(this.imcache)
                npos = this.posMap.Count;
                nchannel = this.channelMap.Count;
                this.imcache = cell(npos,nchannel);
                this.imrange = cell(npos,nchannel);
                this.loadedtps = cell(npos,nchannel);
            else
                cachesize = size(this.imcache);
                if cachesize(1)<this.posMap.Count
                    npos = this.posMap.Count;
                    extracell = cell(npos-cachesize(1),cachesize(2));
                    this.imcache(end+1:npos,:) = extracell;
                    this.imrange(end+1:npos,:) = extracell;
                    this.loadedtps(end+1:npos,:) = extracell;
                end
                if cachesize(2)<this.channelMap.Count
                    nchannel = this.channelMap.Count;
                    extracell = cell(cachesize(1),nchannel-cachesize(2));
                    this.imcache(:,end+1:nchannel) = extracell;
                    this.imrange(:,end+1:nchannel) = extracell;
                    this.loadedtps(:,end+1:nchannel) = extracell;
                end
            end
            
            % Pre-allocate image arrays for each trap:
            nTraps = this.trapMaps{ipos}.Count;
            cachelen = length(this.imcache{ipos,ichannel});
            if isempty(this.imcache{ipos,ichannel})
                this.imcache{ipos,ichannel} = cell(nTraps,1);
                this.imrange{ipos,ichannel} = cell(nTraps,1);
                this.loadedtps{ipos,ichannel} = cell(nTraps,1);
            elseif cachelen<nTraps
                extracell = cell(nTraps-cachelen,1);
                this.imcache{ipos,ichannel}(end+1:nTraps) = extracell;
                this.imrange{ipos,ichannel}(end+1:nTraps) = extracell;
                this.loadedtps{ipos,ichannel}(end+1:nTraps) = extracell;
            end
            for itrap = 1:nTraps
                if isempty(this.imcache{ipos,ichannel}{itrap}) || ...
                        ~isequal(size(this.imcache{ipos,ichannel}{itrap}),this.stackSize)
                    this.imcache{ipos,ichannel}{itrap} = ...
                        zeros(this.stackSize,this.storage_class);
                    this.imrange{ipos,ichannel}{itrap} = [Inf,-Inf];
                    this.loadedtps{ipos,ichannel}{itrap} = ...
                        false(this.stackSize(3),1);
                end
            end
        end
        
        function im = imnorm(this,im,imin,imax,sclass)
            if nargin<5, sclass = this.storage_class; end
            im = (double(im)-imin)/(imax-imin);
            switch sclass
                case 'uint8'
                    im = uint8(double(intmax(sclass))*im);
                case 'uint16'
                    im = uint16(double(intmax(sclass))*im);
            end
        end
        
        function im = imraw(this,im,imin,imax)
            if ~strcmp(this.storage_class,'double')
                im = double(im)/double(intmax(this.storage_class));
            end
            im = imin + im*(imax-imin);
        end
        
        function loadTimepointWithoutChecks(this,ipos,ichannel,channelNum,timepoint)
            %ImageCache.loadTimepointWithoutChecks Quick load a timepoint
            %	Use this function with caution: it assumes that the correct
            %	cTimelapse has already been loaded and requires the
            %	channelNum index for this cTimelapse to be predetermined.
            %	It is intended to be used by internal functions that 
            %   preload these values.
            trapStack = ...
                this.cTimelapse.returnTrapsTimepoint([],...
                this.timepointIndices{this.current_position}(timepoint),...
                channelNum,'max');
            % Loop over the required traps to set the array:
            traps = this.trapMaps{ipos}.keys;
            for trap = [traps{:}]
                itrap = this.trapMaps{ipos}(trap);
                trapim = double(trapStack(:,:,trap));
                % Rescale if required:
                oldrange = this.imrange{ipos,ichannel}{itrap};
                newmin = min(oldrange(1),min(trapim(:)));
                newmax = max(oldrange(2),max(trapim(:)));
                if oldrange(1)>newmin || oldrange(2)<newmax
                    areloaded = this.loadedtps{ipos,ichannel}{itrap};
                    if any(areloaded)
                        previms = this.imraw(this.imcache{ipos,ichannel}{itrap}(:,:,areloaded),...
                            oldrange(1),oldrange(2));
                        this.imcache{ipos,ichannel}{itrap}(:,:,areloaded) = ...
                            this.imnorm(previms,newmin,newmax);
                    end
                    this.imrange{ipos,ichannel}{itrap} = [newmin,newmax];
                end
                this.imcache{ipos,ichannel}{itrap}(:,:,timepoint) = ...
                    this.imnorm(trapim,newmin,newmax);
                this.loadedtps{ipos,ichannel}{itrap}(timepoint) = true;
            end
        end
        
    end
end
