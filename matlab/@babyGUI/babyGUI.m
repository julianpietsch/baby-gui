classdef babyGUI < handle
    %babyGUI Edit segmentation images in a cTimelapse
    %   The babyGUI helps the user to manually curate cell outlines 
    %   and cell tracking
    
    properties
        cExperiment % A cExperiment
        cTimelapse % A cTimelapse
        imcache % An ImageCache object
        
        % Store basic window elements in the frame struct
        frame = struct(...
            'fig',matlab.ui.Figure.empty,... % Handle to the figure
            'imAx',matlab.graphics.axis.Axes.empty,... % Handle to trap image display axis
            'im',matlab.graphics.primitive.Image.empty,... % Handle to Image object
            'slider',matlab.ui.control.UIControl.empty,... % Handle to the timepoint slider control
            'xslider',matlab.ui.control.UIControl.empty,... % Handle to the X slider control
            'yslider',matlab.ui.control.UIControl.empty,... % Handle to the Y slider control
            'track',matlab.graphics.axis.Axes.empty,... % Handle to cell tracking axis
            'pos',matlab.graphics.axis.Axes.empty,... % Handle to pos overview display axis
            'ctrl',matlab.ui.container.Panel.empty); % Handle to control panel container
        
        % Control handles
        ctrls = struct() % structure containing handles to general controls
        trackMarkerLine
        
        % Control whether or not to autosave
        autosave = true
        haschanged = false;
        expthaschanged = false;
    end
    
    properties (Dependent)
        imzoom
        zoomactive
        zoomoffset
        zoomcentre
        currentPos
        currentTrap
        currentTimepoint
        currentChannels
        currentCell
        cellLabels
        cellTracks
        cellLocs
        times
        title % A character array specifying the figure title
        Enable
        layout % either 'vertical' (default), 'horizontal', or 'overview'
        posImFromCache % a boolean specifying whether to generate pos im from ImageCache
        dofocusannot % a boolean specifying whether focus annotation mode is active or not
        srcOffset % [x,y] offset from source image for the displayed tiles
    end
    
    properties (SetAccess=private)
        currentDaughterInd
        births
    end
    
    properties (Access=private)
        % Position properties
        tpInds
        ntimepoints
        ntraps
        timelabel
        posImChannel
        availableChannels
        
        % Image display properties
        windowUpdateRequired = false
        tilesContext % Number of time points to show either side of this one
        srcWidth % Width of a trap/tile in the source image
        srcHeight % Height of a trap/tile in the source image
        tileWidth % Width of a single tile in the displayed tiled image
        tileHeight % Height of a single tile in the dispalyed tiled image
        alphascale = 0.8
        currentColours
        currentPosImColours
        trackColours
        channelNormalisation = struct()
        
        % Intermediate image caches
        bgdIm
        trapIms
        posIm
        posSegIm
        
        progressBar % Instance of Progress class
        
        % State of dependent variables
        currentPos_val
        currentTrap_val
        currentTimepoint_val
        currentChannels_val
        cellLabels_val
        cellTracks_val
        cellLocs_val
        imzoom_val = 1.5
        zoomactive_val
        zoomcentre_val
        times_val
        Enable_val = 'off'
        
        % Variables for draggable contours
        dragPoints = {}
        dragCentre = {}
        currentOutline = {}
        dragUpdate = false
        dragAsRect = false
        cpoint = []
        points = []
        centre = []
        radii = []
        angles = []
        dragtp = 0
        tp_pad = 0
        yoffs = []
        px = []
        py = []
        
        % Save stack variables
        saveStackDir
        
        % Lineage variables
        cachedAreas
        arrowHandles = {}
    end
    
    properties (Constant)
        % Define a constant structure to set default dimensions (in pixels)
        % for window elements:
        dimensions = struct(...
            'lh',15,... % Label height
            'sh',20,... % Selection menu height
            'ch',30,... % Control height
            'sp',3,... % Spacing between controls
            'impad',3,... % Number of pixels worth of padding between images
            'xaxis',30,... % Space to allow for x axes
            'yaxis',30,... % Space to allow for y axes
            'minWidth',400,... % Minimum width of scrolling panel
            'minHeight',300,... % Minimum height of scrolling panel
            'ctrlHeight',370,... % Height of the control panel
            'ctrlWidth',500,... % Width of the control panel
            'posHeight',370,... % Height of the overview panel
            'posWidth',370,... % Width of the overview panel
            'trackHeight',50,... % Height of the cell track plot axis
            'slideHeight',20,... % Height of the timepoint slider control
            'tabWidth',65,... % Width of tab buttons
            'screenHeightPad',70); % Padding to account for OS menus
        
        zoomlevels = [0.5,0.6,0.7,1,1.5,2,2.5,3,4,5];
        
        zoomactivescale = 2;
        
        % Define the edit modes and shortcut keys
        editModeMap = struct(...
            'setactive',{'<no action> / set active (s)','s'},...
            'editoutline',{'Edit outline / set active (e)','e'},...
            'addoutline',{'Add outline / set active (a)','a'},...
            'addcell',{'Add new cell / set active (A)','A'},...
            'removeoutline',{'Remove outline / set active (x)','x'},...
            'removecell',{'Remove cell / set active (X)','X'},...
            'swaptracks',{'Join tracks / set active (q)','q'},...
            'breaktracks',{'Break tracks / set active (d)','d'},...
            'selectcells',{'Select / deselect cells (z)','z'},...
            'setdaughter',{'Set daughter / set active (g)','g'},...
            'unsetdaughter',{'Unset daughter / set active (G)','G'});
        
        % Define the outline editing modes and shortcut keys
        outlineModeMap = struct(...
            'clickmatch',{'Morph to click (1)','1'},...
            'dragpoints',{'Drag contour (2)','2'},...
            'dragrect',{'Drag rectangle (3)','3'});
        
        % Define the outline editing modes and shortcut keys
        addModeMap = struct(...
            'autogen',{'Best guess (b)','b'},...
            'asbud',{'Add as bud (B)','B'},...
            'asrod',{'Add as rod (r)','r'});
        
        % Define map from layout id to label
        layoutMap = struct(...
            'vertical','Vertical layout',...
            'horizontal','Horizontal layout',...
            'overview','Large overview');
        
        trackDisplayMap = struct(...
            'tracks','Standard',...
            'mdarea','Daughter areas',...
            'yloc','Vertical location');
        
        trackSizeMap = struct(...
            'small',{'Small',50},...
            'medium',{'Medium',100},...
            'large',{'Large',200});
        
        nknot_values = [4,6,8,10,12,16];
        
        % The maximum number of colours used to distinguish cells:
        MaxTrackColours = 100;
        
        boxScale = 1.5;
        cpointScale = 3.5;
    end
    
    methods
        function this = babyGUI(cExperiment,varargin)
            ip = inputParser;
            ip.addRequired('cExperiment',@(x) isa(x,'babyExperiment') && isscalar(x));
            ip.addOptional('pos',1,@(x) (isnumeric(x) && isscalar(x)) || isempty(x));
            ip.addParameter('ImageCache',[],@(x) isempty(x) || isa(x,'ImageCache'));
            ip.addParameter('TilesContext',2,@(x) isnumeric(x) && isscalar(x) && round(x)==x && x>=0);
            ip.addParameter('Channels',{},@(x) iscellstr(x) || ...
                (ischar(x) && size(x,1)==1));
            ip.addParameter('CachedOverview',[],@(x) isempty(x) || ...
                (isscalar(x) && islogical(x)));
            ip.parse(cExperiment,varargin{:});
            
            % Save handles to cExperiment and ImageCache
            this.cExperiment = cExperiment;
            this.imcache = ip.Results.ImageCache;
            if isempty(this.imcache)
                if isempty(cExperiment.imcache)
                    % Create a new image cache and assign to cExperiment
                    this.imcache = ImageCache(cExperiment);
                    this.cExperiment.imcache = this.imcache;
                else
                    % Use imcache saved with cExperiment
                    this.imcache = cExperiment.imcache;
                end
            end
            
            % Initialise the image cache
            this.currentPos = ip.Results.pos;
            initChannels = ip.Results.Channels;
            if isempty(initChannels)
                chnames = this.cTimelapse.channelNames;
                isTL = contains(chnames,{'field','dic','tl','phase'},'IgnoreCase',true);
                tlch = chnames(isTL);
                if isempty(tlch)
                    initChannels = chnames{1};
                else
                    initChannels = tlch{1};
                end
            end
            this.currentChannels = initChannels;
            this.posImChannel = find(strcmp(this.cTimelapse.channelNames,this.currentChannels{1}),1);
            this.imcache.addPos(this.currentPos);
            for c=1:length(this.currentChannels)
                this.imcache.addPosChannel(this.currentPos,this.currentChannels{c});
            end
            
            % Initialise properties
            this.srcWidth = this.imcache.stackSize(2);
            this.srcHeight = this.imcache.stackSize(1);
            this.tilesContext = ip.Results.TilesContext;
            this.currentTrap = 1;
            this.currentTimepoint = 1;
            
            this.progressBar = Progress();
            this.progressBar.frame.setLocationRelativeTo([]);
            
            cachedOverview = ip.Results.CachedOverview;
            if isempty(cachedOverview)
                % Attempt to detect whether we can load position images in
                % order to determine default for posImFromCache:
                this.posImFromCache = ...
                    ~isempty(this.imcache.savedir) && ...
                    ~isa(this.cExperiment,'babyExperimentOmero') && ...
                    ~isempty(this.cTimelapse.cTimepoint) && ...
                    ~isempty(this.cTimelapse.cTimepoint(1).filename) && ...
                    exist(fullfile(this.cTimelapse.timelapseDir,...
                    this.cTimelapse.cTimepoint(1).filename{1}),'file')~=2;
            else
                this.posImFromCache = cachedOverview;
            end

            % The first call to UpdateWindowElements also creates the
            % figure and we request auto-adaptation for screen size:
            this.UpdateWindowElements(true);
            
            % Set up the figure title
            if isa(cExperiment,'babyExperimentOmero') && ...
                    ~isnumeric(cExperiment.omeroDs)
                title_info = {char(cExperiment.omeroDs.getName.getValue),...
                    char(cExperiment.OmeroDatabase.getDate(cExperiment.omeroDs))};
            else
                title_info = {cExperiment.saveFolder};
                if ~isempty(cExperiment.metadata) && isfield(cExperiment.metadata,'date')
                    title_info = [title_info,{cExperiment.metadata.date}];
                end
            end
            this.title = strjoin(title_info,' - ');
            
            % Set up close call back
            set(this.frame.fig, 'CloseRequestFcn', @(~,~) delete(this));
            
            % Set up controls in the control panel:
            this.InitialiseControls;
            
            this.refreshColours;
            this.refreshSegIms;
            this.refreshPosImage;
            this.refreshTrackPlot;
            this.Enable = 'on';
            drawnow;
            
            % Turn off a warning that gets produced by splinefit
            warning('off','MATLAB:nargchk:deprecated');
        end
        
        function delete(this)
            this.saveTimelapse;
            this.saveExperiment;
            if ~isempty(this.frame.fig) && isvalid(this.frame.fig)
                delete(this.frame.fig);
            end
            if ~isempty(this.progressBar) && isvalid(this.progressBar)
                this.progressBar.frame.dispose;
            end
        end
        
        function saveTimelapse(this)
            % Automatically called when changing positions, closing the 
            % window, or deleting the handle
            if isempty(this.cTimelapse) || ~this.autosave || ~this.haschanged
                return
            end
            enable = this.Enable;
            if ~strcmp(enable,'off')
                this.Enable = 'off';
                drawnow;
            end
            this.cExperiment.cTimelapse=this.cTimelapse;
            this.cExperiment.saveTimelapseExperiment(this.currentPos_val,false);
            if ~strcmp(enable,'off')
                this.Enable = 'on';
            end
            this.haschanged = false;
        end
        
        function saveExperiment(this)
            % Automatically called when changing positions, closing the 
            % window, or deleting the handle
            if isempty(this.cExperiment) || ~this.autosave || ~this.expthaschanged
                return
            end
            enable = this.Enable;
            if ~strcmp(enable,'off')
                this.Enable = 'off';
                drawnow;
            end
            this.cExperiment.saveExperiment;
            if ~strcmp(enable,'off')
                this.Enable = 'on';
            end
            this.expthaschanged = false;
        end
        
        function refreshColours(this,~,~)
            ncells = size(this.cellTracks,1);
            if ncells==0
                return
            end
            % Colours depend on the editing mode
            modenames = fieldnames(this.editModeMap);
            switch modenames{this.ctrls.editmode.Value}
                case {'setactive','removeoutline','breaktracks'}
                    % Colour by track
                    this.currentColours = this.trackColours(...
                        mod(this.cellLabels{this.currentTrap}-1,this.MaxTrackColours)+1,:);
                case {'editoutline','addoutline','addcell','removecell'}
                    othercolour = [0.1,0.1,0.9]; % all cells are blue...
                    this.currentColours = othercolour(ones(ncells,1),:);
                    this.currentColours(this.currentCell,:) = ...
                        [0.7,0.8,0]; % ...except the active cell
                case 'swaptracks'
                    % 'Swap tracks (l)'
                    othercolour = [0.1,0.1,0.9]; % all cells are blue...
                    this.currentColours = othercolour(ones(ncells,1),:);
                    this.currentColours(this.currentCell,:) = ...
                        [0.7,0.8,0]; % ...except the selected cell
                case 'selectcells'
                    othercolour = [1,0,0]; % all cells are red...
                    this.currentColours = othercolour(ones(ncells,1),:);
                    selected = logical(full(this.cTimelapse.cellsToPlot(...
                        this.currentTrap,this.cellLabels{this.currentTrap})));
                    selcolour = [0,1,0]; % ...except the selected cells
                    this.currentColours(selected,:) = ...
                        selcolour(ones(sum(selected),1),:);
                case {'setdaughter','unsetdaughter'}
                    othercolour = [0,0,1]; % all cells are blue...
                    this.currentColours = othercolour(ones(ncells,1),:);
                    this.currentColours(this.currentCell,:) = ...
                        [1,1,0]; % ...except the active cell
                    % ...and its daughters
                    cLabels = this.cellLabels{this.currentTrap};
                    cellLabel = cLabels(this.currentCell);
                    daughterLabels = find(this.cTimelapse.cellMothers(this.currentTrap,:)==cellLabel);
                    % Filter to current labels
                    daughterLabels = daughterLabels(ismember(daughterLabels,cLabels));
                    daughters = arrayfun(@(x) find(cLabels==x,1),daughterLabels);
                    if ~isempty(daughters)
                        daughtercolour = [0.8,0.5,0];
                        this.currentColours(daughters,:) = ...
                            daughtercolour(ones(length(daughters),1),:);
                    end
            end
        end
        
        function refreshTrackColours(this)
            this.trackColours = jet(this.MaxTrackColours);
            this.trackColours = this.trackColours(randperm(this.MaxTrackColours),:);
        end
        
        function colours = getTrackColours(this)
            colours = this.trackColours;
        end
        
        function refreshCellLabels(this,traps)
            if nargin<2 || isempty(traps)
                traps = 1:this.ntraps;
            end
            if isempty(this.cellLabels_val)
                this.cellLabels_val = cell(1,this.ntraps);
            end
            for t=traps(:)'
                cellsPresent = arrayfun(@(x) logical(x.trapInfo(t).cellsPresent),...
                    this.cTimelapse.cTimepoint(this.tpInds));
                tps = this.tpInds(cellsPresent);
                cLbls = arrayfun(@(x) double(x.trapInfo(t).cellLabel),...
                    this.cTimelapse.cTimepoint(tps),'Uni',0);
                cLbls = unique(horzcat(cLbls{:}));
                this.cellLabels_val{t} = cLbls(cLbls~=0);
            end
        end
        
        function cellArea = getCellArea(this,trap,tp,lbl)
            if isempty(this.cachedAreas)
                this.cachedAreas = cell(this.ntraps,numel(this.tpInds));
            end
            trapInfo = this.cTimelapse.cTimepoint(this.tpInds(tp)).trapInfo(trap);
            ncells = numel(trapInfo.cellLabel);
            if isempty(this.cachedAreas{trap,tp})
                this.cachedAreas{trap,tp} = NaN(ncells,1);
            end
            if numel(this.cachedAreas{trap,tp}) < ncells
                this.cachedAreas{trap,tp}(end+1:ncells) = NaN;
            end
            cInd = find(trapInfo.cellLabel==this.cellLabels{trap}(lbl),1);
            if isempty(cInd)
                cellArea = 0;
                return
            end
            if isnan(this.cachedAreas{trap,tp}(cInd))
                segim = trapInfo.cell(cInd).segmented;
                cmask = imfill(full(segim),'holes');
                this.cachedAreas{trap,tp}(cInd) = sum(cmask(:));
            end
            cellArea = this.cachedAreas{trap,tp}(cInd);
        end
        
        function clearCellAreas(this,trap,tp)
            if nargin<2
                this.cachedAreas = [];
            end
            if ~isempty(this.cachedAreas)
                if nargin<3
                    [this.cachedAreas{trap,:}] = deal([]);
                else
                    this.cachedAreas{trap,tp} = [];
                end
            end
        end
        
        function updateTrapList(this)
            trapText = arrayfun(@(x) sprintf('Trap %u',x),1:this.ntraps,'Uni',0);
            hasCells = ~cellfun(@isempty,this.cellLabels);
            cellsToPlot = full(this.cTimelapse.cellsToPlot);
            hasSelected = arrayfun(@(t) any(cellsToPlot(t,:)),1:length(this.cellLabels));
            trapText(hasCells) = strcat(trapText(hasCells),' [*]');
            trapText(hasSelected) = strcat(trapText(hasSelected),' [+]');
            forCura = this.forCuration;
            forCura = ismember(1:length(this.cellLabels),forCura(:,1));
            trapText(forCura) = strcat(trapText(forCura),' [.]');
            isCura = this.isCurated;
            isCura = ismember(1:length(this.cellLabels),isCura(:,1));
            trapText(isCura) = strcat(trapText(isCura),' [C]');
            val = min(max([1,this.ctrls.traps.Value]),this.ntraps);
            set(this.ctrls.traps,'String',trapText,'Value',val);
        end
        
        %% Get/set functions
        function val = get.currentPos(this)
            val = this.currentPos_val;
        end
        function set.currentPos(this,val)
            if ~(isnumeric(val) && isscalar(val)) && ~isempty(val)
                error('pos must be a scalar numeric');
            end
            load_timelapse = true;
            if isempty(val)
                load_timelapse = false;
                val = this.cExperiment.currentPos;
            end
            if ~ismember(val,1:length(this.cExperiment.dirs))
                error('pos out of range');
            end
            if this.currentPos_val==val, return; end
            % If there is a loaded cTimelapse, make sure we save it before
            % changing positions:
            this.saveTimelapse;
            % Reset normalisation
            this.channelNormalisation = struct();
            if load_timelapse
                % Load timelapse and setup properties
                this.cTimelapse = this.cExperiment.loadCurrentTimelapse(val);
            else
                this.cTimelapse = this.cExperiment.cTimelapse;
            end
            this.imcache.loadTimelapse(val,this.cTimelapse);
            this.imcache.addPos(val);
            if ~isequal([this.srcHeight,this.srcWidth],this.imcache.stackSize(1:2))
                this.srcWidth = this.imcache.stackSize(2);
                this.srcHeight = this.imcache.stackSize(1);
                this.windowUpdateRequired = true;
            end
            this.availableChannels = this.cTimelapse.channelNames(:)';
            if isfield(this.ctrls,'saveChannel')
                svch = this.ctrls.saveChannel;
                if ishandle(svch) && isvalid(svch)
                    newvalue = find(strcmp(this.cTimelapse.channelNames,...
                        svch.String{svch.Value}),1);
                    if isempty(newvalue), newvalue = 1; end
                    set(this.ctrls.saveChannel,'Value',newvalue,...
                        'String',this.cTimelapse.channelNames);
                end
            end
            if any(this.cTimelapse.timepointsProcessed)
                this.tpInds = find(this.cTimelapse.timepointsProcessed);
            else
                this.tpInds = this.cTimelapse.timepointsToProcess;
            end
            this.ntimepoints = length(this.tpInds);
            if this.currentTimepoint > this.ntimepoints
                this.currentTimepoint = this.ntimepoints;
            end
            this.ntraps = length(this.cTimelapse.cTimepoint(this.tpInds(1)).trapLocations);
            
            % Clear the area cache
            this.clearCellAreas();
            
            % Check that cellsToPlot has correct dimensions
            nrows = size(this.cTimelapse.cellsToPlot,1);
            if nrows<this.ntraps
                this.cTimelapse.cellsToPlot(nrows+1:this.ntraps,:) = 0;
            end
            
            % Check that cellMothers has correct dimensions
            nrows = size(this.cTimelapse.cellMothers,1);
            if nrows<this.ntraps
                this.cTimelapse.cellMothers(nrows+1:this.ntraps,:) = 0;
            end
            
            this.cellLabels_val = [];
            this.refreshCellLabels;
            if ~isempty(this.cExperiment.metadata) && ...
                    isfield(this.cExperiment.metadata,'logTimes') && ...
                    all(size(this.cExperiment.metadata.logTimes) ...
                    == [length(this.cExperiment.dirs),length(this.cTimelapse.timepointsProcessed)])
                this.times_val = this.cExperiment.metadata.logTimes(val,this.tpInds);
                if max(this.times_val)>600
                    this.times_val = this.times_val / 60;
                    this.timelabel = 'Time (hours)';
                else
                    this.timelabel = 'Time (min)';
                end
            else
                this.times_val = 1:this.ntimepoints;
                this.timelabel = 'Time point';
            end
            this.currentPos_val = val;
        end
        
        function val = get.currentTrap(this)
            val = this.currentTrap_val;
        end
        function set.currentTrap(this,val)
            if ~(isnumeric(val) && isscalar(val))
                error('trap must be a scalar numeric');
            end
            if ~ismember(val,1:length(this.cTimelapse.cTimepoint(this.tpInds(1)).trapLocations))
                error('trap out of range');
            end
            %if this.currentTrap_val==val, return; end
            this.currentTrap_val = val;
            this.refreshTracks;
            this.refreshTrackColours;
        end
        
        function val = get.cellTracks(this)
            val = this.cellTracks_val;
        end
        
        function val = get.cellLocs(this)
            val = this.cellLocs_val;
        end
        
        function val = get.cellLabels(this)
            val = this.cellLabels_val;
        end
        
        function val = get.times(this)
            val = this.times_val;
        end
        
        function val = get.currentTimepoint(this)
            val = this.currentTimepoint_val;
        end
        function set.currentTimepoint(this,val)
            if ~(isnumeric(val) && isscalar(val))
                error('timepoint must be a scalar numeric');
            end
            if ~ismember(val,1:this.ntimepoints)
                error('timepoint out of range');
            end
            if this.currentTimepoint_val==val, return; end
            this.currentTimepoint_val = val;
            % Ensure that drag tp is in a valid range
            nT = this.tilesContext;
            ctp = this.currentTimepoint;
            this.dragtp = min([this.dragtp,this.ntimepoints-ctp,nT]);
            this.dragtp = max([this.dragtp,1-ctp,-nT]);
            % Update display
            this.resetzoomcentre;
            this.refreshTrapIms;
        end
        
        function val = get.currentChannels(this)
            val = this.currentChannels_val;
        end
        function set.currentChannels(this,val)
            if ischar(val), val = {val}; end
            if ~iscellstr(val)
                error('channels must be specified by name');
            end
            if ~all(ismember(val,this.availableChannels))
                error('some channels could not be identified');
            end
            if length(this.currentChannels_val)==length(val) && ...
                    all(strcmp(this.currentChannels_val,val)), return; end
            this.currentChannels_val = val;
        end
        
        function val = get.currentCell(this)
            val = this.ctrls.cells.Value;
            if val==0, val = []; end
        end
        function set.currentCell(this,val)
            this.resetzoomcentre;
            if val==this.ctrls.cells.Value, return; end
            this.ctrls.cells.Value = val;
        end
        
        function val = get.layout(this)
            keys = fieldnames(this.layoutMap);
            if isfield(this.ctrls,'layout') && isvalid(this.ctrls.layout) ...
                    && ~isempty(this.ctrls.layout.Value)
                val = keys{this.ctrls.layout.Value};
            else
                if this.tilesContext<1
                    val = 'horizontal';
                else
                    val = 'vertical';
                end
            end
        end
        function set.layout(this,val)
            keys = fieldnames(this.layoutMap);
            assert(ismember(val,keys),'"layout" must be one of %s',...
                strjoin(strcat('"',keys,'"'),', '));
            this.ctrls.layout.Value = find(strcmp(keys,val),1);
        end
        
        function val = get.title(this)
            val = this.frame.fig.Name;
        end
        function set.title(this,val)
            this.frame.fig.Name = val;
        end
        
        function val = get.imzoom(this)
            val = this.imzoom_val;
        end
        function set.imzoom(this,val)
            assert(isnumeric(val) && isscalar(val),...
                '"imzoom" must be a scalar numeric');
            if this.imzoom_val==val, return; end
            this.imzoom_val = val;
            this.UpdateWindowElements;
            if ismember(val,this.zoomlevels)
                this.ctrls.zoom.Value = find(this.imzoom==this.zoomlevels,1);
            end
            this.refreshSegIms;
        end
        
        function val = get.zoomactive(this)
            if isempty(this.zoomactive_val)
                this.zoomactive_val = false;
            end
            if isfield(this.ctrls,'zoomactive') && isvalid(this.ctrls.zoomactive)
                val = this.ctrls.zoomactive.Value;
                if this.zoomactive_val~=val
                    this.zoomactive_val = val;
                    % Reset zoomcentre if zoom to active has been toggled
                    this.resetzoomcentre;
                end
            else
                val = this.zoomactive_val;
            end
        end
        function set.zoomactive(this,val)
            assert(islogical(val) && isscalar(val),...
                '"zoomactive" property must be either true or false');
            assert(isfield(this.ctrls,'zoomactive') && isvalid(this.ctrls.zoomactive),...
                '"zoomactive" control is not initialised');
            this.ctrls.zoomactive.Value = val;
            
            this.refreshSegIms;
        end
        
        function val = get.zoomcentre(this)
            % Trigger reset of centre if necessary
            this.zoomactive;
            if isempty(this.zoomcentre_val)
                [val,~,~] = this.getContourParams;
                if isempty(val)
                    val = [this.srcHeight/2,this.srcWidth/2];
                end
                this.zoomcentre_val = val;
            else
                val = this.zoomcentre_val;
            end
        end
        function set.zoomcentre(this,val)
            orig = this.zoomcentre_val;
            assert(isempty(val) || (all(size(val)==[1,2]) && isnumeric(val)),...
                '"zoomcentre" must be a numeric row vector of length 2');
            this.zoomcentre_val = val;
            if ~all(this.zoomcentre_val==orig)
                this.refreshSegIms;
            end
        end
        function resetzoomcentre(this)
            this.zoomcentre_val = [];
        end
        
        function val = get.zoomoffset(this)
            scale = this.zoomactivescale;
            tileSize = [this.tileWidth,this.tileHeight];
            val = round(scale*this.zoomcentre-tileSize/2);
            val(1) = max(val(1),0);
            val(1) = min(val(1),this.srcWidth*scale-this.tileWidth);
            val(2) = max(val(2),0);
            val(2) = min(val(2),this.srcHeight*scale-this.tileHeight);
        end
        
        function val = get.Enable(this)
            val = this.Enable_val;
        end
        function set.Enable(this,val)
            assert(ismember(val,{'off','on'}),...
                '"Enable" must be "off" or "on"');
            ctrlnames = fieldnames(this.ctrls);
            for c=1:length(ctrlnames)
                if isprop(this.ctrls.(ctrlnames{c}),'Enable')
                    this.ctrls.(ctrlnames{c}).Enable = val;
                end
            end
            this.frame.slider.Enable = val;
            if this.tileWidth==this.srcWidth
                this.frame.xslider.Enable = 'off';
            else
                this.frame.xslider.Enable = val;
            end
            if this.tileHeight==this.srcHeight
                this.frame.yslider.Enable = 'off';
            else
                this.frame.yslider.Enable = val;
            end
            this.refreshFocusAnnotStatus; % make sure focus annot is valid
            this.Enable_val = val;
        end
        
        function val = get.posImFromCache(this)
            if isfield(this.ctrls,'posImFromCache')
                val = this.ctrls.posImFromCache.Value;
            else
                val = false;
            end
        end
        
        function set.posImFromCache(this,val)
            assert(islogical(val) && isscalar(val),...
                '"posImFromCache" only takes true or false values');
            this.ctrls.posImFromCache.Value = val;
            if ishandle(this.ctrls.posImFromCache)
                this.refreshPosImage;
            end
        end
        
        function val = get.dofocusannot(this)
            if isfield(this.ctrls,'focusannot')
                val = this.ctrls.focusannot.Value;
            else
                val = false;
            end
        end
        function set.dofocusannot(this,val)
            assert(islogical(val) && isscalar(val),...
                '"dofocusannot" only takes true or false values');
            assert(isfield(this.ctrls,'focusannot'),...
                'controls not yet initialised');
            assert(strcmp(this.ctrls.focusannot.Enable,'on'),...
                'can only set focus annotation when control is enabled');
            this.ctrls.focusannot.Value = val;
            this.refreshSegIms;
        end
        
        function val = get.srcOffset(this)
            val = zeros(1,2);
            xsval = this.frame.xslider.Value;
            ysval = this.frame.yslider.Value;
            if this.tileWidth<this.srcWidth && ~isempty(xsval)
                val(1) = round(xsval*(this.srcWidth-this.tileWidth));
            end
            if this.tileHeight<this.srcHeight && ~isempty(ysval)
                val(2) = round(ysval*(this.srcHeight-this.tileHeight));
            end 
        end
        
        %% Marking functions
        function setMark(this,name,val,trap,tp,tag_on_omero)
            if nargin<6 || isempty(tag_on_omero), tag_on_omero = true; end
            if isempty(this.cExperiment.marks)
                this.cExperiment.marks = struct();
            end
            assert(isstruct(this.cExperiment.marks),...
                'cExperiment.marks has been corrupted; cannot add marks');
            
            pos = this.currentPos;
            
            if ~isfield(this.cExperiment.marks,name)
                this.cExperiment.marks.(name) = ...
                    cell(length(this.cExperiment.dirs),1);
            end
            isMarked = this.cExperiment.marks.(name);
            assert(iscell(isMarked) && length(isMarked)>=pos,...
                'Cannot change mark: cExperiment.marks.%s has been corrupted',...
                name);
            
            isMarked = isMarked{pos};
            if isempty(isMarked)
                isMarked = zeros(0,2,'uint16');
            end
            
            if ismember([trap,tp],isMarked,'rows')
                if ~val
                    isMarked(ismember(isMarked,[trap,tp],'rows'),:) = [];
                    this.cExperiment.marks.(name){pos} = isMarked;
                    this.expthaschanged = true;
                end
            else
                if val
                    isMarked(end+1,:) = [trap,tp];
                    this.cExperiment.marks.(name){pos} = isMarked;
                    this.expthaschanged = true;
                    if isa(this.cExperiment,'babyExperimentOmero') ...
                            && tag_on_omero
                        this.cExperiment.OmeroDatabase.addTag(...
                            name,this.cExperiment.omeroDs,'marktags');
                    end
                end
            end
        end
        
        function val = getMarks(this,name,trap,tp)
            if nargin<3, val = zeros(0,2,'uint16');
            elseif nargin<4, val = [];
            else, val = false;
            end
            
            if isempty(this.cExperiment.marks), return; end
            
            pos = this.currentPos;
            
            if ~isfield(this.cExperiment.marks,name)
                return
            end
            isMarked = this.cExperiment.marks.(name);
            if ~iscell(isMarked) || length(isMarked)<pos
                return
            end
            
            isMarked = isMarked{pos};
            if isempty(isMarked), return; end
            
            if nargin<3, val = isMarked;
            elseif nargin<4
                val = isMarked(ismember(isMarked(:,1),trap),2);
            else, val = ismember([trap,tp],isMarked,'rows');
            end
        end
        
        function markForCuration(this)
            this.setMark('forCuration',true,this.currentTrap,this.currentTimepoint);
            this.updateTrapList;
            this.refreshTrackPlot;
        end
        function unmarkForCuration(this)
            this.setMark('forCuration',false,this.currentTrap,this.currentTimepoint);
            this.updateTrapList;
            this.refreshTrackPlot;
        end
        function val = forCuration(this,trap,tp)
            if nargin<2, val = this.getMarks('forCuration');
            elseif nargin<3, val = this.getMarks('forCuration',trap);
            else, val = this.getMarks('forCuration',trap,tp);
            end
        end
        
        function markAsCurated(this)
            this.setMark('isCurated',true,this.currentTrap,this.currentTimepoint);
            this.updateTrapList;
            this.refreshTrackPlot;
        end
        function unmarkAsCurated(this)
            this.setMark('isCurated',false,this.currentTrap,this.currentTimepoint);
            this.updateTrapList;
            this.refreshTrackPlot;
        end
        function val = isCurated(this,trap,tp)
            if nargin<2, val = this.getMarks('isCurated');
            elseif nargin<3, val = this.getMarks('isCurated',trap);
            else, val = this.getMarks('isCurated',trap,tp);
            end
        end
    end
    
    methods (Access=private)
        function saveOriginalOutlines(this,trap,tp)
            trapInfo = this.cTimelapse.cTimepoint(tp).trapInfo(trap);
            % Check if we need to save the segmentation outline
            if ~isfield(trapInfo,'cellOriginal') ...
                    || ~isstruct(trapInfo.cellOriginal) ...
                    || ~isfield(trapInfo,'cellLabelOriginal')
                if trapInfo.cellsPresent
                    this.cTimelapse.cTimepoint(tp).trapInfo(trap).cellOriginal = trapInfo.cell;
                    this.cTimelapse.cTimepoint(tp).trapInfo(trap).cellLabelOriginal = trapInfo.cellLabel;
                else
                    this.cTimelapse.cTimepoint(tp).trapInfo(trap).cellOriginal = ...
                        this.cTimelapse.cellInfoTemplate([]);
                    this.cTimelapse.cTimepoint(tp).trapInfo(trap).cellLabelOriginal = [];
                end
            end
        end
        
        function markCellOutlineEdit(this,trap,tp,ci)
            modifiedon = now;
            if ~isempty(ci)
                this.cTimelapse.cTimepoint(tp).trapInfo(trap).cell(ci).modifiedon = modifiedon;
                % Reset any assigned autogen attributes
                if isfield(this.cTimelapse.cTimepoint(tp).trapInfo(trap).cell,'autogen')
                    this.cTimelapse.cTimepoint(tp).trapInfo(trap).cell(ci).autogen = false;
                end
            end
            % Reset any cached area calculations
            this.clearCellAreas(trap,tp);
            this.cTimelapse.cTimepoint(tp).trapInfo(trap).modifiedon = modifiedon;
            this.haschanged = true;
        end
        
        function ci = ensureCellIndex(this,trap,tp,cellLabel)
            %babyGUI.ensureCellIndex Ensures cell and template exist
            %   CellIndex = egui.ensureCellIndex(trap,tp,cellLabel)
            %   - ci: index to use into the trapInfo.cell array
            %   - trap: trap index
            %   - tp: time point index (should be obtained from this.tpInds)
            %   - cellLabel: label of the cell
            
            % Ensure original outlines from segmentation are saved before
            % adding any cell indices:
            this.saveOriginalOutlines(trap,tp);
            
            trapInfo = this.cTimelapse.cTimepoint(tp).trapInfo(trap);
            
            if trapInfo.cellsPresent
                ci = find(trapInfo.cellLabel==cellLabel,1);
                if isempty(ci)
                    ci = length(trapInfo.cell)+1;
                    this.cTimelapse.cTimepoint(tp).trapInfo(trap).cellLabel(ci) = cellLabel;
                end
            else
                ci = 1;
                this.cTimelapse.cTimepoint(tp).trapInfo(trap).cell = ...
                    this.cTimelapse.cellInfoTemplate;
                this.cTimelapse.cTimepoint(tp).trapInfo(trap).cellsPresent = true;
                this.cTimelapse.cTimepoint(tp).trapInfo(trap).cellLabel(ci) = cellLabel;
            end
        end
        
        function [zoff,zsc] = getZoomActiveParams(this)
            if this.zoomactive
                zoff = this.zoomoffset; zsc = this.zoomactivescale;
            else
                zoff = this.srcOffset; zsc = 1;
            end
        end
        
        function val = isVisible(this,XY)
            [zoff,zsc] = this.getZoomActiveParams;
            XY = XY - zoff/zsc;
            tW = this.tileWidth/zsc;
            tH = this.tileHeight/zsc;
            val = all(XY>=1,2) & XY(:,1)<=tW & XY(:,2)<=tH;
        end
        
        function val = isDragVisible(this)
            val = false;
            if isempty(this.centre) || isempty(this.radii) || isempty(this.angles)
                return
            end
            N = numel(this.radii);
            XY = this.centre(ones(N,1),:)+this.radii(:,ones(1,2)).*...
                [cos(this.angles),sin(this.angles)];
            val = any(this.isVisible(XY));
        end
        
        function [centre,radii,angles] = getContourParams(this)
            centre = []; radii = []; angles = [];
            ct = this.tpInds(this.currentTimepoint+this.dragtp);
            trapInfo = this.cTimelapse.cTimepoint(ct).trapInfo(this.currentTrap);
            if isempty(this.currentCell), return; end
            if ~trapInfo.cellsPresent, return; end
            cellLabel = this.cellLabels{this.currentTrap}(this.currentCell);
            cc = find(trapInfo.cellLabel==cellLabel,1);
            if isempty(cc), return; end
            centre = reshape(trapInfo.cell(cc).cellCenter,1,[]);
            radii = reshape(trapInfo.cell(cc).cellRadii,[],1);
            angles = reshape(trapInfo.cell(cc).cellAngle,[],1);
        end
        
        function [px,py] = setContourParams(this,centre,radii,angles,trap,tp,ci)
            %babyGUI.setContourParams Set a new cell outline
            %   NB: this function assumes that trapInfo.cell has been
            %   initialised correctly and that the cell index has a
            %   corresponding trapInfo.cellLabel
            if nargin<5, trap = this.currentTrap; end
            if nargin<6, tp = this.tpInds(this.currentTimepoint+this.dragtp); end
            
            % Ensure original outlines from segmentation are saved before
            % editing:
            this.saveOriginalOutlines(trap,tp);
            
            trapInfo = this.cTimelapse.cTimepoint(tp).trapInfo(trap);
            if nargin<7
                assert(trapInfo.cellsPresent,'Cannot set contour when cell not initialised');
                cellLabel = this.cellLabels{this.currentTrap}(this.currentCell);
                ci = find(trapInfo.cellLabel==cellLabel,1);
                assert(~isempty(ci),'Cannot set contour when cell not initialised');
            end
            
            % Generate new segmentation outline
            image_size = [this.srcHeight,this.srcWidth];
            if this.ctrls.cartesianspline.Value
                [px,py] = BABYutil.cartesian_spline_from_radii(...
                    double(radii),double(angles),double(centre),image_size);
            else
                [px,py] = BABYutil.get_full_points_from_radii(...
                    double(radii),double(angles),double(centre),image_size);
            end
            trapInfo.cell(ci).segmented = ...
                sparse(BABYutil.px_py_to_logical(px,py,image_size));
            trapInfo.cell(ci).cellCenter = reshape(centre,1,[]);
            trapInfo.cell(ci).cellRadii = reshape(radii,1,[]);
            trapInfo.cell(ci).cellAngle = reshape(angles,1,[]);
            trapInfo.cell(ci).cellRadius = mean(trapInfo.cell(ci).cellRadii);
            this.cTimelapse.cTimepoint(tp).trapInfo(trap) = trapInfo;
            this.markCellOutlineEdit(trap,tp,ci);
        end
        
        function refreshFocusAnnotStatus(this)
            % Can only select focus annotation if there are context tiles, and if all
            % the brightfield stacks are displayed
            if ~isfield(this.ctrls,'focusannot'), return; end
            enablefocusannot = false;
            if any(strcmpi(this.imcache.channels,'brightfield'))
                bfchans = this.imcache.channels;
                bfpre = 'brightfield_';
                bfchans = bfchans(strncmpi(bfchans,bfpre,numel(bfpre)));
                enablefocusannot = this.tilesContext>0 && (isempty(bfchans) || ...
                    all(ismember(bfchans,this.currentChannels)));
            end
            if enablefocusannot
                this.ctrls.focusannot.Enable = 'on';
            else
                this.ctrls.focusannot.Enable = 'off';
                this.ctrls.focusannot.Value = false;
            end
        end
        
        function isfachan = isFocusAnnotationChannel(this)
            % There should either be a single Brightfield channel, or multiple
            % numbered Brightfield channels corresponding to the stacks
            bfpre = 'brightfield_';
            isfachan = strncmpi(this.currentChannels,bfpre,numel(bfpre));
            if ~any(isfachan)
                isfachan = strcmpi(this.currentChannels,'brightfield');
                assert(any(isfachan),'could not find the brightfield channels');
            end
        end
        
        adjustNormGUI(this)
    end
end
