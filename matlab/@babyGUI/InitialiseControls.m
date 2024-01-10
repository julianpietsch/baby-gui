function InitialiseControls(this)
%InitialiseControls Initialise the controls for a CellAnnotationGUI

%% Retrieve default dimensions
dims = this.dimensions;
sp = dims.sp; % Spacing between controls
lh = dims.lh; % Label height
sh = dims.sh; % Selection menu height
ch = dims.ch; % Control height

%% Overall control panel layout

% Split into two columns
mainbox = uix.HBox('Parent',this.frame.ctrl,'Padding',0,'Spacing',sp);
leftbox = uix.VBox('Parent',mainbox,'Padding',0,'Spacing',sp);
rightbox = uix.VBox('Parent',mainbox,'Padding',0,'Spacing',sp);

%% Left column layout

% Edit mode (add/remove outline; edit outline;
% set active track; track editing; cell selection)
uicontrol('Parent',leftbox,'Style','text','FontWeight','bold',...
    'HorizontalAlignment','left','String','Left/right click actions:');
this.ctrls.editmode = uicontrol('Parent',leftbox,'Style','listbox',...
    'String',struct2cell(this.editModeMap(1)),...
    'Value',1,'Callback',@updateEditMode,'Enable',this.Enable);
modekeys = struct2cell(this.editModeMap(2));

editmodehbox = uix.HBox('Parent',leftbox);
uicontrol('Parent',editmodehbox,'Style','text','FontWeight','bold',...
    'HorizontalAlignment','left','String','Edit mode:');
this.ctrls.outlinemode = uicontrol('Parent',editmodehbox,'Style','popupmenu',...
    'String',struct2cell(this.outlineModeMap(1)),...
    'Value',1,'Callback',@updateEditMode,'Enable',this.Enable);
outlinemodekeys = struct2cell(this.outlineModeMap(2));
set(editmodehbox,'Widths',[-2,-3]);

addmodehbox = uix.HBox('Parent',leftbox);
uicontrol('Parent',addmodehbox,'Style','text','FontWeight','bold',...
    'HorizontalAlignment','left','String','Add mode:');
this.ctrls.addmode = uicontrol('Parent',addmodehbox,'Style','popupmenu',...
    'String',struct2cell(this.addModeMap(1)),...
    'Value',1,'Enable',this.Enable);
addmodekeys = struct2cell(this.addModeMap(2));
set(addmodehbox,'Widths',[-2,-3]);

copyoutlinehbox = uix.HBox('Parent',leftbox);
uicontrol('Parent',copyoutlinehbox,'Style','text','FontWeight','bold',...
    'HorizontalAlignment','left','String','Copy from:');
this.ctrls.copyfromprev = uicontrol('Parent',copyoutlinehbox,...
    'Style','pushbutton','String','prev (C)','Enable',this.Enable,...
    'Callback',@(~,~) copyoutline(-1));
this.ctrls.copyfromnext = uicontrol('Parent',copyoutlinehbox,...
    'Style','pushbutton','String','next (c)','Enable',this.Enable,...
    'Callback',@(~,~) copyoutline(+1));
set(copyoutlinehbox,'Widths',[-4,-3,-3]);

this.ctrls.resetoutline = uicontrol('Parent',leftbox,...
    'Style','pushbutton','String','Reset outline to segmented',...
    'Enable',this.Enable,'Callback',@resetoutline);

viewopthbox = uix.HBox('Parent',leftbox);
this.ctrls.showoutlines = uicontrol('Parent',viewopthbox,'Style','checkbox',...
    'String','Show outlines (f)','Value',true,'Enable',this.Enable,...
    'CallBack',@showoutlines);
this.ctrls.zoomactive = uicontrol('Parent',viewopthbox,'Style','checkbox',...
    'String','Zoom to active (w)','Value',this.zoomactive,...
    'Enable',this.Enable,'Callback',@(~,~) this.refreshSegIms);

nknotsbox = uix.HBox('Parent',leftbox);
uicontrol('Parent',nknotsbox,'Style','text','FontWeight','bold',...
    'HorizontalAlignment','left','String','Number of knots:');
this.ctrls.nknots = uicontrol('Parent',nknotsbox,'Style','popupmenu',...
    'String',[{'default'},arrayfun(@(x) sprintf('%u',x),this.nknot_values,'Uni',0)],...
    'Value',1,'Enable',this.Enable,'Callback',@setnknots);
set(nknotsbox,'Widths',[-3,-2]);

% Save cTimelapse and save cExperiment buttons
savehbox = uix.HBox('Parent',leftbox);
this.ctrls.saveTimelapse = uicontrol('Parent',savehbox,'Style','pushbutton',...
    'String','Save cTimelapse...','Callback',@saveTimelapseClick,'Enable',this.Enable);
this.ctrls.saveExperiment = uicontrol('Parent',savehbox,'Style','pushbutton',...
    'String','Save cExperiment...','Callback',@saveExperimentClick,'Enable',this.Enable);

set(leftbox,'Heights',[lh,-1,sh,sh,0.9*ch,ch,lh,sh,1.1*ch]);

%% Right column tabbed layout

% Use tabs for actions and options
ctrlTabs = uix.TabPanel('Parent',rightbox,'Padding',3);
navCtrls = uix.VBox('Parent',ctrlTabs,'Padding',0,'Spacing',sp);
optionCtrls = uix.VBox('Parent',ctrlTabs,'Padding',0,'Spacing',sp);
annotCtrls = uix.VBox('Parent',ctrlTabs,'Padding',0,'Spacing',sp);
set(ctrlTabs,'TabTitles',{'Navigate','Options','Annotate'},...
    'TabWidth',dims.tabWidth);

%% Navigate tab

% Pos dropdown
this.ctrls.pos = uicontrol('Parent',navCtrls,'Style','popupmenu',...
    'String',this.cExperiment.dirs,'Value',this.currentPos,...
    'Callback',@updatePos,'Enable',this.Enable);

% Trap selection list
this.ctrls.traps = uicontrol('Parent',navCtrls,'Style','listbox',...
    'String',{},'Value',[],'Callback',@updateTrap,'Enable',this.Enable);
updateTrapList;

% Cell selection list
this.ctrls.cells = uicontrol('Parent',navCtrls,'Style','listbox',...
    'String',{},'Value',[],'Callback',@updateCell,'Enable',this.Enable);
updateCellList;

this.ctrls.jumpToLastModified = uicontrol('Parent',navCtrls,...
    'Style','pushbutton','String','Jump to last modified...',...
    'Callback',@jumpToLastModified,'Enable',this.Enable);

% Image cache panel (cache size will be displayed in title)
this.ctrls.cachepanel = uix.Panel('Parent',navCtrls,'Title','Image cache');
updateImageCacheSize;
cachehbox = uix.HBox('Parent',this.ctrls.cachepanel);
this.ctrls.preload = uicontrol('Parent',cachehbox,'Style','pushbutton',...
    'String','Preload...','Callback',@preloadPosition,'Enable',this.Enable);
this.ctrls.clearcache = uicontrol('Parent',cachehbox,'Style','pushbutton',...
    'String','Clear...','Enable',this.Enable,'Callback',@clearcache);

% Local image cache panel
this.ctrls.localcachepanel = uix.Panel('Parent',navCtrls,...
    'Title','Local image cache');
cachehbox = uix.HBox('Parent',this.ctrls.localcachepanel);
this.ctrls.createcache = uicontrol('Parent',cachehbox,...
    'Style','pushbutton','String','Create...',...
    'Callback',@(~,~) this.imcache.createLocalCache(),'Enable',this.Enable);
this.ctrls.loadcache = uicontrol('Parent',cachehbox,...
    'Style','pushbutton','String','Load...',...
    'Callback',@(~,~) this.imcache.openLocalCache(),'Enable',this.Enable);
this.ctrls.savecache = uicontrol('Parent',cachehbox,...
    'Style','pushbutton','String','Save...',...
    'Callback',@saveLocalCache,'Enable',this.Enable);

    function saveLocalCache(~,~)
        this.imcache.saveLocal(this.currentPos,this.currentChannels);
    end

set(navCtrls,'Heights',[sh,-1,-1,ch,lh+ch,lh+ch]);

%% Options tab

% Toggle horizontal layout
this.ctrls.layout = uicontrol('Parent',optionCtrls,'Style','popupmenu',...
    'String',struct2cell(this.layoutMap),'Value',[],...
    'Enable',this.Enable,'Callback',@changeViewMode);
this.layout = this.layout; % set layout to its default value

% Zoom level
this.ctrls.zoom = uicontrol('Parent',optionCtrls,'Style','popupmenu',...
    'String',arrayfun(@(x) sprintf('%0.1fx zoom',x),this.zoomlevels,'Uni',0),...
    'Value',find(this.imzoom==this.zoomlevels,1),...
    'Callback',@updateZoom,'Enable',this.Enable);

% Set number of time points of context
this.ctrls.tpcontext = uicontrol('Parent',optionCtrls,'Style','popupmenu',...
    'String',arrayfun(@(x) sprintf('+%u/-%u context time points',x,x),0:8,'Uni',0),...
    'Value',find(this.tilesContext==0:8),...
    'Callback',@updateTilesContext,'Enable',this.Enable);

channelbox = uix.HBox('Parent',optionCtrls);
this.ctrls.channelSel = uicontrol('Parent',channelbox,'Style','pushbutton',...
    'String','Trap channels...','Callback',@channelSelect,'Enable',this.Enable);
this.ctrls.posChannelSel = uicontrol('Parent',channelbox,'Style','pushbutton',...
    'String','Pos channel...','Callback',@posChannelSelect,'Enable',this.Enable);

normbox = uix.HBox('Parent',optionCtrls);
this.ctrls.autoNormalisation = uicontrol('Parent',normbox,'Style','pushbutton',...
    'String','Auto norm...','Callback',@autoNormalisation,'Enable',this.Enable);
this.ctrls.adjustNormalisation = uicontrol('Parent',normbox,'Style','pushbutton',...
    'String','Adjust norm...','Callback',@adjustNormalisation,'Enable',this.Enable);

% Track display options
trackdisp_panel = uix.Panel('Parent',optionCtrls,...
    'Title','Track display options');
trackdisp_hbox = uix.HBox('Parent',trackdisp_panel);
this.ctrls.trackdisp = uicontrol('Parent',trackdisp_hbox,'Style','popupmenu',...
    'String',struct2cell(this.trackDisplayMap),'Value',1,...
    'Enable',this.Enable,'CallBack',@this.refreshTrackPlot);
this.ctrls.trackdispsz = uicontrol('Parent',trackdisp_hbox,'Style','popupmenu',...
    'String',struct2cell(this.trackSizeMap(1)),'Value',1,...
    'Enable',this.Enable,'CallBack',@changeViewMode);

% Update overview checkbox
this.ctrls.updateoverview = uicontrol('Parent',optionCtrls,'Style','checkbox',...
    'String','Update pos image','Value',true,'Enable',this.Enable);

% Update overview while scrolling checkbox
this.ctrls.scrolloverview = uicontrol('Parent',optionCtrls,'Style','checkbox',...
    'String','Scroll pos image','Value',false,'Enable',this.Enable);

this.ctrls.posImFromCache = uicontrol('Parent',optionCtrls,'Style','checkbox',...
    'String','Make pos image from cache','Enable',this.Enable,...
    'Value',this.ctrls.posImFromCache.Value,'CallBack',@this.refreshPosImage);

this.ctrls.spcreverse = uicontrol('Parent',optionCtrls,'Style','checkbox',...
    'String','Reverse space bar (r)','Value',false,'Enable',this.Enable);

this.ctrls.orphantracks = uicontrol('Parent',optionCtrls,'Style','checkbox',...
    'String','Display orphan tracks (o)','Value',true,'Enable',this.Enable,...
    'CallBack',@this.refreshTrackPlot);

this.ctrls.showlinarrows = uicontrol('Parent',optionCtrls,'Style','checkbox',...
    'String','Display lineage arrows (l)','Value',true,...
    'Enable',this.Enable,'CallBack',@(~,~) this.refreshSegIms());

this.ctrls.refreshtracks = uicontrol('Parent',optionCtrls,'Style','checkbox',...
    'String','Refresh track plot (T)','Value',true,...
    'Enable',this.Enable,'CallBack',@this.refreshTrackPlot);

this.ctrls.cartesianspline = uicontrol('Parent',optionCtrls,...
    'Style','checkbox','Enable',this.Enable,...
    'String','Use cartesian spline','Value',false);

uix.HBox('Parent',optionCtrls);

set(optionCtrls,'Heights',[sh,sh,sh,ch,ch,lh+sh,lh,lh,lh,lh,lh,lh,lh,lh,-1]);

%% Annotate tab

uicontrol('Parent',annotCtrls,'Style','text','FontWeight','bold',...
    'HorizontalAlignment','left','String','Save stacks to:');
this.saveStackDir = BABYutil.uiFile(annotCtrls,'/',pwd,pwd);
this.ctrls.saveStackDir = this.saveStackDir.btn;

savestackbox = uix.HBox('Parent',annotCtrls);
this.ctrls.saveChannel = uicontrol('Parent',savestackbox,'Style','popupmenu',...
    'String',this.cTimelapse.channelNames,'Value',1,'Enable',this.Enable);
this.ctrls.saveStack = uicontrol('Parent',savestackbox,'Style','pushbutton',...
    'String','Save stack...','Callback',@saveStack,'Enable',this.Enable);

forCurationBox = uix.HBox('Parent',annotCtrls);
this.ctrls.markForCuration = uicontrol('Parent',forCurationBox,...
    'Style','pushbutton','String','Mark',...
    'Callback',@(~,~) this.markForCuration);
uicontrol('Parent',forCurationBox,'Style','text','FontWeight','bold',...
    'HorizontalAlignment','center','String','/','Position',[0,0,1,0.8*ch]);
this.ctrls.unmarkForCuration = uicontrol('Parent',forCurationBox,...
    'Style','pushbutton','String','Unmark',...
    'Callback',@(~,~) this.unmarkForCuration);
uicontrol('Parent',forCurationBox,'Style','text','FontWeight','bold',...
    'HorizontalAlignment','left','String','for curation.');
set(forCurationBox,'Widths',[-5,-1,-5,-7]);

asCuratedBox = uix.HBox('Parent',annotCtrls);
this.ctrls.markAsCurated = uicontrol('Parent',asCuratedBox,...
    'Style','pushbutton','String','Mark',...
    'Callback',@(~,~) this.markAsCurated);
uicontrol('Parent',asCuratedBox,'Style','text','FontWeight','bold',...
    'HorizontalAlignment','center','String','/');
this.ctrls.unmarkAsCurated = uicontrol('Parent',asCuratedBox,...
    'Style','pushbutton','String','Unmark',...
    'Callback',@(~,~) this.unmarkAsCurated);
uicontrol('Parent',asCuratedBox,'Style','text','FontWeight','bold',...
    'HorizontalAlignment','left','String','as curated.');
set(asCuratedBox,'Widths',[-5,-1,-5,-7]);

navCuratedBox = uix.HBox('Parent',annotCtrls);
this.ctrls.prevIsCurated = uicontrol('Parent',navCuratedBox,...
    'Style','pushbutton','String','Prev','Callback',@prevIsCurated);
uicontrol('Parent',navCuratedBox,'Style','text','FontWeight','bold',...
    'HorizontalAlignment','center','String','/');
this.ctrls.nextIsCurated = uicontrol('Parent',navCuratedBox,...
    'Style','pushbutton','String','Next','Callback',@nextIsCurated);
uicontrol('Parent',navCuratedBox,'Style','text','FontWeight','bold',...
    'HorizontalAlignment','left','String','curated.');
set(navCuratedBox,'Widths',[-5,-1,-5,-7]);

%%%% babyExperimentSamples not yet available in baby-gui version %%%%%
% this.ctrls.newcEsamples = uicontrol('Parent',annotCtrls,...
%     'Style','pushbutton','String','Create new cEsamples...',...
%     'Callback',@(~,~) babyExperimentSamples);
% 
% this.ctrls.addTocEsamples = uicontrol('Parent',annotCtrls,...
%     'Style','pushbutton','String','Add to cEsamples...',...
%     'Callback',@(~,~) babyExperimentSamples.addExperimentGUI([],this.cExperiment));

this.ctrls.focusannot = uicontrol('Parent',annotCtrls,'Style','checkbox',...
    'String','Focus annotation mode','Value',false,'Enable',this.Enable,...
    'Callback',@(~,~) this.refreshSegIms);
this.ctrls.nextNoFocusAssigned = uicontrol('Parent',annotCtrls,...
    'Style','pushbutton','String','Next with unassigned focus (n)...',...
    'Callback',@jumpToNoFocusAssigned);
this.refreshFocusAnnotStatus;

cellsel_panel = uix.Panel('Parent',annotCtrls,...
    'Title','Cell selection for this trap');
cellselbox = uix.HBox('Parent',cellsel_panel);
this.ctrls.selectAllLineages = uicontrol('Parent',cellselbox,...
    'Style','pushbutton','String','Select all lineages',...
    'Callback',@selectAllLineages);
this.ctrls.unselectAllCells = uicontrol('Parent',cellselbox,...
    'Style','pushbutton','String','Unselect all',...
    'Callback',@unselectAllCells);

uix.HBox('Parent',annotCtrls);

set(annotCtrls,'Heights',[lh,ch,ch,0.8*ch,0.8*ch,0.8*ch,sh,ch,ch+sh,-1]);

%% Slider control
set(this.frame.slider,'Callback',@this.refreshPosImage);
addlistener(this.frame.slider,'Value','PostSet',@slider_update);
addlistener(this.frame.xslider,'Value','PostSet',@(~,~) this.refreshSegIms);
addlistener(this.frame.yslider,'Value','PostSet',@(~,~) this.refreshSegIms);

%% Change mode with keyboard clicks
set(this.frame.fig,'WindowKeyPressFcn',@keypress);

%% Listen to mouse click events

% The following is required to ensure that properties (e.g.,
% ButtonDownFcn) do not get reset with plotting commands
hold(this.frame.imAx,'on');
hold(this.frame.pos,'on');

set(this.frame.im,'ButtonDownFcn',@imclick);
set(this.frame.pos.Children,'ButtonDownFcn',@posImClick);

%% GUI callbacks

    function saveTimelapseClick(~,~)
        this.haschanged = true;
        this.saveTimelapse;
    end

    function saveExperimentClick(~,~)
        this.expthaschanged = true;
        this.saveExperiment;
    end

    function imclick(~,~)
        % Ignore click if inactive
        if strcmp(this.Enable,'off'), return; end
        % Find time point and x,y location of click
        clickloc = get(this.frame.imAx,'CurrentPoint');
        clicktype = get(this.frame.fig,'SelectionType');
        if ~ismember(clicktype,{'normal','alt'}), return; end
        % Perform edits according to mode
        this.EditTracking(clickloc(1,1),clickloc(1,2),strcmp(clicktype,'normal'));
        % Delegate CellLabel/Track refresh to EditTracking for optimisation...
        % this.refreshCellLabels;
        % this.refreshTracks;
        updateTrapList;
        updateCellList;
        updateCell;
    end

    function posImClick(~,~)
        % Ignore click if inactive
        if strcmp(this.Enable,'off'), return; end
        % Find x,y location of click
        clickloc = get(this.frame.pos,'CurrentPoint');
        % Find nearest trap at the current time point
        tpi = this.tpInds(this.currentTimepoint);
        trapLocs = this.cTimelapse.cTimepoint(tpi).trapLocations;
        trapLocs = [[trapLocs.xcenter]',[trapLocs.ycenter]'];
        clickloc = clickloc(ones(size(trapLocs,1),1),1:2);
        trapDist = sum((trapLocs-clickloc).^2,2);
        [~,closestTrap] = min(trapDist);
        this.currentTrap = closestTrap;
        this.ctrls.traps.Value = this.currentTrap;
        this.refreshTrapIms;
        this.refreshPosImage;
        updateCell;
        updateCellList;
    end

    function slider_update(~,~)
        slider_val = get(this.frame.slider,'Value');
        this.currentTimepoint = round(slider_val*(this.ntimepoints-1) + 1);
        this.refreshSegIms;
        if this.ctrls.scrolloverview.Value, this.refreshPosImage(true); end
        this.trackMarkerLine.XData = repmat(this.times(this.currentTimepoint),2,1);
    end

    function setSlider()
        slider_val = (this.currentTimepoint-1)/(this.ntimepoints-1);
        set(this.frame.slider,'Value',slider_val);
    end

    function updateEditMode(~,~)
        % Change seg colours to match mode
        updateCell;
    end

    function keypress(~,event)
        % Ignore keypress if inactive
        if strcmp(this.Enable,'off'), return; end
        if isa(event.Source,'UIControl')
            % Only run default keyboard actions when a control has focus
            return
        end
        switch event.Key
            case 'leftarrow'
                % decrement timepoint
                if this.currentTimepoint > 1
                    this.currentTimepoint = this.currentTimepoint-1;
                    setSlider; slider_update;
                    this.refreshPosImage;
                end
            case 'rightarrow'
                % increment timepoint
                if this.currentTimepoint < this.ntimepoints
                    this.currentTimepoint = this.currentTimepoint+1;
                    setSlider; slider_update;
                    this.refreshPosImage;
                end
            case 'space'
                if this.ctrls.spcreverse.Value
                    % decrement timepoint
                    if this.currentTimepoint > 1
                        this.currentTimepoint = this.currentTimepoint-1;
                        setSlider; slider_update;
                        this.refreshPosImage;
                    end
                else
                    % increment timepoint
                    if this.currentTimepoint < this.ntimepoints
                        this.currentTimepoint = this.currentTimepoint+1;
                        setSlider; slider_update;
                        this.refreshPosImage;
                    end
                end
            case 't'
                if isequal(event.Modifier,{'shift'})
                    % Toggle refresh tracks
                    this.ctrls.refreshtracks.Value = ~this.ctrls.refreshtracks.Value;
                    this.refreshTrackPlot;
                else
                    % increment trap
                    if this.currentTrap < this.ntraps
                        this.currentTrap = this.currentTrap+1;
                        this.ctrls.traps.Value = this.currentTrap;
                        this.refreshTrapIms;
                        this.refreshPosImage;
                        updateCell;
                        updateCellList;
                    end
                end
            case 'f'
                % Toggle view cell outlines
                this.ctrls.showoutlines.Value = ~this.ctrls.showoutlines.Value;
                this.refreshSegIms;
                this.refreshPosImSeg(this.currentTrap);
            case 'w'
                % Toggle zoom to active
                this.zoomactive = ~this.ctrls.zoomactive.Value;
            case 'r'
                if isequal(event.Modifier,{'shift'})
                    % Toggle reverse spacebar
                    this.ctrls.spcreverse.Value = ~this.ctrls.spcreverse.Value;
                else
                    % Add outlines as rods
                    this.ctrls.addmode.Value = find(strcmp(addmodekeys,'r'),1);
                end
            case 'o'
                this.ctrls.orphantracks.Value = ~this.ctrls.orphantracks.Value;
                this.refreshTrackPlot;
            case 'm'
                mInd = find(strcmp(fieldnames(this.trackDisplayMap),'mdarea'));
                if this.ctrls.trackdisp.Value == mInd
                    this.ctrls.trackdisp.Value = 1;
                else
                    this.ctrls.trackdisp.Value = mInd;
                end
                this.refreshTrackPlot;
            case 'l'
                this.ctrls.showlinarrows.Value = ~this.ctrls.showlinarrows.Value;
                this.refreshSegIms;
            case '0'
                if strcmp(event.Character,'0')
                    % Change zoom to minimum
                    this.ctrls.zoom.Value = 1;
                    updateZoom(this.ctrls.zoom);
                end
            case '6'
                this.ctrls.nknots.Value = find(this.nknot_values==6,1)+1;
                setnknots;
            case '8'
                this.ctrls.nknots.Value = find(this.nknot_values==8,1)+1;
                setnknots;
            case 'c'
                if event.Character=='C'
                    copyoutline(-1);
                else
                    copyoutline(+1);
                end
            case 'n'
                jumpToNoFocusAssigned;
            otherwise
                if ismember(event.Character,modekeys)
                    this.ctrls.editmode.Value = find(strcmp(modekeys,event.Character),1);
                    updateEditMode;
                end
                if ismember(event.Character,outlinemodekeys)
                    this.ctrls.outlinemode.Value = find(strcmp(outlinemodekeys,event.Character),1);
                    updateEditMode;
                end
                if ismember(event.Character,addmodekeys)
                    this.ctrls.addmode.Value = find(strcmp(addmodekeys,event.Character),1);
                end
        end
    end

    function showoutlines(~,~)
        this.refreshSegIms;
        this.refreshPosImSeg(this.currentTrap);
    end

    function changeViewMode(~,~)
        this.UpdateWindowElements(true);
        this.refreshSegIms;
        
    end

    function channelSelect(~,~)
        currentChannelInds = ...
            find(ismember(this.availableChannels,this.currentChannels));
        [channels,ok] = listdlg('ListString',this.availableChannels,...
            'SelectionMode','multiple','Name','Pick channels to display',...
            'ListSize',[200,150],'InitialValue',currentChannelInds);
        if ~ok, return; end
        this.currentChannels = this.availableChannels(channels);
        this.UpdateWindowElements;
        this.refreshTrapIms;
        this.refreshSegIms;
        updateImageCacheSize;
    end

    function posChannelSelect(~,~)
        [channels,ok] = listdlg('ListString',this.cTimelapse.channelNames,...
            'SelectionMode','single','Name','Pick a pos channel to display',...
            'ListSize',[200,150],'InitialValue',this.posImChannel);
        if ~ok, return; end
        this.posImChannel = channels;
        this.refreshPosImage;
    end

    function updateZoom(x,~)
        this.imzoom = this.zoomlevels(x.Value);
    end

    function updateTilesContext(x,~)
        this.tilesContext = x.Value-1;
        this.UpdateWindowElements;
        this.refreshTrapIms;
        this.refreshSegIms;
    end

    function preloadPosition(~,~)
        this.Enable = 'off';
        drawnow;
        resetProgressBar;
        this.progressBar.frame.setTitle('Preloading images...');
        nchannels = length(this.currentChannels);
        this.progressBar.push_bar('Channel...',1,nchannels);
        for c=1:nchannels
            try
                this.progressBar.set_val(c);
                this.imcache.loadPosChannel(this.currentPos,...
                    this.currentChannels{c},[],this.progressBar);
            catch err
                errordlg(err.message);
                this.progressBar.pop_bar;
                break
            end
        end
        this.progressBar.pop_bar;
        this.Enable = 'on';
        drawnow;
    end

    function updatePos(x,~)
        this.Enable = 'off';
        drawnow;
        this.saveTimelapse;
        this.currentPos = x.Value;
        this.currentTrap = 1;
        this.currentCell = this.cellLabels{1}(find(this.cellLabels{1},1));
        this.refreshTrapIms;
        this.refreshPosImage;
        updateTrapList;
        updateCellList;
        updateCell;
        updateImageCacheSize;
        this.Enable = 'on';
    end

    function updateTrap(x,~)
        this.currentTrap = x.Value;
        this.refreshTrapIms;
        this.refreshPosImage;
        updateCellList;
        updateCell;
    end

    function updateCell(~,~)
        this.refreshColours;
        this.refreshSegIms;
        this.refreshPosImSeg(this.currentTrap);
        this.refreshTrackPlot;
    end

    function updateImageCacheSize()
        this.ctrls.cachepanel.Title = sprintf('Image cache (%u MB)',...
            round(this.imcache.cachesize/8/1000/1000));
    end

    function clearcache(~,~)
        this.imcache.clearcache;
        this.refreshTrapIms;
        updateImageCacheSize;
    end

    function updateTrapList(~,~)
        this.updateTrapList;
    end

    function updateCellList(~,~)
        clab = this.cellLabels{this.currentTrap};
        cellText = arrayfun(@(x) sprintf('Cell %u',x),clab,'Uni',0);
        val = min(max([1,this.ctrls.cells.Value]),length(clab));
        set(this.ctrls.cells,'String',cellText,'Value',val);
    end

    function resetProgressBar()
        % Reset to top level, which also automatically removes the ticker
        % when the bar count reaches 0:
        while this.progressBar.bar_count > 0
            this.progressBar.pop_bar;
        end
        
        % Just in case the user closed the window manually in
        % the previous run:
        assignin('base', ['prog_terminate',...
            num2str(this.progressBar.window_number)], false);
    end

    function copyoutline(tpoff)
        trap = this.currentTrap;
        ctp = this.currentTimepoint;
        tmpltp = ctp+tpoff;
        if tmpltp<1 || tmpltp>this.ntimepoints
            warning('time point out of range...');
        end
        
        if isempty(this.cellLabels{trap})
            warning('no cells to copy...');
            return
        elseif isempty(this.currentCell)
            warning('a cell needs to be selected...');
            return
        else
            cellLabel = this.cellLabels{trap}(this.currentCell);
        end
        
        % Translate time points to cTimelapse indices
        ctp = this.tpInds(ctp);
        tmpltp = this.tpInds(tmpltp);
        
        % Get template trap and cell info
        trapInfo = this.cTimelapse.cTimepoint(tmpltp).trapInfo(trap);
        if ~trapInfo.cellsPresent
            warning('no outlines to copy...');
            return
        end
        tci = find(trapInfo.cellLabel==cellLabel);
        if length(tci)~=1
            warning('no outline to copy for the active cell...');
            return
        end
        
        needsTrackUpdate = ~ismember(cellLabel,...
            this.cTimelapse.cTimepoint(ctp).trapInfo(trap).cellLabel);
        
        % Use template to set outline for the current cell at ctp:
        ci = this.ensureCellIndex(trap,ctp,cellLabel);
        % this.setContourParams(centre,radii,angles,trap,tp,ci)
        this.setContourParams(trapInfo.cell(tci).cellCenter,...
            trapInfo.cell(tci).cellRadii,...
            trapInfo.cell(tci).cellAngle,trap,ctp,ci);
        
        % Refresh relevant properties
        if needsTrackUpdate
            this.refreshCellLabels;
            this.refreshTracks;
        end
        updateTrapList;
        updateCellList;
        updateCell;
    end

    function resetoutline(~,~)
        trap = this.currentTrap;
        ctp = this.tpInds(this.currentTimepoint);
        trapInfo = this.cTimelapse.cTimepoint(ctp).trapInfo(trap);
        if isempty(this.cellLabels{trap})
            warning('no cells to reset...');
            return
        elseif isempty(this.currentCell)
            warning('a cell needs to be selected...');
            return
        else
            cellLabel = this.cellLabels{trap}(this.currentCell);
        end
        
        if ~(isfield(trapInfo,'cellOriginal') && isfield(trapInfo,'cellLabelOriginal'))
            warning('no original segmentation outline available...');
            return
        end
        
        tci = find(trapInfo.cellLabelOriginal==cellLabel);
        if length(tci)~=1
            warning('no original segmentation outline available for the active cell...');
            return
        end
        
        needsTrackUpdate = ~ismember(cellLabel,trapInfo.cellLabel);
        
        % Use original template to reset outline for current cell:
        ci = this.ensureCellIndex(trap,ctp,cellLabel);
        % this.setContourParams(centre,radii,angles,trap,tp,ci)
        this.setContourParams(trapInfo.cellOriginal(tci).cellCenter,...
            trapInfo.cellOriginal(tci).cellRadii,...
            trapInfo.cellOriginal(tci).cellAngle,trap,ctp,ci);
        
        % Refresh relevant properties
        if needsTrackUpdate
            this.refreshCellLabels;
            this.refreshTracks;
        end
        updateTrapList;
        updateCellList;
        updateCell;
    end

    function setnknots(~,~)
        trap = this.currentTrap;
        tp = this.tpInds(this.currentTimepoint+this.dragtp);
        trapInfo = this.cTimelapse.cTimepoint(tp).trapInfo(trap);
        if isempty(this.cellLabels{trap})
            warning('no cells to reset...');
            return
        elseif isempty(this.currentCell)
            warning('a cell needs to be selected...');
            return
        else
            cellLabel = this.cellLabels{trap}(this.currentCell);
        end
        
        ci = find(trapInfo.cellLabel==cellLabel);
        if length(ci)~=1
            warning('active cell has no outline at drag time point...');
            return
        end
        
        % Use current mask to compensate for shape/orientation
        mask = imfill(full(trapInfo.cell(ci).segmented),'holes');
        rprops = regionprops(mask,'MajorAxisLength','MinorAxisLength','Orientation');
        axmaj = rprops.MajorAxisLength; axmin = rprops.MinorAxisLength;
        ori = pi*rprops.Orientation/180;
        
        % Determine number of knots and new angles
        if this.ctrls.nknots.Value==1
            warning('not implemented...');
            return
        end
        nknots = this.nknot_values(this.ctrls.nknots.Value - 1);
        astep = 2*pi/nknots;
        angles = (0:nknots-1)*astep;
        
        % Roughly compensate for elliptical squeezing by converting
        % parameterised ellipse parameter to true angle. See:
        % https://math.stackexchange.com/a/436125
        angles = angles - atan((axmaj-axmin)*tan(angles)./(axmaj+axmin+tan(angles).^2));
        % Reorient
        angles = mod(-ori+2*pi,2*pi)+angles;
        
        % Determine radii of new knots based on current outline
        if this.ctrls.cartesianspline.Value
            [px,py] = BABYutil.eval_cartesian_spline_from_radii(angles,...
                double(trapInfo.cell(ci).cellRadii),double(trapInfo.cell(ci).cellAngle));
            radii = sqrt(px.^2+py.^2);
        else
            radii = BABYutil.get_radii_from_radii(angles,...
                double(trapInfo.cell(ci).cellRadii),double(trapInfo.cell(ci).cellAngle));
        end
        
        % Ensure any original outlines are saved
        this.saveOriginalOutlines(trap,tp);
        
        this.setContourParams(trapInfo.cell(ci).cellCenter,...
            radii,angles,trap,tp,ci);
        
        % Refresh relevant properties
        updateTrapList;
        updateCellList;
        updateCell;
    end

    function jumpto(trap,tpi,cellLabel)
        this.currentTrap = trap;
        this.ctrls.traps.Value = trap;
        this.dragtp = 0;
        this.currentTimepoint = tpi;
        updateCellList;
        
        if ~isempty(this.ctrls.cells.String)
            if isempty(cellLabel)
                this.currentCell = 1;
            else
                this.currentCell = find(...
                    this.cellLabels{this.currentTrap}==cellLabel,1);
            end
        end
        
        this.refreshTrapIms;
        this.refreshPosImage;
        this.refreshColours;
        setSlider; slider_update;
        this.refreshTrackPlot;
    end

    function jumpToLastModified(~,~)
        % Find last modified trap and time point for this position:
        maxmodifiedon = -Inf;
        latest_tpi = [];
        latest_trap = [];
        for tpi=1:this.ntimepoints
            tp = this.tpInds(tpi);
            trapInfo = this.cTimelapse.cTimepoint(tp).trapInfo;
            if isfield(trapInfo,'modifiedon')
                for trap=1:length(trapInfo)
                    if ~isempty(trapInfo(trap).modifiedon) && ...
                            trapInfo(trap).modifiedon > maxmodifiedon
                        maxmodifiedon = trapInfo(trap).modifiedon;
                        latest_tpi = tpi;
                        latest_trap = trap;
                    end
                end
            end
        end
        
        if isempty(latest_tpi)
            warning('no cells have been modified...');
            return
        end
        
        trapInfo = this.cTimelapse.cTimepoint(this.tpInds(latest_tpi)).trapInfo(latest_trap);
        latest_cell = [];
        if trapInfo.cellsPresent && isfield(trapInfo.cell,'modifiedon')
            for c=1:length(trapInfo.cell)
                if ~isempty(trapInfo.cell(c).modifiedon) && ...
                        trapInfo.cell(c).modifiedon >= maxmodifiedon
                    latest_cell = trapInfo.cellLabel(c);
                end
            end
        end
        
        jumpto(latest_trap,latest_tpi,latest_cell);
    end

    function prevIsCurated(~,~)
        if isempty(this.isCurated)
            warndlg('No traps have been marked as curated...');
            return
        end
        curated = sortrows(double(this.isCurated),[-1,-1]);
        curtrap = this.currentTrap; curtp = this.currentTimepoint;
        prevind = find((curated(:,1)==curtrap & curated(:,2)<curtp) ...
            | curated(:,1)<curtrap,1);
        if ~isempty(prevind)
            jumpto(curated(prevind,1),curated(prevind,2),[]);
        end
    end

    function nextIsCurated(~,~)
        if isempty(this.isCurated)
            warndlg('No traps have been marked as curated...');
            return
        end
        curated = sortrows(double(this.isCurated),[1,1]);
        curtrap = this.currentTrap; curtp = this.currentTimepoint;
        nextind = find((curated(:,1)==curtrap & curated(:,2)>curtp) ...
            | curated(:,1)>curtrap,1);
        if ~isempty(nextind)
            jumpto(curated(nextind,1),curated(nextind,2),[]);
        end
    end
    
    function jumpToNoFocusAssigned(~,~)
        % Find last modified trap and time point for this position:
        for trap=1:this.ntraps
            for tpi=1:this.ntimepoints
                tp = this.tpInds(tpi);
                trapInfo = this.cTimelapse.cTimepoint(tp).trapInfo(trap);
                if ~trapInfo.cellsPresent || isempty(trapInfo.cellLabel)
                    continue
                end
                cellInfo = trapInfo.cell;
                if ~isfield(cellInfo,'focusStack')
                    jumpto(trap,tpi,trapInfo.cellLabel(1));
                    return
                elseif any(cellfun(@isempty,{cellInfo.focusStack}))
                    missing = cellfun(@isempty,{cellInfo.focusStack});
                    jumpto(trap,tpi,trapInfo.cellLabel(find(missing,1)));
                    return
                end
            end
        end
        
        helpdlg('All cells have an assigned focus...');
    end

    function autoNormalisation(~,~)
        c = 1;
        if length(this.currentChannels)>1
            [c,ok] = listdlg('ListString',this.currentChannels,...
                'SelectionMode','single','Name','Pick channel to normalise',...
                'ListSize',[200,150]);
            if ~ok, return; end
        end
        channel = this.currentChannels{c};
        normIm = this.imcache.getTrapTimepoints(this.currentPos,...
            this.currentTrap,channel,this.currentTimepoint);
        normquants = [0.0001,0.9999];
        imrange = quantile(normIm(:),normquants);
        channel = regexprep(channel,'\W','_');
        this.channelNormalisation.(channel) = imrange(1) + ...
            ([0,1]-normquants(1))*diff(imrange)/diff(normquants);
        this.refreshTrapIms;
        this.refreshSegIms;
    end

    function adjustNormalisation(~,~)
        this.adjustNormGUI;
    end

    function saveStack(~,~)
        channel = this.ctrls.saveChannel.String{this.ctrls.saveChannel.Value};
        this.cTimelapse.exportStacks(this.saveStackDir.fileName,...
            channel,this.currentTimepoint,this.currentTrap);
    end

    function selectAllLineages(~,~)
        this.cTimelapse.autoSelectLineages(this.currentTrap);
        this.haschanged = true;
        updateTrapList;
        updateCellList;
        updateCell;
    end

    function unselectAllCells(~,~)
        this.cTimelapse.cellsToPlot(this.currentTrap,:) = 0;
        this.haschanged = true;
        updateTrapList;
        updateCellList;
        updateCell;
    end
end
