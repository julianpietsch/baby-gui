function UpdateWindowElements(this,adapt)
%UpdateWindowElements Create/update figure window, panels and axes

if nargin<2 || isempty(adapt)
    adapt = false;
end
    
%% Retrieve default dimensions
dims = this.dimensions;
sp = dims.sp;
impad = dims.impad;
nchannels = length(this.currentChannels);

% Update track display size according to input box
if isfield(this.ctrls,'trackdispsz') && ~isempty(this.ctrls.trackdispsz)
    trackdispsizes = fieldnames(this.trackSizeMap);
    tdsKey = trackdispsizes{this.ctrls.trackdispsz.Value};
    dims.trackHeight = this.trackSizeMap(2).(tdsKey);
end

%% Calculate figure and component dimensions

% Start by assuming that the tileWidth and tileHeight match the source
this.tileWidth = this.srcWidth;
this.tileHeight = this.srcHeight;

% The image width and height depend on how many tiles of context there are
% and on how many channels have been selected for display:
    function [imW,imH] = getImDims(cntxt)
        imW = impad + (this.tileWidth+impad)*(2*cntxt+1);
        imH = impad + (this.tileHeight+impad)*nchannels;
    end

[imwidth,imheight] = getImDims(this.tilesContext);
imzoomW = round(this.imzoom*imwidth);
imzoomH = round(this.imzoom*imheight);

scrsz = get(0,'ScreenSize');
scrHpad = dims.screenHeightPad;
if ispc, scrHpad = 1.5*scrHpad; end
scrW = round(0.96*scrsz(3)); scrH = scrsz(4)-scrHpad;

% Notes:
% - In all cases the ctrls have a fixed size
% - Flexible spacers (*) are used to centre the pos and trap slider when
%   they are smaller than the minimum width
% - Reduce the default zoom and/or context to fit within screen size:
    function adaptZoom(padW,padH)
        W = this.imzoom*imwidth + padW;
        H = this.imzoom*imheight + padH;
        if ~adapt || (W<=scrW && H<=scrH), return; end
        Zs = this.zoomlevels(this.zoomlevels<this.imzoom);
        widths = Zs*imwidth + padW;
        heights = Zs*imwidth + padH;
        Zs = Zs(widths<=scrW & heights<=scrH);
        if isempty(Zs)
            this.imzoom_val = min(this.zoomlevels);
        else
            this.imzoom_val = max(Zs);
        end
        imzoomW = round(this.imzoom*imwidth);
        imzoomH = round(this.imzoom*imheight);
        if isfield(this.ctrls,'zoom')
            this.ctrls.zoom.Value = find(this.imzoom==this.zoomlevels,1);
        end
    end

    function adaptContext(padW)
        Z = this.imzoom;
        if ~adapt || Z*imwidth+padW<=scrW, return; end
        for cntxt=this.tilesContext:-1:0
            [W,~] = getImDims(cntxt);
            if Z*W+padW<=scrW, break; end
        end
        [imwidth,~] = getImDims(cntxt);
        imzoomW = round(Z*imwidth);
        this.tilesContext = cntxt;
        if isfield(this.ctrls,'tpcontext')
            this.ctrls.tpcontext.Value = find(this.tilesContext==0:4,1);
        end
    end

    function adaptTileDims(padW,padH)
        W = this.imzoom*imwidth + padW;
        H = this.imzoom*imheight + padH;
        if W>scrW
            Nctxt = 2*this.tilesContext+1;
            this.tileWidth = max(floor(((scrW-padW-impad)/Nctxt-impad)/this.imzoom),10);
        end
        if H>scrH
            this.tileHeight = max(floor(((scrH-padH-impad)/nchannels-impad)/this.imzoom),10);
        end
        [imwidth,imheight] = getImDims(this.tilesContext);
        imzoomW = round(this.imzoom*imwidth);
        imzoomH = round(this.imzoom*imheight);
    end

ctrlWidth = dims.ctrlWidth;
ctrlHeight = dims.ctrlHeight;
posAspect = this.cTimelapse.imSize(2)/this.cTimelapse.imSize(1); % height/width

switch this.layout
    case 'vertical'
        % ---------------------
        % |W1  trap slider  W1|
        % ---------------------
        % |         |         |
        % |  ctrls  |W2 pos   |
        % |         |         |
        % ---------------------

        % Fix the height of the pos overview; adjust spacer for aspect and
        % width of trap slider
        posHeight = ctrlHeight;
        posWidth = round(posHeight*posAspect);
        
        % The fixed width depends on the control and pos panels:
        fixedwidth = 3*sp + ctrlWidth + posWidth;
        % The fixed height depends on the control panel, cell track axis 
        % and two slider heights:
        fixedheight = 6*sp + ctrlHeight + dims.trackHeight ...
            + dims.xaxis + 2*dims.slideHeight;
        
        % If necessary, adapt zoom then context
        adaptZoom(2*sp,fixedheight);
        adaptContext(2*sp);
        adaptTileDims(2*sp,fixedheight);
        
        % Determine figure width and height and calculate required spacers:
        imPanelW = max(imzoomW+sp+dims.slideHeight,fixedwidth-2*sp);
        figwidth = 2*sp + imPanelW;
        Wspacer1 = (figwidth-imzoomW-sp-dims.slideHeight)/2;
        Wspacer2 = figwidth - fixedwidth;
        figheight = fixedheight + imzoomH;
    
    case 'horizontal'
        % ---------------------------
        % |       H1      |         |
        % |               |  ctrls  |
        % |               |         |
        % |W trap slider W|---------|
        % |               |    H2   |
        % |               |W2 pos   |
        % |       H1      |         |
        % ---------------------------
        
        % Fix the width of the pos overview; adjust spacer for aspect and
        % height of trap slider
        posWidth = ctrlWidth;
        posHeight = round(posWidth/posAspect);
        Wspacer2 = 0;
        if (posHeight+ctrlHeight)>scrH
            posHeight = scrH-ctrlHeight;
            posWidth = round(posHeight*posAspect);
            Wspacer2 = (ctrlWidth-posWidth)/2;
        end
        
        % The fixed width depends on the control panel and Y slider:
        fixedwidth = 4*sp + ctrlWidth + dims.slideHeight;
        % The fixed height of the left column depends on the cell track 
        % axis and slider heights:
        fixedheightleft = 5*sp + dims.trackHeight + dims.xaxis ...
            + 2*dims.slideHeight;
        % The fixed height of the right column depends on the control and
        % pos panels:
        fixedheightright = 3*sp + ctrlHeight + posHeight;
        
        % If necessary, adapt context then zoom
        adaptContext(fixedwidth);
        adaptZoom(fixedwidth,fixedheightleft);
        adaptTileDims(fixedwidth,fixedheightleft);
        
        % Determine figure width and height and calculate required spacers:
        imPanelW = max(imzoomW+sp+dims.slideHeight,dims.minWidth);
        figwidth = fixedwidth-sp-dims.slideHeight+imPanelW;
        Wspacer = (figwidth-imzoomW-ctrlWidth-dims.slideHeight-2*sp)/2;
        figheight = max(fixedheightleft+imzoomH,fixedheightright);
        Hspacer1 = (figheight - imzoomH - fixedheightleft)/2;
        Hspacer2 = figheight - fixedheightright;
        
    case 'overview'
        % -----------------------------------
        % |        H        |               |
        % |W1 trap slider W1|               |
        % |        H        |               |
        % |-----------------|      pos      |
        % |                 |               |
        % | ctrls         W2|               |
        % |                 |               |
        % -----------------------------------
        
        % Maximise height of pos overview and adjust its width for aspect
        posHeight = scrH;
        posWidth = round(posHeight*posAspect);
        
        % The fixed width depends on the pos panel and Y slider:
        fixedwidth = 4*sp+posWidth+dims.slideHeight;
        % The fixed height of the left column depends on the cell track 
        % axis, sliders and ctrl panel:
        fixedheight = 5*sp + dims.trackHeight + dims.xaxis ...
            + 2*dims.slideHeight + ctrlHeight;
        
        % If necessary, adapt context then zoom
        adaptContext(fixedwidth);
        adaptZoom(fixedwidth,fixedheight);
        adaptTileDims(fixedwidth,fixedheight);
        
        % Recalculate posHeight in case it needs to be larger
        posHeight = max(posHeight,fixedheight+imzoomH-2*sp);
        posWidth = round(posHeight*posAspect);
        
        % Determine figure width and height and calculate required spacers:
        imPanelW = max([imzoomW+sp+dims.slideHeight,dims.minWidth,ctrlWidth]);
        figwidth = 3*sp + posWidth + imPanelW;
        Wspacer1 = (figwidth-imzoomW-posWidth-2*sp-dims.slideHeight)/2;
        Wspacer2 = figwidth - ctrlWidth - posWidth - 2*sp;
        figheight = posHeight + 2*sp;
        Hspacer = (figheight - imzoomH - fixedheight)/2;
        
    otherwise
        error('The "%s" layout has undefined dimension',this.layout);
end

%% Draw/update window elements

% Draw the figure window if it hasn't already been initialised
if isempty(this.frame.fig) || ~isvalid(this.frame.fig)
    this.frame.fig = figure('MenuBar','none','Name','',...
        'Resize','off','Visible','on');
    this.frame.fig.UserData = this;
end
% In any case, update the dimensions and position of the figure
set(this.frame.fig,'Position',...
    [(scrsz(3)-figwidth)/2,(scrsz(4)-figheight)/2,figwidth,figheight]);

% Draw the control panel and position overview axis:
if isempty(this.frame.ctrl) || ~isvalid(this.frame.ctrl)
    this.frame.ctrl = uipanel('Parent',this.frame.fig);
end
if isempty(this.frame.pos) || ~isvalid(this.frame.pos)
    this.frame.pos = axes('Parent',this.frame.fig,'XTick',[],'YTick',[],...
        'XLim',[0.5,this.cTimelapse.imSize(2)+0.5],...
        'YLim',[0.5,this.cTimelapse.imSize(1)+0.5]);
    % Show white image to begin with
    image('CData',ones([this.cTimelapse.imSize,3]),'Parent',this.frame.pos);
end

% Draw the time slider control
if isempty(this.frame.slider) || ~isvalid(this.frame.slider)
    this.frame.slider = uicontrol('Style','slider','Parent',this.frame.fig,...
        'Enable',this.Enable);
end

% Draw the x offset slider control
enable = this.Enable;
if this.tileWidth==this.srcWidth, enable = 'off'; end
if isempty(this.frame.xslider) || ~isvalid(this.frame.xslider)
    this.frame.xslider = uicontrol('Style','slider','Parent',this.frame.fig,...
        'Enable',enable);
else
    set(this.frame.xslider,'Enable',enable);
end

% Draw the y offset slider control
enable = this.Enable;
if this.tileHeight==this.srcHeight, enable = 'off'; end
if isempty(this.frame.yslider) || ~isvalid(this.frame.yslider)
    this.frame.yslider = uicontrol('Style','slider','Parent',this.frame.fig,...
        'Enable',enable);
else
    set(this.frame.yslider,'Enable',enable);
end

% Draw the cell tracking axis (includes an X axis)
if isempty(this.frame.track) || ~isvalid(this.frame.track)
    this.frame.track = axes('Parent',this.frame.fig,'YTick',[]);
end

% Draw the image axis
if isempty(this.frame.imAx) || ~isvalid(this.frame.imAx)
    this.frame.imAx = axes('Parent',this.frame.fig,'XTick',[],'YTick',[]);
    % Show white image to begin with
    this.frame.im = image('CData',ones([imheight,imwidth,3]),'Parent',this.frame.imAx);
end

%% Position elements based on layout mode

hcum = sp; % Used to track the cumulative height of added plots

switch this.layout
    case 'vertical'
        % ---------------------
        % |W1  trap slider  W1|
        % ---------------------
        % |         |         |
        % |  ctrls  |W2 pos   |
        % |         |         |
        % ---------------------
        
        % Position the control panel and position overview axis:
        set(this.frame.ctrl,'Units','pixels','Position',...
            [sp,hcum,ctrlWidth,ctrlHeight]);
        set(this.frame.pos,'Units','pixels','Position',...
            [sp+ctrlWidth+Wspacer2,hcum,posWidth,posHeight]);
        hcum = hcum + ctrlHeight + sp;
        
        % Position the time slider control
        set(this.frame.slider,'Position',...
            [sp,hcum,figwidth-sp,dims.slideHeight]);
        hcum = hcum + dims.slideHeight + sp;
        
        % Position the cell tracking axis (includes an X axis)
        set(this.frame.track,'Units','pixels','Position',...
            [sp,hcum+dims.xaxis,figwidth-sp,dims.trackHeight]);
        hcum = hcum + dims.xaxis + dims.trackHeight + sp;
        
        % Position the x slider control
        set(this.frame.xslider,'Position',...
            [Wspacer1,hcum,imzoomW,dims.slideHeight]);
        hcum = hcum + dims.slideHeight + sp;
        
        % Position the image axis
        set(this.frame.imAx,'Units','pixels',...
            'XLim',[0.5,imwidth+0.5],'YLim',[0.5,imheight+0.5],...
            'Position',[Wspacer1,hcum,imzoomW,imzoomH]);
        
        % Position the y slider control
        set(this.frame.yslider,'Position',...
            [Wspacer1+imzoomW+sp,hcum,dims.slideHeight,imzoomH]);
        
    case 'horizontal'
        % ---------------------------
        % |       H1      |         |
        % |               |  ctrls  |
        % |               |         |
        % |W trap slider W|---------|
        % |               |    H2   |
        % |               |W2 pos   |
        % |       H1      |         |
        % ---------------------------
        
        % Position the time slider control
        set(this.frame.slider,'Position',...
            [sp,hcum,imPanelW,dims.slideHeight]);
        hcum = hcum + dims.slideHeight + sp;
        
        % Position the cell tracking axis (includes an X axis)
        set(this.frame.track,'Units','pixels','Position',...
            [sp,hcum+dims.xaxis,imPanelW,dims.trackHeight]);
        hcum = hcum + dims.xaxis + dims.trackHeight + sp;
        
        % Position the x slider control
        set(this.frame.xslider,'Position',...
            [Wspacer,hcum+Hspacer1,imzoomW,dims.slideHeight]);
        hcum = hcum + Hspacer1 + dims.slideHeight + sp;
        
        % Position the image axis
        set(this.frame.imAx,'Units','pixels',...
            'XLim',[0.5,imwidth+0.5],'YLim',[0.5,imheight+0.5],...
            'Position',[Wspacer,hcum,imzoomW,imzoomH]);
        
        % Position the y slider control
        set(this.frame.yslider,'Position',...
            [Wspacer+imzoomW+sp,hcum,dims.slideHeight,imzoomH]);
        
        % Position the control panel and position overview axis:
        hcum = sp;
        set(this.frame.pos,'Units','pixels','Position',...
            [2*sp+imPanelW+Wspacer2,hcum,posWidth,posHeight]);
        hcum = hcum + posHeight + Hspacer2 + sp;
        
        set(this.frame.ctrl,'Units','pixels','Position',...
            [2*sp+imPanelW,hcum,ctrlWidth,ctrlHeight]);
        
    case 'overview'
        % -----------------------------------
        % |        H        |               |
        % |W1 trap slider W1|               |
        % |        H        |               |
        % |-----------------|      pos      |
        % |                 |               |
        % | ctrls         W2|               |
        % |                 |               |
        % -----------------------------------
        
        % Position the control panel
        set(this.frame.ctrl,'Units','pixels','Position',...
            [sp,hcum,ctrlWidth,ctrlHeight]);
        hcum = hcum + ctrlHeight + sp;
        
        % Position the time slider control
        set(this.frame.slider,'Position',...
            [sp,hcum,imPanelW,dims.slideHeight]);
        hcum = hcum + dims.slideHeight + sp;
        
        % Position the cell tracking axis (includes an X axis)
        set(this.frame.track,'Units','pixels','Position',...
            [sp,hcum+dims.xaxis,imPanelW,dims.trackHeight]);
        hcum = hcum + dims.xaxis + dims.trackHeight + sp;
        
        % Position the x slider control
        set(this.frame.xslider,'Position',...
            [Wspacer1,hcum+Hspacer,imzoomW,dims.slideHeight]);
        hcum = hcum + Hspacer + dims.slideHeight + sp;
        
        % Position the image axis
        set(this.frame.imAx,'Units','pixels',...
            'XLim',[0.5,imwidth+0.5],'YLim',[0.5,imheight+0.5],...
            'Position',[Wspacer1,hcum,imzoomW,imzoomH]);
        
        % Position the y slider control
        set(this.frame.yslider,'Position',...
            [Wspacer1+imzoomW+sp,hcum,dims.slideHeight,imzoomH]);
        
        % Position the position overview axis:
        hcum = sp;
        set(this.frame.pos,'Units','pixels','Position',...
            [2*sp+imPanelW,hcum,posWidth,posHeight]);
        
        % By default, now scroll the position overview image
        if isfield(this.ctrls,'scrolloverview')
            this.ctrls.scrolloverview.Value = true;
        end
        
    otherwise
        error('The "%s" layout is undefined',this.layout);
end

this.refreshFocusAnnotStatus;

% For the background image, use white with light green to highlight current
% timepoint. NB: image axes are flipped relative to typical indexing.
this.bgdIm = ones([imheight,imwidth,3]);
highlightXloc = (dims.impad+this.tileWidth)*this.tilesContext;
highlightSlice = highlightXloc+(1:this.tileWidth+2*dims.impad);
this.bgdIm(:,highlightSlice,[1,3]) = 0.85;

this.windowUpdateRequired = false;

set(this.frame.fig,'Visible','on');
drawnow;

end