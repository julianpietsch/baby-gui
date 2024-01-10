function initialiseControls(this)

%% Initial state of the GUI
posImgFiles = {}; % Start with no image files
posNames = {}; % Cell array specifying pos names for each image file
imgFile = '';
reader = []; % empty when no file loaded
ctrls = struct(); % stores handles to controls if needed

%% Set up the GUI

lh = 15; % Label height in pixels
sh = 20; % Selection menu height
ch = 30; % Control height in pixels
sp = 3; % Spacing in pixels
slideHeight = 20; % Height (or width) for horizontal (or vertical) slider
screenHeightPad = 70; % Padding for OS menus

scrsz = get(0,'ScreenSize');

imgW = 512; imgH = 512; % default image width and height in pixels
rawImg = ones(imgH,imgW);

ctrlW = 300; % control panel width in pixels
minCtrlH = 500; % min control panel height in pixels
ctrlH = minCtrlH; % control panel height in pixels

maxImgW = scrsz(3)-4*sp-ctrlW-slideHeight;
maxImgH = scrsz(4)-3*sp-slideHeight-screenHeightPad;

% Create a dialog for the controls
fig = dialog('Name','Create new cExperiment...',...
    'WindowStyle','normal','Visible','on');

%% Initialise overall layout:

% We want a control panel, image display axes, sliders for xy scroll, and a
% slider for time point scroll
ctrlpanel = uipanel('Parent',fig);
imgax = axes('Parent',fig,'XTick',[],'YTick',[],...
        'XLim',[0.5,imgH+0.5],...
        'YLim',[0.5,imgW+0.5],'YDir','reverse');
image('CData',rawImg(:,:,ones(3,1)),'Parent',imgax);
ctrls.xslider = uicontrol('Style','slider','Parent',fig,'Enable','off');
ctrls.yslider = uicontrol('Style','slider','Parent',fig,'Enable','off');
ctrls.tslider = uicontrol('Style','slider','Parent',fig,'Enable','off');
ctrls.xylabel = uicontrol('Parent',fig,'Style','text',...
            'HorizontalAlignment','center','String','xy');
ctrls.tlabel = uicontrol('Parent',fig,'Style','text',...
            'HorizontalAlignment','center','String','t');
        
% Then use callback to position everything
updateFrame;

%% Slider control
addlistener(ctrls.tslider,'Value','PostSet',@updateTslider);
addlistener(ctrls.xslider,'Value','PostSet',@updateXYslider);
addlistener(ctrls.yslider,'Value','PostSet',@updateXYslider);

%% Set up tabs
tabpanel = uix.TabPanel('Parent',ctrlpanel,'Padding',3);
filebox = uix.VBox('Parent',tabpanel,'Padding',0,'Spacing',sp);
trapbox = uix.VBox('Parent',tabpanel,'Padding',0,'Spacing',sp);
metabox = uix.VBox('Parent',tabpanel,'Padding',0,'Spacing',sp);
cfgbox = uix.VBox('Parent',tabpanel,'Padding',0,'Spacing',sp);
procbox = uix.VBox('Parent',tabpanel,'Padding',0,'Spacing',sp);
set(tabpanel,'TabTitles',{'Files','Traps','Metadata','Settings','Process'},...
    ...% 'TabEnables',[true,false,false,false,false]);
    'TabEnables',{'on','on','on','on','on'});

%% General control panel

ctrlboxHghts = [];

% Save folder
uiLabel(filebox,'Folder to save results to:'); ctrlboxHghts(end+1) = lh;
ctrls.savedir = BABYutil.uiFile(filebox,'/',pwd,pwd); ctrlboxHghts(end+1) = ch;

% Image files
hbx = uix.HBox('Parent',filebox,'Padding',0,'Spacing',sp);
ctrlboxHghts(end+1) = ch;
ctrls.imgfileadd = uicontrol('Parent',hbx,'Style','pushbutton',...
    'String','Add file...','Callback',@addFile);
ctrls.imgfileremove = uicontrol('Parent',hbx,'Style','pushbutton',...
    'String','Remove file','Callback',@removeFile,'Enable','off');

ctrls.imgfilelist = uicontrol('Parent',filebox,'Style','listbox',...
    'String',{},'Value',[]);
ctrlboxHghts(end+1) = -1;

% Stacking
stackingbox = uix.HBox('Parent',filebox,'Padding',0,'Spacing',sp);
ctrls.stackingon = uicontrol('Parent',stackingbox,'Style','checkbox',...
    'String','Stacking axis:','Value',false,'Callback',@toggleStacking);
ctrls.stackingdim = uicontrol('Parent',stackingbox,'Style','popupmenu',...
    'String',{'T','Z','C'},'Enable','off','Callback',@updateReaders);
pps_values = 0;
ctrls.stackingpps = uicontrol('Parent',stackingbox,'Style','popupmenu',...
    'String',{'Position-wise'},'Value',1,'Enable','off');
set(stackingbox,'Widths',[-3,-1,-3]);
ctrlboxHghts(end+1) = 1.1*sh;

% Flip
ctrls.flip = uicontrol('Parent',filebox,'Style','checkbox',...
    'String','Flip image?','Value',false,'Callback',@updateImage);
ctrlboxHghts(end+1) = sh;

% Rotation
uiLabel(filebox,'Rotation:'); ctrlboxHghts(end+1) = lh;
rotLbls = {'None','90','180','270'};
rotVals = [0,1,2,3];
ctrls.rot = uicontrol('Parent',filebox,'Style','popupmenu',...
    'String',rotLbls,'Value',1,'Callback',@updateImage);
ctrlboxHghts(end+1) = sh;

% Position
uiLabel(filebox,'Position:'); ctrlboxHghts(end+1) = lh;
ctrls.poses = uicontrol('Parent',filebox,'Style','popupmenu',...
    'String',{'None'},'Value',1,'Enable','off','Callback',@updatePosition);
ctrlboxHghts(end+1) = sh;

% Channel
uiLabel(filebox,'Channel:'); ctrlboxHghts(end+1) = lh;
ctrls.channel = uicontrol('Parent',filebox,'Style','popupmenu',...
    'String',{'None'},'Value',1,'Enable','off','Callback',@updateRawImage);
ctrlboxHghts(end+1) = sh;

% Z section
uiLabel(filebox,'Z section:'); ctrlboxHghts(end+1) = lh;
ctrls.zsect = uicontrol('Parent',filebox,'Style','popupmenu',...
    'String',{'None'},'Value',1,'Enable','off','Callback',@updateRawImage);
ctrlboxHghts(end+1) = sh;

% Time point
uiLabel(filebox,'Time point:'); ctrlboxHghts(end+1) = lh;
ctrls.tp = uicontrol('Parent',filebox,'Style','popupmenu',...
    'String',{'None'},'Value',1,'Enable','off','Callback',@updateTimePoint);
ctrlboxHghts(end+1) = sh;

% Zoom
uiLabel(filebox,'Zoom:'); ctrlboxHghts(end+1) = lh;
zoomLbls = {'1/4x','1/3x','1/2x','1x','2x','3x','4x'};
zoomMult = cellfun(@(x) eval(x(1:end-1)),zoomLbls);
zoom1 = find(zoomMult==1,1);
ctrls.zoom = uicontrol('Parent',filebox,'Style','popupmenu',...
    'String',zoomLbls,'Value',zoom1,'Enable','off','Callback',@updateImage);
ctrlboxHghts(end+1) = sh;

set(filebox,'Heights',ctrlboxHghts);

%% Trap control panel

trapboxHghts = [];

% Use traps?
ctrls.usetraps = uicontrol('Parent',trapbox,'Style','checkbox',...
    'String','Use traps?','Value',false,'Callback',@toggleTraps);
trapboxHghts(end+1) = sh;

trapchannel = 1;

% Trap area selection
drawbox = uix.HBox('Parent',trapbox,'Padding',0,'Spacing',sp);
trapboxHghts(end+1) = ch;
ctrls.drawTrap = uicontrol('Parent',drawbox,'Style','pushbutton',...
    'String','Draw trap...','Callback',@startTrapDraw,'Enable','off');
ctrls.finishDrawTrap = uicontrol('Parent',drawbox,'Style','pushbutton',...
    'String','...done','Callback',@endTrapDraw,'Enable','off');

trapimpanel = uipanel('Parent',trapbox);
trapboxHghts(end+1) = -1;
trapImg = ones(80,80);
trapimgax = axes('Parent',trapimpanel,'XTick',[],'YTick',[],...
    'XLim',[0.5,size(trapImg,1)+0.5],'YLim',[0.5,size(trapImg,2)+0.5],...
    'YDir','reverse');
image('CData',trapImg(:,:,ones(3,1)),'Parent',trapimgax);

updateTrapImg;

set(trapbox,'Heights',trapboxHghts);

%% Meta control panel

metaboxheights = [];

% Date
hbx = uix.HBox('Parent',metabox,'Padding',0,'Spacing',sp);
uiLabel(hbx,'Date');
ctrls.date = uicontrol('Parent',hbx,'Style','edit');
metaboxheights(end+1) = sh;

% Times/time interval
hbx = uix.HBox('Parent',metabox,'Padding',0,'Spacing',sp);
uiLabel(hbx,'Time interval');
ctrls.date = uicontrol('Parent',hbx,'Style','edit');
metaboxheights(end+1) = sh;

% Pixel size
hbx = uix.HBox('Parent',metabox,'Padding',0,'Spacing',sp);
uiLabel(hbx,'Pixel size');
ctrls.date = uicontrol('Parent',hbx,'Style','edit');
metaboxheights(end+1) = sh;

% Channel names
ctrls.chnames = uitable('Parent',metabox,'Data',cell(0,2),...
    'ColumnName',{'Channel ID','Label'},...
    'ColumnFormat',{'char','char'},...
    'ColumnEditable',[false,true],...
    'RowName',[],'RowStriping','off','CellSelectionCallback',@(~,~) 1);
metaboxheights(end+1) = -1;

% Position names
ctrls.posnames = uitable('Parent',metabox,'Data',cell(0,3),...
    'ColumnName',{'Position','ID','Label'},...
    'ColumnFormat',{'numeric','char','char'},...
    'ColumnEditable',[false,false,true],...
    'RowName',[],'RowStriping','off','CellSelectionCallback',@(~,~) 1);
metaboxheights(end+1) = -1;

% Chamber groups
ctrls.chambergroups = uitable('Parent',metabox,'Data',cell(0,2),...
    'ColumnName',{'Group','Positions'},...
    'ColumnFormat',{'char','char'},...
    'ColumnEditable',[true,false],...
    'RowName',[],'RowStriping','off','CellSelectionCallback',@(~,~) 1);
metaboxheights(end+1) = -1;

hbx = uix.HBox('Parent',metabox,'Padding',0,'Spacing',3);
ctrls.grpaddbtn = uicontrol('Parent',hbx,'Style','pushbutton',...
    'String','Add','Callback',@(~,~) 1);
ctrls.grpremovebtn = uicontrol('Parent',hbx,...
    'Style','pushbutton','String','Remove','Callback',@(~,~) 1);
metaboxheights(end+1) = ch;

set(metabox,'Heights',metaboxheights);

%% Settings control panel

cfgboxheights = [];

% Channel registration (offsets)

% Parameters for image registration (tracking) through time

% Segmentation model to use (and channels)

% Cell tracking method to use (i.e., whether or not to apply mmRetrack)

% Channels to extract/extraction settings

% set(cfgbox,'Heights',cfgboxheights);

%% Processing control panel

procboxheights = [];

uicontrol('Parent',procbox,'Style','pushbutton',...
    'String','Create...','Callback',@checkValidity);
procboxheights(end+1) = 1.5*ch;

uicontrol('Parent',procbox,'Style','pushbutton',...
    'String','Segment...');
procboxheights(end+1) = 1.5*ch;

uicontrol('Parent',procbox,'Style','pushbutton',...
    'String','Curate...');
procboxheights(end+1) = 1.5*ch;

uicontrol('Parent',procbox,'Style','pushbutton',...
    'String','Extract...');
procboxheights(end+1) = 1.5*ch;

% Display a summary?

% Filler box
uix.VBox('Parent',procbox,'Padding',0,'Spacing',0);
procboxheights(end+1) = -1;

set(procbox,'Heights',procboxheights);

%% Run
set(fig,'Visible','on');

uiwait(fig);
if ~isempty(reader), delete(reader); end
if ~isvalid(fig)
    return
end

savedir = ctrls.savedir.fileName;
imgfile = posImgFiles;
dostacking = ctrls.stackingon.Value;
stackingdim = ctrls.stackingdim.String{ctrls.stackingdim.Value};
stackingpps = pps_values(ctrls.stackingpps.Value);
doflip = ctrls.flip.Value;
rot = 90*rotVals(ctrls.rot.Value);
usetraps = ctrls.usetraps.Value;

close(fig);

if dostacking
    this.cExperiment = babyExperimentBioformats(imgfile,savedir,...
        'StackAlong',stackingdim,'PosesPerStack',stackingpps);
else
    this.cExperiment = babyExperimentBioformats(imgfile,savedir);
end
if usetraps
  this.cExperiment.trapTemplates = struct('positiveExamples',trapImg);
end
this.cExperiment.createTimelapsePositions('none','all',[],rot,[],usetraps,doflip);
for p=1:numel(this.cExperiment.dirs)
    cTimelapse = this.cExperiment.loadCurrentTimelapse(p);
    cTimelapse.channelForTrapDetection = ...
        find(strcmp(cTimelapse.channelNames,trapchannel),1);
    this.cExperiment.saveTimelapseExperiment;
end

%% Events
    function addFile(~,~)
        [filenames,pathname] = uigetfile('*.*',...
            'Select image file(s)...','MultiSelect','on');
        if isequal(filenames,0) || isequal(pathname,0)
            return
        end
        filenames = cellfun(@(f) fullfile(pathname,f),...
            cellstr(filenames),'Uniform',false);
        Nadd = numel(filenames);
        newPoses = cell(1,Nadd);
        for f=1:numel(filenames)
            tmpReader = ImageReaderBioformats(filenames{f});
            newPoses{f} = tmpReader.posNames;
            delete(tmpReader);
        end
        posNames(end+1:end+Nadd) = newPoses;
        posImgFiles(end+1:end+Nadd) = filenames;
        
        val = ctrls.imgfilelist.Value;
        if isempty(val) || val > numel(posImgFiles)
            val = 1;
        end
        set(ctrls.imgfilelist,'String',posImgFiles,'Value',val);
        ctrls.imgfileremove.Enable = 'on';
        ctrls.tslider.Enable = 'on';
        updatePositionList;
    end

    function removeFile(~,~)
        val = ctrls.imgfilelist.Value;
        if ~isempty(val)
            posImgFiles(val) = [];
            posNames(val) = [];
            val = max(val - 1,1);
            if isempty(posImgFiles)
                val = [];
                ctrls.imgfileremove.Enable = 'off';
                ctrls.tslider.Enable = 'off';
            end
            set(ctrls.imgfilelist,'String',posImgFiles,'Value',val);
            updatePositionList;
        end
    end

    function toggleStacking(~,~)
        if ctrls.stackingon.Value
            ctrls.stackingdim.Enable = 'on';
            ctrls.stackingpps.Enable = 'on';
        else
            ctrls.stackingdim.Enable = 'off';
            ctrls.stackingpps.Enable = 'off';
        end
        updateReaders;
    end
    
    function updateReaders(~,~)
        nposes_total = sum(cellfun(@numel,posNames));
        if nposes_total > 1
            npos_factors = factor(nposes_total);
            multiples = arrayfun(@(n) prod(nchoosek(npos_factors,n),2),...
                1:numel(npos_factors),'Uniform',false);
            multiples = unique(vertcat(multiples{:}))';
        else
            multiples = [];
        end
        pps_values = [0,multiples];
        ctrls.stackingpps.String = [{'Position-wise'},...
            arrayfun(@(m) sprintf('%u in stack',m),multiples,'Uniform',false)];
    end
    
    function updatePositionList
        prevNames = setdiff(ctrls.poses.String,{'None'});
        prevValue = ctrls.poses.Value;
        if isempty(prevNames), prevValue = []; end
        newNames = [posNames{:}];
        newValue = prevValue;
        setEnable = 'on';
        if isempty(newNames)
            newNames = {'None'};
            newValue = 1;
            setEnable = 'off';
        elseif isempty(prevNames) || newValue > numel(newNames)
            newValue = 1;
        end
        set(ctrls.poses,'String',newNames,'Value',newValue,'Enable',setEnable);
        if strcmp(setEnable,'off')
            if ~isempty(reader), delete(reader); reader = []; end
            updateRawImage;
        else
            updatePosition;
        end
    end
    
    function updateFrame(~,~)
        % Calculate dimensions for frame elements
        ctrlH = max(minCtrlH,imgH+2*sp+2*slideHeight);
        fW = 4*sp+ctrlW+imgW+slideHeight;
        fH = 2*sp+ctrlH;
        fP = [(scrsz(3)-fW)/2,(scrsz(4)-fH)/2,fW,fH];
        Hpad = round((fH-imgH-2*sp-2*slideHeight)/2);
        iP = [3*sp+ctrlW+slideHeight,Hpad+2*slideHeight+2*sp,imgW,imgH];
        
        % Assumes that imgW and imgH have been updated
        set(fig,'Position',fP);
        set(ctrlpanel,'Units','pixels','Position',[sp,sp,ctrlW,ctrlH]);
        set(imgax,'Units','pixels','Position',iP);
        set(ctrls.tslider,'Units','pixels','Position',...
            [iP(1),iP(2)-2*sp-2*slideHeight,imgW,slideHeight]);
        set(ctrls.tlabel,'Units','pixels','Position',...
            [iP(1)-sp-slideHeight,iP(2)-2*sp-2*slideHeight,slideHeight,slideHeight]);
        set(ctrls.xslider,'Units','pixels','Position',...
            [iP(1),iP(2)-sp-slideHeight,imgW,slideHeight]);
        set(ctrls.yslider,'Units','pixels','Position',...
            [iP(1)-sp-slideHeight,iP(2),slideHeight,imgH]);
        set(ctrls.xylabel,'Units','pixels','Position',...
            [iP(1)-sp-slideHeight,iP(2)-sp-slideHeight,slideHeight,slideHeight]);
    end

    function updatePosition(~,~)
        nposes = cellfun(@numel,posNames);
        imglist = repelem(posImgFiles,nposes);
        posindlist = arrayfun(@(n) 1:n,nposes,'Uniform',false);
        posindlist = [posindlist{:}];
        selpos = ctrls.poses.Value;
        if ~strcmp(imglist{selpos},imgFile) && ~isempty(reader)
            delete(reader);
            reader = [];
        end
        if isempty(reader)
            imgFile = imglist{selpos};
            reader = ImageReaderBioformats(imgFile);
        end
        
        % Set reduced zoom if necessary
        zoomVal = zoom1;
        if any(reader.imageSize(1:2)>[maxImgW,maxImgH])
            optScaling = [maxImgW,maxImgH]./reader.imageSize(1:2);
            [~,zoomVal] = find((zoomMult-min(optScaling))<0,1,'last');
            if isempty(zoomVal), zoomVal = 1; end
        end
        set(ctrls.zoom,'Value',zoomVal,'Enable','on');
        
        reader.pos = posindlist(selpos);
        imsz = reader.imageSize;
        
        % Reset channel, Z sect and time point to defaults
        prevChan = 'None';
        if ~isempty(ctrls.channel.String) && ~isempty(ctrls.channel.Value)
            prevChan = ctrls.channel.String{ctrls.channel.Value};
        end
        if sum(strcmp(reader.channels,prevChan))==1
            dfltCh = find(strcmp(reader.channels,prevChan),1);
        else
            dfltCh = find(contains(reader.channels,...
                {'field','phase','dic','tl'},'IgnoreCase',true),1);
            if isempty(dfltCh) && ~isempty(reader.channels), dfltCh = 1; end
        end
        prevZsect = ctrls.zsect.Value;
        if prevZsect<=imsz(3), dfltZ = prevZsect; else, dfltZ = 1; end
        dfltTP = 1;
        
        % Update controls
        set(ctrls.channel,'String',reader.channels,'Value',dfltCh,...
            'Enable','on');
        set(ctrls.zsect,'String',arrayfun(@num2str,1:imsz(3),'uni',0),...
            'Value',dfltZ,'Enable','on');
        set(ctrls.tp,'String',arrayfun(@num2str,1:imsz(5),'uni',0),...
            'Value',dfltTP,'Enable','on');
        
        % Update image
        updateRawImage;
    end

    function updateRawImage(~,~)
        % Load image
        if isempty(reader) || isempty(ctrls.channel.Value)
            rawImg = ones(size(rawImg));
        else
            rawImg = reader.getTimepoint(ctrls.tp.Value,...
                ctrls.zsect.Value,'C',ctrls.channel.Value);
        end
        updateImage;
    end

    function updateImage(~,~)
        % Resize to match zoom
        img = double(rawImg);
        rotVal = rotVals(ctrls.rot.Value);
        if ctrls.flip.Value, img = flipud(img); end
        if rotVal~=0, img = rot90(img,rotVal); end
        img = imresize(img,zoomMult(ctrls.zoom.Value));
        
        % Update size of frame if image has changed size
        if any(size(img)~=[imgH,imgW])
            imgW = size(img,2); imgH = size(img,1);
            updateFrame;
        end
        % Normalise image for display
        imgMin = min(img(:)); imgMax = max(img(:));
        if imgMax-imgMin>0
            img = (img-imgMin)./(imgMax-imgMin);
        end
        img = img(:,:,ones(1,3));
        
        imgax.Children.CData = img;
        set(imgax,'XLim',[0.5,imgW+0.5],'YLim',[0.5,imgH+0.5]);
    end

    function updateTimePoint(~,~)
        ntps = reader.imageSize(5);
        slider_val = (ctrls.tp.Value-1)/(ntps-1);
        set(ctrls.tslider,'Value',slider_val);
        updateRawImage;
    end

    function updateTslider(~,~)
        slider_val = ctrls.tslider.Value;
        ntps = reader.imageSize(5);
        ctrls.tp.Value = round(slider_val*(ntps-1) + 1);
        updateRawImage;
    end

    function updateXYslider(~,~)
    end

    function toggleTraps(~,~)
        endTrapDraw;
        if ctrls.usetraps.Value
            ctrls.drawTrap.Enable = 'on';
            ctrls.finishDrawTrap.Enable = 'off';
        else
            ctrls.drawTrap.Enable = 'off';
            ctrls.finishDrawTrap.Enable = 'off';
            trapImg = ones(80,80,3);
        end
    end

    function startTrapDraw(~,~)
        ctrls.drawTrap.Enable = 'off';
        cPointer = get(fig,'pointer');
        set(fig,'pointer','crosshair','Units','pixels');
        waitforbuttonpress;
        trapRect = rbbox;
        set(fig,'pointer',cPointer);
        set(imgax,'Units','pixels');
        trapRect(1:2) = trapRect(1:2)-imgax.Position(1:2);
        trapRect(2) = imgax.Position(4)-trapRect(2)-trapRect(4);
        trapRect = trapRect/zoomMult(ctrls.zoom.Value);
        trapRect(3:4) = 2*floor(trapRect(3:4)/2)+1;
        trapRect(1:2) = round(trapRect(1:2));
        Xlim = trapRect([1,1])+[0,trapRect(3)-1];
        Ylim = trapRect([2,2])+[0,trapRect(4)-1];
        img = double(rawImg);
        if ctrls.flip.Value, img = flipud(img); end
        rotVal = rotVals(ctrls.rot.Value);
        if rotVal~=0, img = imrotate(img,90*rotVal); end
        trapImg = img(Ylim(1):Ylim(2),Xlim(1):Xlim(2));
        updateTrapImg;
        ctrls.drawTrap.Enable = 'on';
        trapchannel = ctrls.channel.String{ctrls.channel.Value};
    end

    function endTrapDraw(~,~)
        ctrls.drawTrap.Enable = 'on';
        ctrls.finishDrawTrap.Enable = 'off';
    end

    function updateTrapImg(~,~)
        P = trapimpanel.Position;
        tsz = size(trapImg); tsz = flip(tsz(1:2));
        trapScaling = min(P(3:4)./tsz);
        set(trapimgax,'Units','pixels',...
            'Position',[(P(3:4)-trapScaling*tsz)/2,trapScaling*tsz]);
        % Normalise image for display
        img = trapImg;
        imgMin = min(img(:)); imgMax = max(img(:));
        img = (img-imgMin)./(imgMax-imgMin);
        img = img(:,:,ones(1,3));
        
        trapimgax.Children.CData = img;
        set(trapimgax,'XLim',[0.5,size(trapImg,2)+0.5],...
            'YLim',[0.5,size(trapImg,1)+0.5]);
    end

    function checkValidity(~,~)
        if isempty(reader)
            errordlg('A valid image needs to be specified');
            return
        end
        if ~isfolder(ctrls.savedir.fileName)
            errordlg('A valid save folder needs to be specified');
            return
        end
        uiresume(fig);
    end

%% General helper functions

    function uiLabel(parent,label)
        %uiLabel(parent,label)
        %Convenience function for making standard labels
        uicontrol('Parent',parent,'Style','text','FontWeight','bold',...
            'HorizontalAlignment','left','String',label);
    end

end
