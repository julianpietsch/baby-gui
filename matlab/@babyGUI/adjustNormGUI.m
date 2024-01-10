function adjustNormGUI(this)

%% Retrieve default dimensions
dims = this.dimensions;
sp = dims.sp; % Spacing between controls
lh = dims.lh; % Label height
sh = dims.sh; % Selection menu height
ch = 3*dims.ch; % Control height

cw = 200;

%% Initialise dialog and controls
scrsz = get(0,'ScreenSize');
fig = dialog('Name','Adjust normalisation...',...
    'WindowStyle','normal','Visible','on',...
    'Position',[2*scrsz(3)/3,scrsz(4)/3,scrsz(3)/3,ch]);

mainbox = uix.HBox('Parent',fig,'Padding',sp,'Spacing',sp);
chctrl = uicontrol('Parent',mainbox,'Style','listbox',...
    'String',this.currentChannels,'Value',1,'Callback',@channel_update);

sliderbox = uix.VBox('Parent',mainbox,'Padding',0,'Spacing',sp);
uix.VBox('Parent',sliderbox);

minbox = uix.HBox('Parent',sliderbox,'Padding',0,'Spacing',sp);
uicontrol('Parent',minbox,'Style','text',...
    'HorizontalAlignment','left','String','Min:');
minslider = uicontrol('Style','slider','Parent',minbox);
addlistener(minslider,'Value','PostSet',@slider_update);
set(minbox,'Widths',[40,-1]);

maxbox = uix.HBox('Parent',sliderbox,'Padding',0,'Spacing',sp);
uicontrol('Parent',maxbox,'Style','text',...
    'HorizontalAlignment','left','String','Max:');
maxslider = uicontrol('Style','slider','Parent',maxbox);
addlistener(maxslider,'Value','PostSet',@slider_update);
set(maxbox,'Widths',[40,-1]);

uix.VBox('Parent',sliderbox);
set(sliderbox,'Heights',[-1,sh,sh,-1]);

set(mainbox,'Widths',[cw,-1]);

channel_update;

    function channel_update(~,~)
        channel = regexprep(this.currentChannels{chctrl.Value},'\W','_');
        minval = 0; maxval = 1;
        if isfield(this.channelNormalisation,channel)
            minval = this.channelNormalisation.(channel)(1);
            maxval = this.channelNormalisation.(channel)(2);
        end
        set(minslider,'Value',max(minval,0));
        set(maxslider,'Value',min(maxval,1));
    end

    function slider_update(~,~)
        channel = regexprep(this.currentChannels{chctrl.Value},'\W','_');
        this.channelNormalisation.(channel) = ...
            [get(minslider,'Value'),get(maxslider,'Value')];
        this.refreshTrapIms;
        this.refreshSegIms;
    end

end