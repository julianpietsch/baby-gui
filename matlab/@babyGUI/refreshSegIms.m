function refreshSegIms(this)

%% First ensure we update the window elements if necessary
if this.windowUpdateRequired, this.UpdateWindowElements; end

%% Assign local dimension variables
pad = this.dimensions.impad;
tW = this.tileWidth;
tH = this.tileHeight;
nT = this.tilesContext;
nC = length(this.currentChannels);
zoff = this.zoomoffset;
srcOff = this.srcOffset;

ctp = this.currentTimepoint;
tps = max(1,ctp-nT):min(this.ntimepoints,ctp+nT);

% Ensure dragtp is correctly bounded
this.dragtp = min([this.dragtp,this.ntimepoints-ctp,nT]);
this.dragtp = max([this.dragtp,1-ctp,-nT]);

modenames = fieldnames(this.editModeMap);
outlinemodenames = fieldnames(this.outlineModeMap);
editoutline = strcmp('editoutline',modenames{this.ctrls.editmode.Value}) && ...
    ~strcmp('clickmatch',outlinemodenames{this.ctrls.outlinemode.Value});

cc = this.currentCell;

% Initialise from background image
imtemp = this.bgdIm;

%% Precalculate dimensions for focus annotation tiles if active
dofa = this.dofocusannot;
if dofa
    isfachan = this.isFocusAnnotationChannel;
    
    nS = sum(isfachan); % number of focus stacks
    shW = round(tW/6); % half width
    sW = shW*2+1; % width of stack sample around current cell
    shH = round(tH/6); % half height
    sH = shH*2+1; % height of stack sample around current cell
    mH = round(sH/2); % height of mid-stack markers
    faW = 2*pad+sW;
    faH = 4*nS*pad+(nS+1)*mH+nS*sH;
    faOff = round((size(imtemp,1)-faH)/2);
    assert(faH<=size(imtemp,1) && faOff>=0,...
        'Focus annotation array does not fit in tile display...');
    
    ct = this.tpInds(ctp);
    trapInfo = this.cTimelapse.cTimepoint(ct).trapInfo(this.currentTrap);
    centre = round([tW/2,tH/2]); % centre of trap if no cell selected
    activeS = NaN;
    if ~isempty(cc) && trapInfo.cellsPresent 
        cellLabel = this.cellLabels{this.currentTrap}(cc);
        ccli = find(trapInfo.cellLabel==cellLabel,1);
        if ~isempty(ccli)
            centre = round(reshape(trapInfo.cell(ccli).cellCenter,1,[]));
            if isfield(trapInfo.cell,'focusStack') && ...
                    ~isempty(trapInfo.cell(ccli).focusStack)
                activeS = trapInfo.cell(ccli).focusStack;
            end
        end
    end
end

%% Overlay segmentation masks on trapIm and insert into background image
trap = this.currentTrap;
cLabels = this.cellLabels{trap};
cellLabelMap = zeros(max(cLabels),1);
cellLabelMap(cLabels) = 1:numel(cLabels);
cellMothers = full(this.cTimelapse.cellMothers(trap,cLabels));
linarrows = struct('start',{},'end',{},'colour',{});
dola = this.ctrls.showlinarrows.Value;
for t=1:numel(tps)
    tp=tps(t); tt=t;
    if dofa && tp==ctp-1, continue; end
    if dofa && tp<ctp-1, tt=t+1; end
    trapInfo = this.cTimelapse.cTimepoint(this.tpInds(tps(tt))).trapInfo(trap);
    clabs = [];
    if trapInfo.cellsPresent, clabs = trapInfo.cellLabel; end
    colinds = cellLabelMap(clabs);
    for ch=1:nC
        % Position in the final image that this will be placed
        xind = pad+(tp-ctp+nT)*(tW+pad);
        yind = pad+(ch-1)*(tH+pad);
        dozoom = this.zoomactive && tp==(ctp+this.dragtp);
        
        tIm = this.trapIms{ch}(:,:,tt);
        fIm = tIm(:,:,ones(3,1)); % convert to RGB
        if this.ctrls.showoutlines.Value
            for col=1:3
                sIm = tIm;
                for c=1:length(colinds)
                    isdrag = editoutline && isequal(cc,colinds(c)) ...
                        && tp==(ctp+this.dragtp);
                    if dofa
                        isdrag = false;
                        if isequal(cc,colinds(c))
                            % Active cell is coloured according to distance
                            % to focus slice
                            colval = 0.8;
                            if ~isnan(activeS) && isfachan(ch)
                                s = sum(isfachan(1:ch));
                                colval = 0.8-0.5*(1/(1+20*(abs(s-activeS)/nS)^4));
                            end
                            if col==2, colval=1; end
                        elseif isfield(trapInfo.cell,'focusStack') ...
                                && ~isempty(trapInfo.cell(c).focusStack)
                            if col==3, colval=1; else, colval=0.8; end
                        else
                            if col==1, colval=1; else, colval=0.8; end
                        end
                    else
                        colval = this.currentColours(colinds(c),col);
                    end
                    if ~isdrag && ~(dofa && ctp~=tp)
                        if c>numel(trapInfo.cell)
                            warning('bad trapInfo');
                        end
                        segmask = logical(trapInfo.cell(c).segmented);
                        if isfield(trapInfo.cell,'autogen') ...
                                && isequal(trapInfo.cell(c).autogen,true) ...
                                && ~dofa
%                             colval = rgb2hsv(colval);
%                             colval(2) = 0.5 + 0.5*colval(2);
%                             colval = hsv2rgb(colval);
                            sIm(segmask) = (0.35+0.65*colval)*sqrt(sIm(segmask));
                        else
                            sIm(segmask) = (0.5*colval) + 0.5*sqrt(sIm(segmask));
                        end
                    end
                end
                fIm(:,:,col) = sIm;
            end
            
            if dola && ~dozoom
                for c=1:length(colinds)
                    mlbl = cellMothers(cellLabelMap(clabs(c)));
                    if ismember(mlbl,clabs)
                        % Pre-calculate lineage arrow properties if requested,
                        % present and not in a zoomed view:
                        mi = find(clabs==mlbl,1);
                        % Offset lineage arrow according to srcOffset and 
                        % check if it points to a cell outside of the 
                        % viewable area:
                        sXY = trapInfo.cell(mi).cellCenter-this.srcOffset;
                        eXY = trapInfo.cell(c).cellCenter-this.srcOffset;
                        if any([sXY,eXY]<1) || any([sXY(1),eXY(1)]>tW) ...
                                || any([sXY(2),eXY(2)]>tH)
                            continue
                        end
                        m = cellLabelMap(mlbl);
                        if this.currentDaughterInd(m,tps(tt)) ~= colinds(c)
                            continue
                        end
                        ind = numel(linarrows) + 1;
                        offset = [xind,yind];
                        linarrows(ind).start = offset + sXY;
                        linarrows(ind).end = offset + eXY;
                        linarrows(ind).colour = this.currentColours(m,:);
                        linarrows(ind).active = this.currentDaughterInd(m,tps(tt)) == colinds(c);
                    end
                end
            end
        end
        
        if dozoom
            fIm = imresize(fIm,this.zoomactivescale);
            fIm = fIm(zoff(2)+(1:tH),zoff(1)+(1:tW),:);
        else
            fIm = fIm(srcOff(2)+(1:tH),srcOff(1)+(1:tW),:);
        end
        
        imtemp(yind+1:yind+tH,xind+1:xind+tW,:) = fIm;
    end
end

%% Add in focus annotation images if mode is active
if dofa
    xind = pad+(nT-1)*(tW+pad)+round(tW/2);
    xi = [centre(1)+(-shW:shW); 1:sW]; % source and target x indices
    xi = xi(:,xi(1,:)>0 & xi(1,:)<=tW);
    yi = [centre(2)+(-shH:shH); 1:sH]; % source and target y indices
    yi = yi(:,yi(1,:)>0 & yi(1,:)<=tH);
    
    tt = find(tps==ctp,1); % index of the active timepoint
    fachan = find(isfachan);
    
    % Set the colour for the first mid-stack marker
    imtemp(faOff+(1:mH),xind+(1:faW),:) = 0.8;
    if activeS<0.75
        imtemp(faOff+(1:mH),xind+(1:faW),2) = 1;
    end
    for s=1:numel(fachan)
        ch = fachan(s);
        
        tIm = 0.8*ones(sW,sH);
        tIm(yi(2,:),xi(2,:)) = this.trapIms{ch}(yi(1,:),xi(1,:),tt);
        tIm = tIm(:,:,ones(3,1)); % convert to RGB
        
        yind = faOff+mH+pad+(s-1)*(sH+mH+4*pad);
        
        % Set the border around this stack
        imtemp(yind+(1:sH+pad*2),xind+(1:faW),:) = 0.8;
        if (activeS>=s-0.25 && activeS<s+0.25)
            imtemp(yind+(1:sH+pad*2),xind+(1:faW),2) = 1;
        end
        
        % Set the stack image
        imtemp(yind+pad+(1:sH),xind+pad+(1:sW),:) = tIm;
        
        % Set the colour for the next mid-stack marker
        imtemp(yind+3*pad+sH+(1:mH),xind+(1:faW),:) = 0.8;
        if (s==numel(fachan) && activeS>=s+0.25) || ...
                (activeS>=s+0.25 && activeS<s+0.75)
            imtemp(yind+3*pad+sH+(1:mH),xind+(1:faW),2) = 1;
        end
    end
end

%% Update in GUI
this.frame.im.CData = imtemp;

% Draw lineage arrows if mode is active
narrows = numel(linarrows);
nhandles = numel(this.arrowHandles);
if narrows < nhandles
    for h=narrows+1:nhandles
        delete(this.arrowHandles{h});
    end
    this.arrowHandles = this.arrowHandles(1:narrows);
    nhandles = numel(this.arrowHandles);
end
if dola
    [this.arrowHandles{nhandles+1:narrows}] = deal([]);
    for l=1:numel(linarrows)
        la = linarrows(l);
        h = this.arrowHandles{l};
        if isempty(h) || ~isvalid(h)
            h = drawArrow(this.frame.imAx,...
                la.start,la.end,'MaxHeadSize',10);
            this.arrowHandles{l} = h;
        end
        if la.active, lw = 2; else, lw = 0.5; end
        set(h,'XData',la.start(1),'YData',la.start(2),...
            'UData',la.end(1)-la.start(1),'VData',la.end(2)-la.start(2),...
            'Color',la.colour,'LineWidth',lw);
    end
end

this.updateDraggable;
end

function h = drawArrow(ax,s,e,varargin)
h = quiver(ax,s(1),s(2),e(1)-s(1),e(2)-s(2),0,varargin{:});
end