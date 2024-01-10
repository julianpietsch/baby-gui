function refreshPosImSeg(this,trap,~)

if ~this.ctrls.updateoverview.Value, return; end

if ~this.ctrls.showoutlines.Value
    % Following needs to be RGB for shared colormap/grayscale in single
    % figure:
    this.frame.pos.Children.CData = this.posIm;
    return
end

% Overlay segmentation masks on pos image
imW = size(this.posIm,2);
imH = size(this.posIm,1);
if this.cTimelapse.trapsPresent
    bbW = this.cTimelapse.cTrapSize.bb_width;
    bbH = this.cTimelapse.cTrapSize.bb_height;
else
    bbW = (imW-1)/2;
    bbH = (imH-1)/2;
end

ctp = this.tpInds(this.currentTimepoint);
trapInfoArray = this.cTimelapse.cTimepoint(ctp).trapInfo;
trapLocArray = this.cTimelapse.cTimepoint(ctp).trapLocations;

if nargin>1 && isnumeric(trap) && isscalar(trap) && trap>=1 && ...
        trap<=length(trapInfoArray) && round(trap)==trap
    traps = trap;
    imtemp = this.posSegIm;
else
    traps = 1:length(trapInfoArray);
    imtemp = this.posIm;
end

% Ensure that cellsToPlot is large enough
maxlabel = max(cellfun(@(x) max([x(:);0]),this.cellLabels));
if size(this.cTimelapse.cellsToPlot,2)<maxlabel
    this.cTimelapse.cellsToPlot(:,end+1:maxlabel) = false;
end

for ti=traps
    trapInfo = trapInfoArray(ti);
    trapLoc = trapLocArray(ti);
    
    % Determine cell labels
    clabs = [];
    if trapInfo.cellsPresent, clabs = trapInfo.cellLabel; end
    
    % Determine trap location
    if this.cTimelapse.trapsPresent
        y = round(trapLoc.ycenter);
        x = round(trapLoc.xcenter);
    else
        y = (imH+1)/2;
        x = (imW+1)/2;
    end
    yinds = y-bbH:y+bbH;
    yindsTrap = 1:2*bbH+1;
    yindsTrap = yindsTrap(yinds>=1 & yinds<=imH);
    yinds = yinds(yinds>=1 & yinds<=imH);
    xinds = x-bbW:x+bbW;
    xindsTrap = 1:2*bbW+1;
    xindsTrap = xindsTrap(xinds>=1 & xinds<=imW);
    xinds = xinds(xinds>=1 & xinds<=imW);
    
    if isfield(trapInfo.cell,'segmented')
        segims = arrayfun(@(x) logical(x.segmented(yindsTrap,xindsTrap)),...
            trapInfo.cell,'UniformOutput',false);
    else
        segims = trapInfo.segmented(yindsTrap,xindsTrap);
    end
    
    % refreshPosImage sets this.posIm to be HxWx3
    tIm = this.posIm(yinds,xinds,:);
    for col=1:2
        sIm = tIm(:,:,col);
        for c=1:length(clabs)
            % cell labels of 0 may arise in old expts due to a baby bug:
            if clabs(c) == 0, continue; end
            if (this.cTimelapse.cellsToPlot(ti,clabs(c)) && col==2) || ...
                    (~this.cTimelapse.cellsToPlot(ti,clabs(c)) && col==1)
                sIm(segims{c}) = sqrt(sIm(logical(segims{c})));
            end
        end
        tIm(:,:,col) = sIm;
    end
    imtemp(yinds,xinds,:) = tIm;
end

this.posSegIm = imtemp;

% Following needs to be RGB for shared colormap/grayscale in single
% figure:
this.frame.pos.Children.CData = imtemp;

end