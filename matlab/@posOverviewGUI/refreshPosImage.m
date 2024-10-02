function refreshPosImage(this)

tp = this.ctp;

if ~this.loadedtps(tp), this.addTP(tp); end

alpha = this.alphascale;
posim = alpha*double(this.imcache(:,:,tp))/255;

imW = size(posim,2); imH = size(posim,1);

if this.cTimelapse.trapsPresent
    bbW = this.cTimelapse.cTrapSize.bb_width;
    bbH = this.cTimelapse.cTrapSize.bb_height;

    trapLocs = this.cTimelapse.cTimepoint(tp).trapLocations;
    nottraps = true(size(posim));
    for ti=1:length(trapLocs)
        xloc = round(trapLocs(ti).xcenter);
        yloc = round(trapLocs(ti).ycenter);
        trapXinds = max(xloc-bbW,1):min(xloc+bbW,imW);
        trapYinds = max(yloc-bbH,1):min(yloc+bbH,imH);
        nottraps(trapYinds,trapXinds) = false;
    end
    
    posim(nottraps) = posim(nottraps)*this.darkscale;
end

posim = posim(:,:,ones(3,1));

if ~this.showcells
    this.imgax.Children.CData = posim;
    return
end

% Overlay segmentation masks on pos image
trapInfoArray = this.cTimelapse.cTimepoint(tp).trapInfo;
trapLocArray = this.cTimelapse.cTimepoint(tp).trapLocations;

imtemp = posim;

cellCounter = 0;

for ti=1:length(trapInfoArray)
    trapInfo = trapInfoArray(ti);
    trapLoc = trapLocArray(ti);
    
    % Determine cell labels
    clabs = [];
    if trapInfo.cellsPresent, clabs = trapInfo.cellLabel; end
    
    if this.cTimelapse.trapsPresent
        % Determine trap location
        y = round(trapLoc.ycenter);
        x = round(trapLoc.xcenter);
        yinds = y-bbH:y+bbH;
        yindsTrap = 1:2*bbH+1;
        yindsTrap = yindsTrap(yinds>=1 & yinds<=imH);
        yinds = yinds(yinds>=1 & yinds<=imH);
        xinds = x-bbW:x+bbW;
        xindsTrap = 1:2*bbW+1;
        xindsTrap = xindsTrap(xinds>=1 & xinds<=imW);
        xinds = xinds(xinds>=1 & xinds<=imW);
    else
        yinds = 1:imH;
        xinds = 1:imW;
        yindsTrap = 1:imH;
        xindsTrap = 1:imW;
    end
    
    segims = arrayfun(@(x) logical(x.segmented(yindsTrap,xindsTrap)),...
        trapInfo.cell,'UniformOutput',false);
    
    cellinds = 1:length(clabs);
    trapColours = this.trackColours(...
        mod(cellCounter+clabs-1,this.MaxTrackColours)+1,:);
    
    if this.fillcells
        for c=cellinds
            segims{c} = imfill(full(segims{c}),'holes');
        end
    end

    % refreshPosImage sets this.posIm to be HxWx3
    tIm = posim(yinds,xinds,:);
    for col=1:3
        sIm = tIm(:,:,col);
        for c=cellinds
            sIm(segims{c}) = alpha*trapColours(c,col) + ...
                (1-alpha)*sqrt(sIm(segims{c}));
        end
        tIm(:,:,col) = sIm;
    end
    
    % cellCounter = cellCounter + length(clabs);
    
    imtemp(yinds,xinds,:) = tIm;
end

this.imgax.Children.CData = imtemp;

end