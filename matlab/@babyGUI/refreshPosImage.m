function refreshPosImage(this,fromCache,~)

if ~this.ctrls.updateoverview.Value, return; end

if nargin<2 || ~islogical(fromCache) || ~isscalar(fromCache)
    fromCache = this.posImFromCache;
end

tp = this.tpInds(this.currentTimepoint);

if fromCache
    trapims = this.imcache.getTimepointTraps(...
        this.currentPos,this.currentTimepoint,...
        this.cTimelapse.channelNames{this.posImChannel},...
        [],'raw');
    posim = this.cTimelapse.returnWholeTrapImage(...
        trapims,tp);
    image_vals = posim(posim(:)~=0);
    posim_range = num2cell(quantile(image_vals,[0.0001,0.9999]));
    [posim_min,posim_max] = deal(posim_range{:});
else
    posim = this.cTimelapse.returnSingleTimepoint(...
        tp,this.posImChannel,'max');
    posim_range = num2cell(quantile(posim(:),[0.0001,0.9999]));
    [posim_min,posim_max] = deal(posim_range{:});
end

posim = (posim - posim_min)/(posim_max-posim_min);
posim(posim(:)<0) = 0;
posim(posim(:)>1) = 1;
posim = this.alphascale*posim;

if this.cTimelapse.trapsPresent
    halfW = round((this.srcWidth - 1)/2);
    halfH = round((this.srcHeight - 1)/2);
    maxW = size(posim,2); maxH = size(posim,1);
    trapLocs = this.cTimelapse.cTimepoint(tp).trapLocations;
    nottraps = true(size(posim));
    curTrapOutline = false(size(posim));
    for ti=1:length(trapLocs)
        xloc = round(trapLocs(ti).xcenter);
        yloc = round(trapLocs(ti).ycenter);
        trapXinds = max(xloc-halfW,1):min(xloc+halfW,maxW);
        trapYinds = max(yloc-halfH,1):min(yloc+halfH,maxH);
        nottraps(trapYinds,trapXinds) = false;
        if ti==this.currentTrap
            outerXinds = max(xloc-halfW-3,1):min(xloc+halfW+3,maxW);
            outerYinds = max(yloc-halfH-3,1):min(yloc+halfH+3,maxH);
            curTrapOutline(outerYinds,outerXinds) = true;
            curTrapOutline(trapYinds,trapXinds) = false;
        end
    end
    dark_scaling = 0.7;
    posim(nottraps) = posim(nottraps)*dark_scaling;
    posimcurtrap = posim;
    posimcurtrap(curTrapOutline) = 0.4+0.6*posimcurtrap(curTrapOutline)/dark_scaling;
    posim = cat(3,posimcurtrap,posimcurtrap,posim);
else
    posim = posim(:,:,ones(3,1));
end

this.posIm = posim;
this.refreshPosImSeg;

end