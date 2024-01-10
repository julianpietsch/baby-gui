function extractMotherDaughterArea(cTimelapse,calcEllipse)

if nargin<2 || isempty(calcEllipse), calcEllipse = false; end

timepoints = find(cTimelapse.timepointsProcessed);
assert(~isempty(timepoints),'timepoints have not been segmented yet');

if isempty(cTimelapse.extractedData)
    cTimelapse.extractedData = struct();
end

assert(isstruct(cTimelapse.extractedData),'extractedData is corrupt');

if all(isfield(cTimelapse.extractedData,{'trapNum','cellNum'})) ...
        && numel(cTimelapse.extractedData(1).trapNum)==numel(cTimelapse.extractedData(1).cellNum)
    % If extractedData has already been filled, overwrite existing
    trapNums = cTimelapse.extractedData(1).trapNum(:);
    cellNums = cTimelapse.extractedData(1).cellNum(:);
else
    if any(isfield(cTimelapse.extractedData,{'trapNum','cellNum'}))
        error('cTimelapse.extractedData structure is corrupt');
    end 
    % If extractedData is empty, then determine selected cells
    [trapNums,cellNums] = find(cTimelapse.cellsToPlot);
    trapNums = trapNums(:); cellNums = cellNums(:);
    %reorder so cells in the same trap contiguous, nicer for viewing later.
    [trapNums,I] = sort(trapNums); % NB: sort is stable, so cells should be ordered
    cellNums =cellNums(I);
    cTimelapse.extractedData.trapNum = trapNums;
    cTimelapse.extractedData.cellNum = cellNums;
end

% NB: to properly align missing time points between positions, we need to
% place results at the actual time point in the array, not the index into
% `timepoints`.
ncells = numel(trapNums);
ntps = timepoints(end);
mAreas = -ones(ncells,ntps,'int16');
dAreas = -ones(ncells,ntps,'int16');
currentDaughters = zeros(ncells,ntps,'uint16');
births = false(ncells,ntps);
mEllipseMajor = zeros(ncells,ntps);
mEllipseMinor = zeros(ncells,ntps);
dEllipseMajor = zeros(ncells,ntps);
dEllipseMinor = zeros(ncells,ntps);
mVols = zeros(ncells,ntps);
dVols = zeros(ncells,ntps);

traps = unique(trapNums);
for trap=traps(:)'
    inds = find(trapNums==trap);
    sel = cellNums(inds);
    
    trapInfo = arrayfun(@(x) x.trapInfo(trap),...
        cTimelapse.cTimepoint(timepoints),'uni',0);
    
    cellLabels = cellfun(@(x) x.cellLabel(:)',trapInfo,'uni',0);
    cellLabels = setdiff(unique([cellLabels{:}]),0);
    cellLabelMap = zeros(max(cellLabels),1);
    cellLabelMap(cellLabels) = 1:numel(cellLabels);
    
    isConsumed = false(numel(cellLabels),1);
    for c=1:numel(sel)
        ind = inds(c);
        dLabels = find(cTimelapse.cellMothers(trap,:)==sel(c));
        currentDaughter = 0;
        prevDArea = 0;
        for tt=1:numel(trapInfo)
            tp = timepoints(tt);
            clabs = trapInfo{tt}.cellLabel;
            mi = find(clabs==sel(c),1);
            if isempty(mi), continue; end
            mfill = connect_and_fill(trapInfo{tt}.cell(mi).segmented);
            mArea = sum(mfill(:));
            dNew = dLabels(ismember(dLabels(:),clabs) ...
                & ~isConsumed(cellLabelMap(dLabels)));
            if ~isempty(dNew)
                currentDaughter = dNew(1);
                isConsumed(cellLabelMap(currentDaughter)) = true;
                prevDArea = 0;
                births(ind,tp) = true;
            end
            
            di = find(clabs==currentDaughter,1);
            if isempty(di)
                dArea = prevDArea;
                dLabel = 0;
            else
                dfill = connect_and_fill(trapInfo{tt}.cell(di).segmented);
                dArea = sum(dfill(:));
                dLabel = clabs(di);
            end
            
            prevDArea = dArea;
            
            mAreas(ind,tp) = mArea;
            dAreas(ind,tp) = dArea;
            currentDaughters(ind,tp) = dLabel;
            
            if calcEllipse
                rp = regionprops(mfill,{'MajorAxisLength','MinorAxisLength'});
                mEllipseMajor(ind,tp) = rp.MajorAxisLength;
                mEllipseMinor(ind,tp) = rp.MinorAxisLength;
                if ~isempty(di)
                    rp = regionprops(dfill,{'MajoraxisLength','MinoraxisLength'});
                    dEllipseMajor(ind,tp) = rp.MajorAxisLength;
                    dEllipseMinor(ind,tp) = rp.MinorAxisLength;
                end
            elseif isfield(trapInfo{tt},'ellipseDims') ...
                    && size(trapInfo{tt}.ellipseDims,1)==numel(clabs)
                mEllipseMajor(ind,tp) = trapInfo{tt}.ellipseDims(mi,1);
                mEllipseMinor(ind,tp) = trapInfo{tt}.ellipseDims(mi,2);
                if ~isempty(di)
                    dEllipseMajor(ind,tp) = trapInfo{tt}.ellipseDims(di,1);
                    dEllipseMajor(ind,tp) = trapInfo{tt}.ellipseDims(di,2);
                end
            end
            
            if isfield(trapInfo{tt},'cellVolumes') ...
                    && numel(trapInfo{tt}.cellVolumes)==numel(clabs)
                mVols(ind,tp) = trapInfo{tt}.cellVolumes(mi);
                if ~isempty(di)
                    dVols(ind,tp) = trapInfo{tt}.cellVolumes(di);
                end
            end
        end
    end
end

cTimelapse.extractedData(1).area = mAreas;
cTimelapse.extractedData(1).daughterArea = dAreas;
cTimelapse.extractedData(1).daughterLabel = currentDaughters;
cTimelapse.extractedData(1).births = sparse(births);
cTimelapse.extractedData(1).mothers = sparse(cTimelapse.cellMothers(...
    sub2ind(size(cTimelapse.cellMothers),trapNums,cellNums)));

if any(mEllipseMajor(:)) || any(mEllipseMinor(:)) ...
        || any(dEllipseMajor(:)) || any(dEllipseMinor(:))
    cTimelapse.extractedData(1).ellipseMajor = sparse(mEllipseMajor);
    cTimelapse.extractedData(1).ellipseMinor = sparse(mEllipseMinor);
    cTimelapse.extractedData(1).daughterEllipseMajor = sparse(dEllipseMajor);
    cTimelapse.extractedData(1).daughterEllipseMinor = sparse(dEllipseMinor);
end

if any(mVols(:)) || any(dVols(:))
    cTimelapse.extractedData(1).volume = sparse(mVols);
    cTimelapse.extractedData(1).daughterVolume = sparse(dVols);
end

    function fillIm = connect_and_fill(edgeIm)
        fillIm = imfill(full(edgeIm),'holes');
    end
end