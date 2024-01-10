function trackCellsOverlapDivide(cTimelapse,timepoints,traps,overwrite)
if nargin<2 || isempty(timepoints)
    timepoints = cTimelapse.timepointsToProcess;
    timepoints = timepoints(cTimelapse.timepointsProcessed(timepoints));
end
if nargin<3 || isempty(traps)
    traps = 1:numel(cTimelapse.cTimepoint(timepoints(1)).trapInfo);
end
if nargin<4 || isempty(overwrite)
    overwrite = false;
end

look_back_tps = 4;

for trap=traps
    max_label = 0;
    % mothers = zeros(1,100); % TODO: make this work for subset timepoints
    cTimelapse.cellMothers(trap,:) = 0;
    mergedMDs = zeros(1,size(cTimelapse.cellMothers,2));
    cellMasks_prev = cell(1,0);
    cellAreas_prev = cell(1,0);
    for t=1:numel(timepoints)
        tp = timepoints(t);
        trapInfo = cTimelapse.cTimepoint(tp).trapInfo(trap);
        
        ncells = numel(trapInfo.cell);
        cellMasks = cell(1,ncells);
        cellAreas = NaN(1,ncells);
        for c=1:ncells
            cellMasks{c} = imfill(full(trapInfo.cell(c).segmented),'holes');
            cellAreas(c) = sum(cellMasks{c},'all');
        end
        
        if t==1 && ~overwrite && numel(trapInfo.cellLabel)==ncells
            cellMasks_prev = {cellMasks};
            cellAreas_prev = {cellAreas};
            max_label = max(trapInfo.cellLabel);
            continue
        end
        
        if size(cTimelapse.cellMothers,2)<max_label+ncells
            cTimelapse.cellMothers(:,end+100) = 0;
            mergedMDs(end+100) = 0;
        end
        cellCenters = vertcat(trapInfo.cell.cellCenter);
        cellLabels = zeros(1,ncells);
        for t_prev=1:look_back_tps
            if t-t_prev<1, continue; end
            tp_prev = timepoints(t-t_prev);
            trapInfo_prev = cTimelapse.cTimepoint(tp_prev).trapInfo(trap);
            cellImage_prev = zeros(cTimelapse.trapImSize);
            masks_prev = cellMasks_prev{end-t_prev+1};
            areas_prev = cellAreas_prev{end-t_prev+1};
            for c=1:numel(trapInfo_prev.cell)
                cellImage_prev(masks_prev{c}) = c;
            end
            for c=1:ncells
                if cellLabels(c)>0, continue; end
                pxIds = cellImage_prev(cellMasks{c});
                if ~any(pxIds>0), continue; end
                [candIds,~,uidx] = unique(pxIds(pxIds>0));
                candIdCounts = accumarray(uidx,1)';
                IoUs = candIdCounts./(cellAreas(c)+areas_prev(candIds)-candIdCounts);
                [~,bestidxs] = sort(IoUs,'descend');
                if numel(bestidxs)>1 && IoUs(bestidxs(2))>0.1
                    % Check if these cells are a merging M-D pair
                    bestlbls = double(trapInfo_prev.cellLabel(candIds(bestidxs)));
                    if cTimelapse.cellMothers(trap,bestlbls(1))==bestlbls(2)
                        mergedMDs(bestlbls(2)) = bestlbls(1);
                        cellLabels(c) = bestlbls(2);
                        continue
                    elseif cTimelapse.cellMothers(trap,bestlbls(2))==bestlbls(1)
                        mergedMDs(bestlbls(1)) = bestlbls(2);
                        cellLabels(c) = bestlbls(1);
                        continue
                    end
                end 
                bestId = candIds(bestidxs(1));
                bestLabel = double(trapInfo_prev.cellLabel(bestId));
                IoU = IoUs(bestidxs(1));
                if  IoU > 0.5
                    cellLabels(c) = bestLabel;
                elseif IoU > 0.1
                    cellLabels(c) = -bestLabel;
                end
            end
            [uniqLabels,~,uidx] = unique(abs(cellLabels));
            uniqLabelCounts = accumarray(uidx,1);
            for c=find(uniqLabelCounts==2 & uniqLabels(:)>0)'
                pairIdx = find(abs(cellLabels)==uniqLabels(c));
                [~,motherIdx] = min(cellCenters(pairIdx,2));
                mLabel = uniqLabels(c);
                if mergedMDs(mLabel)>0
                    dLabel = mergedMDs(mLabel);
                    mergedMDs(mLabel) = 0;
                else
                    max_label = max_label + 1;
                    dLabel = max_label;
                end
                cellLabels(pairIdx) = dLabel;
                cellLabels(pairIdx(motherIdx)) = mLabel;
                cTimelapse.cellMothers(trap,dLabel) = mLabel;
            end
            cellLabels(cellLabels<0) = 0;
        end
        if t<=look_back_tps, sidx = 1; else, sidx = 2; end
        cellMasks_prev = [cellMasks_prev(sidx:end),{cellMasks}];
        cellAreas_prev = [cellAreas_prev(sidx:end),{cellAreas}];
        unassigned = cellLabels==0;
        N_unassigned = sum(unassigned);
        cellLabels(unassigned) = double(max_label) + (1:N_unassigned);
        max_label = max_label + N_unassigned;
        cTimelapse.cTimepoint(tp).trapInfo(trap).cellLabel = uint16(cellLabels);
    end
end