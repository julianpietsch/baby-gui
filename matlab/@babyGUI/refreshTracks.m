function refreshTracks(this)
%refreshTracks Refresh tracks and lineage tree for the current trap
trap = this.currentTrap;
cLabels = this.cellLabels{trap};
ncells = numel(cLabels);
cellLabelMap = zeros(max(cLabels),1);
cellLabelMap(cLabels) = 1:ncells;

% First refresh tracks
this.cellTracks_val = false(ncells,this.ntimepoints);
this.cellLocs_val = NaN(ncells,2,this.ntimepoints);
for t=1:length(this.tpInds)
    tp = this.tpInds(t);
    trapInfo = this.cTimelapse.cTimepoint(tp).trapInfo(trap);
    if ~trapInfo.cellsPresent, continue; end
    this.cellTracks_val(:,t) = ismember(cLabels,trapInfo.cellLabel);
    this.cellLocs_val(cellLabelMap(trapInfo.cellLabel),:,t) = ...
        vertcat(trapInfo.cell.cellCenter);
end
this.cellLocs_val = permute(this.cellLocs_val,[3,1,2]);

% Then reset any missing cellsToPlot
sel = find(this.cTimelapse.cellsToPlot(trap,:)>0);
missing_sel = sel(~ismember(sel,cLabels));
this.cTimelapse.cellsToPlot(trap,missing_sel) = 0;

% Then reset any cellLabels assigned to missing mothers or daughters
daughters = find(this.cTimelapse.cellMothers(trap,:) > 0);
mothers = this.cTimelapse.cellMothers(trap,daughters);
missingMDs = daughters(~ismember(daughters,cLabels)|~ismember(mothers,cLabels));
this.cTimelapse.cellMothers(trap,missingMDs) = 0;

% Then refresh lineage tree
this.currentDaughterInd = NaN(ncells,this.ntimepoints);
this.births = false(ncells,this.ntimepoints);
this.cTimelapse.cellMothers(:,end+1:max(cLabels)) = 0;
cellMothers = full(this.cTimelapse.cellMothers(trap,cLabels));
mothers = setdiff(unique(cellMothers),0);
currentDaughters = zeros(1,max(mothers));
isConsumed = false(ncells,1);
for t=1:numel(this.tpInds)
    tp = this.tpInds(t);
    trapInfo = this.cTimelapse.cTimepoint(tp).trapInfo(trap);
    if ~trapInfo.cellsPresent, continue; end
    clabs = trapInfo.cellLabel;
    cinds = cellLabelMap(clabs);
    currentMothers = cellMothers(cinds);
    hasmother = currentMothers > 0 & ~isnan(currentMothers);
    % if sum(hasmother)>2
    %     disp('has mother');
    % end
    cinds = cinds(hasmother); clabs = clabs(hasmother);
    currentMothers = currentMothers(hasmother);
    dNew = ~isConsumed(cinds);
    if any(dNew)
        dNew = find(dNew);
        [~,IA,~] = unique(currentMothers(dNew));
        dNew = dNew(IA);
        for d=dNew(:)'
            currentDaughters(currentMothers(d)) = clabs(d);
            isConsumed(cinds(d)) = true;
            this.births(cellLabelMap(currentMothers(d)),t) = true;
        end
    end
    cDpresent = find(currentDaughters(currentMothers)==clabs);
    for d=cDpresent(:)'
        this.currentDaughterInd(cellLabelMap(currentMothers(d)),t) = cinds(d);
    end
end

% for c=1:numel(mothers)
%     tic;
%     mLabel = mothers(c);
%     mi = cellLabelMap(mLabel);
%     dLabels = cLabels(cellMothers==mLabel);
%     currentDaughter = 0;
%     toc
%     tic;
%     isdlabel = false(ncells,1);
%     isdlabel(dLabels) = true;
%     for t=1:numel(this.tpInds)
%         tp = this.tpInds(t);
%         clabs = cellLabelMap(this.cTimelapse.cTimepoint(tp).trapInfo(trap).cellLabel);
%         dNew = isdlabel(clabs) & ~isConsumed(clabs);
%     end
%     toc
%     for t=1:numel(this.tpInds)
%         tp = this.tpInds(t);
%         clabs = this.cTimelapse.cTimepoint(tp).trapInfo(trap).cellLabel;
%         cinds = cellLabelMap(clabs);
%         dNew = isdlabel(cinds) & ~isConsumed(cinds);
%         if any(dNew)
%             currentDaughter = cinds(find(dNew,1));
%         dNew = dLabels(ismember(dLabels(:),clabs) ...
%             & ~isConsumed(cellLabelMap(dLabels)));
%         if ~isempty(dNew)
%             currentDaughter = dNew(1);
%             isConsumed(cellLabelMap(currentDaughter)) = true;
%             this.births(mi,t) = true;
%         end
%         if currentDaughter ~= 0 && ismember(currentDaughter,clabs)
%             this.currentDaughterInd(mi,t) = cellLabelMap(currentDaughter);
%         end
%     end
% end
end