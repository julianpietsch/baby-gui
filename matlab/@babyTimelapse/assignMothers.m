function assignMothers(cTimelapse,traps,mintps)

timepoints = find(cTimelapse.timepointsProcessed);
assert(~isempty(timepoints),'timepoints have not been segmented yet');

assert(any(arrayfun(@(x) isfield(x.trapInfo,'baProbs'),...
    cTimelapse.cTimepoint(timepoints))),'not segmented with baby');

if nargin<2 || isempty(traps)
    traps = 1:length(cTimelapse.cTimepoint(timepoints(1)).trapInfo);
end

if nargin<3 || isempty(mintps)
    mintps = 3;
end

cTimelapse.cellMothers(traps,:) = 0;

for trap=traps(:)'
    trapInfo = arrayfun(@(x) x.trapInfo(trap),...
        cTimelapse.cTimepoint(timepoints),'uni',0);
    trapInfo = [trapInfo{:}];
    
    % Eliminate any cells whose duration is shorter than mintps
    cellLabels = setdiff(unique([trapInfo.cellLabel]),0);
    cellLabelMap = zeros(max(cellLabels),1);
    cellLabelMap(cellLabels) = 1:numel(cellLabels);
    cellDurations = zeros(numel(cellLabels),numel(trapInfo));
    for ti=1:numel(trapInfo)
        clabs = trapInfo(ti).cellLabel;
        if ~trapInfo(ti).cellsPresent || isempty(clabs) ...
                || (length(clabs)==1 && ismember(0, clabs))
            continue;
        end
        clabs = clabs(clabs > 0);
        ci = cellLabelMap(clabs);
        cellDurations(ci,ti) = 1;
    end
    cellDurations = sum(cellDurations,2);
    
    cellLabels = cellLabels(cellDurations>=mintps);
    cellLabelMap = zeros(max(cellLabels),1);
    cellLabelMap(cellLabels) = 1:numel(cellLabels);
    
    cumAssign = zeros(numel(cellLabels));
    cumIsMother = zeros(numel(cellLabels),1);
    cumIsDaughter = zeros(1,numel(cellLabels));
    for ti=1:length(trapInfo)
        clabs = trapInfo(ti).cellLabel;
        validlabs = ismember(clabs,cellLabels);
        clabs = clabs(validlabs);
        if ~trapInfo(ti).cellsPresent || isempty(clabs), continue; end
        baProbs = trapInfo(ti).baProbs(validlabs,validlabs);
        ci = cellLabelMap(clabs);
        cumIsMother(ci) = max(cumIsMother(ci), sum(baProbs,2));
        cumIsDaughter(ci) = max(cumIsDaughter(ci), max(baProbs));
        cumAssign(ci,ci) = cumAssign(ci,ci) + ...
            baProbs .* (1-cumIsMother(ci,ones(numel(ci),1))');
    end
    assignMat = exp(cumAssign);
    assignNorm = sum(assignMat,1);
    assignMat = assignMat ./ assignNorm(ones(size(assignMat,1),1),:);
    [~,assignments] = max(assignMat,[],1);
    isdaughter = cumIsDaughter > 0.5;
    if any(isdaughter)
        cTimelapse.cellMothers(trap,cellLabels(isdaughter)) = ...
            cellLabels(assignments(isdaughter));
    end
end

end