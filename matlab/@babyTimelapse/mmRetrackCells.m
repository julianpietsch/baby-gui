function mmRetrackCells(cTimelapse,timepoints,traps,overwritefirst,inverted,gr_init,varargin)
if nargin<2 || isempty(timepoints)
    if isempty(cTimelapse.timepointsToProcess)
        timepoints = 1:length(cTimelapse.cTimepoint);
    else
        timepoints = cTimelapse.timepointsToProcess;
    end
end

if nargin<3 || isempty(traps)
    traps = 1:length(cTimelapse.cTimepoint(timepoints(1)).trapInfo);
end

if nargin<4 || isempty(overwritefirst)
    % Continue tracking from first time point by default, that is, do not 
    % amend labels for first tp
    overwritefirst = false;
end

if nargin<5 || isempty(inverted)
    inverted = false;
end

tp1 = timepoints(1);
othertps = find(cellfun(@numel,{cTimelapse.cTimepoint.trapInfo}) ...
    == numel(cTimelapse.cTimepoint(tp1).trapInfo));
othertps = othertps(~ismember(othertps,timepoints));
% If not overwriting first tp, then it needs to be included for determining
% the next new_lbl:
if ~overwritefirst, othertps = [othertps,tp1]; end

if nargin<6 || isempty(gr_init)
    gr_init = 0.02;
end

ip = inputParser;
ip.addParameter('tol',[0.2,2],@(x) isvector(x) && isnumeric(x) && numel(x)==2 && all(x>0));
ip.addParameter('min_area',50,@(x) isscalar(x) && isnumeric(x));
ip.addParameter('max_gr',0.5,@(x) isscalar(x) && isnumeric(x));
ip.addParameter('Nlag',20,@(x) isscalar(x) && isnumeric(x) && round(x)==x && x>0);
ip.parse(varargin{:});

tol = ip.Results.tol;
max_gr = ip.Results.max_gr;
min_area = ip.Results.min_area;
min_area = min_area*0.065/cTimelapse.pixelSize;
Nlag = ip.Results.Nlag;

% Organise input data for parallelisation
trapdata = cell(1,numel(traps));
mothers = cell(1,numel(traps));
new_lbls = cell(1,numel(traps));
for ti=1:numel(traps)
    trap=traps(ti);
    
    % Determine the next new label from all pre-existing labels in this
    % trap that will not be overwritten
    new_lbl = max(arrayfun(@(x) max([x.trapInfo(trap).cellLabel(:);uint16(0)]),...
        cTimelapse.cTimepoint(othertps)))+1;
    if isempty(new_lbl), new_lbl = uint16(1); end
    new_lbls{ti} = new_lbl;

    % Reset mothers for any cells that we could create
    cTimelapse.cellMothers(trap,new_lbl:end) = 0;

    trapdata{ti} = arrayfun(@(x) x.trapInfo(trap),...
        cTimelapse.cTimepoint(timepoints),'Uniform',false);
    mothers{ti} = cTimelapse.cellMothers(trap,:);
end

cellLabels = cell(1,numel(traps));
parfor ti=1:numel(traps)
    [cellLabels{ti},mothers{ti}] = track_cells(trapdata{ti},mothers{ti},...
        new_lbls{ti},overwritefirst,inverted,gr_init,...
        tol,max_gr,min_area,Nlag);
end

maxcells = max(cellfun(@numel,mothers));
if size(cTimelapse.cellMothers,2)<maxcells
    cTimelapse.cellMothers(:,end+1:maxcells) = 0;
end
for ti=1:numel(traps)
    trap = traps(ti);
    cTimelapse.cellMothers(trap,1:numel(mothers{ti})) = mothers{ti};
    for t=1:numel(timepoints)
        tp = timepoints(t);
        cTimelapse.cTimepoint(tp).trapInfo(trap).cellLabel = cellLabels{ti}{t};
    end
end
    
end

function [ylocs,cAreas,cLengths,lb,ub]=get_features(trapInfo,inverted)
ncells = numel(trapInfo.cell);
ylocs = NaN(1,ncells);
cAreas = NaN(1,ncells);
cLengths = NaN(1,ncells);
lb = NaN(1,ncells);
ub = NaN(1,ncells);
for c=1:numel(trapInfo.cell)
    seg = imfill(full(trapInfo.cell(c).segmented),'holes');
    if ~any(seg(:))
        ylocs(c) = Inf;
        cAreas(c) = 0;
        cLengths(c) = 0;
        lb(c) = Inf;
        ub(c) = Inf;
        continue
    end
    rp = regionprops(seg,'Area','Centroid','MajorAxisLength','BoundingBox');
    if numel(rp)>1
        % Pick region with largest area if there are multiple
        [~,largest_area] = max([rp.Area]);
        rp = rp(largest_area);
    end
    ylocs(c) = rp.Centroid(2);
    cAreas(c) = rp.Area;
    cLengths(c) = rp.MajorAxisLength;
    if inverted
        ub(c) = -rp.BoundingBox(2);
        lb(c) = -sum(rp.BoundingBox([2,4]));
    else
        lb(c) = rp.BoundingBox(2);
        ub(c) = sum(rp.BoundingBox([2,4]));
    end
end
if inverted, ylocs = -ylocs; end
end

function [cellLabels,mothers]=track_cells(trapInfos,mothers,new_lbl,...
    overwritefirst,inverted,gr_init,tol,max_gr,min_area,Nlag)

cellLabels = cell(1,numel(trapInfos));

% Precalculate the tolerance bounds
ltm = diff(log(tol))/2; ltb = ltm-log(tol(2));

% Obtain features for first time point to initialise state
trapInfo = trapInfos{1};
[ylocs,cAreas,cLengths,lb,ub]=get_features(trapInfo,inverted);

% Label in order of y index starting from base of well (descending)
[~,yOrder] = sort(-ylocs);
% Ignore any outlines under the minimum size threshold
yOrder = yOrder(cAreas(yOrder)>=min_area);

% Initialise starting state
l_state = cLengths(yOrder);
lb_state = lb(yOrder);
ub_state = ub(yOrder);
if isempty(ub_state)
    lb0_state = 0;
else
    lb0_state = ub_state(1);
end
grs = gr_init(ones(1,numel(yOrder)));
merge_err_state = zeros(1,numel(yOrder));

% Reinitialise cellLabels if requested
if overwritefirst
    Nc = numel(trapInfo.cell);
    cLabels = zeros(1,Nc,'uint16');
    cLabels(yOrder) = new_lbl+uint16(0:numel(yOrder)-1);
    new_lbl = new_lbl + numel(yOrder);
    debris = cLabels==0;
    cLabels(debris) = new_lbl+uint16(0:sum(debris)-1);
    new_lbl = new_lbl+sum(debris);
else
    cLabels = trapInfo.cellLabel;
end
cellLabels{1} = cLabels;

lbl_state = cLabels(yOrder);
debris = cAreas<min_area;
Nd = sum(debris);
debris_state_labels = NaN(1,Nd+20);
debris_state_ylocs = NaN(1,Nd+20);
debris_state_labels(1:Nd) = cLabels(debris);
debris_state_ylocs(1:Nd) = ylocs(debris);

% Loop over all subsequent time points
for t=2:numel(trapInfos)
    % Obtain features for this time point
    trapInfo = trapInfos{t};
    [ylocs,cAreas,cLengths,lb,ub]=get_features(trapInfo,inverted);
    
    % Work in order of y position
    [~,yOrder] = sort(-ylocs);
    % Ignore any outlines under the minimum size threshold
    yOrder = yOrder(cAreas(yOrder)>=min_area);
    
    % If there are no outlines (or none over the minimum size
    % threshold), simply assign any small objects as new cells/debris
    % and skip to next time point:
    if isempty(yOrder)
        if trapInfo.cellsPresent
            ndebris = numel(trapInfo.cell);
            cellLabels{t} = new_lbl+uint16(0:ndebris-1);
            new_lbl = new_lbl + ndebris;
        end
        continue
    end
    
    cLengths = cLengths(yOrder);
    lb = lb(yOrder); ub = ub(yOrder);
    Ns = numel(l_state); Nt = numel(yOrder);
    % Max size of proposal is sum of cells in source and target lists,
    % since all have a chance to disappear or appear:
    Nmax = Ns + Nt;
    
    % Generate weight vectors for the previous and current time points
    % from the lengths...
    Wl = [[l_state(:);Inf*ones(Nmax-Ns,1)],...
        [cLengths(:);Inf*ones(Nmax-Nt,1)]];
    % ...and neighbour separation of each cell
    Wn = inf(Nmax,2);
    if Ns>0
        Wn(1:Ns,1) = ([lb0_state,lb_state(1:end-1)]...
            -[ub_state(2:end),lb_state(end)]);
    end
    if Nt>0
        Wn(1:Nt,2) = [lb0_state,lb(1:end-1)]-[ub(2:end),lb(end)];
    end
    % Wmats = {Wl,Wn};
    Wmats = {Wl};
    grs = min(grs,max_gr);
    grspad = [grs(:);inf(Nmax-Ns,1)];
    % Wl_pred = exp(grspad).*Wl(:,1);
    
    % Start with single basic proposal
    proposals = {[(1:Nmax)'*[1,1],ones(Nmax,2)]};
    for r=1:Nmax
        % Loop over all existing proposals and generate candidate
        % children proposals for current 'row'
        new_proposals = cell(1,numel(proposals));
        for pp=1:numel(proposals) % parent proposal
            parent_proposal = proposals{pp};
            rVals = [r,r];
            si = find(parent_proposal(:,1)==rVals(1),1);
            ti = find(parent_proposal(:,2)==rVals(2),1);
            
            % Always add proposals for cell to appear/disappear:
            fixed_proposals = cell(1,2);
            % For appearance, shift start labels up by one and fill top
            % with NaN so we can count number of appearances:
            if ~isempty(ti) && ti <= Nt
                proposal = parent_proposal;
                if ~isempty(si)
                    proposal(si:end,1) = [proposal(si+1:end,1);NaN];
                end
                % Finally, set weight type for the appearing cell to be the
                % special 0 case:
                proposal(ti,4) = 0;
                fixed_proposals{1} = proposal;
            end
            % For disappearance, shift current labels up by one,
            % filling with NaN and setting weight type as for appearance:
            if ~isempty(si) && si <= Ns
                proposal = parent_proposal;
                if ~isempty(ti)
                    proposal(ti:end,2) = [proposal(ti+1:end,2);NaN];
                end
                proposal(si,3) = 0;
                fixed_proposals{2} = proposal;
            end
            fixed_proposals = fixed_proposals(~cellfun(@isempty,fixed_proposals));
            
            if isempty(si) || isempty(ti) || si > Ns || ti > Nt
                new_proposals{pp} = [{parent_proposal},fixed_proposals];
                continue;
            end
            
            rInds = [si;ti];
            grmat = [grspad,-grs(si)*ones(Nmax,1)];
            candidates = cell(1,numel(Wmats));
            errors = cell(1,numel(Wmats));
            for w=1:numel(Wmats)
                W = Wmats{w}; W(1:si-1,1) = 0; W(1:ti-1,2) = 0;
                rref = flip(W(sub2ind(size(W),rInds,[1;2])));
                rref = rref(:)';
                err = log(cumsum(W))+grmat-ones(Nmax,1)*log(rref);
                % W_cum = cumsum(W);
                % err = [abs(W_cum(:,1)-W_cum(ti,2)),...
                %     abs(W_cum(:,2)-grs(si)*W_cum(si,1))]<tol_init;
                errInds = find(abs(err(:)+ltb)<ltm);
                [I,J] = ind2sub(size(err),errInds);
                err = err(errInds);
                errors{w} = abs(err(:));
                candidates{w} = [I,J];
            end
            % Sort candidates according to score and eliminate
            % duplicates in favour of those with the lowest scores
            wtypes = cellfun(@(w,c) w*ones(size(c,1),1),...
                num2cell(1:numel(candidates)),candidates,'uni',0);
            candidates = vertcat(candidates{:});
            wtypes = vertcat(wtypes{:});
            [~,err_order] = sort(vertcat(errors{:}));
            assert(isvector(err_order));
            candidates = candidates(err_order,:);
            wtypes = wtypes(err_order);
            % `unique` picks first occurrence (lowest error) by default
            [candidates,uInds] = unique(candidates,'rows');
            wtypes = wtypes(uInds);
            % Eliminate invalid candidates
            valid = candidates(:,1)>=rInds(candidates(:,2)) & ...
                ~all(candidates==[rInds(2),2],2);
            candidates = candidates(valid,:);
            wtypes = wtypes(valid);
            Ncandidates = size(candidates,1);
            new_proposals{pp} = [cell(1,Ncandidates),fixed_proposals];
            for p=1:Ncandidates
                pr = candidates(p,1); pc = candidates(p,2);
                proposal = parent_proposal;
                % Shift labels down by proposed jump...
                proposal(pr+1:end,pc) = ...
                    proposal(rInds(pc)+1:end-pr+rInds(pc),pc);
                % ...and repeat label according to proposed jump
                proposal(rInds(pc):pr,pc) = rVals(pc);
                % Update weight scheme selected for this proposal
                proposal(si,3) = wtypes(p); proposal(ti,4) = wtypes(p);
                proposal(rInds(pc):pr,pc+2) = wtypes(p);
                new_proposals{pp}{p} = proposal;
            end
        end
        new_proposals = [new_proposals{:}];
        if isempty(new_proposals), break; end
        proposals = new_proposals;
        scores = Inf(1,numel(proposals));
        for p=1:numel(proposals)
            proposal = proposals{p};
            % Determine weight matrix according to proposal
            wp = proposal(:,3:4);
            W = Wmats{1};
            % Occassionally a proposal ends up bigger in size than
            % weights... not sure how this happens, but we will score
            % these poorly:
            if ~isequal(size(wp),size(W)), continue; end
            for w=2:numel(Wmats)
                W(wp==w) = Wmats{w}(wp==w);
            end
            W(:,1) = exp(grspad).*W(:,1); % predicted expansion
            % Determine scores up to current 'row'
            ps = proposal(proposal(:,1)<=rVals(1),1);
            pt = proposal(proposal(:,2)<=rVals(2),2);
            score = abs(...
                log(accumarray(ps,W(1:numel(ps),1),[Nmax,1])) ...
                -log(accumarray(pt,W(1:numel(pt),2),[Nmax,1])));
            % Set penalty score for each appear/disappear event:
            score(accumarray(pt,proposal(1:numel(pt),4)==0,[Nmax,1])>0) = ltm;
            score(accumarray(ps,proposal(1:numel(ps),3)==0,[Nmax,1])>0) = ltm;
            % Adjust the cost for the last cell in the trap if it is
            % disappearing or appearing:
            if numel(ps) >= Ns && Ns > 0 && proposal(Ns,3) == 0
                score(ps(Ns)) = 0.5 * ltm;
            end
            if numel(pt) >= Nt && Nt > 0 && proposal(Nt,4) == 0
                score(pt(Nt)) = 0.5 * ltm;
            end
            % Determine score as an average over finite terms
            scores(p) = mean(score(isfinite(score)));
        end
        [~,order] = sort(scores);
        proposals = proposals(order(1:min(numel(order),10)));
    end
    % Best proposal is the first in the list
    proposal = proposals{1};
    % Determine weight matrix according to this proposal
    wp = proposal(:,3:4);
    W = Wmats{1};
    for w=2:numel(Wmats)
        W(wp==w) = Wmats{w}(wp==w);
    end
    
    % Update Ns to reflect any appear events in this proposal:
    Ns = Ns + sum(proposal(:,4)==0);
    % Nt = Nt + sum(proposal(:,3)==0);
    
    % We need to label all cells
    cLabels = zeros(1,numel(trapInfo.cell),'uint16');
    % And correct the growth rate estimates
    grs_next = NaN(1,Nt);
    % And carry over any merge errors
    merge_errs = zeros(1,Nt);
    
    % First try to match state labels to the current time point
    for r=1:Ns
        si = find(proposal(:,1)==r); Nsi = numel(si);
        ti = find(proposal(:,2)==r); Nti = numel(ti);
        if any(si>Ns) || any(ti>Nt) || any(si>numel(lbl_state)), break; end
        if Nsi == 1 && Nti == 1
            % Single cell tracks to single cell
            cLabels(yOrder(ti)) = lbl_state(si);
            cLengths(ti) = W(ti,2);
            gr = log(W(ti,2))-log(W(si,1));
            grs_next(ti) = (Nlag*grs(si)+max(gr,0))/(Nlag+1);
            merge_errs(ti) = merge_err_state(si);
        elseif Nsi == 1 && Nti == 2
            % Single cell divides into two
            cLabels(yOrder(ti(1))) = lbl_state(si);
            if merge_err_state(si) > 0
                % This cell previously merged, so keep previous label
                sister_label = merge_err_state(si);
                % NB: we now stop carrying over the merge error
                % NB: we do not alter the lineage tree
            else
                % Label the other daughter as new- the usual case
                sister_label = new_lbl; new_lbl = new_lbl + 1;
                % Update the lineage tree
                mothers(sister_label) = lbl_state(si);
                % NB: merge errors are not carried over in this case
            end
            cLabels(yOrder(ti(2))) = sister_label;
            cLengths(ti) = W(ti,2);
            gr = log(sum(W(ti,2)))-log(W(si,1));
            grs_next(ti) = (Nlag*grs(si)+max(gr,0))/(Nlag+1);
        elseif Nsi == 2 && Nti == 1
            % Two cells merge into one- this is a segmentation error,
            % but try to handle gracefully:
            cLabels(yOrder(ti)) = lbl_state(si(1));
            % Remember the absorbed cell in case division occurs later
            merge_errs(ti) = lbl_state(si(2));
            % Update length with proposed weight as per normal...
            cLengths(ti) = W(ti,2);
            % ...but just carry over the growth rate of cell 1
            grs_next(ti) = grs(si(1));
        else
            % If ambiguous, simply mark cells as new tracks...
            cLabels(yOrder(ti)) = new_lbl+uint16(0:Nti-1);
            new_lbl = new_lbl + Nti;
            % ...and carry over an averaged state
            if isempty(si)
                grs_next(ti) = gr_init;
            else
                grs_next(ti) = mean(grs(si));
            end
            % NB: any merge errors are reset
        end
        
    end
    
    % Assign cells under min_area to 'debris' tracks, which essentially
    % find closest yloc from previous time points:
    debris = find(cAreas<min_area);
    assert(all(cLabels(debris)==0));
    Nd = numel(debris);
    if Nd>numel(debris_state_labels)
        debris_state_labels(end+1:end+50) = NaN;
        debris_state_ylocs(end+1:end+50) = NaN;
    end
    debris_ylocs = ylocs(debris);
    debris_dist = abs(debris_ylocs(:)-debris_state_ylocs);
    debris_ind = NaN(1,Nd);
    for d=1:Nd
        [d_dist,ind] = nanmin(debris_dist(:));
        if isnan(d_dist), break; end
        [I,J] = ind2sub(size(debris_dist),ind);
        debris_ind(I) = J;
        debris_dist(I,:) = NaN;
        debris_dist(:,J) = NaN;
    end
    d_missed = isnan(debris_ind); Ndm = sum(d_missed);
    if Ndm>0
        dm_inds = find(isnan(debris_state_labels));
        dm_inds = dm_inds(1:Ndm);
        debris_ind(d_missed) = dm_inds;
        debris_state_labels(dm_inds) = new_lbl+uint16(0:Ndm-1);
        new_lbl = new_lbl + Ndm;
    end
    cLabels(debris) = debris_state_labels(debris_ind);
    
    % Mark any remaining unbalanced/extra cells as new tracks
    extra = cLabels==0;
    cLabels(extra) = new_lbl+uint16(0:sum(extra)-1);
    new_lbl = new_lbl + sum(extra);
    grs_next(isnan(grs_next)) = gr_init;
    
    % Update the state
    l_state = cLengths;
    ub_state = ub; lb_state = lb;
    lb0_state = (Nlag*lb0_state+ub(1))/(Nlag+1);
    grs = grs_next;
    lbl_state = cLabels(yOrder);
    merge_err_state = merge_errs;
    debris_state_ylocs(debris_ind) = debris_ylocs;
    
    % Finally update the cTimelapse
    cellLabels{t} = cLabels;
end
end