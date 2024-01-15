function EditTracking(this,xloc,yloc,leftclick)
% EditTracking(this,xloc,yloc)
%
% does various things based on the edit mode
%

%% Convert click coords to coords relative to the appropriate trap image
pad = this.dimensions.impad;
tW = this.tileWidth;
tH = this.tileHeight;
nT = this.tilesContext;
ctp = this.currentTimepoint;

% Adjustments to align pointer tip with pixel
xloc = xloc-0.6; yloc = yloc+0.6;

column = floor((xloc-pad)/(tW+pad)); % zero-indexed column
row = floor((yloc-pad)/(tH+pad)); % zero-indexed row
xind = pad+column*(tW+pad);
yind = pad+row*(tH+pad);
timepoint = ctp-nT+column;

if column<0 || row<0 || xloc-xind>=tW
    % clicked in padding (for debugging)
    % fprintf('clicked in padding tp=%u, ch=%u (%0.1f,%0.1f)\n',...
    %     column,row,xind,yind);
    return
end

dofa = this.dofocusannot;
if dofa
    if timepoint==ctp-1
        trap = this.currentTrap;
        if isempty(this.cellLabels{trap}) || isempty(this.currentCell)
            % Invalid if no cells in trap or no active cell selected
            return
        end
        cellLabel = this.cellLabels{trap}(this.currentCell);
        tp = this.tpInds(ctp);
        clbls = this.cTimelapse.cTimepoint(tp).trapInfo(trap).cellLabel;
        % Invalid if active cell not present at current time point:
        if ~ismember(cellLabel,clbls), return; end
        
        nS = sum(this.isFocusAnnotationChannel); % number of focus stacks
        shH = round(tH/6); % half height
        sH = shH*2+1; % height of stack sample around current cell
        mH = round(sH/2); % height of mid-stack markers
        faH = 4*nS*pad+(nS+1)*mH+nS*sH;
        faOff = round((size(this.bgdIm,1)-faH)/2);
        stack = (yloc-faOff-mH-pad)/(sH+mH+4*pad)+1;
        yind = (sH+mH+4*pad)*(stack-floor(stack));
        stack = floor(stack);
        if stack<1
            stack = 0.5;
        elseif stack>nS
            stack = nS+0.5;
        else
            if (yind>sH+2*pad && yind<sH+3*pad) || ...
                    (yind>sH+3*pad+mH)
                return
            elseif yind>=sH+3*pad
                stack = stack+0.5;
            end
        end
        
        ci = find(cellLabel==clbls,1);
        this.cTimelapse.cTimepoint(tp).trapInfo(trap).cell(ci).focusStack = stack;
        this.haschanged = true;
        missing = cellfun(@isempty,...
            {this.cTimelapse.cTimepoint(tp).trapInfo(trap).cell.focusStack});
        if any(missing)
            % Move on to the next cell that has no assigned focus
            this.currentCell = ...
                find(ismember(this.cellLabels{trap},clbls(missing)),1);
        end
        return
    elseif timepoint~=ctp
        % All timepoints but current are inactive
        return
    end
end

if yloc-yind>=tH || timepoint<1 || timepoint>this.ntimepoints
    % clicked in padding (for debugging)
    fprintf('clicked in padding tp=%u, ch=%u (%0.1f,%0.1f)\n',...
        column,row,xind,yind);
    return
end
    
% Convert time point to the cTimelapse version:
guitp = timepoint;
timepoint = this.tpInds(timepoint);

Cx=xloc-xind;
Cy=yloc-yind;
cellPt=[Cx Cy];

if this.zoomactive
    cellPt = (cellPt + this.zoomoffset)/this.zoomactivescale;
else
    cellPt = cellPt + this.srcOffset;
end

% For debugging:
% fprintf('  selected timepoint = %i\n',timepoint);
% fprintf('  Cx=%0.1f Cy=%0.1f  \n\n',Cx,Cy);

%% Do something based on the selected mode

trap = this.currentTrap;
trapInfo = this.cTimelapse.cTimepoint(timepoint).trapInfo(trap);
if isempty(this.cellLabels{trap})
    cellLabel = 0;
elseif isempty(this.currentCell)
    error('a cell needs to be selected...');
else
    cellLabel = this.cellLabels{trap}(this.currentCell);
end

if dofa
    foundCell = nearestCell;
    if ~isempty(foundCell)
        this.currentCell = foundCell;
    end
    return
end

TrapInfoMissing = [];

modenames = fieldnames(this.editModeMap);
modeval = modenames{this.ctrls.editmode.Value};
switch modeval
    case 'setactive'
        if ~leftclick
            % Select new cell
            foundCell = nearestCell;
            if ~isempty(foundCell)
                this.currentCell = foundCell;
                % Switch to edit outline mode
                this.ctrls.editmode.Value = find(strcmp(modenames,'editoutline'),1);
                this.dragtp = column-nT;
            end
        end
        
    case {'addoutline','addcell'}
        if leftclick
            if cellLabel==0, cellLabel = 1; end % create the first cell
            if strcmp(modeval,'addcell')
                % Add outline with a new cell label:
                cellLabel = max([this.cellLabels{trap}(:);0])+1;
            end
            
            % Abort if an outline has already been added for active cell
            % unless it is an auto-generated outline
            addmodenames = fieldnames(this.addModeMap);
            addmodeval = addmodenames{this.ctrls.addmode.Value};
            isautogen = strcmp(addmodeval,'autogen') && ...
                trapInfo_isautogen(trapInfo,cellLabel);
            if ismember(cellLabel,trapInfo.cellLabel) && ~isautogen
                return
            end
            
            logmsg(this.cTimelapse,...
                'Add new cell outline at (%0.0f,%0.0f) in trap %d at time point %d with cell label %d',...
                Cx,Cy,trap,timepoint,cellLabel);
            switch addmodeval
                case 'autogen'
                    % Auto generation of outlines works by interpolating
                    % from the nearest outline for the current cell
                    guiLabel = find(this.cellLabels{trap} == cellLabel,1);
                    ctrack = [];
                    if ~isempty(guiLabel)
                        % Determine this cell track in cTimelapse time points
                        ctrack = false(numel(this.cTimelapse.cTimepoint),1);
                        ctrack(this.tpInds(this.cellTracks(guiLabel,:))) = true;
                    end
                    if strcmp(modeval,'addcell') || ~any(ctrack)
                        % this.ensureCellIndex also ensures that trapInfo.cell struct
                        % is initialised correctly; NCI = new cell index
                        NCI = this.ensureCellIndex(trap,timepoint,cellLabel);
                        % Add a mid-sized cell at current time point
                        % this.setContourParams(centre,radii,angles,trap,tp,ci)
                        this.setContourParams(round(cellPt),8*ones(1,6),...
                            (0:5)*pi/3,trap,timepoint,NCI);
                    else
                        % Determine autogen state for this cell for all tps
                        isautogen = false(size(ctrack));
                        isautogen(ctrack) = arrayfun(...
                            @(x) trapInfo_isautogen(x.trapInfo(trap),cellLabel), ...
                            this.cTimelapse.cTimepoint(ctrack));
                        
                        % Make sure that the current outline is either
                        % absent or autogenerated:
                        assert(~ctrack(timepoint) || isautogen(timepoint));
                        validtrack = ctrack & ~isautogen;
                        
                        runAutoGen(validtrack,cellLabel);
                    end
                case 'asbud'
                    % this.ensureCellIndex also ensures that trapInfo.cell struct
                    % is initialised correctly
                    NCI = this.ensureCellIndex(trap,timepoint,cellLabel); % NCI = new cell index
                    % Add a small round cell of fixed size at the clicked point
                    % this.setContourParams(centre,radii,angles,trap,tp,ci)
                    this.setContourParams(round(cellPt),4*ones(1,6),...
                        (0:5)*pi/3,trap,timepoint,NCI);
                case 'asrod'
                    % this.ensureCellIndex also ensures that trapInfo.cell struct
                    % is initialised correctly
                    NCI = this.ensureCellIndex(trap,timepoint,cellLabel); % NCI = new cell index
                    % Add a vertical rod of fixed size at the clicked point
                    % NB: works best if spline mode is cartesian
                    % Determine node locations first in cartesian coordinates
                    sc_angles = [pi,pi/2,0];
                    xc = 3*[-1,cos(sc_angles),1,flip(cos(-sc_angles))];
                    yc = 3*[0,sin(sc_angles)+2,0,flip(sin(-sc_angles))-2];
                    dev = [-1,1i]*[xc;yc];
                    % Convert to radial coordinates:
                    % this.setContourParams(centre,radii,angles,trap,tp,ci)
                    this.setContourParams(round(cellPt),abs(dev),...
                        pi-angle(dev),trap,timepoint,NCI);
            end
            this.refreshCellLabels(this.currentTrap);
            this.refreshTracks;
        else %rightclick
            % Select new cell
            foundCell = nearestCell;
            if ~isempty(foundCell)
                this.currentCell = foundCell;
            end
        end
        
    case 'removeoutline'
        if leftclick
            % Ensure original outlines from segmentation are saved before
            % editing:
            this.saveOriginalOutlines(trap,timepoint);
            logmsg(this.cTimelapse,'Remove cell outline at (%0.0f,%0.0f) in trap %d at time point %d',Cx,Cy,trap,timepoint);
            cell_to_remove_index = this.cTimelapse.ReturnNearestCellCentre(timepoint,trap,round(cellPt));
            this.cTimelapse.removeCell(timepoint,trap,cell_to_remove_index);
            this.markCellOutlineEdit(trap,timepoint,[]);
            this.refreshCellLabels(this.currentTrap);
            this.refreshTracks;
        else %rightclick
            % Select new cell
            foundCell = nearestCell;
            if ~isempty(foundCell)
                this.currentCell = foundCell;
                % Switch to add outline mode
                this.ctrls.editmode.Value = find(strcmp(modenames,'addoutline'),1);
            end
        end
        
    case 'removecell'
        if leftclick
            % Find clicked cell
            [foundCell,foundCellInd] = nearestCell;
            delLabel = trapInfo.cellLabel(foundCellInd);
            if isempty(foundCell) || foundCell ~= this.currentCell || delLabel ~= cellLabel
                uiwait(errordlg('Only the active cell can be removed. Click on an active outline to remove that cell.',...
                    'Error removing cell'));
                return
            end
            TrapInfoMissing = [];
            TPS = this.cTimelapse.timepointsToProcess;
            for TP = min(TPS):max(TPS)
                if ~isempty(this.cTimelapse.cTimepoint(TP).trapInfo)
                    TPLabels = this.cTimelapse.cTimepoint(TP).trapInfo(trap).cellLabel;
                    delInd = find(TPLabels == delLabel);
                    if ~isempty(delInd)
                        % Ensure original outlines from segmentation
                        % are saved before removing any:
                        this.saveOriginalOutlines(trap,TP);
                        this.cTimelapse.removeCell(TP,trap,delInd);
                        this.markCellOutlineEdit(trap,TP,[]);
                    end
                else
                    TrapInfoMissing = [TrapInfoMissing TP];
                end
            end
            
            % Update the log
            logmsg(this.cTimelapse,'Remove cell with label %d in trap %d',delLabel,trap);
            this.haschanged = true;
            this.refreshCellLabels(this.currentTrap);
            this.refreshTracks;
        else %rightclick
            % Select new cell
            foundCell = nearestCell;
            if ~isempty(foundCell)
                this.currentCell = foundCell;
                % Switch to add outline mode
                this.ctrls.editmode.Value = find(strcmp(modenames,'addcell'),1);
            end
        end
        
    case 'editoutline'
        if leftclick
            % Current time point handled by createDraggable in draggable
            % outline editing modes:
            outlinemodenames = fieldnames(this.outlineModeMap);
            outlinemode = outlinemodenames{this.ctrls.outlinemode.Value};
            if timepoint~=ctp || strcmp('clickmatch',outlinemode)
                % change active contour result
                CellIndex = trapInfo.cellLabel == cellLabel;
                if any(CellIndex)
                    radii = trapInfo.cell(CellIndex).cellRadii;
                    
                    % sometimes the centre gets saved as an integer for some reason
                    centre = double(trapInfo.cell(CellIndex).cellCenter);
                    angles = trapInfo.cell(CellIndex).cellAngle;
                    radii =  BABYutil.edit_radii_from_point(...
                        cellPt,centre,radii',angles');
                    
                    % this.setContourParams(centre,radii,angles,trap,tp,ci)
                    this.setContourParams(centre,radii,angles,trap,timepoint,CellIndex);
                end
            end
        else %rightclick
            % Select new cell
            foundCell = nearestCell;
            if ~isempty(foundCell)
                this.currentCell = foundCell;
                this.dragtp = column-nT;
            end
        end
        
    case 'breaktracks'
        % Select new cell
        [foundCell,foundCellInd] = nearestCell;
        if ~isempty(foundCell)
            if leftclick
                % break tracking from the new cell so that the cell is
                % seperated from that tp forwards.
                oldLabel = trapInfo.cellLabel(foundCellInd);
                TrapInfoMissing = [];
                % new highest label - if new cells need to be added.
                newCellLabel = this.cTimelapse.returnMaxCellLabel(trap)+1;
                % Retrieve lineage annotations that need to be checked for
                % reassignment
                daughters = find(full(this.cTimelapse.cellMothers(trap,:))==oldLabel);
                d_tp1 = NaN(size(daughters));
                for TP = min(this.cTimelapse.timepointsToProcess):max(this.cTimelapse.timepointsToProcess)
                    if isempty(this.cTimelapse.cTimepoint(TP).trapInfo)
                        TrapInfoMissing = [TrapInfoMissing TP];
                        continue
                    end
                    
                    TPLabels = this.cTimelapse.cTimepoint(TP).trapInfo(trap).cellLabel;
                    d_appear = isnan(d_tp1) & ismember(daughters,TPLabels);
                    d_tp1(d_appear) = TP;
                    
                    if TP<timepoint, continue; end
                    
                    if any(d_appear)
                        this.cTimelapse.cellMothers(trap,daughters(d_appear)) = newCellLabel;
                    end
                    
                    if any(TPLabels == oldLabel)
                        this.cTimelapse.cTimepoint(TP).trapInfo(trap).cellLabel(TPLabels==oldLabel) = newCellLabel;
                    end
                    
                    if isfield(this.cTimelapse.cTimepoint(TP).trapInfo(trap),'cellLabelOriginal')
                        % Modify the original labels to stay in sync
                        % with current labels
                        TPLabels = this.cTimelapse.cTimepoint(TP).trapInfo(trap).cellLabelOriginal;
                        
                        if any(TPLabels == oldLabel)
                            this.cTimelapse.cTimepoint(TP).trapInfo(trap).cellLabelOriginal(TPLabels==oldLabel) = newCellLabel;
                        end
                    end
                end
                
                % Update the log
                logmsg(this.cTimelapse,'Modified cell with label %d in trap %d at timepoint %d to have label %d',...
                    oldLabel,trap,timepoint,newCellLabel);
                this.haschanged = true;
                this.refreshCellLabels(this.currentTrap);
                this.refreshTracks;
            else % right click
                % set active
                this.currentCell = foundCell;
                % Switch to swap tracks mode
                this.ctrls.editmode.Value = find(strcmp(modenames,'swaptracks'),1);
            end
        end
        
    case 'swaptracks'
        % change tracking
        [foundCell,foundCellInd] = nearestCell;
        if ~isempty(foundCell)
            if leftclick
                oldLabel = trapInfo.cellLabel(foundCellInd);
                % Abort if clicked cell is same as active
                if oldLabel == cellLabel, return; end
                TrapInfoMissing = [];
                % new highest label - if new cells need to be added.
                newCellLabel = this.cTimelapse.returnMaxCellLabel(trap)+1;
                % Retrieve lineage annotations that need to be checked for
                % reassignment (both clicked cell and active cell)
                d_clicked = find(full(this.cTimelapse.cellMothers(trap,:))==oldLabel);
                dc_tp1 = NaN(size(d_clicked));
                d_active = find(full(this.cTimelapse.cellMothers(trap,:))==cellLabel);
                da_tp1 = NaN(size(d_active));
                for TP = min(this.cTimelapse.timepointsToProcess):max(this.cTimelapse.timepointsToProcess)
                    if isempty(this.cTimelapse.cTimepoint(TP).trapInfo)
                        TrapInfoMissing = [TrapInfoMissing TP];
                        continue
                    end
                    
                    TPLabels = this.cTimelapse.cTimepoint(TP).trapInfo(trap).cellLabel;
                    
                    dc_appear = isnan(dc_tp1) & ismember(d_clicked,TPLabels);
                    dc_tp1(dc_appear) = TP;
                    da_appear = isnan(da_tp1) & ismember(d_active,TPLabels);
                    da_tp1(da_appear) = TP;
                    
                    if TP<timepoint, continue; end
                    
                    if any(dc_appear)
                        this.cTimelapse.cellMothers(trap,d_clicked(dc_appear)) = cellLabel;
                    end
                    if any(da_appear)
                        this.cTimelapse.cellMothers(trap,d_active(da_appear)) = newCellLabel;
                    end
                    
                    if any(TPLabels == oldLabel)
                        this.cTimelapse.cTimepoint(TP).trapInfo(trap).cellLabel(TPLabels==oldLabel) = cellLabel;
                    end
                    
                    if any(TPLabels == cellLabel)
                        this.cTimelapse.cTimepoint(TP).trapInfo(trap).cellLabel(TPLabels==cellLabel) = newCellLabel;
                    end
                    
                    if isfield(this.cTimelapse.cTimepoint(TP).trapInfo(trap),'cellLabelOriginal')
                        % Modify the original labels to stay in sync
                        % with current labels
                        TPLabels = this.cTimelapse.cTimepoint(TP).trapInfo(trap).cellLabelOriginal;
                        
                        if any(TPLabels == oldLabel)
                            this.cTimelapse.cTimepoint(TP).trapInfo(trap).cellLabelOriginal(TPLabels==oldLabel) = cellLabel;
                        end
                        
                        if any(TPLabels == cellLabel)
                            this.cTimelapse.cTimepoint(TP).trapInfo(trap).cellLabelOriginal(TPLabels==cellLabel) = newCellLabel;
                        end
                    end
                end
                
                % Update the log
                logmsg(this.cTimelapse,'Modified cell with label %d in trap %d at timepoint %d to have label %d',...
                    oldLabel,trap,timepoint,cellLabel);
                this.haschanged = true;
                this.refreshCellLabels(this.currentTrap);
                this.refreshTracks;
            else %rightclick
                % set active
                this.currentCell = foundCell;
            end
        end
        
    case 'selectcells'
        [~,foundCellInd] = nearestCell;
        if ~isempty(foundCellInd)
            this.cTimelapse.cellsToPlot(trap,trapInfo.cellLabel(foundCellInd)) = leftclick;
            this.haschanged = true;
        end
        
    case 'setdaughter'
        [foundCell,foundCellInd] = nearestCell;
        if ~isempty(foundCellInd)
            if leftclick
                this.cTimelapse.cellMothers(trap,trapInfo.cellLabel(foundCellInd)) = cellLabel;
                this.haschanged = true;
                this.refreshCellLabels(this.currentTrap);
                this.refreshTracks;
            else
                %rightclick
                % set active
                this.currentCell = foundCell;
            end
        end
    case 'unsetdaughter'
        [foundCell,foundCellInd] = nearestCell;
        if ~isempty(foundCellInd)
            if leftclick
                this.cTimelapse.cellMothers(trap,trapInfo.cellLabel(foundCellInd)) = 0;
                this.haschanged = true;
                this.refreshCellLabels(this.currentTrap);
                this.refreshTracks;
            else
                %rightclick
                % set active
                this.currentCell = foundCell;
            end
        end
end

if ~isempty(TrapInfoMissing)
    fprintf('\n \nWARNING!! Trap info missing for the following timepoints:')
    display(TrapInfoMissing)
    fprintf('\n \n')
end

    function [cellNum,cellInd] = nearestCell()
        cellNum = [];
        cellInd = this.cTimelapse.ReturnNearestCellCentre(timepoint,trap,cellPt);
        if ~isempty(cellInd)
            cellNum = find(this.cellLabels{trap}==trapInfo.cellLabel(cellInd),1);
        end
    end

    function isautogen = trapInfo_isautogen(trapInfo,cellLabel)
        CI = find(trapInfo.cellLabel==cellLabel,1);
        isautogen = ~isempty(CI) && trapInfo.cellsPresent && ...
            isfield(trapInfo.cell,'autogen') && ...
            isequal(trapInfo.cell(CI).autogen,true);
    end

    function acparams = getACparams(tp,cellLabel)
        tInf = this.cTimelapse.cTimepoint(tp).trapInfo(trap);
        acparams = struct();
        assert(tInf.cellsPresent);
        CI = find(tInf.cellLabel==cellLabel,1);
        assert(~isempty(CI));
        acparams.centre = double(tInf.cell(CI).cellCenter(:)');
        acparams.radii = double(tInf.cell(CI).cellRadii(:));
        acparams.angles = double(tInf.cell(CI).cellAngle(:));
    end

    function runAutoGen(validtrack,cellLabel)
        % Find nearest non-autogenerated outlines to use as templates
        lTPoffset = find(validtrack(timepoint-1:-1:1),1);
        rTPoffset = find(validtrack(timepoint+1:end),1);
        
        % At least one of the template time points should be non-empty:
        assert(~isempty(lTPoffset) || ~isempty(rTPoffset));
        
        % Retrieve info for templates
        if ~isempty(lTPoffset)
            lTP = timepoint-lTPoffset;
            lInfo = getACparams(lTP,cellLabel);
            lInfo.tp = lTP;
        end
        
        if ~isempty(rTPoffset)
            rTP = timepoint+rTPoffset;
            rInfo = getACparams(rTP,cellLabel);
            rInfo.tp = rTP;
        end
        
        % Extrapolate if either of the template points is absent
        if isempty(lTPoffset) || isempty(rTPoffset)
            cellNum = find(this.cellLabels{trap}==cellLabel,1);
            tpIndMap = zeros(numel(this.cTimelapse.cTimepoint),1);
            tpIndMap(this.tpInds) = 1:numel(this.tpInds);
            assert(all(tpIndMap(validtrack)>0));
            trackarea = NaN(size(validtrack));
            trackarea(validtrack) = arrayfun(@(tp) ...
                this.getCellArea(trap,tp,cellNum),...
                tpIndMap(validtrack));
            trackvol = 4/3*pi*sqrt(trackarea/pi).^3;
            % Strip NaN values
            isvalid = ~isnan(trackvol);
            trackvol = trackvol(find(isvalid,1):find(isvalid,1,'last'));
            if isempty(rTPoffset), trackvol = flip(trackvol); end
            volgrowth = diff(trackvol(:));
            % Estimate volume expansion from a weighted expectation with
            % the exponential distribution:
            maxTP = min(15,numel(volgrowth));
            weights = exppdf(0:maxTP-1,5);
            volgrowthest = sum(weights(:).*volgrowth(1:maxTP),'omitnan') ...
                /sum(weights(~isnan(volgrowth(1:maxTP))));
            
            if isempty(lTPoffset)
                lTP = timepoint-1;
                lInfo = struct();
                lInfo.tp = timepoint;
                lInfo.centre = round(cellPt);
                if isnan(volgrowthest) || volgrowthest<=0
                    r_est = 2;
                else
                    v_est = max(trackvol(1)-rTPoffset*volgrowthest,10);
                    r_est = (v_est/trackvol(1))^(1/3)*mean(rInfo.radii);
                    r_est = max(r_est,0.7);
                end
                lInfo.radii = r_est*ones(4,1);
                lInfo.angles = (0:3)'*pi/2;
            else
                assert(isempty(rTPoffset)); % This should be true by outer if clause
                rTP = timepoint+1;
                rInfo = lInfo; % extrapolate from leftTP
                rInfo.tp = timepoint;
                rInfo.centre = round(cellPt);
                % Since trackvol is flipped, we expect negative volgrowthest
                if isnan(volgrowthest) || volgrowthest>=0
                    r_ratio = 1;
                else
                    v_est = trackvol(1)-lTPoffset*volgrowthest;
                    r_ratio = min((v_est/trackvol(1))^(1/3),2);
                end
                rInfo.radii = r_ratio*rInfo.radii;
            end
        end
        
        % Ensure that left and right have same number of knots/angles:
        p_angles = lInfo.angles; % take number of knots from left
        rInfo.radii = BABYutil.get_radii_from_radii(...
            p_angles,rInfo.radii,rInfo.angles);
        rInfo.angles = p_angles;
        
        % Interpolate between left and right template time points
        % Linear vol interp: v = v_l + (v_r-v_l)*(t-t_l)/(t_r-t_l)
        newTPs = lTP+1:rTP-1;
        tpfrac = (newTPs-lInfo.tp)/(rInfo.tp-lInfo.tp);
        lVol = 4/3*pi*(lInfo.radii(:)).^3;
        rVol = 4/3*pi*(rInfo.radii(:)).^3;
        trep = ones(numel(newTPs),1);
        krep = ones(numel(p_angles),1);
        p_vols = lVol(:,trep) + tpfrac(krep,:) .* (rVol(:,trep)-lVol(:,trep));
        p_radii = (p_vols/pi*3/4).^(1/3);
        p_centre = lInfo.centre(trep,:) + tpfrac([1,1],:)' ...
            .* (rInfo.centre(trep,:)-lInfo.centre(trep,:));
        
        % Generate new interpolated segmentation masks
        for t=1:numel(newTPs)
            % this.ensureCellIndex also ensures that trapInfo.cell struct
            % is initialised correctly; CI = cell index
            CI = this.ensureCellIndex(trap,newTPs(t),cellLabel);
            % this.setContourParams(centre,radii,angles,trap,tp,ci)
            this.setContourParams(p_centre(t,:),p_radii(:,t),...
                p_angles,trap,newTPs(t),CI);
            this.cTimelapse.cTimepoint(newTPs(t)).trapInfo(trap).cell(CI).autogen = true;
        end
    end
end
