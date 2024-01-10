function segout = babySegment(cTimelapse,varargin)

ip = inputParser;
isnat = @(x) all(isnumeric(x)) && all(round(x)==x) && all(x>=0);
ip.addParameter('timepoints',[],@(x) isempty(x) || (isvector(x) && isnat(x)));
ip.addParameter('traps',[],@(x) isempty(x) || (isvector(x) && isnat(x)));
ip.addParameter('showGUI',true,@(x) islogical(x) && isscalar(x));
ip.addParameter('keeperrors',false,@(x) islogical(x) && isscalar(x));
ip.addParameter('keep_baprobs',false,@(x) islogical(x) && isscalar(x));
ip.addParameter('assign_mothers',true,@(x) islogical(x) && isscalar(x));
ip.addParameter('refine_outlines',true,@(x) islogical(x) && isscalar(x));
ip.addParameter('with_volumes',true,@(x) islogical(x) && isscalar(x));
ip.addParameter('verbose',false,@(x) islogical(x) && isscalar(x));
ip.parse(varargin{:});

traps = ip.Results.traps;
timepoints = ip.Results.timepoints;
showGUI = ip.Results.showGUI;
keeperrors = ip.Results.keeperrors;
keep_baprobs = ip.Results.keep_baprobs;
assign_mothers = ip.Results.assign_mothers;
refine_outlines = ip.Results.refine_outlines;
with_volumes = ip.Results.with_volumes;
verbose = ip.Results.verbose;

if isempty(timepoints)
    if isempty(cTimelapse.timepointsToProcess)
        timepoints = 1:length(cTimelapse.cTimepoint);
    else
        timepoints = cTimelapse.timepointsToProcess;
    end
end

if isempty(traps)
    traps = 1:length(cTimelapse.cTimepoint(timepoints(1)).trapInfo);
end

bb = cTimelapse.babyBrain;
assert(~strcmp(bb.status,'asleep'),...
    'A valid BabyBrain configuration needs to be specified');

% Ensure that server is available by ensuring session initialisation here
sessionid = bb.get_session;

% Wait for session
for i=1:100
    sessioninfo = cTimelapse.babyBrain.activeSessions;
    % By default, Matlab concatenates output into struct array if all sessions
    % have the same fields. Normalise to cell array:
    if isstruct(sessioninfo), sessioninfo = num2cell(sessioninfo); end
    if iscell(sessioninfo)
        sessioninfo = sessioninfo(cellfun(@(x) isstruct(x) && ...
            all(isfield(x,{'id','runner'})), sessioninfo));
        sessionids = cellfun(@(x) x.id,sessioninfo,'uni',0);
        sessioninfo = sessioninfo(strcmp(sessionids,sessionid));
        if numel(sessioninfo)>=1 && strcmp(sessioninfo{1}.runner,'ready')
            break
        end
    end
    if i==3, fprintf('Waiting for baby server to load'); end
    if mod(i,5) == 3
        fprintf('.');
    end
    pause(0.5);
end
fprintf('\n');

bbchs = bb.channel;
if ~iscell(bbchs), bbchs = cellstr(bbchs); end
channelNums = cellfun(@(ch) find(strcmpi(cTimelapse.channelNames,ch),1),bbchs,'uni',0);
assert(~any(cellfun(@isempty,channelNums)),...
    'this cTimelapse does not have the "%s" channel',strjoin(bbchs,'+'));
channelNums = cell2mat(channelNums);
nch = numel(channelNums);
assert(numel(bb.nstacks)==nch,...
    '"nstacks" of BabyBrain config does not match number of channels');
zs = bb.zstacks;
if isempty(zs)
    zs = cell(1,nch);
elseif ~iscell(zs)
    zs = {zs};
end
assert(numel(zs)==nch,...
    '"zstacks" of BabyBrain config does not match number of channels');

warning('off','MATLAB:nargchk:deprecated');

% Preload raw images to determine bit depth and confirm Z stacks
bitdepth = '';
nz = zeros(1,nch);
for c=1:nch
    tp_img = cTimelapse.returnSingleTimepointRaw(timepoints(1),channelNums(c));

    if isempty(bitdepth)
        bitdepth = class(tp_img);
    else
        assert(strcmp(bitdepth,class(tp_img)),...
            'All channels must have same bit depth');
    end
    assert(isa(tp_img,'uint8') || isa(tp_img,'uint16'),...
        'Raw images must be either 8- or 16-bit unsigned integers');
    if isa(tp_img,'uint8'), bitfunc = @uint8; else, bitfunc = @uint16; end

    nz(c) = size(tp_img,3);
    if ~isempty(zs{c})
        assert(all(zs{c}<=nz(c)),'Configured Z stacks are out of range');
        nz(c) = numel(zs{c});
    end
    assert(nz(c)==bb.nstacks(c),...
        'Number of Z stacks does not match configured BabyBrain input');
end

if showGUI, pgui = posOverviewGUI(cTimelapse); end

% baby will still be 'thinking' if the segmentation terminates early
cTimelapse.babyBrain.set_status('thinking');

servererror = false(1,numel(timepoints));
stack_size = [];

if assign_mothers
    cTimelapse.cellMothers(traps,:) = 0;
end

for t=0:numel(timepoints)
    % Processing and image retrieval are scheduled so that the image
    % for the next time point can be downloaded while the current image is
    % being processed, i.e.:
    % t = 0. download 1, queue 1
    % t = 1. download 2, retrieve 1, queue 2, process 1
    % t = 2. download 3, retrieve 2, queue 3, process 2
    % ...
    % t = N-1. download N, retrieve N-1, queue N, process N-1
    % t = N. retrieve N, process N
    
    % Download image for next tp
    if t<numel(timepoints)
        if verbose, t_start = tic; end
        next_tp = timepoints(t+1);
        
        tp_img = cell(1,nch);
        trap_stack = cell(1,nch);
        for c=1:nch
            % Process to get correct alignment (i.e., rotation/offset) with 
            % trap definitions and then convert back to raw bit depth
            tp_img{c} = cTimelapse.returnSingleTimepoint(next_tp,channelNums(c),'stack');
            tp_img{c} = bitfunc(round(tp_img{c}));
            if ~isempty(zs{c}), tp_img{c} = tp_img{c}(:,:,zs{c}); end
            assert(nz(c)==size(tp_img{c},3),'Time points have differing Z depths!');
        
            % Need to returnTrapsFromImage separately for each channel so
            % that the padding value is also channel-specific
            tsz = cTimelapse.trapImSize;
            trap_stack{c} = zeros(tsz(1),tsz(2),numel(traps),nz(c),bitdepth);
            for s=1:nz(c)
                trap_stack{c}(:,:,:,s) = ...
                    cTimelapse.returnTrapsFromImage(tp_img{c}(:,:,s),next_tp,traps);
            end
            % Dimension order is traps * x * y * z
            trap_stack{c} = permute(trap_stack{c},[3,1,2,4]);
            if isempty(stack_size)
                stack_size = arrayfun(@(x) size(trap_stack{c},x),1:3);
            else
                assert(isequal(stack_size,arrayfun(@(x) size(trap_stack{c},x),1:3)),...
                    'Time points and/or channels have differing numbers of traps!');
            end
        end
        trap_stack = cat(4,trap_stack{:});

        if showGUI
            if ~isvalid(pgui), error('terminated by user...'); end
            pgui.addTP(next_tp,tp_img{1}(:,:,1));
        end
        if verbose, fprintf('Image download took %.3fs\n',toc(t_start)); end
    end
    
    % Update logs
    if t>0
        tp = timepoints(t);
        
        % Trigger the TimepointChanged event for babyLogging
        babyLogging.changeTimepoint(cTimelapse,t);
        
        % Also log memory usage every 100 time points
        if mod(t,100)==0
            if showGUI
                logmsg(cTimelapse,'cTimelapse: %s; posOverviewGUI: %s\n',...
                    cTimelapse.getMemUsage,pgui.getMemUsage);
            else
                logmsg(cTimelapse,'cTimelapse: %s\n',cTimelapse.getMemUsage);
            end
        end
    end
        
    % Retrieve segmented output
    if t>0 && ~servererror(t)
        if verbose, t_start = tic; end
        try
            segout = bb.get_segmentation;
        catch err
            warning(err.identifier,'There was a server error: %s',err.message);
            servererror(t) = true;
            cTimelapse.cTimepoint(tp).server_error = true;
        end
        if verbose, fprintf('Segmentation retrieval took %.3fs\n',toc(t_start)); end
    end

    % Queue the next image for processing
    if t<numel(timepoints)
        if verbose, t_start = tic; end
        try
            bb.queue_image(trap_stack, cTimelapse.pixelSize,...
                assign_mothers && t == (numel(timepoints)-1), ...
                keep_baprobs, refine_outlines, with_volumes);
        catch err
            warning(err.identifier,'There was a server error: %s',err.message);
            servererror(t+1) = true;
            cTimelapse.cTimepoint(next_tp).server_error = true;
        end
        if verbose, fprintf('Image queuing took %.3fs\n',toc(t_start)); end
    end
    
    % Organise and store segmented output
    if t>0 && ~servererror(t)
        if verbose, t_start = tic; end
        assert(length(segout)==stack_size(1),...
            'Server returned the wrong number of traps...');
        trapSize = stack_size(2:3);
        
        % Normalise segout to a cell array:
        if isstruct(segout), segout = num2cell(segout); end
        
        if keep_baprobs
            cTimelapse.cTimepoint(tp).trapInfo(1).baProbs = [];
        end
        
        for tt=1:length(segout)
            segInf = segout{tt};
            trapNum = traps(tt);
            
            trapInf = cTimelapse.trapInfoTemplate;
            if keep_baprobs || isfield(cTimelapse.cTimepoint(tp).trapInfo,'baProbs')
                trapInf.baProbs = [];
            end
            
            if isfield(segInf,'error')
                trapInf.error = segInf.error;
            end
            
            if ~isempty(segInf.centres)
                trapInf.cellsPresent = true;
                ncells = size(segInf.centres,1);
                
                if isempty(segInf.radii)
                    segInf.radii = cell(ncells,1);
                end
                if isnumeric(segInf.radii)
                    segInf.radii = num2cell(segInf.radii,2);
                end
                
                if isempty(segInf.angles)
                    segInf.angles = cell(ncells,1);
                end
                if isnumeric(segInf.angles)
                    segInf.angles = num2cell(segInf.angles,2);
                end
                
                if keep_baprobs && isfield(segInf,'p_bud_assign')
                    trapInf.baProbs = single(segInf.p_bud_assign);
                end
                
                goodcells = true(ncells,1);
                for c=1:ncells
                    % NB: need to swap rows and columns for Matlab coords
                    centre = flip(segInf.centres(c,:)+1);
                    radii = reshape(segInf.radii{c},1,[]);
                    angles = pi/2-reshape(segInf.angles{c},1,[]);
                    if any(isnan(radii)) || any(isnan(angles)) || isempty(radii) || isempty(angles)
                        goodcells(c) = false;
                    else
                        try
                            %TODO: hacky code to determine if the spline is
                            %fit in radial or cartesian coordinates; this
                            %should instead be queried from the baby
                            %server for the given model
                            if startsWith(bb.modelset,'mm') || startsWith(bb.modelset,'ecoli_mothermachine')
                                [px,py] = BABYutil.cartesian_spline_from_radii(...
                                    double(radii),double(angles),double(centre),trapSize);
                                segoutline = ...
                                    BABYutil.px_py_to_logical(px,py,trapSize);
                            else
                                segoutline = BABYutil.get_outline_from_radii(...
                                    radii,angles,centre,trapSize);
                            end
                        catch
                            warning('outline could not be generated');
                            goodcells(c) = false;
                        end
                    end
                    if ~goodcells(c) && keeperrors
                        trapInf.cell(c).outline_error = true;
                        segoutline = false(trapSize);
                    end
                    if goodcells(c) || keeperrors
                        trapInf.cell(c).cellCenter = single(centre);
                        trapInf.cell(c).cellRadius = single(mean(radii));
                        trapInf.cell(c).cellRadii = single(radii);
                        trapInf.cell(c).cellAngle = single(angles);
                        trapInf.cell(c).segmented = sparse(segoutline);
                    end
                end
                
                if isfield(segInf,'cell_label')
                    seg_cell_label = segInf.cell_label;
                else
                    seg_cell_label = [];
                end
                
                save_ellipse = isfield(segInf,'ellipse_dims');
                if save_ellipse, ellipse_dims = segInf.ellipse_dims; end
                save_vols = isfield(segInf,'volumes');
                if save_vols, cell_volumes = segInf.volumes; end
                
                if keeperrors
                    trapInf.cellLabel = uint16(seg_cell_label(:)');
                    if save_ellipse, trapInf.ellipseDims = ellipse_dims; end
                    if save_vols, trapInf.cellVolumes = cell_volumes; end
                else
                    if any(goodcells)
                        trapInf.cell = trapInf.cell(goodcells);
                        nvalid = min([numel(seg_cell_label),numel(goodcells)]);
                        maxlbl = max(seg_cell_label);
                        seg_cell_label = seg_cell_label(1:nvalid);
                        seg_cell_label = seg_cell_label(goodcells(1:nvalid));
                        cLabel = zeros(1,numel(trapInf.cell));
                        nvalid = numel(seg_cell_label);
                        cLabel(1:nvalid) = seg_cell_label;
                        unLabelled = nvalid+1:numel(cLabel);
                        if ~isempty(unLabelled)
                            cLabel(unLabelled) = maxlbl + (1:numel(unLabelled));
                        end
                        trapInf.cellLabel = uint16(cLabel);
                        if save_ellipse
                            if isempty(ellipse_dims) || size(ellipse_dims,2) ~= 2
                                trapInf.ellipseDims = NaN(numel(trapInf.cell),2);
                            else
                                nvalid = min([size(ellipse_dims,1),numel(goodcells)]);
                                ellipse_dims = ellipse_dims(1:nvalid,:);
                                ellipse_dims = ellipse_dims(goodcells(1:nvalid),:);
                                ellipse_dims(nvalid+1:numel(trapInf.cell),:) = NaN;
                                trapInf.ellipseDims = ellipse_dims;
                            end
                        end
                        if save_vols
                            nvalid = min([numel(cell_volumes),numel(goodcells)]);
                            cell_volumes = cell_volumes(1:nvalid);
                            cell_volumes = cell_volumes(goodcells(1:nvalid));
                            cell_volumes(nvalid+1:numel(trapInf.cell)) = NaN;
                            trapInf.cellVolumes = cell_volumes;
                        end
                    else
                        trapInf.cellLabel = uint16([]);
                    end
                end
            end
            
            if assign_mothers && isfield(segInf,'mother_assign') ...
                    && ~isempty(segInf.mother_assign)
                m_assign = segInf.mother_assign;
                try
                    cTimelapse.cellMothers(trapNum,1:numel(m_assign)) = m_assign;
                catch err
                    fprintf('cellMothers assignment error caught for trap:\n');
                    disp(trapNum);
                    disp('Value of m_assign:');
                    disp(m_assign);
                    fprintf('m_assign class: %s; numel: %u\nError:\n',...
                        class(m_assign), numel(m_assign));
                    disp(err);
                    for s=1:numel(err.stack)
                        disp(err.stack(s));
                    end
                end
            end
            
            fnames = fieldnames(trapInf);
            for f=1:numel(fnames)
                cTimelapse.cTimepoint(tp).trapInfo(trapNum).(fnames{f}) = trapInf.(fnames{f});
            end
        end
        
        if isfield(cTimelapse.cTimepoint,'server_error')
            % In case there was an error in a previous run, reset the error
            % flag for this time point.
            cTimelapse.cTimepoint(tp).server_error = false;
        end
        cTimelapse.timepointsProcessed(tp) = true;
        if verbose, fprintf('Segmentation storage took %.3fs\n',toc(t_start)); end
    end
    
    if t>0 && showGUI
        if ~isvalid(pgui), error('terminated by user...'); end
        pgui.ctp = tp; drawnow;
    end
end

% Baby will be happy if it has successfully completed segmentation
cTimelapse.babyBrain.set_status('happy');

if showGUI, delete(pgui); end

end
