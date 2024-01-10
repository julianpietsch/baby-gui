function exportCurated(cExperiment,varargin)
%EXPORTCURATED Export curated traps from an babyExperiment object
%
%   EXPORTCURATED(CEXPERIMENT) exports all trap/time-point pairs that have
%   been marked as 'isCurated' in the 'marks' property of an
%   babyExperiment object CEXPERIMENT, or alternatively all traps and 
%   time points for an babyExperimentSamples object CEXPERIMENT. A 
%   folder selection dialog will appear to choose the folder into which the
%   brightfield and segmentation outline images will be saved. Brightfield
%   Z-sections and outlines from multiple cells will be tiled in images
%   named according to their originating position, trap and time point
%   (e.g., 'pos001_trap001_tp0001_brightfield.png').
%
%   EXPORTCURATED(...,'Folder',FOLDER) specifies the folder FOLDER into
%   which images should be saved rather than prompting with a dialog.
%
%   EXPORTCURATED(...,'Marks',MARKTYPE) specifies that the trap/time-point
%   pairs marked as MARKTYPE in the 'marks' property of CEXPERIMENT should
%   be exported rather than the default type 'isCurated'. If 'all' is
%   specified, then all traps and time points will be exported. If
%   'triplets' is specified then POSES, TRAPS and TIMEPOINTS should all be
%   equal length arrays where each index specifies a
%   position/trap/time-point triplet to export (see below).
%   
%   EXPORTCURATED(...,'Channel',CHANNEL) specifies the (case-insensitive)
%   channel name CHANNEL that should be exported instead of 'brightfield'.
%
%   EXPORTCURATED(...,'poses',POSES) filters the exported trap/time-point
%   pairs to only those whose position index is in the array POSES. If a 
%   MARKTYPE of 'triplets' has been specified, then this should either be a
%   scalar index, which will be repeated to match the lengths of TRAPS 
%   and/or TIMEPOINTS, or an array of indices specifying the position for
%   each position/trap/time-point triplet.
%
%   EXPORTCURATED(...,'traps',TRAPS) filters the exported trap/time-point
%   pairs to only those whose trap index is in the array TRAPS. If a
%   MARKTYPE of 'triplets' has been specified, then this should either be a
%   scalar index, which will be repeated to match the lengths of POSES 
%   and/or TIMEPOINTS, or an array of indices specifying the trap for
%   each position/trap/time-point triplet.
%
%   EXPORTCURATED(...,'timepoints',TIMEPOINTS) filters the exported 
%   trap/time-point pairs to only those whose time-point index is in the 
%   array TIMEPOINTS. If a MARKTYPE of 'triplets' has been specified, then 
%   this should either be a scalar index, which will be repeated to match 
%   the lengths of POSES and/or TRAPS, or an array of indices specifying 
%   the time-point for each position/trap/time-point triplet.
%
%   EXPORTCURATED(...,'exportFocusStacks',true) will additionally save any
%   annotated Z-section focus for each cell outline in the meta-data of
%   the exported segmentation outline images (that is, those with suffix
%   '_segoutline.png').
%
%   EXPORTCURATED(...,'weight',WEIGHT) will add the numeric scalar WEIGHT
%   as a property in the meta-data of the exported segmentation outline
%   images.
%
%   EXPORTCURATED(...,'pixel_size',PIXELSIZE) will add the numeric scalar 
%   PIXELSIZE as a property in the meta-data of the exported images.

ip = inputParser;
ip.addRequired('cExperiment',@(x) isa(x,'babyExperiment') && isscalar(x));
ip.addParameter('Folder',[],@(x) isempty(x) || isfolder(x));
ip.addParameter('Marks',[],@(x) isempty(x) || (ischar(x) && size(x,1)==1));
ip.addParameter('Channel','brightfield',@(x) isempty(x) ||...
    (ischar(x) && size(x,1)==1) || iscellstr(x) || isstring(x));
isnat = @(x) all(isnumeric(x)) && all(round(x)==x) && all(x>=0);
ip.addParameter('poses',[],@(x) isempty(x) || (isvector(x) && isnat(x)));
ip.addParameter('traps',[],@(x) isempty(x) || (isvector(x) && isnat(x)));
ip.addParameter('timepoints',[],@(x) isempty(x) || (isvector(x) && isnat(x)));
ip.addParameter('exportFocusStacks',false,@(x) islogical(x) && isscalar(x));
ip.addParameter('weight',[],@(x) isempty(x) || (isscalar(x) && isnumeric(x)));
ip.addParameter('pixel_size',[],@(x) isempty(x) || (isscalar(x) && isnumeric(x)));
ip.parse(cExperiment,varargin{:});

folder = ip.Results.Folder;
marktype = ip.Results.Marks;
channels = ip.Results.Channel;
poses = ip.Results.poses;
trapFilter = ip.Results.traps;
tpFilter = ip.Results.timepoints;
exportFocus = ip.Results.exportFocusStacks;
weight = ip.Results.weight;
pixel_size = ip.Results.pixel_size;

if ~isempty(channels)
    % Ensure that any char/string inputs are converted to cellstr format
    channels = cellstr(channels);
    for ch=1:numel(channels)
        chNum = find(strcmpi(cExperiment.channelNames,channels{ch}),1);
        assert(~isempty(chNum),...
            'experiment does not appear to have %s stacks',channels{ch});
        % get capitalisation from cExperiment:
        channels{ch} = cExperiment.channelNames{chNum};
    end
end

if isempty(poses)
    poses = 1:numel(cExperiment.dirs);
end

assert(all(ismember(poses,1:numel(cExperiment.dirs))),...
    'specified poses are outside the possible range');

isSamples = isa(cExperiment,'babyExperimentSamples');

if isempty(marktype)
    if isSamples
        marktype = 'all';
    else
        marktype = 'isCurated';
    end
end

if ~ismember(marktype,{'all','triplets'})
    assert(isstruct(cExperiment.marks),...
        'cExperiment has a corrupted marks property');
    assert(isfield(cExperiment.marks,marktype),...
        '"%s" marks have not been defined for this cExperiment',marktype);
    
    marks = cExperiment.marks.(marktype);

    assert(iscell(marks) && length(marks)==length(cExperiment.dirs),...
        '"%s" marks are corrupted',marktype);

    goodposes = cellfun(@(x) isempty(x) || (isnumeric(x) && size(x,2)==2),marks);
    assert(all(goodposes),'corrupted trap/tp pairs in positions: %s',...
        strjoin(cExperiment.dirs(~goodposes),', '));

    poses = intersect(poses,...
        find(~cellfun(@(x) isempty(x) || size(x,1)==0,marks)));
elseif strcmp(marktype,'triplets')
    ntriplets = max([numel(poses),numel(trapFilter),numel(tpFilter)]);
    if isscalar(poses), poses = poses(ones(ntriplets,1)); end
    if isscalar(trapFilter), trapFilter = trapFilter(ones(ntriplets,1)); end
    if isscalar(tpFilter), tpFilter = tpFilter(ones(ntriplets,1)); end
    assert(numel(poses)==ntriplets && numel(trapFilter)==ntriplets ...
        && numel(tpFilter)==ntriplets,...
        '"poses", "traps" and "timepoints" must all be of equal length when specifying triplets');
    pos_trap_tp = [poses(:),trapFilter(:),tpFilter(:)];
    poses = unique(poses);
    trapFilter = [];
    tpFilter = [];
end

use_ds = isa(cExperiment,'babyExperimentOmero') && ...
    ~(strcmp(marktype,'all') && isempty(trapFilter));
if use_ds
    ds = OmeroDataset(cExperiment.omeroDs,'meta',struct(...
        'acqfiledata',cExperiment.metadata.acq,...
        'logfiledata',rmfield(cExperiment.metadata,'acq')));
end

if isempty(folder)
    uiInfo = 'specify folder in which to export images';
    fprintf('%s\n',uiInfo);
    folder = uigetdir('',uiInfo);
end

assert(isfolder(folder),'specified folder could not be found');

for pos=poses(:)'
    cTimelapse = cExperiment.returnTimelapse(pos);
    chNums = zeros(size(channels));
    for ch=1:numel(channels)
        chNums(ch) = find(strcmpi(cTimelapse.channelNames,channels{ch}),1);
        assert(~isempty(chNums),'channel with %s stacks missing in pos %u',...
            channels{ch},pos);
    end

    fprintf('Position %s',cExperiment.dirs{pos});
    
    if any(cTimelapse.timepointsProcessed)
        tpInds = find(cTimelapse.timepointsProcessed);
    else
        tpInds = cTimelapse.timepointsToProcess;
    end
    ntimepoints = numel(tpInds);
    
    if strcmp(marktype,'all')
        traps = 1:numel(cTimelapse.cTimepoint(tpInds(1)).trapInfo);
        tps = 1:ntimepoints;
        [trapV,tpV] = meshgrid(traps,tps);
        trap_tp_pairs = [trapV(:),tpV(:)];
        alltptraps = {traps};
        alltptraps = alltptraps(ones(ntimepoints,1));
    else
        if strcmp(marktype,'triplets')
            posMask = pos_trap_tp(:,1)==pos;
            trap_tp_pairs = pos_trap_tp(posMask,[2,3]);
        else
            trap_tp_pairs = marks{pos};
        end
        traps = unique(trap_tp_pairs(:,1));
        tps = unique(trap_tp_pairs(:,2));
        alltptraps = arrayfun(@(tp) ...
            trap_tp_pairs(trap_tp_pairs(:,2)==tp,1),tps,'uni',0);
    end
    if ~isempty(tpFilter)
        [~,goodtps,~] = intersect(tpInds(tps),tpFilter);
        tps = tps(goodtps);
        alltptraps = alltptraps(goodtps);
        trap_tp_pairs = trap_tp_pairs(ismember(trap_tp_pairs(:,2),tps),:);
    end
    if ~isempty(trapFilter)
        traps = intersect(traps,trapFilter);
        alltptraps = cellfun(@(t) intersect(t,trapFilter),...
            alltptraps,'uni',0);
        trap_tp_pairs = trap_tp_pairs(ismember(trap_tp_pairs(:,1),traps),:);
    end
    
    % Uniquely allocate daughters for each tp
    l = struct(); % lineage info
    for ti=1:numel(traps)
        trap = traps(ti);
        cellMothers = cTimelapse.cellMothers(trap,:);
        mothers = setdiff(full(unique(cellMothers)),0);
        nmothers = numel(mothers);
        l(ti).mothers = mothers;
        l(ti).currentdaughter = zeros(nmothers,ntimepoints,'uint16');
        
        trapInfo = arrayfun(@(x) x.trapInfo(trap),...
            cTimelapse.cTimepoint(tpInds),'uni',0);
        cellLabels = cellfun(@(x) x.cellLabel(:)',trapInfo,'uni',0);
        cellLabels = setdiff(unique([cellLabels{:}]),0);
        cellLabelMap = zeros(max(cellLabels),1);
        cellLabelMap(cellLabels) = 1:numel(cellLabels);
        
        isConsumed = false(numel(cellLabels),1);
        for m=1:nmothers
            dLabels = find(cellMothers==mothers(m));
            currentDaughter = 0;
            for tt=1:ntimepoints
                clabs = trapInfo{tt}.cellLabel;
                mi = find(clabs==mothers(m),1);
                if isempty(mi), continue; end
                dNew = dLabels(ismember(dLabels(:),clabs) ...
                    & ~isConsumed(cellLabelMap(dLabels)));
                if ~isempty(dNew)
                    currentDaughter = dNew(1);
                    isConsumed(cellLabelMap(currentDaughter)) = true;
                end
                di = find(clabs==currentDaughter,1);
                if isempty(di)
                    dLabel = 0;
                else
                    dLabel = clabs(di);
                end
                l(ti).currentdaughter(m,tt) = dLabel;
            end
        end
    end
    
    % Choose image loading method based on whether an OmeroDataset is
    % available or not:
    if use_ds
        trapIndMap = zeros(max(traps),1);
        trapIndMap(traps) = 1:numel(traps);
        tpIndMap = zeros(max(tpInds),1);
        tpIndMap(tpInds) = 1:numel(tpInds);
        
        batch_size = 40; % process pairs in batches
        npairs = size(trap_tp_pairs,1);
        nbatches = ceil(npairs/batch_size);
        for b=1:nbatches
            fprintf('.');
            if mod(b,40)==0, fprintf('\n'); end
            inds = (b-1)*batch_size+1:min(npairs,b*batch_size);
            % TODO: the following does not correctly account for
            % channel offset
            trapIms = cell(size(chNums));
            for ch=1:numel(chNums)
                trapIms{ch} = ds.returnSparseTrapTPs(cTimelapse,...
                    trap_tp_pairs(inds,1),trap_tp_pairs(inds,2),...
                    cTimelapse.channelNames{chNums(ch)});
            end
            cat(3,trapIms{:});
            for i=1:numel(inds)
                ind = inds(i);
                trap = trap_tp_pairs(ind,1);
                tp = trap_tp_pairs(ind,2);
                if ~isempty(chNums), trapIm = trapIms(:,:,:,i);
                else, trapIm = []; end
                exportTrapTP(pos,trap,tp,tpIndMap(tp),...
                    trapIm,l(trapIndMap(trap)));
            end
        end
        fprintf('\n');
    else
        rawStack = cell(size(chNums));
        for tt=1:numel(tps)
            tp=tps(tt);
            fprintf('.');
            if mod(tt,40)==0, fprintf('\n'); end

            tpStack = cell(numel(chNums),1);
            for ch=1:numel(chNums)
                tpStack{ch} = cTimelapse.returnSingleTimepoint(tpInds(tp),chNums(ch),'stack');
                if isempty(rawStack{ch})
                    rawStack{ch} = cTimelapse.returnSingleTimepointRaw(tpInds(tp),chNums(ch));
                    if cTimelapse.image_rotation ~= 0
                        if mod(cTimelapse.image_rotation,90) == 0
                            rawStack{ch} = rot90(rawStack{ch},cTimelapse.image_rotation / 90);
                        end
                    end
                    assert(isequal(size(rawStack{ch}),size(tpStack{ch})),...
                        'image preprocessing should not change image size');
                    rawQs = quantile(double(rawStack{ch}(:)),0:0.1:1);
                    adjQs = quantile(tpStack{ch}(:),0:0.1:1);
                    if sum(rawQs-adjQs)<0.05*(rawQs(end)-rawQs(1))
                        warning('image preprocessing slightly changes intensity scale');
                    elseif ~isequal(rawQs,adjQs)
                        error('image preprocessing should not change intensity scale');
                    end
                end
            end
            tpStack = cat(3,tpStack{:});

            tptraps = alltptraps{tt};
            if ~isempty(chNums)
                trapIms = zeros([cTimelapse.trapImSize,numel(tptraps),...
                    size(tpStack,3)],'like',rawStack{1});
                for z=1:size(tpStack,3)
                    trapIms(:,:,:,z) = cTimelapse.returnTrapsFromImage(...
                        tpStack(:,:,z),tpInds(tp),tptraps);
                end
                trapIms = permute(trapIms,[1,2,4,3]);
            end

            for ti=1:numel(tptraps)
                trap = tptraps(ti);
                if ~isempty(chNums), trapIm = trapIms(:,:,:,ti);
                else, trapIm = []; end
                exportTrapTP(pos,trap,tpInds(tp),tp,trapIm,l(ti));
            end
        end
        fprintf('\n');
    end
end

    function exportTrapTP(pos,trap,tp,tpInd,trapIm,lInfo)
        if isSamples
            ntpp = size(cTimelapse.ds_traps,1);
            spi = ceil(trap/ntpp);
            sti = mod(trap-1,ntpp)+1;
            
            % Fill description structure
            desc = struct();
            desc.experimentID = cTimelapse.ds.id;
            desc.position = cTimelapse.ds_poses(spi);
            desc.trap = cTimelapse.ds_traps(sti,spi);
            desc.tp = cTimelapse.ds_tps(tpInd,sti,spi);
            
            subfolder = cExperiment.dirs{pos};
            if ~isfolder(fullfile(folder,subfolder))
                [success,msg,~] = mkdir(folder,subfolder);
                assert(success,msg);
            end
            
            fname_base = fullfile(folder,subfolder,...
                sprintf('pos%03u_trap%03u_tp%04u_',desc.position,...
                desc.trap,desc.tp));
        else
            % Fill description structure
            desc = struct();
            desc.experimentID = cExperiment.id;
            desc.position = pos;
            desc.trap = trap;
            desc.tp = tp;
            
            fname_base = fullfile(folder,...
                sprintf('pos%03u_trap%03u_tp%04u_',pos,trap,tp));
        end
        
        if ~isempty(pixel_size), desc.pixel_size = pixel_size; end
        
        % Output channel stack
        if ~isempty(channels)
            desc.channel = strjoin(channels,'_');
            save_tiled(trapIm,...
                strcat(fname_base,desc.channel,'.png'),desc);
        end

        % Output stack of cell outlines
        desc.channel = 'segoutlines';
        if ~isempty(weight), desc.weight = weight; end
        tInfo = cTimelapse.cTimepoint(tp).trapInfo(trap);
        if tInfo.cellsPresent && ~isempty(tInfo.cellLabel) ...
                && length(tInfo.cellLabel)==length(tInfo.cell)
            tiled_img = arrayfun(@(x) full(logical(x.segmented)),...
                tInfo.cell,'uni',0);
            tiled_img = cat(3,tiled_img{:});
            
            desc.cellLabels = tInfo.cellLabel;
            if exportFocus
                assert(isfield(tInfo.cell,'focusStack') ...
                    && ~any(cellfun(@isempty,{tInfo.cell.focusStack})),...
                    'Focus annotation missing in pos %u, trap %u, time point %u',...
                    pos,trap,tp);
                desc.focusStack = [tInfo.cell.focusStack];
            end
            
            % For each mother, specify the bud label if one is present
            desc.buds = zeros(1,numel(desc.cellLabels));
            for c=find(ismember(lInfo.mothers,tInfo.cellLabel))
                d = lInfo.currentdaughter(c,tpInd);
                if d>0 && ismember(d,tInfo.cellLabel)
                    motherInd = find(tInfo.cellLabel==lInfo.mothers(c));
                    if numel(motherInd)~=1
                        warning('more than one mother label found in trapInfo');
                    end
                    desc.buds(motherInd(1)) = d;
                end
            end
        else
            if tInfo.cellsPresent ...
                    || (~isempty(tInfo.cellLabel) && tInfo.cellLabel(1)>0) ...
                    || (~isempty(tInfo.cell) && ~isempty(tInfo.cell.cellCenter))
                warning('cells may be present in trap %u, tp %u but trapInfo is corrupted',...
                    trap,tp);
            end
            tiled_img = false(cTimelapse.trapImSize);
            desc.cellLabels = [];
            desc.buds = [];
        end
        save_tiled(tiled_img,strcat(fname_base,'segoutlines.png'),desc);
    end
end
