function [M,D] = grsFromCellInf(cellInf,varargin)
%GRSFROMCELLINF Growth rate estimation from cellInf or cellResults
%
%   GRSFROMCELLINF(CELLINF) estimates growth rates for both mothers and
%   their most current daughter from CELLINF, which can either be a
%   cExperiment.cellInf, cTimelapse.extractedData or a cellResults object.
%
%   GRSFROMCELLINF(...,'cellfilt',CELLFILT) applies the filter CELLFILT to
%   select the cells considered mothers.
%
%   GRSFROMCELLINF(...,'timepoints',TPS) filters all data (including 
%   extracted) to the time points TPS.
%
%   GRSFROMCELLINF(...,'vol_type',VOLTYPE) sets the volume estimation
%   method as one of 'conical', 'conical-corrected' (default), 'fromarea',
%   'ellipsoid', 'area', 'length' or 'rod'. The conical methods use the 
%   extracted volume directly. 'fromarea' estimates the volume from the 
%   area  assuming a spherical shape. The ellipsoid uses the extracted 
%   major and minor ellipse radii to calculate an ellipsoid volume. The
%   area and length methods use area or major axis length directly (with
%   dimensions to match). 'rod' uses major and minor axis lengths to 
%   compute the volume of a cylinder with spherical caps of radius half the
%   minor axis length.
%
%   GRSFROMCELLINF(...,'gr_method',GRMETHOD) uses the method GRMETHOD for
%   smoothing and estimating the growth rate. Available methods are
%   'sgolay', 'slmfixed' (default) and 'slmfree'.
%
%   GRSFROMCELLINF(...,'log_gr',true) specifies that growth rate estimation
%   should be performed on the log-volume, thus obtaining specific growth
%   rates.
%
%   GRSFROMCELLINF(...,'pixel_size',PXSIZE) specifies the pixel size of the
%   data (default = 0.263 um/pixel). The Prime95b camera with 60x zoom 
%   seems to have pixel size 0.182. If PXSIZE is specified in um, then 
%   volumes (and derived growth rates) will have units of um^3.
%
%   GRSFROMCELLINF(...,'full_daughters',true) specifies that daughters for 
%   the selected mothers have been included in the extracted data, so each 
%   daughter should be fit in full (i.e., not just for the currently active
%   daughter).
%
%   GRSFROMCELLINF(...,'usefilters',false) will not apply any track
%   filtering or joining (only applicable when 'full_daughters' is true).
%
%   GRSFROMCELLINF(...,'smallVolThresh',V,'smallGrowthThresh',G) specify 
%   the threshold growth rate G, measured as (max(volume)-min(volume))/dt
%   for each track, and threshold max volume V, below which tracks will be
%   discarded as non-growing debris.
%
%   GRSFROMCELLINF(...,'joinVolThresh',V) specifies the difference in
%   volume V between endpoints of contiguous tracks, below which those
%   tracks will be joined.
%
%   GRSFROMCELLINF(...,'shortTPthresh',T) specifies the minimum number of
%   time points duration a track must have before it is rejected.
%
%   GRSFROMCELLINF(...,'extract',PROPNAMES) specifies additional properties
%   to extract from the CELLINF and include in the output mother and
%   daughter arrays. PROPNAMES should be a cell string (cell array of
%   character vectors) specifying the field names from the CELLINF to
%   include, and optionally also prefixed by the channel name to use. E.g.,
%   the 'nucEstConv' fields for the GFP and mCherry channels could be
%   specified as: {'GFP.nucEstConv','mCherry.nucEstConv'}
%
%   GRSFROMCELLING(...,'birthHandling',METHOD) specifies how to treat
%   volumes across birth events in the calculation of growth rate. When
%   METHOD is 'none' (default) no special treatment is applied. When METHOD
%   is 'split', growth rates are fit independently for each cell division
%   cycle (i.e., between births). When METHOD is 'cumulative', then volumes
%   are recalculated as the cumulative sum over differences; across birth
%   events, the difference is between the sum of the mother and daughter 
%   after birth and the volume of the mother cell before birth (i.e., 
%   the birth event is handled like a mother splitting into two daughters).

vol_types = {'conical','conical-corrected','fromarea','ellipsoid',...
    'length','area','rod'};
grfuns = struct('sgolay',@grSGolay,'slmfixed',@grSLM,'slmfree',@grSLMfree);
mOpt = struct('fP',3,'fW',11,'sW',7,'knotfit','fixed');
% TODO: make daughters increasing only over the interval until their first
% daughter
dOpt = struct('fP',4,'fW',11,'sW',5,'increasing','on');
birthHandlingTypes = {'none','split','cumulative'};

ip = inputParser;
isindex = @(x) isvector(x) && isnumeric(x) && all(x>0) && all(round(x)==x);
ip.addParameter('cellfilt',[],@(x) isempty(x) || isindex(x) || islogical(x));
ip.addParameter('timepoints',[],@(x) isempty(x) || isindex(x) || islogical(x));
ip.addParameter('vol_type','conical-corrected',@(x) ismember(x,vol_types));
ip.addParameter('gr_method','slmfixed',@(x) ismember(x,fieldnames(grfuns)));
ip.addParameter('log_gr',false,@(x) isscalar(x) && islogical(x));
ip.addParameter('full_daughters',false,@(x) isscalar(x) && islogical(x));
ip.addParameter('pixel_size',0.263,@(x) isscalar(x) && isnumeric(x));
ip.addParameter('extract',{},@(x) iscellstr(x) || isstring(x));
ip.addParameter('mW',[],@(x) isempty(x) || isindex(x));
ip.addParameter('dW',[],@(x) isempty(x) || isindex(x));
ip.addParameter('usefilters',true,@(x) isscalar(x) && islogical(x));
ip.addParameter('smallVolThresh',7,@(x) isscalar(x) && isnumeric(x));
ip.addParameter('smallGrowthThresh',10,@(x) isscalar(x) && isnumeric(x));
ip.addParameter('joinVolThresh',7,@(x) isscalar(x) && isnumeric(x));
ip.addParameter('shortTPthresh',5,@(x) isscalar(x) && isnumeric(x));
ip.addParameter('birthHandling','none',@(x) ismember(x,birthHandlingTypes));
ip.parse(varargin{:});

cellfilt = ip.Results.cellfilt;
timepoints = ip.Results.timepoints;
vol_type = ip.Results.vol_type;
gr_method = ip.Results.gr_method;
log_gr = ip.Results.log_gr;
full_daughters = ip.Results.full_daughters;
pixel_size = ip.Results.pixel_size;
extract = cellstr(ip.Results.extract);
mW = ip.Results.mW;
dW = ip.Results.dW;
usefilters = ip.Results.usefilters;
smallVolThresh = ip.Results.smallVolThresh;
smallGrowthThresh = ip.Results.smallGrowthThresh;
joinVolThresh = ip.Results.joinVolThresh;
shortTPthresh = ip.Results.shortTPthresh;
birthHandling = ip.Results.birthHandling;

if ~isempty(mW)
    mOpt.fW = mW; mOpt.sW = mW;
    if mW<=mOpt.fP, mOpt.fP = mW-1; end
end
if ~isempty(dW)
    dOpt.fW = dW; dOpt.sW = dW;
    if dW<=dOpt.fP, dOpt.fP = dW-1; end
end

cr = cellInf; % save the original input for extraction
if isa(cr,'cellResults')
    cellInf = cr.cellInf.general;
    cellInf.posNum = cr.cellPosNum;
    cellInf.trapNum = cr.cellTrapNum;
    cellInf.cellNum = cr.cellNum;
    if isempty(cellfilt) && cr.filtercells
        cellfilt = cr.cellfilter;
    end
    if isempty(timepoints) && cr.filtertimes
        timepoints = cr.timefilter;
    end
end

cellInf = cellInf(1); % only ever need the first index for outline info

if isempty(cellfilt)
    cellfilt = true(size(cellInf.cellNum));
end

if ~isfield(cellInf,'posNum')
    % This typically occurs when this is run on cTimelapse.extractedData
    cellInf.posNum = ones(size(cellInf.cellNum));
end

if isempty(timepoints)
    timepoints = true(1,size(cellInf.times,2));
end

times = cellInf.times(cellfilt,timepoints)/60; % times in hours
births = full(cellInf.births(cellfilt,timepoints)>0);

% Determine indices to daughters
mInds = 1:numel(cellInf.posNum);
mInds = mInds(cellfilt);
ncells = numel(mInds);
dInds = cell(ncells,1);
dFilt = cell(ncells,1);
dNum = double(cellInf.daughterLabel(cellfilt,timepoints));
dNum(dNum==0) = NaN;
ndcells = 0;
for m=1:ncells
    mInd = mInds(m);
    dl = unique(dNum(m,~isnan(dNum(m,:))),'stable');
    if numel(dl) ~= sum(births(m,:))
        warning('daughter labels and births do not match for pos %u, trap %u, cell %u',...
            cellInf.posNum(mInd),cellInf.trapNum(mInd),cellInf.cellNum(mInd));
    end
    if full_daughters
        trapFilt = (cellInf.posNum(:) == cellInf.posNum(mInd)) & ...
            (cellInf.trapNum(:) == cellInf.trapNum(mInd));
        dFilt{m} = NaN(numel(dl),1);
        for d=1:numel(dl)
            dInd = find(trapFilt & (cellInf.cellNum(:) == dl(d)),1);
            if isempty(dInd)
                warning('daughter %u missing for pos %u, trap %u, cell %u',...
                    dl(d),cellInf.posNum(mInd),cellInf.trapNum(mInd),cellInf.cellNum(mInd));
            else
                dFilt{m}(d) = dInd;
            end
        end
        dl = dl(~isnan(dFilt{m}));
        dFilt{m} = dFilt{m}(~isnan(dFilt{m}));
    else
        dFilt{m} = mInd*ones(numel(dl),1);
    end
    dInds{m} = ndcells+(1:numel(dl));
    ndcells = ndcells+numel(dl);
end
dFilt = vertcat(dFilt{:});
dFilt = dFilt(~isnan(dFilt));

% Retrieve volumes
switch vol_type
    case 'conical'
        mVols = full(cellInf.volume(cellfilt,timepoints));
        mVols(mVols==0) = NaN;
        if full_daughters
            dVols = full(cellInf.volume(dFilt,timepoints));
        else
            dVols = full(cellInf.daughterVolume(dFilt,timepoints));
        end
        dVols(dVols==0) = NaN;
    case 'conical-corrected'
        mAreas = full(double(cellInf.area(cellfilt,timepoints)));
        mAreas(mAreas<0) = NaN;
        mVols = full(cellInf.volume(cellfilt,timepoints));
        mVols(mVols==0) = NaN;
        mVols = mVols + mAreas;
        if full_daughters
            dAreas = full(double(cellInf.area(dFilt,timepoints)));
            dVols = full(cellInf.volume(dFilt,timepoints));
        else
            dAreas = full(double(cellInf.daughterArea(dFilt,timepoints)));
            dVols = full(cellInf.daughterVolume(dFilt,timepoints));
        end
        dAreas(dAreas<0) = NaN;
        % dVols(dVols==0) = NaN;
        % NB: the above allows for cases where uncorrected conical 
        % volume -> 0 (i.e., for outlines that are < 9 pixels in size)
        dVols = dVols + dAreas;
    case 'fromarea'
        mAreas = full(double(cellInf.area(cellfilt,timepoints)));
        mAreas(mAreas<0) = NaN;
        mVols = 4/3*pi*sqrt(mAreas/pi).^3;
        if full_daughters
            dAreas = full(double(cellInf.area(dFilt,timepoints)));
        else
            dAreas = full(double(cellInf.daughterArea(dFilt,timepoints)));
        end
        dAreas(dAreas<0) = NaN;
        dVols = 4/3*pi*sqrt(dAreas/pi).^3;
    case 'ellipsoid'
        mEmaj = full(cellInf.ellipseMajor(cellfilt,timepoints)); mEmaj(mEmaj==0) = NaN;
        mEmin = full(cellInf.ellipseMinor(cellfilt,timepoints)); mEmin(mEmin==0) = NaN;
        mVols = pi/6*mEmaj.*mEmin.^2;
        if full_daughters
            dEmaj = full(cellInf.ellipseMajor(dFilt,timepoints));
            dEmin = full(cellInf.ellipseMinor(dFilt,timepoints));
        else
            dEmaj = full(cellInf.daughterEllipseMajor(dFilt,timepoints));
            dEmin = full(cellInf.daughterEllipseMinor(dFilt,timepoints));
        end
        dEmaj(dEmaj==0) = NaN;
        dEmin(dEmin==0) = NaN;
        dVols = pi/6*dEmaj.*dEmin.^2;
    case 'length'
        mVols = full(cellInf.ellipseMajor(cellfilt,timepoints)); mVols(mVols==0) = NaN;
        if full_daughters
            dVols = full(cellInf.ellipseMajor(dFilt,timepoints));
        else
            dVols = full(cellInf.daughterEllipseMajor(dFilt,timepoints));
        end
        dVols(dVols==0) = NaN;
        mVols = mVols / pixel_size^2; % anticipates volume conversion below
        dVols = dVols / pixel_size^2; % anticipates volume conversion below
    case 'area'
        mVols = full(double(cellInf.area(cellfilt,timepoints)));
        mVols(mVols<=0) = NaN;
        if full_daughters
            dVols = full(double(cellInf.area(dFilt,timepoints)));
        else
            dVols = full(double(cellInf.daughterArea(dFilt,timepoints)));
        end
        dVols(dVols<=0) = NaN;
        mVols = mVols / pixel_size; % anticipates volume conversion below
        dVols = dVols / pixel_size; % anticipates volume conversion below
    case 'rod'
        mEmaj = full(cellInf.ellipseMajor(cellfilt,timepoints)); mEmaj(mEmaj==0) = NaN;
        mEmin = full(cellInf.ellipseMinor(cellfilt,timepoints)); mEmin(mEmin==0) = NaN;
        mVols = pi*(mEmin/2).^2.*(mEmaj-mEmin)+4/3*pi*(mEmin/2).^3;
        if full_daughters
            dEmaj = full(cellInf.ellipseMajor(dFilt,timepoints));
            dEmin = full(cellInf.ellipseMinor(dFilt,timepoints));
        else
            dEmaj = full(cellInf.daughterEllipseMajor(dFilt,timepoints));
            dEmin = full(cellInf.daughterEllipseMinor(dFilt,timepoints));
        end
        dEmaj(dEmaj==0) = NaN;
        dEmin(dEmin==0) = NaN;
        dVols = pi*(dEmin/2).^2.*(dEmaj-dEmin)+4/3*pi*(dEmin/2).^3;
    otherwise
        error('unrecognised "vol_type"');
end

mVols = pixel_size^3 * mVols;
dVols = pixel_size^3 * dVols;

% Collect any requested fields to extract except those handled specially
extract = setdiff(extract,{'posNum','trapNum','cellNum','times','mothers'});
notformothers = {'births','daughterLabel'};
forextraction = struct('M',struct(),'D',struct());
if isa(cr,'cellResults')
    chnames = cr.channels;
elseif isfield(cr,'extractionParameters') ...
        && isfield(cr(1).extractionParameters.functionParameters,'channels')
    chnames = cr(1).extractionParameters.functionParameters.channels;
else
    chnames = strcat('channel',arrayfun(@num2str,1:numel(cr),'uni',0));
end
chkeys = regexprep(chnames,'[^a-zA-Z0-9_]','');
assert(numel(chkeys)==numel(unique(chkeys)),...
    'Channel names has duplicates after removing invalid characters');
infsz = size(cellInf.times);
for f=1:numel(extract)
    fname = extract{f};
    if numel(strfind(fname,'.'))==1
        C = strsplit(fname,'.');
        fname = C{2};
        chname = C{1};
        if isempty(regexp(chname,'^\d+$','once'))
            chind = find(strcmp(chnames,chname));
        else
            chind = str2double(chname);
            chname = chnames{chind};
        end
        outname = strjoin([chkeys(chind),C(2)],'_');
    else
        chname = '1';
        chind = 1;
        outname = fname;
    end
    
    if isa(cr,'cellResults')
        if chind==1 && isfield(cr.cellInf.general,fname)
            var = cr.cellInf.general.(fname);
        else
            var = cr.cellInf.(chname).(fname);
        end
    else
        var = cr(chind).(fname);
    end
    if isvector(var) && numel(var)==infsz(1)
        mvar = var(cellfilt);
        dvar = var(dFilt);
    elseif isvector(var) && numel(var)==infsz(2)
        mvar = var(timepoints); % save only for mothers
        dvar = [];
    elseif isequal(size(var),infsz)
        mvar = var(cellfilt,timepoints);
        dvar = var(dFilt,timepoints);
    else
        mvar = var; % save only for mothers
        dvar = [];
    end
    if issparse(mvar)
        mvar = full(mvar); dvar = full(dvar);
        if ~islogical(mvar)
            mvar(mvar==0) = NaN; dvar(dvar==0) = NaN;
        end
    end
    if isinteger(mvar)
        mvar = double(mvar); dvar = double(dvar);
    end
    if ~ismember(fname,notformothers)
        forextraction.M.(outname) = mvar;
    end
    forextraction.D.(outname) = dvar;
end

[ncells,ntps] = size(mVols);
M = struct(); D = cell(ncells,1);
M.times = times;
M.vol = mVols;
M.births = births;
M.daughterLabel = dNum;
M.posNum = cellInf.posNum(cellfilt);
M.trapNum = cellInf.trapNum(cellfilt);
M.cellNum = cellInf.cellNum(cellfilt);
if isfield(cellInf,'mothers')
    M.mothers = cellInf.mothers(cellfilt);
end

% Add any requested fields to extract
fnames = fieldnames(forextraction.M);
for f=1:numel(fnames)
    M.(fnames{f}) = forextraction.M.(fnames{f});
end
fnames = fieldnames(forextraction.D);
for c=1:ncells
    nB = numel(dInds{c});
    btps = find(births(c,:));
    fextra = [fnames,repmat({{}},size(fnames))]';
    D{c} = struct('times',{},'vol',{},'svol',{},'grs',{},'d2v',{},fextra{:});
    for b=1:nB
        if full_daughters
            sT = 1; eT = ntps;
        else
            sT = btps(b); % start tp
            % end tp is either next birth or end of time course if no next
            if b == numel(btps), eT = ntps; else, eT = btps(b+1)-1; end
            D{c}(b).tps = sT:eT;
        end
        
        D{c}(b).times = times(c,sT:eT);
        D{c}(b).vol = dVols(dInds{c}(b),sT:eT);
        
        % Add any requested fields to extract
        for f=1:numel(fnames)
            dvar = forextraction.D.(fnames{f});
            if isempty(dvar), continue; end
            if isvector(dvar)
                D{c}(b).(fnames{f}) = dvar(dInds{c}(b));
            else
                D{c}(b).(fnames{f}) = dvar(dInds{c}(b),sT:eT);
            end
        end
    end
end

% Perform additional filtering if all daughter tracks are present
if full_daughters && usefilters
    % NB: don't need to filter on dVols or dFilt. Can filter dInds instead-
    % this specifies the daughter index into the daughter arrays produced
    % with dFilt (i.e., dVols and extracted dvars below).
    for m=1:numel(D)
        dvars = D{m}; m_di = dInds{m};
        if isempty(dvars) || isempty(fieldnames(dvars)), continue; end
        % Filter out missing tracks
        validtracks = any(~isnan(vertcat(dvars.vol)),2);
        dvars = dvars(validtracks); m_di = m_di(validtracks);
        
        % Filter out small non-growing tracks
        nd = numel(dvars); m_dVols = vertcat(dvars.vol);
        isPresent = ~isnan(m_dVols);
        tp_start = arrayfun(@(d) find(isPresent(d,:),1),1:nd);
        tp_end = arrayfun(@(d) find(isPresent(d,:),1,'last'),1:nd);
        [dvmin,dvmin_ind] = min(m_dVols,[],2,'omitnan');
        [dvmax,dvmax_ind] = max(m_dVols,[],2,'omitnan');
        duration = times(m,tp_end)-times(m,tp_start);
        m_dGRest = sign(dvmax_ind-dvmin_ind).*(dvmax-dvmin)./duration(:);
        smalltracks = max(m_dVols,[],2,'omitnan')<smallVolThresh ...
            & m_dGRest<smallGrowthThresh;
        dvars = dvars(~smalltracks); m_di = m_di(~smalltracks);
        m_dGRest = m_dGRest(~smalltracks,:);
        
        % Join daughter tracks that are contiguous and within a volume 
        % threshold of each other
        nd = numel(dvars); m_dVols = vertcat(dvars.vol);
        isPresent = ~isnan(m_dVols);
        tp_start = arrayfun(@(d) find(isPresent(d,:),1),1:nd);
        tp_end = arrayfun(@(d) find(isPresent(d,:),1,'last'),1:nd);
        [tp_start,t_order] = sort(tp_start); tp_end = tp_end(t_order);
        dvars = dvars(t_order); m_di = m_di(t_order);
        % Determine all possible pairs of valid contiguous tracks
        d_contig = arrayfun(@(s) find(tp_end==s-1),tp_start,'uni',0);
        dc_stats = struct('dn',{},'dp',{},'dv',{}); p = 0;
        for dn=1:numel(d_contig)
            for dp=d_contig{dn}
                % Extrapolate average track growth rates forward and backward
                dt = times(m,tp_start(dn))-times(m,tp_end(dp));
                v_dp = m_dVols(dp,tp_end(dp));
                v_dn = m_dVols(dn,tp_start(dn));
                e_dp = v_dp+m_dGRest(dp)*dt;
                e_dn = v_dn-m_dGRest(dn)*dt;
                dv = min(abs(e_dp-v_dn),abs(e_dn-v_dp));
                if dv<joinVolThresh
                    p = p+1;
                    dc_stats(p).dn = dn;
                    dc_stats(p).dp = dp;
                    dc_stats(p).dv = dv;
                end
            end
        end
        d_joined = [];
        for d=1:numel(dc_stats)
            if numel(dc_stats)==0, break; end
            [~,i] = min([dc_stats.dv]);
            dp = dc_stats(i).dp; dn = dc_stats(i).dn;
            fnames = fieldnames(dvars);
            for f=1:numel(fnames)
                fname = fnames{f};
                if numel(dvars(dn).(fname))~=ntps || ...
                        numel(dvars(dp).(fname))~=ntps
                    continue;
                end
                tps = tp_start(dp):tp_end(dp);
                dvars(dn).(fname)(tps) = dvars(dp).(fname)(tps);
            end
            d_joined(end+1) = dp; %#ok<AGROW>
            dc_stats([dc_stats.dp]==dp | [dc_stats.dn]==dn) = [];
        end
        dvars(d_joined) = []; m_di(d_joined) = [];

        % Discard tracks that are short
        shorttracks = sum(~isnan(vertcat(dvars.vol)),2)<shortTPthresh;
        dvars = dvars(~shorttracks); m_di = m_di(~shorttracks);
        D{m} = dvars;
        
        % Need to regenerate births and daughterLabel for the mothers
        isPresent = ~isnan(vertcat(dvars.vol)); nd=numel(dvars);
        tp_start = arrayfun(@(d) find(isPresent(d,:),1),1:nd);
        tp_end = arrayfun(@(d) find(isPresent(d,:),1,'last'),1:nd);
        M.births(m,:) = false;
        M.births(m,tp_start) = true;
        m_dnum = cellInf.cellNum(dFilt(m_di));
        M.daughterLabel(m,:) = NaN;
        for d=1:nd
            if d==nd, tpe = tp_end(d);
            else, tpe = min(tp_end(d),tp_start(d+1)); end
            M.daughterLabel(m,tp_start(d):tpe) = m_dnum(d);
        end
    end
end

% Calculate growth rates
M.svol = NaN(ncells,ntps);
M.grs = NaN(ncells,ntps);
M.d2v = NaN(ncells,ntps);
grfun = grfuns.(gr_method);
for c=1:ncells
    % Estimate mother growth rate
    Vm = M.vol(c,:);
    if log_gr, Vm = log(Vm); end
    switch birthHandling
        case 'none'
            [M.svol(c,:),M.grs(c,:),M.d2v(c,:)] = grfun(times(c,:),Vm,mOpt);
        case 'split'
            btps = find(M.births(c,:));
            stp = 1;
            for btp=[btps,ntps]
                tps = stp:btp-1;
                [M.svol(c,tps),M.grs(c,tps),M.d2v(c,tps)] = ...
                    grfun(times(c,tps),Vm(tps),mOpt);
                stp = btp;
            end
        case 'cumulative'
            % NB: daughter IDs are found from daughterLabel, so cannot
            % account for two daughters appearing at same time point! I.e.,
            % births is the best reference for birth tps. Furthermore, if
            % full_daughters=false, then the daughter tracks are already
            % truncated to validtps.
            validtps = find(~isnan(Vm),1):find(~isnan(Vm),1,'last');
            btps = []; Vd = [];
            if ~isempty(fieldnames(D{c}))
                if full_daughters
                    btps = arrayfun(@(x) find(~isnan(x.vol),1),D{c});
                    Vd = arrayfun(@(x,btp) x.vol(btp),D{c},btps);
                else
                    btps = find(M.births(c,:));
                    Vd = arrayfun(@(x) x.vol(1),D{c});
                end
            end
            validbtps = ismember(btps,validtps(2:end));
            Vd = Vd(validbtps);
            btps = btps(validbtps);
            % convert btps to validtps indexing
            btps_valid = btps-validtps(1)+1;
            cvol = diff(Vm(validtps));
            if log_gr
                cvol(btps_valid-1) = log(exp(Vm(btps))+Vd)-Vm(btps-1);
            else
                cvol(btps_valid-1) = Vm(btps)+Vd-Vm(btps-1);
            end
            cvol = cumsum([Vm(validtps(1)),cvol],'omitnan');
            if numel(cvol) ~= numel(validtps)
                error('');
            end
            Vm(validtps) = cvol;
            [M.svol(c,:),M.grs(c,:),M.d2v(c,:)] = grfun(times(c,:),Vm,mOpt);
        otherwise
            error('the specified "birthHandling" method is not yet implemented');
    end
    
    % Estimate growth rates for each daughter individually
    if isempty(fieldnames(D{c})), continue; end
    for b=1:numel(D{c})
        Vd = D{c}(b).vol;
        if log_gr, Vd = log(Vd); end
        
        [D{c}(b).svol,D{c}(b).grs,D{c}(b).d2v] = ...
            grfun(D{c}(b).times,Vd,dOpt);
    end
end

end

function [sv,gr,d2v] = grSGolay(T,V,opt)
fP = opt.fP; fW = opt.fW;
sv = NaN(size(V)); gr = NaN(size(V)); d2v = NaN(size(V));
[V,min_tp,max_tp] = getContinuousTrack(V);
ntps = max_tp-min_tp+1;
if ntps<2, return; end
V = V(min_tp:max_tp);
T = T(min_tp:max_tp);
dt = mean(diff(T));
if ntps>=fW
    sv(min_tp:max_tp) = savitzkyGolayFilt(V,fP,0,fW);
    gr(min_tp:max_tp) = savitzkyGolayFilt(V,fP,1,fW)/(-dt);
    d2v(min_tp:max_tp) = savitzkyGolayFilt(V,fP,2,fW)/(dt^2);
else
    % Perform linear regression
    T = [ones(ntps,1),T(:)];
    b = T\V(:);
    sv(min_tp:max_tp) = T*b;
    gr(min_tp:max_tp) = b(2);
    d2v(min_tp:max_tp) = 0;
end
end

function [clean,tp_min,tp_max] = getContinuousTrack(raw,T)
if nargin<2 || isempty(T), T = 1:numel(raw); end

clean = raw; tp_min=1; tp_max=0;

% Find end points of track
mask = ~isnan(raw) & ~isinf(raw); nvalid = sum(mask); 
if nvalid==0, return; end
tp_min = find(mask,1); tp_max = find(mask,1,'last');
if nvalid==1, return; end

% Fill in any interior NaN/Inf values using linear interpolation
clean(tp_min:tp_max) = interp1(T(mask),raw(mask),T(tp_min:tp_max),'linear');
end

function [sv,gr,d2v] = grSLM(T,V,opt)
opt.knotfit = 'fixed';
[sv,gr,d2v] = grSLMfree(T,V,opt);
end

function [sv,gr,d2v] = grSLMfree(T,V,opt)
sW = opt.sW;
if isfield(opt,'knotfit'), kF = opt.knotfit; else, kF = 'free'; end
if isfield(opt,'increasing'), incr = opt.increasing; else, incr = 'off'; end
sv = NaN(size(V)); gr = NaN(size(V)); d2v = NaN(size(V));
mask = ~isnan(V) & ~isinf(V);
if sum(mask)<2, return; end
min_tp = find(mask,1); max_tp = find(mask,1,'last');
nvalid = max_tp-min_tp+1;
nknots = max(round(nvalid/sW),2);
if nknots>2, knotfit = kF; else, knotfit = 'fixed'; end
slm = slmengine(T(mask),V(mask),'increasing',incr,...
    'knots',nknots,'interiorknots',knotfit,...
    'endconditions','natural');
sv(min_tp:max_tp) = slmeval(T(min_tp:max_tp),slm,0);
gr(min_tp:max_tp) = slmeval(T(min_tp:max_tp),slm,1);
d2v(min_tp:max_tp) = slmeval(T(min_tp:max_tp),slm,2);
end