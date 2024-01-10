function extractedData = fastGRextraction(cTimelapse,varargin)
%FASTGREXTRACTION Estimate growth rates directly from BABY output
%
%   E = cTimelapse.FASTGREXTRACTION returns the struct E containing growth 
%   data for a selection of mothers and corresponding daughters at each
%   processed time point of the cTimelapse. For NCELL mothers and NTPS time
%   points, the fields in E are:
%       - 'trapNum' the trap index for each mother (vector of length NCELL)
%       - 'cellNum' the within-trap cell label for each mother (vector of
%         length NCELL)
%       - 'daughterLabel' the within-trap cell label for the current 
%         daughter of each mother at each time point (NCELL x NTPS matrix)
%       - 'volume' the estimated volume of each mother at each time point
%         (NCELL x NTPS matrix; in um^3)
%       - 'daughterVolume' the estimated volume of the mother's current 
%         daughter at each time point (NCELL x NTPS matrix; in um^3)
%       - 'GR' the estimated growth rate (in um^3/hour) of the 
%         mother volume (NCELL x NTPS matrix)
%       - 'daughterGR' the estimated growth rate (in um^3/hour) of the 
%         daughter volumes (NCELL x NTPS matrix)
%       - 'logGR' the specific growth rate of the mothers (in hour^{-1}) 
%         (NCELL x NTPS matrix)
%       - 'daughterLogGR' the specific growth rate of the daughters 
%         (in hour^{-1}) (NCELL x NTPS matrix)
%       - 'births' a logical NCELL x NTPS matrix marking the time points at
%          which a new bud is first detected
%   
%   NB: unlike cTimelapse.extractedData, none of the above are sparse and
%   any time points at which a cell is not present or for which a value is 
%   unavailable are left as NaN.
%
%   E = cTimelapse.FASTGREXTRACTION(TPS) limits the extraction to the time 
%   points TPS specified as indices into cTimelapse.cTimepoint.
%
%   E = cTimelapse.FASTGREXTRACTION(...,'selTimepoints',SELTPS) specifies
%   the time-points for which auto-selection metrics will apply (i.e., F).
%   
%   E = cTimelapse.FASTGREXTRACTION(...,'fraction',F) specifies the
%   fraction F of time points over which a mother must have an assigned 
%   daughter in order to be selected (default F=0.8).
%
%   E = cTimelapse.FASTGREXTRACTION(...,'daughter_fraction',DF) specifies
%   the fraction DF of a mothers time-lapse over which daughters must be
%   assigned.
% 
%   E = cTimelapse.FASTGREXTRACTION(...,'volEst',V) specifies the volume
%   estimation method V as one of:
%       - 'baby' the volume estimate returned by BABY (at the time of 
%          writing, a conical estimation method; the default);
%       - 'ellipsoid' uses the major and minor axes to calculate the volume
%          of a prolate spheroid;
%       - 'area' uses the mask area to estimate a radius and calculate the
%          volume of a sphere;
%
%   -- PARAMETERS PASSED TO GRSFROMCELLINF --
%
%   E = cTimelapse.FASTGREXTRACTION(...,'gr_method',GRMETHOD) uses the 
%   method GRMETHOD for smoothing and estimating the growth rate. Available
%   methods are 'sgolay' (default), 'slmfixed' and 'slmfree'.
%
%   E = cTimelapse.FASTGREXTRACTION(...,'log_gr',true) specifies that 
%   growth rate estimation should be performed on the log-volume, thus 
%   obtaining specific growth rates.
%
%   E = cTimelapse.FASTGREXTRACTION(...,'dW',DW,'mW',MW) specify the 
%   smoothing windows DW for daughter cells or MW for mother cells.
%   window DW of the Savitzky-Golay smoothing filter that gets used to 
%   determine growth rates for daughter cells. 'mW' is the smoothing
%   window for mother growth rates. Defaults are DW=11, MW=11.
%

ip = inputParser;
isindex = @(x) isvector(x) && isnumeric(x) && all(x>0) && all(round(x)==x);
ip.addOptional('timepoints',[],@(x) isempty(x) || isindex(x));
ip.addParameter('selTimepoints',[],@(x) isempty(x) || isindex(x));
ip.addParameter('fraction',0.8,@(x) isscalar(x) && isnumeric(x) && x>0 && x<=1);
ip.addParameter('daughter_fraction',0.8,@(x) isscalar(x) && isnumeric(x) && x>0 && x<=1);
ip.addParameter('volEst','baby',@(x) ismember(x,{'baby','ellipsoid','area'}));
ip.addParameter('gr_method','sgolay');
ip.addParameter('pixel_size',0.182,@(x) isempty(x) || (isscalar(x) && isnumeric(x)));
ip.KeepUnmatched = true;
ip.parse(varargin{:});

timepoints = ip.Results.timepoints;
seltps = ip.Results.selTimepoints;
fraction = ip.Results.fraction;
daughter_fraction = ip.Results.daughter_fraction;
volEst = ip.Results.volEst;
gr_method = ip.Results.gr_method;
pixel_size = ip.Results.pixel_size;
unmatched = [fieldnames(ip.Unmatched),struct2cell(ip.Unmatched)]';

if isempty(timepoints)
    timepoints = find(cTimelapse.timepointsProcessed);
end
if isempty(seltps), seltps = timepoints; end
% selTimepoints must be a subset of timepoints
seltps = intersect(seltps,timepoints);
tps_map = NaN(1,max(timepoints));
tps_map(timepoints) = 1:numel(timepoints);
seltpi = tps_map(seltps);

switch volEst
    case 'baby'
        volFun = @getBabyVolume;
    case 'ellipsoid'
        volFun = @getEllipsoidVolume;
    case 'area'
        volFun = @getAreaVolume;
end

% Collect all cell labels and autoselect
ntraps = numel(cTimelapse.cTimepoint(timepoints(1)).trapInfo);
cellLabels = cell(1,ntraps);
selCellLabels = cell(1,ntraps);
for t=1:ntraps
    all_cLbls = arrayfun(@(x) double(x.trapInfo(t).cellLabel),...
        cTimelapse.cTimepoint(timepoints),'UniformOutput',false);
    % For selection, limit collected cell labels to just those present
    % during the time points of interest (i.e., selTimepoints):
    cLbls = horzcat(all_cLbls{seltpi}); cLbls = cLbls(cLbls~=0);
    [cLbls,~,cLbls_index] = unique(cLbls);
    % Save all unique cell labels for later reference
    all_cLbls = horzcat(all_cLbls{:}); all_cLbls = all_cLbls(all_cLbls~=0);
    all_cLbls = unique(all_cLbls); cellLabels{t} = all_cLbls;
    % Select mother cells according to fraction of time points present in 
    % the selTimepoints array:
    mothers = cTimelapse.cellMothers(t,:);
    mLabels = intersect(full(unique(mothers(mothers>0))),all_cLbls);
    cLbl_counts = accumarray(cLbls_index,1);
    selMothers = cLbls(ismember(cLbls(:),mLabels) & ...
        (cLbl_counts(:)/numel(seltps))>=fraction);
    % Select all descendents of the selected mothers
    selLineage = selMothers;
    while true
        selDaughters = intersect(find(ismember(mothers,selLineage)),all_cLbls);
        if isempty(selDaughters) || all(ismember(selDaughters,selLineage))
            break
        end
        selLineage = union(selLineage,selDaughters);
    end
    selCellLabels{t} = selLineage;
end

% Collect cell IDs for selected cells
cellNums = [selCellLabels{:}];
trapNums = repelem(1:ntraps,cellfun(@numel,selCellLabels));

ncells = numel(trapNums);
ntps = numel(timepoints);
if isempty(pixel_size), pixel_size = cTimelapse.pixelSize; end
dt = cTimelapse.metadata.acq.times.interval/60; % in minutes
cTimes = (0:cTimelapse.metadata.acq.times.ntimepoints-1)*dt;
cTimes = cTimes(timepoints); % subset to timepoints in case of seg errors

% Initialise storage arrays
daughterNum = nan(ncells,ntps);
births = false(ncells,ntps);
vols = nan(ncells,ntps);

traps = unique(trapNums);
for trap=traps(:)'
    inds = find(trapNums==trap);
    sel = cellNums(inds);
    all_cLbls = cellLabels{trap};
    indsMap = NaN(max(all_cLbls),1); indsMap(sel) = inds;
    for tt=1:ntps
        trapInfo = cTimelapse.cTimepoint(timepoints(tt)).trapInfo(trap);
        clabs = trapInfo.cellLabel;
        clab_inds = indsMap(clabs);
        valid = find(~isnan(clab_inds));
        for i=1:numel(valid)
            ci = valid(i);
            vols(clab_inds(ci),tt) = volFun(trapInfo,ci);
        end
    end
    % Collate births data for any mothers in sel
    mothers = cTimelapse.cellMothers(trap,:);
    mLabels = sel(ismember(sel,mothers));
    mInds = indsMap(mLabels);
    for m=1:numel(mLabels)
        ind = mInds(m);
        
        % Find all daughters and the tps for which they are present
        dLabels = find(mothers==mLabels(m));
        dLabels = intersect(dLabels,all_cLbls); % ensure only valid dLabels
        dIsPresent = ~isnan(vols(indsMap(dLabels),:));
        dfilt = any(dIsPresent,2);
        dLabels = dLabels(dfilt);
        dIsPresent = dIsPresent(dfilt,:);
        nd = numel(dLabels);
            
        % Flatten the daughter metrics
        starttp = arrayfun(@(d) find(dIsPresent(d,:),1),1:nd);
        % Order first by start tp then by length (for any coincident)
        [~,dOrder] = sortrows([starttp(:),sum(dIsPresent,2)]);
        dIsPresent = dIsPresent(dOrder,:);
        dLabels = dLabels(dOrder);
        starttp_prev = 0;
        for d=1:nd
            valid = dIsPresent(d,:);
            % The following handles the case of coincident daughters:
            valid(1:starttp_prev) = false;
            if ~any(valid)
                warning('unresolved coincident daughters!');
                continue
            end
            starttp(d) = find(valid,1); starttp_prev = starttp(d);
            daughterNum(ind,valid) = dLabels(d);
        end
        births(ind,starttp) = true;
    end
end

inputData = struct();
inputData.cellNum = cellNums;
inputData.trapNum = trapNums;
inputData.times = cTimes(ones(size(cellNums)),:);
inputData.daughterLabel = daughterNum;
inputData.births = births;
inputData.volume = vols;

% Filter gr estimates just to those mother cells that are present for the 
% specified fraction of the time course (daughters are collated
% independently):
cellfilt = any(births,2) ...
    & sum(~isnan(vols(:,seltpi)),2)/numel(seltpi)>=fraction;
expt = grData(inputData,...
    'cellfilt',cellfilt,'vol_type','conical','gr_method',gr_method,...
    'full_daughters',true,'pixel_size',pixel_size,unmatched{:});

extractedData = expt.M;
extractedData.daughterVolume = expt.flatten_dVols;
extractedData.daughterGR = expt.flatten_dGRs;

% Final autoselect after births filtering performed by grsFromCellInf:
mIsPresent = ~isnan(extractedData.vol);
dIsPresent = ~isnan(extractedData.daughterLabel);
cellfilt = any(extractedData.births,2) ...
    & sum(dIsPresent,2)./sum(mIsPresent,2)>=daughter_fraction;
fnames = fieldnames(extractedData);
for f=1:numel(fnames)
    fname = fnames{f};
    if ismember(fname,{'posNum','trapNum','cellNum'})
        extractedData.(fname) = extractedData.(fname)(cellfilt);
    else
        extractedData.(fname) = extractedData.(fname)(cellfilt,:);
    end
end
end

function vol = getBabyVolume(trapInfo,ci)
vol = NaN;
if isfield(trapInfo,'cellVolumes') && ~isempty(trapInfo.cellVolumes)
    vol = trapInfo.cellVolumes(ci);
end
end

function vol = getEllipsoidVolume(trapInfo,ci)
vol = NaN;
if isfield(trapInfo,'ellipseDims') && ~isempty(trapInfo.ellipseDims)
    vol = pi/6*prod(trapInfo.ellipseDims(ci,[1,2,2]));
end
end

function vol = getAreaVolume(trapInfo,ci)
mask = imfill(full(trapInfo.cell(ci).segmented),'holes');
vol = 4/3*pi*sqrt(sum(mask(:))/pi).^3;
end