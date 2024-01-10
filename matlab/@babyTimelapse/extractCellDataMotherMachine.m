function extractCellDataMotherMachine(cTimelapse)
% extractCellDataMotherMachine(cTimelapse)
%
% extraction function to extract common cell measurements from experiments
% where cells are vertically stacked in a mother machine well.
% The primary advantage of this extraction method is an enhanced estimation
% of background fluorescence.
%
% Uses the parameters structure
%
%       cTimelapse.extractionParameters.functionParameters
%
% with fields:
% channels                  : cell string specifying channels to extract
% 
% stacks                    : cell array of stack indices to include in
%                             projection
%
% type                      : string max,min,sum or std - specifies how to handle
%                             channels referring to z stacks. Stack is retrieved as
%                             stack, turned into a double, and the following
%                             appropriate statistical procedure applied along the
%                             3rd axis.
%                                    max  - takes maximum z projection
%                                    min  - takes maximum z projection
%                                    sum  - takes the sum of z projection
%                                    std  - takes standard deviation of z
%                                           projection
%
% ffcorr                    : logical array specifying whether to do flat
%                             field correction for corresponding channel
%
% cStats                    : experimentSampleStats object
%
% Calculates a raft of statistics, hopefully mostly self explanatory, about 
% the cell pixels.
%
% size and position information common to all cells is stored in the first
% entry of extractedData only.
%
% ===== STATITICS EXTRACTED =======
%
% The following is a list of the the extracted statistics with some
% description. Those marked (legacy) have ben left for legacy reasons
% should probably be ignored.
%
% To clarify, the following nomenclature will be used:
% cell pixels - the set of all pixels that make up the cell (i.e. the
%               filled in outline of the cell)
% cell membrane pixels - the pixels that make up the outline of the cell.
% reduced cell pixels - all pixels contained in the cell outline after it
%                       has been eroded by 2 pixels.
%
% 'radius' - the approximate radius of the cell calculated as the mean of
%            the radii used to make the active contour outline.
% 'area' - number of pixels which make up the cell.
% 'mean' - the mean of the cell pixels.
% 'median' - median of the cell pixels.
% 'max5px' - the mean of the five brightest cell pixels.
% 'max2p5pc' - the mean of the brightest 2.5% of cell pixels.
% 'bgdvert' - Background calculated as a function of vertical location.
% 'imBackground' - median of all pixels which are not part of any cell.
% 'xloc' - x location in the trap. (i.e.relative to trap centre).
% 'yloc' - y location in the trap. (i.e.relative to trap centre).
% 'trapNum' - list of trap indices for each cell extracted.
% 'cellNum' - list of cell labels for each cell extracted.
% 'extractionParameters' - the parameters passed to this function.
% 'posNum' - (appended by compileCellInformation) position index the
%            position in which each cell was found.
% 'date' - populated by compileCellInformation.
% 'times' - populated by compileCellInformation.
% 'timepointsNotProcessed' - populated by compileCellInformation.
% 'annotations' - populated by compileCellInformation.

parameters = cTimelapse.extractionParameters.functionParameters;

channels = parameters.channels;
assert(iscellstr(channels),'"channels" must be a cell string');
assert(all(ismember(channels,cTimelapse.channelNames)),...
    'some "channels" are not in the cTimelapse');

if ~isfield(parameters,'stacks') || isempty(parameters.stacks)
    stacks = cell(size(channels));
else
    stacks = parameters.stacks;
end
assert(iscell(stacks) && all(cellfun(@(x) isnumeric(x) || isempty(x),stacks)),...
    '"stacks" must be a cell array of indices');
if numel(stacks)==1, stacks = stacks(ones(size(channels))); end
assert(numel(stacks)==numel(channels),'length of "stacks" must match "channels"');

type = parameters.type;
if ischar(type), type = {type}; end
assert(iscellstr(type),'"type" must be a cell string');
if numel(type)==1, type = type(ones(size(channels))); end
assert(numel(type)==numel(channels),'length of "type" must match "channels"');

if ~isfield(parameters,'nerode') || isempty(parameters.nerode)
    parameters.nerode = 2;
end
nerode = parameters.nerode;

if ~isfield(parameters,'ndilate') || isempty(parameters.ndilate)
    parameters.ndilate = 2;
end
ndilate = parameters.ndilate;

hascStats = isfield(parameters,'cStats');
if ~isfield(parameters,'ffcorr') || isempty(parameters.ffcorr)
    if hascStats
        ffcorr = true(size(channels));
    else
        ffcorr = false(size(channels));
    end
else
    ffcorr = parameters.ffcorr;
end
assert(islogical(ffcorr),'"ffcorr" must be a logical vector');
if numel(ffcorr)==1, ffcorr = ffcorr(ones(size(channels))); end
assert(numel(ffcorr)==numel(channels),'length of "ffcorr" must match "channels"');

if ~isfield(parameters,'maskonly') || isempty(parameters.maskonly)
    maskonly = false(size(channels));
else
    maskonly = parameters.maskonly;
end
if any(maskonly & ~ffcorr)
    error('"maskonly" must be specified for channels with ffcorr=true');
end

if any(ffcorr)
    assert(hascStats,'must specify valid "cStats" if any "ffcorr" are true');
    cStats = parameters.cStats;
    ffChans = unique(channels(ffcorr));
    ffIms = struct();
    if isa(cStats,'experimentSampleStats')
        assert(all(ismember(ffChans,cStats.channelNames)),...
            'some requested channels missing in the "cStats" object');
        for c=1:numel(ffChans)
            cStats.channel = ffChans{c};
            ffIms.(ffChans{c}) = struct();
            ffIms.(ffChans{c}).darkfield = cStats.darkfield;
            ffIms.(ffChans{c}).flatfield = cStats.flatfield;
            if any(maskonly)
                ffIms.(ffChans{c}).apermask = cStats.aperture.mask;
            end
        end
    elseif isstruct(cStats)
        assert(all(isfield(cStats,ffChans)),...
            'some requested channels missing in the "cStats" struct');
        assert(all(structfun(@(x) isstruct(x) && ...
            all(isfield(x,{'darkfield','flatfield'})),cStats)),...
            '"cStats" struct has badly formatted channel fields');
        if any(maskonly)
            assert(all(structfun(@(x) isfield(x,'apermask'),cStats)),...
                '"cStats" struct is missing "apermask" fields');
        end
        ffIms = cStats;
    else
        error('unrecognised type of "cStats"');
    end
end

numCells=sum(cTimelapse.cellsToPlot(:));
[trap, cells]=find(cTimelapse.cellsToPlot);

% NB: since cTimelapse.pixelSize is currently tied to rescaling of images,
% this is a temporary fix to allow specification of the true pixel size.
px_size = cTimelapse.pixelSize;
if isfield(parameters,'pixelSize')
    px_size = parameters.pixelSize;
end

% For 3D convolution across stacks, need ratio of Z pixels to X pixels
zx_ratio = [];
if recurse_fields(cTimelapse.metadata,{'acq','zsections','spacing'}) ...
        && isscalar(cTimelapse.metadata.acq.zsections.spacing) ...
        && isscalar(cTimelapse.pixelSize)
    % Express ratio in pixels: number of Z pixels per X pixel
    % Therefore need to express in consistent metric (i.e., pixels/um)
    % So take (1/Z_px_len)/(1/X_px_len):
    zx_ratio = cTimelapse.pixelSize / cTimelapse.metadata.acq.zsections.spacing;
end

% annoying necessary line to deal with find behaviour. If find is applied
% to a row vector, it returns a row vector. If it is applied to a matrix it
% returns a column vector.
trap = trap(:);
cells = cells(:);

%reorder so cells in the same trap contiguous, nicer for viewing later.
[trap,I] = sort(trap);
cells =cells(I);

uniqueChannels = unique(channels);
channelNums = cellfun(@(x) find(strcmp(cTimelapse.channelNames,x),1),uniqueChannels);
timepoints = find(cTimelapse.timepointsProcessed);

% First delete any timer objects previously associated with this function
%oldTimers = timerfindall('TimerFcn',@getImageData);
%delete(oldTimers);
% Should be able to use the above argument, but seems to have broken in
% Matlab 2020a... So search manually:
oldTimers = timerfindall;
ourTimers = arrayfun(@(x) isequal(x.TimerFcn,@getImageData),oldTimers);
delete(oldTimers(ourTimers));

% Create a timer object to asynchronously load image data
loadImageTimer = timer('TimerFcn',@getImageData,'StartDelay',0.001);
timerdata = struct();
timerdata.channelNums = channelNums;
timerdata.timepointStacks = cell(size(channelNums));
timerdata.timepoint = timepoints(1);
loadImageTimer.UserData = timerdata;
start(loadImageTimer); % trigger loading of first set of images

imQs = (0:0.025:1)';

%preallocate cellInf
extractedData = struct();
ntps = numel(cTimelapse.timepointsProcessed);
ntraps = numel(cTimelapse.cTimepoint(timepoints(1)).trapInfo);
for chi=1:numel(channels)
    extractedData(chi).imQuantiles=zeros(numel(imQs),ntps);
    extractedData(chi).nonTrapBackground=zeros(1,ntps);
    
    % Trap-level parameters
    extractedData(chi).trapX=sparse(zeros(ntraps,ntps));
    extractedData(chi).trapY=sparse(zeros(ntraps,ntps));
    extractedData(chi).trapNCells=sparse(zeros(ntraps,ntps));
    extractedData(chi).trapBgdMean=sparse(zeros(ntraps,ntps));
    extractedData(chi).trapBgdMedian=sparse(zeros(ntraps,ntps));
    extractedData(chi).trapBgdArea=sparse(zeros(ntraps,ntps));
    
    extractedData(chi).radius=sparse(zeros(numCells,ntps));
    extractedData(chi).area=sparse(zeros(numCells,ntps));
    
    extractedData(chi).mean=sparse(zeros(numCells,ntps));
    extractedData(chi).median=sparse(zeros(numCells,ntps));
    extractedData(chi).innerMean=sparse(zeros(numCells,ntps));
    extractedData(chi).innerMedian=sparse(zeros(numCells,ntps));
    
    extractedData(chi).meanBlur=sparse(zeros(numCells,ntps));
    extractedData(chi).innerMeanBlur=sparse(zeros(numCells,ntps));
    
    extractedData(chi).meanLog=sparse(zeros(numCells,ntps));
    extractedData(chi).meanLogBgdSub=sparse(zeros(numCells,ntps));
    extractedData(chi).meanLogBlurBgdSub=sparse(zeros(numCells,ntps));
    if isfield(parameters,'experiment_background')
        extractedData(chi).meanLogBlurExptBgd=sparse(zeros(numCells,ntps));
    end
    if all(isfield(parameters,{'xyshift_model','xyshift_coefs'}))
        extractedData(chi).meanShift=sparse(zeros(numCells,ntps));
        extractedData(chi).meanLogShiftBgdSub=sparse(zeros(numCells,ntps));
        extractedData(chi).meanBlurShift=sparse(zeros(numCells,ntps));
        extractedData(chi).meanLogBlurShiftBgdSub=sparse(zeros(numCells,ntps));
        if isfield(parameters,'experiment_background')
            extractedData(chi).meanLogBlurShiftExptBgd=sparse(zeros(numCells,ntps));
        end
    end
    
    extractedData(chi).bgdvert=sparse(zeros(numCells,ntps));
    extractedData(chi).imBackground=sparse(zeros(numCells,ntps));
    
    extractedData(chi).xloc=sparse(zeros(numCells,ntps));
    extractedData(chi).yloc=sparse(zeros(numCells,ntps));
    
    extractedData(chi).trapNum = trap';
    extractedData(chi).cellNum = cells';
end

se_erode=strel('disk',nerode);
se_dilate=strel('disk',ndilate);

bbH = cTimelapse.cTrapSize.bb_height;
bbW = cTimelapse.cTrapSize.bb_width;

% rearranged so that each cell is filled in at the same time at each
% timepoint, this allows parfor use and reduces the extraction time
% significantly
for t=1:numel(timepoints)
    timepoint = timepoints(t);
    
    % Trigger the TimepointChanged event for babyLogging
    babyLogging.changeTimepoint(cTimelapse,t);
    % fprintf('.');
    % if mod(t,50)==0, fprintf('\n'); end
    
    % Retrieve trap stack images from timer
    wait(loadImageTimer);
    timerdata = loadImageTimer.UserData;
    rawTimepointStacks = timerdata.timepointStacks;
    if t<numel(timepoints)
        % Start downloading images for next time point
        timerdata.timepoint = timepoints(t+1);
        loadImageTimer.UserData = timerdata;
        start(loadImageTimer);
    end
            
    trapInfo = cTimelapse.cTimepoint(timepoint).trapInfo;
    
    cellLocAll = false([cTimelapse.trapImSize,numel(trap)]);
    cellLocAllInner = false([cTimelapse.trapImSize,numel(trap)]);
    for c=1:numel(trap)
        cellInd = find(trapInfo(trap(c)).cellLabel==cells(c));
        if ~isempty(cellInd)
            loc = double(trapInfo(trap(c)).cell(cellInd).cellCenter);
            if ~isempty(loc)
                cellLoc = full(trapInfo(trap(c)).cell(cellInd).segmented);
                cellLoc = imfill(cellLoc,'holes');
                cellLocAll(:,:,c) = cellLoc;
                cellLocAllInner(:,:,c) = imerode(cellLoc,se_erode);
            end
        end
    end
    
    trapNums = unique(trap);
    cellLocTrapCellsBkg = false([cTimelapse.trapImSize,numel(trapNums)]);
    cellLocTrapCellsBkgOuter = false([cTimelapse.trapImSize,numel(trapNums)]);
    for tn=1:numel(trapNums)
        trapNum = trapNums(tn);
        if trapInfo(trapNum).cellsPresent
            tSeg = cellLocTrapCellsBkg(:,:,tn);
            for cellInd=1:numel(trapInfo(trapNum).cellLabel)
                tSeg = tSeg | full(trapInfo(trapNum).cell(cellInd).segmented);
            end
            tSeg = imfill(tSeg,'holes');
            cellLocTrapCellsBkg(:,:,tn) = tSeg;
            cellLocTrapCellsBkgOuter(:,:,tn) = imdilate(tSeg,se_dilate);
        end
    end
    cellLocAllCellsBkg = any(cellLocTrapCellsBkg,3);
    
    trapInfo=cTimelapse.cTimepoint(timepoint).trapInfo;
    for chi=1:numel(channels) % channel index
        chname = channels{chi};
        rawchi = find(strcmp(uniqueChannels,chname),1);
        
        % Get the correct channel images
        tpStack = rawTimepointStacks{rawchi};
        channelNum = channelNums(rawchi);
        
        %if empty do nothing
        if all(tpStack(:)==0), continue; end
        
        % Apply flat field correction if requested
        if ffcorr(chi) && ismember(chname,ffChans)
            ff = ffIms.(chname);
            if maskonly(chi)
                if size(ff.apermask,3)==1
                    apermask = ff.apermask(:,:,ones(size(tpStack,3),1));
                end
                apermask = flipFFimg(apermask,channelNum);
                tpStack = double(tpStack);
                tpStack(~apermask) = NaN;
            else
                if size(ff.darkfield,3)==1
                    darkfield = ff.darkfield(:,:,ones(size(tpStack,3),1));
                else
                    darkfield = ff.darkfield;
                end
                if size(ff.flatfield,3)==1
                    flatfield = ff.flatfield(:,:,ones(size(tpStack,3),1));
                else
                    flatfield = ff.flatfield;
                end
                darkfield = flipFFimg(darkfield,channelNum);
                flatfield = flipFFimg(flatfield,channelNum);
                tpStack = (double(tpStack)-double(darkfield))./double(flatfield);
            end
        end
        
        tpStack = cTimelapse.applyStandardImageTransformation(tpStack,channelNum);
        % Ensure we are working with double
        tpStack = double(tpStack);
        
        % Select a subset of the stacks if specified
        zInds = stacks{chi};
        if isempty(zInds), zInds = 1:size(tpStack,3); end
        tpStack = tpStack(:,:,zInds);
        
        % Apply Z-projection method as specified
        switch type{chi}
            case 'max'
                tpImg = max(tpStack,[],3);
            case 'min'
                tpImg = min(tpStack,[],3);
            case 'mean'
                tpImg = mean(tpStack,3);
            case 'std'
                tpImg = std(tpStack,[],3);
            case 'sum'
                tpImg = sum(tpStack,3);
        end
        
        % Save whole image quantiles
        tpImg_blur = imgaussfilt(tpImg,2);
        extractedData(chi).imQuantiles(:,timepoint) = ...
            quantile(tpImg_blur(:),imQs);
        
        % Calculate background excluding all selected traps
        tpBgdImg = tpImg;
        trapLocations = cTimelapse.cTimepoint(timepoint).trapLocations;
        for ti=1:numel(trapLocations)
            trap_y = round(trapLocations(ti).ycenter);
            trap_x = round(trapLocations(ti).xcenter);
            yL = max(trap_y-bbH,1); yU = min(trap_y+bbH,size(tpBgdImg,1));
            xL = max(trap_x-bbW,1); xU = min(trap_x+bbW,size(tpBgdImg,2));
            tpBgdImg(yL:yU,xL:xU) = NaN;
        end
        extractedData(chi).nonTrapBackground(timepoint) = nanmedian(tpBgdImg(:));
        
        % Split into traps
        trapStack = cTimelapse.returnTrapsFromImage(tpImg,timepoint);
        assert(size(trapStack,3) == ntraps,...
            'Number of traps cannot change between time points!');
        
        % Blurred image
        tpImg_blur = imgaussfilt(tpImg,1);
        trapStack_blur = cTimelapse.returnTrapsFromImage(tpImg_blur,timepoint);
        
        % Extract trap-level parameters
        for ti=1:ntraps
            if chi==1
                trapLoc = cTimelapse.cTimepoint(timepoint).trapLocations(ti);
                extractedData(chi).trapX(ti,timepoint)=trapLoc.xcenter;
                extractedData(chi).trapY(ti,timepoint)=trapLoc.ycenter;
                if trapInfo(ti).cellsPresent
                    if numel(trapInfo(ti).cellLabel)~=numel(trapInfo(ti).cell)
                        warning('cell labels do not match cell outlines at time point %u, trap %u',...
                            timepoint,ti);
                    end
                    extractedData(chi).trapNCells(ti,timepoint)=...
                        numel(trapInfo(ti).cellLabel);
                end
            end
            trapImage = trapStack(:,:,ti);
            trapAllCellBkgFL = trapImage(cellLocAllCellsBkg);
            extractedData(chi).trapBgdMean(ti,timepoint)=mean(trapAllCellBkgFL);
            extractedData(chi).trapBgdMedian(ti,timepoint)=median(trapAllCellBkgFL);
            extractedData(chi).trapBgdArea(ti,timepoint)=numel(trapAllCellBkgFL);
        end
        
        for c=1:length(trap)
            % Skip this cell if not present at this time point
            cellInd = find(trapInfo(trap(c)).cellLabel==cells(c));
            if isempty(cellInd), continue; end
            
            % Skip this cell if segmentation mask is empty
            cellLoc = cellLocAll(:,:,c)>0;
            if ~sum(cellLoc(:)), continue; end
            
            trapImage = trapStack(:,:,trap(c));
            trapImage_blur = trapStack_blur(:,:,trap(c));
            
            %vector of cell pixels
            cellFL=trapImage(cellLoc);
            nCellFL = numel(cellFL);
            
            extractedData(chi).mean(c,timepoint)=mean(cellFL(:));
            extractedData(chi).median(c,timepoint)=median(cellFL(:));
            extractedData(chi).meanLog(c,timepoint)=mean(log(cellFL(:)));
            
            cellInnerFL=trapImage(logical(cellLocAllInner(:,:,c)));
            extractedData(chi).innerMean(c,timepoint)=mean(cellInnerFL(:));
            extractedData(chi).innerMedian(c,timepoint)=median(cellInnerFL(:));
            
            cellFL_blur = trapImage_blur(cellLoc);
            extractedData(chi).meanBlur(c,timepoint)=mean(cellFL_blur(:));
            cellInnerFL_blur = trapImage_blur(logical(cellLocAllInner(:,:,c)));
            extractedData(chi).innerMeanBlur(c,timepoint)=mean(cellInnerFL_blur(:));
            
            bgd = nanmedian(tpBgdImg(:));
            extractedData(chi).meanLogBgdSub(c,timepoint) = ...
                mean(log(cellFL(cellFL>bgd)-bgd));
            extractedData(chi).meanLogBlurBgdSub(c,timepoint) = ...
                mean(log(cellFL_blur(cellFL_blur>bgd)-bgd));
            if isfield(parameters,'experiment_background')
                ebgd = parameters.experiment_background(chi);
                extractedData(chi).meanLogBlurExptBgd(c,timepoint) = ...
                    mean(log(cellFL_blur(cellFL_blur>ebgd)-ebgd));
            end
            
            if all(isfield(parameters,{'xyshift_coefs','xyshift_model'}))
                trapLoc = cTimelapse.cTimepoint(timepoint).trapLocations(trap(c));
                cc = trapInfo(trap(c)).cell(cellInd).cellCenter;
                trapOffset = [trapLoc.ycenter-bbH-1,trapLoc.xcenter-bbW-1];
                cc_pos = trapOffset+flip(cc(:)');
                designmat = x2fx(cc_pos,parameters.xyshift_model);
                xyshift = NaN(1,2);
                for target=1:2
                    xyshift(target) = designmat * parameters.xyshift_coefs{chi,target};
                end
                pad = ceil(max(abs(xyshift)))+3;
                [~,xL] = max(cellLoc,[],1); xL = min(xL(any(cellLoc,1)));
                [~,xU] = max(flipud(cellLoc),[],1);
                xU = size(cellLoc,1)-min(xU(any(cellLoc,1)))+1;
                [~,yL] = max(cellLoc,[],2); yL = min(yL(any(cellLoc,2)));
                [~,yU] = max(fliplr(cellLoc),[],2);
                yU = size(cellLoc,2)-min(yU(any(cellLoc,2)))+1;
                trapOffset = round(trapOffset);
                xL_pos = max(xL+trapOffset(1)-pad,1);
                xU_pos = min(xU+trapOffset(1)+pad,size(tpImg_blur,1));
                yL_pos = max(yL+trapOffset(2)-pad,1);
                yU_pos = min(yU+trapOffset(2)+pad,size(tpImg_blur,2));
                
                % Shift raw image
                shiftimg = imtranslate(tpImg(xL_pos:xU_pos,yL_pos:yU_pos),xyshift);
                shiftimg = shiftimg(xL+trapOffset(1)-xL_pos+(1:xU-xL+1),...
                    yL+trapOffset(2)-yL_pos+(1:yU-yL+1));
                cellFL_shift = shiftimg(cellLoc(xL:xU,yL:yU));
                extractedData(chi).meanShift(c,timepoint) = ...
                    mean(cellFL_shift(:));
                extractedData(chi).meanLogShiftBgdSub(c,timepoint) = ...
                    mean(log(cellFL_shift(cellFL_shift>bgd)-bgd));
                
                % Shift blurred image
                shiftimg = imtranslate(tpImg_blur(xL_pos:xU_pos,yL_pos:yU_pos),xyshift);
                shiftimg = shiftimg(xL+trapOffset(1)-xL_pos+(1:xU-xL+1),...
                    yL+trapOffset(2)-yL_pos+(1:yU-yL+1));
                cellFL_shift = shiftimg(cellLoc(xL:xU,yL:yU));
                extractedData(chi).meanBlurShift(c,timepoint) = ...
                    mean(cellFL_shift(:));
                extractedData(chi).meanLogBlurShiftBgdSub(c,timepoint) = ...
                    mean(log(cellFL_shift(cellFL_shift>bgd)-bgd));
                if isfield(parameters,'experiment_background')
                    ebgd = parameters.experiment_background(chi);
                    extractedData(chi).meanLogBlurShiftExptBgd(c,timepoint) = ...
                        mean(log(cellFL_shift(cellFL_shift>ebgd)-ebgd));
                end
            end
            
            
%             if strcmp(chname,'CFP'), bgd = 110; end
%             if strcmp(chname,'YFP'), bgd = 117; end
%             cellFL_blur_bgdsub = cellFL_blur - bgd;
%             valid = cellFL_blur_bgdsub > 0;
%             if any(valid)
%                 extractedData(chi).meanLogBlurBgdSub(c,timepoint) = ...
%                     mean(log(cellFL_blur_bgdsub(valid)));
%             end
            
            bgdLoc = ~cellLocTrapCellsBkgOuter(:,:,find(trapNums==trap(c),1));
            bgdImg = trapImage;
            bgdImg(~bgdLoc) = NaN;
            vbgd = medfilt1(median(bgdImg,2,'omitnan'),11,'omitnan','truncate');
            xyloc = trapInfo(trap(c)).cell(cellInd).cellCenter;
            extractedData(chi).bgdvert(c,timepoint)=vbgd(round(xyloc(2)));
            
            bgdFL = trapImage(bgdLoc);
            bgdFL = bgdFL(~isnan(bgdFL(:)));
            if isempty(bgdFL), bgdFL = trapImage; end
            extractedData(chi).imBackground(c,timepoint)=median(bgdFL(:));
            
            % information common to all channels (basically
            % shape information) is stored only in the
            % channel 1 structure.
            if chi==1
                extractedData(chi).area(c,timepoint)=nCellFL;
                extractedData(chi).radius(c,timepoint)=trapInfo(trap(c)).cell(cellInd).cellRadius;
                xyloc = trapInfo(trap(c)).cell(cellInd).cellCenter;
                extractedData(chi).xloc(c,timepoint)=xyloc(1);
                extractedData(chi).yloc(c,timepoint)=xyloc(2);
            end
        end
    end
end

% Clean up timer since it persists beyond function termination
delete(loadImageTimer);

cTimelapse.extractedData = extractedData;

fprintf('\n');

    function getImageData(timerobj,~)
        tdata = timerobj.UserData;
        chNums = tdata.channelNums;
        tp = tdata.timepoint;
        for cn=1:numel(chNums)
            chnum = chNums(cn);
            tdata.timepointStacks{cn} = ...
                cTimelapse.returnSingleTimepointRaw(tp,chnum);
        end
        timerobj.UserData = tdata;
    end

    function img = flipFFimg(img,chnum)
        % Since the FF is determined directly from Omero, we need to
        % apply the same transformations to any FF images that are applied
        % by the returnSingleTimepointRaw function
        
        %Images are flipped in both directions on Omero upload - so if the
        %data was segmented from a folder before being uploaded/converted
        %then it should be flipped to ensure data is consistent with any
        %segmentation results
        if strcmp(cTimelapse.segmentationSource,'Folder')
            img=rot90(img);
            img=flipud(img);
        end
        
        % A bug in the microscope code meant that some channels were saved
        % flipped relative to brightfield. The flipchannels property can be
        % used to specify which channels should be flipped to match the
        % orientation of the brightfield image.
        if ~isempty(cTimelapse.flipchannels) && ...
                length(cTimelapse.flipchannels)>chnum && ...
                cTimelapse.flipchannels(chnum)
            img=flipud(img);
        end
    end
end

%%
    function iscomplete = recurse_fields(x,f)
        iscomplete = isempty(f) || ...
            (isfield(x,f{1}) && recurse_fields(x.(f{1}),f(2:end)));
    end
