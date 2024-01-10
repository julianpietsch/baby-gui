function extractCellDataNucLocVariants(cTimelapse)
% extractCellDataNucLocVariants(cTimelapse)
%
% Extracts for variants of the original nucConvEst algorithm
% Also estimates size of gaussian filter that should be used
% Based on extractCellDataStacksParfor
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
% The standard data extraction method. Calculates a raft of statistics,
% hopefully mostly self explanatory, about the cell pixels. Also allows for
% nuclear localisation measurement by specifying a nuclearMarkerChannel :
% a channel of the timelapse which is a nuclear marker. The max z
% projection of this channel is then used to define nuclear pixels and the
% mean of the maxAllowedOverlap brightest pixels in this nuclear channel is
% used as a candidate set of nuclear pixels. The maxPixOverlap brightest
% pixels in the measurement channel of these candidate pixels are
% determined to be the nuclear ones and the mean of these divided by the
% median of the whole cell are stored as nuclearTagLoc field of
% extractedData.
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
        && isscalar(px_size)
    % Express ratio in pixels: number of Z pixels per X pixel
    % Therefore need to express in consistent metric (i.e., pixels/um)
    % So take (1/Z_px_len)/(1/X_px_len):
    zx_ratio = cTimelapse.pixelSize / cTimelapse.metadata.acq.zsections.spacing;
end

gaussianLocParams = struct('mp',struct(),'dw',struct());

gaussianLocParams.mp.optBgdQ = 0.24;
gaussianLocParams.mp.alpha = 0.85;
gaussianLocParams.mp.nucFrac = 0.34;
gaussianLocParams.mp.nucArea = 4.3; % µm^2
if isfield(parameters,'gaussianLocParams')
    fnames = fieldnames(parameters.gaussianLocParams);
    for f=1:numel(fnames)
        gaussianLocParams.mp.(fnames{f}) = parameters.gaussianLocParams.(fnames{f});
    end
end

gaussianLocParams.dw.optBgdQ = 0.24;
gaussianLocParams.dw.alpha = 0.85;
gaussianLocParams.dw.nucFrac = 0.28;
gaussianLocParams.dw.nucArea = 4.0; % µm^2
if isfield(parameters,'gaussianLocDWParams')
    fnames = fieldnames(parameters.gaussianLocDWParams);
    for f=1:numel(fnames)
        gaussianLocParams.dw.(fnames{f}) = parameters.gaussianLocDWParams.(fnames{f});
    end
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
trapHw = cTimelapse.cTrapSize.bb_width;
trapHh = cTimelapse.cTrapSize.bb_height;

% First delete any timer objects previously associated with this function
oldTimers = timerfindall('TimerFcn',@getImageData);
delete(oldTimers);

% Create a timer object to asynchronously load image data
loadImageTimer = timer('TimerFcn',@getImageData,'StartDelay',0.001);
loadImageTimer.UserData = struct();
loadImageTimer.UserData.channelNums = channelNums;
loadImageTimer.UserData.timepointStacks = cell(size(channelNums));
loadImageTimer.UserData.timepoint = timepoints(1);
start(loadImageTimer); % trigger loading of first set of images

imQs = (0:0.025:1)';

%preallocate cellInf
extractedData = struct();
ntps = numel(cTimelapse.timepointsProcessed);
for chi=1:numel(channels)
    extractedData(chi).imQuantiles=zeros(numel(imQs),ntps);
    extractedData(chi).radius=sparse(zeros(numCells,ntps));
    extractedData(chi).area=sparse(zeros(numCells,ntps));
    
    extractedData(chi).mean=sparse(zeros(numCells,ntps));
    extractedData(chi).median=sparse(zeros(numCells,ntps));
    extractedData(chi).max5px=sparse(zeros(numCells,ntps));
    extractedData(chi).max2p5pc=sparse(zeros(numCells,ntps));
    extractedData(chi).imBackground=sparse(zeros(numCells,ntps));
    
    extractedData(chi).smallPeakConv=sparse(zeros(numCells,ntps));
    extractedData(chi).nucEstConv=sparse(zeros(numCells,ntps));
    extractedData(chi).nucEstConvDW=sparse(zeros(numCells,ntps));
    extractedData(chi).nucEstConv3=sparse(zeros(numCells,ntps));
    
    extractedData(chi).gaussianLoc = sparse(zeros(numCells,ntps));
    extractedData(chi).gaussianLocUnscaled = sparse(zeros(numCells,ntps));
    extractedData(chi).gaussianLocDW = sparse(zeros(numCells,ntps));
    extractedData(chi).gaussianLocDWUnscaled = sparse(zeros(numCells,ntps));
    
    extractedData(chi).gaussianLocBgd = sparse(zeros(numCells,ntps));
    extractedData(chi).gaussianLocBgdQDev = sparse(zeros(numCells,ntps));
    extractedData(chi).gaussianLocDWBgd = sparse(zeros(numCells,ntps));
    extractedData(chi).gaussianLocDWBgdQDev = sparse(zeros(numCells,ntps));
    
    extractedData(chi).ksOptBgd = sparse(zeros(numCells,ntps));
    extractedData(chi).ksOptBgdQuantile = sparse(zeros(numCells,ntps));
    extractedData(chi).alphaRadius = sparse(zeros(numCells,ntps));
    extractedData(chi).alphaRadiusArea = sparse(zeros(numCells,ntps));
    
    extractedData(chi).ksOptBgdDW = sparse(zeros(numCells,ntps));
    extractedData(chi).ksOptBgdQuantileDW = sparse(zeros(numCells,ntps));
    extractedData(chi).alphaRadiusDW = sparse(zeros(numCells,ntps));
    extractedData(chi).alphaRadiusAreaDW = sparse(zeros(numCells,ntps));
    
    extractedData(chi).xloc=sparse(zeros(numCells,ntps));
    extractedData(chi).yloc=sparse(zeros(numCells,ntps));
    
    extractedData(chi).xpos=sparse(zeros(numCells,ntps));
    extractedData(chi).ypos=sparse(zeros(numCells,ntps));
    
    extractedData(chi).trapNum = trap';
    extractedData(chi).cellNum = cells';
end

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
    rawTimepointStacks = loadImageTimer.UserData.timepointStacks;
    if t<numel(timepoints)
        % Start downloading images for next time point
        loadImageTimer.UserData.timepoint = timepoints(t+1);
        start(loadImageTimer);
    end
            
    trapInfo = cTimelapse.cTimepoint(timepoint).trapInfo;
    
    cellLocAll = false([cTimelapse.trapImSize,numel(trap)]);
    for c=1:numel(trap)
        cellInd = find(trapInfo(trap(c)).cellLabel==cells(c));
        if ~isempty(cellInd)
            loc = double(trapInfo(trap(c)).cell(cellInd).cellCenter);
            if ~isempty(loc)
                cellLocAll(:,:,c) = ...
                    full(trapInfo(trap(c)).cell(cellInd).segmented);
            end
        end
    end
    
    trapNums = unique(trap);
    cellLocAllCellsBkg = false([cTimelapse.trapImSize,numel(trapNums)]);
    for tn=1:numel(trapNums)
        trapNum = trapNums(tn);
        if trapInfo(trapNum).cellsPresent
            tSeg = cellLocAllCellsBkg(:,:,tn);
            for cellInd=1:numel(trapInfo(trapNum).cellLabel)
                tSeg = tSeg | full(trapInfo(trapNum).cell(cellInd).segmented);
            end
            cellLocAllCellsBkg(:,:,tn) = tSeg;
        end
    end
    
    % Parallelisation of the morphological operations is not helpful...
    %tic;
    for c=1:size(cellLocAll,3)
    %parfor (c=1:size(cellLocAll,3), 4)
        cellLocAll(:,:,c) = imfill(cellLocAll(:,:,c),'holes');
    end
    for tn=1:size(cellLocAllCellsBkg,3)
    %parfor (tn=1:size(cellLocAllCellsBkg,3), 4)
        cellLocAllCellsBkg(:,:,tn)=imfill(cellLocAllCellsBkg(:,:,tn),'holes');
    end
    %disp(toc);

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
        
        % Split Z-projected image into traps
        trapStack = cTimelapse.returnTrapsFromImage(tpImg,timepoint);
        
        % Also split raw sections into traps
        trapSectionStack = []; nz = size(tpStack,3);
        if nz>1
            trapSectionStack = cell(1,nz);
            for z=1:nz
                trapSectionStack{z} = cTimelapse.returnTrapsFromImage(tpStack(:,:,z),timepoint);
            end
            trapSectionStack = cat(4,trapSectionStack{:});
            trapSectionStack = permute(trapSectionStack,[1,2,4,3]);
        end
        
        for c=1:length(trap)
            % Skip this cell if not present at this time point
            cellInd = find(trapInfo(trap(c)).cellLabel==cells(c));
            if isempty(cellInd), continue; end
            
            % Skip this cell if segmentation mask is empty
            cellLoc = cellLocAll(:,:,c)>0;
            if ~sum(cellLoc(:)), continue; end
            
            trapImage = trapStack(:,:,trap(c));
            
            %vector of cell pixels
            cellFL=trapImage(cellLoc);
            nCellFL = numel(cellFL);
            
            %below is the function to extract the fluorescence information
            %from the cells. Change to mean/median FL etc
            flsorted=sort(cellFL(:),'descend');
            
            ratioOverlap=nCellFL*.025;
            numberOverlapPixels = min(ceil(ratioOverlap),nCellFL);
            
            extractedData(chi).max5px(c,timepoint)=mean(flsorted(1:min(5,nCellFL)));
            extractedData(chi).max2p5pc(c,timepoint)=mean(flsorted(1:numberOverlapPixels));
            extractedData(chi).mean(c,timepoint)=mean(cellFL(:));
            extractedData(chi).median(c,timepoint)=median(cellFL(:));
            
            convMatrix=zeros(3,3);
            convMatrix(2,:)=1;
            convMatrix(:,2)=1;
            convMatrix=imresize(convMatrix,ratioOverlap/5);
            
            flPeak=conv2(double(trapImage),convMatrix,'same');
            flPeak=flPeak(cellLoc);
            extractedData(chi).smallPeakConv(c,timepoint)=max(flPeak(:));
            
            bgdLoc = ~cellLocAllCellsBkg(:,:,find(trapNums==trap(c),1));
            bgdFL = trapImage(bgdLoc);
            bgdFL = bgdFL(~isnan(bgdFL(:)));
            if isempty(bgdFL), bgdFL = trapImage; end
            extractedData(chi).imBackground(c,timepoint)=median(bgdFL(:));
            
            % Original nucEstConv
            % Estimate concentration in nucleus
            alpha = 0.95;
            approxNucRadius = sqrt(0.085*nCellFL/pi);
            sdEst = approxNucRadius/sqrt(chi2inv(alpha,2));
            nucFiltHw = ceil(2*approxNucRadius);
            nucFilt = fspecial('gaussian',(2*nucFiltHw+1)*ones(1,2),sdEst);
            
            cellImg = trapImage-median(cellFL(:));
            cellImg(~cellLoc) = 0;
            
            nucConv = conv2(cellImg,nucFilt,'same');
            extractedData(chi).nucEstConv(c,timepoint) = ...
                max(nucConv(:))/sum(nucFilt(:).^2)*alpha ...
                /(pi*sdEst^2*chi2inv(alpha,2));
            
            trapSections = [];
            cellLocDW = cellLoc;
            if ~isempty(trapSectionStack)
                trapSections = trapSectionStack(:,:,:,trap(c));
                cellLocDW = cellLoc(:,:,ones(1,nz));
            end
            
            % nucConvEstDW: depth-wise version of original nucEstConv
            if ~isempty(trapSections)
                nucConvMaxDW = NaN(1,nz);
                for z=1:nz
                    zImg = trapSections(:,:,z);
                    zImg = zImg-median(zImg(cellLoc)); zImg(~cellLoc) = 0;
                    nucConvDW = conv2(cImgDW,nucFilt,'same');
                    nucConvMaxDW(z) = max(nucConvDW(:));
                end
                extractedData(chi).nucEstConvDW(c,timepoint) = ...
                    max(nucConvMaxDW)/sum(nucFilt(:).^2)*alpha ...
                    /(pi*sdEst^2*chi2inv(alpha,2));
            end
            
            if ~isempty(trapSections) && ~isempty(zx_ratio)
                % Estimate nuclear concentration in 3D
                nucFiltHwZ = ceil(2*approxNucRadius*zx_ratio);
                nucFilt3 = fspecial('gaussian',...
                    [(2*nucFiltHw+1)*ones(1,2),(2*nucFiltHwZ+1)],...
                    [sdEst*ones(1,2),sdEst*zx_ratio]);
                cellImg3 = tpStack-median(tpStack(:));
                cellImg3(~cellLocDW) = 0;
                convNorm = conv3(cellLocDW,nucFilt3,'same');
                nucConv3 = conv3(cellImg3,nucFilt3,'same')./convNorm;
                extractedData(chi).nucEstConv3(c,timepoint) = ...
                    max(nucConv3(:))/sum(nucFilt3(:).^2) ...
                    *alpha/(4/3*pi*approxNucRadius^3);
            end
            
            % Revised Gaussian localisation measure
            % Convolution with image is normalised by convolution with mask
            % Radius has been optimised such that, after background 
            % subtraction, the correct fraction of fluorescence is 
            % contained within the corresponding circular disk (and
            % background subtraction is optimised to match distribution of
            % disk-external pixels with a normal distribution using the K-S
            % statistic).
            nucPropRange = [0.010,2.5,0.085];
            nucRadiusRange = sqrt(nucPropRange*nFl/pi);
            testRadii = linspace(max(nucRadiusRange(1),1),min(nucRadiusRange(2),20),30);
            
            % Find location of nucleus with standard metric...
            sdLoc = nucRadiusRange(3)/sqrt(chi2inv(0.95,2));
            nucFilt = fspecial('gaussian',(2*ceil(2*nucRadiusRange(3))+1)*ones(1,2),sdLoc);
            % ...for MP
            cImgMP = proj_fl_img-median(proj_fl_img(cmask));
            cImgMP(~cmask) = 0;
            nucConv = conv2(cImgMP,nucFilt,'same')./conv2(cmask,nucFilt,'same');
            [~,maxInd] = max(nucConv(:));
            [xl,yl] = ind2sub(size(cImgMP),maxInd);
            nucDistMP = sum([mX(:)-xl,mY(:)-yl].^2,2);
            cImgMP = proj_fl_img; cImgMP(~cmask) = 0;
            
            % ...for DW
            maxInds = NaN(1,size(fl_img,3));
            nucConvMaxDW = NaN(1,size(fl_img,3));
            for z=1:size(fl_img,3)
                cImgDW = fl_img(:,:,z);
                cImgDW = cImgDW-median(cImgDW(cmask)); cImgDW(~cmask) = 0;
                nucConvDW = conv2(cImgDW,nucFilt,'same')./conv2(cmask,nucFilt,'same');
                [nucConvMaxDW(z),maxInds(z)] = max(nucConvDW(:));
            end
            [~,maxZ] = max(nucConvMaxDW);
            [xldw,yldw] = ind2sub(size(cImgDW),maxInds(maxZ));
            nucDistDW = sum([mX(:)-xldw,mY(:)-yldw].^2,2);
            cImgDW = fl_img(:,:,maxZ); cImgDW(~cmask) = 0;
            
            % Optimise K-S statistic
            nucDists = {nucDistMP,nucDistDW};
            cImgs = {cImgMP,cImgDW};
            for zt=1:numel(ztypes)
                ztype = ztypes{zt};
                nucDist = nucDists{zt};
                cImg = cImgs{zt};
                ks_Fl = NaN(1,numel(testRadii));
                ks_Bg = NaN(1,numel(testRadii));
                for r=1:numel(testRadii)
                    cFl = sort(cImg(cmask(:) & nucDist<=testRadii(r)^2))';
                    nFl = numel(cFl); F_cFl = (0:nFl-1)/(nFl-1);
                    cBg = sort(cImg(cmask(:) & nucDist>testRadii(r)^2))';
                    nBg = numel(cBg); F_cBg = (0:nBg-1)/(nBg-1);
                    if nBg == 0 || nFl == 0, continue; end
                    ks_Fl(r) = max(abs(F_cFl-normcdf(cFl,median(cFl),mad(cFl)/norminv(3/4))));
                    ks_Bg(r) = max(abs(F_cBg-normcdf(cBg,median(cBg),mad(cBg)/norminv(3/4))));
                end
                [~,rInd] = min(ks_Fl+ks_Bg);
                optR = testRadii(rInd);
                optBgd = median(cImg(cmask(:) & nucDist>optR^2));
                optFgd = median(cImg(cmask(:) & nucDist<=optR^2));
                optAlpha = sum(cImg(cmask(:) & nucDist<=optR^2)-optBgd)...
                    /sum(cImg(cmask(:))-optBgd);
                nNuc = sum(cmask(:) & nucDist<=optR^2);
                nFl = sum(cmask(:));
                F_cFl = (0:nFl-1)/(nFl-1);
                cFl = sort(cImg(cmask(:)))';
                [~,qInd] = min(abs(cFl-optBgd));
                ks_opt.(ztype){t}{trap}.ncfrac(c) = pi*optR^2/nFl;
                ks_opt.(ztype){t}{trap}.alpha(c) = optAlpha;
                ks_opt.(ztype){t}{trap}.bgdq(c) = F_cFl(qInd);
                ks_opt.(ztype){t}{trap}.bgd(c) = optBgd;
                ks_opt.(ztype){t}{trap}.fgd(c) = optFgd;
                ks_opt.(ztype){t}{trap}.Ncyt(c) = nFl;
                ks_opt.(ztype){t}{trap}.Nnuc(c) = nNuc;
                
                % Optimise background value using K-S statistic as 
                % calculated above, and then pick the radius that produces
                % the desired alpha value
                [~,rInd] = min(ks_Bg);
                bgdR = testRadii(rInd);
                bgd = median(cImg(cmask(:) & nucDist>bgdR^2));
                nuc_fl_frac = NaN(1,numel(testRadii));
                for r=1:numel(testRadii)
                    nuc_fl_frac(r) = ...
                        sum(cImg(cmask(:) & nucDist<=testRadii(r)^2)-bgd)...
                        /sum(cImg(cmask(:))-bgd);
                end
                aInd = find(nuc_fl_frac-alpha>0,1)-1;
                if isempty(aInd) || aInd<=0, continue; end
                aInds = aInd+(0:1);
                alphaRadius = testRadii(aInd)+diff(testRadii(aInds))...
                    *(alpha-nuc_fl_frac(aInd))/diff(nuc_fl_frac(aInds));
                [~,qInd] = min(abs(cFl-bgd));
                ks_opt.(ztype){t}{trap}.optBgd(c) = bgd;
                ks_opt.(ztype){t}{trap}.optBgdQ(c) = F_cFl(qInd);
                ks_opt.(ztype){t}{trap}.alphaRadius(c) = alphaRadius;
            end
            
            extractedData(chi).gaussianLoc(c,timepoint) = NaN;
            extractedData(chi).gaussianLocBgd(c,timepoint) = NaN;
            extractedData(chi).gaussianLocBgdQDev(c,timepoint) = NaN;
            extractedData(chi).gaussianLocDW(c,timepoint) = NaN;
            extractedData(chi).gaussianLocDWBgd(c,timepoint) = NaN;
            extractedData(chi).gaussianLocDWBgdQDev(c,timepoint) = NaN;
            
            extractedData(chi).gaussianLocUnscaled(c,timepoint) = NaN;
            extractedData(chi).gaussianLocDWUnscaled(c,timepoint) = NaN;
            
            extractedData(chi).ksOptBgd(c,timepoint) = NaN;
            extractedData(chi).ksOptBgdQuantile(c,timepoint) = NaN;
            extractedData(chi).alphaRadius(c,timepoint) = NaN;
            extractedData(chi).alphaRadiusArea(c,timepoint) = NaN;
            
            extractedData(chi).ksOptBgdDW(c,timepoint) = NaN;
            extractedData(chi).ksOptBgdQuantileDW(c,timepoint) = NaN;
            extractedData(chi).alphaRadiusDW(c,timepoint) = NaN;
            extractedData(chi).alphaRadiusAreaDW(c,timepoint) = NaN;
    
            % information common to all channels (basically
            % shape information) is stored only in the
            % channel 1 structure.
            if chi==1
                extractedData(chi).area(c,timepoint)=nCellFL;
                extractedData(chi).radius(c,timepoint)=trapInfo(trap(c)).cell(cellInd).cellRadius;
                xyloc = trapInfo(trap(c)).cell(cellInd).cellCenter(1);
                extractedData(chi).xloc(c,timepoint)=xyloc(1);
                extractedData(chi).yloc(c,timepoint)=xyloc(2);
                trapLoc = cTimelapse.cTimepoint(timepoint).trapLocations(trap(c));
                extractedData(chi).xpos(c,timepoint)=trapLoc.xcenter-trapHw-1+xyloc(1);
                extractedData(chi).ypos(c,timepoint)=trapLoc.ycenter-trapHh-1+xyloc(2);
            end
        end
    end
end

% Clean up timer since it persists beyond function termination
delete(loadImageTimer);

cTimelapse.extractedData = extractedData;

fprintf('\n');

    function getImageData(timerobj,~)
        chNums = timerobj.UserData.channelNums;
        tp = timerobj.UserData.timepoint;
        for cn=1:numel(chNums)
            chnum = chNums(cn);
            timerobj.UserData.timepointStacks{cn} = ...
                cTimelapse.returnSingleTimepointRaw(tp,chnum);
        end
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

function iscomplete = recurse_fields(x,f)
iscomplete = isempty(f) || ...
    (isfield(x,f{1}) && recurse_fields(x.(f{1}),f(2:end)));
end
