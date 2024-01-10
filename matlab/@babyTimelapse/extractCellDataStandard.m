function extractCellDataStandard(cTimelapse)
% extractCellDataStandard(cTimelapse)
%
% standard extraction function to extract all the commonly used
% measurement for the cells. Uses the parameters structure
%
%       cTimelapse.extractionParameters.functionParameters
%
% with fields:
% channels                  : array of channels to extract or the string 'all', in
%                             which case all channels are extracted
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
% nuclearMarkerChannel    : If a nuclear tag is used, the number of that channel.
%                           if NaN then this is ignored and nuclear
%                           localisation set to zero.
% maxPixOverlap           : number of nuclear pixels used to calculate
%                           nuclearTagLoc. Also used for max5 and membrane max5. 
% 
% maxAllowedOverlap       : number of candidate nuclear pixels
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
% cell Halo pixels - the pixels obtained when the cell is dilated by 3 and
%                    the cell area then removed. Basically, a ring of width
%                    3 pixels around the cell.
%
% 'radius' - the approximate radius of the cell calculated as the mean of
%            the radii used to make the active contour outline.
% 'area' - number of pixels which make up the cell.
% 'eccentricity' - eccentricity of the cell taken from regionprops
% 'mean' - the mean of the cell pixels.
% 'median' - median of the cell pixels.
% 'max5' - the mean of the brightest maxPixOverlap of cell pixels.
% 'std' - std of the cell pixels.
% 'proteinLocalization' - the 'max5' field divided by the 'median' field.
%                         Commonly used for nuclear localisation. 
% 'min' - the value of the dimmest cell pixel.
% 'imBackground' - median of all pixels which are not part of any cell.
% 'membraneMax5' - mean of the brightest 2.5% of membrane pixels.
% 'membraneMedian' - median of membrane pixels.
% 'nuclearTagLoc' - If a nuclear channel is specified, this is obtained by
%                   taking the brightest MaxAllowedOverlap pixels in the
%                   of the cell pixels in the nuclear channel as candidate
%                   pixels. It then takes the mean of the brightest
%                   maxPixOverlap pixels in the image channel amongst the
%                   candidate pixels. This, divided by the median of all
%                   the cell pixels in the image channel, is the
%                   'nuclearTagLoc'
% 'smallmean' - the mean of the reduced cell pixels.
% 'smallmedian' - the median of the reduced cell pixels.
% 'smallmax5' - To obtain this data field the image is convolved with a
%               a small disc of approximately the size of 2.5% of the cell
%               pixels. 'smallmax5' is the maximum pixel within the cell
%               outline in this resulting image. If bright pixels are
%               clumped together, this should be approximately the same as
%               'max5'.
% 'distToNuc' - (legacy) field populated by extractNucAreaFL method.
% 'nucArea' - (legacy) field populated by extractNucAreaFL method.
% 'segmentedRadius' - radius calculated from 'area' assuming a circular
%                     cell
% 'xloc' - x location in the trap. (i.e.relative to trap centre).
% 'yloc' - y location in the trap. (i.e.relative to trap centre).
% 'cellHaloMedian' - median of the cell Halo pixels.
% 'cellHaloMean' - mean of the cell Halo pixels.
% 'cellHaloQ95' - 9th percentile of the cell Halo pixels.
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

type = parameters.type;
channels = parameters.channels;
nuclearMarkerChannel = parameters.nuclearMarkerChannel;
maxPixOverlap = parameters.maxPixOverlap;
maxAllowedOverlap = parameters.maxAllowedOverlap;

%number of candidate pixels should not be larger than number of finally
%allowed centre pixels.
if maxAllowedOverlap<maxPixOverlap
    maxAllowedOverlap = maxPixOverlap;
end


if strcmp(channels,'all')
    channels = 1:length(cTimelapse.channelNames);
end

if ~ismember(nuclearMarkerChannel,channels)
    nuclearMarkerChannel = NaN;
end

numCells=sum(cTimelapse.cellsToPlot(:));
[trap, cells]=find(cTimelapse.cellsToPlot);

% annoying necessary line to deal with find behaviour. If find is applied
% to a row vector, it returns a row vector. If it is applied to a matrix it
% returns a column vector.
trap = trap(:);
cells = cells(:);

%reorder so cells in the same trap contiguous, nicer for viewing later.
[trap,I] = sort(trap);
cells =cells(I);

%preallocate cellInf
for channel=1:length(channels)
    
    extractedData(channel).mean=sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    extractedData(channel).median=sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    extractedData(channel).max5=sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    extractedData(channel).std=sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    extractedData(channel).smallmean=sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    extractedData(channel).smallmedian=sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    extractedData(channel).smallmax5=sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    extractedData(channel).min=sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    extractedData(channel).imBackground=sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    extractedData(channel).area=sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    extractedData(channel).radius=sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    extractedData(channel).distToNuc=sparse(numCells,length(cTimelapse.timepointsProcessed));
    extractedData(channel).radiusFL=sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    extractedData(channel).segmentedRadius=sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    extractedData(channel).xloc=sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    extractedData(channel).yloc=sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    
    extractedData(channel).membraneMax5=sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    extractedData(channel).membraneMedian=sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    extractedData(channel).nuclearTagLoc=sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    
    extractedData(channel).cellHaloMedian = sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    extractedData(channel).cellHaloMean = sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    extractedData(channel).cellHaloQ95 = sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    extractedData(channel).imBackgroundDistant = sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));

    extractedData(channel).trapNum = trap';
    extractedData(channel).cellNum = cells';
    
end

% cell array of images for each channel extracted at each timpoint
tpStacks = cell(size(channels));

for timepoint=find(cTimelapse.timepointsProcessed)
    %disp(['Timepoint Number ',int2str(timepoint)]);
    % Trigger the TimepointChanged event for babyLogging
    babyLogging.changeTimepoint(cTimelapse,timepoint);
    
    for channel=1:length(channels)
        channel_number = channels(channel);
        
        %switch to doube so that mathematical operations are as expected.
        tpStack=double(cTimelapse.returnSingleTimepoint(timepoint,channel_number,'stack'));
        
        % if the channels is the nuclear marker channel, populate the
        % nuclearStack array with max projection, which is taken as a
        % marker of nuclearity.
        if channel_number == nuclearMarkerChannel
            nuclearStack = cTimelapse.returnTrapsFromImage(max(tpStack,[],3),timepoint);
        end
        
        switch type
            case 'max'
                tpStack = max(tpStack,[],3);
            case 'min'
                tpStack = min(tpStack,[],3);
            case 'mean'
                tpStack = mean(tpStack,3);
            case 'std'
                tpStack = std(tpStack,[],3);
            case 'sum'
                tpStack = sum(tpStack,3);
        end
        
        tpStacks{channel} = cTimelapse.returnTrapsFromImage(tpStack,timepoint);
        
    end
    
    for channel=1:length(channels)
        channel_number = channels(channel);
        
        if ~all(tpStacks{channel}(:)==0)
            %if empty do nothing
            
            trapInfo=cTimelapse.cTimepoint(timepoint).trapInfo;
            uniqueTraps=unique(trap);
            
            for j=1:length(uniqueTraps)
                currTrap=uniqueTraps(j);
                cellsCurrTrap=cells(trap==currTrap);
                
                %protect code from cases where no cells have been found and
                %the cell structure is weird and weak
                if ~trapInfo(currTrap).cellsPresent
                    cellsCurrTrap = [];
                end
                
                trapImage=tpStacks{channel}(:,:,currTrap);
                
                if ~isnan(nuclearMarkerChannel)
                    nuclearTrapImage = nuclearStack(:,:,currTrap);
                end
                
                
                for cellIndex=1:length(cellsCurrTrap)
                    currCell=cellsCurrTrap(cellIndex);
                    temp_loc=find(trapInfo(currTrap).cellLabel==currCell);
                    
                    % give the row in cellInf in which data for this cell
                    % shoud be inserted.
                    dataInd = find(trap==currTrap & cells == currCell);
                    if ~isempty(temp_loc)
                        if length({trapInfo(currTrap).cell(temp_loc).segmented})~=1
                            error('Oh no!!');
                        end
                        seg_areas=full(trapInfo(currTrap).cell(temp_loc).segmented);
                        
                        cellLoc=zeros(size(seg_areas));
                        loc=double(cTimelapse.cTimepoint(timepoint).trapInfo(currTrap).cell(temp_loc).cellCenter);
                        if ~isempty(loc)
                            cellLoc=imfill(seg_areas(:,:,1),'holes');
                        end
                        %logical of cell pixels
                        cellLoc=cellLoc>0;
                        
                        %logical of membrane pixels
                        membraneLoc = seg_areas >0;
                        
                        %vector of cell pixels
                        cellFL=trapImage(cellLoc);
                        
                        %vector of cell membrane pixels
                        membraneFL = trapImage(membraneLoc);
                        
                        %below is the function to extract the fluorescence information
                        %from the cells. Change to mean/median FL etc
                        flsorted=sort(cellFL(:),'descend');
                        mflsorted=sort(membraneFL(:),'descend');
                        
                        numberOverlapPixels = min(maxPixOverlap,length(cellFL));
                        
                        %nuclear extraction
                        if ~isnan(nuclearMarkerChannel)
                            cellFLnuclear = nuclearTrapImage(cellLoc);
                            
                            [~, indMax5nuclear]=sort(cellFLnuclear(:),'descend'); %nuclear tag pixels sorted
                            
                            allowedNuclearPixels=false(size(cellFLnuclear));
                            
                            numberAllowedPixels = min(maxAllowedOverlap,length(indMax5nuclear));
                            allowedNuclearPixels(indMax5nuclear(1:numberAllowedPixels))=1;
                            nucleusValHOG=sort(cellFL(allowedNuclearPixels),'descend');
                            
                            nuclocalization =  mean(nucleusValHOG(1:numberOverlapPixels))/median(cellFL(:));
                            
                            extractedData(channel).nuclearTagLoc(dataInd,timepoint)=nuclocalization;
                            %end nuclear extraction
                            
                        end
                        
                        
                        convMatrix=zeros(3,3);
                        convMatrix(2,:)=1;
                        convMatrix(:,2)=1;
                        
                        extractedData(channel).max5(dataInd,timepoint)=mean(flsorted(1:numberOverlapPixels));
                        extractedData(channel).mean(dataInd,timepoint)=mean(cellFL(:));
                        extractedData(channel).median(dataInd,timepoint)=median(cellFL(:));
                        extractedData(channel).std(dataInd,timepoint)=std(cellFL(:));
                        extractedData(channel).min(dataInd,timepoint)=min(cellFL(:));
                        
                        extractedData(channel).membraneMedian(dataInd, timepoint)=median(membraneFL(:));
                        extractedData(channel).membraneMax5(dataInd, timepoint)=mean(mflsorted(1:numberOverlapPixels));
                        
                        cellLocSmall=imerode(cellLoc,strel('disk',2));
                        
                        cellFLsmall=trapImage(cellLocSmall);
                        
                        flPeak=conv2(double(trapImage),convMatrix);
                        flPeak=flPeak(cellLoc);
                        
                        extractedData(channel).smallmax5(dataInd,timepoint)=max(flPeak(:));
                        extractedData(channel).smallmean(dataInd,timepoint)=mean(cellFLsmall(:));
                        extractedData(channel).smallmedian(dataInd,timepoint)=median(cellFLsmall(:));
                        
                        % Quantify the effect of neighbouring cells, by
                        % calculating statistics for the cell halo:
                        cellHaloLoc = imdilate(cellLoc,strel('disk',3));
                        cellHaloLoc(cellLoc) = false;
                        cellHaloFL = trapImage(cellHaloLoc);
                        
                        extractedData(channel).cellHaloMedian(dataInd,timepoint) = median(cellHaloFL(:));
                        extractedData(channel).cellHaloMean(dataInd,timepoint) = mean(cellHaloFL(:));
                        extractedData(channel).cellHaloQ95(dataInd,timepoint) = quantile(cellHaloFL(:),0.95);
                        
                        seg_areas=zeros(size(trapInfo(currTrap).cell(1).segmented));
                        for allCells=1:length(trapInfo(currTrap).cellLabel)
                            
                            tSeg=full(trapInfo(currTrap).cell(allCells).segmented);
                            seg_areas=seg_areas|tSeg;
                            
                            loc=double(cTimelapse.cTimepoint(timepoint).trapInfo(currTrap).cell(allCells).cellCenter);
                            if ~isempty(loc)
                                seg_areas=imfill(seg_areas(:,:,1),'holes');
                            end
                        end
                        
                        % Calculate image background for regions distant to
                        % cell locations by dilating seg_areas to avoid 
                        % fluorescence 'halo' around cells:
                        dilated_seg_areas = imdilate(seg_areas,strel('disk',3));
                        bgd = trapImage(~dilated_seg_areas);
                        bgd = bgd(~isnan(bgd(:)));
                        if isempty(bgd)
                            % Trap is overrun with cells, there is no
                            % sensible image background measure
                            extractedData(channel).imBackgroundDistant(dataInd,timepoint)=NaN;
                        else
                            extractedData(channel).imBackgroundDistant(dataInd,timepoint)=...
                                median(bgd(:));
                        end
                        
                        % Standard image background calculation
                        seg_areas=~seg_areas;
                        
                        bkg=trapImage(seg_areas);
                        bkg=bkg(~isnan(bkg(:)));
                        if isempty(bkg)
                            bkg=trapImage;
                        end
                        extractedData(channel).imBackground(dataInd,timepoint)=median(bkg(:));
                        
                        % information common to all channels (basically
                        % shape information) is stored only in the
                        % channel 1 structure.
                        if channel==1
                            extractedData(channel).area(dataInd,timepoint)=length(cellFL);
                            extractedData(channel).radius(dataInd,timepoint)= trapInfo(currTrap).cell(temp_loc).cellRadius;
                            
                            %radiusFL populated by extractSegAreaFl
                            %method.
                            if isfield(trapInfo(currTrap).cell(temp_loc),'cellRadiusFL');
                                extractedData(channel).radiusFL(dataInd,timepoint)= trapInfo(currTrap).cell(temp_loc).cellRadiusFL;%trapInfo(currTrap).cell(temp_loc).cellRadius;
                            end
                            extractedData(channel).segmentedRadius(dataInd,timepoint)= sqrt(sum(cellLoc(:))/pi);%trapInfo(currTrap).cell(temp_loc).cellRadius;
                            
                            % nucArea populated by extractNucAreaFL method.
                            if isfield(trapInfo(currTrap).cell(temp_loc),'nucArea');
                                if isempty(trapInfo(currTrap).cell(temp_loc).nucArea)
                                    extractedData(channel).nucArea(dataInd,timepoint)=NaN;
                                    extractedData(channel).distToNuc(dataInd,timepoint)=NaN;
                                else
                                    extractedData(channel).nucArea(dataInd,timepoint)=trapInfo(currTrap).cell(temp_loc).nucArea;
                                    extractedData(channel).distToNuc(dataInd,timepoint)=trapInfo(currTrap).cell(temp_loc).distToNuc;
                                end
                            end
                            
                            extractedData(channel).xloc(dataInd,timepoint)= trapInfo(currTrap).cell(temp_loc).cellCenter(1);
                            extractedData(channel).yloc(dataInd,timepoint)=trapInfo(currTrap).cell(temp_loc).cellCenter(2);
                            
                        end
                    end
                end
            end
        end
    end
end
cTimelapse.extractedData=extractedData;
end
