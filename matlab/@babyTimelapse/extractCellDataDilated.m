function extractCellDataDilated(cTimelapse)
%extractCellDataDilated Custom extraction method to extract data for 
%   dilated segmentation images
%   Based on "extractCellDataStandardParfor"
%
% Uses parfor loop to speed up the extraction.
%
% THE MAX5 IS THE MEAN OF THE 5 BRIGHTEST PIXELS.
% The MAX25P IS THE BRIGHTEST 2.5% OF THE PIXELS OF THE CELL. THIS ALLOWS
% FOR SCALING WITH CELL SIZE SO YOU DON'T GET DIFFERENT RELATIVE 
% LOCALIZATION SCORES DEPENDING ON CELL SIZE.
%
% Uses the parameters structure
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
% The standard data extraction method. Calculates a raft of statistics,
% hopefully mostly self explanatory, about the cell pixels. 
%
% Size and position information common to all cells is stored in the first
% entry of extractedData only.

parameters = cTimelapse.extractionParameters.functionParameters;

% Standard parameters
type = parameters.type;
channels = parameters.channels;

% Custom parameters for NucLocMeasures
if isfield(parameters,'timepoints')
    timepoints = parameters.timepoints;
    timepoints = timepoints(cTimelapse.timepointsProcessed(timepoints));
else
    timepoints = find(cTimelapse.timepointsProcessed);
end

if strcmp(channels,'all')
    channels = 1:length(cTimelapse.channelNames);
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

se1=strel('disk',1);
se2=strel('disk',2);

extractedData = struct();

%preallocate cellInf
for channel=1:length(channels)
    extractedData(channel).radius=sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    extractedData(channel).area=sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    extractedData(channel).eccentricity=sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    extractedData(channel).segmentedRadius=sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    extractedData(channel).xloc=sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    extractedData(channel).yloc=sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    
    extractedData(channel).mean=sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    extractedData(channel).median=sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    extractedData(channel).max5=sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    extractedData(channel).max25p=sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    extractedData(channel).std=sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    extractedData(channel).min=sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    extractedData(channel).imBackground=sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    extractedData(channel).imBackgroundOuter=sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    
    extractedData(channel).outerArea=sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    extractedData(channel).outerMax5=sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    extractedData(channel).outerMedian=sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    extractedData(channel).outerMean=sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    
    extractedData(channel).membraneArea=sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    extractedData(channel).membraneMax5=sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    extractedData(channel).membraneMedian=sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    extractedData(channel).membraneMean=sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    
    extractedData(channel).trapNum = trap';
    extractedData(channel).cellNum = cells';
    
end

% cell array of images for each channel extracted at each timpoint
tpStacks = cell(size(channels));

% rearranged so that each cell is filled in at the same time at each
% timepoint, this allows parfor use and reduces the extraction time
% significantly
for timepoint=timepoints
    disp(['Timepoint Number ',int2str(timepoint)]);
    
    for channel=1:length(channels)
        channel_number = channels(channel);
        
        %switch to doube so that mathematical operations are as expected.
        tpStack=cTimelapse.returnSingleTimepoint(timepoint,channel_number,'stack');
        
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
    
    trapInfo=cTimelapse.cTimepoint(timepoint).trapInfo;
    
    seg_areas=tpStacks{channel}(:,:,1);
    cellLocAll=false([size(seg_areas) length(trap)]);
    outerLocAll=false([size(seg_areas) length(trap)]);
    for allIndex=1:length(trap)
        temp_loc=find(trapInfo(trap(allIndex)).cellLabel==cells(allIndex));
        if ~isempty(temp_loc)
            seg_areas=full(trapInfo(trap(allIndex)).cell(temp_loc).segmented);
            loc=double(trapInfo(trap(allIndex)).cell(temp_loc).cellCenter);
            if ~isempty(loc)
                cellLocAll(:,:,allIndex)=seg_areas(:,:,1);
            end
        end
    end
    
    cellLocAllCellsBkg=false([size(seg_areas) length(trap)]);
    cellLocAllCellsOuterBkg=false([size(seg_areas) length(trap)]);
    for allTrapIndex=1:length(trap)
        if trapInfo(trap(allIndex)).cellsPresent
            for allCells=1:length(trapInfo(trap(allIndex)).cellLabel)
                tSeg=full(trapInfo(trap(allIndex)).cell(allCells).segmented);
                cellLocAllCellsBkg(:,:,allTrapIndex)=cellLocAllCellsBkg(:,:,allTrapIndex)|tSeg;
                cellLocAllCellsOuterBkg(:,:,allTrapIndex)=cellLocAllCellsBkg(:,:,allTrapIndex)|tSeg;
            end
        end
    end
    
    t=size(cellLocAllCellsBkg,3);
    parfor allIndex=1:size(cellLocAll,3)
        cellLocAll(:,:,allIndex)=imfill(cellLocAll(:,:,allIndex),'holes');
        outerLocAll(:,:,allIndex)=imdilate(cellLocAll(:,:,allIndex),se2);
        if allIndex<=t
            cellLocAllCellsBkg(:,:,allIndex)=imfill(cellLocAllCellsBkg(:,:,allIndex),'holes');
            cellLocAllCellsOuterBkg(:,:,allIndex)=cellLocAllCellsBkg(:,:,allIndex);
            cellLocAllCellsOuterBkg(:,:,allIndex)=imdilate(cellLocAllCellsOuterBkg(:,:,allIndex),se2);
        end
    end
    
    trapInfo=cTimelapse.cTimepoint(timepoint).trapInfo;
    for channel=1:length(channels)
        
        if ~all(tpStacks{channel}(:)==0)
            %if empty do nothing
            
            tpImCh=tpStacks{channel};
            
            for allIndex =1:length(trap)
                trapImage=tpImCh(:,:,trap(allIndex));
                temp_loc=find(trapInfo(trap(allIndex)).cellLabel==cells(allIndex));
                cellLoc=cellLocAll(:,:,allIndex)>0;
                % give the row in cellInf in which data for this cell
                % shoud be inserted.
                if ~isempty(temp_loc) && sum(cellLoc(:))
                    
                    %vector of cell pixels
                    cellFL=trapImage(cellLoc);
                    
                    %below is the function to extract the fluorescence information
                    %from the cells. Change to mean/median FL etc
                    flsorted=sort(cellFL(:),'descend');
                    
                    %                     numberOverlapPixels = min(maxPixOverlap,length(cellFL));
                    ratioOverlap=ceil(length(cellFL(:))*.025);
                    numberOverlapPixels = min(ratioOverlap,length(cellFL));
                    
                    nmax5 = min(5,length(flsorted));
                    extractedData(channel).max5(allIndex,timepoint)=mean(flsorted(1:nmax5));
                    extractedData(channel).max25p(allIndex,timepoint)=mean(flsorted(1:numberOverlapPixels));
                    extractedData(channel).mean(allIndex,timepoint)=mean(cellFL(:));
                    extractedData(channel).median(allIndex,timepoint)=median(cellFL(:));
                    extractedData(channel).std(allIndex,timepoint)=std(cellFL(:));
                    extractedData(channel).min(allIndex,timepoint)=min(cellFL(:));
                    
                    % Calculate stats for expanded cell
                    outer_area=outerLocAll(:,:,allIndex)>0;
                    outerFL=trapImage(outer_area);
                    flsorted=sort(outerFL(:),'descend');
                    nmax5 = min(5,length(flsorted));
                    extractedData(channel).outerMax5(allIndex,timepoint)=mean(flsorted(1:nmax5));
                    extractedData(channel).outerMean(allIndex,timepoint)=mean(outerFL(:));
                    extractedData(channel).outerMedian(allIndex,timepoint)=median(outerFL(:));
                    extractedData(channel).outerArea(allIndex,timepoint)=length(outerFL);
                    
                    % Calculate stats for expanded cell minus original 
                    membrane_area=outerLocAll(:,:,allIndex)>0 & ~cellLoc;
                    membraneFL=trapImage(membrane_area);
                    flsorted=sort(membraneFL(:),'descend');
                    nmax5 = min(5,length(flsorted));
                    extractedData(channel).membraneMax5(allIndex,timepoint)=mean(flsorted(1:nmax5));
                    extractedData(channel).membraneMean(allIndex,timepoint)=mean(membraneFL(:));
                    extractedData(channel).membraneMedian(allIndex,timepoint)=median(membraneFL(:));
                    extractedData(channel).membraneArea(allIndex,timepoint)=length(membraneFL);
                    
                    seg_areas=~cellLocAllCellsBkg(:,:,allIndex);
                    
                    bkg=trapImage(seg_areas);
                    bkg=bkg(~isnan(bkg(:)));
                    if isempty(bkg)
                        bkg=trapImage;
                    end
                    extractedData(channel).imBackground(allIndex,timepoint)=median(bkg(:));
                    
                    seg_areas=~cellLocAllCellsOuterBkg(:,:,allIndex);
                    
                    bkg=trapImage(seg_areas);
                    bkg=bkg(~isnan(bkg(:)));
                    if isempty(bkg)
                        bkg=trapImage;
                    end
                    extractedData(channel).imBackgroundOuter(allIndex,timepoint)=median(bkg(:));
                    
                    % information common to all channels (basically
                    % shape information) is stored only in the
                    % channel 1 structure.
                    if channel==1
                        extractedData(channel).area(allIndex,timepoint)=length(cellFL);
                        extractedData(channel).radius(allIndex,timepoint)=trapInfo(trap(allIndex)).cell(temp_loc).cellRadius;
                        tP=regionprops(cellLoc,'Eccentricity');
                        extractedData(channel).eccentricity(allIndex,timepoint)=tP.Eccentricity;
                        
                        extractedData(channel).segmentedRadius(allIndex,timepoint)=sqrt(sum(cellLoc(:))/pi);
                        
                        extractedData(channel).xloc(allIndex,timepoint)=trapInfo(trap(allIndex)).cell(temp_loc).cellCenter(1);
                        extractedData(channel).yloc(allIndex,timepoint)=trapInfo(trap(allIndex)).cell(temp_loc).cellCenter(2);
                    end
                    
                end
                
            end
        end
    end
end

cTimelapse.extractedData = extractedData;
end
