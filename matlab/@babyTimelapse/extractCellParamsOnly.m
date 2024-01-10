function extractCellParamsOnly(cTimelapse)
% extractCellParamsOnly(cTimelapse)
%
% very similar to extractData, just cut back. Extracts only those results
% accessible with no image loading and stores them cTimelapse.extractedData
% (which is then a 1x1 structure array).
% 
% lack of image loading nd processing makes it much faster.
%
% written not work if there is only one timepoint

numCells=sum(cTimelapse.cellsToPlot(:));
[trap,cells]=find(cTimelapse.cellsToPlot);
trap = trap(:); cells = cells(:);
%reorder so cells in the same trap contiguous, nicer for viewing later.
[trap,I] = sort(trap); % NB: sort is stable, so cells should be ordered
cells =cells(I);

if isempty(cTimelapse.timepointsProcessed) || length(cTimelapse.timepointsProcessed)==1
    tempSize=[cTimelapse.cTimepoint.trapInfo];
    cTimelapse.timepointsProcessed=ones(1,length(tempSize)/length(cTimelapse.cTimepoint(1).trapInfo));
    if length(cTimelapse.timepointsProcessed)==1
        cTimelapse.timepointsProcessed=0;
    end
end

extractedData = struct();
for channel=1
    extractedData(channel).mean=sparse(numCells,length(cTimelapse.timepointsToProcess));
    extractedData(channel).median=sparse(numCells,length(cTimelapse.timepointsToProcess));
    extractedData(channel).max5=sparse(numCells,length(cTimelapse.timepointsToProcess));
    extractedData(channel).std=sparse(numCells,length(cTimelapse.timepointsToProcess));
    
    extractedData(channel).smallmean=sparse(numCells,length(cTimelapse.timepointsToProcess));
    extractedData(channel).smallmedian=sparse(numCells,length(cTimelapse.timepointsToProcess));
    extractedData(channel).smallmax5=sparse(numCells,length(cTimelapse.timepointsToProcess));
    extractedData(channel).min=sparse(numCells,length(cTimelapse.timepointsToProcess));
    extractedData(channel).imBackground=sparse(numCells,length(cTimelapse.timepointsToProcess));

    extractedData(channel).distToNuc=sparse(numCells,length(cTimelapse.timepointsToProcess));
    extractedData(channel).nucArea=sparse(numCells,length(cTimelapse.timepointsToProcess));
    extractedData(channel).radius=sparse(numCells,length(cTimelapse.timepointsToProcess));
    extractedData(channel).radiusAC=sparse(zeros(numCells,length(cTimelapse.timepointsProcessed)));
    extractedData(channel).xloc=sparse(numCells,length(cTimelapse.timepointsToProcess));
    extractedData(channel).yloc=sparse(numCells,length(cTimelapse.timepointsToProcess));
    
    extractedData(channel).trapNum=trap;
    extractedData(channel).cellNum=cells;
    
    for timepoint=1:length(cTimelapse.timepointsToProcess)
        if cTimelapse.timepointsProcessed(timepoint)
            trapInfo=cTimelapse.cTimepoint(timepoint).trapInfo;
            for j=1:length(extractedData(channel).cellNum)
                currCell=extractedData(channel).cellNum(j);
                currTrap=extractedData(channel).trapNum(j);
                
                temp_loc=find(trapInfo(currTrap).cellLabel==currCell);
                if temp_loc 
                    extractedData(channel).radius(j,timepoint)=trapInfo(currTrap).cell(temp_loc).cellRadius;
                    extractedData(channel).xloc(j,timepoint)=trapInfo(currTrap).cell(temp_loc).cellCenter(1);
                    extractedData(channel).yloc(j,timepoint)=trapInfo(currTrap).cell(temp_loc).cellCenter(2);

                    if isfield(trapInfo(currTrap).cell(temp_loc),'nucArea')
                        if isempty(trapInfo(currTrap).cell(temp_loc).nucArea)
                            extractedData(channel).nucArea(j,timepoint)=NaN;
                            extractedData(channel).distToNuc(j,timepoint)=NaN;
                        else
                            extractedData(channel).nucArea(j,timepoint)=trapInfo(currTrap).cell(temp_loc).nucArea;
                            extractedData(channel).distToNuc(j,timepoint)=trapInfo(currTrap).cell(temp_loc).distToNuc;
                        end
                    end
                    if isfield(trapInfo(currTrap).cell(temp_loc),'radiusAC')
                        extractedData(channel).radiusAC(j,timepoint)=trapInfo(currTrap).cell(temp_loc).radiusAC;
                    end
                end
            end
        end
    end
end
cTimelapse.extractedData=extractedData;