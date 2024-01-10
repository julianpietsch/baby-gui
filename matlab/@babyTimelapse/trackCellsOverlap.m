function trackCellsOverlap(cTimelapse,cellMovementThresh)
%trackCells: Tracks individual cells over multiple images, and finds trapinfo
%---------------------------------------------------------
%       Function of babyTimelapse class to assign labels to the cells in
%       each timepoint, and so track the individual cells between the 
%       timepoints. Also calculates if cells are present in each trap if
%       traps are present.
%
%       Calculates the distance between all cell centers between one timeframe ant
%       the next and uses this, weighted for factors like cell size etc.,
%       to find the most likely candidate for the same cell. Each cell is
%       then assigned a label.
%
%IMPORTANT VARS:    cTimelapse.cTimepoint(i).cellLabel
%                   1xN double vector
%                     cellLabel is a vector which contains labels for all
%                     of the cells being tracked. 
%                     The label itself will be a number, corresponding to
%                     the cell's position in the first timepoint. The
%                     cell's number in the current timepoint i is given by
%                     the position in the vector. This gives a way to link
%                     a given cell in all timepoints, by matching the
%                     labels.
%
%INPUT:             cTimelapse
%                   babyTimelapse
%                     babyTimelapse object with all cells already segmented
%
%                   cellMovementThresh
%                   Double
%                     Number to indicate how far (in pixels) a cell can
%                     move before it should be considered a new cell. 
%                     Optional argument, if left blank it is assigned by
%                     UI.
if nargin<2
    prompt = {'Max change in position and radius before a cell is classified as a new cell'};
    dlg_title = 'Tracking Threshold';
    num_lines = 1;
    def = {'5'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    cellMovementThresh=str2double(answer{1});
end

%%
%Identify the mother index ... cells that are closest to the center of the
%trap and most likely to be the cells of interest. The tracking will be
%more lenient on these cells.
% motherIndex=cTimelapse.findMotherIndex;

%%
if isempty(cTimelapse.timepointsProcessed)
    tempSize=[cTimelapse.cTimepoint.trapInfo];
    cTimelapse.timepointsProcessed=ones(1,length(tempSize)/length(cTimelapse.cTimepoint(1).trapInfo));
end

trapMaxCell=[];
for timepoint=cTimelapse.timepointsToProcess
    if cTimelapse.timepointsProcessed(timepoint)
        % Trigger the TimepointChanged event for babyLogging
        babyLogging.changeTimepoint(cTimelapse,timepoint);
        
        %Shuffles data around, trapinfo is data of current timepoint, 
        %trapinfom1 is data of previous timepoint,
        %trapinfom2 is data of timepoint before that.
        if length(cTimelapse.timepointsToProcess)>2 && timepoint>cTimelapse.timepointsToProcess(2) 
            trapInfom2=trapInfom1;
        end
        if timepoint>cTimelapse.timepointsToProcess(1)
            trapInfom1=trapInfo;
        end
        trapInfo=cTimelapse.cTimepoint(timepoint).trapInfo;
        
        %this is to correct for a bug in some old timelapses that I was
        %processing, shouldn't be needed generally. For when a timepoint
        %wasn't processed but the tp before and after was
        if length(cTimelapse.timepointsToProcess)>1 && timepoint>cTimelapse.timepointsToProcess(2) && ~cTimelapse.timepointsProcessed(timepoint-1)
            trapInfom2=trapInfom1;
        end
        
        currTPCells=[];
        for trapInd=1:length(trapInfo)
            %fill in the current tp cells so that they can be used in the
            %future for the overlap labelling
            for cellInd=1:length(trapInfo(trapInd).cell)
                tempCell=full(trapInfo(trapInd).cell(cellInd).segmented);
                tempCell=imfill(tempCell,'holes');
                currTPCells{trapInd}(:,:,cellInd)=tempCell;
            end
            
            if timepoint==cTimelapse.timepointsToProcess(1)
                if trapInfo(trapInd).cellsPresent
                    len=length(trapInfo(trapInd).cell);
                    trapInfo(trapInd).cellLabel=1:length(trapInfo(trapInd).cell);
                else
                    len=0;
                    trapInfo(trapInd).cellLabel=0;
                    trapInfo(trapInd).cellLabel=0;
                    trapInfo(trapInd).cell(1).cellCenter=[];
                    trapInfo(trapInd).cell(1).cellRadius=[];
                end
                trapMaxCell(trapInd)=len;
            elseif ~trapInfo(trapInd).cellsPresent
                trapInfo(trapInd).cellLabel=0;
                trapInfo(trapInd).cell(1).cellCenter=[];
                trapInfo(trapInd).cell(1).cellRadius=[];
            else
                trapInfo(trapInd).cellLabel=zeros(1,length(trapInfo(trapInd).cell));
                
                cellOverlap=[]; %curr Cell row, test prev TP 
                for cellInd=1:length(trapInfo(trapInd).cell)
                    tCell=currTPCells{trapInd}(:,:,cellInd);
                    try
                        for prevCellInd=1:size(prevTPCells{trapInd},3)
                        prevTpCell=prevTPCells{trapInd}(:,:,prevCellInd);
                        unionCell=tCell | prevTpCell;
                        intersectCell=(tCell & prevTpCell); %| (~tCell & gtCell);
                        cellOverlap(prevCellInd,cellInd)=sum(intersectCell(:))/sum(unionCell(:));
                        end
                    catch
                        b=1
                    end
                end
                
                dist=cellOverlap;
                index=1;
                noLabel=ones(1,size(dist,2));
                if all(size(dist)>0)
                    for i=1:size(dist,2)
                        [val loc]=max(dist(:));
                        [row col]=ind2sub(size(dist),loc);

                        if val>cellMovementThresh
                            %cell number update
                            try
                                temp_val=trapInfom1(trapInd).cellLabel(row);
                                trapInfo(trapInd).cellLabel(1,col)=temp_val;
                                dist(:,col)=0;
                                dist(row,:)=0;
                                noLabel(col)=0;
                                index=index+1;
                            catch
                                b=1
                            end
                        end
                    end
                    
                end
                
%                 for all cells that are "new" cells to the image, update them
%                 and the maxCell value
                if ~trapInfo(trapInd).cellsPresent
                    unlabelledCellNum=0;
                elseif trapInfo(trapInd).cellsPresent && isempty(trapInfo(trapInd).cell(1).cellCenter)
                    unlabelledCellNum=0;
                    trapInfo(trapInd).cellsPresent=0;
                else
                    unlabelledCellNum=length(trapInfo(trapInd).cell)-sum(trapInfo(trapInd).cellLabel>0);
                end
 
                
                if unlabelledCellNum>0
                    locsUnlabelled=find(trapInfo(trapInd).cellLabel==0);
                    trapInfo(trapInd).cellLabel(locsUnlabelled(1:unlabelledCellNum))=trapMaxCell(trapInd)+1:trapMaxCell(trapInd)+unlabelledCellNum;
                    trapMaxCell(trapInd)=trapMaxCell(trapInd)+unlabelledCellNum;
                end
            end
            
        end
        prevTPCells=currTPCells;
        
%         trapMaxCell=zeros(1,length(cTimelapse.cTimepoint(1).trapInfo));


%         for trap=1:length(cTimelapse.cTimepoint(cTimelapse.timepointsToProcess(1)).trapInfo)
%             if timepoint==cTimelapse.timepointsToProcess(1)
%                 if trapInfo(trap).cellsPresent
%                     len=length(trapInfo(trap).cell);
%                     trapInfo(trap).cellLabel=1:length(trapInfo(trap).cell);
%                 else
%                     len=0;
%                     trapInfo(trap).cellLabel=0;
%                     trapInfo(trap).cellLabel=0;
%                     trapInfo(trap).cell(1).cellCenter=[];
%                     trapInfo(trap).cell(1).cellRadius=[];
%                 end
%                 trapMaxCell(trap)=len;            
%             elseif ~trapInfo(trap).cellsPresent
%                 trapInfo(trap).cellLabel=0;
%                 trapInfo(trap).cell(1).cellCenter=[];
%                 trapInfo(trap).cell(1).cellRadius=[];
% 
% 
%             else
%                 trapInfo(trap).cellLabel=zeros(1,length(trapInfo(trap).cell));
%                 circen=[trapInfo(trap).cell(:).cellCenter];
%                 circen=reshape(circen,2,length(circen)/2)';
%                 cirrad=[trapInfo(trap).cell(:).cellRadius]';
%                 pt2=[circen cirrad];
%                 
%                 circen=[trapInfom1(trap).cell(:).cellCenter];
%                 circen=reshape(circen,2,length(circen)/2)';
%                 cirrad=[trapInfom1(trap).cell(:).cellRadius]';
%                 pt1=[circen cirrad];
%                 
%                 if timepoint>cTimelapse.timepointsToProcess(2)
%                     circen=[trapInfom2(trap).cell(:).cellCenter];
%                     circen=reshape(circen,2,length(circen)/2)';
%                     cirrad=[trapInfom2(trap).cell(:).cellRadius]';
%                     pt3=[circen cirrad];
%                 else
%                     pt3=ones(1,3)*Inf;
%                 end
%                 pt1=double(pt1);pt2=double(pt2);pt3=double(pt3);
%                 if isempty(pt1)
%                     pt1=ones(1,3)*Inf;
%                 end
%                 if isempty(pt2) && timepoint>cTimelapse.timepointsToProcess(1)
%                     pt2=ones(1,3)*Inf;
%                 end
%                 if isempty(pt3) && timepoint>cTimelapse.timepointsToProcess(2)
%                     pt3=ones(1,3)*Inf;
%                 end
% %                 dist=pdist2(pt1,pt2,'euclidean');
%                 dist=cTimelapse.alternativeDist(pt1,pt2);
%                 
%                 index=1;
%                 noLabel=ones(1,size(dist,2));
%                 if all(size(dist)>0)
%                     for i=1:size(dist,2)
%                         [val loc]=min(dist(:));
%                         [row col]=ind2sub(size(dist),loc);
% 
%                         if val<cellMovementThresh
%                             %cell number update
%                             temp_val=trapInfom1(trap).cellLabel(row);
%                             trapInfo(trap).cellLabel(1,col)=temp_val;
%                             dist(:,col)=Inf;
%                             dist(row,:)=Inf;
%                             noLabel(col)=0;
%                             index=index+1;
%                         end
%                     end
%                     
%                 end
%                 
% %                 for all cells that are "new" cells to the image, update them
% %                 and the maxCell value
%                 if ~trapInfo(trap).cellsPresent
%                     unlabelledCellNum=0;
%                 elseif trapInfo(trap).cellsPresent && isempty(trapInfo(trap).cell(1).cellCenter)
%                     unlabelledCellNum=0;
%                     trapInfo(trap).cellsPresent=0;
%                 else
%                     unlabelledCellNum=length(trapInfo(trap).cell)-sum(trapInfo(trap).cellLabel>0);
%                 end
%  
%                 
%                 if unlabelledCellNum>0
%                     locsUnlabelled=find(trapInfo(trap).cellLabel==0);
%                     trapInfo(trap).cellLabel(locsUnlabelled(1:unlabelledCellNum))=trapMaxCell(trap)+1:trapMaxCell(trap)+unlabelledCellNum;
%                     trapMaxCell(trap)=trapMaxCell(trap)+unlabelledCellNum;
%                 end
%             end
%         end
        cTimelapse.cTimepoint(timepoint).trapInfo=trapInfo;
    end
end

cTimelapse.cTimepoint(1).trapMaxCell=trapMaxCell;
trap=1:length(cTimelapse.cTimepoint(1).trapInfo);
end
