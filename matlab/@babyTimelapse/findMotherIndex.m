function motherIndex=findMotherIndex(cTimelapse,centre_method,offset)
% motherIndex=findMotherIndex(cTimelapse,centre_method,offset)
%Identify the mother index ... cells that are closest to the center of the
%trap, but based on their location on trapInfo.cell
%1) use the distance to the "center"
%2) if closest cell is <1/4 of trap dimensions away from it, it is the
%mother
%
% cTimelapse    -   an object of the babyTimelapse class
% centre_method -   method for finding the centre of the trap. defaults to
%                   cell_centre:
%                   options:
%                       'cell_centre'  - uses the centre of the cells in the
%                                        trap.
%                       'refined_trap' - uses the centre of the refined
%                                        trap pixel outline.
% offset        -   if refined trap methd, and arbitary offset added in the
% x direction (i.e. along the length of the trap).

%identify the center of the trap by finding themode of the x and y
%locations

if cTimelapse.trapsPresent

if nargin<2 || isempty(centre_method)
    centre_method = 'cell_centre';
end


if nargin<3 || isempty(offset)
    offset = -3;
end

xloc=zeros(1,1e5);
yloc=zeros(1,1e5);
xlocM=zeros(100,length(cTimelapse.timepointsProcessed));
ylocM=zeros(100,length(cTimelapse.timepointsProcessed));

traps = 1:length(cTimelapse.cTimepoint(cTimelapse.timepointsToProcess(1)).trapInfo);

switch centre_method
    case 'cell_centre'
ind=1;
for timepoint=find(cTimelapse.timepointsProcessed)
    indM=1;
    if cTimelapse.timepointsProcessed(timepoint)
        trapInfo=cTimelapse.cTimepoint(timepoint).trapInfo;
        for trap=1:length(cTimelapse.cTimepoint(cTimelapse.timepointsToProcess(1)).trapInfo)
            if trapInfo(trap).cellsPresent
                circen=[trapInfo(trap).cell(:).cellCenter];
                circen=reshape(circen,2,length(circen)/2)';
                xloc(ind:ind+size(circen,1)-1)=circen(:,1);
                xlocM(indM:indM+size(circen,1)-1,timepoint)=circen(:,1);
                yloc(ind:ind+size(circen,1)-1)=circen(:,2);
                ylocM(indM:indM+size(circen,1)-1,timepoint)=circen(:,2);
                ind=ind+size(circen,1);
                indM=indM+size(circen,1);
            end
        end
    end
end
cellPres=xloc>0;
trapCenterX=median(xloc(cellPres));
trapCenterY=median(yloc(cellPres));
ylocM(ylocM==0)=NaN; xlocM(xlocM==0)=NaN;
trapCenterXTime=smooth(nanmedian(xlocM),size(xlocM,2)/2,'rlowess');
trapCenterYTime=smooth(nanmedian(ylocM),size(xlocM,2)/2,'rlowess');

trapCenterXTime = repmat(trapCenterXTime,1,length(traps));
trapCenterYTime = repmat(trapCenterYTime,1,length(traps));

    case 'refined_trap'
        
        trapCenterXTime = zeros(length(cTimelapse.cTimepoint),length(traps));
        trapCenterYTime = trapCenterXTime;
        for timepoint=cTimelapse.timepointsToProcess
            indM=1;
            if cTimelapse.timepointsProcessed(timepoint)
                trapInfo=cTimelapse.cTimepoint(timepoint).trapInfo;
                for trap=1:length(trapInfo)
                        
                        % find centre as the mean of the trap pixels with
                        % an arbitrary x offset of -3;
                        [I,J] = find(trapInfo(trap).refinedTrapPixelsInner);
                        trapCenterXTime(timepoint,trap) = mean(J)+offset;
                        trapCenterYTime(timepoint,trap) = mean(I);
                        
                    
                end
            end
        end

        
    otherwise
        error('invalid mother identification method')

end
%debug for old cTimelapses without the cTrapSize parameter
if isempty(cTimelapse.cTrapSize)
    cTimelapse.cTrapSize.bb_height=40;
    cTimelapse.cTrapSize.bb_width=40;
end
    % if the closest cell is within a 1/4 of the frame from the center of the
% trap, that is the mother

cutoff=ceil(cTimelapse.cTrapSize.bb_height/4);
motherIndex=[];

for timepoint=1:length(cTimelapse.timepointsProcessed)
    if cTimelapse.timepointsProcessed(timepoint)
        % Trigger the TimepointChanged event for babyLogging
        babyLogging.changeTimepoint(cTimelapse,timepoint);
        
        trapInfo=cTimelapse.cTimepoint(timepoint).trapInfo;
        
        
        for trap=traps
            if trapInfo(trap).cellsPresent
                %Below is if you want to use the median cell location as a function
                %of time, rather than just a single point. This is needed if the
                %timelapse is really long, and the cells grow a lot over time.
                pt1=[trapCenterXTime(timepoint,trap) trapCenterYTime(timepoint,trap)];
                pt1=double(pt1);
        
                circen=[trapInfo(trap).cell(:).cellCenter];
                circen=reshape(circen,2,length(circen)/2)';
                pt2=[circen];
                pt2=double(pt2);
                
                dist=pdist2(pt1,pt2);
                [val, ind]=min(dist);
                if val<cutoff
                    motherIndex(trap,timepoint)=ind;
                else
                    motherIndex(trap,timepoint)=0;
                end
            else
                motherIndex(trap,timepoint)=0;
            end
        end
    end
end


else
    motherIndex=zeros(1,length(cTimelapse.timepointsProcessed));
end
%%
% During the tracking step, deal with mothers differently


% Once cells have been labelled, go back throug hand create a list of
% mothers
