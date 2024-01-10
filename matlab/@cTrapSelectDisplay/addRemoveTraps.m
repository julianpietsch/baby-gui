function addRemoveTraps(cDisplay)
% addRemoveTraps(cDisplay)
%
% function add and remove traps from cTimelapse. If SelectionType is 'alt'
% this is a right click and the trap is removed, if not then it is a left
% click and it is added. In both cases the identifyTrapLocationsSingleTP
% method of babyTimelapse is used to fix the location of the traps and
% produce the overlap image.

% get location of click
cp=get(cDisplay.axesHandle,'CurrentPoint');
cp=round(cp);
Cx=cp(1,1);
Cy=cp(1,2);

% switch on click type ('alt' for right click)
if strcmp(get(gcbf,'SelectionType'),'alt')
    pts=[];
    pts(:,1)=[cDisplay.trapLocations.xcenter];
    pts(:,2)=[cDisplay.trapLocations.ycenter];
    
    
    trapPt=[Cx Cy];
    D = pdist2(pts,trapPt,'euclidean');
    [~, loc]=min(D);
    
    cDisplay.trapLocations(loc)=[];
    
    % refine trapLocations and set data in cDisplay.cTimelapse
    cDisplay.trapLocations = cDisplay.cTimelapse.identifyTrapLocationsSingleTP(cDisplay.timepoint,cDisplay.trapLocations,'none',cDisplay.cc);
    
    cDisplay.setImage;
    
    % Update the log
    logmsg(cDisplay.cTimelapse,'Remove trap at %s', num2str([Cx,Cy]));
else %not right click.
    cDisplay.trapLocations(end+1).xcenter=Cx;
    cDisplay.trapLocations(end).ycenter=Cy;
    
    % refine trapLocations and set data in cDisplay.cTimelapse
    [cDisplay.trapLocations]=cDisplay.cTimelapse.identifyTrapLocationsSingleTP(cDisplay.timepoint,cDisplay.trapLocations,length(cDisplay.trapLocations),cDisplay.cc);
    
    cDisplay.setImage;
    
    % Update the log
    logmsg(cDisplay.cTimelapse,'Add trap at %s',num2str([Cx,Cy]));
end
end
