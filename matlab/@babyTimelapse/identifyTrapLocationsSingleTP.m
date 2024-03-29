function [trapLocations] = identifyTrapLocationsSingleTP(cTimelapse,timepoint,trapLocations,trapLocationsToCheck,trap_prediction_image)
% [trapLocations] = identifyTrapLocationsSingleTP(cTimelapse,timepoint,trapLocations,trapLocationsToCheck,trap_prediction_image)
%
% cTimelapse                -   babyTimelapse object.
% timepoint                 -   timepoint at which to identify traps
% trapLocations             -   structure array with fields xcenter,ycenter
%                               holding the locations of traps. If this is
%                               empty, it is produced from
%                               trap_prediction_image.
% trapLocationsToCheck      -   array of indices of traps that should be
%                               shifted to their nearest maximum in the
%                               trap_prediction_image. If 'none', this is
%                               not performed.
% trap_prediction_image     -   image of the size of BABYTIMELAPSE.IMSIZE
%                               with high values at likely trap centres.
%                               If not provided, it is generated by
%                               BABYTIMELAPSE.GENERATETRAPLOCATIONSPREDICTIONIMAGE
%
% trapLocations         -   structure array with fields xcenter,ycenter
%                           holding the locations of traps.
%
% Identifies trapLocations from the trap_prediction_image and uses it to
% populate cTimelapse.cTimepoint(timepoint).trapLocations - a structure
% array of trapLocations. If trapLocationsToCheck is specified, these
% locations are shifted to their nearest maximum in the
% trap_prediction_image. This is useful in GUI trap selection. 
%
% See also, CTRAPSELECTDISPLAY ,
% BABYTIMELAPSE.GENERATETRAPLOCATIONSPREDICTIONIMAGE.



if nargin<4 || isempty(trapLocationsToCheck)
    trapLocationsToCheck='none'; %traps to put through the 'find nearest best point and set trap location to that' mill. if string 'none' does none of them.
end

if nargin<5 || isempty(trap_prediction_image)
    channel = cTimelapse.trapTemplateChannel;
    trap_prediction_image = generateTrapLocationPredictionImage(cTimelapse,timepoint,channel);
end

cTrap = cTimelapse.cTrapSize;

if nargin<3 || isempty(trapLocations)
    [trapLocations] = predictTrapLocations(cTrap,trap_prediction_image); 
end

[trapLocations]=updateTrapLocations(cTrap,trapLocations,trapLocationsToCheck,trap_prediction_image);

cTimelapse.cTimepoint(timepoint).trapLocations=trapLocations;

trapInfo = cTimelapse.trapInfoTemplate;

cTimelapse.cTimepoint(timepoint).trapInfo=trapInfo ;
cTimelapse.cTimepoint(timepoint).trapInfo(1:length(trapLocations))= trapInfo;
end

function [trapLocations]=predictTrapLocations(cTrap,trap_prediction_image)
% [trapLocations]=predictTrapLocations(cTrap,trap_prediction)
%
% cTrap                     -   cTrap structure with size of trap image.
% trap_prediction_image     -   image of the size of cTimelapse.imSize.
%                               Should be high at locations likely to be
%                               the centre of the trap and low otherwise.
%
% trapLocations             -   structure array as found in
%                               cTimelapse.cTimepoint (fields are xcenter,ycenter)
%
% iteratively goes through pixels in trap_prediction_image, finding the
% maximum, making it a trapLocation, and ruling out the area within 1.2
% trapwidths of it. 
% Applies a threshold of either 0.05 or 5* the standard deviation -
% whichever is higher.

[max_dynamic, imax] = max(trap_prediction_image(:));

val_thresh = 0.05;

% Start by assuming the image_thresh is the 95th percentile. We will assume
% that we have at least two traps in the image. We will update the image
% threshold based on the variability in the maxima.
image_thresh = quantile(trap_prediction_image(:),0.95);

val_thresh = max(val_thresh,image_thresh);

trap_index=1;

trapLocations = struct('xcenter',[],'ycenter',[]);

ccBounding=1.2;
bY=floor(cTrap.bb_height*ccBounding);
bX=floor(cTrap.bb_width*ccBounding);
    
blotting_image = zeros(2*[bY,bX] + 1);
allmaxima = nan(1,40);
allmaxima(trap_index) = max_dynamic;

while max_dynamic> val_thresh
    [ypeak, xpeak] = ind2sub(size(trap_prediction_image),imax(1));
    trap_prediction_image = BABYutil.PutSubStack(...
        trap_prediction_image,[ypeak,xpeak],{blotting_image});
    
    trapLocations(trap_index).xcenter = xpeak;
    trapLocations(trap_index).ycenter = ypeak;
    
    if sum(trap_prediction_image(:)>0) < numel(blotting_image)
        % If there are fewer non-zero pixels left than the area of the
        % blotting images, we should cancel
        break
    end
    trap_index=trap_index+1;
    [max_dynamic, imax] = max(trap_prediction_image(:));
    
    % Update threshold
    if numel(allmaxima) < trap_index
        allmaxima(end+1:end+40) = NaN;
    end
    allmaxima(trap_index) = max_dynamic;
    val_thresh = max(nanmedian(allmaxima) - 2*iqr(allmaxima),image_thresh);
end

% no traps found
if trap_index==1
    trapLocations = [];
end

end


function [trapLocations] = updateTrapLocations(cTrap,trapLocations,trapLocationsToCheck,trap_prediction_image)
% [trapLocations] = updateTrapLocations(image,cTrap,trapLocations,trapLocationsToCheck,trap_prediction_image)
% 
% shifts those locations marked as 'ToCheck' to their nearest maximum in
% the trap_prediction_image.
%
% cTrap                     -   cTrap structure with size of trap image.
% trapLocations             -   structure array of trapLocation as stored
%                               by BABYTIMELAPSE
% trapLocationsToCheck      -   array of indices of traps that should be
%                               shifted to their nearest maximum.
% trap_prediction_image     -   image of the size of BABYTIMELAPSE.IMSIZE
%                               with high values at likely trap centres. 

if nargin<4 || isempty(trapLocationsToCheck)
    trapLocations = 1:length(trapLocationsToCheck);
elseif strcmp(trapLocationsToCheck,'none')
    trapLocationsToCheck = [];
end

sub_image_extent = round([cTrap.bb_height/3, cTrap.bb_width/3]);


for i=1:length(trapLocations)

    xcurrent=trapLocations(i).xcenter;
    ycurrent=trapLocations(i).ycenter;
    
    if ismember(i,trapLocationsToCheck) 
        %if this is one of the trap locations to check, make it's location
        %the maximum cross correlation value within a third of a trap width
        %of the point selected for the trap.
        temp_im = BABYutil.GetSubStack(trap_prediction_image,...
            [ycurrent,xcurrent],2*sub_image_extent +1,-Inf);
        temp_im = temp_im{1};
        
        [~, maxloc]=max(temp_im(:));
        [ypeak, xpeak] = ind2sub(size(temp_im),maxloc);

        xcenter=(xcurrent+xpeak-sub_image_extent(2)-1);
        ycenter=(ycurrent+ypeak-sub_image_extent(1)-1);

        trapLocations(i).xcenter=xcenter;
        trapLocations(i).ycenter=ycenter;

    end
    
end


end


