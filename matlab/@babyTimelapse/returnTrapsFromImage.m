function trapsTimepoint=returnTrapsFromImage(cTimelapse,image,timepoint,traps,padding_value)
% trapsTimepoint=returnTrapsFromImage(cTimelapse,image,timepoint,traps)
%
% fairly niche function. Takes an image that is either a strip of trap
% sized images or a full image the size of imSize and extracts the traps
% from it. Checks for one or the other by comparing size(image,1) to
% cTimelapse.cTrap.bb_height*2 + 1
%
% cTimelapse    :   object of the babyTimelapse class
% image         :   an image, either a full image of size imSize or a strip
%                   of trap images.
% timepoint     :   the timepoint at which to extract trap images. Only
%                   used if using the 'full image' mode
% traps         :   the indices of the traps to extract. If in 'full image'
%                   mode this the trap positions are used to make the
%                   extraction. If in 'strip image' mode then they are
%                   counted along the width from the left hand side.
%                   Defaults to all traps if empty.
% padding_value :   Used to pad output when it moves out of image.
%                   defaults to falsefor logical and image mean for all
%                   others.
%
% if there are not traps present in this timelapse it just gives the image
% back.

if nargin<4
    traps=1:length(cTimelapse.cTimepoint(timepoint).trapInfo);
end

if nargin<5 || isempty(padding_value)
    if islogical(image)
        padding_value=0;
    else
        padding_value=mean(image(:));
    end
end

if cTimelapse.trapsPresent
    
    cTrap=cTimelapse.cTrapSize;
    bb=max([cTrap.bb_width cTrap.bb_height])+100;

    bb_image=padarray(image,[bb bb],padding_value);
    %if the traps have been converted to be flat in a single image
    if size(image,1)==(cTrap.bb_height*2+1)
        trapsTimepoint=zeros(2*cTrap.bb_height+1,2*cTrap.bb_width+1,length(traps),'double');
        for j=1:length(traps)
            x=(j-1)*(cTrap.bb_width*2+1)+bb;
            temp_im=bb_image(1+bb:bb+2*cTrap.bb_height+1,x+1:x+2*cTrap.bb_width+1);
            trapsTimepoint(:,:,j)=temp_im;
        end
        %else, if the image is the size of the image acquired by the camera
    elseif isequal(size(image),cTimelapse.imSize)
        trapsTimepoint=zeros(2*cTrap.bb_height+1,2*cTrap.bb_width+1,length(traps),'double');
        for j=1:length(traps)
            y=round(cTimelapse.cTimepoint(timepoint).trapLocations(traps(j)).ycenter + bb);
            x=round(cTimelapse.cTimepoint(timepoint).trapLocations(traps(j)).xcenter + bb);
            temp_im=bb_image(y-cTrap.bb_height:y+cTrap.bb_height,x-cTrap.bb_width:x+cTrap.bb_width);
            trapsTimepoint(:,:,j)=temp_im;
        end
    else
        disp(size(image))
        disp(cTimelapse.imSize)
        disp(cTimelapse.cTrapSize)
        error('image must be either a horizontal strip of single trap images or a full size image')
    end
    
    
else
    trapsTimepoint = image;
end
