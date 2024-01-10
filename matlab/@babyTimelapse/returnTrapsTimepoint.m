function trapsTimepoint=returnTrapsTimepoint(cTimelapse,traps,timepoint,channel,type)
%trapsTimepoint=returnTrapsTimepoint(cTimelapse,traps,timepoint,channel,type)
%
% If there are traps in the timelapse, this returns a zstack containing the
% images of the set of traps indicated at the timepoint indicated. If there
% are no traps in the timelapse however, it return the entire frame in a
% single 2D image.
%
% pads image with mean value, so pixels outside the image that are still in
% the trap region will have a value of the mean of the whole image.
%
% if the requested channel is not present at the requested timepoint then
% an image stack of the appropriate size but all zeros is returned.
%
% should really use the returnTrapsFromImage - too lazy.

if nargin<3
    timepoint=1;
end

if nargin<2||isempty(traps)
    traps=1:length(cTimelapse.cTimepoint(timepoint).trapInfo);
end



if nargin<4
    channel=1;
end

if nargin<5
    type='max';
end

if cTimelapse.trapsPresent 
    
    cTrap=cTimelapse.cTrapSize;
    image=cTimelapse.returnSingleTimepoint(timepoint,channel,type);
    bb=max([cTrap.bb_width cTrap.bb_height])+100;
    bb_image=padarray(image,[bb bb],1.2*mean(image(:)));
    
    trapsTimepoint=zeros(2*cTrap.bb_height+1,2*cTrap.bb_width+1,length(traps),'like',image);
    for j=1:length(traps)
        y=round(cTimelapse.cTimepoint(timepoint).trapLocations(traps(j)).ycenter + bb);
        x=round(cTimelapse.cTimepoint(timepoint).trapLocations(traps(j)).xcenter + bb);
        temp_im=bb_image(y-cTrap.bb_height:y+cTrap.bb_height,x-cTrap.bb_width:x+cTrap.bb_width);
        trapsTimepoint(:,:,j)=temp_im;
        
    end
    
else
    trapsTimepoint=cTimelapse.returnSingleTimepoint(timepoint,channel,type);
end

