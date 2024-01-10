function trap_image = returnWholeTrapImage(cTimelapse,trapOutline,timepoint,traps)
% trapOutline = returnWholeTrapImage(cTimelapse,trapOutline,timepoint,traps(optional))
%
% returns an image (double) of the size cTimelapse.imSize with a trap outline at
% every trap location. trapOutline must be either an nxmxz array of trap
% outlines (1 for each trap) or a single trap outline.
%
% WARNING: currently this only gives trap pixels for the traps being
% tracked, so any non-tracked traps will not be blotted out.

cTrap=cTimelapse.cTrapSize;
bb=max([cTrap.bb_width cTrap.bb_height])+100;
trap_image=zeros(cTimelapse.imSize + bb);

if nargin<4 || isempty(traps)
    traps = 1:length(cTimelapse.cTimepoint(timepoint).trapInfo);
end

if ~ismember(size(trapOutline,3),[1,length(traps)])
    error('size(trapOutline,3) should be length(traps) or 1')
end

for trapi=traps
    trap = traps(trapi);
    y=round(cTimelapse.cTimepoint(timepoint).trapLocations(trap).ycenter + bb);
    x=round(cTimelapse.cTimepoint(timepoint).trapLocations(trap).xcenter + bb);
    
    % put the trap in the appropriate place in the image. min in
    % trapOutline index ensures that if only 1 trapOutline is provided it
    % is used for all the traps.
    trap_image(y-cTrap.bb_height:y+cTrap.bb_height,x-cTrap.bb_width:x+cTrap.bb_width) = trapOutline(:,:,min([length(traps),trapi]));
end
trap_image = trap_image((1:cTimelapse.imSize(1))+bb,(1:cTimelapse.imSize(2))+bb);

end
