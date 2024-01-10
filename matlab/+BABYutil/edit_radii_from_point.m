function [radii] = edit_radii_from_point(point,center,radii,angles)
% function [radii] = edit_radii_from_point(center,radii,angles,in)
% provide a point selected in the image, and all the relevant info about
% the cell and it updates the radii.
%expets an

if size(angles,1) <size(angles,2)
    angles = angles';
end

xnew = point(1)-center(1);
ynew = point(2)-center(2);
[Rnew,angle_new] = ACBackGroundFunctions.xy_to_radial(xnew,ynew);

[~,minindex] = min(abs([angles;(2*pi)] - angle_new));

if minindex==(length(radii)+1)
    minindex=1;
end

radii(minindex) = Rnew;

end