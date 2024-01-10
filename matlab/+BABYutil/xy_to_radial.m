function [R,angle] = xy_to_radial(x,y)
%function [R,angle] = xy_to_radial(x,y)

%converts to a radius and an angle between 0 and 2pi
R = sqrt((x.^2) + (y.^2));
angle = atan(y./x);

angle = angle+(pi*(x<0));

angle = mod(angle,2*pi);


end