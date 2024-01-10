function [px,py] = get_points_from_radii(radii,angles,center,point_number,image_size)
%function [px,py] = get_points_from_radii(radii,angles,center,point_number,image_size)

% function to take a set of radii,angles,a center and return an unbroken
% edge of the cell with no repeats.

% radii        -   vector of radii around the cell
% angles       -   angles to the x axis at which these radii are given (clockwise
%                  is positive)
% point_number -   number of points to find and return.
% image_size   -   size of the image in which the points should be confined.

% px           -   x coordinates of resultant end points.
% py           -   y coordinates of resultant end points.


if point_number<2
    error('need at least 2 radial points')
end


%order the angles vector (may not be necessary)
[angles,indices_angles] = sort(angles,1);
radii = radii(indices_angles);

%construct spline using file exchange function 'splinefit'
r_spline = splinefit([angles; 2*pi+angles(1)],[radii;radii(1)],[angles; 2*pi+angles(1)],'p');

steps = linspace(0,2*pi,point_number+1)';
steps = steps(1:point_number,1);
radii_full = ppval(r_spline,steps);

%convert radial coords to x y coords
px = round(center(1)+radii_full.*cos(steps));%radial cords
py = round(center(2)+radii_full.*sin(steps));

%check they are sensible
px(px<1) = 1;
px(px>image_size(2)) = image_size(2);

py(py<1) = 1;
py(py>image_size(1)) = image_size(1);



end