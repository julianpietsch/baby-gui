function [px,py] = cartesian_spline_from_radii(radii,angles,center,image_size)
%function [px,py] = cartesian_spline_from_radii(radii,angles,center,image_size)

% function to take a set of radii,angles,a center and return an unbroken
% edge of the cell with no repeats.

% radii        -   vector of radii around the cell
% angles       -   angles to the x axis at which these radii are given (clockwise
%                  is positive)
% image_size   -   size of the image in which the points should be confined.

% px           -   x coordinates of resultant end points.
% py           -   y coordinates of resultant end points.

%NB: no longer ordering radii and vectors by angle. We now assume that the
%ordering is intentional

% Make the boundaries periodic and use a parametric variable
angles_loop = [angles(:); angles(1)];
radii_loop = [radii(:); radii(1)];
t_ = linspace(0,2*pi,numel(angles_loop))';
%construct x and y splines using file exchange function 'splinefit'
x_spline = splinefit(t_,radii_loop.*cos(angles_loop),t_,'p');
y_spline = splinefit(t_,radii_loop.*sin(angles_loop),t_,'p');

% Estimate required sampling density from lengths of piecewise linear
% segments
xy_loop = [radii_loop.*cos(angles_loop), radii_loop.*sin(angles_loop)];
linperim = sum(sqrt(sum(diff(xy_loop).^2,2)));
n_steps = round(2.5*linperim);
steps = linspace(0,2*pi,n_steps);

%Calculate spline for the dense parametric representation
px = round(center(1)+ppval(x_spline,steps));
py = round(center(2)+ppval(y_spline,steps));

%check they are sensible
px(px<1) = 1;
px(px>image_size(2)) = image_size(2);

py(py<1) = 1;
py(py>image_size(1)) = image_size(1);

%remove repeats (i.e. pixels that do not differ from their neighbours)
pxpy = unique([px(:),py(:)],'rows','stable');
px = pxpy(:,1);
py = pxpy(:,2);

end