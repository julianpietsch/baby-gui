function [px,py] = eval_cartesian_spline_from_radii(eval_angles,radii,angles,center)
%function [px,py] = eval_cartesian_spline_from_radii(newangles,radii,angles,center)
%
% function to form a cartesian spline with knots defined by a set of radii,
% angles and a center, and evaluate this spline at arbitrary angles.
% eval_angles  -   vector of (circular) angles to evaluate the spline at
% radii        -   vector of radii around the cell
% angles       -   angles to the x axis at which these radii are given (clockwise
%                  is positive)
%
% Returns
% px           -   x coordinates of resultant end points.
% py           -   y coordinates of resultant end points.

if nargin<4 || isempty(center), center = [0,0]; end

%NB: no longer ordering radii and vectors by angle. We now assume that the
%ordering is intentional

% Make the boundaries periodic and use a parametric variable
angles_loop = [angles(:); angles(1)];
radii_loop = [radii(:); radii(1)];
t_ = linspace(0,2*pi,numel(angles_loop))';
%construct x and y splines using file exchange function 'splinefit'
x_spline = splinefit(t_,radii_loop.*cos(angles_loop),t_,'p');
y_spline = splinefit(t_,radii_loop.*sin(angles_loop),t_,'p');

% Evaluate splines at requested angles
eval_angles = mod(eval_angles,2*pi);
px = center(1)+ppval(x_spline,eval_angles);
py = center(2)+ppval(y_spline,eval_angles);
end