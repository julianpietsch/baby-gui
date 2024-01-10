function newRadii = get_radii_from_radii(newAngles,radii,angles)
%GET_RADII_FROM_RADII Return interpolated radii from a radial spline
%
%   NEWRADII = GET_RADII_FROM_RADII(NEWANGLES,RADII,ANGLES) returns radii
%   NEWRADII at the angles NEWANGLES obtained by interpolating a radial 
%   spline defined by RADII and ANGLES.

if numel(angles) ~= numel(radii)
    error('number of angles and radii defining the radial spline must match');
end

angles = angles(:); radii = radii(:);

%order the angles vector (may not be necessary)
[angles,indices_angles] = sort(angles,1);
radii = radii(indices_angles);

%construct spline using file exchange function 'splinefit'
r_spline = splinefit([angles; 2*pi+angles(1)],[radii;radii(1)],[angles; 2*pi+angles(1)],'p');
% Convert angles to domain of fitted spline knots
newAngles = angles(1)+mod(newAngles-angles(1),2*pi);
newRadii = ppval(r_spline,newAngles(:));

end
