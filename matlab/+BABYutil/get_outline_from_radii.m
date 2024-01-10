function [ outline ] = get_outline_from_radii( radii,angles,center,image_size )
%[ outline ] = get_outline_from_radii( radii,angles,center,image_size )
%
% makes a logical array outline of the cell from the inputs (those used in
% active contour method).
%
% See also, ACBACKGROUNDFUNCTIONS.GET_FULL_POINTS_FROM_RADII,
% ACBACKGROUNDFUNCTIONS.PX_PY_TO_LOGICAL


[px,py] = ACBackGroundFunctions.get_full_points_from_radii(radii,angles,center,image_size);
outline = ACBackGroundFunctions.px_py_to_logical( px,py,image_size );

end

