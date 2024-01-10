function image = min_max_normalise( image )
% image = min_max_normalise( image )
% subtracts the minimum, divides by the new maximum (if it is non zero);

image = image-min(image(:));
m = max(image(:));
if m>0;
    image = image/m;
end
end

