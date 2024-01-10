function show_image = get_cell_image(image,size_subimage,centerStack,m)
% function show_image = get_cell_image(image,size_subimage,centerStack,m)
%
% get a stack of size_subimage by size_sub_image chunks of the image
% centered on centers in the centre vector.
% size_subimage should be odd and can be a [vertical_size horiontal_size]
% sub image.
% m is the value with which to padarray if necessary - defaults to false or
% image median depending on class of image.
% centerStack = [x's  y's]
%
% See also, BABYUTIL.GETSUBSTACK

if numel(size_subimage) ==1
    size_subimage = [1 1] * size_subimage;
end

if nargin<4 || isempty(m)
    info = whos('image');
    if strcmp(info.class,'logical')
        m = false;
        show_image = false(size_subimage(1),size_subimage(2),size(centerStack,1));
    else
        m = median(image(:));
        show_image = zeros(size_subimage(1),size_subimage(2),size(centerStack,1));
    end
end

image_cell = BABYutil.GetSubStack(image,round(fliplr(centerStack)),size_subimage,m);


if size(centerStack,1)>1
    show_image = m*ones([size_subimage size(centerStack,1)]);
    
    for i=1:length(image_cell)
        
        show_image(:,:,i) = image_cell{i};
        
    end
else
    show_image = image_cell{1};
end

%
% image = padarray(image,((size_subimage-1)/2),m);
%
%
% %gets 30 by 30 square centered on 'center' in the original image
% for i=1:size(centerStack,1)
%
%     show_image(:,:,i) = image(round(centerStack(i,2))+(0:(size_subimage(1)-1))',round(centerStack(i,1))+(0:(size_subimage(2)-1))');
%
%
% end

end
