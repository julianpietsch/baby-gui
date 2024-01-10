function LogicalPoints = px_py_to_logical( px,py,ImageSize )
% LogicalPoints = px_py_to_logical( px,py,ImageSize )
% 
% takes set of points given by their x,y location(x is positive in right
% direction, y is positive in down direction) and makes a logical of size
% image size with true at each location (px(i),py(i)). Image size is given
% in the matrix convention [height,width]

        LogicalPoints = false(ImageSize);
        LogicalPoints(py + (px-1)*ImageSize(1)) = true;
        
end

