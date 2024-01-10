function Stack =  PutSubStack( Stack,Centres,SubStacks)
% [ SubStackCell ] = PutSubStack( Stack,Centres,SubStacks )
% puts SubStack{i} at in Stack at location Centres(i,:).
% centres must have the appropriate dimension for stack.



for centrei = 1:size(Centres,1)
    
    StackSize = size(Stack);
    if length(SubStacks)==1
        SubStack = SubStacks{1};
    else
        SubStack = SubStacks{centrei};
    end
    
    SizeOfSubStack = size(SubStack);
    
    if length(SizeOfSubStack)<length(StackSize)
        
        SizeOfSubStack((length(SizeOfSubStack)+1):length(StackSize)) = 1;
        
    end
    
    extent = ceil(SizeOfSubStack/2);
    range = cell(length(SizeOfSubStack),1);
    corresponding_range = range;
    
    for dimensioni =1:length(SizeOfSubStack)
        
        range{dimensioni} = Centres(centrei,dimensioni) - extent(dimensioni)  + (1:SizeOfSubStack(dimensioni));
        %centres for odd sizes
        %for even sizes(not recommended) stack is up down from 'true'
        %fractional centre.

        
        corresponding_range{dimensioni} = 1:SizeOfSubStack(dimensioni);
        
        good_values = range{dimensioni}>=1 & range{dimensioni}<=StackSize(dimensioni);
        
        range{dimensioni} = range{dimensioni}(good_values);
        corresponding_range{dimensioni} = corresponding_range{dimensioni}(good_values);
    end
    
    if length(StackSize)==3
        Stack(range{1},range{2},range{3}) = ...
            SubStack(corresponding_range{1},corresponding_range{2},corresponding_range{3});
    elseif length(StackSize) ==2
         Stack(range{1},range{2}) = ...
            SubStack(corresponding_range{1},corresponding_range{2});
    end
    
end

end
