function [ SubStackCell ] = GetSubStack( Stack,Centres,SizeOfSubStack,varargin )
% [ SubStackCell ] = GetSubStack( Stack,Centres,SizeOfSubStack,varargin )
% Returns sections of a NxMxL matrix defined by Centres and SizeOfSubStack.
% Centres given as matrix dimensions [Y X Z]
% intend it to be able to do this for various different cases but let's stick
% with just the simplest for now. 
% If the size of substack is not odd then the centre will be
% taken to be 1 less in every direction.

SubStackCell = {};
StackSize = size(Stack);
if length(SizeOfSubStack)<length(StackSize)
    
    SizeOfSubStack((length(SizeOfSubStack)+1):length(StackSize)) = StackSize((length(SizeOfSubStack)+1):length(StackSize));
    
end

if nargin<4
    PadImage = mean(Stack(:));
else
    PadImage = varargin{1};
end


    %Stack = padarray(Stack,SizeOfSubStack,PadImage);

    extent = ceil(SizeOfSubStack/2);
    range = cell(length(SizeOfSubStack),1);
    corresponding_range = range;
    
    for centrei = 1:size(Centres,1)
        

        for dimensioni =1:length(SizeOfSubStack)
            if dimensioni<=size(Centres,2)
                range{dimensioni} = Centres(centrei,dimensioni) - extent(dimensioni)  + (1:SizeOfSubStack(dimensioni));
                %centres for odd sizes
                %for even sizes(not recommended) stack is up down from 'true'
                %fractional centre.
            else
                
                range{dimensioni} = ceil(StackSize(dimensioni)/2) - extent(dimensioni) + (1:SizeOfSubStack(dimensioni)) ;
            end
                
            corresponding_range{dimensioni} = 1:SizeOfSubStack(dimensioni);
        
            good_values = range{dimensioni}>=1 & range{dimensioni}<=StackSize(dimensioni);
            
            range{dimensioni} = range{dimensioni}(good_values);
            corresponding_range{dimensioni} = corresponding_range{dimensioni}(good_values);
        end
        
        if length(StackSize)==3
            SubStackCell{centrei} = PadImage*ones(SizeOfSubStack);
            SubStackCell{centrei}(corresponding_range{1},corresponding_range{2},corresponding_range{3}) = Stack(range{1},range{2},range{3});
        elseif length(StackSize) ==2
            SubStackCell{centrei} = PadImage*ones(SizeOfSubStack);
            SubStackCell{centrei}(corresponding_range{1},corresponding_range{2}) = Stack(range{1},range{2});
        end
        
    end
    
end
