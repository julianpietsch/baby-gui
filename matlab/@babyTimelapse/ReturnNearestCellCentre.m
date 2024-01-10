function CellNumNearestCell = ReturnNearestCellCentre(cTimelapse,timepoint,trap,pt,thresh)
% ReturnNearestCellCentre: Get nearest cell to a point on the image
% --------------------------------------------------------------
% Method written by Elco Bakker to find the nearest cell in a trap to a given
% pt (i.e. the cell with the center nearest the point pt) and return its cell
% number (its index in the trapinfo(trap).cell array.
% 
% SYNTAX: 
% 
%   CellNumNearestCell = ReturnNearestCellCentre(cTimelapse,timepoint,trap,pt,thresh)
% 
% INPUT
% 
%   cTimelapse    -       object of cTimelape class
%   timepoint     -       integer indicating timepoint of interest
%   trap          -       integer indicating trap of interest
%   pt            -       1 by 2 double of the form [x y] indicating the
%                        point of reference to which to find the closest
%                        cell
%   thresh        -       (optional) a maximum allowed distance for a cell to
%                        be returned.
% 
% 
% OUTPUT
% 
%   CellNumNearestCell        -      index of cell closest to the point pt.
%                                    Empty if there are no cells in the trap
%                                    or if no cells are closer than thresh.


if nargin<5
    thresh = inf;
end



CellNumNearestCell = [];

if cTimelapse.cTimepoint(timepoint).trapInfo(trap).cellsPresent
    
    circen=[cTimelapse.cTimepoint(timepoint).trapInfo(trap).cell(:).cellCenter];
    circen=reshape(circen,2,length(circen)/2)';
    pts=double(circen);
    
    aPointMatrix = repmat(pt,size(pts,1),1);
    D = (sum(((aPointMatrix-pts).^2), 2)).^0.5;
    [MinVal, CellNumNearestCell]=min(D);
    MinVal = MinVal(1);
    CellNumNearestCell = CellNumNearestCell(1);
    if MinVal > thresh
        CellNumNearestCell = [];
    end
        
    

end

end
