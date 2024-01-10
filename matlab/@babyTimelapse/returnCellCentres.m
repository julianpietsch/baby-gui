function [ cell_centres ] = returnCellCentres( cTimelapse,trap_index,timepoint,cell_index )
% [ cell_centres ] = returnCellCentres( cTimelapse,trap_index,timepoint,cell_index )
%
%returns the RELATIVE(relative to the trap centre) position of
%the cells in the image. can handle an CellIndex as an array,
%in which case returns column of form [x's   y's].

if nargin<4 || isempty(cell_index)
    cell_index = 1:length(cTimelapse.cTimepoint(timepoint).trapInfo(trap_index).cell);
end

cell_centres =  reshape([cTimelapse.cTimepoint(timepoint).trapInfo(trap_index).cell(cell_index).cellCenter],2,[])';


end

