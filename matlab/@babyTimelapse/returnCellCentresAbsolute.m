function [cell_centres] = returnCellCentresAbsolute( cTimelapse,trap_index,timepoint,cell_index )
% cell_centres] = returnCellCentresAbsolute( cTimelapse,trap_index,timepoint,cell_index )
%
%returns the ABSOLUTE position (as double) of the cells in the image.
%can handle an CellIndex as an array, in which case returns
%column of form [x's   y's].

if nargin<4 || isempty(cell_index)
    cell_index = 1:length(cTimelapse.cTimepoint(timepoint).trapInfo(trap_index).cell);
end

cell_centres =  cTimelapse.returnCellCentres(trap_index,timepoint,cell_index );

if cTimelapse.trapsPresent && ~isempty(cell_centres)
    
    cell_centres = cell_centres + ...
        repmat(cTimelapse.returnTrapCentre(trap_index,timepoint) - [cTimelapse.cTrapSize.bb_width cTimelapse.cTrapSize.bb_height],length(cell_index),1) - 1;
    % the -1 on this might seem strange but think about it. if the
    % cell center relative is [1 1] it should be at the first square
    % [1 1] square of the trap image so then you should add [0 0] to
    % the first entry of the trap image which is
    %    trap_centre - bb_width
    % so then for [1 1] the answer [1 1] + trap_centre - bb_width - [1 1]
    
end

cell_centres = double(cell_centres);



end



