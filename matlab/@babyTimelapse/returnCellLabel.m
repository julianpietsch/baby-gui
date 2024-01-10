function cell_label = returnCellLabel(cTimelapse,trap_index,time_point,cell_indices)
% cell_label = returnCellLabel(cTimelapse,trap_index,time_point,cell_indices)
%
% return the labels of cells at timepoint time_point and trap index
% trap_index.
% If cell_indices is empty/absent it returns all.
cell_label = cTimelapse.cTimepoint(time_point).trapInfo(trap_index).cellLabel;
   
if nargin>=3 && ~isempty(cell_indices)
    
    cell_label = cell_label(cell_indices);
    
end

end


