function importResultFromImage(cTimelapse,segmentation_mask,tp,dilate_length)
% importResultFromImage(cTimelapse,segmentation_mask,tp)
% takes a segmentation_mask and puts the data into the cTimelapse
% data_structure. Useful for comparing with other softwares or integrating
% them into our data.
% 
% segmentation_mask     - segmentation to be imported. Can be either an
%                         nxnx1 array in which case each area with a given
%                         number is taken to be a cell. NOTE, this means if
%                         they are all labelled 1 you should run bwlabel
%                         first. The number on the cell is taken to be it's
%                         tracking number.
% tp                    - the time point to which the segmentation_mask
%                         corresponds
% dilate_length         - dilate cell outlines by this amount.


trap_seg_mask_stack = cTimelapse.returnTrapsFromImage(segmentation_mask,tp,cTimelapse.defaultTrapIndices(tp),0);
trap_indices = cTimelapse.defaultTrapIndices(tp);

for TI = 1:length(trap_indices)
    trap = trap_indices(TI);
    trap_seg_mask = trap_seg_mask_stack(:,:,TI);
    trap_outline = false(size(trap_seg_mask,1),size(trap_seg_mask,2));
    importResultFromTrapImage(cTimelapse,trap_seg_mask,trap,tp,trap_outline,dilate_length);

end

cTimelapse.cTimepoint(1).trapMaxCell = cTimelapse.returnMaxCellLabel;
end



