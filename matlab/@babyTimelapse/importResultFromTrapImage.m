function importResultFromTrapImage(cTimelapse,trap_seg_mask,trap,tp,trap_outline,dilate_length)
% importResultFromTrapImage(cTimelapse,segmentation_mask,trap,TP,trap_outline,dilate_length)
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
% trap                  - trap index to which trap_seg_mask applies
% tp                    - the time point to which the trap_seg_mask
%                         corresponds
% trap_outline          - a logical array of feature pixels. Any cell who's
%                         centre is in a True pixel will be ignored.
%                         Designed for traps but can be generally used to
%                         define a background.
% dilate_length         - dilate cells by this much before import.
%                         sometimes useful. negative =>erode.

if any(size(trap_seg_mask)~= cTimelapse.trapImSize)
    tsize =cTimelapse.trapImSize;
    error('trap seg mask of size [%d,%d]. ctimelapse expects a trap image of size [%d,%d]',size(trap_seg_mask_stack,1),size(trap_seg_mask_stack,2),tsize(1),tsize(2));
end


if nargin<5 || isempty(trap_outline) 
    trap_outline = false(size(trap_seg_mask));
end

if nargin<6 || isempty(dilate_length) || dilate_length==0
    dilate_length = 0;
else
    dilate_strel = strel('disk',abs(dilate_length));
end

trapInfo = cTimelapse.trapInfoTemplate;
cell_indices = unique(trap_seg_mask(trap_seg_mask>0));
num_cells = length(cell_indices);
% put remaining cells in trapInfo
if num_cells>0
    n = 0; %cells position in cell struct
    for CI = 1:num_cells
        cell_num = cell_indices(CI);
        cell_im = trap_seg_mask==cell_num;
        if dilate_length >0
            cell_im = imdilate(cell_im,dilate_strel,'same');
        elseif dilate_length<0
            cell_im = imerode(cell_im,dilate_strel,'same');
        end
        
        cellStruct = cTimelapse.cellInfoTemplate;
        
        [Y,X] = find(cell_im);
        
        cell_centre_Y = round(mean(Y));
        cell_centre_X = round(mean(X));
        
        % if the centre of a cell is in the trap, ignore
        if trap_outline(cell_centre_Y,cell_centre_X)>0
            continue
        end
        
        trapInfo.cellsPresent = true;
        n = n+1;
        cellStruct.cellCenter = [cell_centre_Y cell_centre_X];
        cellStruct.cellRadius = sqrt(sum(cell_im(:))/pi);
        
        segmented = false(size(cell_im));
        boundary = bwboundaries(cell_im);
        segmented(boundary{1}(:,1)+(size(segmented,1)*(boundary{1}(:,2)-1))) = true;
        cellStruct.segmented = sparse(segmented);
        
        trapInfo.cell(n) = cellStruct;
        trapInfo.cellLabel(n) = cell_num;
    end
end

cTimelapse.cTimepoint(tp).trapInfo(trap) = trapInfo;

end
