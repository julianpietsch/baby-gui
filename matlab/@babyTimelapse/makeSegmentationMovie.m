function makeSegmentationMovie(cTimelapse,channel,location,name,if_tracked,timepoints,traps)
% makeSegmentationMovie(channel,location,if_tracked,traps)
% export a whole pile of png's that can be fused into a movie of the cell
% segmentation result. If numerous traps given, makes a big fused stack.

if nargin<2 || isempty(channel)
    channel = 1;
end
channel = channel(1);
if nargin<3 || isempty(location)
    location = uigetdir([],'please select a location for storing the exported images');
end

if nargin<4 || isempty(name)
    name = inputdlg('please provide an experiment name to be used in the files','position name',1,{'name'});
    name = name{1};
end

if nargin<5 || isempty(if_tracked)
    if_tracked = true;
end

if nargin<6 || isempty(timepoints)
    timepoints = cTimelapse.timepointsToProcess;
end
if nargin<7 || isempty(traps)
    traps = cTimelapse.defaultTrapIndices;
end

rng(1);

max_all_cell_labels = max(cTimelapse.returnMaxCellLabel);


cmap = hsv(max_all_cell_labels);
cmap = cmap(randperm(max_all_cell_labels),:);

% so as not to disrupt random calls in the future.
rng('shuffle')

for TPi = 1:length(timepoints)
    TP = timepoints(TPi);
    
    trap_image_stack = cTimelapse.returnTrapsTimepoint(traps,TP,channel);
    seg_res_image_stack = cTimelapse.returnTrapsSegResTimepoint(traps,TP);
    
    % mage tiled images
    mega_trap_image = MakeMegaImage(trap_image_stack);
    mega_trap_image = SwainImageTransforms.min_max_normalise(mega_trap_image);
    mega_seg_res_image = MakeMegaImage(seg_res_image_stack);
    
    % add colours on outlines
    mega_trap_image = repmat(mega_trap_image,[1,1,3]);
    colour_seg_res_image = colour_image_from_indices(mega_seg_res_image,cmap);
    mega_seg_res_image =repmat(mega_seg_res_image,[1,1,3]);
    mega_trap_image(mega_seg_res_image>0) = colour_seg_res_image(mega_seg_res_image>0);
    
    mega_trap_image = uint8(mega_trap_image*254);
    imwrite(mega_trap_image,fullfile(location,sprintf([name 'ch%2.2d_tp%6.6d.png'],channel,TP)));
    
end



end

function colour_image = colour_image_from_indices(index_image,colour_map)

colour_image = zeros(numel(index_image),3);
colour_image(index_image(:)>0,:) = colour_map(index_image(index_image(:)>0),:);
colour_image = reshape(colour_image,size(index_image,1),size(index_image,2),[]);
end

