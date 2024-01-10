function img = applyStandardImageTransformation(cTimelapse,img,channel)

%used for padding data
medVal=mean(img(:));

stack_depth = size(img,3);

% apply background correction. 
% This should really be called 'flat field correction', since it .multiplies
% the image. If backgroundOffset is not empty for this channel, this is
% first subtracted, the image is then .multiplied by BackgroundCorrection,
% and backgroundOffset is then added back. This is to preserve
% combarability with images where the background is not subtracted. 
if isprop(cTimelapse,'BackgroundCorrection') && size(cTimelapse.BackgroundCorrection,2)>=channel && ~isempty(cTimelapse.BackgroundCorrection{channel})
    bgdScaling = cTimelapse.BackgroundCorrection{channel}(:,:,ones(size(img,3),1));
    if isprop(cTimelapse,'BackgroundOffset') && length(cTimelapse.BackgroundOffset)>=channel && ~isempty(cTimelapse.BackgroundOffset{channel})
        bgdOffset = cTimelapse.BackgroundOffset{channel};
        img = bgdOffset + (img - bgdOffset).*bgdScaling;
    else
        img = img.*bgdScaling;
    end
end



% if scaledImSize (the size of the final image before rotation) and
% rawImSize (the size of the loaded image) are different, then rescale
if any(cTimelapse.scaledImSize ~= cTimelapse.rawImSize)
    new_im = zeros([cTimelapse.scaledImSize stack_depth]);
    for si = 1:stack_depth
        new_im(:,:,si) = imresize(img,cTimelapse.scaledImSize);
    end
    img = new_im;
    clear new_im
end


if isequal(cTimelapse.image_flipud,true)
    img = flipud(img);
end

% rotate image. It is first padded to try and prevent zeros occuring in the
% final image as an artefact of the padding.
if cTimelapse.image_rotation~=0
    bbN = ceil(0.5*max(cTimelapse.scaledImSize)); tIm=[];
    for slicei = 1:stack_depth
        tpImtemp = padarray(img(:,:,slicei),[bbN bbN],medVal,'both');
        tpImtemp = imrotate(tpImtemp,cTimelapse.image_rotation,'bilinear','loose');
        tIm(:,:,slicei) = tpImtemp(bbN+1:end-bbN,bbN+1:end-bbN);
        if slicei==1 && stack_depth>1
            tIm(:,:,2:stack_depth) = 0;
        end
    end
    img = tIm;
    clear tIm
    clear tpImtemp
end

% shift image by 'offset' to make is align with other channels.
if isprop(cTimelapse,'offset') && size(cTimelapse.offset,1)>=channel && any(cTimelapse.offset(channel,:)~=0)
    boundaries = fliplr(cTimelapse.offset(channel,:));
    lower_boundaries = abs(boundaries) + boundaries +1;
    upper_boundaries = [size(img,1) size(img,2)] + boundaries + abs(boundaries);
    img = padarray(img,[abs(boundaries) 0],medVal);
    img = img(lower_boundaries(1):upper_boundaries(1),lower_boundaries(2):upper_boundaries(2),:);
end

% populate this if 
if isempty(cTimelapse.imSize)
    cTimelapse.imSize = [size(img,1),size(img,2)];
end

end