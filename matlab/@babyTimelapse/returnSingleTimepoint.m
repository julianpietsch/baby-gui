function timepointIm=returnSingleTimepoint(cTimelapse,timepoint,channel,type)
% timepointIm=returnSingleTimepoint(cTimelapse,timepoint,channel,type)
%
% Returns the image for a particular channel at a particular timepoint,
% applying some corrections defined by BABYTIMELAPSE properties.
% 
% timepoint   :   number indicating the desired timepoint (will access this
%                 element of the cTimepoint array).
% channel     :   number (default 1) indicating which of the channels in
%                 cTimelapse.channelNames to return.
%
% type        :   string (default 'max'). if more than one file is loaded
%                 (e.g. in the case of a z stack) how to handle this stack.
%                 either 'max','sum','min' or 'stack' - applying
%                 concomitant operation to the stack or returning the whole
%                 stack in the case of 'stack'. This is done after the
%                 image is converted to a double.
%
% Loads the image (or images in the case of a 'stack' channel) from file
% using the RETURNSINGLETIMEPOINTRAW method, which returns a stack of
% images.
%
% This stack is then treated accoding to the 'type' argument, with
% either a projection ('min','max' etc.) being made or the whole stack
% being returned ('stack'). The default is max.
% A number of operations are then applied to the image in the following
% order.
%
% background correction - if the field cTimelapse.BackgroundCorrection{channel}
%                         is not empty it is taken to be a flat field
%                         correction and is appplied to the image by
%                         element wise multiplication.
%                         if BackgroundOffset is also populated then this
%                         is first subtracted, then added back after the
%                         flat field mutliplication. This was found to
%                         prevent the flat field correction increasing the
%                         noisy camera background in the case of low image
%                         values.
%
%
% rotation - if cTimelapse.image_rotation is not zero the image is rotated
%            by this amount (in degrees) using imrotate function. Any extra
%            elements added are padded with the mean value of the image.
%            This will change the size of the image if the angle is not a
%            multiple of 90.
%
% offset - if the channel row of the array cTimelapse.offset is not zero
%          then the image is shifted by this amount (offset is specified as
%          [x_offset y_offset]). Allows images from different channels to
%          be registered properly. Any extra values are padded by the
%          mean value of the image.
%
% These corrections are applied in this order.
%
% If there is no filename matching the channel at the timepoint requested
% an image of the appropriate size of all zeros is returned and a warning
% displayed.
%
% If the channel has been loaded into memory
% (BABYTIMELAPSE.LOADCHANNELINTOMEMORY) then no image is loaded and this
% image, which will already have been corrected, is used instead.
%
% See also BABYTIMELAPSE.RETURNSINGLETIMEPOINTRAW, BABYTIMELAPSE.LOADCHANNELINTOMEMORY

if nargin<3 || isempty(channel)
    channel=1;
end

if nargin<4 || isempty(type)
    type='max';
end

% if the correct channel has been preloaded, use this image.
if isfield(cTimelapse.temporaryImageStorage,'channel') && channel == cTimelapse.temporaryImageStorage.channel
    timepointIm = cTimelapse.temporaryImageStorage.images(:,:,timepoint);
    timepointIm = double(timepointIm);
    return
end


timepointIm = returnSingleTimepointRaw(cTimelapse,timepoint,channel);

%necessary for background correction and just makes life easier
timepointIm = double(timepointIm);

%change if want things other than maximum projection
switch type
    case 'min'
        timepointIm=min(timepointIm,[],3);
    case 'max'
        timepointIm=max(timepointIm,[],3);
    case 'stack'
    case 'sum'
        timepointIm=sum(timepointIm,3);
end

timepointIm = cTimelapse.applyStandardImageTransformation(timepointIm,channel);

end
