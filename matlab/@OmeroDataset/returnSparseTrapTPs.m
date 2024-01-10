function [trapStack] = returnSparseTrapTPs(this,cTimelapse,traps,tps,channel)
%RETURNSPARSETRAPTPS Efficiently retrieve images for a trap/tp sequence
%   TRAPSTACK=DS.RETURNSPARSETRAPTPS(CTIMELAPSE,TRAPS,TPS,CHANNEL) returns 
%   a NX x NY x NZ x Ntraps array of images from a babyTimelapse object 
%   CTIMELAPSE corresponding to the trap-time point pairs specified in the 
%   (equal-length) arrays TRAPS and TPS for the channel name CHANNEL 
%   specified as a char vector.
%   
%   TODO: does not currently account for specification of channel offset

ntraps = numel(traps);
assert(ntraps==numel(tps),...
    '"traps" and "tps" must be of equal length specifying trap/tp pairs');

this.ensure_session;
this.pos = cTimelapse.metadata.posname;
imsize = this.imageSize;

% Determine the channels and sections we need
assert(ismember(channel,cTimelapse.channelNames),...
    'the specified channel "%s" is not present in the cTimelapse',channel);
ci = find(arrayfun(@(x) strncmpi(x,channel,length(x)),this.channelNames));
assert(~isempty(ci),'the specified channel "%s" was not on Omero',channel);
if ismember(channel,this.channelNames)
    % Channel may have Z-sections
    if this.hasZstacks(ci), zsect = []; nz = imsize(3);
    else, zsect = 1; nz = 1; end
else
    % Channel name specifies Z-section as trailing digits
    zsect = str2double(channel(numel(this.channelNames{ci})+2:end));
    nz = 1;
end

% A number of transformations were potentially applied to the
% image before retrieving the traps. We need to reverse these
% to identify trap location in the raw image, and then reapply to the trap
% images for correct orientation.

% We need to undo image rotation to identify trap location in the raw image
image_rotation = cTimelapse.image_rotation;
assert(mod(image_rotation,90)==0,...
    'cannot handle image rotations that are not multiples of 90 degrees');
N_rot90 = floor(image_rotation/90);
% NB: imrotate applies rotations in (y,x) space, so to undo that rotation
% in (x,y) space, we simply need to apply the same rotation directly. The
% image rotation is applied after any rescaling of the raw image, so the
% mid-point may need to be treated differently than for operations on raw
invrot = image_rotation;
scaled_rotmat = [cosd(invrot) -sind(invrot); sind(invrot) cosd(invrot)];
scaled_mid = (cTimelapse.imSize(:)+1)/2;

% Remaining transformations originally applied before image scaling
rotmat = eye(2); % identity
mid = (reshape(imsize(1:2),[],1)+1)/2;
bbsize = struct2array(cTimelapse.cTrapSize);
swapbb = mod(image_rotation,180)~=0;

% A bug in the microscope code meant that some channels were saved
% flipped relative to brightfield. The flipchannels property can be
% used to specify which channels should be flipped to match the
% orientation of the brightfield image.
flipchannels = cTimelapse.flipchannels;
cTci = find(strcmp(cTimelapse.channelNames,channel),1);
doflip = ~isempty(flipchannels) && ...
    length(flipchannels)>cTci && ...
    flipchannels(cTci);
if doflip
    % undo vertical flip
    rotmat = [-1,0;0,1]*rotmat;
end

%Images are flipped in both directions on Omero upload - so if the
%data was segmented from a folder before being uploaded/converted
%then it should be flipped to ensure data is consistent with any
%segmentation results
isfolderexpt = strcmp(cTimelapse.segmentationSource,'Folder');
if isfolderexpt
    % undo vertical flip
    rotmat = [-1,0;0,1]*rotmat;
    % undo 90 degree rotation
    rotmat = [0,1;-1,0]*rotmat;
    swapbb = ~swapbb;
end

if swapbb, bbsize = flip(bbsize); end

% Use the same pixel store to get all tiles for this position
[store,pixels] = this.database.Session.getRawPixelsStore(this.image);
try
    type = char(pixels.getPixelsType().getValue().getValue());
    if strcmp(type,'float'), type = 'single'; end
    trapSize = [bbsize*2+1,nz,ntraps];
    trapStack = zeros(trapSize,type);
    for ti=1:ntraps
        trap = traps(ti);
        tp = tps(ti);
        loc = struct2array(cTimelapse.cTimepoint(tp).trapLocations(trap));
        loc = round(scaled_rotmat*(loc(:)-scaled_mid)+scaled_mid);
        loc = round(rotmat*(double(loc(:))-mid)+mid);
        lloc = loc(:)-bbsize(:); uloc = loc(:)+bbsize(:);
        im = this.getHypercube(...
            'X',[max(1,lloc(1)),min(imsize(1),uloc(1))],...
            'Y',[max(1,lloc(2)),min(imsize(2),uloc(2))],...
            'Z',zsect,'C',ci,'T',tp,'pixelstore',{store,pixels});
        if size(im,1)~=trapSize(1) || size(im,2)~=trapSize(2)
            % If trap goes outside image border, need to pad
            imraw = im;
            im = ones(trapSize(1:3),'like',imraw)*median(imraw(:));
            im((1:size(imraw,1))+max(0,1-lloc(2)),...
                (1:size(imraw,2))+max(0,1-lloc(1)),:) = imraw;
        end
        if isfolderexpt
            % apply 90 degree rotation and then vertical flip
            im = rot90(im);
            im = flipud(im);
        end
        if doflip, im = flipud(im); end
        % Apply rotation originally specified by cTimelapse
        im = rot90(im,N_rot90);
        trapStack(:,:,:,ti) = median(im(:));
        trapStack(1:size(im,1),1:size(im,2),:,ti) = im;
    end
catch err
    store.close(); % make sure we clean up on error as well
    rethrow(err);
end
store.close();
            
end
