function initializeImageProperties(cTimelapse,image,image_rotation,trapsPresent,pixel_size,image_flipud)
% INITIALIZEIMAGEPROPERTIES subfunction of LOADTIMELAPSE that sets up
% rotation,magnification, etc. 
% Put in a separate function so it could be shared by BABYTIMELAPSE and
% BABYTIMELAPSEOMERO. 
%
% INPUTS
% cTimelapse         -  object of the babyTimelapse class
% image              -  an image from the timelapse (for display) - of
%                       class uint8/16 or double.
% 
% all other inputs taken verbatim from BABYTIMELAPSE.LOADTIMELAPSE.
%
% populates the rawImSize and imSize properties
%
% NOTE: This code needs to be able to handle a raw (i.e. uint8/16) image
%
% See also, BABYTIMELAPSE.LOADTIMELAPSE

image = double(image(:,:,1));
cTimelapse.rawImSize=size(image);
cTimelapse.imSize = cTimelapse.rawImSize;
shown_figure = false;
%
if nargin<4 || isempty(trapsPresent)
    if ~shown_figure
        h = figure;
        imshow(image,[]);
        shown_figure = true;
    end
    prompt = {'Are traps present in this Timelapse?'};
    dlg_title = 'Traps Present';
    answer = questdlg(prompt,dlg_title,'Yes','No','Yes');
    if strcmp(answer,'Yes')
        cTimelapse.trapsPresent=true;
    else
        cTimelapse.trapsPresent=false;
    end
else
    cTimelapse.trapsPresent=trapsPresent;
end

% only asks for rotation of traps are present, otherwise it is needless.
if nargin<3 || isempty(image_rotation)
    if cTimelapse.trapsPresent
        if ~shown_figure
            h = figure;
            imshow(image,[]);
            shown_figure = true;
        end
        prompt = {'Enter the rotation (in degrees counter-clockwise) required to orient opening of traps to the left'};
        dlg_title = 'Rotation';
        num_lines = 1;
        def = {'0'};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        cTimelapse.image_rotation=str2double(answer{1});
    else
        cTimelapse.image_rotation=0;
    end
else
    cTimelapse.image_rotation=image_rotation;
end
if nargin<5 || isempty(pixel_size)
    if ~shown_figure
        h = figure;
        imshow(image,[]);
        shown_figure = true;
    end
    prompt = {'Enter the size of the pixels in this image in micrometers (for swainlab microscopes at 60x magnification this is 0.263 micrometers)'};
    dlg_title = 'Pixel Size';
    num_lines = 1;
    def = {'0.263'};
    answer = inputdlg(prompt,dlg_title,num_lines,def,struct('Interpreter','tex'));
    cTimelapse.pixelSize=str2double(answer{1});
else
    cTimelapse.pixelSize=pixel_size;
end
if nargin<6 || isempty(image_flipud)
    cTimelapse.image_flipud = false;
else
    cTimelapse.image_flipud = image_flipud;
end

if shown_figure
    close(h);
end

end
