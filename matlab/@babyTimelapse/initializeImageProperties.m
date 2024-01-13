function initializeImageProperties(cTimelapse,image,varargin)
% INITIALIZEIMAGEPROPERTIES subfunction of LOADTIMELAPSE that sets up
% rotation, magnification, etc.
% 
% Put in a separate function so it could be shared by children of 
% BABYTIMELAPSE. 
%
% INPUTS
% cTimelapse         -  object of the babyTimelapse class
% image              -  an image from the timelapse (for display) - of
%                       class uint8/16 or double.
% 
% 
%
% populates the rawImSize, imSize, pixelSize, image_rotation, image_flipud,
% trapsPresent, trapTemplates and trapDetectionChannel properties
%
% NOTE: This code needs to be able to handle a raw (i.e. uint8/16) image
%
% See also, BABYTIMELAPSE.LOADTIMELAPSE

ip = inputParser;
ip.addParameter('PixelSize',[],@(x) isempty(x) || ...
    (isscalar(x) && isnumeric(x)));
ip.addParameter('ImageRotation',[],@(x) isempty(x) || ...
    (isscalar(x) && isnumeric(x)));
ip.addParameter('FlipImage',false,@(x) isscalar(x) && islogical(x));
ip.addParameter('TrapTemplate',[],@(x) isempty(x) || isnumeric(x) || ...
    (isstruct(x) && isfield(x,'positiveExamples') && isnumeric(x.positiveExamples)));
ip.addParameter('TrapTemplateChannel',[],@(x) isempty(x) || ...
    (isscalar(x) && isnumeric(x)) || (ischar(x) && isrow(x)));
ip.addParameter('TrapsPresent',[],@(x) isempty(x) || ...
    (isscalar(x) && islogical(x)));
ip.parse(varargin{:});

pixel_size = ip.Results.PixelSize;
image_rotation = ip.Results.ImageRotation;
cTimelapse.image_flipud = ip.Results.FlipImage;
trapsPresent = ip.Results.TrapsPresent;
trapTemplate = ip.Results.TrapTemplate;
trapTemplateChannel = ip.Results.TrapTemplateChannel;

if isempty(trapsPresent) && ~isempty(trapTemplate)
    trapsPresent = true;
end

image = double(image(:,:,1));
cTimelapse.rawImSize=size(image);
cTimelapse.imSize = cTimelapse.rawImSize;
shown_figure = false;

if isempty(trapsPresent)
    if ~shown_figure
        h = figure;
        imshow(image,[]);
        shown_figure = true;
    end
    prompt = {'Are traps present in this Timelapse?'};
    dlg_title = 'Traps Present';
    answer = questdlg(prompt,dlg_title,'Yes','No','Yes');
    if strcmp(answer,'Yes')
        trapsPresent=true;
    else
        trapsPresent=false;
    end
end
cTimelapse.trapsPresent = trapsPresent;

if ~isempty(trapTemplate)
    trapTemplates = struct();
    if isstruct(trapTemplate)
        trapTemplates.positiveExamples = double(trapTemplate.positiveExamples);
        if isfield(trapTemplate,'negativeExamples')
            trapTemplates.negativeExamples = double(trapTemplate.negativeExamples);
        end
    else
        trapTemplates.positiveExamples = double(trapTemplate);
    end
    cTimelapse.trapTemplates = trapTemplates;
end

if ~isempty(trapTemplateChannel)
    if ischar(trapTemplateChannel)
        trapTemplateChannel = find(strcmp(cTimelapse.channelNames,...
            trapTemplateChannel),1);
    end
    cTimelapse.trapTemplateChannel = trapTemplateChannel;
end

% only asks for rotation if traps are present, otherwise it is needless.
if isempty(image_rotation)
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
        image_rotation=str2double(answer{1});
    else
        image_rotation=0;
    end
end
cTimelapse.image_rotation = image_rotation;

if isempty(pixel_size)
    if ~shown_figure
        h = figure;
        imshow(image,[]);
        shown_figure = true;
    end
    prompt = {'Enter the size of the pixels in this image in micrometers (for swainlab microscopes at 60x magnification this is 0.182 micrometers)'};
    dlg_title = 'Pixel Size';
    num_lines = 1;
    def = {'0.182'};
    answer = inputdlg(prompt,dlg_title,num_lines,def,struct('Interpreter','tex'));
    pixel_size=str2double(answer{1});
end
cTimelapse.pixelSize = pixel_size;

if shown_figure
    close(h);
end

end
