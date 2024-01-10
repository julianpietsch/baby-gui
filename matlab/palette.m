function [ colour ] = palette(name)
%palette Return a predefined colour by name
%   name: either a string defining a colour name or an integer index
%   to specify default palette. Specify 'all' to return a cell array of
%   all possible colour names.

%% Define colours in a map
% Make the colours variable persistent so that it only needs to be defined
% on the first call to this function
persistent colours;

if isempty(colours)
    colours = containers.Map('KeyType','char','ValueType','any');

    %% Pastels
    colours('snow') = addcolour([255,250,250]);
    colours('snow2') = addcolour([238,233,233]);
    colours('snow3') = addcolour([205,201,201]);
    colours('snow4') = addcolour([139,137,137]);
    colours('ghostwhite') = addcolour([248,248,255]);
    colours('whitesmoke') = addcolour([245,245,245]);
    colours('gainsboro') = addcolour([220,220,220]);
    colours('floralwhite') = addcolour([255,250,240]);
    colours('oldlace') = addcolour([253,245,230]);
    colours('linen') = addcolour([240,240,230]);
    colours('antiquewhite') = addcolour([250,235,215]);
    colours('antiquewhite2') = addcolour([238,223,204]);
    colours('antiquewhite3') = addcolour([205,192,176]);
    colours('antiquewhite4') = addcolour([139,131,120]);
    colours('papayawhip') = addcolour([255,239,213]);
    colours('blanchedalmond') = addcolour([255,235,205]);
    colours('bisque') = addcolour([255,228,196]);
    colours('bisque2') = addcolour([238,213,183]);
    colours('bisque3') = addcolour([205,183,158]);
    colours('bisque4') = addcolour([139,125,107]);
    colours('peachpuff') = addcolour([255,218,185]);
    colours('peachpuff2') = addcolour([238,203,173]);
    colours('peachpuff3') = addcolour([205,175,149]);
    colours('peachpuff4') = addcolour([139,119,101]);
    colours('navajowhite') = addcolour([255,222,173]);
    colours('moccasin') = addcolour([255,228,181]);
    colours('cornsilk') = addcolour([255,248,220]);
    colours('cornsilk2') = addcolour([238,232,205]);
    colours('cornsilk3') = addcolour([205,200,177]);
    colours('cornsilk4') = addcolour([139,136,120]);
    colours('ivory') = addcolour([255,255,240]);
    colours('ivory2') = addcolour([238,238,224]);
    colours('ivory3') = addcolour([205,205,193]);
    colours('ivory4') = addcolour([139,139,131]);
    colours('lemonchiffon') = addcolour([255,250,205]);
    colours('seashell') = addcolour([255,245,238]);
    colours('seashell2') = addcolour([238,229,222]);
    colours('seashell3') = addcolour([205,197,191]);
    colours('seashell4') = addcolour([139,134,130]);
    colours('honeydew') = addcolour([240,255,240]);
    colours('honeydew2') = addcolour([244,238,224]);
    colours('honeydew3') = addcolour([193,205,193]);
    colours('honeydew4') = addcolour([131,139,131]);
    colours('mintcream') = addcolour([245,255,250]);
    colours('azure') = addcolour([240,255,255]);
    colours('aliceblue') = addcolour([240,248,255]);
    colours('lavender') = addcolour([230,230,250]);
    colours('lavenderblush') = addcolour([255,240,245]);
    colours('mistyrose') = addcolour([255,228,225]);
    colours('white') = addcolour([255,255,255]);

    %% Greys
    colours('black') = addcolour([0,0,0]);
    colours('darkslategrey') = addcolour([49,79,79]);
    colours('dimgrey') = addcolour([105,105,105]);
    colours('slategrey') = addcolour([112,138,144]);
    colours('lightslategrey') = addcolour([119,136,153]);
    colours('darkgrey') = addcolour([130,130,130]);
    colours('grey') = addcolour([190,190,190]);
    colours('lightgrey') = addcolour([211,211,211]);

    %% Blues
    colours('midnightblue') = addcolour([25,25,112]);
    colours('navy') = addcolour([0,0,128]);
    colours('cornflowerblue') = addcolour([100,149,237]);
    colours('darkslateblue') = addcolour([72,61,139]);
    colours('slateblue') = addcolour([106,90,205]);
    colours('mediumslateblue') = addcolour([123,104,238]);
    colours('lightslateblue') = addcolour([132,112,255]);
    colours('mediumblue') = addcolour([0,0,205]);
    colours('royalblue') = addcolour([65,105,225]);
    colours('blue') = addcolour([0,0,255]);
    colours('dodgerblue') = addcolour([30,144,255]);
    colours('deepskyblue') = addcolour([0,191,255]);
    colours('skyblue') = addcolour([135,206,250]);
    colours('lightskyblue') = addcolour([135,206,250]);
    colours('steelblue') = addcolour([70,130,180]);
    colours('lightsteelblue') = addcolour([176,196,222]);
    colours('lightblue') = addcolour([173,216,230]);
    colours('powderblue') = addcolour([176,224,230]);
    colours('paleturquoise') = addcolour([175,238,238]);
    colours('darkturquoise') = addcolour([0,206,209]);
    colours('mediumturquoise') = addcolour([72,209,204]);
    colours('turquoise') = addcolour([64,224,208]);
    colours('cyan') = addcolour([0,255,255]);
    colours('lightcyan') = addcolour([224,255,255]);
    colours('cadetblue') = addcolour([95,158,160]);

    %% Greens
    colours('mediumaquamarine') = addcolour([102,205,170]);
    colours('aquamarine') = addcolour([127,255,212]);
    colours('darkgreen') = addcolour([0,100,0]);
    colours('darkolivegreen') = addcolour([85,107,47]);
    colours('darkseagreen') = addcolour([143,188,143]);
    colours('seagreen') = addcolour([46,139,87]);
    colours('mediumseagreen') = addcolour([60,179,113]);
    colours('lightseagreen') = addcolour([32,178,170]);
    colours('palegreen') = addcolour([152,251,152]);
    colours('springgreen') = addcolour([0,255,127]);
    colours('lawngreen') = addcolour([124,252,0]);
    colours('chartreuse') = addcolour([127,255,0]);
    colours('mediumspringgreen') = addcolour([0,250,154]);
    colours('greenyellow') = addcolour([173,255,47]);
    colours('limegreen') = addcolour([50,205,50]);
    colours('yellowgreen') = addcolour([154,205,50]);
    colours('forestgreen') = addcolour([34,139,34]);
    colours('olivedrab') = addcolour([107,142,35]);
    colours('darkkhaki') = addcolour([189,183,107]);
    colours('khaki') = addcolour([240,230,140]);
    colours('green') = addcolour([0,255,0]);

    %% Yellows
    colours('palegoldenrod') = addcolour([238,232,170]);
    colours('lightgoldenrodyellow') = addcolour([250,250,210]);
    colours('lightyellow') = addcolour([255,255,224]);
    colours('yellow') = addcolour([255,255,0]);
    colours('gold') = addcolour([255,215,0]);
    colours('lightgoldenrod') = addcolour([238,221,130]);
    colours('goldenrod') = addcolour([218,165,32]);
    colours('darkgoldenrod') = addcolour([184,134,11]);

    %% Browns
    colours('rosybrown') = addcolour([188,143,143]);
    colours('indianred') = addcolour([205,92,92]);
    colours('saddlebrown') = addcolour([139,69,19]);
    colours('sienna') = addcolour([160,82,45]);
    colours('peru') = addcolour([205,133,63]);
    colours('burlywood') = addcolour([222,184,135]);
    colours('beige') = addcolour([245,245,220]);
    colours('wheat') = addcolour([245,222,179]);
    colours('sandybrown') = addcolour([244,164,96]);
    colours('tan') = addcolour([210,180,140]);
    colours('chocolate') = addcolour([210,105,30]);
    colours('firebrick') = addcolour([178,34,34]);
    colours('brown') = addcolour([165,42,42]);

    %% Oranges
    colours('darksalmon') = addcolour([233,150,122]);
    colours('salmon') = addcolour([250,128,114]);
    colours('lightsalmon') = addcolour([255,160,122]);
    colours('orange') = addcolour([255,165,0]);
    colours('darkorange') = addcolour([255,140,0]);
    colours('coral') = addcolour([255,127,80]);
    colours('lightcoral') = addcolour([240,128,128]);
    colours('tomato') = addcolour([255,99,71]);
    colours('orangered') = addcolour([255,69,0]);
    colours('red') = addcolour([255,0,0]);

    %% Pinks/violets
    colours('hotpink') = addcolour([255,105,180]);
    colours('deeppink') = addcolour([255,20,147]);
    colours('pink') = addcolour([255,192,203]);
    colours('lightpink') = addcolour([255,182,193]);
    colours('palevioletred') = addcolour([219,112,147]);
    colours('maroon') = addcolour([176,48,96]);
    colours('mediumvioletred') = addcolour([199,21,133]);
    colours('violetred') = addcolour([208,32,144]);
    colours('violet') = addcolour([238,130,238]);
    colours('plum') = addcolour([221,160,221]);
    colours('orchid') = addcolour([218,112,214]);
    colours('mediumorchid') = addcolour([186,85,211]);
    colours('darkorchid') = addcolour([153,50,204]);
    colours('darkviolet') = addcolour([148,0,211]);
    colours('blueviolet') = addcolour([138,43,226]);
    colours('purple') = addcolour([160,32,240]);
    colours('mediumpurple') = addcolour([147,112,219]);
    colours('thistle') = addcolour([216,191,216]);
    
    %% MI paper colours
    colours('happy') = addcolour([112,143,195]);
    colours('oxidative') = addcolour([230,103,125]);
    colours('osmotic') = addcolour([248,190,107]);
    colours('carbon') = addcolour([151,202,117]);
    
    %% Wong (2011) Nat. Methods (colour blind safe palette)
    colours('wong1') = addcolour([0,0,0]); % black
    colours('wong2') = addcolour([230,159,0]); % orange
    colours('wong3') = addcolour([86,180,233]); % sky blue
    colours('wong4') = addcolour([0,158,115]); % bluish green
    colours('wong5') = addcolour([240,228,66]); % yellow
    colours('wong6') = addcolour([0,114,178]); % blue
    colours('wong7') = addcolour([213,94,0]); % vermillion
    colours('wong8') = addcolour([204,121,167]); % reddish purple

    %% IBM design library (colour blind safe palette)
    colours('ibm1') = addcolour([100,143,255]);
    colours('ibm2') = addcolour([120,94,240]);
    colours('ibm3') = addcolour([220,38,127]);
    colours('ibm4') = addcolour([254,97,0]);
    colours('ibm5') = addcolour([255,176,0]);
    
    %% Paul Tol (https://personal.sron.nl/~pault/)
    % Vibrant palette
    colours('vibrant_blue') = addcolour([0,119,187]);
    colours('vibrant_cyan') = addcolour([51,187,238]);
    colours('vibrant_teal') = addcolour([0,153,136]);
    colours('vibrant_orange') = addcolour([238,119,51]);
    colours('vibrant_red') = addcolour([204,51,17]);
    colours('vibrant_magenta') = addcolour([238,51,119]);
    colours('vibrant_grey') = addcolour([187,187,187]);

    % Muted palette
    colours('muted_indigo') = addcolour([51,34,136]);
    colours('muted_cyan') = addcolour([136,204,238]);
    colours('muted_teal') = addcolour([68,170,153]);
    colours('muted_green') = addcolour([17,119,51]);
    colours('muted_olive') = addcolour([153,153,51]);
    colours('muted_sand') = addcolour([221,204,119]);
    colours('muted_rose') = addcolour([204,102,119]);
    colours('muted_wine') = addcolour([136,34,85]);
    colours('muted_purple') = addcolour([170,68,153]);
    colours('muted_palegrey') = addcolour([221,221,221]);

    % Medium contrast palette
    colours('medium_light_yellow') = addcolour([238,204,102]);
    colours('medium_light_red') = addcolour([238,153,170]);
    colours('medium_light_blue') = addcolour([102,153,204]);
    colours('medium_dark_yellow') = addcolour([153,119,0]);
    colours('medium_dark_red') = addcolour([153,68,85]);
    colours('medium_dark_blue') = addcolour([0,68,136]);
    
    % Bright colour palette (default)
    colours('bright_blue') = addcolour([68,119,170]);
    colours('bright_cyan') = addcolour([102,204,238]);
    colours('bright_green') = addcolour([34,136,51]);
    colours('bright_yellow') = addcolour([204,187,68]);
    colours('bright_red') = addcolour([238,102,119]);
    colours('bright_purple') = addcolour([170,51,119]);
    colours('bright_grey') = addcolour([187,187,187]);
    
end

if all(strcmp(name,'all'))
    colour = keys(colours);
    % order the colours by the ordering index:
    [~,colOrder] = sort(cellfun(@(x) x(4), values(colours)));
    colour = colour(colOrder);
    return
end

%% Set a default palette
default_palette = {'bright_blue','bright_red','bright_yellow',...
    'bright_cyan','bright_purple','bright_green','bright_grey'};

% For numeric arguments, pick a colour from the default palette, looping if
% the number exceeds the length of the palette:
if isnumeric(name)
    if size(name,2)==3
        % Assume that a correctly-formatted RGB triplet array has been provided
        colour = name;
        return
    elseif length(name)==1
        % Pick colour from default palette
        name = default_palette{mod(floor(name)-1,length(default_palette))+1};
    else
        error('Can only specify palette indices one at a time...');
    end
end

% If a cell array has been provided, recursively call palette for each cell
if iscell(name)
    colour = zeros(length(name),3);
    for i=1:length(name)
        colour(i,:) = palette(name{i});
    end
    return
end

if ~isKey(colours,name)
    error(['The specified colour "',name,'" is not in the palette']);
end

colour = colours(name);
colour = colour(1:3)./255; % Ignore the last ordering index

end

function cvec = addcolour(rgb)
%addcolour Add an index to store the order that the colours were added
persistent i
if isempty(i)
    i = 1;
end
cvec = [rgb,i];
i = i+1;
end
