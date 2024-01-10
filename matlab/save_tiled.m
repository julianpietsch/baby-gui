function save_tiled(imgstack,filename,desc,layout)
if nargin<3 || isempty(desc)
    desc = struct();
end

if nargin<4 || isempty(layout)
    layout = 'row';
end

imgsize = arrayfun(@(x) size(imgstack,x),1:3);

desc.tilesize = imgsize(1:2);
desc.ntiles = imgsize(3);

switch layout
    case 'row'
        desc.layout = [1,imgsize(3)];
        tiled_img = num2cell(imgstack,[1,2]);
        tiled_img = horzcat(tiled_img{:});
    
    case 'column'
        desc.layout = [imgsize(3),1];
        tiled_img = num2cell(imgstack,[1,2]);
        tiled_img = vertcat(tiled_img{:});

    otherwise
        error('unrecognised tile layout specification');
end

imwrite(tiled_img,filename,'png','Description',savejson('',desc,'Compact',1));

end