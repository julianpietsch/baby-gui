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

if isequal(layout,'row')
    desc.layout = [1,imgsize(3)];
    tiled_img = num2cell(imgstack,[1,2]);
    tiled_img = horzcat(tiled_img{:});
elseif isequal(layout,'column')
    desc.layout = [imgsize(3),1];
    tiled_img = num2cell(imgstack,[1,2]);
    tiled_img = vertcat(tiled_img{:});
else
    desc.layout = layout;
    tiled_img = num2cell(imgstack,[1,2]);
    tiled_img = reshape(tiled_img,layout(1),[]);
    for r=1:size(tiled_img,1)
        tiled_img{r,1} = [tiled_img{r,:}];
    end
    tiled_img = vertcat(tiled_img{:,1});
end

imwrite(tiled_img,filename,'png','Description',savejson('',desc,'Compact',1));

end