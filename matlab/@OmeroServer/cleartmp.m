function cleartmp(obj)
% cleartmp Clear all temporary files/folders in CachePath
sp = Nursery.active;
foldersToClear = [{'tmp'},...
    sp.localOmeroCache.dirname(~sp.localOmeroCache.cExperiment)];
% Convert to full paths, and add the DownloadPath (which may be identical
% to the 'tmp' path):
foldersToClear = unique([fullfile(obj.cachePath,foldersToClear),...
    {obj.downloadPath}]);

%TODO: check whether or not to clear any downloaded images

for f=1:length(foldersToClear)
    folder = foldersToClear{f};
    if isdir(folder)
        folderContents = dir(folder);
        % Remove '..' and '.' from directory listing
        folderContents = folderContents(~ismember({folderContents.name},{'.','..'}));
        
        % Delete all files
        ffiles = {folderContents(~[folderContents.isdir]).name};
        ffiles = fullfile(folder,ffiles);
        if ~isempty(ffiles), delete(ffiles{:}); end
        
        % Delete all subdirectories
        fdirs = {folderContents([folderContents.isdir]).name};
        for d=1:length(fdirs)
            rmdir(fullfile(folder,fdirs{d}),'s');
        end        
    end
end
end
