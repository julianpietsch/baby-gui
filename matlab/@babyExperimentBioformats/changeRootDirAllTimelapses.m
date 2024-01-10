function changeRootDirAllTimelapses(cExperiment,poses,newSharedFolder)
% changeRootDirAllTimelapses(cExperiment,poses,newSharedFolder)
%
% Updates the timelapseDir property of every cTimelapse covered by the
% cExperiment. NB: does not update if all files in cExperiment.rootFolder
% are found.

if nargin<2 || isempty(poses)
    poses=1:numel(cExperiment.dirs);
end

if nargin<3 || isempty(newSharedFolder)
    newSharedFolder = [];
end

rootFolders = cellstr(cExperiment.rootFolder);
if numel(rootFolders) == numel(cExperiment.dirs)
    rootFolders = rootFolders(poses);
end
if ~all(cellfun(@isfolder,rootFolders)) ...
        && any(cellfun(@(x) exist(x,'file')~=2,rootFolders))
    % Attempt to find a shared folder
    splitfolders = cellfun(@(f) strsplit(f,{'\\','/'}),rootFolders,'Uniform',false);
    maxsplits = min(cellfun(@numel,splitfolders));
    sharedsplits = cellfun(@(f) f(1:maxsplits),splitfolders,'Uniform',false);
    sharedsplits = vertcat(sharedsplits{:});
    all_matching = arrayfun(@(f) ...
        all(strcmp(sharedsplits(:,f),sharedsplits(1,f))),1:maxsplits);
    i_shared = find(all_matching,1,'last');
    if isempty(i_shared)
        error('Cannot find a shared root directory between files!');
    end
    fprintf('Shared root folder is currently set as:\n  %s\n',...
        fullfile(splitfolders{1}{1:i_shared}));
    if numel(splitfolders)==1
        i_shared = i_shared - 1;
    end
    relpaths = cellfun(@(f) fullfile(f{i_shared+1:end}),splitfolders,'Uniform',false);
    if isempty(newSharedFolder)
        fprintf('Please pick a new shared root folder...\n');
        newSharedFolder = uigetdir('',...
            'Please pick a new shared root folder...');
        if isequal(newSharedFolder,0), return; end % cancelled
    end
    newPaths = cellfun(@(f) fullfile(newSharedFolder,f),relpaths,'Uniform',false);
    assert(all(cellfun(@(f) exist(f,'file')==2,newPaths)),...
        'At least some of the files could not be found in that root folder.');
    fprintf('Saving changes to timelapses');
    cExperiment.rootFolder = newPaths;
    for p=poses
        fprintf('.');
        cTimelapse = cExperiment.loadCurrentTimelapse(p);
        fi = find(strcmp(rootFolders,cTimelapse.timelapseDir),1);
        assert(~isempty(fi),...
            'timelapse files differ from those in cExperiment!');
        cTimelapse.timelapseDir = newPaths{fi};
        cExperiment.saveTimelapseExperiment([],false);
    end
    cExperiment.saveExperiment;
    fprintf('\nShared root folder has been updated to:\n  %s\n',newSharedFolder);
else
    fprintf('All files were found.\n');
end
