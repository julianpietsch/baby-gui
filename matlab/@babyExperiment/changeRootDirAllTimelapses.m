function changeRootDirAllTimelapses(cExperiment,dirsToSearch,newRootFolder)
% changeRootDirAllTimelapses(cExperiment,dirsToSearch,newRootFolder)
%
% changes the timelapseDir property of every cTimelapse covered by the
% cExperiment.

if nargin<2 || isempty(dirsToSearch)
    dirsToSearch=1:length(cExperiment.dirs);
end

if nargin<3 || isempty(newRootFolder)

%fix for windows file systems:
oldRootFolder = cExperiment.rootFolder;

if ischar(oldRootFolder)
    oldRootFolder = regexprep(oldRootFolder,'\\','/');
end

fprintf('Select the correct folder for: \n %s \n',oldRootFolder);
h = helpdlg('This is the new root folder containing all images of the timelapse objects');
waitfor(h);

newRootFolder=uigetdir(pwd,sprintf('Select the correct folder for: %s',oldRootFolder));
if newRootFolder==0
    
    fprintf('\n\nchangeRootDirAll cancelled\n\n')
    return
    
end
end
% cExperiment.saveFolder=uigetdir(pwd,['Select the folder where you want to save the timelapses: ',cExperiment.rootFolder]);

cExperiment.rootFolder=newRootFolder;
fprintf('Saving changes to timelapses');
for i=1:length(dirsToSearch)
    fprintf('.');
    posIndex=dirsToSearch(i);
    cTimelapse = cExperiment.loadCurrentTimelapse(posIndex);
    newDir=fullfile(newRootFolder,cExperiment.dirs{posIndex});
    cTimelapse.timelapseDir=newDir;
    cExperiment.saveTimelapseExperiment([],false)
end
cExperiment.rootFolder=newRootFolder;
cExperiment.saveExperiment;
fprintf('\n\n')
