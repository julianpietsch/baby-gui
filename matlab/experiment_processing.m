%% -- EXPERIMENT INITIALISATION -- %%

% Want MM traps with mother cell at bottom
cExperiment = experimentInitGUI;
cExperiment.identifyTrapsTimelapses([],false);
% Since traps are thin we need to specify magnitude of 'huge_move' param
cExperiment.trackTrapsInTime([],[],[],[],200,[150,150]);

%% BabyBrain setup for phase contrast
cExperiment.babyBrain.map_channel('phase','channel1');
cExperiment.babyBrain.config.camera = 'mmsCMOS';
cExperiment.saveExperiment;

%% BabyBrain setup for FRET segmentation
cExperiment.babyBrain.config.camera = 'mmEMCCD';
cExperiment.babyBrain.config.channel = {'CFP','YFP'};
cExperiment.saveExperiment;

%% Set extract parameters
available_extract_channels = cExperiment.channelNames(:)';
[extract_channels, ok] = listdlg('Name','Channels for extraction',...
    'ListString',available_extract_channels,'SelectionMode','multiple');
if ok
    extract_channels = available_extract_channels(extract_channels);
    extraction_parameters = struct();
    extraction_parameters.extractFunction = @extractCellDataMotherMachine;
    extraction_parameters.functionParameters = struct();
    extraction_parameters.functionParameters.channels = extract_channels;
    extraction_parameters.functionParameters.type = 'max';
    if exist('poses','var')~=1, poses = 1:numel(cExperiment.dirs); end
    cExperiment.setExtractParameters(poses,extraction_parameters);
end

%% -- EXPERIMENT PROCESSING -- %%

%% Run segmentation and cell tracking
cExperiment.babySegment([],'refine_outlines',false);
cExperiment.mmRetrackCells;
%cExperiment.mmRetrackCells([],[],[],true,true,0);

%% Track inverted (pos, trap, tp, replace first tp, invert, gr_init = 0)
cExperiment.mmRetrackCells([],[],[],true,true,0);

%% Curate and select tracks with GUI
% Autoselection within traps is available in the annotation tab
bgui = babyGUI(cExperiment);

%% Extract data
if exist('poses','var')~=1, poses = 1:numel(cExperiment.dirs); end
cExperiment.extractCellInformation(poses,false);
cExperiment.babyExtractLineage(poses,false,[],true); % recalculate ellipse params

%% Define chamber groups
cExperiment.chamberMap.all = 1:numel(cExperiment.dirs);
cExperiment.saveExperiment;

%% Compile chamber groups to h5 files
cExperiment.compileChambers('cellfilt',@(x) sum(x(1).radius>0,2)>4,...
    'vol_type','length','log_gr',true,'gr_method','slmfixed',...
    'mW',15,'dW',15,'birthHandling','cumulative','usefilters',false,...
    'save_D',false);

%% -- GENERAL TASKS -- %%

%% Load experiment
cExperiment = babyExperiment.loadFrom;

%% Inspect with BABY GUI
bgui = babyGUI(cExperiment);

%% Copy experiment
cExperiment.copyExperiment;

%% Check that correct image file is being used
rootFolders = cellstr(cExperiment.rootFolder);
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
    fprintf('Please pick a new shared root folder...\n');
    newRoot = uigetdir('','Please pick a new shared root folder...');
    if ~isequal(newRoot,0)
        newPaths = cellfun(@(f) fullfile(newRoot,f),relpaths,'Uniform',false);
        assert(all(cellfun(@(f) exist(f,'file')==2,newPaths)),...
            'At least some of the files could not be found in that root folder.');
        fprintf('Saving changes to timelapses');
        cExperiment.rootFolder = newPaths;
        for p=1:numel(cExperiment.dirs)
            fprintf('.');
            cTimelapse = cExperiment.loadCurrentTimelapse(p);
            fi = find(strcmp(rootFolders,cTimelapse.timelapseDir),1);
            assert(~isempty(fi),...
                'timelapse files differ from those in cExperiment!');
            cTimelapse.timelapseDir = newPaths{fi};
            cExperiment.saveTimelapseExperiment([],false);
        end
        cExperiment.saveExperiment;
        fprintf('\nShared root folder has been updated to:\n  %s\n',newRoot);
    end
else
    fprintf('All files were found.\n');
end
