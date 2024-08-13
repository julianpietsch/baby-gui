function cExperiment = loadFrom(cExperiment_filepath,update_savedir,copy_disco)
%BABYEXPERIMENT.LOADFROM Load a cExperiment
%
%   cExperiment = BABYEXPERIMENT.LOADFROM() opens a dialog box to
%   specify a cExperiment.mat file to load.
%
%   cExperiment = BABYEXPERIMENT.LOADFROM(FILEPATH) loads the
%   cExperiment.mat file specified in FILEPATH. If FILEPATH is empty, then
%   a dialog box to pick the file will still be opened.
%
%   By default, if the "save folder" saved with the cExperiment is
%   different to that of the specified FILEPATH, a second dialog will open
%   asking if the user wants to update the save folder. The dialog box can
%   be avoided by specifying BABYEXPERIMENT.LOADFROM(FILEPATH,'Yes') to
%   force the save folder on the cExperiment to be updated, or 
%   BABYEXPERIMENT.LOADFROM(FILEPATH,'No') to leave the save folder as
%   the version saved in the file.

if nargin<1 || isempty(cExperiment_filepath)
    fprintf('Select the saved cExperiment you want to load...\n');
    [FileName,PathName] = uigetfile('*.mat',...
        'Name of previously created cExperiment variable');
    if isequal(FileName,0) || isequal(PathName,0), return; end
    cExperiment_filepath = fullfile(PathName,FileName);
end

assert(exist(cExperiment_filepath,'file')==2,...
    'The specified file could not be found');

l1 = load(cExperiment_filepath);

was_disco = isprop(l1.cExperiment,'recovered_data') && ...
        exist('disco_baby_map','file') == 2;
if was_disco
    % Partially backwards compatible map for old DISCO-style cExperiments. 
    % This can only be called if the 'baby-gui-migration' repo is on the 
    % path and only saves trap template data from cellVision, so there will
    % be data loss for DISCO experiments.
    new_class_name = disco_baby_map(class(l1.cExperiment));
    recovered_data = l1.cExperiment.recovered_data;
    if isfield(l1,'cCellVision') && isprop(l1.cCellVision,'recovered_data')
        cTrap = l1.cCellVision.recovered_data.cTrap;
        trapTemplates = cTrap.trap1;
        if ~isempty(cTrap.trap2) && ~all(cTrap.trap1==cTrap.trap2,'all')
            trapTemplates = cat(3,trapTemplates,cTrap.trap2);
        end
        recovered_data.trapTemplates = ...
            struct('positiveExamples',trapTemplates);
    end
    assert(startsWith(new_class_name,'babyExperiment') ...
        && exist(new_class_name,'class')==8,'bad load_class_name!');
    l1.cExperiment = eval([new_class_name,'(true)']);
    l1.cExperiment.copyprops(recovered_data);
end

assert(isfield(l1,'cExperiment') && isa(l1.cExperiment,'babyExperiment'),...
    'That file does not contain a cExperiment');
cExperiment=l1.cExperiment;

cExperiment_path = fileparts(cExperiment_filepath);

if nargin<2 || isempty(update_savedir)
    update_savedir = 'No';
    if ~any(strcmp({cExperiment_path,fullfile(pwd,cExperiment_path)},...
            cExperiment.saveFolder))
        update_savedir = questdlg('The save folder from which this file was loaded does not match the save location of the cExperiment. Would you like to make them match? (If you have no idea what this means, press ''yes'')','change saveFolder');
    end
end
switch update_savedir
    case 'Yes'
        cExperiment.saveFolder = cExperiment_path;
        if ~was_disco
            cExperiment.saveExperiment;
        end
    case 'Cancel'
        error('User cancelled loading cExperiment');
end

if was_disco
    if nargin<2, copy_disco = ''; end
    if isempty(copy_disco)
        copy_disco = questdlg(...
            ['You are loading an old DISCO experiment in the new BABY format.', ...
            'Deprecated data will be lost. Do you want to save in a new ', ...
            'location?'],'Deprecation warning');
    end
    switch copy_disco
        case 'Yes'
            cExperiment.copyExperiment;
        case 'No'
            if strcmp(update_savedir,'Yes')
                cExperiment.saveExperiment;
            end
        case 'Cancel'
            error('User cancelled loading cExperiment');
        otherwise
            error('Unrecognised value for "copy_disco" arg');
    end
end

end
