function cTimelapse = loadFrom(cTimelapse_filepath)
%BABYTIMELAPSE.LOADFROM Load a cExperiment
%
%   cTimelapse = BABYTIMELAPSE.LOADFROM() opens a dialog box to
%   specify a cTimelapse.mat file to load.
%
%   cTimelapse = BABYTIMELAPSE.LOADFROM(FILEPATH) loads the
%   cTimelapse.mat file specified in FILEPATH. If FILEPATH is empty, then
%   a dialog box to pick the file will still be opened.

if nargin<1 || isempty(cTimelapse_filepath)
    fprintf('Select the saved cTimelapse you want to load...\n');
    [FileName,PathName] = uigetfile('*.mat',...
        'Name of previously created cTimelapse variable');
    if isequal(FileName,0) || isequal(PathName,0), return; end
    cTimelapse_filepath = fullfile(PathName,FileName);
end

assert(exist(cTimelapse_filepath,'file')==2,...
    'The specified file could not be found');

l1 = load(cTimelapse_filepath);
cTimelapse = l1.cTimelapse;

was_disco = isprop(cTimelapse,'recovered_data') && ...
    exist('disco_baby_map','file') == 2;
if was_disco
    % Partially backwards compatible map for old DISCO-style cTimelapses. 
    % This can only be called if the 'baby-gui-migration' repo is on the 
    % path and only saves trap template data from cellVision, so there will
    % be data loss for DISCO experiments.
    new_class_name = disco_baby_map(class(cTimelapse));
    recovered_data = cTimelapse.recovered_data;
    assert(startsWith(new_class_name,'babyTimelapse') ...
        && exist(new_class_name,'class')==8,'bad load_class_name!');
    cTimelapse = eval([new_class_name,'([],true)']);
    cTimelapse.copyprops(recovered_data);
end

end
