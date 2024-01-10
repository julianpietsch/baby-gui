function saveTimelapseExperiment(cExperiment,currentPos, saveCE)
% saveTimelapseExperiment(cExperiment,currentPos, saveCE)
% 
% saves cExperiment.cTimelapse to:
%   [cExperiment.saveFolder filesep cExperiment.dirs{currentPos},'cTimelapse']
%
% also saves the cExperiment to:
%       [cExperiment.saveFolder filesep 'cExperiment.mat']
%
% If currentPos is not provided, cExperiment.currentPos (populated when
% BABYEXPERIMENT.LOADCURRENTTIMELAPSE is called) is used. It will be
% empty if babyExperiment.TimelapseTraps has been replaced by a
% non-identical object (see BABYEXPERIMENT.SET.CTIMELAPSE)
% 
% Third input is boolean - saveCE: logical - if true,
% save the cExperiment file as well as the timelapse. Defaults to false.
%
% See also, BABYEXPERIMENT.LOADCURRENTTIMELAPSE
if nargin<2 || isempty(currentPos)
    currentPos = cExperiment.currentPos;
end

if nargin<3
    saveCE=false;
end

cTimelapse = cExperiment.cTimelapse;
cTimelapseFilename=[cExperiment.saveFolder filesep cExperiment.dirs{currentPos},'cTimelapse'];
save(cTimelapseFilename,'cTimelapse');

if saveCE
    cExperiment.saveExperiment; 
end

end
