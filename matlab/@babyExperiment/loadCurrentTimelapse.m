function cTimelapse=loadCurrentTimelapse(cExperiment,positionToLoad)
%cTimelapse=loadCurrentTimelapse(cExperiment,positionToLoad)
%
% loads a cTimelapse. PositionToLoad should be a single number
% indicating which position to load. Note, positionToLoad indicated index
% in cExperiment.dirs to load, so depending on the ordering of the
% directories in dirs cExperiment.loadCurrentTimelapse(2) will not
% necessarily load the cTimlapse related to directory pos2, and will in
% general load pos10 - his is due to alphabetic ordering.
%
% now just a wrapper for returnTimelapse plus population cTimelapse field
% of cExperiment and currentPos
%
% good if you want to cycle through timelapse loading and saving.
%
% See also, BABYEXPERIMENT.RETURNTIMELAPSE
cTimelapse=cExperiment.returnTimelapse(positionToLoad);

cExperiment.cTimelapse=cTimelapse;

cExperiment.currentPos = positionToLoad;
end    
