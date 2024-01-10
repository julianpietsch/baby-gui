function saveTimelapse(cExperiment,currentPos)
% saveTimelapse(cExperiment,currentPos)
%
% now just a wrapper for BABYEXPERIMENT.SAVETIMELAPSEEXPERIMENT  with
% save experiment set to false.
% only kept for legacy reasons.

cExperiment.saveTimelapseExperiment(currentPos,false);


end
