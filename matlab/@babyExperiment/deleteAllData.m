function deleteAllData(cExperiment,positionsToDelete)
% deleteAllData(cExperiment,positionsToDelete)
%
% Delete the extracted data (i.e. cellInf and extractedData) from
% cExperiment and each of the cTimelapse objects. Usually done to make
% loading faster and save space when reprocessing.


if nargin<2
    positionsToDelete=1:length(cExperiment.dirs);
end

for i=1:length(positionsToDelete)
    currentPos=positionsToDelete(i);
    cTimelapse=cExperiment.loadCurrentTimelapse(currentPos);
    cTimelapse.extractedData = [];
    cExperiment.cTimelapse = cTimelapse;
    cExperiment.saveTimelapseExperiment;   
    clear cTimelapse;
end

cExperiment.cellInf = [];

cExperiment.saveExperiment;

end
