function cropTimepoints(cExperiment,positionsToCrop)
%cropTimepoints(cExperiment,positionsToCrop)
%Completely removes a set of timepoints from each cTimelapse in the
%cExperiment. Timepoints to keep are selected by GUI.

num_lines=1;
dlg_title = 'Tp to crop?';
prompt = {['This removes the timepoints from the timelapses completely. If you end up wanting to use the timepoints that you have cropped, you will need to create a new cTimelapse file.' ...
    'Enter the first timepoint that would like to keep.'], 'And the last timepoint that you would like to keep.'};    
answer = inputdlg(prompt,dlg_title,num_lines);
startTP=str2double(answer{1});
endTP=str2double(answer{2});

if nargin<2
    positionsToCrop=1:length(cExperiment.dirs);
end

for i=1:length(positionsToCrop)
    currentPos=positionsToCrop(i);
    cExperiment.cTimelapse=cExperiment.returnTimelapse(currentPos);
    cExperiment.cTimelapse.cTimepoint=cExperiment.cTimelapse.cTimepoint(startTP:endTP);
    cExperiment.cTimelapse.timepointsProcessed = cExperiment.cTimelapse.timepointsProcessed(startTP:endTP);
    cExperiment.cTimelapse.timepointsToProcess(~ismember(cExperiment.cTimelapse.timepointsToProcess,startTP:endTP)) = [];
    cExperiment.saveTimelapseExperiment(currentPos);   
    clear cTimelapse;
end
