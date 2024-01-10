function selectTPToProcess(cExperiment,positionsToSet,tpToProcess)
% selectTPToProcess(cExperiment,positionsToSet,tpToProcess)
%
% tpToprocess  -  an array of timepoints to process. If empty selected by
%                 GUI. Best if it is continuous (i.e. x:y)
%
% sets the timepointsToProcess field of both cExperiment object and all its
% children cTimelapse objects. also sets their timepointsProcessed field to
% false for timepoints outside the range of timepointsToProcess. This is
% done because subsequent processing steps will only be applied to the
% timepoints listed as timepointsToProcess.

cTimelapse=cExperiment.returnTimelapse(length(cExperiment.dirs));
%Load last one in case it has fewer timepoints as sometimes happens in an interrupted experiment.

if nargin<2
    positionsToSet=1:length(cExperiment.dirs);
end

if nargin<3 || isempty(tpToProcess)    
    if ~isempty(cExperiment.timepointsToProcess)
        startTP=(cExperiment.timepointsToProcess(1));
        endTP=cExperiment.timepointsToProcess(end);
    else
        startTP=1;
        endTP=length(cTimelapse.cTimepoint);
        
    end
    
    num_lines=1;clear prompt; clear def;
    prompt(1) = {'In contrast with cropTP, this keeps the timepoints associated with the timelapse, but just notes that they should not be processed. This is generally prefferable to cropping. Starting Timepoint'};
    prompt(2) = {['Ending Timepoint (max TP ',num2str(length(cTimelapse.cTimepoint)),')']};
    dlg_title = 'Tp To Process';
    def(1) = {num2str(startTP)};
    def(2) = {num2str(endTP)};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    
    startTP = str2double(answer{1});
    endTP = str2double(answer{2});
    
    cExperiment.timepointsToProcess=startTP:endTP;
    
else
    cExperiment.timepointsToProcess= tpToProcess;
    startTP = min(tpToProcess);
    endTP = max(tpToProcess);
end

waitbar_h = waitbar(0);
for i=1:length(positionsToSet)
    currentPos=positionsToSet(i);
    cTimelapse=cExperiment.loadCurrentTimelapse(currentPos);
    cTimelapse.timepointsToProcess = cExperiment.timepointsToProcess;
    % set any elements of the timepointsProcessed field outside the range
    % of timepoints to be processed to false. and ensure its length is
    % endTP.
    
    % if new timepoint to process have been added - set these to false
    cTimelapse.timepointsProcessed(end+1:endTP) = false;
    
    % if some timpoints are no longer being processed - remove them from
    % the list of processed/unprocessed timepoints.
    cTimelapse.timepointsProcessed(endTP+1:end) = [];
    
    % this line is to mark timepoints as not processed (i.e. false) if they
    % are not in timepointsToProcess (normally is no longer processing
    % earlier timpoints).
    cTimelapse.timepointsProcessed(~ismember(1:length(cTimelapse.timepointsProcessed),cTimelapse.timepointsToProcess)) = false;
    cExperiment.cTimelapse = cTimelapse;
    cExperiment.saveTimelapseExperiment([],false);
    
    waitbar(i/length(positionsToSet),waitbar_h,sprintf('timepoints set for position %d of %d ...',i,length(positionsToSet)));
end
cExperiment.saveExperiment();
close(waitbar_h)