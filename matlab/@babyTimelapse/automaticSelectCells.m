function automaticSelectCells(cTimelapse,params)
% automaticSelectCells(cTimelapse,params)
% ----------------------------------------------------------
% automatically select cells to populate cTimelapse.cellsToPlot based on a
% number of criterion encapsulated in params (which are chosen by GUI if
% not provided). cellsToPlot is a sparse matrix with a 1 at 
% (trapNUM,cellLabel) when a cell in trap trapNUM and with label
% cellLabel meets the criteria specified by params.
%
% param fields:
%   params.fraction      -   fraction cells must be present
%     params.duration      -   duration for which cell must be present (same as above but number)
%     params.framesToCheck      -   cells must appear in first X frames
%     params.framesToCheckEnd      -   cells must be present in last X frames
%     params.maximumNumberOfCells      -   limit to this number of cells (handy if curating).
%     params.includeDaughters    - also select daughters of selected cells (default is false)
%
% used for selecting cells to extract data for and also in the
% combineTracklets code. 


if isempty(cTimelapse.timepointsProcessed)
    tempSize=[cTimelapse.cTimepoint.trapInfo];
    cTimelapse.timepointsProcessed=ones(1,length(tempSize)/length(cTimelapse.cTimepoint(cTimelapse.timepointsToProcess(1)).trapInfo));
end

if nargin<2
    params.fraction=.8; %fraction of timelapse length that cells must be present or
    params.duration=5; %number of frames cells must be present
%     params.cellsToCheck=4;
    params.framesToCheck=find(cTimelapse.timepointsProcessed,1,'last');
    params.framesToCheckEnd=find(cTimelapse.timepointsProcessed,1,'first');
    params.maximumNumberOfCells = Inf;
    params.includeDaughters = false;
    
    num_lines=1;clear prompt; clear def;
    prompt(1) = {'Fraction of whole timelapse a cell must be present'};
    prompt(2) = {'OR - number of frames a cell must be present'};
    prompt(3) = {'Cell must appear in the first X frames'};
    prompt(4) = {'Cell must be present after frame X'};
    prompt(5) = {'Select a maximum of X cells (useful if you want to check cells and not spend ages)'};
    prompt(6) = {'Include daughters of selected? (specify "true" if so)'};

    dlg_title = 'Tracklet params';    
    def(1) = {num2str(params.fraction)};
    def(2) = {num2str(params.duration)};
    def(3) = {num2str(params.framesToCheck)};
    def(4) = {num2str(params.framesToCheckEnd)};
    def(5) = {num2str(params.maximumNumberOfCells)};
    if params.includeDaughters, def(6) = {'true'}; else, def(6) = {'false'}; end
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    params.fraction=str2double(answer{1}); %fraction cells must be present
    params.duration=str2double(answer{2}); %duration for which cell must be present (same as above but number)
    params.framesToCheck=str2double(answer{3}); % cells must appear in first X frames
    params.framesToCheckEnd=str2double(answer{4}); % cells must be present in last X frames
    params.maximumNumberOfCells = str2double(answer{5}); % limit to this number of cells (handy if curating).
    params.includeDaughters = strcmpi(answer{6},'true'); % include any daughters of the selected cells
end

cTimelapse.cellsToPlot(:)=0;
if ~isfield(params, 'maximumNumberOfCells')
    params.maximumNumberOfCells = Inf;
end

if ~isfield(params,'includeDaughters')
    params.includeDaughters = false;
end

%keep track of number of cells found to ensure you don't go over
%maxmimumNumberOfCells.
cellsLeft = params.maximumNumberOfCells;

cTimelapse.cellsToPlot=zeros(length(cTimelapse.cTimepoint(cTimelapse.timepointsToProcess(1)).trapInfo),5000);
cTimepoint=cTimelapse.cTimepoint;
for trap=1:length(cTimelapse.cTimepoint(cTimelapse.timepointsToProcess(1)).trapInfo)
    %disp(['Trap Number ' int2str(trap)]); % SUPERSEDED BY babyLogging
    cellLabels=zeros(1,100*sum(cTimelapse.timepointsProcessed));% list of each occurrence of a cellLabel upto timepoint params.framesToCheck
    cellLabelsEnd=zeros(1,100*sum(cTimelapse.timepointsProcessed));% list of each occurrence of a cellLabel after timepoint params.framesToCheckEnd
    cellsSeen=zeros(1,100*sum(cTimelapse.timepointsProcessed));% list of each occurrence of a cellLabel before timepoint params.framesToCheckEnd
    index=0;
    for timepoint=cTimelapse.timepointsToProcess
        if cTimelapse.timepointsProcessed(timepoint)
            tempLabels=cTimepoint(timepoint).trapInfo(trap).cellLabel;
            cellLabels(1,index+1:index+length(tempLabels))=tempLabels;
            if timepoint<=params.framesToCheck
                cellsSeen(1,index+1:index+length(tempLabels))=tempLabels;
            end
            if timepoint>=params.framesToCheckEnd
                cellLabelsEnd(1,index+1:index+length(tempLabels))=tempLabels;
            end
            index=index+length(tempLabels);
        end
    end
    cellLabels=cellLabels(1:index);
    cellLabelsEnd(cellLabelsEnd==0)=[];
    cellsSeen(cellsSeen==0) = [];
    
    % tabulate the number of times (and therefore timepoints) each
    % cell label occurs. n(1) is an artifact of the hist call and is
    % removed.
    n=hist(cellLabels,0:max(cellLabels));
    if ~isempty(n)
        n(1)=[];
    end
    cellsSeen = unique(cellsSeen);
    cellLabelsEnd = unique(cellLabelsEnd);
    locs=find(n>=sum(cTimelapse.timepointsProcessed)*params.fraction | n>=params.duration);
    
    if ~isempty(cellsSeen) && ~isempty(locs)
        % pick out those cells (locs) that are in both cellSeen (those
        % cells in the first n timepoints) and cellLabelEnd (those in the
        % last m timepoints).
        locs=locs(ismember(locs,intersect(cellsSeen,cellLabelsEnd)));
        if length(locs)>cellsLeft
            locs = locs(1:cellsLeft);
            cellsLeft = 0;
        else
            cellsLeft = max((cellsLeft - length(locs)),0);        
        end
        if ~isempty(locs)
            if params.includeDaughters
                % Also select any cells that are daughters
                daughters = find(ismember(cTimelapse.cellMothers(trap,:),locs));
                locs = unique([locs(:);daughters(:)]);
            end
            for cellsForPlot=1:length(locs)
                cTimelapse.cellsToPlot(trap,locs(cellsForPlot))=1;
            end
        end
    end
end
