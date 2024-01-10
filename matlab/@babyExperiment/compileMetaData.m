function compileMetaData(cExperiment,extractedTimepoints,progress_bar)
%babyExperiment.compileMetaData Compile meta data into cellInf
%   This function adds the meta data to the cellInf array (only for the
%   first channel), including the generation of a times field that provides
%   a time of acquisition for each cell at each time point.

% First check that this cExperiment has been compiled:
if ~isprop(cExperiment,'cellInf') || isempty(cExperiment.cellInf)
    error('The cellInf hasn''t been compiled for this cExperiment.');
end

use_progress = true;
close_progress = false;
pop_progress_bar = false;

if nargin<3, progress_bar = true; end
if ~isa(progress_bar,'Progress')
    if islogical(progress_bar) && ~progress_bar
        use_progress = false;
    else
        % Initialise a progress bar
        progress_bar = Progress();
        % Centre the dialog box
        progress_bar.frame.setLocationRelativeTo([]);
        % Set the progress bar title
        progress_bar.frame.setTitle('Compiling meta data...');
        close_progress = true;
    end
end

% Check that the log file has been parsed for this cExperiment (and not
% just parsed with the 'meta_only' flag):
if isempty(cExperiment.posTimes)
    cExperiment.parseLogFile([],progress_bar);
end

% If we still don't have meta data, silently fail
if isempty(cExperiment.posTimes)
    if close_progress
        progress_bar.frame.dispose;
    end
    return
end

approxNtimepoints = size(cExperiment.posTimes,2);

%% Determine the timepoints for which extraction was performed if not provided:

% assumes each position provides at least on cell and throws error otherwise - potentially
% problematic.
extractedPositions = sort(unique(cExperiment.cellInf(1).posNum));
if nargin<2 || isempty(extractedTimepoints)
    extractedTimepoints = false(max(extractedPositions(:)),approxNtimepoints);
    
    try
        % First run a test call to see if the cExperiment can access cTimelapse
        cExperiment.returnTimelapse(extractedPositions(1));
        
        if use_progress
            progress_bar.push_bar('Determining extracted timepoints...',1,...
                length(extractedPositions));
            pop_progress_bar = true;
        end
        
        for i=1:length(extractedPositions)
            if use_progress, progress_bar.set_val(i); end
            ipos=extractedPositions(i);
            cTimelapse = cExperiment.returnTimelapse(ipos);
            extractedTimepoints(ipos,1:length(cTimelapse.timepointsProcessed)) = ...
                cTimelapse.timepointsProcessed;
            delete(cTimelapse);
        end
        
    catch
        warning('Cannot load cTimelapse for this cExperiment. Using cExperiment timepointsToProcess.');
        extractedTimepoints(:,cExperiment.timepointsToProcess) = true;
    end
else
    if size(extractedTimepoints,1) < max(extractedPositions(:))
        warning('"extractedTimepoints" has too few positions. Padding with false.');
        extractedTimepointsArg = extractedTimepoints;
        extractedTimepoints = false(max(extractedPositions(:)),...
            size(extractedTimepointsArg,2));
        extractedTimepoints(1:size(extractedTimepointsArg,1),:) = ...
            extractedTimepointsArg;
    end
end

%% Refactor the arrays to the same format as cellInf arrays:

% Truncate arrays to have the same time points as cellInf
% As per compileCellInformation, use the first field that is not considered
% special to determine the 'data_template' size:
fields_treated_special = {'posNum','trapNum','cellNum','extractionParameters'};
cellInfFields = fieldnames(cExperiment.cellInf);
template_field = find(~ismember(cellInfFields,fields_treated_special),1);
Ntimepoints = approxNtimepoints; % The number of time points in log file
if ~isempty(template_field)
    Ntimepoints = size(cExperiment.cellInf(1).(cellInfFields{template_field}),2);
end
if size(extractedTimepoints,2) < Ntimepoints
    warning('"extractedTimepoints" has too few time points. Padding with false.');
    extractedTimepoints = [extractedTimepoints,...
        zeros(size(extractedTimepoints,1),Ntimepoints-size(extractedTimepoints,2))];
else
    extractedTimepoints = extractedTimepoints(:,1:Ntimepoints);
end
timesInMinutes = cExperiment.posTimes;
if size(timesInMinutes,2) < Ntimepoints
    warning('The log file has too few times. Padding with zeros.');
    timesInMinutes = [timesInMinutes,...
        zeros(size(timesInMinutes,1),Ntimepoints-size(timesInMinutes,2))];
else
    timesInMinutes = timesInMinutes(:,1:Ntimepoints);
end

% Set times and extractedTimepoints based on each cell's position
positions = cExperiment.cellInf(1).posNum;
extractedTimepoints = extractedTimepoints(positions,:);
timesInMinutes = timesInMinutes(positions,:);

%% Update cExperiment:
cExperiment.cellInf(1).extractedTimepoints = extractedTimepoints;
cExperiment.cellInf(1).date = cExperiment.metadata.date;
cExperiment.cellInf(1).times = timesInMinutes;
cExperiment.cellInf(1).timepointsNotProcessed = sparse(~extractedTimepoints);
cExperiment.cellInf(1).annotations = cExperiment.metadata;

% Clean up progress bar:
if pop_progress_bar
    progress_bar.pop_bar; % finished reading all timepoints
end
if close_progress
    progress_bar.frame.dispose;
end

end
