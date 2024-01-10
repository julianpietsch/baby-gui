function extractCellInformation(cExperiment,positionsToExtract,doParameterGUI,extractParameters,cStats)
% extractCellInformation(cExperiment,positionsToExtract,doParameterGUI,extractParameters)
%
% function to extract data for each position.
%
% cExperiment           :   object of the babyExperiment class
% positionsToExtract    :   list of position to extract data from. Also the
%                           list of positions for which the extraction
%                           parameters will be set if they are changed.
%                           Defaults to all positionstracked in the
%                           experiment (so cExperiment.posTracked)
% doParameterGUI        :   boolean of whether to open an editing GUI that allows you
%                           to set various extraction parameters. Defaults to true. 
% extractParameters     :   a parameter structure that will be saved to:
%                                cTimelapse.extractionParameters
%                           for each of the positions in
%                           positionsToExtract. For form of the structure
%                           see babyTimelapse.defaultExtractionParameters
%
% Basically just call cTimelapse.extractCellData for each timepoint.
if nargin<2 || isempty(positionsToExtract)
    positionsToExtract=find(cExperiment.posTracked);
end

if nargin<3 || isempty(doParameterGUI)
    doParameterGUI = true;
end

if nargin<4, extractParameters = []; end

if nargin<5, cStats = []; end

if doParameterGUI
    if isempty(extractParameters)
        cTimelapse = cExperiment.loadCurrentTimelapse(positionsToExtract(1));
        extractParameters = cTimelapse.extractionParameters;
    end
    extractParameters = cExperiment.guiSetExtractParameters(extractParameters);
end

% only set the Extract Parameters if either parameters have been provided
% or the parameter setting GUI was invoked. Was felt best to set it all in
% one go - less time efficient but if it gets stopped the appropriate
% extractionParameters are saved.
if doParameterGUI || ~isempty(extractParameters)
    cExperiment.setExtractParameters(positionsToExtract,extractParameters)
end

% Start logging protocol
cExperiment.logger.start_protocol('extracting cell information',length(positionsToExtract));
try
    
for i=1:length(positionsToExtract)
    experimentPos=positionsToExtract(i);
    cTimelapse=cExperiment.returnTimelapse(experimentPos);
    cTimelapse.extractCellData(cStats);
    cExperiment.cTimelapse=cTimelapse;
    cExperiment.saveTimelapseExperiment(experimentPos,false);
end

% Parse the full log file as part of the extraction to obtain the actual
% times of acquisition for each time point:
if isempty(cExperiment.logger.shouldLog) || ~cExperiment.logger.shouldLog || ~cExperiment.logger.use_gui
    progress_bar = false;
else
    progress_bar = cExperiment.logger.progress_bar;
end

%cExperiment.parseLogFile([],progress_bar);
cExperiment.saveExperiment;

% Finish logging protocol
cExperiment.logger.complete_protocol;
catch err
    cExperiment.logger.protocol_error;
    rethrow(err);
end

end
