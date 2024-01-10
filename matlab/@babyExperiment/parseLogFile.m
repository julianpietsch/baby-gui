function fail_flag = parseLogFile(cExperiment,logFile,progress_bar)
%parseLogFile Parse experiment's log file to obtain meta data
%   cExperiment: this object
%   logFile (optional): if specified, that log file is used, otherwise
%       a good guess is made by searching 'cExperiment.rootFolder'
%   progress_bar (optional): specify either a parent 'Progress' object in
%       which to create a child progress bar, or 'meta_only' to process 
%       only the header and Acq file and skip processing the timepoints, or
%       a logical false value to disable the GUI.
%
%   fail_flag : boolean
%       an indicator of success. Is set to true if either the log file or
%       the acq file could not be found. Consider refining.

fail_flag = false;

if nargin<2 || isempty(logFile)
    if isa(cExperiment,'babyExperimentOmero') && ...
            isa(cExperiment.OmeroDatabase,'OmeroDatabase') && ...
            cExperiment.OmeroDatabase.sessionActive
        % Ensure all log files (including the acq file) are already
        % downloaded to the saveFolder (NB: downloadFiles will not
        % re-download if they are already present and up-to-date):
        filePaths = cExperiment.OmeroDatabase.downloadFiles(...
            cExperiment.omeroDs,[],cExperiment.saveFolder);
        logFile = filePaths(1).log;
    else
        logDirs = {cExperiment.rootFolder,cExperiment.saveFolder};
        for d=1:length(logDirs)
            logFile = dir(fullfile(logDirs{d},'*log.txt'));
            logFile = logFile(~strcmp({logFile.name},'cExperiment_log.txt'));
            % Break this loop as soon as we find a suitable candidate
            if ~isempty(logFile), break; end
        end
        if isempty(logFile)
            warning('The log file could not be found. Skipping parseLogFile...');
            fail_flag = true;
            return
        else
            if length(logFile)>1
                warning(['More than one log file available in "%s". ',...
                    'Using first found...'],logDirs{d});
            end
            logFile = fullfile(logDirs{d},logFile(1).name);
        end
    end
end

close_progress = false;
pop_progress_bar = false;
meta_only = false;

if nargin<3 
    progress_bar = true;
end

if ~isa(progress_bar,'Progress')
    if ischar(progress_bar) && strcmp(progress_bar,'meta_only')
        progress_bar = [];
        meta_only = true;
    elseif islogical(progress_bar) && ~progress_bar
        progress_bar = [];
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

[logdir,~,~] = fileparts(logFile);

%% Initialise default header variables and regular expressions

microscope = ''; acqfile = ''; acqdir = '';
comments = ''; project = ''; tags = {};

% Specify the regular expressions to test:
headerRegExp = {'^Microscope name is:',...
    '^Acquisition settings are saved in:','^Experiment details:',...
    '^Omero project:','^Omero tags:','^Experiment started at:'};
parseMicroscope = '^Microscope name is: (.*)$';
headerExpNames = {'microscope','acqfile','comments',...
    'omeroproj','omerotags','endheader'};
header = true; acqfileline = false;
commentlines = 0; projectline = false; tagline = false;

%% Initialise time course variables and regular expressions

% Pre-allocate arrays by guessing sizes
approxNtimepoints = 100;
approxNpositions = length(cExperiment.dirs);
%NB: times gets an updated size if we make it past the header lines
times = zeros(approxNtimepoints,approxNpositions);
positionStrs = cell(approxNpositions,1);
positionExposure = struct();
pumpSwitches = zeros(0,1);
acq = {}; % Gets filled with structure parsed from Acq file
timepoint = 0;
npos = 0;

parseRegExpFast = {'^------Time point_','^Position:',...
    '^Channel:','^Exposure time:','^Switching pumps at'};
parseRegExp = {'^------Time point_(\d)+-','^Position:(\d+),(.*)$',...
    '^Channel:(.*) set at:(.*)$','^Exposure time:(\d+)ms$',...
    '^Switching pumps at (.*)[.] Fast'};
parseExpNames = {'timepoint','position','channel','exposure','pumpswitch'};

dateFormat = 'dd-mmm-yyyy HH:MM:SS';

%% Parse the log file line by line:
fprintf('Parsing log file:\n%s\n',logFile);
fid = fopen(logFile);
tline = fgets(fid);

while ischar(tline)
    tline = strtrim(tline);
    
    if header
        if acqfileline
            [acqdir,filename,ext] = fileparts(tline);
            acqfile = [filename,ext];
            % Attempt to find and parse the Acq file:
            if exist(fullfile(logdir,acqfile),'file')
                acq = parseAcqFile(fullfile(logdir,acqfile));
            else
                fail_flag = true;
            end
            acqfileline = false;
        end
        if projectline
            project = tline; projectline = false;
        end
        if tagline
            tags = regexp(tline,'([^,]+),','tokens');
            tags = [tags{:}]; tagline = false;
        end
        
        % Parse the current header line:
        matchline = regexp(tline,headerRegExp,'once','start');
        matches = ~cellfun(@isempty,matchline);
        if any(matches)
            lineType = headerExpNames(matches);
            switch lineType{1}
                case 'microscope'
                    parsedline = regexp(tline,parseMicroscope,'once','tokens');
                    microscope = parsedline{1};
                case 'acqfile'
                    acqfileline = true;
                case 'comments'
                    acqfileline = false; commentlines = 1;
                case 'omeroproj'
                    commentlines = 0; projectline = true;
                case 'omerotags'
                    projectline = false; tagline = true;
                case 'endheader'
                    header = false; acqfileline = false; tagline = false;
                    commentlines = 0; projectline = false;
                    % Attempt to update the approximated number of timepoints 
                    if ~isempty(acq) && ~isempty(acq.times) && ~isempty(acq.times.ntimepoints)
                      approxNtimepoints = acq.times.ntimepoints;
                      times = zeros(approxNtimepoints,approxNpositions);
                    elseif ~isempty(cExperiment.cTimelapse)
                      approxNtimepoints = length(cExperiment.cTimelapse.cTimepoint);
                      times = zeros(approxNtimepoints,approxNpositions);
                    elseif ~isempty(cExperiment.cellInf) && isfield(cExperiment.cellInf,'xloc')
                      approxNtimepoints = size(cExperiment.cellInf(1).xloc,2);
                      times = zeros(approxNtimepoints,approxNpositions);
                    end
            end
        end
        
        if commentlines>0
            if commentlines==2
                comments = tline;
            elseif commentlines>2
                comments = strjoin({comments,tline},'\n');
            end
            commentlines = commentlines + 1;
        end
    end
    
    % Parse the current line:
    matchline = regexp(tline,parseRegExpFast,'once','start');
    matches = ~cellfun(@isempty,matchline);
    
    if any(matches)
        lineType = parseExpNames(matches);
        lineRegExp = parseRegExp(matches);
        switch lineType{1}
            case 'timepoint'
                parsedline = regexp(tline,lineRegExp{1},'once','tokens');
                timepoint = str2double(parsedline{1});
                % For 'meta_only' flag, parse only to the second timepoint
                if meta_only && timepoint==2
                    break
                end
                fprintf('.');
                if mod(timepoint,60)==0
                    fprintf('\n');
                end
                % If we see a timepoint, begin progress loop over timepoints:
                if ~pop_progress_bar && ~isempty(progress_bar)
                    progress_bar.push_bar('Parsing log file...',1,approxNtimepoints);
                    pop_progress_bar = true;
                end
                if ~isempty(progress_bar)
                    progress_bar.set_val(timepoint);
                end
            case 'position'
                parsedline = regexp(tline,lineRegExp{1},'once','tokens');
                position = str2double(parsedline{1});
                npos = max([position,npos]);
                positionStrs{position} = parsedline{2};
            case 'channel'
                parsedline = regexp(tline,lineRegExp{1},'once','tokens');
                channelName = parsedline{1};
                if ~isfield(positionExposure,channelName)
                    positionExposure.(channelName) = zeros(approxNpositions,1);
                end
                if timepoint>size(times,1) || position>size(times,2) || ...
                        times(timepoint,position)==0
                    times(timepoint,position) = datenum(parsedline{2},dateFormat);
                end
            case 'exposure'
                if position>length(positionExposure.(channelName)) || ...
                        positionExposure.(channelName)(position) == 0
                    parsedline = regexp(tline,lineRegExp{1},'once','tokens');
                    positionExposure.(channelName)(position) = ...
                        str2double(parsedline{1});
                end
            case 'pumpswitch'
                parsedline = regexp(tline,lineRegExp{1},'once','tokens');
                pumpSwitches(end+1) = datenum(parsedline{1},dateFormat);
        end
    end
    
    tline = fgets(fid);
end

fprintf('done.\n');

fclose(fid);

if timepoint==0
    warning('The specified log file did not contain any time points. "metadata" property not set.');
    if pop_progress_bar
        progress_bar.pop_bar; % finished reading all timepoints
    end
    if close_progress
        progress_bar.frame.dispose;
    end
    return
end

%% Save the date
date = datestr(times(1,1),'yyyy_mm_dd');

if ~meta_only
    
    %% Convert the pump switching times to minutes:
    if ~isempty(pumpSwitches)
        relPumpSwitches = pumpSwitches - times(1,1);
        pumpSwitches = zeros(size(relPumpSwitches));
        for i = 1:length(relPumpSwitches)
            pumpSwitches(i) = sum(datevec(relPumpSwitches(i)) .* [0 0 24*60 60 1 1/60]);
        end
    end
    
    %% Convert timepoint times to minutes
    times = times - times(1,1);
    timesInMinutes = zeros(size(times));
    for i = 1:size(times,1)
        for j = 1:size(times,2)
            timesInMinutes(i,j) = sum(datevec(times(i,j)) .* [0 0 24*60 60 1 1/60]);
        end
    end
    
    %% Subset the timepoints and positions to those recorded in log file
    
    % The final value of timepoint is the last in the log file, and npos is the
    % maximum position number seen in the log file:
    timesInMinutes = timesInMinutes(1:timepoint,1:npos)';
    positionStrs = positionStrs(1:npos);
    channels = fieldnames(positionExposure);
    for i=1:length(channels)
        positionExposure.(channels{i}) = positionExposure.(channels{i})(1:npos);
    end

end

%% Update cExperiment:
experiment = regexp(acqfile,'^(.*)Acq\.txt$','tokens','once');
username = regexp(acqdir,'Swain Lab[/\\]([^/\\]+)[/\\]RAW DATA','tokens','once');

annotations = struct();
annotations.experiment = experiment{1};
annotations.username = username{1};
annotations.microscope = microscope;
annotations.logFile = logFile;
annotations.acqFile = fullfile(acqdir,acqfile);
annotations.comments = comments;
annotations.project = project;
annotations.tags = tags;
annotations.date = date;
if ~meta_only
    annotations.pumpSwitches = pumpSwitches;
    annotations.logTimes = timesInMinutes;
end
annotations.logPosNames = positionStrs;
annotations.logExposureTimes = positionExposure;




annotations.acq = acq;

fields_to_replace = {'pumpSwitches','logTimes','logPosNames','logExposureTimes'};
if isempty(cExperiment.metadata)
    cExperiment.metadata = annotations;
else
    fields_to_update = fieldnames(annotations);
    for j=1:length(fields_to_update)
        field = fields_to_update{j};
        if ismember(field,fields_to_replace) || ...
                ~isfield(cExperiment.metadata,field) || ...
                isempty(cExperiment.metadata.(field))
            cExperiment.metadata.(field) = annotations.(field);
        end
    end
end

% Clean up progress bar:
if pop_progress_bar
    progress_bar.pop_bar; % finished reading all timepoints
end
if close_progress
    progress_bar.frame.dispose;
end

end


%% parseAcqFile function

function acqannot = parseAcqFile(acqfile)
%parseAcqFile Parse experiment's Acq file to obtain experiment settings
%   acqfile: the name of the Acq file to be parsed

%% Initialise default variables and regular expressions

channelTable = table({},[],[],[],[],[],[],[],'VariableNames',...
    {'names','exposure','skip','zsect','start','mode','gain','voltage'});
zsectStruct = struct('sections',{},'spacing',{},'PFSon',{},'anyz',{},...
    'drift',{},'method',{});
timeStruct = struct('istimelapse',{},'interval',{},'ntimepoints',{},...
    'totalduration',{});
pointsTable = table({},[],[],[],[],[],'VariableNames',...
    {'name','xpos','ypos','zpos','PFSoffset','group'});
npumps = 0;
pumpStartTable = table({},[],[],{},[],{},'VariableNames',...
    {'port','diameter','rate','direction','isrunning','contents'});
switchParams = struct('volume',[],'rate',[],'nchanges',[],...
    'switchtimes',[],'switchto',[],'switchfrom',[],...
    'pumpflow',[]);

% Specify default values for each item in the data rows:
channelDefaults = [{''},repmat({NaN},1,7)];
zsectDefaults = repmat({NaN},1,6);
timeDefaults = repmat({NaN},1,4);
pointsDefaults = [{''},repmat({NaN},1,5)];
pumpStartDefaults = {'',NaN,NaN,'',NaN,''};

% Specify the regular expressions to test:
sectionRegExp = {'^Channels:$','^Channel name,',...
    '^Z_sectioning:$','^Sections,Spacing,',...
    '^Time_settings:$','^Points:$','^Position name, X position',...
    '^Flow_control:$','^Syringe pump details:',...
    '^Pump states at beginning of experiment:$','^Pump port, Diameter,',...
    '^Dynamic flow details:$','^Number of pump changes:',...
    '^Switching parameters:','^Infuse/withdraw volumes:$',...
    '^Infuse/withdraw rates:$','^Times:','^Switched to:',...
    '^Switched from:','^Flow post switch:$'};
sectionNames = {'channels','channels','zsect','zsect','times',...
    'positions','positions','flow','npumps','pumpstart','pumpstart',...
    'pumpstartend','nchanges','switchparams','switchvol','switchrate',...
    'switchtimes','switchto','switchfrom','switchflow'};

npumpsRegExp = '^Syringe pump details: (\d+) pumps.$';
nchangesRegExp = '^Number of pump changes:(\d+)$';
switchparamsRegExp = '^Switching parameters:(\d+),(\d+)$';
switchtimesRegExp = '^Times:(\d+.*)$';
switchtoRegExp = '^Switched to:(\d+.*)$';
switchfromRegExp = '^Switched from:(\d+.*)$';

activeSection = '';
sectionLineNum = 0;

%% Parse the Acq file line by line:
fprintf('...and also parsing Acq file:\n%s\n',acqfile);
fid = fopen(acqfile);
tline = fgets(fid);

while ischar(tline)
    % Trim unwanted whitespace and check that the line isn't empty:
    tline = strtrim(tline);
    if isempty(tline)
        tline = fgets(fid);
        continue
    end
    sectionLineNum = sectionLineNum + 1;
    
    % Parse the current line:
    matchline = regexp(tline,sectionRegExp,'once','start');
    matches = ~cellfun(@isempty,matchline);
    
    if any(matches)
        activeSection = sectionNames(matches);
        activeSection = activeSection{1};
        sectionLineNum = 0;
        
        switch activeSection
            case 'npumps'
                parsedline = regexp(tline,npumpsRegExp,'once','tokens');
                npumps = str2double(parsedline{1});
                activeSection = 'pumpstart'; % Included since there is not always a header
            case 'nchanges'
                parsedline = regexp(tline,nchangesRegExp,'once','tokens');
                if ~isempty(parsedline)
                    switchParams.nchanges = str2double(parsedline{1});
                end
            case 'switchparams'
                parsedline = regexp(tline,switchparamsRegExp,'once','tokens');
                if ~isempty(parsedline)
                    switchParams.volume = str2double(parsedline{1});
                    switchParams.rate = str2double(parsedline{2});
                end
            case 'switchtimes'
                parsedline = regexp(tline,switchtimesRegExp,'once','tokens');
                if ~isempty(parsedline)
                    switchParams.switchtimes = parseNumList(parsedline{1});
                end
            case 'switchto'
                parsedline = regexp(tline,switchtoRegExp,'once','tokens');
                if ~isempty(parsedline)
                    switchParams.switchto = parseNumList(parsedline{1});
                end
            case 'switchfrom'
                parsedline = regexp(tline,switchfromRegExp,'once','tokens');
                if ~isempty(parsedline)
                    switchParams.switchfrom = parseNumList(parsedline{1});
                end
        end
    else
        switch activeSection
            case 'channels'
                channelTable = addRowFromLine(channelTable,tline,channelDefaults);
                % Add columns to positions list and defaults:
                channelName = channelTable.names{end};
                count = sum(ismember(channelTable.names,channelName));
                if count>1
                    channelName = sprintf('%s_%u',channelName,count);
                end
                pointsTable = horzcat(pointsTable,table([],...
                    'VariableNames',{channelName}));
                pointsDefaults = [pointsDefaults,{NaN}];
            case 'zsect'
                zsectStruct = ...
                    cell2struct(parseListWithDefaults(tline,zsectDefaults),...
                    fieldnames(zsectStruct),2);
            case 'times'
                timeStruct = ...
                    cell2struct(parseListWithDefaults(tline,timeDefaults),...
                    fieldnames(timeStruct),2);
            case 'positions'
                pointsTable = addRowFromLine(pointsTable,tline,pointsDefaults);
            case 'pumpstart'
                pumpStartTable = addRowFromLine(pumpStartTable,tline,...
                    pumpStartDefaults);
            case 'nchanges'
                switchParams.nchanges = str2double(tline);
            case 'switchvol'
                switchParams.volume = parseNumList(tline);
            case 'switchrate'
                switchParams.rate = parseNumList(tline);
            case 'switchtimes'
                switchParams.switchtimes = parseNumList(tline);
            case 'switchto'
                switchParams.switchto = parseNumList(tline);
            case 'switchfrom'
                switchParams.switchfrom = parseNumList(tline);
            case 'switchflow'
                switchParams.pumpflow(sectionLineNum,:) = parseNumList(tline);
        end
    end
    
    tline = fgets(fid);
end

fclose(fid);

%% Compile final output structure:
acqannot = struct();
acqannot.channels = channelTable;
acqannot.zsections = zsectStruct;
acqannot.times = timeStruct;
acqannot.positions = pointsTable;
acqannot.npumps = npumps;
acqannot.pumpInit = pumpStartTable;
acqannot.switchParams = switchParams;

    function numOut = parseNumList(numStrIn)
        numOut = cellfun(@str2double,strsplit(numStrIn,','));
    end

    function listout = parseListWithDefaults(rowstr,defaults)
        listout = strsplit(rowstr,',');
        listout = cellfun(@strtrim,listout,'UniformOutput',false);
        numtokens = cellfun(@isnumeric,defaults);
        if length(listout)<length(defaults)
            listout = [listout,defaults(length(listout)+1:end)];
            numtokens(length(listout)+1:end) = false;
        elseif length(listout)>length(defaults)
            listout = listout(1:length(defaults));
            warning('More elements than expected in a row of the acquisition file');
        end
        listout(numtokens) = ...
            cellfun(@str2double,listout(numtokens),'UniformOutput',false);
    end

    function tblout = addRowFromLine(tblin,rowstr,defaults)
        rowlist = parseListWithDefaults(rowstr,defaults);
        tblrow = cell2table(rowlist,'VariableNames',tblin.Properties.VariableNames);
        tblout = vertcat(tblin,tblrow);
    end
end
