function compileCellInformation(cExperiment,positionsToExtract,force,writelog)
% compileCellInformation(cExperiment,positionsToExtract,force)
%
% compiles all the data extracted in each of the positions into cellInf
% fiels of cExperiment.
%
% positionsToExtract    :   indexes of positions to extract. defaults to
%                           all positions.
% force                 :   boolean, defaults to false. If true, data is
%                           compiled even if the extraction parameters for
%                           each of the individual timelapses done match.
%                           Otherwise it will skip positions who's
%                           extractedData.extractionParameters field do not
%                           match those of the first timelapse.

if nargin<2 || isempty(positionsToExtract)
    positionsToExtract=find(cExperiment.posTracked);
end

if nargin<3 || isempty(force)
    
    force = false;
    
end

if nargin<4 || isempty(writelog)
    writelog = true;
end

% Start logging protocol
cExperiment.logger.start_protocol('compiling cell information',length(positionsToExtract));
try

cTimelapse=cExperiment.returnTimelapse(positionsToExtract(1));

cExperiment.cellInf=cTimelapse.extractedData;

cExperiment.cellInf(1).posNum=positionsToExtract(1)*ones(size(cExperiment.cellInf(1).trapNum));
if numel(cExperiment.cellInf)>1
    [cExperiment.cellInf(2:end).trapNum] = deal([]);
    [cExperiment.cellInf(2:end).cellNum] = deal([]);
end

% Track extracted timepoints for compilation of times from meta data:
extractedTimepoints = false(length(cExperiment.dirs),...
    length(cTimelapse.timepointsProcessed));
extractedTimepoints(positionsToExtract(1),:) = cTimelapse.timepointsProcessed;

% list of fields that are not identically sized arrays
fields_treated_special = {'posNum','trapNum','cellNum','mothers',...
    'extractionParameters','imQuantiles','nonTrapBackground'};
trapwise_fields = {'trapX','trapY','trapNCells',...
    'trapBgdMean','trapBgdMedian','trapBgdArea'};

if isfield(cExperiment.cellInf,'imQuantiles')
    % Preallocate imQuantiles for all positions
    for chi=1:numel(cExperiment.cellInf)
        cExperiment.cellInf(chi).imQuantiles(:,:,2:numel(positionsToExtract)) = NaN;
    end
end

if isfield(cExperiment.cellInf,'nonTrapBackground')
    % Preallocate nonTrapBackground for all positions
    for chi=1:numel(cExperiment.cellInf)
        cExperiment.cellInf(chi).nonTrapBackground(2:numel(positionsToExtract),:) = NaN;
    end
end

% size of chunks to preallocate.
tempLen=50e3;
field_names = fieldnames(cExperiment.cellInf);
fields_to_treat = setdiff(field_names,union(fields_treated_special,trapwise_fields));

trapwise_field_values = cell(numel(positionsToExtract),1);
trapwise_field_values{1} = rmfield(cTimelapse.extractedData,...
    setdiff(fieldnames(cTimelapse.extractedData),trapwise_fields));

data_template = [];
index = length(cExperiment.cellInf(1).posNum);
for posi=2:length(positionsToExtract)

    if index>= length(cExperiment.cellInf(1).posNum) || posi==2
        %preallocate more space
        for chi = 1:length(cExperiment.cellInf)
            for fi = 1:length(fields_to_treat)
                fn = fields_to_treat{fi};
                if issparse(cExperiment.cellInf(chi).(fn))    
                    if isempty(data_template)
                        %timepoints in the data
                        tps = size(cExperiment.cellInf(chi).(fn),2);
                        data_template = spalloc(tempLen,tps,0.5*tempLen*tps);
                        %assumes only 50 percent of timepoints will be present on
                        %average.
                    end
                    cExperiment.cellInf(chi).(fn)((index+1):(index+tempLen),:)=data_template;
                elseif isa(cExperiment.cellInf(chi).(fn),'int16')
                    tps = size(cExperiment.cellInf(chi).(fn),2);
                    cExperiment.cellInf(chi).(fn)((index+1):(index+tempLen),:)=-ones(tempLen,tps,'int16');
                elseif isa(cExperiment.cellInf(chi).(fn),'uint8')
                    tps = size(cExperiment.cellInf(chi).(fn),2);
                    cExperiment.cellInf(chi).(fn)((index+1):(index+tempLen),:)=zeros(tempLen,tps,'uint8');
                elseif isa(cExperiment.cellInf(chi).(fn),'uint16')
                    tps = size(cExperiment.cellInf(chi).(fn),2);
                    cExperiment.cellInf(chi).(fn)((index+1):(index+tempLen),:)=zeros(tempLen,tps,'uint16');
                elseif ~isempty(cExperiment.cellInf(chi).(fn))
                    error('unrecognised cellInf type');
                end
            end
        end
        cExperiment.cellInf(1).trapNum((index+1):(index+tempLen)) = zeros(1,tempLen);
        cExperiment.cellInf(1).cellNum((index+1):(index+tempLen)) = zeros(1,tempLen);
        cExperiment.cellInf(1).posNum((index+1):(index+tempLen)) = zeros(1,tempLen);
        if ismember('mothers',field_names)
            cExperiment.cellInf(1).mothers((index+1):(index+tempLen)) = zeros(1,tempLen);
        end
    end
    
    pos = positionsToExtract(posi);
    cTimelapse=cExperiment.returnTimelapse(pos);

    if ~isequaln(cTimelapse.extractedData(1).extractionParameters,cExperiment.cellInf(1).extractionParameters) && ~force
        logmsg(cExperiment,'Not compiling data %d, extraction parameters do not match.',pos);
        continue
    end
    
    if isempty(cTimelapse.extractedData(1).cellNum)
        continue
    end
    
    trapwise_field_values{posi} = rmfield(cTimelapse.extractedData,...
        setdiff(fieldnames(cTimelapse.extractedData),trapwise_fields));
    
    % Handle the case where different positions have different numbers of
    % extracted timepoints (assumes all timepointsProcessed are left aligned)
    ntps = numel(cTimelapse.timepointsProcessed);
    if  ntps > size(extractedTimepoints,2)
        extractedTimepoints(:,end+1:ntps) = false;
    end
    extractedTimepoints(pos,1:ntps) = cTimelapse.timepointsProcessed;
    extractedTimepoints(pos,:) = ...
        extractedTimepoints(positionsToExtract(1),:) & extractedTimepoints(pos,:);
    tps = find(extractedTimepoints(pos,:));
    
    num_cells = length(cTimelapse.extractedData(1).cellNum);
    for chi = 1:length(cExperiment.cellInf)
        for fi = 1:length(fields_to_treat)
            fn = fields_to_treat{fi};
            if isempty(cTimelapse.extractedData(chi).(fn)), continue; end
            cExperiment.cellInf(chi).(fn)((index+1):(index+num_cells),tps)=cTimelapse.extractedData(chi).(fn)(:,tps);
        end
        if isfield(cTimelapse.extractedData(chi),'imQuantiles') && ...
                isfield(cExperiment.cellInf(chi),'imQuantiles')
            cExperiment.cellInf(chi).imQuantiles(:,tps,posi) = ...
                cTimelapse.extractedData(chi).imQuantiles(:,tps);
        end
        if isfield(cTimelapse.extractedData(chi),'nonTrapBackground') && ...
                isfield(cExperiment.cellInf(chi),'nonTrapBackground')
            cExperiment.cellInf(chi).nonTrapBackground(posi,tps) = ...
                cTimelapse.extractedData(chi).nonTrapBackground(tps);
        end
    end
    cExperiment.cellInf(1).trapNum((index+1):(index+num_cells)) = cTimelapse.extractedData(1).trapNum;
    cExperiment.cellInf(1).cellNum((index+1):(index+num_cells)) = cTimelapse.extractedData(1).cellNum;
    cExperiment.cellInf(1).posNum((index+1):(index+num_cells)) = pos*ones(1,num_cells);
    if ismember('mothers',field_names)
        cExperiment.cellInf(1).mothers((index+1):(index+num_cells)) = cTimelapse.extractedData(1).mothers;
    end
    index = index + num_cells;

    cExperiment.cTimelapse = [];

end

%remove left over zeros from preallocation.
for chi=1:length(cExperiment.cellInf)
    for fi = 1:length(fields_to_treat)
        fn = fields_to_treat{fi};
        if isempty(cExperiment.cellInf(chi).(fn)), continue; end
        cExperiment.cellInf(chi).(fn)((index+1):end,:)=[];
    end
end
cExperiment.cellInf(1).trapNum((index+1):end) =[];
cExperiment.cellInf(1).cellNum((index+1):end) = [];
cExperiment.cellInf(1).posNum((index+1):end) = [];
if ismember('mothers',field_names)
    cExperiment.cellInf(1).mothers((index+1):end) = [];
end

% Compile trap-wise values
if all(cellfun(@(x) isfield(x,'trapNCells'),trapwise_field_values))
    npositraps = cellfun(@(x) size(x(1).trapNCells,1),trapwise_field_values);
    twTrapNum = arrayfun(@(n) 1:n,npositraps,'Uniform',false);
    cExperiment.cellInf(1).trapWiseTrapNum = [twTrapNum{:}]';
    twPosNum = arrayfun(...
        @(n,p) p*ones(1,n),npositraps(:),positionsToExtract(:),'Uniform',false);
    cExperiment.cellInf(1).trapWisePosNum = [twPosNum{:}]';
    common_fields = fieldnames(trapwise_field_values{1});
    for p=2:numel(trapwise_field_values)
        common_fields = intersect(common_fields,...
            fieldnames(trapwise_field_values{p}));
    end
    for f=1:numel(common_fields)
        tw_field = common_fields{f};
        for chi=1:numel(cExperiment.cellInf)
            fvals = cellfun(@(x) full(x(chi).(tw_field)),...
                trapwise_field_values,'uni',0);
            cExperiment.cellInf(chi).(tw_field) = vertcat(fvals{:});
        end
    end
end

if force
    [cExperiment.cellInf(:).extractionParameters] = deal([]);
end

% Compile meta data into the cellInf:

if cExperiment.logger.use_gui
    cExperiment.compileMetaData(extractedTimepoints,cExperiment.logger.progress_bar);
else
    cExperiment.compileMetaData(extractedTimepoints,false);
end

cExperiment.saveExperiment();

if isa(cExperiment,'babyExperimentOmero')
    %Tag dataset in Omero for archiving - we assume that if the data are worth
    %compiling then the experiment should be archived.
    %If is an babyExperimentOmero object - tag for archiving
    try
        addArchiveTag(cExperiment.id);
        if writelog
            writeOmeroLog('Added archive tag',['Experiment '
                cExperiment.metadata.experiment ' tagged for automated archiving: babyExperiment.compileCellInformation']);
        end
    catch err
        warning ('Error adding archive tag');
        if writelog
            writeOmeroLog('Error adding archive tag',['Archive tagging failed for experiment ' cExperiment.metadata.experiment ': ' err.message]);
        end
    end
end

% Finish logging protocol
cExperiment.logger.complete_protocol;
catch err
    cExperiment.logger.protocol_error;
    rethrow(err);
end

end
