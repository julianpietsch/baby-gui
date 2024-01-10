function compileChambers(cExperiment,varargin)
%BABYEXPERIMENT.COMPILECHAMBERS Save extracted data to grData files
%
%   Each group of positions specified in the cExperiment.chamberMap struct
%   is compiled and saved as a grData object.
%   Additional parameters are passed to the grData constructor.
%   Some behaviours are changed:
%   - all extracted fields are included by default
%   - the 'cellfilt' arg should be specified as a function handle that
%     accepts a cellInf as input and returns a valid cell filter as output
%   - full_daughters is true by default
%   - pixel_size is taken from cExperiment.pixelSize
ip = inputParser;
ip.addParameter('cellfilt',[],@(x) isempty(x) || isa(x,'function_handle'));
ip.addParameter('extract','all',@(x) isempty(x) || isequal(x,'all') || ...
    (iscellstr(x) || isstring(x)));
ip.addParameter('full_daughters',true,@(x) isscalar(x) && islogical(x));
ip.KeepUnmatched = true;
ip.parse(varargin{:});

cellfilt_fun = ip.Results.cellfilt;
extract = ip.Results.extract;
full_daughters = ip.Results.full_daughters;
unmatched = [fieldnames(ip.Unmatched),struct2cell(ip.Unmatched)]';

cnames = fieldnames(cExperiment.chamberMap);
cellInf_original = cExperiment.cellInf;
for f=1:numel(cnames)
    cname = cnames{f};
    poses = cExperiment.chamberMap.(cname);
    cExperiment.compileCellInformation(poses);
    cellInf = cExperiment.cellInf;
    if isfield(cellInf,'extractionParameters') ...
            && isfield(cellInf(1).extractionParameters.functionParameters,'channels')
        chnames = cellInf(1).extractionParameters.functionParameters.channels;
    else
        chnames = strcat('channel',arrayfun(@num2str,1:numel(cellInf),'Uniform',false));
    end
    to_extract = extract;
    if isequal(extract,'all')
        if numel(cellInf)>1
            fnames = fieldnames(cellInf);
            fempty = arrayfun(@(s) structfun(@(x) ...
                isempty(x) || (isnumeric(x) && all(x(:)==0)),s),...
                cellInf,'Uniform',false);
            fempty = [fempty{:}];
            isgeneral = all(fempty(:,2:end),2);
            generalFields = fnames(isgeneral);
            channelFields = fnames(~isgeneral);
        else
            channelFields = {};
            generalFields = fieldnames(cellInf);
        end
        to_extract = cell(1+numel(chnames),1);
        to_extract{1} = generalFields;
        for ch=1:numel(chnames)
            to_extract{1+ch} = strcat(chnames{ch},'.',channelFields);
        end
        to_extract = vertcat(to_extract{:});
    end
    cellfilt = [];
    if isa(cellfilt_fun,'function_handle')
        cellfilt = cellfilt_fun(cellInf);
    end
    
    fprintf('Calculating growth rates for chamber %s...',cname);
    chamberResults = grData(cellInf,'cellfilt',cellfilt,...
        'extract',to_extract,'full_daughters',full_daughters,...
        'pixel_size',cExperiment.pixelSize,unmatched{:});
    fprintf('done.\n');
    
    fprintf('Saving data for chamber %s...',cname);
    chamberResults.filename = fullfile(cExperiment.saveFolder,...
        sprintf('chamberResults_%s.h5',cname));
    chamberResults.saveData;
    fprintf('done.\n');
end

cExperiment.cellInf = cellInf_original;
cExperiment.saveExperiment;
