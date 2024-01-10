function extractionParameters = guiSetExtractParameters( cExperiment,extractionParameters )
%extractParameters = guiSetExtractParameters( cExperiment,extractParameters )
%
% a place for all the random GUI interfaces that have accumulated around
% extractParameters. If you want some sort of dialog box for your
% parameters extraction method, put it here in the switch/case function
%
% cExperiment              :   object of the babyExperiment class
% extractionParameters     :   parameters structure of the form:
%                              extractFunction    : function handle for function used in extraction
%                              functionParameters : structure of parameters taken by that function
%
% See also, BABYTIMELAPSE.EXTRACTCELLDATA,
% BABYTIMELAPSE.EXTRACTCELLDATASTANDARD,
% BABYEXPERIMENT.SETEXTRACTPARAMETERS

if nargin<2 || isempty(extractionParameters)
    
    extractionParameters = babyTimelapse.defaultExtractParameters;
end

% Though there is now no alternative to these 3 methods, i have left in the
% if statement so that it is clear how one might change the extraction
% process. 

% These 3 extraction methods share the same gui.
if isequal(extractionParameters.extractFunction,@extractCellDataStandardParfor) ...
        || isequal(extractionParameters.extractFunction,@extractCellDataStandard)...
        || isequal(extractionParameters.extractFunction,@extractCellParamsOnly)...
        
    functionParameters = extractionParameters.functionParameters;
    
    list = {'max','mean','std','sum','basic'};
    dlg_title = 'What to extract?';
    prompt = {['All Params using max projection (max), std (std), mean (mean) or sum (sum) of stacks; or basic (basic)' ...
        ' the basic measure only compiles the x, y locations of cells along with the estimated radius so it is much faster, but less informative.'],'','',''};
    answer = listdlg('PromptString',prompt,'Name',dlg_title,'ListString',list,'SelectionMode','single',...
        'ListSize',[300 100]);
    
    type=list{answer};
    switch type
        case {'max','std','mean','sum'}
            % bit of an ugly hack that should be sorted out, but means that
            % if someone selects basic they can later change their
            % selection to standard (i.e. basic extraction will run through
            % the gui).
            if isequal(extractionParameters.extractFunction,@extractCellParamsOnly)
                extractionParameters.extractFunction = @extractCellDataStandardParfor;
                functionParameters = struct;
            end
            functionParameters.type = type;
        case {'basic'}
            % if the extraction is basic there are not parameters, it
            % simply extracts the cell size and positions. 
            extractionParameters.extractFunction = @extractCellParamsOnly;
            functionParameters = [];
    end
    
    if ~strcmp(type,'basic')
        
        channel_list = cExperiment.channelNames;
        
        dlg_title = 'Which channels to extract?';
        prompt = {['please select the channels for which you would like to extract data'],'',''};
        answer = listdlg('PromptString',prompt,'Name',dlg_title,'ListString',channel_list,'SelectionMode','multiple',...
            'ListSize',[300 100]);
        functionParameters.channels = answer;
        
        % this is only used by the Standard method, not the
        % StandardParforMethod.
        if isequal(extractionParameters.extractFunction,@extractCellDataStandard)
            % this rather ugly call
            settings_dlg_struct = struct(...
                'title', 'nuclear label?',...
                'Description','If one of the channels is a nuclear label, please specify it here. This must be on of your extraction channels. If you have no particular marker please select '' not applicable '' ',...
                'forgotten_field_1',struct('entry_name',{{'nuclear tag field','nuclearChannel'}},'entry_value',{vertcat({'not applicable'},channel_list(functionParameters.channels))}),...
                'forgotten_field_2',struct('entry_name',{{'number of candidate nuclear pixels','maxAllowedOverlap'}},'entry_value',{25}),...
                'forgotten_field_3',struct('entry_name',{{'number of final nuclear pixels','maxPixOverlap'}},'entry_value',{5})...
                );
            
            answer_struct = settingsdlg(settings_dlg_struct);
            
            functionParameters.maxAllowedOverlap = answer_struct.maxAllowedOverlap;
            functionParameters.maxPixOverlap = answer_struct.maxPixOverlap;
            if strcmp(answer_struct.nuclearChannel,'not applicable')
                functionParameters.nuclearMarkerChannel = NaN;
            else
                functionParameters.nuclearMarkerChannel = find(strcmp(answer_struct.nuclearChannel,channel_list));
            end
        end
    end
    
    extractionParameters.functionParameters = functionParameters;
    
  
end

end

