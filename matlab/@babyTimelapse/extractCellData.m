function extractCellData(cTimelapse,cStats)
% extractCellData(cTimelapse)
%
% wrapper function that performs the extraction according to
% cTimelapse.extractionParameters. 
% Simply applies the function:
%   cTimelapse.extractionParameters.extractFunction
% which will usually make use of the parameters:
%   cTimelapse.extractionParameters.functionParameters
%
% sets the extractionParameters to be a field in the extractedData.

if nargin<2, cStats = []; end
do_cStats = ~isempty(cStats) && isa(cStats,'experimentSampleStats');

if do_cStats
    cTimelapse.extractionParameters.functionParameters.cStats = cStats;
end
    
cTimelapse.extractionParameters.extractFunction(cTimelapse)

if do_cStats
    cTimelapse.extractionParameters.functionParameters.cStats = [];
end

cTimelapse.extractedData(1).extractionParameters = cTimelapse.extractionParameters;

end

