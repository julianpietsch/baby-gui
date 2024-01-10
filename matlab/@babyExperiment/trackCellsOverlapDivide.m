function trackCellsOverlapDivide(cExperiment,poses,varargin)
%trackCellsOverlapDivide Track cells according to overlaps with previous tp
%
% 
if nargin<2 || isempty(poses)
    poses = 1:numel(cExperiment.dirs);
end

% Start logging protocol
cExperiment.logger.start_protocol('tracking cells by overlaps',length(poses));
try

for posi = 1:numel(poses)
    pos = poses(posi);
    cExperiment.loadCurrentTimelapse(pos);
    cExperiment.cTimelapse.trackCellsOverlapDivide(varargin{:});
    cExperiment.saveTimelapseExperiment([],false);
    cExperiment.posTrapsTracked(pos) = true;
end

cExperiment.saveExperiment;

% Finish logging protocol
cExperiment.logger.complete_protocol;

catch err
    cExperiment.logger.protocol_error;
    rethrow(err);
end

end