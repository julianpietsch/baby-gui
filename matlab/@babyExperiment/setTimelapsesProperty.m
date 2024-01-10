function setTimelapsesProperty(cExperiment,poses,property,value)
% setTimelapsesProperty(cExperiment,poses,property,value)
% 
% set the property 'property' to value 'value' in all the timelapses
% specified by poses.
%
% rather a dangerous function.
if isempty(poses)
    poses = 1:numel(cExperiment.dirs);
end

for posi = poses
    pos = posi;
    cExperiment.loadCurrentTimelapse(pos);
    cExperiment.cTimelapse.(property) = value;
    cExperiment.saveTimelapseExperiment([],false);
end

end