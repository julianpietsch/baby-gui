function new_location = copyExperiment(cExperiment,new_location)
% new_location = copyExperiment(cExperiment,new_location)
%
% copies all the cTimelapse files to a new locaiton, changing the
% saveFolder of cExperiment so that it won't overwrite the original.
    

if nargin<2 || isempty(new_location)
    fprintf('\n\n   please select a location to save the copy of the experiment info     \n\n')
    new_location = uigetdir(cExperiment.saveFolder); 
end

for diri=1:length(cExperiment.dirs)
    
    copyfile([cExperiment.saveFolder filesep cExperiment.dirs{diri} 'cTimelapse.mat'],...
        [new_location filesep cExperiment.dirs{diri} 'cTimelapse.mat']);
    
end

% Log the copy command to the old log:
old_location = cExperiment.saveFolder;
logmsg(cExperiment,'\n=====================');
logmsg(cExperiment,'Experiment copied from:\n\t%s\n\tto:\n\t%s',old_location,new_location);
logmsg(cExperiment,'---------------------');

% Update the saveFolder and save the experiment:
cExperiment.saveFolder = new_location;
cExperiment.saveExperiment;

% Log the copy command to a new log file in the new location:
logmsg(cExperiment,'\n=====================');
logmsg(cExperiment,'Experiment copied from:\n\t%s\n\tto:\n\t%s',old_location,new_location);
logmsg(cExperiment,'---------------------');

end

