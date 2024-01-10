function saveExperiment(cExperiment)
% SAVEEXPERIMENT(cExperiment)
% save the experiment in the saveFolder.

save([cExperiment.saveFolder filesep 'cExperiment.mat'],'cExperiment');

end


