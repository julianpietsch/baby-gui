function cTimelapse=saveTimelapse(cTimelapse)
%Save function for Omero timelapses - can be used when the relevant cExperiment
%is not available.

server = cExperiment.server;
ds = cExperiment.dataset.dataset;

fa=getFileAnnotations(server.session, cTimelapse.fileAnnotation_id);
file=fa.getFile;
fullPath=fullfile(char(file.getPath.getValue), char(file.getName.getValue));
description = ['cTimelapse file uploaded by @babyExperimentOmero.saveTimelapse' datestr(now)];

save(fullPath,'cTimelapse');

% Update the files on Omero:
server.updateFile(ds,fullPath,'dsFiles',fa,'description',description);
end
