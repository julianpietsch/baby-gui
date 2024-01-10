function fileAnnotation = updateFile(this,dataset,file,varargin)
%UPDATEFILE Update existing file(s) or upload new file(s) to Omero
%   UPDATEFILE(OBJ,DATASET,FILE) upload the file(s) specified by FILE (as a
%   char or cellstr of file paths) to the DATASET specified as an Omero 
%   Dataset object or its numeric ID.
%
%   FA = UPDATEFILE(OBJ,...,'Description',DESC) specifies an optional
%   description that will be set for new files only (not currently updated)
%
%   FA = UPDATEFILE(OBJ,...,'dsFiles',FILEANNOTATIONS) can be specified to
%   optionally include pre-obtained FILEANNOTATIONS for the specified
%   dataset to save an additional network call.
%
%   FA = UPDATEFILE(OBJ,...,'Dummy',true) can be specified to ensure that 
%   fileAnnotations exist for the specified files. If the file(s) already
%   exist, the existing annotations are returned; if not present on the
%   server, this function will create 'Dummy' files to produce 
%   fileAnnotations that can later be updated.

% Ensure that file is always a cellstr
if ischar(file), file = {file}; end

ip = inputParser;
ip.addRequired('dataset',@(x) isnumeric(x) || isa(x,'omero.model.DatasetI'));
ip.addRequired('file',@(x) iscellstr(x))
ip.addParameter('dsFiles',[],@(x) isempty(x) || ...
    isa(x,'omero.model.FileAnnotationI[]') || ...
    isa(x,'omero.model.FileAnnotationI'));
ip.addParameter('Description','',@(x) ischar(x) || iscellstr);
ip.addParameter('Dummy',false,@(x) islogical(x) && isscalar(x));
ip.parse(dataset,file,varargin{:});

dsFiles = ip.Results.dsFiles;
description = ip.Results.Description;
dummy = ip.Results.Dummy;

if ischar(description), description = {description}; end
if length(description)~=length(file)
    error('The number of descriptions provided must match the number of files');
end

% Unless this is a Dummy call, then make sure that all files exist:
fileExists = cellfun(@(f) exist(f,'file')==2,file);
if ~dummy && ~all(fileExists)
    error('The file(s) %s are missing. Specify "Dummy" if you want to create placeholders.',...
        strjoin(strcat('"',file(~fileExists),'"'),', '));
end

if ~isa(dsFiles,'omero.model.FileAnnotationI[]') ...
        && ~isa(dsFiles,'omero.model.FileAnnotationI')
    dsFiles = getDatasetFileAnnotations(this.session,dataset);
end

faNames = arrayfun(@(x) char(x.getFile().getName().getValue()),dsFiles,'Uni',0);

fileAnnotation = dsFiles([]); % Start with an empty fileAnnotation list

for f=1:length(file)
    [~,filename,ext] = fileparts(file{f});
    filename = [filename,ext];

    ind = find(strcmp(faNames,filename));

    if isempty(ind)
        % File has not been uploaded before
        if dummy
            if exist(file{f},'file')==2
                % Create a backup copy of this file (being careful not to
                % overwrite any existing files):
                backup_index = 0;
                backup_file = sprintf('%s.%u.tmp',file{f},backup_index);
                while exist(backup_file,'file')==2
                    backup_index = backup_index + 1;
                    backup_file = sprintf('%s.%u.tmp',file{f},backup_index);
                end
                movefile(file{f},backup_file);
                try
                    % Create a dummy file using the original file name:
                    fh = fopen(file{f},'w'); fclose(fh);
                    fileAnnotation(f) = create_on_omero(file{f},description{f});
                catch err
                    movefile(backup_file,file{f});
                    rethrow(err);
                end
                movefile(backup_file,file{f});
            else
                % File does not exist, so temporarily create one:
                fh = fopen(file{f},'w');
                try
                    fclose(fh);
                    fileAnnotation(f) = create_on_omero(file{f},description{f});
                catch err
                    delete(file{f});
                    rethrow(err);
                end
                delete(file{f});
            end
        else
            disp(['Uploading "' filename '" to Omero Database...']);
            fileAnnotation(f) = create_on_omero(file{f},description{f});
        end
        
    else
        if length(ind)>1
            warning('More than one file matching "%s" was found on Omero',filename);
            ind = ind(1);
        end

        %There is a recorded file annotation attached to the dataset representing
        %the cExperiment.
        %Get annotation and check if the file is different to the one on omero:
        fileAnnotation(f) = dsFiles(ind);
        if ~dummy
            omeroHash = char(fileAnnotation(f).getFile().getHash().getValue());
            localHash = DataHash(file{f},struct('Input','file','Method','SHA-1'));
            if ~strcmp(omeroHash,localHash)
                % File has changed, so update it:
                disp(['Updating "' filename '" on Omero Database...']);
                updateFileAnnotation(this.session, fileAnnotation(f), file{f});
            end
        end
    end
end

if length(fileAnnotation)==1
    fileAnnotation = fileAnnotation(1);
end

%% Helper functions
    function fa = create_on_omero(f,d)
        % Create a new fileAnnotation (and OriginalFile):
        fa = writeFileAnnotation(this.session,f,'description',d);
        % Link the fileAnnotation to the current dataset:
        link = omero.model.DatasetAnnotationLinkI;
        link.setParent(dataset);
        link.setChild(fa);
        this.session.getUpdateService().saveAndReturnObject(link);
    end
end