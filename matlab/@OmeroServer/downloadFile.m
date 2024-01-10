function [filePath,fileAnnotation,updated]=downloadFile(this,dataset,fileName,targetFolder)
%downloadFile Downloads files attached to an Omero dataset
%   downloadFile(dataset,fileName) downloads the files specified by the
%   fileName argument to a folder in obj.CachePath that is named 
%   according to the dataset ID.
%
%   - dataset: the integer ID of an Omero dataset, or any Omero object that
%       can have file annotations.
%   - fileName: either a char/cellstr specifying the file names to download
%       from the dataset, an empty array to specify download of all files
%       associated with the dataset, or an array of FileAnnotation objects 
%       specifying the files to download. Note that in the latter case, no 
%       validation is done to ensure that the FileAnnotation objects belong
%       to the specified dataset.
%   - targetFolder (optional): specify this to override the default folder
%       into which the file(s) will be downloaded. Note that in this case,
%       the directory must already exist. Alternatively, specify 'tmp' to 
%       save to a temporary folder within the obj.CachePath folder.
%
%   The function returns path(s) pointing to the files that were downloaded
%   and also the fileAnnotation(s). Note that if more than one file is
%   specified, the outputs will be a cellstr and a FileAnnotation array.
%   
%   If the files already exist, the hash contained in the FileAnnotation
%   will be used to check if the file needs to be updated. If so, the
%   existing file will be replaced.

if nargin<3
    fileName = {};
end

if ischar(fileName), fileName = {fileName}; end
if iscellstr(fileName) || isempty(fileName)
    % We need to obtain file annotations for the specified files
    if isnumeric(dataset) && ~isempty(dataset)
        fileAnnotation = getDatasetFileAnnotations(this.session,dataset);
    else
        otypes = struct2table(getObjectTypes());
        otype = cellfun(@(x) ~isempty(strfind(class(dataset),x)),otypes.class);
        if sum(otype)~=1
            error('Invalid Omero object supplied');
        end
        %THIS APPEARS TO ONLY DOWNLOAD THE LOG AND ACQ FILES FOR AT LEAST
        %SOME DATASETS - NEED TO FIX
        fileAnnotation = getObjectAnnotations(this.session,'file',...
            otypes.name{otype},dataset);
    end
    
    if ~isempty(fileName)
        faNames = arrayfun(@(x) char(x.getFile.getName.getValue),...
            fileAnnotation,'Uni',0);
        if ~all(ismember(fileName,faNames))
            error('The %s file(s) could not be found in the specified dataset',...
                strjoin(strcat('"',fileName(~ismember(fileName,faNames)),'"'),', '));
        end
        faIndices = find(ismember(faNames,fileName));
        faNames = faNames(faIndices);
        [ufaNames,uInds] = unique(faNames);
        if numel(ufaNames)~=numel(faNames)
            warning('Dataset has multiple files with the same name. Downloading first found.');
        end
        fileName = fileAnnotation(faIndices(uInds));
    else
        fileName = fileAnnotation;
    end
end

if ~isa(fileName,'omero.model.FileAnnotationI[]') && ...
        ~isa(fileName,'omero.model.FileAnnotationI')
    error('The "fileName" argument must be a char, cellstr or FileAnnotation array');
end

fileAnnotation = fileName;

% Get the faNames in case a FileAnnotation array was passed in:
faNames = arrayfun(@(x) char(x.getFile.getName.getValue),...
    fileAnnotation,'Uni',0);

if nargin<4, targetFolder = ''; end

if isempty(targetFolder) && ~isnumeric(dataset) && ~isa(dataset,'omero.model.DatasetI')
    targetFolder = 'tmp';
end

if isempty(targetFolder) && ~isempty(dataset) ...
        && (isnumeric(dataset) || isa(dataset,'omero.model.DatasetI'))
    targetFolder = this.downloadDir(dataset);
elseif isempty(targetFolder) || strcmp(targetFolder,'tmp')
	% Files not downloaded from datasets automatically go into the tmp dir
    % unless otherwise specified:
    if ~isdir(this.cachePath)
        error('The "OmeroTemp" dir could not be found.');
    end
    targetFolder = fullfile(this.cachePath,'tmp');
    if ~isdir(targetFolder)
        % The temporary folder does not yet exist, so create it
        mkdir(this.cachePath,'tmp');
    end
else
    if ~isdir(targetFolder)
        %error('The specified targetFolder "%s" does not exist.',targetFolder);
        mkdir(targetFolder);
    end
end

filePath = fullfile(targetFolder,faNames); % cellstr compatible

updated = true(size(filePath));

for f=1:length(fileAnnotation)
    file = filePath{f};
    fa = fileAnnotation(f);
    
    if exist(file,'file')==2
        % If this file already exists, compute a hash to see if we need to
        % update it:
        localHash = DataHash(file,struct('Input','file','Method','SHA-1'));
        omeroHash = char(fa.getFile.getHash.getValue);
        if strcmp(omeroHash,localHash)
            % File has not changed, so skip download:
            updated(f) = false;
            continue
        end
    end
    
%     fileSize=fileAnnotations(n).getFile.getSize.getValue;
%     try
%         if fileSize<67100000%This seems to be the current size limit for downloading files.
%             getFileAnnotationContent(obj.session, fileAnnotations(n), fullPath);
%         else
%             errFilename=fileAnnotations(n).getFile.getName.getValue;
%             errordlg(['File ' char(errFilename) 'is too large to download. Please download from insight client or copy from microscope computer and copy manually into ' obj.CachePath]);
%             disp(['File ' char(errFilename) 'is too large to download. Please download from insight client or copy from microscope computer and copy manually into ' obj.CachePath]);
%             uiwait;
%         end
%     catch err
%         error('File download failed');
%     end
    
    disp(['Downloading ' faNames{f} ' to ' targetFolder]);
    getFileAnnotationContent(this.session, fa, file);
end

if length(fileAnnotation)==1, fileAnnotation = fileAnnotation(1); end
if length(filePath)==1, filePath = filePath{1}; end 
