function timepointIm = returnSingleTimepointRaw(cTimelapse,timepoint,channel)
% timepointIm = returnSingleTimepointRaw(cTimelapse,timepoint,channel)
%
% returns the raw image stack before corrections (rescaling, shifting,
% background correction,projection) are made. Intended as a background
% function, generally better to use BABYTIMELAPSE.RETURNSINGLETIMEPOINT
%
% Important to note that in the case of numerous 'stack' channels this will
% return a stack, as oppose to BABYTIMELAPSE.RETURNSINGLETIMEPOINT which
% generally returns a projected image (unless the 'stack' 'type' is
% chosen).
%
% INPUTS
% timepoint   :   number indicating the desired timepoint (will access this
%                 element of the cTimepoint array).
% channel     :   number (default 1) indicating which of the channels in
%                 cTimelapse.channelNames to use to identify appropriate
%                 files to load.
% 
% A timepoint and channel are specified, and the channel index is used to
% get a channel string from the cTimelapse.channelNames field. all
% filenames associated with that timepoint (i.e. in cTimepoint.filename
% cell array) that contain the channel string are identified. If there is
% only one this is loaded, whereas if there are more than one they are all
% loaded and put into a stack.
%
% The function also takes a number of other liberties. If the timelapseDir
% is set to 'ignore' it takes the filenames to be absolute (this can be
% done using the method makeFileNamesAbsolute). If this is not the case, it
% constructs the file name from timelapseDir and filename{i}. In doing so
% it takes the liberty of resaving the file name with only the relative
% path (i.e. it throws away anything behind the last / or \).
%
% before loading the existence of the putative file is checked, and if it
% doesn't exist the user is asked to specify a new file location. This is
% done using a special GUI that rechecks whether the file has been found
% every 5 seconds, this is done in case the disk was just temporarily
% disconnected. For cExperiments this is best done using the
% changeRootDirsAll method. This is not done if the timelapseDir is 'ignore'
% and filenames are absolute, since in this case it is very hard to work
% out where the files are (since they could be in different folders).
%
% See also BABYTIMELAPSE.RETURNSINGLETIMEPOINT,
% BABYTIMELAPSE.MAKEFILENAMESABSOLUTE,
% BABYEXPERIMENT.CHANGEROOTDIRALLTIMELAPSES


fileNum = regexp(cTimelapse.cTimepoint(timepoint).filename,cTimelapse.channelNames{channel},'match');

% if the channelName matches numerous parts of the filename, want it to
% pull out only those files for which it has the maximum number of matches.
% This can be important if, for example, you give the name GFP_ and the
% filenames are something like:
%    { gal1-GFP_000002_GFP_.png , gal1-GFP_000002_DIC_.png}
loc = cellfun('length',fileNum);
threshM = max([max(loc) 1]);
loc = loc>=threshM;

ind = find(loc);

if any(loc)
    
    % if the timelapseDir is 'ignore', the filenames are absolute. If not, they
    % are first checked and any file separators removed. This may occur if
    % files are made 'ignore' is used and then switched back.
    if ~strcmp(cTimelapse.timelapseDir,'ignore')
        for i = ind
            filename = cTimelapse.cTimepoint(timepoint).filename{i};
            locSlash = regexp(filename,'[\\|/]','start');
            if ~isempty(locSlash)
                filename = filename(locSlash(end)+1:end);
                cTimelapse.cTimepoint(timepoint).filename{i} = filename;
            end
        end
    end

    ind = find(loc);
    for i = 1:length(ind)
        filename = cTimelapse.cTimepoint(timepoint).filename{ind(i)};
        if strcmp(cTimelapse.timelapseDir,'ignore')
            ffile = filename;
        else
            ffile = fullfile(cTimelapse.timelapseDir,filename);
        end
        
        %this code checks if the file is still where it is suppose
        %to be and resets the timelapseDir by user interface if it
        %is not.
        while ~exist(ffile,'file')
            if ~strcmp(cTimelapse.timelapseDir,'ignore')
                select_new_directory = selectGUIWithWait(5,...
                    sprintf('The file \n\n %s \n\n could not be found. Would you like to select a new location for the images in this timelapse?',ffile),...
                    'select new directory',...
                    'yes',...
                    'retry old location');
                if select_new_directory == 1
                    fprintf('Select the correct folder for: \n %s \n',cTimelapse.timelapseDir);
                    folder = uigetdir(pwd,['Select the correct folder for: ',cTimelapse.timelapseDir]);
                    cTimelapse.timelapseDir = folder;
                    
                    filename = cTimelapse.cTimepoint(timepoint).filename{ind(i)};
                    ffile = fullfile(cTimelapse.timelapseDir,filename);
                end
            else
                fprintf('\n file not found:\n\n %s \n\n',ffile)
                error('your files are absolute (timelapseDir is ignore) and they have also not been found, need to sort out carefully.')
            end
        end
        
        %look for TIF at end of filename and change load method
        %appropriately.
        if ~isempty(regexp(ffile,'TIF$','once'))
            timepointIm(:,:,i) = imread(ffile,'Index',i);
        else
            timepointIm(:,:,i) = imread(ffile);
        end
        %if it is a stack, preallocate. Done in this strange way to
        %preserve data type without making single slices unduly
        %slow.
        if i==1 && isempty(cTimelapse.rawImSize)
            cTimelapse.rawImSize = size(timepointIm);
        end
        if i==1 && sum(loc)>1
            timepointIm(:,:,2:sum(loc)) = 0;
        end
        
        
    end
    
    
    
else
    if ~isempty(cTimelapse.rawImSize)
        timepointIm = zeros(cTimelapse.rawImSize);
    else
        filename = cTimelapse.cTimepoint(timepoint).filename{1};
        if strcmp(cTimelapse.timelapseDir,'ignore')
            ffile = filename;
        else
            ffile = fullfile(cTimelapse.timelapseDir,filename);
        end
        timepointIm = imread(ffile);
        timepointIm(:,:) = 0;
        cTimelapse.rawImSize = size(timepointIm);
    end
    disp('There is no data in this channel at this timepoint');
    return
end




end
