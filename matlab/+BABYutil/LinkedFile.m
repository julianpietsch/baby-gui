classdef LinkedFile < handle
    properties (Dependent)
        filename
        relname
        basedir
        exists
        filter
    end
    
    properties (Access=protected,Transient)
        current_relname
        current_basedir
    end
    
    properties (Access=protected)
        relnamelist
        basedirlist
        filter_val = '*.mat'
    end
    
    properties (Constant)
        % Error message types that get raised by this class
        errBadParams = 'BABY:LinkedFile:badParams';
        errMissingFile = 'BABY:LinkedFile:missingFile';
        errUserCancel = 'BABY:LinkedFile:userCancel';
        errConflicting = 'BABY:LinkedFile:fileConflict';
        errNotFound = 'BABY:LinkedFile:targetNotFound';
    end
    
    events
        LinksUpdated % Event triggered when links are modified; used to notify parent objects to re-save
    end
    
    methods
        function this = LinkedFile(relname,basedir,varargin)
            ip = inputParser;
            ip.addRequired('relname',@(x) ischar(x) && isvector(x));
            ip.addRequired('basedir',@(x) ischar(x) && isvector(x));
            ip.addParameter('Filter','*.mat',@(x) ischar(x) && isvector(x));
            ip.parse(relname,basedir,varargin{:});
            
            this.filter_val = ip.Results.Filter;
            
            fullname = fullfile(basedir,relname);
            if ~isdir(basedir)
                error(this.errMissingFile,...
                    'The specified base directory "%s" does not exist.',basedir);
            end
            if exist(fullname,'file')~=2 && ~isdir(fullname)
                error(this.errMissingFile,...
                    'The specified file/directory "%s" does not exist.',fullname);
            end
            this.relnamelist = {relname};
            this.basedirlist = {basedir};
            this.current_relname = relname;
            this.current_basedir = basedir;
        end
        
        function val = get.relname(this)
            if isempty(this.current_relname)
                this.filename; % Runs search for existing basedir/relnames
            end
            val = this.current_relname;
        end
        
        function val = get.basedir(this)
            if isempty(this.current_basedir)
                this.filename; % Runs search for existing basedir/relnames
            end
            val = this.current_basedir;
        end
        
        function updateLink(this,relname,directory)
            %LinkedFile.updateLink Update both name and directory of target
            %   LinkedFile.updateLink(relname,directory) 
            %   - relname: either a char specifying the location of the 
            %   target relative to the specified directory, or another 
            %   LinkedFile object
            %   - directory: the root of the directory tree containing the
            %   linked file. Required unless a LinkedFile object is
            %   provided as the first argument.
            
            if isa(relname,'BABYutil.LinkedFile')
                link = relname;
                relname = link.relname;
                directory = link.basedir;
            else
                if nargin<3 || isempty(directory) || isempty(relname) || ...
                        ~ischar(directory) || ~ischar(relname)
                    error(this.errBadParams,...
                        '"LinkedFile.updateLink" requires both directory and relname as arguments.');
                end
            end
            
            % Arguments must be one-dimensional
            directory = directory(:)';
            relname = relname(:)';
            
            fullname = fullfile(directory,relname);
            if ~isdir(directory)
                error(this.errMissingFile,...
                    'The specified base directory "%s" does not exist.',directory);
            end
            if exist(fullname,'file')~=2 && ~isdir(fullname)
                error(this.errMissingFile,...
                    'The specified file "%s" does not exist.',fullname);
            end
            
            % Set update to true if the lists change:
            update = false;
            if ~ismember(directory,this.basedirlist)
                this.basedirlist{end+1} = directory;
                update = true;
            end
            if ~ismember(relname,this.relnamelist)
                this.relnamelist{end+1} = relname;
                update = true;
            end
            if update
                % The link has been updated, so trigger LinksUpdated event:
                notify(this,'LinksUpdated');
            end
            this.current_relname = relname;
            this.current_basedir = directory;
        end
        
        function updateRelName(this,relname)
            if nargin<2 || isempty(relname) || ~ischar(relname)
                error(this.errBadParams,...
                    '"LinkedFile.updateRelName" requires a character string as argument.');
            end
            relname = relname(:)'; % Needs to be one-dimensional
            
            fileExists = cellfun(@(x) exist(fullfile(x,relname),'file')==2,...
                this.basedirlist);
            fileExists = fileExists | cellfun(@(x) isdir(fullfile(x,relname)),...
                this.basedirlist);
            if ~any(fileExists)
                error(this.errMissingFile,...
                    'None of the registered base directories (%s) contain the linked file "%s"',...
                    strjoin(strcat('"',this.basedirlist,'"'),', '),relname);
            end
            
            if ~ismember(relname,this.relnamelist)
                this.relnamelist{end+1} = relname;
                % The link has been updated, so trigger LinksUpdated event:
                notify(this,'LinksUpdated');
            end
            
            this.current_basedir = this.basedirlist{find(fileExists,1)};
            this.current_relname = relname;
        end
        
        function updateBaseDir(this,directory)
            if nargin<2 || isempty(directory) || ~ischar(directory)
                error(this.errBadParams,...
                    '"LinkedFile.updateBaseDir" requires a character string as argument.');
            end
            directory=directory(:)'; % Needs to be one-dimensional
            
            if ~isdir(directory)
                error(this.errMissingFile,...
                    'The specified directory "%s" does not exist.',directory);
            end
            
            fileExists = this.relNamesExist(directory);
            if ~any(fileExists)
                error(this.errMissingFile,...
                    'None of the registered file names (%s) could be found in "%s".',...
                    strjoin(strcat('"',this.relnamelist,'"'),', '),directory);
            end
            
            if ~ismember(directory,this.basedirlist)
                this.basedirlist{end+1} = directory;
                % The link has been updated, so trigger LinksUpdated event:
                notify(this,'LinksUpdated');
            end
            
            this.current_basedir = directory;
            this.current_relname = this.relnamelist{find(fileExists,1)};
        end
        
        function val = get.filename(this)
            if isempty(this.current_basedir) || isempty(this.current_relname) || ...
                    exist(fullfile(this.current_basedir,this.current_relname),'file')~=2 || ...
                    isdir(fullfile(this.current_basedir,this.current_relname))
                try
                    [dirindex,nameindex] = this.findTarget;
                    this.current_basedir = this.basedirlist{dirindex};
                    this.current_relname = this.relnamelist{nameindex};
                catch err
                    switch err.identifier
                        case this.errNotFound
                            % Could not find file, so prompt user to add one:
                            this.uiAddFile(err.message);
                        case this.errConflicting
                            % Conflicting relnames found, so prompt user to
                            % resolve:
                            this.uiResolveConflicts;
                        otherwise
                            err.rethrow;
                    end
                end
            end
            
            % Unless additional errors were thrown, the basedir and relname
            % should now have been set to an existing target:
            val = fullfile(this.basedir,this.relname);
        end
        
        function val = get.exists(this)
            % Attempt to find target and return false if it fails
            val = true;
            if isempty(this.current_basedir) || isempty(this.current_relname) || ...
                    exist(fullfile(this.current_basedir,this.current_relname),'file')~=2 || ...
                    ~isdir(fullfile(this.current_basedir,this.current_relname))
                try
                    [dirindex,nameindex] = this.findTarget;
                    this.current_basedir = this.basedirlist{dirindex};
                    this.current_relname = this.relnamelist{nameindex};
                catch err
                    if strcmp(err.identifier,this.errNotFound)
                        val = false;
                    else
                        err.rethrow;
                    end
                end
            end
        end
        
        function val = get.filter(this)
            val = this.filter_val;
        end
        
        function [newrelname,newbasedir] = uiAddFile(this,errmsg)
            %uiAddFile Prompt the user to update the base dir or file
            %   This function updates the object, but also returns the new
            %   relname and basedir if the caller wants those values.
            if nargin<2 || isempty(errmsg)
                errmsg = '';
            end
            if ~ischar(errmsg)
                error(this.errBadParams,...
                    'The "msg" argument must be a char');
            end
            
            % Define title and available answers for question dialog
            title = 'Update Linked File';
            answers = {'Update Base Directory','Update File Name','Cancel'};
            
            % Initialise scope of output from uigetfile dialog:
            newfile = '';
            newpath = '';
            
            % Run an infinite loop until the correct file has been found or
            % the user cancels:
            while true
                msg = 'Would you like to:';
                if ~isempty(errmsg)
                    msg = strjoin({errmsg,msg},'\n');
                end
                answer = questdlg(msg,title,answers{1},answers{2},answers{3},answers{1});
                switch answer
                    case answers{1}
                        % Update the base directory
                        newdir = uigetdir('',answers{1});
                        if newdir==0
                            userCancel;
                        else
                            % If there is no error, then directory is
                            % guaranteed to exist. First check if a new
                            % file has also been specified:
                            if ~isempty(newpath)
                                reldir = relativepath(newpath,newdir);
                                if isempty(regexp(reldir,'^\./','once'))
                                    errmsg = 'The specified directory is not a base relative to the new target.';
                                    continue
                                end
                                newrelname = fullfile(reldir,newfile);
                                newbasedir = newdir;
                                this.updateLink(newrelname,newbasedir);
                                return
                            end
                            % Otherwise try setting the base dir for
                            % existing relnames:
                            try
                                newbasedir = newdir;
                                this.updateBaseDir(newbasedir);
                                newrelname = this.relname;
                                return
                            catch err
                                switch err.identifier
                                    case this.errMissingFile
                                        errmsg = err.message;
                                        continue
                                    case this.errConflicting
                                        this.uiResolveConflicts;
                                        continue
                                    otherwise
                                        err.rethrow;
                                end
                            end
                        end
                    case answers{2}
                        % Update the filename
                        if strcmp(this.filter,'/')
                            newfile = uigetdir('',answers{2});
                            if ~isequal(newfile,0)
                                [newpath,newfile,~] = fileparts(newfile);
                            end
                        else
                            [newfile,newpath] = uigetfile(this.filter,answers{2});
                        end
                        if isequal(newfile,0), userCancel;
                        else
                            % If there is no error, then file is
                            % guaranteed to exist. Try to determine a
                            % relative location using existing base dirs:
                            reldirs = cellfun(@(x) relativepath(newpath,x),...
                                this.basedirlist,'UniformOutput',false);
                            isbase = cellfun(@(x) ~isempty(regexp(x,'^\./','once')),reldirs);
                            if ~any(isbase)
                                msg = sprintf(['No base directories exist for that file. ',...
                                    'Do you want to add "%s" as a new base directory?'],newpath);
                                answer = questdlg(msg,'Add new base dir?',...
                                    'Yes','No','Yes');
                                if strcmp(answer,'Yes')
                                    newrelname = newfile;
                                    newbasedir = newpath;
                                    this.updateLink(newrelname,newbasedir);
                                    return
                                end
                                errmsg = sprintf('Choose a base dir for "%s"?',...
                                    fullfile(newpath,newfile));
                                continue
                            end
                            newrelname = fullfile(reldirs{find(isbase,1)},newfile);
                            this.updateRelName(newrelname);
                            newbasedir = this.basedir;
                            return
                        end
                    otherwise
                        userCancel;
                end
            end
            
            function userCancel()
                % User cancelled so raise an error
                if isempty(errmsg)
                    error(this.errUserCancel,'Updating linked file was cancelled.');
                else
                    error(this.errUserCancel,errmsg);
                end
            end
        end
        
        function val = isequal(this,link)
            assert(all(size(this)==size(link)) || length(this)==1 || length(link)==1,...
                this.errBadParams,'Incompatible dimensions for input arrays');
            % Replicate the arrays to the same size if possible:
            if length(link)==1
                link = link(ones(size(this)));
            elseif length(this)==1
                this = this(ones(size(link)));
            end
            val = false(size(this));
            for l=1:length(this)
                if isa(link(l),'BABYutil.LinkedFile')
                    val(l) = any(ismember(this(l).relnamelist,link(l).relnamelist)) ...
                        && any(ismember(this(l).basedirlist,link(l).basedirlist));
                end
            end
        end
    end
    
    methods (Access=private)
        function mask = baseDirsExist(this)
            mask = cellfun(@isdir,this.basedirlist);
        end
        
        function mask = relNamesExist(this,directory)
            mask = cellfun(@(x) exist(fullfile(directory,x),'file')==2,...
                this.relnamelist);
            mask = mask | cellfun(@(x) isdir(fullfile(directory,x)),...
                this.relnamelist);
            if sum(mask)>1
                error(this.errConflicting,...
                    'The specfied base directory "%s" contains conflicting links.',...
                    directory);
            end
        end
        
        function [dirindex,nameindex] = findTarget(this)
            for dirindex=find(this.baseDirsExist)
                nameindex = find(this.relNamesExist(this.basedirlist{dirindex}),1);
                if ~isempty(nameindex)
                    return
                end
            end
            error(this.errNotFound,...
                'The target (%s) could not be found in the registered directories (%s).',...
                strjoin(strcat('"',this.relnamelist,'"'),', '),...
                strjoin(strcat('"',this.basedirlist,'"'),', '));
        end
        
        function uiResolveConflicts(~)
            errordlg('There are conflicts in the LinkedFile object, but the UI hasn''t been written yet...');
        end
    end
end
