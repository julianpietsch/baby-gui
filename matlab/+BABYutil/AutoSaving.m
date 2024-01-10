classdef AutoSaving < handle
    properties (Access=protected)
        savefile % A LinkedFile object that tracks file location
        objtype % The class of the object
    end
    
    properties (Access=private,Transient)
        staticLoad = false;
        linkListener
    end
    
    properties (Constant)
        % Error message types that get raised by this class
        AutoSaveErrors = struct(...
            'BadParams','BABY:AutoSaving:badParams',...
            'MissingDir','BABY:AutoSaving:missingDir',...
            'FileExists','BABY:AutoSaving:fileExists');
    end
    
    methods (Access=protected)
        function this = AutoSaving(objtype,relname,basedir,force)
            %AutoSaving Instantiate a new AutoSaving object
            %   This class is designed to be abstract, so should only be
            %   instantiated by children.
            %   AutoSaving(objtype,relname,basedir,force)
            %   - objtype: class of object being loaded (to be specified by
            %   subclasses)
            %   - relname: save file name relative to "basedir"
            %   - basedir: root of the directory tree within which this
            %   class will be stored
            %   - force: whether or not to force overwriting an existing
            %   file (as specified by basedir/relname).
            
            if nargin<2
                relname = '';
            end
            
            if nargin<3 || isempty(relname) || isempty(basedir)
                fprintf('Please specify a save location for the new %s object.\n',...
                    objtype);
                if isempty(relname)
                    relname = sprintf('my%s.mat',objtype);
                end
                [relname,basedir] = uiputfile('*.mat',...
                    sprintf('Save new %s object as...',objtype),relname);
                if isequal(relname,0) || isequal(basedir,0)
                    error(BABYutil.AutoSaving.AutoSaveErrors.BadParams,...
                        'A save location must be specified to create new %s objects.',...
                        objtype);
                end
                force = true;
            elseif nargin<4 || isempty(force)
                force = false;
            end
            
            % Check that the basedir exists
            if ~isdir(basedir)
                error(BABYutil.AutoSaving.AutoSaveErrors.MissingDir,...
                    'The directory "%s" does not exist.',basedir);
            end
            
            this.objtype = objtype;
            fn = fullfile(basedir,relname);
            
            % Check if the file already exists
            if exist(fn,'file')==2 && ~force
                error(BABYutil.AutoSaving.AutoSaveErrors.FileExists,...
                    'The file "%s" already exists. Use "force" argument to override.',fn);
            end
            
            % First generate a dummy save file to link to
            dummyvar = {}; %#ok<NASGU>
            save(fn,'dummyvar');
            
            % Link to the dummy file and do first autosave
            this.savefile = BABYutil.LinkedFile(relname,basedir);
            this.autosave;
            
            % Finally, add a listener to detect any changes to the list of
            % linked files:
            this.linkListener = addlistener(this.savefile,'LinksUpdated',...
                @this.autosave);
        end
        
        function makeSaveVar(this,varName)
            assignin('caller',varName,this);
        end
    end
    
    methods
        function autosave(this,~,~)
            %AutoSaving.autosave Save this object to its stored location
            %   This function admits two arguments so that it can be
            %   triggered by events.
            fn = this.savefile.filename;
            % Use the file name as the variable name:
            [~,varName,~] = fileparts(fn);
            varName = regexprep(varName,'\W','_');
            this.makeSaveVar(varName);
            save(fn,varName);
        end
    end
    
    methods (Static)
        function obj = loadobj(s)
            if BABYutil.AutoSaving.static_loader
                obj = s;
            else
                error('This object is best loaded using its static "load" function.');
            end
        end
        
        function this = load(relname,basedir,objtype)
            %AutoSaving.LoadAutoSaved Load an AutoSaving object
            %   Children of AutoSaving should overload this function and
            %   specify the objtype
            
            % Use LinkedFile class to handle file existence and updating:
            flink = BABYutil.LinkedFile(relname,basedir);
            
            % Loop until an object of the correct type has been loaded or
            % an error has been thrown
            this = {};
            while ~isa(this,objtype)
                if ~isempty(this)
                    % The loaded object is of the wrong class: use the
                    % LinkedFile method to get user to change target:
                    flink.uiAddFile(sprintf(...
                        'The linked data in "%s" is not of class "%s"',...
                        flink.filename,objtype));
                end
                
                BABYutil.AutoSaving.static_loader(true);
                loaded_struct = load(flink.filename,'-mat');
                BABYutil.AutoSaving.static_loader(false);
                vars = fieldnames(loaded_struct);
                % Take only the first loaded variable in the .mat file
                this = loaded_struct.(vars{1});
                
                % Avoid possible infinite while loops:
                if isempty(this)
                    this = NaN;
                end
            end
            
            % Update the list of basedirs on the loaded object:
            this.savefile.updateLink(flink);
            
            % Finally, add a listener to detect any changes to the list of
            % linked files:
            this.linkListener = addlistener(this.savefile,'LinksUpdated',...
                @this.autosave);
        end
    end

    methods (Static,Access=protected)
        function out = static_loader(in)
            persistent load_state;
            if isempty(load_state)
                load_state=false;
            end
            if nargin
               load_state = in;
            else
               out = load_state;
            end
        end
    end

end