classdef Nursery < BABYutil.AutoSaving
    %Nursery Manages setting and retrieval of preferences
    %   This class provides global access (via its static functions) to
    %   preferences for BABY and its GUIs
    
    properties (Dependent)
        basedir
        filename
        localOmeroCache % Table specifying experiments to retain downloads for
    end
    
    properties (Transient,Access=private)
        couch_listener % Event handle to listen for changes to CouchServer
        omerodb_cache % OmeroDatabase object
    end
    
    properties(Access=private)
        couch_server
        username_val
        omero_url
        omero_username
        omero_password
        omero_port
        Omeropwd_val
        couchUsername_val
        couchPwd_val

        localOmeroCache_val = struct2table(struct(...
            'dirname',{{}},'cExperiment',[],'rawimages',[],...
            'previewimages',[]))
    end
    
    properties (Constant,Access=private)
        myclass = 'Nursery';
    end
    
    properties (Constant)
        % Error message types that get raised by this class
        Errors = struct(...
            'BadParams','BABY:Nursery:badParams',...
            'NotActive','BABY:Nursery:notActive',...
            'WrongType','BABY:Nursery:wrongType',...
            'DoesNotExist','BABY:Nursery:DoesNotExist',...
            'AlreadyExists','BABY:Nursery:alreadyExists');
        defaultOmeroCacheRow = struct2table(struct(...
            'dirname',{''},'cExperiment',true,'rawimages',false,...
            'previewimages',false),'AsArray',true)
    end
    
    methods
        function this = Nursery(varargin)
            % Initialise the parent AutoSaving class
            this@BABYutil.AutoSaving(Nursery.myclass,...
                varargin{:});
                
            % Create a new Nursery object
            this.autosave;
            
            % Set this object as the new active Nursery:
            Nursery.active(this);
        end
        
        function val = get.basedir(this)
            % The base directory for the Nursery should be the
            % directory in which this file has been saved
            val = fileparts(this.savefile.filename);
        end
        
        function val = get.filename(this)
            % The location where these preferences are saved
            val = this.savefile.filename;
        end
        
        function val = get.localOmeroCache(this)
            val = this.localOmeroCache_val;
        end
        
        function set.localOmeroCache(this,val)
            if ~istable(val) || ~all(strcmp(val.Properties.VariableNames,...
                    this.defaultOmeroCacheRow.Properties.VariableNames))
                error(['The "localOmeroCache" must be a table with ',...
                    'variable names matching Nursery.defaultOmeroCacheRow']);
            end
            this.localOmeroCache_val = val;
            this.autosave;
        end
        
        function clearcaches(this)
            if ~isempty(this.slib_cache)
                cellfun(@delete,this.slib_cache.values);
                this.slib_cache.remove(this.slib_cache.keys);
            end
            if ~isempty(this.mlib_cache)
                cellfun(@delete,this.mlib_cache.values);
                this.mlib_cache.remove(this.mlib_cache.keys);
            end
        end
    end
    
    methods (Static)
        function this = load(relname,basedir)
            this = load@BABYutil.AutoSaving(...
                relname,basedir,Nursery.myclass);
        end
        
        function this = loadobj(obj)
            % Set this object to be the active one
            Nursery.active(obj);
            this = obj;
        end
        
        function val = active(obj)
            persistent active_preferences;
            if nargin<1 || isempty(obj)
                if isempty(active_preferences)
                    % First check to see if we can find a Nursery
                    % object in the global namespace
                    PrefsVars = evalin('base','whos');
                    PrefsVars = PrefsVars(strcmp({PrefsVars.class},'Nursery'));
                    if(length(PrefsVars)==1)
                        active_preferences = evalin('base',PrefsVars.name);
                    else
                        error(Nursery.Errors.NotActive,...
                            'Saved preferences have not been loaded yet.');
                    end
                end
                val = active_preferences;
            else
                if ~isa(obj,'Nursery')
                    error(Nursery.Errors.WrongType,...
                        'Can only set active preferences to "Nursery" objects.');
                end
                active_preferences = obj;
            end
        end
        
        function val = username()
            this = Nursery.active;
            if isempty(this.username_val)
                error('Username has not been set');
            end
            val = this.username_val;
        end
        
        function val = Omeropwd()
            this = Nursery.active;
            if isempty(this.Omeropwd_val)
                val=passcode;
                this.setOmeroPwd(val);
            end
            val = this.Omeropwd_val;
        end     
        
        function val = couchPwd()
            this = Nursery.active;
            if isempty(this.couchPwd_val)
                val=passcode;
                this.setCouchPwd(val);
            end
            val = this.couchPwd_val;
        end 
               
         function val = couchUsername()
            this = Nursery.active;
            if isempty(this.couchUsername_val)
                prompt = 'couchDb username';
                dlgtitle = 'Enter username for couchDb';
                dims = [1 55];
                definput={getenv('username')};
                answer = inputdlg(prompt,dlgtitle,dims,definput);
                val=answer{:};
                this.setCouchUsername(val);
            end
            val = this.couchUsername_val;
        end 
        
        function setUsername(val)
            this = Nursery.active;
            if isempty(val)
                error('Username cannot be empty');
            end
            if ~isempty(this.username_val)
                warning('Changing username from "%s" to "%s"',...
                    this.username_val,val);
            end
            this.username_val = val;
            this.autosave;
        end
        
        function setOmeroPwd(val)
            this = Nursery.active;
            if isempty(val)
                error('Omero details cannot be empty');
            end
            if ~isempty(this.username_val)
                warning('Changing Omero from "%s" to "%s"',...
                    this.Omeropwd_val,val);
            end
            this.Omeropwd_val = val;
            this.autosave;
        end
        
        function setCouchPwd(val)
            this = Nursery.active;
            if isempty(val)
                error('couchDb password cannot be empty');
            end
            if ~isempty(this.couchPwd_val)
                warning('Changing couchDb password from "%s" to "%s"',...
                    this.couchPwd_val,val);
            end
            this.couchPwd_val = val;
            this.autosave;
        end
        
        function setCouchUsername(val)
            this = Nursery.active;
            if isempty(val)
                error('couchDb user name cannot be empty');
            end
            if ~isempty(this.couchUsername_val)
                warning('Changing couchDb username from "%s" to "%s"',...
                    this.couchUsername_val,val);
            end
            this.couchUsername_val = val;
            this.autosave;
        end
        
        
        function val = csrv()
            this = Nursery.active;
            if isempty(this.couch_server)
                this.couch_server = CouchServer;
                this.autosave;
            end
            if isempty(this.couch_listener)
                this.couch_listener = ...
                    addlistener(this.couch_server,'Updated',@this.autosave);
            end
            val = this.couch_server;
        end
        
        function setOmeroConfig(url,username,password,port)
            this = Nursery.active;
            this.omero_url = url;
            this.omero_username = username;
            this.omero_password = password;
            this.omero_port = port;
            this.autosave;
        end
        
        function val = omerodb()
            this = Nursery.active;
            if isempty(this.omerodb_cache) || ~isvalid(this.omerodb_cache)
                this.omerodb_cache = OmeroServer(this.omero_url,...
                    this.omero_username,this.omero_password,this.omero_port);
            end
            val = this.omerodb_cache;
        end
    end
    
    methods (Static,Access=private)
        function confirmReplaceDlg(name)
            msg = sprintf(...
                'Are you sure you want to replace the %s?',name);
            answer = questdlg(msg,'Replace?','Yes','No','No');
            if ~strcmp(answer,'Yes')
                error(MetaInfo.Errors.AlreadyExists,...
                    'The %s has already been set in Nursery.',name);
            end
        end
    end
end