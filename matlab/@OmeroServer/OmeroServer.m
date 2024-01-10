classdef OmeroServer < handle
    %Class of objects holding information about the contents of the Omero
    %database. Used to coordinate such information between uploading
    %scripts, manual entry to the database and microscope acquisition
    %software
    properties
        cachePath % Full path to folder in which locally cached dataset files are organised
        downloadPath % Full path to temporary folder to which files will be downloaded
    end
    
    properties (Transient)
        client = []
        session = []
        omTimer = []
    end
    
    properties (Transient,Access=private)
        username
        password
        url
        port
    end
    
    methods
        function this = OmeroServer(url,username,password,port)
            
            %First ensure all necessary scripts are are on the Matlab path
            this.preparePath
            
            % Determine the local username and set a userdir accordingly:
            if ismac
                [~, localuser] = system('whoami');
                localuser=localuser(1:end-1); % Remove newline character
                userdir = fullfile('/Users',localuser,'Documents');
            elseif ispc
                userdir = fullfile('C:\Users',getenv('USERNAME'));
            elseif isunix
                userdir = '~/Documents';
            else
                error('unknown system encountered');
            end
            
            % Define default path to organised cache of files stored locally
            this.cachePath = fullfile(userdir,'OmeroTemp');
            
            % Define the default download path to be a sub directory of the
            % cache (note that a call to obj.cleartmp clears this below):
            this.downloadPath = fullfile(this.cachePath,'tmp');
            
            if nargin<4 || isempty(port)
                port = 4064;
            end
            
            this.url = url;
            this.username = username;
            this.password = password;
            this.port = port;
            
            % Clear any temporary directories
            this.cleartmp;
            
            this.login;
        end
        
        function active = sessionActive(this)
            active = false;
            try
                % If this doesn't give an error then the session is active
                this.session.getConfigService;
            catch
                % active is false
                return
            end
            % Session must be active
            active = true;
        end
        
        function logout(this)
            if ~isempty(this.omTimer)
                if isvalid(this.omTimer), stop(this.omTimer); end
                delete(this.omTimer);
            end
            if ~isempty(this.client)
                fprintf('...closing Omero session.\n');
                this.client.closeSession;
            end
        end
        
        function delete(this)
            %Destructor function for OmeroDatabase objects - used to close the
            %session preventing excessive use of server resources
            if ~isempty(this.omTimer) && isvalid(this.omTimer)
                stop(this.omTimer);
                delete(this.omTimer);
            end
            if ~isempty(this.client)
                fprintf('...closing Omero session.\n');
                this.client.closeSession;
            end
        end
    end
    
    methods (Static)
        preparePath
    end
end