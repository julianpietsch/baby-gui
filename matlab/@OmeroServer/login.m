function login(this)
    %Populates the Session and Client fields of an OmeroDatabase object by
    %logging in to the database
    sesh=true;
    try
        %If this line doesn't give an error then there is an already active session - don't do anything       
        this.session.getConfigService;    
    catch
        % Close any sessions that may have expired
        if ~isempty(this.omTimer)
            if isvalid(this.omTimer), stop(this.omTimer); end
            delete(this.omTimer);
        end
        if ~isempty(this.client), this.client.closeSession(); end
        
        if exist('omero.client','class') == 0
            this.client = loadOmero(this.url,this.port);
        else
            this.client = omero.client(this.url,this.port);
        end
        
        % Obtain a fresh session and ping to keep alive
        fprintf('Starting Omero session...\n');
        this.session= this.client.createSession(this.username, this.password);
        this.omTimer = omeroKeepAlive(this.client);
        sesh=false;
    end
        
    if sesh
        fprintf('Omero session is already active.\n');
    end
end