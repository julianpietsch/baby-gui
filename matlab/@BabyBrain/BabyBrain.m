classdef BabyBrain < handle
    properties (Dependent)
        config
        url
        version
        modelset
        zstacks
    end
    
    properties (SetAccess=private)
        channel
        nstacks
        modelsets
        all_modelsets % includes modelsets that do not match cExperiment
        status
        last_status_update
    end
    
    properties (Access=private)
        config_val = BabyConfig
        channel_nzstacks
        ismatching % logical specifying model sets that exactly match cExperiment
        mset_channels
        mset_nstacks
        channelMap
    end
    
    properties (Access=private,Transient)
        last_update
        session_id
        session_modelset
    end
    
    methods
        function this = BabyBrain(varargin)
            if nargin>0, this.config = BabyConfig(varargin{:}); end
            this.set_status('asleep');
        end
        
        function updateServerDetails(this,force)
            if nargin<2 || isempty(force), force = false; end
            if ~isempty(this.last_update) && toc(this.last_update)<30 && ~force
                % Don't update server details unless more than 30s have
                % elapsed since the last update...
                return
            end
                        
            version_info = webread(this.url);
            assert(strcmp(regexprep(this.version,'\.\d+$',''),...
                regexprep(version_info.baby,'\.\d+$','')),...
                'BabyBrain:VersionMismatch',...
                'Configuration version does not match server version');
            
            msets = webread(strcat(this.url,'models?meta=true'));
            this.all_modelsets = msets;
            msets_ch = struct2cell(structfun(@(x) x.channel,msets,'uni',0));
            msets_nz = struct2cell(structfun(@(x) x.z_sections,msets,'uni',0));
            this.modelsets = fieldnames(msets);
            this.mset_channels = msets_ch;
            this.mset_nstacks = msets_nz;
            
            if ~isempty(this.channel_nzstacks)
                % Filter out model sets that are incompatible with
                % experiment:
                valid_mset = false(size(this.modelsets));
                matching_mset = false(size(this.modelsets));
                this.channelMap_defaults;
                channels = this.channel_nzstacks.keys;
                nzstacks = cell2mat(this.channel_nzstacks.values);
                for m=1:numel(msets_ch)
                    mch = msets_ch{m};
                    mnz = msets_nz{m};
                    mchs = strsplit(mch,'+');
                    if isscalar(mnz) && numel(mchs)>1
                        mnz = mnz(ones(size(mchs)));
                        msets_nz{m} = mnz;
                    end
                    chs = cell(size(mchs));
                    hasZ = cell(size(mchs));
                    matchZ = cell(size(mchs));
                    for c=1:numel(mchs)
                        if this.channelMap.isKey(mchs{c})
                            mchs{c} = union(mchs{c},this.channelMap(mchs{c}));
                        end
                        isch = contains(channels,mchs{c},'IgnoreCase',true);
                        chs{c} = channels(isch);
                        nz = nzstacks(isch);
                        hasZ{c} = any(nz>=mnz(c));
                        matchZ{c} = any(nz==mnz(c));
                    end
                    msets_ch{m} = chs;
                    valid_mset(m) = all([hasZ{:}]);
                    matching_mset(m) = all([matchZ{:}]);
                end
                this.ismatching = matching_mset(valid_mset);
                this.modelsets = this.modelsets(valid_mset);
                this.mset_channels = msets_ch(valid_mset);
                this.mset_nstacks = msets_nz(valid_mset);
            else
                this.mset_channels = cellfun(@(x) {{x}},msets_ch,'uni',0);
            end
            
            this.last_update = tic;
        end
        
        function sessions = activeSessions(this)
             sessions = webread(strcat(this.url,'sessions'));
        end
        
        function val = get.config(this)
            val = this.config_val;
        end
        function set.config(this,val)
            assert(isa(val,'BabyConfig'),'BabyBrain:BadConfig',...
                '"config" must be an instance of "BabyConfig"');
            
            if isempty(this.config_val), this.config_val = BabyConfig; end
            
            modelparams = {'modelset_filter','camera','channel','zoom',...
                'nstacks'};
            checkparams = [modelparams,{'zstacks','url','modelset'}];
            ismodified = cellfun(@(x)...
                ~isequal(this.config_val.(x),val.(x)),checkparams);
            modified = checkparams(ismodified);
            ismodified = cell2struct(num2cell(ismodified(:)),checkparams);
            
            % If anything has been modified, then we need to reset session
            if ~isempty(modified)
                this.session_id = [];
            end
            
            % If the url has changed make sure the relevant
            % server-dependent properties will get refreshed
            if ismodified.url
                this.last_update = [];
                % Need a temporary config to change the url
                current_config = this.config_val;
                this.config_val = BabyConfig('url',val.url);
                % Force an update since URL has changed
                this.updateServerDetails(true);
                this.config_val = current_config;
            end
            
            mset = this.config_val.modelset;
            
            if any(ismember(modelparams,modified)) || ...
                    numel(val.zstacks)~=numel(this.config_val.zstacks)
                ch = val.channel;
                if ~isempty(this.channel_nzstacks) && ~isempty(ch)
                    assert(all(ismember(ch,this.channel_nzstacks.keys)),...
                        'BabyBrain:InvalidChannel',...
                        'Experiment does not have the "%s" channel(s)',...
                        strjoin(val.channel,'+'));
                end
                % NB: If channel_nzstacks has been specified, then the available
                % modelsets will also be appropriately filtered by
                % this.updateServerDetails
                
                stackmax = Inf;
                if ~isempty(this.channel_nzstacks)
                    if isempty(ch)
                        stackmax = max(cell2mat(this.channel_nzstacks.values));
                    else
                        stackmax = cellfun(@(c) this.channel_nzstacks(c),ch);
                    end
                end
                
                zs = val.zstacks;
                if ~isempty(zs)
                    if isscalar(stackmax)
                        stackmax = stackmax(ones(size(zs)));
                    end
                    assert(all(cellfun(@(s,m) all(numel(s)<=m),zs,num2cell(stackmax))),...
                        'BabyBrain:InvalidZstacks',...
                        '"zstacks" specifies stack(s) that are out of range');
                end
                
                nz = val.nstacks;
                if ~isempty(nz)
                    assert(all(nz<=stackmax),'BabyBrain:InvalidNstacks',...
                        '"nstacks" exceeds the number of available stacks');
                elseif ~isempty(zs), nz = cellfun(@numel,zs);
                end
                
                this.updateServerDetails;
                
                cam = val.camera; zm = val.zoom;
                msf = val.modelset_filter;
                msets = this.modelsets;
                valid_mset = true(size(msets));
                
                % Apply modelset_filter if specified
                if ~isempty(msf)
                    valid_mset = valid_mset & ...
                        contains(msets,msf,'IgnoreCase',true);
                end
                
                % Apply channel filter if specified
                if ~isempty(ch)
                    valid_mset = valid_mset & ...
                        cellfun(@(m) numel(m)==numel(ch) && ...
                        all(cellfun(@(c) any(ismember(ch,c)),m)), ...
                        this.mset_channels);
                end
                
                if ~isempty(nz)
                    valid_mset = valid_mset & ...
                        cellfun(@(m) all(m==nz),this.mset_nstacks);
                end
                
                params = {};
                
                % Apply param filters if specified
                if ~isempty(cam)
                    params = [params,{cam}];
                    valid_mset = valid_mset & cellfun(@(x) ...
                        contains(this.all_modelsets.(x).camera,cam,'IgnoreCase',true),...
                        msets);
                end
                if ~isempty(zm)
                    params = [params,{zm}];
                    valid_mset = valid_mset & cellfun(@(x) strcmp(...
                        sprintf('%ux',this.all_modelsets.(x).optical_zoom),...
                        zm),msets);
                end
                
                % Check that there is a valid modelset
                if ~any(valid_mset)
                    if ~isempty(ch), params = [params,strjoin(ch,'+')]; end
                    if ~isempty(msf), params = [{msf},params]; end
                    error('BabyBrain:NoMatches',...
                        'No model sets could be found matching: %s',...
                        strjoin(strcat('"',params,'"'),', '));
                end
                
                % Ensure that "zstacks" is specified if there are fewer stacks
                % than in the experiment
                if ~isempty(this.channel_nzstacks) && isempty(zs)
                    valid_mset = valid_mset & this.ismatching;
                    assert(any(valid_mset),'BabyBrain:NeedsZstacks',...
                        '"zstacks" must be specified');
                end
                
                % Pick first matching model set (i.e., order of definition of
                % the model sets determines the default):
                msetind = find(valid_mset,1);
                mset = msets{msetind};
                this.nstacks = this.mset_nstacks{msetind};
                mchs = this.mset_channels{msetind};
                this.channel = cell(size(mchs));
                for c=1:numel(mchs)
                    mch = mchs{c};
                    i = 1; % pick first available by default
                    if isempty(ch) && ~isempty(this.channel_nzstacks)
                        % if known, prefer channel with the closest number
                        % of zstacks
                        nzstacks = cellfun(@(m) this.channel_nzstacks(m),mch);
                        [~,i] = min(abs(nzstacks-this.nstacks(c)));
                    elseif ~isempty(ch)
                        % if channel filter specified then use that
                        i = find(ismember(mch,ch),1);
                    end
                    this.channel{c} = mch{i};
                end
                if isscalar(this.channel), this.channel = this.channel{1}; end
            else
                % No filters were changed
                
                if ismodified.modelset
                    % Only allow silent modelset changes to/from empty
                    % This check is designed to prevent the case where a
                    % new configuration is being copied that silently
                    % changes to a different default model
                    assert(isempty(this.config_val.modelset) || ...
                        isempty(val.modelset),'BabyBrain:NoSilent',...
                        'The "modelset" can only be silently changed to/from empty');
                    mset = val.modelset;
                end
            end
            
            if ~isempty(mset)
                % Ensure updateServerDetails has been called
                this.updateServerDetails;
                % Ensure the selected modelset is available on the server
                assert(ismember(mset,this.modelsets),'BabyBrain:ModelUnavailable',...
                    'The model "%s" is unavailable on the current BABY server',mset);
            end
            
            if isempty(mset)
                % BabyBrain is asleep until it has a valid modelset
                this.set_status('asleep');
            elseif ~isequal(this.config_val.modelset,mset)
                 % BabyBrain is awake if modelset has changed
                this.set_status('awake');
            end
            
            val.modelset = mset;
            this.config_val = val;
        end
        
        function map_channel(this,model_channel,expt_channel)
            iskey = @(x) (ischar(x) && isrow(x)) || (isstring(x) && isscalar(x));
            assert(iskey(model_channel),...
                '"model_channel" must be a single text key');
            model_channel = char(model_channel);
            if iskey(expt_channel), expt_channel = cellstr(expt_channel); end
            if isstring(expt_channel), expt_channel = cellstr(expt_channel); end
            assert(iscell(expt_channel),...
                '"expt_channel" must be a single text key or cellstr of keys');
            assert(all(cellfun(iskey,expt_channel)),...
                '"expt_channel" must be a single text key or cellstr of keys');
            if isempty(this.channelMap)
                this.channelMap = containers.Map;
            end
            if ~this.channelMap.isKey(model_channel)
                this.channelMap(model_channel) = {};
            end
            this.channelMap(model_channel) = union(...
                this.channelMap(model_channel),expt_channel);
        end
        
        function copy_channel_map(this,that)
            this.channelMap = that.channelMap;
        end
        
        function val = get.url(this)
            val = this.config.url;
        end
        
        function val = get.version(this)
            val = this.config.version;
        end
        
        function val = get.modelset(this)
            val = this.config.modelset;
        end
        
        function val = get.zstacks(this)
            val = this.config.zstacks;
        end
        
        function sessionid = get_session(this)
            assert(~isempty(this.modelset),'BabyBrain:NoModelset',...
                'A "modelset" must be specified to obtain a session');
            resp = webread(strcat(this.url,'session/',...
                regexprep(this.modelset,'_','-')));
            sessionid = resp.sessionid;
            this.session_id = sessionid;
        end
        
        function status = queue_image(this,img,pixel_size,assign_mothers,...
                keep_baprobs,refine_outlines,with_volumes)
            if nargin<3, pixel_size = []; end
            if nargin<4 || isempty(assign_mothers)
                assign_mothers = false;
            end
            if nargin<5 || isempty(keep_baprobs)
                keep_baprobs = false;
            end
            if nargin<6 || isempty(refine_outlines)
                refine_outlines = false;
            end
            if nargin<7 || isempty(with_volumes)
                with_volumes = false;
            end
            
            assert(isnumeric(img) && ...
                ismember(class(img),{'double','uint8','uint16'}),...
                'BabyBrain:BadImgType',...
                'Image must be numeric double, uint8 or uint16');
            assert(ndims(img)<=4,'BabyBrain:BadImgShape',...
                'Image must have four or fewer dimensions');
            
            if isempty(this.session_id), this.get_session; end
            
            [N,W,H,Z] = size(img);
            if ~isempty(this.nstacks)
                assert(sum(this.nstacks)==Z,'BabyBrain:BadImgShape',...
                    'Image has incorrect number of stacks for the current modelset');
            end
            
            if isa(img,'double')
                assert(all(img(:)>=0) && all(img(:)<=1),...
                    'BabyBrain:BadImgRange',...
                    'Images of type "double" must lie between 0 and 1');
                img = uint16(round(img*(2^16-1)));
            end
            
            if isa(img,'uint16'), bitdepth = '16';
            elseif isa(img,'uint8'), bitdepth = '8';
            else
                error('BabyBrain:BadImgType',...
                    'Image must be numeric double, uint8 or uint16');
            end
            
            % Construct request to send image for segmentation
            crlf = sprintf('\r\n');
            % NB: the boundary must not be present in image encoding
            bdry = '----BabyFormBoundary';
            options = weboptions('CharacterEncoding','ISO-8859-1',...
                'MediaType',strcat('multipart/form-data; boundary=',bdry));
            bdry_tmplt = @(name) ['--',bdry,crlf,...
                sprintf('Content-disposition: form-data; name="%s"',name),crlf];
            
            data = [bdry_tmplt('dims'),crlf];
            data = [data,savejson('',[N,W,H,Z],'Compact',1),crlf];
            data = [data,bdry_tmplt('bitdepth'),crlf,bitdepth,crlf];
            data = [data,bdry_tmplt('img'),crlf];
            data = [data,char(typecast(img(:)','uint8')),crlf];
            if ~isempty(pixel_size)
                data = [data,bdry_tmplt('pixel_size'),crlf,...
                    num2str(pixel_size),crlf];
            end
            if assign_mothers
                data = [data,bdry_tmplt('assign_mothers'),crlf,'true',crlf];
            end
            if keep_baprobs
                data = [data,bdry_tmplt('return_baprobs'),crlf,'true',crlf];
            end
            if refine_outlines
                data = [data,bdry_tmplt('refine_outlines'),crlf,'true',crlf];
            end
            if with_volumes
                data = [data,bdry_tmplt('with_volumes'),crlf,'true',crlf];
            end
            data = [data,'--',bdry,'--',crlf];
            
            status = webwrite(strcat(this.url,'segment?sessionid=',...
                this.session_id),data,options);
        end
        
        function segout = get_segmentation(this)
            assert(~isempty(this.session_id),'BabyBrain:NoSession',...
                'A session is required before retrieving segmentation tasks');
            segout = webread(strcat(this.url,'segment?sessionid=',...
                this.session_id),weboptions('Timeout',600));
        end
    end
    
    methods (Access={?babyTimelapse,?BabyBrain})
        function set_status(this,val)
            assert(ismember(val,{'asleep','awake','thinking','happy'}));
            this.status = val;
            this.last_status_update = datestr(now);
        end
    end
    
    methods (Access={?babyExperiment,?babyTimelapse})
        function set_channel_nzstacks(this,nzstacks)
            assert(isa(nzstacks,'containers.Map'),'BabyBrain:BadNZstacks',...
                '"nzstacks" must be specified as a "containers.Map"');
            assert(all(cellfun(@(x) isscalar(x) && isnumeric(x) ...
                && x>0 && round(x)==x,nzstacks.values)),'BabyBrain:BadNZstacks',...
                '"nzstacks" for each channel must be scalar indices');
            
            this.channel_nzstacks = nzstacks;
            
            % Ensure that valid model sets get refreshed
            this.last_update = [];
        end
    end
    
    methods (Access=private)
        function channelMap_defaults(this)
            this.map_channel('phase',{'tl'});
        end
    end
end
