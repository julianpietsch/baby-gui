classdef BabyConfig % Deliberately NOT a handle for copying between cTimelapses
    properties (SetAccess={?BabyBrain})
        modelset % the selected model
    end
    
    properties (Dependent)
        modelset_filter % specify to select a particular model
        camera % specify to filter models by camera
        channel % specify to filter models by channel
        zoom % specify to filter models by magnification/lens
        nstacks % specify to filter models by number of Z stacks
        zstacks % specify Z stack indices to use if fewer stacks in model than experiment
        url
    end
    
    properties (SetAccess=private)
        version = 'v0.30.1'
    end
    
    properties (Access=private)
        modelset_filter_val
        camera_val
        channel_val
        zoom_val
        nstacks_val
        zstacks_val
        url_val
    end
    
    properties (Constant,Access=private)
        default_url = 'http://127.0.0.1:5101/'
    end
    
    methods
        function this = BabyConfig(varargin)
            if nargin>0, this = this.with(varargin{:}); end
        end
        
        function this = with(this,varargin)
            %BabyConfig.with Convenience function for adjusting parameters
            %
            %   This function also performs basic type validation for each
            %   of the parameters. Validation against a given experiment is
            %   performed by an instance of BabyBrain.
            
            isindex = @(x) isnumeric(x) && all(x>0 & round(x)==x);
            
            ip = inputParser;
            ip.addParameter('modelset_filter',[],@(x) isempty(x) || ...
                (isrow(x) && ischar(x)));
            ip.addParameter('camera',[],@(x) isempty(x) || ...
                (isrow(x) && ischar(x)));
            ip.addParameter('zoom',[],@(x) isempty(x) || ...
                (isscalar(x) && isindex(x)) || (isrow(x) && ischar(x)));
            ip.addParameter('channel',[],@(x) isempty(x) || ...
                (isrow(x) && ischar(x)) || iscellstr(x) || isstring(x));
            ip.addParameter('nstacks',[],@(x) isempty(x) || ...
                (isvector(x) && isindex(x)));
            ip.addParameter('zstacks',[],@(x) isempty(x) || ...
                (isvector(x) && isindex(x)) || ...
                (iscell(x) && all(@(y) isvector(y) && isindex(y),x)));
            ip.addParameter('url',this.default_url,...
                @(x) isvector(x) && ischar(x) && ...
                ~isempty(regexp(x,'^https?://[-/.:~_a-zA-Z0-9]+/$','once')));
            ip.parse(varargin{:});
            
            % Validate compatibility of channels, nstacks and zstacks
            chname = this.channel_val;
            if ~ismember('channel',ip.UsingDefaults)
                chname = ip.Results.channel;
                if ischar(chname) || isstring(chname)
                    chname = cellstr(chname);
                end
            end
            nch = numel(chname);
            nz = this.nstacks_val;
            if ~ismember('nstacks',ip.UsingDefaults)
                nz = ip.Results.nstacks;
            end
            zs = this.zstacks_val;
            if ~ismember('zstacks',ip.UsingDefaults)
                zs = ip.Results.zstacks;
                if ~iscell(zs) && ~isempty(zs), zs = {zs}; end
            end
            if ~isempty(chname)
                if ~isempty(nz)
                    assert(numel(nz) == nch,...
                        'BabyConfig:InvalidNstacks',...
                        'length of "nstacks" must match number of channels');
                end
                if ~isempty(zs)
                    assert(numel(zs) == nch,...
                        'BabyConfig:InvalidZstacks',...
                        'length of "zstacks" must match number of channels');
                end
            end
            if ~isempty(nz) && ~isempty(zs)
                assert(all(cellfun(@numel,zs)==nz),'BabyConfig:InvalidStacks',...
                    'all "nstacks" must match lengths of "zstacks"');
            end
            
            % Validation complete: update all modified variables
            this.channel_val = chname;
            this.nstacks_val = nz;
            this.zstacks_val = zs;
            
            if ~ismember('modelset_filter',ip.UsingDefaults)
                this.modelset_filter_val = ip.Results.modelset_filter;
            end
            if ~ismember('camera',ip.UsingDefaults)
                this.camera_val = ip.Results.camera;
            end
            if ~ismember('zoom',ip.UsingDefaults)
                zm = ip.Results.zoom;
                if ~isempty(zm) && isnumeric(zm)
                    zm = sprintf('%ux',zm);
                end
                this.zoom_val = zm;
            end
            if ~ismember('url',ip.UsingDefaults)
                this.url_val = ip.Results.url;
            end
        end
        
        function val = get.modelset_filter(this)
            val = this.modelset_filter_val;
        end
        function this = set.modelset_filter(this,val)
            this = this.with('modelset_filter',val);
        end
        
        function val = get.camera(this), val = this.camera_val; end
        function this = set.camera(this,val)
            this = this.with('camera',val);
        end
        
        function val = get.channel(this), val = this.channel_val; end
        function this = set.channel(this,val)
            this = this.with('channel',val);
        end
        
        function val = get.zoom(this), val = this.zoom_val; end
        function this = set.zoom(this,val)
            this = this.with('zoom',val);
        end
        
        function val = get.nstacks(this), val = this.nstacks_val; end
        function this = set.nstacks(this,val)
            this = this.with('nstacks',val);
        end
        
        function val = get.zstacks(this), val = this.zstacks_val; end
        function this = set.zstacks(this,val)
            this = this.with('zstacks',val);
        end
        
        function val = get.url(this)
            if isempty(this.url_val), val = this.default_url;
            else, val = this.url_val; end
        end
        function this = set.url(this,val)
            this = this.with('url',val);
        end
    end
end