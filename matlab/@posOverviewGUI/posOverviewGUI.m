classdef posOverviewGUI < MemAware
    properties
        cTimelapse % A cTimelapse
        fig
        imgax
        slider
        alphascale = 0.8;
        darkscale = 0.7;
        showcells = true;
        fillcells = false;
    end
    
    properties (Dependent)
        ctp
        channel
    end
    
    properties (Access={?MemAware,?posOverviewGUI})
        loadedtps
        imcache
        trackColours
        
        ntimepoints
        ctp_val = 1
        channel_val
    end
    
    properties (Constant)
        MaxTrackColours = 100;
        maxcachetps = 100;
    end
    
    methods
        function this = posOverviewGUI(cTimelapse,channel)
            this.cTimelapse = cTimelapse;
            
            if nargin<2 || isempty(channel)
                chnames = cTimelapse.channelNames;
                isslice = ~cellfun(@isempty,regexp(chnames,'^\w+_\d*$'));
                channel = find(isslice,1);
                if isempty(channel), channel = 1; end
            end
            
            this.channel_val = channel;
            
            figtitle = 'Segmenting...';
            if isstruct(cTimelapse.metadata) && isfield(cTimelapse.metadata,'date')
                figtitle = strcat(figtitle,' - ',cTimelapse.metadata.date);
            end
            if isstruct(cTimelapse.metadata) && isfield(cTimelapse.metadata,'experiment')
                figtitle = strcat(figtitle,' - ',cTimelapse.metadata.experiment);
            end
            
            this.fig = figure('MenuBar','none','Name',figtitle,'Visible','on');
            this.fig.UserData = this;
            
            scrsz = get(0,'ScreenSize');
            scrW = scrsz(3); scrH = scrsz(4);
            imsize = this.cTimelapse.imSize;
            sp = 2;
            
            aspect = imsize(2)/imsize(1);
            figheight = 0.8*scrH;
            figwidth = aspect*figheight;

            % Update the dimensions and position of the figure
            set(this.fig,'Position',...
                [(scrW-figwidth)/2,(scrH-figheight)/2,figwidth,figheight]);

            this.imgax = axes('Parent',this.fig,'XTick',[],'YTick',[],...
                'XLim',[0.5,this.cTimelapse.imSize(2)+0.5],...
                'YLim',[0.5,this.cTimelapse.imSize(1)+0.5],'YDir','reverse',...
                'Units','pixels','Position',[sp,sp,figwidth-2*sp,figheight-2*sp]);
            
            % Show white image to begin with
            image('CData',ones([imsize,3]),'Parent',this.imgax);
            
            % Set up cache
            this.ntimepoints = length(cTimelapse.cTimepoint);
            this.loadedtps = false(this.ntimepoints,1);
            this.imcache = zeros([imsize,this.ntimepoints],'uint8');
            
            % Set up a palette of colours for cells
            this.trackColours = jet(this.MaxTrackColours);
            this.trackColours = this.trackColours(randperm(this.MaxTrackColours),:);
        end
        
        function delete(this)
            if ~isempty(this.fig) && isvalid(this.fig)
                delete(this.fig);
            end
        end
        
        function addTP(this,tp,posim)
            assert(isscalar(tp) && isnumeric(tp),...
                '"tp" must be a single time point index');
            if this.loadedtps(tp), return; end
            
            if nargin<3 || isempty(posim)
                posim = this.cTimelapse.returnSingleTimepoint(...
                    tp,this.channel,'max');
            end
            
            assert(size(posim,1)==size(this.imcache,1) && ...
                size(posim,2)==size(this.imcache,2),...
                'Dimensions of posim need to match those of the cTimelapse');
            
            % Ensure that the posim is a double
            posim = double(posim);
            
            posim_range = num2cell(quantile(posim(:),[0.0001,0.9999]));
            [posim_min,posim_max] = deal(posim_range{:});
            
            posim = (double(posim) - double(posim_min))/double(posim_max-posim_min);
            posim(posim(:)<0) = 0;
            posim(posim(:)>1) = 1;
            
            this.imcache(:,:,tp) = uint8(round(255*posim));
            this.loadedtps(tp) = true;
        end
        
        refreshPosImage(this)
        
        function val = get.ctp(this)
            val = this.ctp_val;
        end
        function set.ctp(this,val)
            this.ctp_val = val;
            this.refreshPosImage;
        end
        
        function val = get.channel(this)
            val = this.channel_val;
        end
        function set.channel(this,val)
            this.channel_val = val;
            this.imcache = zeros([this.cTimelapse.imSize,this.ntimepoints],'uint8');
            this.loadedtps = false(this.ntimepoints,1);
            this.refreshPosImage;
        end
    end
end