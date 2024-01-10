classdef ImageReaderBioformats < ImageReader
    properties (Dependent, SetAccess=private)
        npos
        posNames
        nchannels
        channels
        imageSize
        pixelSize
        timeInterval
        times
        meta
        instruments
        rawReader
        rawMetaData
    end
    
    properties (Dependent)
        pos
        posName
    end
    
    properties (Transient, Access=private)
        r % bfReader
        m % meta
        pos_val
        meta_val = struct()
        instruments_val
        rawMetaData_val
    end
    
    methods
        function this = ImageReaderBioformats(filename)
            status = bfCheckJavaPath();
            assert(status, ['Missing Bio-Formats library. Either add bioformats_package.jar '...
                'to the static Java path or add it to the Matlab path.']);
            this.r = loci.formats.Memoizer(bfGetReader());
            this.r.setId(filename);
            this.m = this.r.getMetadataStore();
        end
        
        function delete(this)
            if ~isempty(this.r) ...
                    && (isa(this.r,'loci.formats.ChannelSeparator') ...
                    || isa(this.r,'loci.formats.Memoizer')) ...
                    && ~isempty(this.r.getCurrentFile())
                this.r.close();
            end
        end
        
        function refresh(this)
            this.meta_val = struct();
            this.instruments_val = {};
        end
        
        function val = get.meta(this)
            val = this.meta_val;
        end
        
        function val = get.npos(this)
            if ~isfield(this.meta_val,'npos')
                assert(this.r.getSeriesCount()==this.m.getImageCount(),...
                    'Metadata is corrupted');
                this.meta_val.npos = this.r.getSeriesCount();
            end
            val = this.meta_val.npos;
        end
        
        function val = get.posNames(this)
            if ~isfield(this.meta_val,'posnames')
                this.meta_val.posnames = cell(1,this.npos);
                for p=1:this.npos
                    this.meta_val.posnames{p} = char(this.m.getImageName(p-1));
                end
            end
            val = this.meta_val.posnames;
        end
        
        function val = get.pos(this)
            val = this.r.getSeries()+1;
        end
        
        function set.pos(this,val)
            this.r.setSeries(val-1);
        end
        
        function val = get.posName(this)
            val = this.posNames{this.pos};
        end
        
        function set.posName(this,val)
            assert(sum(strcmp(val,this.posNames))==1,...
                'posName not found or not unique');
            this.pos = find(strcmp(this.posNames,val),1);
        end
        
        function val = get.instruments(this)
            if isempty(this.instruments_val)
                this.instruments_val = cell(1,this.m.getInstrumentCount());
                for i=1:numel(this.instruments_val)
                    this.instruments_val{i} = char(this.m.getInstrumentID(i-1));
                end
            end
            val = this.instruments_val;
        end
        
        function val = get.channels(this)
            if ~isfield(this.meta_val,'channels')
                chnames = cell(1,this.npos);
                for p=1:this.npos
                    % getChannelCount(imageIndex)
                    nCh = this.m.getChannelCount(p-1);
                    chnames{p} = cell(1,nCh);
                    for c=1:nCh
                        % getChannelName(imageIndex,channelIndex)
                        chnames{p}{c} = char(this.m.getChannelName(p-1,c-1));
                        if isempty(chnames{p}{c})
                            chnames{p}{c} = sprintf('channel%u',c);
                        end
                    end
                end
                this.meta_val.channels = chnames;
            end
            val = this.meta_val.channels{this.pos};
        end
        
        function val = get.nchannels(this)
            val = this.r.getSizeC();
        end
        
        function val = get.imageSize(this)
            val = [this.r.getSizeX,this.r.getSizeY,this.r.getSizeZ,...
                this.r.getSizeC,this.r.getSizeT];
        end
        
        function val = get.pixelSize(this)
            iid = this.getInstrumentID(this.pos);
            Xpxsz = this.m.getPixelsPhysicalSizeX(iid);
            if ~isempty(Xpxsz)
                Xpxsz = double(Xpxsz.value());
            end
            Ypxsz = this.m.getPixelsPhysicalSizeY(iid);
            if ~isempty(Ypxsz)
                Ypxsz = double(Ypxsz.value());
            end
            val = Xpxsz;
            return
            if isequal(Xpxsz,Ypxsz)
                val = Xpxsz;
            else
                val = [Xpxsz,Ypxsz];
            end
        end
        
        function val = get.rawMetaData(this)
            if isempty(this.rawMetaData_val)
                mraw = this.r.getGlobalMetadata;
                mdata = containers.Map;
                k_iter = mraw.keys;
                for k=1:mraw.size
                    key = k_iter.next;
                    if isempty(key), break; end
                    mdata(key) = mraw.get(key);
                end
                this.rawMetaData_val = mdata;
            end
            val = this.rawMetaData_val;
        end
        
        function val = get.timeInterval(this)
            keyprefs = {'dAvgPeriodDiff','dPeriod','Frame Interval'};
            isms = [1,1,0];
            kInd = find(ismember(keyprefs,this.rawMetaData.keys),1);
            if isempty(kInd)
                warning('interval cannot be determined, assuming 5 mins...');
                val = 5 * 60;
                return
            end
            val = this.rawMetaData(keyprefs{kInd});
            if ischar(val)
                val = str2double(val);
            elseif isa(val,'ome.units.quantity.Time')
                unit = char(val.unit.getSymbol);
                val = double(val.value);
                switch unit
                    case 's'
                    case 'm'
                        val = val*60;
                    case 'h'
                        val = val*60*60;
                    otherwise
                        error('unknown time interval unit');
                end
            elseif ~isnumeric(val)
                error('unrecognised interval specification format');
            end
            if isms(kInd)
                val = val/1000;
            end
        end
        
        function val = get.times(this)
            if ~isfield(this.meta_val,'times')
                this.meta_val.times = cell(1,this.npos);
            end
            p = this.pos;
            if isempty(this.meta_val.times{p})
                ntps = this.r.getSizeT;
                t = NaN(1,ntps);
                for tp=1:ntps
                    iPlane = this.r.getIndex(0,0,tp-1);
                    try
                        t(tp) = double(this.m.getPlaneDeltaT(p-1,iPlane).value());
                    catch
                        t(tp) = this.timeInterval*(tp-1);
                    end
                end
                this.meta_val.times{p} = t;
            end
            val = this.meta_val.times{this.pos};
        end
        
        function val = get.rawReader(this)
            val = this.r;
        end
        
        function val = getTimepoint(this,varargin)
            ip = inputParser;
            ip.addOptional('T',1,@(x) isscalar(x) && isnumeric(x));
            ip.addOptional('Z',[],@(x) isempty(x) || ...
                (isscalar(x) && isnumeric(x)));
            ip.addParameter('C',[],@(x) isempty(x) || ...
                (isscalar(x) && isnumeric(x)) || ...
                (isrow(x) && ischar(x)));
            ip.parse(varargin{:});
            
            T = ip.Results.T;
            Z = ip.Results.Z;
            C = ip.Results.C;
            
            if ischar(C)
                chmatch = contains(this.channels,C,'IgnoreCase',true);
                assert(sum(chmatch)==1,'ambiguous channel specification');
                C = find(chmatch,1);
            end
            
            if isempty(Z)
                Z = 1:this.r.getSizeZ;
            end
            if isempty(C)
                C = 1:this.r.getSizeC;
            end
            
            switch this.r.getBitsPerPixel()
                case 8
                    pxClass = 'uint8';
                case 16
                    pxClass = 'uint16';
                case 12
                    % bfGetPlane returns a uint16 in this case
                    pxClass = 'uint16';
                case 14
                    % bfGetPlane returns a uint16 in this case
                    pxClass = 'uint16';
                otherwise
                    pxClass = 'double';
            end
                
            val = zeros(this.r.getSizeY,this.r.getSizeX,numel(Z),numel(C),pxClass);
            for c=1:numel(C)
                for z=1:numel(Z)
                    iPlane = this.r.getIndex(Z(z)-1,C(c)-1,T-1)+1;
                    val(:,:,z,c) = bfGetPlane(this.r,iPlane);
                end
            end
        end
    end
    
    methods (Access=private)
        function val = getInstrumentID(this,pos)
            if ~isempty(this.instruments)
                val = find(strcmp(this.instruments,...
                    char(this.m.getImageInstrumentRef(pos-1))),1)-1;
            else
                val = 0;
            end
        end
    end
end