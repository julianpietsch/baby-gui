classdef OmeroDataset < handle
    properties
        meta
        use_omero_cache = false
    end
    
    properties (Dependent)
        filenames
        cExperiments
        database
        dataset
        name
        id
        omero_tags
        images
        posNames
        pos
        posNum
        image
        imageSize
        channelNames
        hasZstacks
    end
    
    properties (Transient, Access=private)
        database_val
        dataset_val
        images_val
        imageSizes_val
        channelNames_val
        hasZstacks_val
    end
    
    properties (Access=private)
        datasetId
        pos_val
    end
    
    methods
        function this = OmeroDataset(exptid, varargin)
            ip = inputParser;
            ip.addRequired('exptid',@(x) ...
                (isnumeric(x) && round(x) == x && isscalar(x)) || ...
                (ischar(x) && isvector(x)) || ...
                (isa(x,'omero.model.DatasetI') && isscalar(x)));
            ip.addParameter('meta',cell2struct({},{}),@(x) ...
                isempty(x) || (isstruct(x) && isscalar(x)));
            ip.addParameter('meta_only',true,@(x) isscalar(x) && islogical(x));
            ip.addParameter('use_omero_cache',false,@(x) isscalar(x) && islogical(x));
            ip.parse(exptid, varargin{:});
            
            this.ensure_session;
            
            if isnumeric(exptid)
                this.datasetId = exptid;
            elseif ischar(exptid)
                datasets = this.database.getDsListFromTag(exptid);
                if isempty(datasets)
                    error('No datasets found matching that search term');
                elseif numel(datasets) > 1
                    error('Multiple datasets found matching that search term');
                end
                this.dataset_val = datasets(1);
                this.datasetId = datasets(1).getId().getValue();
            else
                % the argument is a datasetId
                this.datasetId = exptid.getId().getValue();
                this.dataset_val = exptid;
            end
            
            % Ensure that the dataset has been properly loaded
            if ~isa(this.dataset, 'omero.model.DatasetI')
                error('the dataset could not be loaded');
            end
            
            this.use_omero_cache = ip.Results.use_omero_cache;
            this.meta = ip.Results.meta;
            
            if isempty(this.meta)
                this.parseMeta(ip.Results.meta_only);
            end
        end
        
        function parseMeta(this, meta_only)
            if nargin<2, meta_only = true; end
            
            this.ensure_session;
            session = this.database.session;
            
            this.meta = struct();
            this.meta.dsId = this.dataset.getId().getValue();
            this.meta.name = char(this.dataset.getName().getValue());
            
            % Get the fileAnnotations for this dataset. Note that the 
            % following Omero function accepts either numeric IDs or 
            % Dataset objects:
            fileAnnotations = getDatasetFileAnnotations(session, this.dataset);
            faNames = arrayfun(@(x) char(x.getFile().getName().getValue()),...
                fileAnnotations, 'UniformOutput', false);
            
            % Does this dataset have cExperiments
            cExps = ~cellfun(@isempty, regexp(faNames,'^cExperiment_.*\.mat$'));
            cExpNames = faNames(cExps);
            this.meta.cExperiments = regexprep(cExpNames,'^cExperiment_(.*)\.mat$','$1');
            
            if this.use_omero_cache
                downloadArgs = {};
            else
                downloadArgs = {'tmp'};
            end
            
            % Get and parse the log file if it is present
            faLog = fileAnnotations(~cellfun(@isempty,regexp(faNames,'log\.txt$')));
            if isempty(faLog)
                this.meta.logfiledata = 'missing';
            else
                if length(faLog)>1
                    this.meta.multiplelogfiles = true;
                    faLog = faLog(1);
                end
                logfile = this.database.downloadFile(this.dataset,faLog,downloadArgs{:});
                try
                    this.meta.logfiledata = parseLogFile(logfile, meta_only);
                catch err
                    this.meta.logfiledata = 'parse_error';
                    warning(err.identifier,...
                        'error occurred when parsing log file: %s',err.message);
                end
                if ~this.use_omero_cache
                    delete(logfile);
                end
            end
            
            % Get and parse the acq file if it is present
            faAcq = fileAnnotations(~cellfun(@isempty,regexp(faNames,'Acq\.txt$')));
            if isempty(faAcq)
                this.meta.acqfiledata = 'missing';
            else
                if length(faAcq)>1
                    this.meta.multipleacqfiles = true;
                    faAcq = faAcq(1);
                end
                acqfile = this.database.downloadFile(this.dataset,faAcq,downloadArgs{:});
                try
                    this.meta.acqfiledata = parseAcqFile(acqfile);
                catch err
                    this.meta.acqfiledata = 'parse_error';
                    warning(err.identifier,...
                        'error occurred when parsing acq file: %s',err.message);
                end
                if ~this.use_omero_cache
                    delete(acqfile);
                end
            end
        end
        
        function hypercube = getHypercube(this,varargin)
            %GETHYPERCUBE Returns a 5D stack of images as XYZCT
            %
            %   GETHYPERCUBE() returns the entire hypercube for this
            %   position
            %   GETHYPERCUBE(hsize,offset,step) returns a hypercube of
            %   length hsize (XYZCT vector; use NaN to specify default),
            %   starting at offset (XYZCT vector as for hsize), and with 
            %   steps along each dimension of length step (XYZCT vector as 
            %   for hsize)
            
            im = this.image;
            imsize = this.imageSize;
            
            ip = inputParser();
            validvec = @(x) isnumeric(x) && numel(x)==numel(imsize) ...
                && all(round(x)==x|isnan(x)) && all(x>0|isnan(x));
            ip.addOptional('hsize',imsize,validvec);
            ip.addOptional('offset',ones(1,numel(imsize)),validvec);
            ip.addOptional('step',ones(1,numel(imsize)),validvec);
            validrange = @(x) isnumeric(x) && numel(x)<4 && ...
                all(round(x)==x) && all(x>0);
            ip.addParameter('X',[],validrange);
            ip.addParameter('Y',[],validrange);
            ip.addParameter('Z',[],validrange);
            ip.addParameter('C',[],validrange);
            ip.addParameter('T',[],validrange);
            ip.addParameter('pixelstore',{},@(x) isempty(x) || ...
                (iscell(x) && numel(x)==2));
            ip.parse(varargin{:});
            
            hsize = ip.Results.hsize;
            offset = ip.Results.offset;
            step = ip.Results.step;
            
            % Fill in NaNs with defaults
            hsize(isnan(hsize)) = imsize(isnan(hsize));
            offset(isnan(offset)) = 1;
            step(isnan(step)) = 1;
            
            % Determine the extent in the full array considering step size
            % NB: the offset and asize determine the size of the hypercube
            % and step changes the sampling frequency within that
            % hypercube, so the dimensions of the returned data scale with
            % step size:
            asize = (hsize-1).*step+1; % 
            
            setIndexes('X',1);
            setIndexes('Y',2);
            setIndexes('Z',3);
            setIndexes('C',4);
            setIndexes('T',5);
            hsize = floor((asize-1)./step)+1;
            
            j_asize = toJavaList(asize,'java.lang.Integer');
            j_offset = toJavaList(offset-1,'java.lang.Integer');
            j_step = toJavaList(step,'java.lang.Integer');
            
            if isempty(ip.Results.pixelstore)
                this.ensure_session;
                [store, pixels] = this.database.session.getRawPixelsStore(im);
                try
                    hcube = store.getHypercube(j_offset, j_asize, j_step);
                catch err
                    store.close(); % make sure we clean up on error as well
                    rethrow(err);
                end
                store.close();
            else
                [store,pixels] = deal(ip.Results.pixelstore{:});
                hcube = store.getHypercube(j_offset, j_asize, j_step);
            end
            type = char(pixels.getPixelsType().getValue().getValue());
            
            % Convert raw binary output into correct type
            % NB: for some reason getHypercube returns the bytes in reverse
            % order, so need to flip here and unflip below
            if strcmp(type,'float'), type = 'single'; end
            hcube = typecast(flip(hcube), type);
            
            hypercube = permute(reshape(flip(hcube), hsize), [2,1,3:5]);
            
            function setIndexes(dn,di)
                dv = ip.Results.(dn);
                if ~isempty(dv)
                    offset(di) = dv(1);
                    if numel(dv)==1, asize(di) = 1; return; end
                    if numel(dv)>1, asize(di) = dv(2)-dv(1)+1; end
                    if numel(dv)>2, step(di) = dv(3); end
                end
            end
        end
        
        function imstruct = getSampleImages(this,nsamples,varargin)
            ip = inputParser();
            ip.addRequired('nsamples',@(x) isscalar(x) && isnumeric(x) ...
                && round(x)==x && x>0);
            validrange = @(x) isnumeric(x) && numel(x)<4 && ...
                all(round(x)==x) && all(x>0);
            ip.addParameter('X',[],validrange);
            ip.addParameter('Y',[],validrange);
            ip.addParameter('Z',[],validrange);
            chnames = this.channelNames;
            ip.addParameter('Channels',chnames,@(x) iscellstr(x) ...
                && all(ismember(x,chnames)));
            ip.parse(nsamples,varargin{:});
            
            usechannel = ismember(this.channelNames,ip.Results.Channels);
            
            % Calculate step size in time dimension
            ntps = this.imageSize(5);
            tstep = min(floor((ntps-1)/(nsamples-1)),ntps+1);
            
            imstruct = struct();
            [store,pixels] = this.database.session.getRawPixelsStore(this.image);
            try
                for zstacks=0:1
                    % treat channels with or without z stacks differently
                    if zstacks
                        indsToDo = find(this.hasZstacks & usechannel);
                        zrange = ip.Results.Z;
                    else
                        indsToDo = find(~this.hasZstacks & usechannel);
                        zrange = 1;
                    end
                    
                    while ~isempty(indsToDo)
                        if length(indsToDo)==1
                            chInds = indsToDo;
                            chStep = 1;
                        else
                            chMax = max(indsToDo); chMin = min(indsToDo);
                            testChSteps = 1:(chMax-chMin);
                            testInds = cell(size(testChSteps));
                            for c=testChSteps
                                tInds = chMin:c:chMax;
                                badInd = find([~ismember(tInds,indsToDo),true],1);
                                testInds{c} = tInds(1:(badInd-1));
                            end
                            [~,chStep] = max(cellfun(@length,testInds));
                            chInds = testInds{chStep};
                        end
                        
                        hcube = this.getHypercube('X',ip.Results.X,...
                            'Y',ip.Results.Y,'Z',zrange,...
                            'C',[min(chInds),max(chInds),chStep],...
                            'T',[1,ntps,tstep],'pixelstore',{store,pixels});
                        hcube = permute(hcube,[1:3,5,4]);
                        for c=1:length(chInds)
                            imstruct.(this.channelNames{chInds(c)}) = ...
                                hcube(:,:,:,:,c);
                        end
                        
                        indsToDo = setdiff(indsToDo,chInds);
                    end
                end
            catch err
                store.close(); % make sure we clean up on error as well
                rethrow(err);
            end
            store.close();
        end
        
        function val = get.database(this)
            if isempty(this.database_val)
                this.database_val = Nursery.omerodb;
            end
            val = this.database_val;
        end
        
        function val = get.dataset(this)
            if isempty(this.dataset_val)
                if ~isnumeric(this.datasetId) || ~isscalar(this.datasetId)
                    error('value of datasetId has been corrupted');
                end
                this.ensure_session;
                this.dataset_val = this.database.session.getDatasets(this.datasetId);
            end
            val = this.dataset_val;
        end
        
        function val = get.name(this)
            val = char(this.dataset.getName().getValue());
        end
        
        function val = get.id(this)
            val = '';
            acqFile = OmeroDataset.get_from_tree(this.meta,...
                'logfiledata','acqFile');
            if ~isempty(acqFile) && ischar(acqFile) && isvector(acqFile)
                val = babyExperiment.parseAcqFileIntoID(acqFile);
            end
        end
        
        function val = get.filenames(this)
            fileAnnotations = getDatasetFileAnnotations(...
                this.database.session, this.dataset);
            val = arrayfun(@(x) char(x.getFile().getName().getValue()),...
                fileAnnotations, 'UniformOutput', false);
        end
        
        function val = get.cExperiments(this)
            val = {};
            
            try
                % Attempt to get the most up-to-date list
                this.ensure_session;
                % Get the fileAnnotations for this dataset. Note that the 
                % following Omero function accepts either numeric IDs or 
                % Dataset objects:
                fileAnnotations = getDatasetFileAnnotations(...
                    this.database.session, this.dataset);
                faNames = arrayfun(@(x) char(x.getFile().getName().getValue()),...
                    fileAnnotations, 'UniformOutput', false);

                % Does this dataset have cExperiments
                cExps = ~cellfun(@isempty, regexp(faNames,'^cExperiment_.*\.mat$'));
                cExpNames = faNames(cExps);
                this.meta.cExperiments = regexprep(cExpNames,'^cExperiment_(.*)\.mat$','$1');
            catch err
                if ~strcmp(err.identifier,'put:noconnection:error:here')
                    rethrow(err);
                end
            end
            
            if ~isempty(this.meta) && isstruct(this.meta) ...
                    && isfield(this.meta,'cExperiments')
                val = this.meta.cExperiments;
            end
        end
        
        function val = get.omero_tags(this)
            this.ensure_session;
            tags = getDatasetTagAnnotations(this.database.session, this.dataset);
            val = cell2struct([...
                tagsiter(@(x) x.getId(),false),...
                tagsiter(@(x) x.getDescription(),true),...
                tagsiter(@(x) x.getTextValue(),true)],...
                {'id','description','name'},2);
            
            function proplist = tagsiter(fn,ischar)
                proplist = cell(numel(tags),1);
                for t=1:numel(tags)
                    prop = fn(tags(t));
                    if ~isempty(prop), prop = prop.getValue(); end
                    if ischar, prop = char(prop); end
                    proplist{t} = prop;
                end
            end
        end
        
        function val = get.images(this)
            if isempty(this.images_val)
                if ~isnumeric(this.datasetId) || ~isscalar(this.datasetId)
                    error('value of datasetId has been corrupted');
                end
                this.ensure_session;
                imgs = this.database.session.getImages('dataset', this.datasetId);
                img_names = arrayfun(@(x) char(x.getName().getValue()), ...
                    imgs, 'UniformOutput', false);
                % Make all img names valid field names for a struct:
                img_names = regexprep(img_names,'^_','x_','once');
                img_names = regexprep(img_names,'^([0-9])','x_$1','once');
                this.images_val = cell2struct(...
                    arrayfun(@(x) x, imgs, 'uni', 0), img_names);
            end
            val = this.images_val;
        end
        
        function val = get.posNames(this)
            val = {};
            if ~isempty(this.meta)
                if isfield(this.meta, 'logfiledata')
                    logdata = this.meta.logfiledata;
                    if ~isempty(logdata) && isstruct(logdata) ...
                            && isfield(logdata, 'logPosNames')
                        lpns = logdata.logPosNames;
                        if iscell(lpns)
                            lpns = lpns(~cellfun(@isempty, lpns));
                            if iscellstr(lpns)
                                val = lpns;
                            end
                        end
                    end
                end
                
                if isempty(val) && isfield(this.meta, 'acqfiledata')
                    acqdata = this.meta.acqfiledata;
                    if ~isempty(acqdata) && isstruct(acqdata) ...
                            && isfield(acqdata, 'positions') ...
                            && isstruct(acqdata.positions) ...
                            && isfield(acqdata.positions, 'name') ...
                            && iscellstr(acqdata.positions.name)
                        val = acqdata.positions.name;
                    end
                end
            end
            
            % Make all pos names valid field names for a struct:
            val = regexprep(val,'^_','x_','once');
            val = regexprep(val,'^([0-9])','x_$1','once');
            
            imNames = fieldnames(this.images);
            if isempty(val) || any(~ismember(imNames, val))
                val = imNames;
            end
            val = intersect(val, imNames);
        end
        
        function val = get.pos(this)
            pnames = this.posNames;
            if isempty(this.pos_val) || ~ismember(this.pos_val,pnames)
                if isempty(pnames), this.pos_val = '';
                else, this.pos_val = pnames{1}; end
            end
            val = this.pos_val;
        end
        function set.pos(this,val)
            assert(ischar(val) && isvector(val) && ismember(val,this.posNames),...
                '"pos" must be one of the "posNames"');
            this.pos_val = val;
        end
        
        function val = get.posNum(this)
            if isempty(this.pos)
                val = [];
            else
                val = find(strcmp(this.posNames,this.pos),1);
            end
        end
        function set.posNum(this,val)
            pnames = this.posNames;
            assert(isnumeric(val) && isscalar(val) && round(val)==val ...
                && val>0 && val<=length(pnames),...
                '"posNum" must be an index into "posNames"');
            this.pos_val = pnames{val};
        end 
        
        function val = get.image(this)
            val = this.images.(this.pos);
        end    
        
        function val = get.imageSize(this)
            if isempty(this.imageSizes_val)
                this.imageSizes_val = structfun(@getSize, this.images, ...
                    'UniformOutput', false);
            end
            val = this.imageSizes_val.(this.pos);
            
            function xyzct = getSize(img)
                pixels = img.getPrimaryPixels();
                xyzct = [pixels.getSizeX().getValue(), pixels.getSizeY().getValue(), ...
                    pixels.getSizeZ().getValue(), pixels.getSizeC().getValue(), ...
                    pixels.getSizeT().getValue()];
            end
        end
        
        function val = get.channelNames(this)
            if ~isempty(this.channelNames_val)
                val = this.channelNames_val.(this.pos);
                return
            end
            
            this.ensure_session;
            session = this.database.session;
            vals = struct();
            this.imageSize; % ensure that this.imageSizes_val is initialised
            pnames = this.posNames;
            for p=1:length(pnames)
                pname = pnames{p};
                channels = loadChannels(session,this.images.(pname));
                assert(length(channels)==this.imageSizes_val.(pname)(4),...
                    'number of channels does not match image size');
                vals.(pname) = cell(length(channels),1);
                for c=1:length(channels)
                    % First attempt to get channel name directly from
                    % Omero:
                    channel = channels(c).getLogicalChannel();
                    if ~isempty(channel.getName())
                        vals.(pname){c} = char(channel.getName().getValue());
                        continue
                    end
                    
                    % Otherwise try to get from log file
                    etimes = OmeroDataset.get_from_tree(this.meta,...
                        'logfiledata','logExposureTimes');
                    if isscalar(etimes) && isstruct(etimes)
                        lchnames = fieldnames(etimes);
                        lchnames = lchnames(structfun(@(c) c(p)>0,etimes));
                        if c<=length(lchnames)
                            vals.(pname){c} = lchnames{c};
                            continue
                        end
                    end
                    
                    % Lastly attempt to get from acq file
                    achannels = OmeroDataset.get_from_tree(this.meta,...
                        'acqfiledata','channels');
                    if OmeroDataset.hascol(achannels,'names')
                        achnames = achannels.names;
                        aetimes = OmeroDataset.get_from_tree(this.meta,...
                            'acqfiledata','positions');
                        achnames = achnames(cellfun(@(c) ...
                            ~OmeroDataset.hascol(aetimes,c) ...
                            || aetimes.(c)(p)>0,achnames));
                        if c<=length(achnames)
                            vals.(pname){c} = achnames{c};
                            continue
                        end
                    end
                end
            end
            % Cache for reuse
            this.channelNames_val = vals;
            val = vals.(this.pos);
        end
        
        function val = get.hasZstacks(this)
            if isempty(this.hasZstacks_val)
                this.hasZstacks_val = struct();
            end
            if isfield(this.hasZstacks_val,this.pos)
                % Value has been cached
                val = this.hasZstacks_val.(this.pos);
                return
            end
            
            % First attempt to determine if this channel uses z stacks from
            % the channels table in the acq file
            chtbl = OmeroDataset.get_from_tree(this.meta,...
                'acqfiledata','channels');
            chnames = this.channelNames;
            if ~isempty(chtbl) && OmeroDataset.hascol(chtbl,'zsect') ...
                    && OmeroDataset.hascol(chtbl,'names') ...
                    && all(ismember(chnames,chtbl.names))
                val = logical(chtbl.zsect(ismember(chtbl.names,chnames)));
                % Cache for reuse
                this.hasZstacks_val.(this.pos) = val;
                return
            end
            
            % NB: z stack info is not present in the log file
            
            % Infer stacking based on whether z>1 slices sum to zero or not
            imsize = this.imageSize;
            if imsize(3)<=1
                val = false(length(chnames),1);
            else
                Nx2 = NaN(1,2);
                tpstep = floor((imsize(5)-1)/2);
                hcube = this.getHypercube([Nx2,imsize(3)-1,NaN,3],...
                    [Nx2,2,Nx2],[NaN(1,4),tpstep]);
                val = sum(reshape(permute(hcube,[4,1:3,5]),imsize(4),[]),2)>0;
            end
            
            % Cache for reuse
            this.hasZstacks_val.(this.pos) = val;
        end
        
        function ensure_session(this)
            try
                this.database.session.keepAlive([]);
            catch
                this.database.login;
            end
        end
    end
    
    methods (Static)
        function val = get_from_tree(s, varargin)
            val = [];
            lvl = s;
            for l=1:length(varargin)
                if ~isstruct(lvl) || ~isscalar(lvl), return; end
                index = varargin{l};
                if ~isfield(lvl, index), return; end
                lvl = lvl.(index);
            end
            val = lvl;
        end
        
        function val = hascol(tbl,colname)
            val = istable(tbl) && ...
                ismember(colname,tbl.Properties.VariableNames);
        end
    end
end
