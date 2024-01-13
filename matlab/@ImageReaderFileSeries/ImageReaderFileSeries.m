classdef ImageReaderFileSeries < ImageReader
    %ImageReaderFileSeries Read time lapses stored as file series
    %
    %   ImageReaderFileSeries(...,'pixelSize',PXSZ) specifies a pixel size
    %   of PXSZ microns for the data set.
    %
    %   ImageReaderFileSeries(...,'timeInterval',INT) specifies the
    %   expected or average time in minutes that has elapsed between
    %   consecutive frames.
    %
    %   ImageReaderFileSeries(...,'times',TIMES) specifies a vector of
    %   times in minutes for each frame in the time lapse. Alternatively,
    %   it can be a matrix of times with a row for each point/position and
    %   a column for each frame in the time-lapse data set.

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
        root
        files
    end
    
    properties (Dependent)
        pos
        posName
    end
    
    properties (Transient, Access=private)
        root_val
        files_val
        fileP_val
        fileC_val
        fileZ_val
        fileT_val
        pos_xy
        pos_class
        posC_map
        posZ_map
        posT_map
        pos_file_map
        pos_val
        meta_val = struct()
    end
    
    properties (Constant)
        templates = struct(...
            'posdir_T_channel_Z','<pos>/[^.]*\D<t>_<channel>_<z>',...
            'prefixed','<pos>{t<t>;c<c>;z<z>}')
        var_templates = {...
            '*','.*';...
            '/','[/\\\\]';...
            '(','\\(';...
            ')','\\)';...
            '<t>','(?<t>\\d+)';...
            '<c>','(?<c>\\d+)';...
            '<z>','(?<z>\\d*)';... % star allows for Z sectioning off case
            '<p>','(?<p>\\d+)';...
            '<pos>','(?<pos>[^./\\\\]+)';...
            '<channel>','(?<channel>[^./\\\\]+)'}
    end
    
    methods
        function this = ImageReaderFileSeries(root,varargin)
            
            if nargin<1, return; end
            
            ip = inputParser;
            ischarstr = @(x) ischar(x) && isrow(x);
            isstrlike = @(x) isstring(x) || ischarstr(x) || ...
                (iscell(x) && all(ischarstr,x));
            ip.addRequired('root',@(x) isfolder(x));
            ip.addParameter('ext',{'png','tif','jpg'},isstrlike);
            ip.addParameter('template','posdir_T_channel_Z',ischarstr);
            ip.addParameter('pixelSize',[],@(x) isempty(x) || ...
                (isscalar(x) && isnumeric(x)));
            ip.addParameter('timeInterval',[],@(x) isempty(x) || ...
                (isvector(x) && isnumeric(x)));
            ip.addParameter('times',[],@(x) isempty(x) || ...
                ((isvector(x) || ismatrix(x)) && isnumeric(x)));
            %ip.addParameter('dim_order','tcz',ischarstr);
            %ip.addParameter('dim_prefix',true,@(x) islogical(x) && isscalar(x));
            %ip.addParameter('pos_by_dir',true,@(x) islogical(x) && isscalar(x));
            %ip.addParameter('named_channels',true,@(x) islogical(x) && isscalar(x));
            ip.parse(root,varargin{:});
            
            ext = cellstr(ip.Results.ext);
            template = ip.Results.template;
            % dim_order = ip.Results.dimorder;
            % pos_by_dir = ip.Results.posbydir;
            
            % Get absolute path to root
            d = dir(root);
            root = d(1).folder;
            this.root_val = root;
            
            % Search for all image files
            imgfiles = cell(1,numel(ext));
            for e=1:numel(ext)
                imgfiles{e} = dir(fullfile(root,['**/*.',ext{e}]));
            end
            imgfiles = [imgfiles{:}];
            
            sidx = numel(root)+2;
            imgnames = cell(1,numel(imgfiles));
            imgnames_ext = cell(1,numel(imgfiles));
            for i=1:numel(imgfiles)
                f = imgfiles(i);
                % Remove the extension
                eidx = find(f.name=='.',1,'last')-1;
                % Concatenate with relative folder location
                imgnames{i} = [f.folder(sidx:end),filesep,f.name(1:eidx)];
                imgnames_ext{i} = [f.folder(sidx:end),filesep,f.name];
            end
            
            if isfield(this.templates,template)
                template = this.templates.(template);
            end
            % Format the template into a regex
            % First find the permutation invariant block (if any)
            [permS,permE] = regexp(template,'\{(.+;).+\}','start','end');
            if numel(permS)==1
                elements = strip(strsplit(template(permS+1:permE-1),';'));
                assert(numel(elements)<5,...
                    'Permutation invariant blocks may only contain up to 4 elements');
                P = perms(1:numel(elements));
                template_permuted = cell(1,size(P,1));
                for p=1:numel(template_permuted)
                    template_permuted{p} = [template(1:permS-1),...
                        strjoin(elements(P(p,:)),'.*'),template(permE+1:end)];
                end
                template = template_permuted;
            elseif numel(permS)>1
                warning('Only one permutation invariant block allowed- template will be treated literally');
            end
            % Apply rules to translate to regex
            template = regexprep(template,...
                ImageReaderFileSeries.var_templates(:,1),...
                ImageReaderFileSeries.var_templates(:,2));
            
            matches = regexp(imgnames,template,'once','names');
            % Eliminate any files that do not match the template
            mfilt = ~cellfun(@isempty,matches);
            this.files_val = imgnames_ext(mfilt);
            matches = matches(mfilt);
            assert(all(cellfun(@numel,matches)==1));
            matches = [matches{:}];
            
            % First check for positions
            imgPosNum = [];
            if isfield(matches,'p')
                imgPosNum = str2double({matches.p});
                assert(all(~isnan(imgPosNum) & ...
                    round(imgPosNum)==imgPosNum & imgPosNum>=0),...
                    'position numbers must be valid non-negative integers');
            end
            % Adjust any zero-indexing to one-indexing
            if any(imgPosNum==0), imgPosNum = imgPosNum + 1; end
            
            imgPosName = {};
            if isfield(matches,'pos')
                imgPosName = {matches.pos};
            end
            if ~isempty(imgPosNum) && ~isempty(imgPosName)
                % Ensure a one-to-one correspondence with pos numbers
                posmap = arrayfun(@(x) unique(imgPosName(imgPosNum==x)), ...
                    unique(imgPosNum),'uni',0);
                assert(all(arrayfun(@numel,posmap)==1) && ...
                    numel(unique([posmap{:}]))==numel(posmap), ...
                    'pos names must correspond with pos numbers');
            end
            
            if isempty(imgPosNum) && isempty(imgPosName)
                imgPosNum = ones(size(matches));
            end
            
            if isempty(imgPosName)
                [posNums,~,this.fileP_val] = unique(imgPosNum);
                this.meta_val.posNames = ...
                    arrayfun(@(p) sprintf('pos%03d',p),posNums,'Uniform',false);
            else
                % We have position names and possibly also numbers.
                [posNames,firstOccurs,posIdx] = unique(imgPosName);
                if ~isempty(imgPosNum)
                    % Use corresponding numbers to sort posNames
                    [~,sortIdx] = sort(imgPosNum(firstOccurs));
                    posNames = posNames(sortIdx);
                    [~,revSortIdx] = sort(sortIdx);
                    posIdx = revSortIdx(posIdx);
                end
                this.fileP_val = posIdx;
                this.meta_val.posNames = posNames;
            end
            
            % Then check for channels
            imgChNum = [];
            if isfield(matches,'c')
                imgChNum = str2double({matches.c});
                assert(all(~isnan(imgChNum) & ...
                    round(imgChNum)==imgChNum & imgChNum>=0),...
                    'channel numbers must be valid non-negative integers');
            end
            % Adjust any zero-indexing to one-indexing
            if any(imgChNum==0), imgChNum = imgChNum + 1; end
            
            imgChName = {};
            if isfield(matches,'channel')
                imgChName = {matches.channel};
            end
            if ~isempty(imgChNum) && ~isempty(imgChName)
                % Ensure a one-to-one correspondence with pos numbers
                chmap = arrayfun(@(x) unique(imgChName(imgChNum==x)), ...
                    unique(imgChNum),'uni',0);
                assert(all(arrayfun(@numel,chmap)==1) && ...
                    numel(unique([chmap{:}]))==numel(chmap), ...
                    'channel names must correspond with channel numbers');
            end
            
            if isempty(imgChNum) && isempty(imgChName)
                imgChNum = ones(size(matches));
            end
            
            if isempty(imgChName)
                [chNums,~,this.fileC_val] = unique(imgChNum);
                this.meta_val.channelNames = ...
                    arrayfun(@(c) sprintf('channel%02d',c),chNums,'Uniform',false);
            else
                % We have channel names and possibly also numbers.
                [chNames,firstOccurs,chIdx] = unique(imgChName);
                if ~isempty(imgChNum)
                    % Use corresponding numbers to sort channel names
                    [~,sortIdx] = sort(imgChNum(firstOccurs));
                    chNames = chNames(sortIdx);
                    [~,revSortIdx] = sort(sortIdx);
                    chIdx = revSortIdx(chIdx);
                end
                this.fileC_val = chIdx;
                this.meta_val.channelNames = chNames;
            end
            this.fileC_val = this.fileC_val(:);
            
            % Check for Z index
            imgZ = [];
            if isfield(matches,'z')
                imgZ = {matches.z};
                % If Z is empty, then we assume Z sectioning is off for
                % this channel and set Z = 1
                imgZ(cellfun(@isempty,imgZ)) = {'1'};
                imgZ = str2double(imgZ);
                assert(all(~isnan(imgZ) & ...
                    round(imgZ)==imgZ & imgZ>=0),...
                    'Z indices must be valid non-negative integers');
            end
            % Adjust any zero-indexing to one-indexing
            if any(imgZ==0), imgZ = imgZ + 1; end
            if isempty(imgZ)
                imgZ = ones(size(matches));
            end
            this.fileZ_val = imgZ(:);
            
            % Check for T index
            imgT = [];
            if isfield(matches,'t')
                imgT = str2double({matches.t});
                assert(all(~isnan(imgT) & ...
                    round(imgT)==imgT & imgT>=0),...
                    'time point indices must be valid non-negative integers');
            end
            % Adjust any zero-indexing to one-indexing
            if any(imgT==0), imgT = imgT + 1; end
            if isempty(imgT)
                imgT = ones(size(matches));
            end
            this.fileT_val = imgT(:);
            
            if ~isempty(ip.Results.pixelSize)
                this.meta_val.pixelSize = ip.Results.pixelSize;
            end
            
            if ~isempty(ip.Results.times)
                posTimes = ip.Results.times;
                if isvector(posTimes)
                    posTimes = posTimes(:)';
                    posTimes = posTimes(ones(this.npos,1),:);
                end
                assert(size(posTimes,1)==this.npos,...
                    'Number of rows in times matrix must match number of positions');
                assert(size(posTimes,2)>=max(imgT),...
                    'Number of columns in times matrix must be at least as large as the maximum time point')
                this.meta_val.posTimes = posTimes;
            end
            if ~isempty(ip.Results.timeInterval)
                posTimeIntervals = ip.Results.timeInterval;
                if isscalar(posTimeIntervals)
                    posTimeIntervals = posTimeIntervals(ones(this.npos,1));
                end
                assert(numel(posTimeIntervals)==this.npos,...
                    'Number time intervals must match number of positions');
                this.meta_val.posTimeIntervals = posTimeIntervals;
            end
            
            % Start at position 1
            this.pos = 1;
        end
        
        function val = get.meta(this)
            val = this.meta_val;
        end
        
        function val = get.root(this)
            val = this.root_val;
        end
        
        function val = get.files(this)
            val = this.files_val;
        end
        
        function val = get.npos(this)
            val = numel(this.posNames);
        end
        
        function val = get.posNames(this)
            val = this.meta_val.posNames;
        end
        
        function val = get.pos(this)
            val = this.pos_val;
        end
        
        function set.pos(this,val)
            assert(ismember(val,1:this.npos),...
                '"pos" must be specified as an index into posNames');
            
            posmask = this.fileP_val == val;
            refimg = imread(fullfile(this.root,this.files{find(posmask,1)}));
            this.pos_xy = size(refimg);
            this.pos_class = class(refimg);
            
            posZ_vals = unique(this.fileZ_val(posmask));
            this.posZ_map = zeros(max(posZ_vals),1);
            this.posZ_map(posZ_vals) = 1:numel(posZ_vals);
            posC_vals = unique(this.fileC_val(posmask));
            this.posC_map = zeros(max(posC_vals),1);
            this.posC_map(posC_vals) = 1:numel(posC_vals);
            posT_vals = unique(this.fileT_val(posmask));
            this.posT_map = zeros(max(posT_vals),1);
            this.posT_map(posT_vals) = 1:numel(posT_vals);
            
            this.pos_file_map = zeros(numel(posZ_vals),...
                numel(posC_vals),numel(posT_vals));
            this.pos_file_map(sub2ind(size(this.pos_file_map),...
                this.posZ_map(this.fileZ_val(posmask)),...
                this.posC_map(this.fileC_val(posmask)),...
                this.posT_map(this.fileT_val(posmask)))) = find(posmask);
            
            this.pos_val = val;
        end
        
        function val = get.posName(this)
            val = this.posNames{this.pos};
        end
        
        function set.posName(this,val)
            assert(sum(strcmp(val,this.posNames))==1,...
                'posName not found');
            this.pos = find(strcmp(this.posNames,val),1);
        end
        
        function val = get.channels(this)
            val = this.meta_val.channelNames(this.posC_map>0);
        end
        
        function val = get.nchannels(this)
            val = numel(this.channels);
        end
        
        function val = get.imageSize(this)
            val = [this.pos_xy,numel(this.posZ_map),this.nchannels,...
                numel(this.posT_map)];
        end
        
        function val = get.pixelSize(this)
            val = [];
            if isfield(this.meta_val,'pixelSize')
                val = this.meta_val.pixelSize;
            end
        end
        
        function val = get.timeInterval(this)
            if isfield(this.meta,'posTimeIntervals')
                val = this.meta.posTimeIntervals(this.pos);
            elseif isfield(this.meta,'posTimes')
                val = mean(diff(this.meta.posTimes(this.pos,:)));
            else
                val = NaN;
            end
        end
        
        function val = get.times(this)
            if isfield(this.meta,'posTimes')
                val = this.meta.posTimes(this.pos,:);
            elseif isfield(this.meta,'posTimeIntervals')
                val = this.timeInterval * (0:this.imageSize(5)-1);
            else
                val = NaN(1,this.imageSize(5));
            end
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
            
            sz = this.imageSize;
            
            if isempty(Z)
                Z = 1:sz(3);
            end
            if isempty(C)
                C = 1:sz(4);
            end
            
            assert(T>0 && T<=sz(5),'Time T is out of range');
            assert(all(C>0 & C<=sz(4)),'At least one channel is out of range');
            assert(all(Z>0 & Z<=sz(3)),'At least one Z section is out of range');
                
            val = zeros(sz(1),sz(2),numel(Z),numel(C),this.pos_class);
            for c=1:numel(C)
                for z=1:numel(Z)
                    posZ = this.posZ_map(Z(z));
                    if posZ<=0, continue; end
                    posT = this.posT_map(T);
                    if posT<=0, continue; end
                    fInd = this.pos_file_map(posZ,C(c),posT);
                    if fInd<=0, continue; end
                    val(:,:,z,c) = imread(fullfile(this.root,this.files{fInd}));
                end
            end
        end
    end
end