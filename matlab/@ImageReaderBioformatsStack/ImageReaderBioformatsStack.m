classdef ImageReaderBioformatsStack < ImageReader
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
    end
    
    properties (Dependent)
        pos
        posName
    end
     
    properties (Access=private)
        stack_along
        num_per_stack
    end
    
    properties (Access=private,Transient)
        readers
        imageSizes
        allChannelNames
        pos_val = 1
    end
    
    methods
        function this = ImageReaderBioformatsStack(filenames,varargin)
            %IMAGEREADERBIOFORMATSSTACK Stacks Bioformats files/positions
            %
            %   IMAGEREADERBIOFORMATSSTACK(FILENAMES) stacks each of the
            %   image files FILENAMES (cellstr or string) position-wise
            %   along the time dimension.
            %
            %   IMAGEREADERBIOFORMATSSTACK(...,'StackAlong',D) stacks
            %   alternatively along dimension D, which can be either 'T',
            %   'Z' or 'C'.
            %
            %   IMAGEREADERBIOFORMATSSTACK(...,'PosesPerStack',N) for N > 0 
            %   creates new positions from a single image file by stacking
            %   groups of N consecutive positions along dimension D. The
            %   number of positions in the file should be a multiple of N.
            
            ip = inputParser;
            ischarlike = @(x) ischar(x) || iscellstr(x) || isstring(x);
            isposint = @(x) isnumeric(x) && all(round(x(:))==x(:)) && all(x(:)>=0);
            ip.addRequired('filenames',@(x) ischarlike(x));
            ip.addOptional('poses',[],@(x) isempty(x) || (isposint(x) && isrow(x)) || ...
                (iscell(x) && all(cellfun(@(y) isposint(y) && isrow(y),x))));
            ip.addParameter('StackAlong','T',@(x) ismember(x,{'T','Z','C'}));
            ip.addParameter('PosesPerStack',0,@(x) isposint(x) && isscalar(x));
            ip.parse(filenames,varargin{:});
            
            filenames = cellstr(filenames);
            poses = ip.Results.poses;
            if ~iscell(poses), poses = {poses}; end
            if isscalar(poses)
                poses = poses(ones(1,numel(filenames)));
            end
            assert(numel(filenames)==numel(poses),...
                'Number of "poses" vectors must equal the number of "filenames"');
            
            poses_per_stack = ip.Results.PosesPerStack;
            if poses_per_stack > 0
                assert(numel(filenames)==1,...
                    'If NumPerStack>0, then only one file can be specified');
            end
            
            this.num_per_stack = poses_per_stack;
            dim_map = struct('T',5,'C',4,'Z',3);
            this.stack_along = dim_map.(ip.Results.StackAlong);
            
            this.readers = cell(1,numel(filenames));
            for f=1:numel(filenames)
                this.readers{f} = ImageReaderBioformats(filenames{f});
            end
            r1 = this.readers{1};
            assert(all(cellfun(@(x) all(x.pixelSize==r1.pixelSize),this.readers)),...
                'All image files must have the same pixel size');
            if poses_per_stack > 0
                assert(mod(r1.npos,poses_per_stack)==0,...
                    'If NumPerStack>0, then number of positions must be a multiple');
            else
                assert(all(cellfun(@(x) x.npos==r1.npos,this.readers)),...
                    'Every file must have the same number of positions');
            end
            imageSizes = cell(1,numel(filenames));
            for f=1:numel(filenames)
                r = this.readers{f};
                imageSizes{f} = zeros(r.npos,5);
                for p=1:r.npos
                    r.pos = p;
                    imageSizes{f}(p,:) = r.imageSize;
                end
            end
            this.imageSizes = imageSizes;
            imageSizeCheck = vertcat(imageSizes{:});
            % Allow for differing numbers of channels (these can be filled
            % with zeros)
            imageSizeCheck(:,[4,this.stack_along]) = [];
            assert(all(all(imageSizeCheck==imageSizeCheck(1,:))),...
                'All image dimensions except "StackAlong" must match');
        end
        
        function val = get.npos(this)
            val = this.readers{1}.npos;
            if this.num_per_stack > 0
                val = val / this.num_per_stack;
            end
        end
        
        function val = get.posNames(this)
            val = this.readers{1}.posNames;
            if this.num_per_stack > 0
                val = val(1:this.num_per_stack:this.readers{1}.npos);
            end
        end
        
        function val = get.pos(this)
            val = this.pos_val;
        end
        
        function set.pos(this,val)
            isposint = @(x) isnumeric(x) && all(round(x(:))==x(:)) && all(x(:)>=0);
            if ~isposint(val) || val<1 || val>this.npos
                error('"pos" must be an integer between 1 and %u',this.npos);
            end
            this.pos_val = val;
        end
        
        function val = get.posName(this)
            val = this.posNames{this.pos};
        end
        
        function set.posName(this,val)
            p = find(strcmp(this.posNames,val),1);
            if isempty(p)
                error('That position could not be found');
            end
            this.pos = p;
        end
        
        function val = get.channels(this)
            if isempty(this.allChannelNames)
                chnames = cell(1,numel(this.readers));
                for i=1:numel(this.readers)
                    r = this.readers{i};
                    chnames = cell(1,r.npos);
                    for p=1:r.npos
                        r.pos = p;
                        chnames{i}{p} = r.channels(:);
                    end
                    chnames{i} = vertcat(chnames{i}{:});
                end
                chnames = vertcat(chnames{:});
                this.allChannelNames = unique(chnames,'stable')';
            end
            val = this.allChannelNames;
        end
        
        function val = get.nchannels(this)
            val = numel(this.channels);
        end
        
        function val = get.imageSize(this)
            poslist = this.getCurrentPosList;
            imszs = zeros(size(poslist,1),5);
            for i=1:size(poslist,1)
                poslist{i,1}.pos = poslist{i,2};
                imszs(i,:) = poslist{i,1}.imageSize;
            end
            val = imszs(1,:);
            val(this.stack_along) = sum(imszs(:,this.stack_along));
        end
        
        function val = get.pixelSize(this)
            val = this.readers{1}.pixelSize;
        end
        
        function val = get.timeInterval(this)
            poslist = this.getCurrentPosList;
            val = poslist{1,1}.timeInterval;
        end
        
        function val = get.times(this)
            poslist = this.getCurrentPosList;
            val = poslist{1,1}.times;
        end
        
        function val = get.meta(this)
            val = cellfun(@(x) x.meta,this.readers,'Uniform',false);
            try
                val = [val{:}];
            catch
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
            
            poslist = this.getCurrentPosList;
            imszs = zeros(size(poslist,1),5);
            for i=1:size(poslist,1)
                poslist{i,1}.pos = poslist{i,2};
                imszs(i,:) = poslist{i,1}.imageSize;
            end
            imsz = imszs(1,:);
            imsz(this.stack_along) = sum(imszs(:,this.stack_along));
            if ~isempty(Z) && (Z<1 || Z>imsz(3))
                error('Z must lie between 1 and %u',imsz(3));
            end
            if ~isempty(C) && (C<1 || C>imsz(4))
                error('C must lie between 1 and %u',imsz(4));
            end
            if T<1 || T>imsz(5)
                error('T must lie between 1 and %u',imsz(5));
            end
            
            index = cell(1,5);
            index(3:5) = {Z,C,T};
            if isempty(index{this.stack_along})
                val = cell(1,size(poslist,1));
                for i=1:size(poslist,1)
                    poslist{i,1}.pos = poslist{i,2};
                    val{i} = poslist{i,1}.getTimepoint(T,Z,'C',C);
                end
                val = cat(this.stack_along,val{:});
            else
                stack_lims = cumsum(imszs(:,this.stack_along));
                i = find(index{this.stack_along}<=stack_lims,1);
                stack_lims = [0,stack_lims];
                index{this.stack_along} = index{this.stack_along} - stack_lims(i);
                poslist{i,1}.pos = poslist{i,2};
                val = poslist{i,1}.getTimepoint(index{5},index{3},'C',index{4});
            end
        end
    end
    
    methods (Access=private)
        function val = getCurrentPosList(this)
            if this.num_per_stack > 0
                val = cell(this.num_per_stack,2);
                for i=1:this.num_per_stack
                    val{i,1} = this.readers{1};
                    val{i,2} = (this.pos-1)*this.num_per_stack+i;
                end
            else
                val = cell(numel(this.readers),2);
                for i=1:numel(this.readers)
                    val{i,1} = this.readers{i};
                    val{i,2} = this.pos;
                end
            end
        end
    end
end