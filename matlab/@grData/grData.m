classdef grData < handle
    properties
        M = struct()
        D = {}
        args
        meta = struct()
    end
    
    properties (Dependent)
        filename
        Ml % cells with expanded lineage
        mother_inds
    end
    
    properties (Access=private)
        filename_var
        save_D = true
        log_gr
        shortTPthresh
    end
    
    properties (Access=private,Transient)
        Ml_val = []
        mother_inds_val = []
    end
    
    methods
        function this = grData(cellInf,varargin)
            if nargin<1 || isempty(cellInf)
                return
            end
            % Save args, but ensure that any sparse arguments are made full
            this.args = varargin;
            sparse_args = cellfun(@issparse,this.args);
            this.args(sparse_args) = cellfun(@full,this.args(sparse_args),'Uniform',false);
            
            ip = inputParser;
            ip.addParameter('save_D',true,@(x) isscalar(x) && islogical(x));
            ip.addParameter('log_gr',false,@(x) isscalar(x) && islogical(x));
            ip.addParameter('shortTPthresh',5,@(x) isscalar(x) && isnumeric(x));
            ip.addParameter('extract',{},@(x) iscellstr(x) || isstring(x));
            ip.KeepUnmatched = true;
            ip.parse(varargin{:});
            
            this.save_D = ip.Results.save_D;
            this.log_gr = ip.Results.log_gr;
            this.shortTPthresh = ip.Results.shortTPthresh;
            extract = ip.Results.extract;
            unmatched = [fieldnames(ip.Unmatched),struct2cell(ip.Unmatched)]';
            
            % Retain any non-numeric fields as meta data
            fnames = fieldnames(cellInf);
            is_num_field = structfun(@(x) isnumeric(x) || islogical(x),cellInf(1));
            this.meta = rmfield(cellInf(1),fnames(is_num_field));
            % and do not propagate to grsFromCellInf if specified
            extract = setdiff(extract,fnames(~is_num_field));
            
            % The extractionProperties field of the cellInf has a function
            % handle that we need to convert to text
            if isfield(this.meta,'extractionParameters') ...
                    && isfield(this.meta.extractionParameters,'extractFunction') ...
                    && isa(this.meta.extractionParameters.extractFunction,'function_handle')
                this.meta.extractionParameters.extractFunction = ...
                    func2str(this.meta.extractionParameters.extractFunction);
            end
            
            [this.M,this.D] = grsFromCellInf(cellInf,...
                'log_gr',this.log_gr,'shortTPthresh',this.shortTPthresh,...
                'extract',extract,unmatched{:});
        end
        
        function saveData(this,deflate)
            fname = this.filename;
            if nargin<2 || isempty(deflate), deflate = 4; end
            
            % Always write to a new file
            if exist(fname,'file') == 2, delete(fname); end
            
            % Save args and meta
            %h5writeatt(fname,'/','args',jsonencode(this.args));
            %h5writeatt(fname,'/','meta',jsonencode(this.meta));
            argsjson = uint16(jsonencode(this.args));
            sz = size(argsjson);
            h5create(fname,'/args',sz,'Deflate',4,...
                'ChunkSize',[min(sz(1),50),min(sz(2),50)],...
                'Datatype',class(argsjson));
            h5write(fname,'/args',argsjson);
            metajson = uint16(jsonencode(this.meta));
            sz = size(metajson);
            h5create(fname,'/meta',sz,'Deflate',4,...
                'ChunkSize',[min(sz(1),50),min(sz(2),50)],...
                'Datatype',class(argsjson));
            h5write(fname,'/meta',metajson);
            
            % Write mother data
            mId = '/M/';
            fields = fieldnames(this.M);
            for f=1:numel(fields)
                field = fields{f};
                var = full(this.M.(field));
                if islogical(var), var = uint8(var); end
                sz = size(var);
                chnksz = arrayfun(@(x) min(x,50),sz);
                h5create(fname,strcat(mId,field),sz,'Datatype',class(var),...
                    'Deflate',deflate,'ChunkSize',chnksz);
                h5write(fname,strcat(mId,field),var);
            end
            
            if this.save_D
                % Write daughter data
                for m=1:numel(this.D)
                    dId = sprintf('/D/%u/',m);
                    if isempty(this.D{m}), continue; end
                    fields = fieldnames(this.D{m});
                    for f=1:numel(fields)
                        field = fields{f};
                        var = full(vertcat(this.D{m}.(field)));
                        if islogical(var), var = uint8(var); end
                        if isempty(var), continue; end
                        sz = size(var);
                        h5create(fname,strcat(dId,field),sz,'Deflate',deflate,...
                            'ChunkSize',[min(sz(1),50),min(sz(2),50)]);
                        h5write(fname,strcat(dId,field),var);
                    end
                end
            end
        end
        
        function flat_dVols = flatten_dVols(this,nmin)
            if nargin<2 || isempty(nmin), nmin = this.shortTPthresh; end
            flat_dVols = NaN(numel(this.D),...
                size(this.D{find(cellfun(@(x) isfield(x,'vol'),this.D),1)}(1).vol,2));
            for m=1:numel(this.D)
                if ~isfield(this.D{m},'vol'), continue; end
                for d=1:numel(this.D{m})
                    dv = this.D{m}(d).vol;
                    valid = ~isnan(dv);
                    if sum(valid)<nmin, continue; end
                    flat_dVols(m,valid) = dv(valid);
                end
            end
        end
        
        function flat_dGRs = flatten_dGRs(this,nmin)
            if nargin<2 || isempty(nmin), nmin = this.shortTPthresh; end
            islog = this.log_gr;
            flat_dGRs = NaN(numel(this.D),...
                size(this.D{find(cellfun(@(x)...
                ~isempty(x) && isfield(x,'svol'),this.D),1)}(1).svol,2));
            for m=1:numel(this.D)
                if ~isfield(this.D{m},'svol') || ~isfield(this.D{m},'grs'), continue; end
                for d=1:numel(this.D{m})
                    if islog
                        gr = exp(this.D{m}(d).svol).*this.D{m}(d).grs;
                    else
                        gr = this.D{m}(d).grs;
                    end
                    valid = ~isnan(gr);
                    if sum(valid)<nmin, continue; end
                    flat_dGRs(m,valid) = gr(valid);
                end
            end
        end
        
        function var = get.filename(this)
            if isempty(this.filename_var)
                error('"filename" has not been set');
            end
            var = this.filename_var;
        end
        
        function set.filename(this,var)
            assert(isvector(var) && ischar(var) && endsWith(var,{'.h5','.hdf5'}),...
                '"filename" must be a string with extension ".h5" or ".hdf5"');
            this.filename_var = var;
        end
        
        function val = get.Ml(this)
            if ~isempty(this.Ml_val)
                val = this.Ml_val;
                return
            end
            if ~isfield(this.M,'mothers')
                error('"mothers" needs to be included in extraction');
            end
            mothers = this.M.mothers;
            mothers(mothers<=0) = NaN;
            posNum = this.M.posNum;
            trapNum = this.M.trapNum;
            cellNum = this.M.cellNum;
            motherids = [posNum(:),trapNum(:),mothers(:)];
            if any(~ismember(motherids(~isnan(mothers),:),...
                    [posNum(:),trapNum(:),cellNum(:)],'rows'))
                warning('Mothers are missing. All their descendents will be NaN.');
            end
            % Generate indices of entire lineage for each cell
            mothermap = NaN(max(cellNum),max(trapNum),max(posNum));
            mothermap(sub2ind(size(mothermap),cellNum,trapNum,posNum)) = mothers;
            indmap = NaN(max(cellNum),max(trapNum),max(posNum));
            indmap(sub2ind(size(indmap),cellNum,trapNum,posNum)) = 1:numel(cellNum);
            lineages = cell(numel(cellNum),1);
            for c=1:numel(cellNum)
                linlocal = recurseLineage(cellNum(c),...
                    mothermap(:,trapNum(c),posNum(c)));
                linlocal(linlocal>size(indmap,1)) = NaN;
                linfilt = ~isnan(linlocal);
                linlocal(linfilt) = indmap(linlocal(linfilt),trapNum(c),posNum(c))';
                lineages{c} = linlocal;
            end
            
            % Expand lineage for each row of each matrix field in M
            fnames = fieldnames(this.M);
            fsz = size(this.M.vol); % field guaranteed by grsFromCellInf
            isvalid = ~isnan(this.M.vol);
            fnames = fnames(structfun(@(x) isequal(size(x),fsz),this.M));
            val = struct();
            for f=1:numel(fnames)
                fname = fnames{f};
                srcmat = this.M.(fname);
                if isa(srcmat,'double')
                    outmat = NaN(size(srcmat));
                else
                    outmat = zeros(size(srcmat),'like',srcmat);
                end
                for c=1:numel(cellNum)
                    if any(isnan(lineages{c}))
                        continue;
                    end
                    for li=1:numel(lineages{c})
                        l = lineages{c}(li);
                        stp = find(isvalid(l,:),1); % start tp
                        outmat(c,stp:end) = srcmat(l,stp:end);
                        if strcmp(fname,'births') && li>1
                            % This cell has a mother that becomes a sister and
                            % in a lineage tree, we should record the shared
                            % birth point:
                            outmat(c,stp) = true;
                        end
                    end
                end
                val.(fname) = outmat;
            end
            
            % Generate a sisterLabel array that is similar to
            % daughterLabel, but will not be NaN for the nominal
            % daughter as soon as the sister pair is born:
            srcmat = val.daughterLabel;
            outmat = NaN(size(srcmat));
            rawbirths = this.M.births;
            rawbirths(:,end+1) = true; % ensures at least one birth below
            for c=1:numel(cellNum)
                if any(isnan(lineages{c})), continue; end
                for li=1:numel(lineages{c})
                    l = lineages{c}(li);
                    stp = find(isvalid(l,:),1);
                    outmat(c,stp:end) = srcmat(l,stp:end);
                    if li>1
                        % This cell has a mother that becomes a sister
                        etp = find(isvalid(l,:),1,'last');
                        next_btp = find(rawbirths(l,:),1);
                        sistps = stp:min(next_btp-1,etp);
                        outmat(c,sistps) = cellNum(lineages{c}(li-1));
                    end
                end
            end
            val.sisterLabel = outmat;
            
            this.Ml_val = val;
            
            function lineage = recurseLineage(clbl,mothers,lineage)
                if nargin<3, lineage = []; end
                if ismember(clbl,lineage)
                    warning('a lineage has a cycle- is a mother to itself');
                    lineage = [NaN,lineage];
                    return
                end
                lineage = [clbl,lineage];
                if clbl > numel(mothers) || isnan(mothers(clbl))
                    return
                else
                    lineage = recurseLineage(mothers(clbl),mothers,lineage);
                end
            end
        end
        
        function val = get.mother_inds(this)
            if isempty(this.mother_inds_val)
                assert(numel(this.M.posNum)<2^16-1,...
                    'Too many cells; increased bit depth required...');
                pnum = this.M.posNum(:); tnum = this.M.trapNum(:);
                cnum = this.M.cellNum(:); mnum = this.M.mothers(:);
                cell_map = zeros(max(pnum),max(tnum),max(cnum),'uint16');
                cell_map(sub2ind(size(cell_map),pnum,tnum,cnum)) = 1:numel(pnum);
                this.mother_inds_val = NaN(1,numel(pnum));
                hasm = ~isnan(mnum) & mnum > 0 & mnum <= size(cell_map,3);
                this.mother_inds_val(hasm) = cell_map(sub2ind(size(cell_map),...
                    pnum(hasm),tnum(hasm),mnum(hasm)));
                if any(this.mother_inds_val(hasm) == 0)
                    warning('some mothers are missing...');
                end
            end
            val = this.mother_inds_val;
        end
    end
    
    methods (Static)
        function this = loadFrom(fname)
            this = grData;
            this.filename = fname;
            
            % Read args and meta
            hinf = h5info(fname,'/');
            attrs = {};
            if ~isempty(hinf.Attributes)
                attrs = {hinf.Attributes.Name};
            end
            datasets = {};
            if ~isempty(hinf.Datasets)
                datasets = {hinf.Datasets.Name};
            end
            groups = {};
            if ~isempty(hinf.Groups)
                groups = {hinf.Groups.Name};
            end
            
            % TODO: the arguments here should be loaded to set some private
            % variables (e.g., save_D = true, log_gr, shortTPthresh)
            if ismember('args',attrs)
                args = h5readatt(fname,'/','args');
                if iscell(args) && numel(args) == 1, args = args{1}; end
                this.args = jsondecode(args);
            elseif ismember('args',datasets)
                this.args = jsondecode(char(h5read(fname,'/args')));
            end
            
            if ismember('meta',attrs)
                meta = h5readatt(fname,'/','meta');
                if iscell(meta) && numel(meta) == 1, meta = meta{1}; end
                this.meta = jsondecode(meta);
            elseif ismember('meta',datasets)
                this.meta = jsondecode(char(h5read(fname,'/meta')));
            end
            
            % Read mother data
            hinf_m = hinf.Groups(strcmp(groups,'/M'));
            dnames = {hinf_m.Datasets.Name};
            for d=1:numel(dnames)
                dname = dnames{d};
                this.M.(dname) = h5read(fname,strcat('/M/',dname));
                if isa(this.M.(dname),'uint8') && max(this.M.(dname)(:))<=1
                    this.M.(dname) = logical(this.M.(dname));
                end
            end
            
            % Read daughter data
            if ~ismember('/D',groups)
                this.save_D = false;
                return
            end
            hinf_d = hinf.Groups(strcmp(groups,'/D'));
            mIds = cellfun(@str2double,...
                regexprep({hinf_d.Groups.Name},'^/D/(\d+)$','$1'));
            this.D = cell(size(this.M.vol,1),1);
            for m=mIds
                this.D{m} = struct();
                dId = sprintf('/D/%u',m);
                hinf_d_m = hinf_d.Groups(strcmp({hinf_d.Groups.Name},dId));
                dnames = {hinf_d_m.Datasets.Name};
                for n=1:numel(dnames)
                    dname = dnames{n};
                    ds = h5read(fname,strcat(dId,'/',dname));
                    for d=1:size(ds,1)
                        this.D{m}(d).(dname) = ds(d,:);
                    end
                end
            end
        end
    end
end
