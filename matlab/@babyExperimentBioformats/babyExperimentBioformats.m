classdef babyExperimentBioformats < babyExperiment
    % BABYEXPERIMENTBIOFORMATS a subclass of babyExperiment that
    % manages cExperiments that access data from files loaded by the
    % Bio-Formats library. Requires the OME Bio-Formats Matlab Toolbox,
    % which can be downloaded from:
    % https://www.openmicroscopy.org/bio-formats/downloads/
    %
    % See also BABYEXPERIMENT, BABYTIMELAPSEBIOFORMATS
    
    properties
        posImgFiles
        posImgIndices
        stackingArgs
    end
    
    methods
        function cExperiment = babyExperimentBioformats(imgfile,savefolder,varargin)
            %BABYEXPERIMENTBIOFORMATS Create a new cExperiment
            %   - imgfile (optional): specify a file(s) that can be loaded by 
            %       BioFormats.
            %   - savefolder (optional): specify a folder in which to save
            %       the cExperiment and cTimelapses.
            
            % call super class constructor such that it does not initialise
            % anything.
            cExperiment@babyExperiment(true);
            
            % If initialised with single input true, make an empty object
            % (used in load function)
            if nargin==1 && islogical(imgfile) && imgfile
                % if folder is true, cExperiment returned bare for load
                % function.
                return
            end
            
            ip = inputParser;
            ischarlike = @(x) ischar(x) || iscellstr(x) || isstring(x);
            isposint = @(x) isnumeric(x) && all(round(x(:))==x(:)) && all(x(:)>=0);
            ip.addRequired('imgfile',@(x) ischarlike(x));
            ip.addRequired('savefolder',@(x) ischar(x) && isrow(x));
            ip.addParameter('StackAlong',[],@(x) isempty(x) || ismember(x,{'T','Z','C'}));
            ip.addParameter('PosesPerStack',[],@(x) isempty(x) || (isposint(x) && isscalar(x)));
            ip.parse(imgfile,savefolder,varargin{:});
            
            cExperiment.rootFolder = imgfile;
            cExperiment.saveFolder = savefolder;
            
            % Default to enabling logging:
            cExperiment.shouldLog=true;
            % NB: constructor no longer needs to instantiate the logger
            % property, see babyExperiment
            
            %Record the user who is creating the cExperiment
            if ispc
                cExperiment.creator=getenv('USERNAME');
            else
                [~, cExperiment.creator] = system('whoami');
            end
            %Initialize records of positions segmented and tracked
            cExperiment.posTrapsTracked=0;
            cExperiment.posSegmented=0;
            cExperiment.posTracked=0;
            
            cExperiment.cellsToPlot=cell(1);
            
            stack_along = ip.Results.StackAlong;
            poses_per_stack = ip.Results.PosesPerStack;
            is_stacked = false;
            if ~isempty(stack_along) || ~isempty(poses_per_stack)
                is_stacked = true;
                if isempty(stack_along)
                    stack_along = 'T';
                end
                if isempty(poses_per_stack)
                    poses_per_stack = 0;
                end
                cExperiment.stackingArgs = struct('StackAlong',stack_along,...
                    'PosesPerStack',poses_per_stack);
            end
            
            % Ensure that imgfile is a cellstr
            imgfile = cellstr(imgfile);
            if is_stacked
                imgfile = {imgfile};
            end
            
            % Collect meta data from all provided input files
            pixelSize = []; % pixel size needs to be consistent
            baseChannelNames = {}; % collect union of all channel names
            nposes = NaN(numel(imgfile),1);
            posNames = cell(numel(imgfile),1);
            posTimes = cell(numel(imgfile),1);
            maxnz = 1;
            for f=1:numel(imgfile)
                % Start a temporary reader to get meta data
                if is_stacked
                    reader = ImageReaderBioformatsStack(imgfile{f},...
                        'StackAlong',stack_along,'PosesPerStack',poses_per_stack);
                else
                    reader = ImageReaderBioformats(imgfile{f});
                end
                baseChannelNames = union(baseChannelNames,reader.channels);
                maxnz = max(maxnz,reader.imageSize(3));
                if isempty(pixelSize)
                    pixelSize = reader.pixelSize;
                else
                    assert(isequal(pixelSize,reader.pixelSize),...
                        'Pixel sizes must not differ between image files');
                end
                nposes(f) = reader.npos;
                posNames{f} = reader.posNames;
                posTimes{f} = cell(1,nposes(f));
                for p=1:nposes(f)
                    reader.pos = p;
                    posTimes{f}{p} = reader.times/60; % times in minutes
                end
            end
            cExperiment.posImgFiles = imgfile(repelem(1:numel(imgfile),nposes));
            inds = arrayfun(@(n) 1:n,nposes,'Uniform',false);
            cExperiment.posImgIndices = [inds{:}];
            nposes = numel(cExperiment.posImgFiles);
            cExperiment.dirs = arrayfun(@(x) sprintf('pos%03u',x),...
                1:nposes,'Uniform',false);
            
            channelNames = cell(numel(baseChannelNames),1);
            for c=1:numel(channelNames)
                channelNames{c} = cell(maxnz,1);
                for z=1:maxnz
                    channelNames{c}{z} = sprintf('%s_%03u',baseChannelNames{c},z);
                end
            end
            
            cExperiment.channelNames = vertcat(baseChannelNames{:},channelNames{:});
            
            % Parse the microscope acquisition metadata and populate the 
            % metadata field of cExperiment
            % The following is the minimum detail required for a 
            % cExperiment with a single position:
            cExperiment.metadata = struct();
            if is_stacked, imgfile = imgfile{1}; end
            cExperiment.metadata.acqFile = strjoin(imgfile,';'); % For ID
            cExperiment.metadata.acq = struct(...
                'channels',struct('names',{cExperiment.channelNames},...
                'zsect',zeros(size(cExperiment.channelNames))),...
                'zsections', struct('sections',1),'positions',zeros(0,1));
            
            cExperiment.metadata.originalPosNames = [posNames{:}];
            
            posTimes = [posTimes{:}];
            nT = max(cellfun(@numel,posTimes));
            cExperiment.metadata.posTimes = NaN(nposes,nT);
            for p=1:nposes
                cExperiment.metadata.posTimes(p,1:numel(posTimes{p})) = posTimes{p};
            end
            
            % See if we can determine a date from file name
            date_matches = regexp(cExperiment.metadata.acqFile,{...
                '\D(\d{4})(\d{2})(\d{2})\D','\D(\d{4})\D(\d{2})\D(\d{2})\D'},'tokens');
            is_matching = ~cellfun(@isempty,date_matches);
            if any(is_matching)
                date_match = date_matches{find(is_matching,1)};
                cExperiment.metadata.date = strjoin(date_match{1},'-');
            else
                date_match = inputdlg({'Enter date of acquistion in YYYY-MM-DD format:'},'Enter date');
                assert(~isempty(date_match),'a date must be specified to continue');
                cExperiment.metadata.date = date_match{1};
            end
            assert(ischar(cExperiment.metadata.date) && isrow(cExperiment.metadata.date),...
                'a date char vector must be specified to continue');
            assert(~isempty(regexp(cExperiment.metadata.date,'\d{4}-\d{2}-\d{2}','once')),...
                'a valid date must be specified to continue');
            
        end
    end
    
    methods (Access={?babyExperiment,?babyExperimentOmero,...
            ?OmeroDatabase,?babyExperimentSamples,...
            ?babyExperimentBioformats,?babyExperimentFileSeries,...
            ?babyExperimentOmero})
        function propNames = copyprops(cExperiment,TemplateExperiment,omit)
            %COPYPROPS Copy all properties from a cExperiment into this one
            %   This function can copy both public and private properties.
            %   Use OMIT to specify a cellstr of properties that will not
            %   be copied. This function gets used in the loadobj method
            %   and also by the convertSegmented method of the 
            %   OmeroDatabase class.
            
            if nargin<3 || isempty(omit), omit = {}; end
            if ~iscellstr(omit)
                error('The "omit" argument must be a cellstr.');
            end
            
            copied = copyprops@babyExperiment(cExperiment,TemplateExperiment,omit);
            
            % Only populate copyable fields occuring in both this object
            % and the template object:
            propNames = intersect(...
                getCopyableProperties(cExperiment,'babyExperimentBioformats'),...
                getCopyableProperties(TemplateExperiment,'babyExperimentBioformats'));
            % Omit requested properties:
            propNames = setdiff(propNames,omit);
            
            % Omit properties copied by parent:
            copyNames = setdiff(propNames,copied);
            
            % Copy all properties/fields to this cExperiment:
            for f = 1:numel(copyNames)
                cExperiment.(copyNames{f}) = TemplateExperiment.(copyNames{f});
            end
        end
    end
    
    methods (Static)
        function cExperiment = loadobj(load_structure)
            % Store the class name in the load_structure before passing to
            % parent load routine
            if isstruct(load_structure)
                load_structure.load_class_name = 'babyExperimentBioformats';
            end
            cExperiment = loadobj@babyExperiment(load_structure);
        end
    end
end
