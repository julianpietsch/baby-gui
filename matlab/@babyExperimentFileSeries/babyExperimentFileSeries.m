classdef babyExperimentFileSeries < babyExperiment
    % BABYEXPERIMENTFILESERIES a subclass of babyExperiment to
    % manage cExperiments built from series of image files.
    %
    % See also BABYEXPERIMENT, BABYTIMELAPSEFILESERIES
    
    properties
        readerArgs
    end
    
    methods
        function cExperiment = babyExperimentFileSeries(rootfolder,savefolder,varargin)
            %BABYEXPERIMENTBIOFORMATS Create a new cExperiment
            %   - rootfolder (optional): specify a file(s) that can be loaded by 
            %       BioFormats.
            %   - savefolder (optional): specify a folder in which to save
            %       the cExperiment and cTimelapses.
            
            % call super class constructor such that it does not initialise
            % anything.
            cExperiment@babyExperiment(true);
            
            % If initialised with single input true, make an empty object
            % (used in load function)
            if nargin==1 && islogical(rootfolder) && rootfolder
                % if folder is true, cExperiment returned bare for load
                % function.
                return
            end
            
            ip = inputParser;
            ip.addRequired('rootfolder',@(x) ischar(x) && isrow(x));
            ip.addRequired('savefolder',@(x) ischar(x) && isrow(x));
            ip.addParameter('date','',@(x) ischar(x) && isrow(x));
            ip.KeepUnmatched = true;
            ip.parse(rootfolder,savefolder,varargin{:});
            
            cExperiment.readerArgs = [fieldnames(ip.Unmatched),struct2cell(ip.Unmatched)]';
            reader = ImageReaderFileSeries(rootfolder,cExperiment.readerArgs{:});
            cExperiment.rootFolder = rootfolder;
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
            
            cExperiment.dirs = reader.posNames;
            
            baseChannelNames = reader.meta.channelNames;
            maxnz = 1;
            for p=1:reader.npos, maxnz = max(maxnz,reader.imageSize(3)); end
            channelNames = cell(numel(baseChannelNames),1);
            for c=1:numel(channelNames)
                channelNames{c} = cell(maxnz,1);
                for z=1:maxnz
                    channelNames{c}{z} = sprintf('%s_%03u',baseChannelNames{c},z);
                end
            end
            cExperiment.channelNames = vertcat(baseChannelNames{:},channelNames{:});
            
            % Populate the metadata field of cExperiment
            % The following is the minimum detail required for a 
            % cExperiment with a single position:
            cExperiment.metadata = struct();
            cExperiment.metadata.acqFile = cExperiment.rootFolder; % For ID
            cExperiment.metadata.acq = struct(...
                'channels',struct('names',{cExperiment.channelNames},...
                'zsect',(maxnz>1)*ones(size(cExperiment.channelNames))),...
                'zsections',struct('sections',maxnz),'positions',zeros(0,1));
            
            if isfield(reader.meta,'posTimes')
                % Want times in minutes
                cExperiment.metadata.posTimes = reader.meta.posTimes/60;
            elseif isfield(reader.meta,'posTimeIntervals')
                posTimes = cell(reader.npos,1);
                for p=1:reader.npos
                    reader.pos = p;
                    posTimes{p} = reader.times;
                end
                maxtps = max(cellfun(@numel,posTimes));
                posTimes = cellfun(@(x) padarray(x,maxtps-numel(x),NaN),...
                    posTimes,'uni',0);
                cExperiment.metadata.posTimes = vertcat(posTimes{:});
            end
            
            if isempty(ip.Results.date)
                % See if we can determine a date from root folder name
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
            else
                cExperiment.metadata.date = ip.Results.date;
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
                getCopyableProperties(cExperiment,'babyExperimentFileSeries'),...
                getCopyableProperties(TemplateExperiment,'babyExperimentFileSeries'));
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
                load_structure.load_class_name = 'babyExperimentFileSeries';
            end
            cExperiment = loadobj@babyExperiment(load_structure);
        end
    end

    methods (Access=protected)
        function cTimelapse = newTimelapse(cExperiment,pos)
            cTimelapse = babyTimelapseFileSeries(...
                cExperiment.rootFolder,pos,...
                'channelNames',cExperiment.channelNames,...
                cExperiment.readerArgs{:});
        end
    end
end
