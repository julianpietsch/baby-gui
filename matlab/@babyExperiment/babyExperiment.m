classdef babyExperiment < MemAware
    % class for organising and numerous babyTimelapse objects, one for each
    % position in the experiment. Mostly used to apply identical processing
    % steps to each position, organise loading and saving of them, and
    % compile the data from all the separate positions in one location.
    % Indiviual babyTimelapse objects are created from the folders in the
    % rootFolder and saved, with the babyExperiment object, in the
    % saveFolder.
    % Various subclasses handle alternative load/save/access methods,
    % including use of the OME Bioformats library, OMERO database and file
    % series defined by regular expressions.
    %
    % See also BABYCREATEGUI, BABYTIMELAPSE, BABYEXPERIMENTOMERO,
    % BABYEXPERIMENTBIOFORMATS, BABYEXPERIMENTFILESERIES
    
    properties
        rootFolder %folder where images are. When images are held in an Omero database (babyExperimentOmero subclass) this property is the suffix defining the filename: cExperiment_SUFFIX.mat
        creator %string, the user who created this object(obtained by getenv('USERNAME'))
        saveFolder %folder to save the timelapse objects
        dirs % cell array of directories in rootFolder. 
        posTrapsTracked % logical array of positions with traps tracked
        posSegmented % logical array of position already segmented
        posTracked % logical of positions with cells tracked
        cellsToPlot % Doesn't actually seem to be used anywhere
        metadata % structure of meta data filled by babyExperiment.parseLogFile
        marks % structure containing labelled timepoints for each position and trap
        
        cellAutoSelectParams
        cellInf % cell data compuled from extractedData in each of the individual timelapseTrap objects
        chamberMap = struct() % fields specify named groups of positions for final results
        shouldLog %a parameter that tells the logger whether it should do things 

        % The following all match their equivalents in babyTimelapse and
        % are populated when calling createTimelapsePositions
        pixelSize
        trapsPresent % whether or not images should be split into traps
        trapTemplates % templates and associated info for identifying traps; empty if there are no traps
        trapTemplateChannel % channel index for the trap template
        image_flipud
        image_rotation
        timepointsToProcess
        trackTrapsOverwrite
        channelNames %this has the list of the channel names  
    end
    
    properties (Transient)
        % Transient properties won't be saved
        cTimelapse % populated when loadCurrentTimelapse is used, and the cTimelapse saved when saveCurrentTimelapse is called.
        currentPos % populated when loadCurrentTimelapse is called. Defaultfor where to then save the timelapse.
    end
    
    properties (Dependent)
        id % A unique ID that links this experiment to Omero (filled by babyExperiment.parseLogFile); cannot be set
        logger % handle to an babyLogging object to keep a log
        imcache % handle to an ImageCache object for caching images in memory and on disk
        nzstacks % a Map specifying the number of z stacks for each channel
        posTimes % an nposes * ntimepoints array specifying times in minutes for each position
        babyBrain % options for and status of segmentation by the baby
    end
    
    properties (Access=protected)
        id_val = '' % This should never be updated if non-empty
        imcache_val
        babyBrain_val
    end
    
    properties (Transient,Access=private)
        logger_val
    end
    
    properties (Hidden=true)
        % these properties are not visible to the user
        clearOldTrapInfo % if this is true, when reselecting the taps through IdentifyTrapsTimelapses it will clear trapInfo first.
    end
    
    events
        PositionChanged
        LogMsg
    end
    
    methods
        
        function cExperiment=babyExperiment(rootFolder,saveFolder)
            %cExperiment=babyExperiment(rootFolder,saveFolder)
            % 
            % INPUTS 
            % rootFolder   -  EITHER :
            %                   - a boolean (true) indicating that a bare
            %                     babyExperiment object should be
            %                     created (used in loading and subclasses)
            %                 OR:
            %                   - a string with the full path to the root
            %                     folder (i.e. folder where all the
            %                     postition folders are). 
            %                 If empty, rootFolder is selected by user
            %                 input.
            % saveFolder   -  string. Full path to the folder where the
            %                 babyExperiment object and created
            %                 babyTimelapse objects should be saved.
            
            
            % Initialise source (root folder) and save paths
            if nargin<1
                fprintf('\n   Select the Root of a single experimental set containing folders of multiple positions \n');
                rootFolder=uigetdir(pwd,'Select the Root of a single experimental set containing folders of multiple positions');
            elseif islogical(rootFolder) && rootFolder
                % if folder is true, cExperiment returned bare for load
                % function.
                return
            end
            
            % Default to enabling logging:
            cExperiment.shouldLog=true;
            
            if nargin<2
                fprintf('\n   Select the folder where data should be saved \n');
                saveFolder=uigetdir(rootFolder,'Select the folder where data should be saved');
                if isempty(saveFolder)
                    fprintf('\n\n   No folder selected, no cExperiment created\n\n')
                    return
                end
            end            

            cExperiment.rootFolder=rootFolder;
            cExperiment.saveFolder=saveFolder;
            
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
            %Define the source folders (or Omero image names) for each
            %position

            tempdir=dir(cExperiment.rootFolder);

            cExperiment.dirs=cell(1);
            % created dirs - the list of positions - as the directories in
            % the rootFolder
                index=1;
                for i=1:length(tempdir)
                    if tempdir(i).isdir
                        if ~strcmp(tempdir(i).name(1),'.')
                            cExperiment.dirs{index}=tempdir(i).name;
                            index=index+1;
                        end

                    end
                end

            
            cExperiment.cellsToPlot=cell(1);
            
            
            %Parse the microscope acquisition metadata and attach the 
            %structure to the cExperiment object - this populates the 
            %metadata field of cExperiment. Only the meta data is collected
            %at this stage; the full log file can be parsed at extraction
            %since this can take an annoyingly long time with lots of
            %positions/timepoints...
            cExperiment.parseLogFile([],'meta_only');
            
        end
        
        function set.cTimelapse(cExperiment,cTimelapse)
            
            %though this is bad code, since both properties are transient
            %it shouldn't be a problem.
            if ~isequal(cExperiment.cTimelapse,cTimelapse)
                cExperiment.currentPos = [];
            end
            cExperiment.cTimelapse = cTimelapse;
        end
        
        function val = get.id(cExperiment)
            if isempty(cExperiment.id_val)
                % The ID will only be missing if the latest parseLogFile 
                % function has not yet been run on the experiment
                if isempty(cExperiment.metadata) || ~isfield(cExperiment.metadata,'acqFile')
                    fail_flag = cExperiment.parseLogFile([],'meta_only');
                    if fail_flag
                        val = '';
                    else
                        val = cExperiment.parseAcqFileIntoID(cExperiment.metadata.acqFile);
                        cExperiment.id_val = val;
                    end
                else
                    val = cExperiment.parseAcqFileIntoID(cExperiment.metadata.acqFile);
                    cExperiment.id_val = val;
                    %cExperiment.saveExperiment; % Save the new id into the cExperiment
                end
            else
                val = cExperiment.id_val;
            end
        end
        
        function val = get.logger(cExperiment)
            if isempty(cExperiment.logger_val)
                % Create a new logger to log changes for this cExperiment:
                cExperiment.logger_val = babyLogging(cExperiment);
            end
            % Always ensure that the 'shouldLog' state of the logger 
            % matches that of the cExperiment:
            cExperiment.logger_val.shouldLog = cExperiment.shouldLog;
            val = cExperiment.logger_val;
        end
        
        function val = get.imcache(cExperiment)
            val = cExperiment.imcache_val;
            if isscalar(val) && isa(val, 'ImageCache')
                % Ensure that the ImageCache always has an up-to-date
                % handle to this cExperiment
                val.cExperiment = cExperiment;
            else
                val = [];
            end
        end
        
        function set.imcache(cExperiment,val)
            assert(isa(val,'ImageCache') && isscalar(val),...
                'The "imcache" must be a valid ImageCache object');
            assert(strcmp(val.cExperiment.id,cExperiment.id),...
                'The "imcache" is for a different cExperiment');
            cExperiment.imcache_val = val;
        end
        
        function val = get.nzstacks(cExperiment)
            val = containers.Map;
            meta = cExperiment.metadata;
            if isempty(meta) || ~isfield(meta,'acq') || ...
                    ~all(isfield(meta.acq,{'channels',...
                    'zsections','positions'}))
                return;
            end
            
            acq = meta.acq;
            chnames = acq.channels.names;
            
            nz = acq.zsections.sections;
            for c=1:numel(chnames)
                chname = chnames{c};
                if acq.channels.zsect(strcmp(acq.channels.names,chname))
                    val(chname) = nz;
                else
                    val(chname) = 1;
                end
            end
        end
        
        function val = get.posTimes(cExperiment)
            val = [];
            if ~isstruct(cExperiment.metadata), return; end
            if isfield(cExperiment.metadata,'posTimes')
                val = cExperiment.metadata.posTimes;
            elseif isfield(cExperiment.metadata,'logTimes')
                val = cExperiment.metadata.logTimes;
                posNames = cExperiment.dirs;
                %If statement avoids an error when running single position experiments.
                if ~isempty(posNames) && isfield(cExperiment.metadata,'logPosNames') ...
                        && isequal(sort(posNames(:)),sort(cExperiment.metadata.logPosNames(:)))
                    posIndices = cellfun(@(s) find(strcmp(cExperiment.metadata.logPosNames,s),1),posNames,'Uniform',false);
                    posIndices=[posIndices{:}];
                    val = val(posIndices,:);
                end
            end
            if size(val,1)==1
                % Copy times for each position if only one reference given
                nposes = max(1,numel(cExperiment.dirs));
                val = val(ones(nposes,1),:);
            end
        end
        
        function val = get.babyBrain(cExperiment)
            if isempty(cExperiment.babyBrain_val)
                cExperiment.babyBrain_val = BabyBrain;
                cExperiment.babyBrain_val.set_channel_nzstacks(...
                    cExperiment.nzstacks);
            end
            val = cExperiment.babyBrain_val;
        end
        
        function set.babyBrain(cExperiment,val)
            assert(isa(val,'BabyBrain'),...
                '"babyBrain" must be specified as a "BabyBrain" object');
            % Copy only the configuration
            cExperiment.babyBrain.channelMap = val.channelMap;
            cExperiment.babyBrain.updateServerDetails(true);
            cExperiment.babyBrain.config = val.config;
        end
        
        function resetBrain(cExperiment)
            cExperiment.babyBrain_val = [];
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
            %   be copied. This function gets used in the loadobj method,
            %   by OmeroDatabase.convertSegmented and by 
            %   babyExperimentOmero.convertToFolderExperiment
            
            if nargin<3 || isempty(omit), omit = {}; end
            if ~iscellstr(omit)
                error('The "omit" argument must be a cellstr.');
            end
            
            % Only populate copyable fields occuring in both this object
            % and the template object:
            propNames = intersect(...
                getCopyableProperties(cExperiment,'babyExperiment'),...
                getCopyableProperties(TemplateExperiment,'babyExperiment'));
            % Omit requested properties
            propNames = setdiff(propNames,omit);
            
            % Copy all properties/fields to this cExperiment:
            for f = 1:numel(propNames)
                cExperiment.(propNames{f}) = TemplateExperiment.(propNames{f});
            end
        end
    end
    
    methods(Static)
        function val = parseAcqFileIntoID(acqfile)
            val = acqfile;
            
            re_root = '^.*AcquisitionData(?<mic>[^/\\]+)[/\\]'; % Microscope ID
            re_id = ['Swain Lab[/\\](?<user>[^/\\]+)[/\\]RAW DATA[/\\]',... % User ID
                '(?<year>\d+)[/\\](?<month>\w+)[/\\](?<day>\d+)',... % Date
                '[^/\\]*[/\\](?<name>.*)[/\\].*$']; % Experiment name
            
            % Handle old Batman experiments with a different root
            if isempty(regexp(acqfile,[re_root,re_id],'once')) ...
                    && ~isempty(regexp(acqfile,re_id,'once'))
                % Assume the microscope ID part is to blame
                exptdate = datenum(regexprep(acqfile,['^.*',re_id],...
                    '$<year>-$<month>-$<day>'),'yyyy-mmm-dd');
                if exptdate < datenum('2015-05-05','yyyy-mm-dd')
                    acqfile = regexprep(acqfile,...
                        '^(?<root>.*)(?<dirtree>[/\\]Swain Lab[/\\].*)$',...
                        'AcquisitionDataBatman$<dirtree>');
                end
            end
            
            if ~isempty(regexp(acqfile,[re_root,re_id], 'once'))
                val = regexprep(acqfile,[re_root,re_id], ...
                    '$<user>_$<mic>_$<year>_$<month>_$<day>_$<name>'); % Replacement string
            end
        end
        
        function cExperiment = loadobj(load_structure)
            % cExperiment = loadobj(load_structure)
            % load function to help maintain back compatability and take
            % care of fiddly loading behaviour
            
            if isstruct(load_structure)
                % The 'load_structure' argument could be for this class or any 
                % of its children. We work this out based on the 
                % 'load_class_name' property set by the children.
                if isfield(load_structure,'load_class_name')
                    lcname = load_structure.load_class_name;
                    assert(startsWith(lcname,'babyExperiment') ...
                        && exist(lcname,'class')==8,'bad load_class_name!');
                    cExperiment = eval([lcname,'(true)']);
                    load_structure = rmfield(load_structure,'load_class_name');
                else
                    % Default to instantiating according to this class
                    cExperiment = babyExperiment(true);
                end
            else
                cExperiment = eval([class(load_structure),'(true)']);
            end
            
            cExperiment.copyprops(load_structure);
            
            %% addtional stuff for back compatability etc.
            
            % Warn the user when loading if shouldLog is false
            if ~cExperiment.shouldLog
                warning('Logging is not active. To change, set "cExperiment.shouldLog = true"');
            end
            
            % back compatibility to put channel names into cExperiment
            % channelNames
            if isempty(cExperiment.channelNames) && ~isempty(cExperiment.dirs)
                cTimelapse = cExperiment.loadCurrentTimelapse(1);
                cExperiment.channelNames = cTimelapse.channelNames;
            end
            
            if isempty(cExperiment.clearOldTrapInfo)
                cExperiment.clearOldTrapInfo = false(size(cExperiment.dirs));
            end
            
            if isempty(cExperiment.posSegmented)
                cExperiment.posSegmented = false(size(cExperiment.dirs));
            end
            
            % This property was introduced after posSegmented. Here it is 
            % retrospectively assumed that if positions have been 
            % segmented, then the traps have also been tracked.
            if isempty(cExperiment.posTrapsTracked)
                cExperiment.posTrapsTracked = cExperiment.posSegmented;
            end
            
        end
                
        cExperiment = loadFrom(cExperiment_filepath,answer)
    end

    methods (Access=protected)
        function cTimelapse = newTimelapse(cExperiment,pos)
            cTimelapse = babyTimelapse(fullfile(cExperiment.rootFolder,cExperiment.dirs{pos}));
        end
    end
end

