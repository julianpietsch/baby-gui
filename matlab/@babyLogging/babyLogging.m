classdef babyLogging<handle
    %babyLogging Class to track progress of segmentation
    %   Use this class to update a log file and progress bars when
    %   iterating through time-consuming procedures
    
    properties
        cExperiment
        
        start_time = [] % Also used to identify if a protocol is running
        protocol_name = 'creating experiment' % default first protocol name
        protocol_args = ''
        position = [] % Index of current position, empty if not initialised
        npos = [] % Total number of positions to iterate over
        timepoint = [] % Index of current timepoint, empty if not initialised
        ntimes = [] % Total number of time points to iterate over
        trap = [] % Index of current trap, empty if not initialised
        ntraps = [] % Total number of traps to iterate over
        
        window_closed = false
        cancel = false
    end
    
    properties (Dependent)
        file_name
        file_dir
        
        progress_bar
        shouldLog
        use_gui
        stop_on_error
    end
    
    properties (Transient,Access=private)
        progress_bar_obj = []
        shouldLog_val
        use_gui_val
        stop_on_error_val = true
        file_handle = []
        file_handle_dir = ''
    end
    
    properties (Transient)
        % Handles to all of the listeners are stored here
        listenPositionChanged
        listenTimepointChanged
        listenTrapChanged
        listenExptLogMsg
        listenPosLogMsg
    end
    
    methods
        function this = babyLogging(cExperiment,shouldLog)
            %Constructor Construct a babyLogging object
            %   Pass in an babyExperiment object to intialise
            if nargin<2 || isempty(shouldLog)
                shouldLog = true;
            end
            
            % Save a handle to the cExperiment, since it will be used to
            % update the log file name and find the loaded timeLapse:
            this.cExperiment = cExperiment;
            
            % Also start listening for events on the cExperiment:
            this.listenPositionChanged = ...
                addlistener(cExperiment,'PositionChanged',@this.update_pos);
            this.listenExptLogMsg = ...
                addlistener(cExperiment,'LogMsg',@this.log_message);
            
            % Initialise the log file properties
            this.shouldLog = shouldLog;
        end
        
        function val = get.file_name(this)
            if isa(this.cExperiment,'babyExperimentOmero')
                val = ['cExperiment_log_' this.cExperiment.rootFolder '.txt'];
            else
                val = 'cExperiment_log_.txt';
            end
        end
            
        function val = get.file_dir(this)
            % Always return the most up-to-date saveFolder
            val = this.cExperiment.saveFolder;
        end
        
        function val = get.progress_bar(this)
            % Only initialise the progress bar if it hasn't already been
            % initialised:
            if isempty(this.progress_bar_obj)
                % Initialise the progress bar class
                this.progress_bar_obj = Progress('0;'); % Dummy callback to create 'finalise' button
                this.progress_bar_obj.finalise_button.setLabel('Cancel...');
                button_handle = handle(this.progress_bar_obj.finalise_button,'callbackproperties');
                set(button_handle, 'MouseClickedCallback', {@this.cancel_protocol});
                % Centre the dialog box
                this.progress_bar_obj.frame.setLocationRelativeTo([]);
            end
            val = this.progress_bar_obj;
        end
        
        function val = get.shouldLog(this)
            if isempty(this.shouldLog_val)
                % Unless explicitly set to false, default is true
                this.shouldLog_val = true;
            end
            val = this.shouldLog_val;
        end
        function set.shouldLog(this,val)
            if isempty(val)
                % Unless explicitly set to false, default is true
                val = true;
            else
                val = logical(val(1));
            end
            this.shouldLog_val = val;
            if ~val
                % Logging has been turned off so safely delete the GUI 
                % components if they have been initialised and close the
                % file handle if it is still open:
                this.close_gui;
                this.close_logfile;
            end
        end
        
        function val = get.use_gui(this)
            if isempty(this.use_gui_val)
                this.use_gui_val = ~strcmp(...
                    java.lang.System.getProperty('java.awt.headless'),'true');
            end
            val = this.use_gui_val;
        end
        function set.use_gui(this,val)
            this.use_gui_val = val;
            if ~val
                % GUI has been turned off so safely delete the GUI 
                % components if they have been initialised:
                this.close_gui;
            end
        end
        
        function val = get.stop_on_error(this)
            val = this.stop_on_error_val;
        end
        function set.stop_on_error(this,val)
            this.stop_on_error = val;
        end
        
        function add_arg(this,name,value)
            if this.shouldLog
                value = this.flatten_struct('',value);
                this.protocol_args = [this.protocol_args,name,': ',...
                    strrep(value,'\','\\'),'\n'];
            end
        end
        
        function start_protocol(this,name,npos)
            % Call this function before a protocol starts that requires
            % logging. Optionally supply the number of positions that will
            % be iterated over.
            if this.shouldLog
                if ~isempty(this.start_time)
                    % Another protocol may have terminated prematurely, so
                    % force a reset:
                    this.reset;
                end
                if nargin<3
                    this.npos = length(this.cExperiment.dirs);
                else
                    this.npos = npos;
                end
                this.start_time = datetime('now');
                this.protocol_name = name;
                if this.use_gui
                    this.progress_bar.frame.setTitle(name);
                end
                % The following should automatically open the log file if it
                % isn't already open:
                if isempty(this.protocol_args)
                    this.append(['\n\n=====================\n',...
                        datestr(now),'\tStart ',name,'...\n']);
                else
                    this.append(['\n\n=====================\n',...
                        datestr(now),'\tStart ',name,' using parameters:\n',...
                        this.protocol_args]);
                end
            end
        end
        
        function timestring = time_taken(this)
            % Return a string reporting the time since the protocol was
            % started
            if this.shouldLog
                if ~isempty(this.start_time)
                    runtime = diff([this.start_time;datetime('now')]);
                    if runtime<minutes(1)
                        timestring = sprintf('%0.0f secs',seconds(runtime));
                    elseif runtime<hours(1)
                        timestring = sprintf('%0.1f mins',minutes(runtime));
                    else
                        timestring = sprintf('%0.1f hours',hours(runtime));
                    end
                else
                    timestring = '0 secs';
                end
            end
        end
        
        function complete_protocol(this)
            % Call this function after a protocol that has been logging has
            % now finished. Can be safely called even if protocol was never
            % started.
            if this.shouldLog
                if ~isempty(this.start_time)
                    this.append([datestr(now),'\tSuccessfully completed ',...
                        this.protocol_name,' in ',this.time_taken,...
                        '.\n---------------------\n']);
                    if isa(this.cExperiment,'babyExperimentOmero')
                        this.cExperiment.uploadLogFile;
                    end
                end
                this.reset;
            end
        end
        
        function protocol_error(this)
            % Call this function whenever an error has been thrown
            % mid-protocol. Logging terminates with an error message and
            % the state is reset.
            if this.shouldLog
                if this.cancel
                    if ~isempty(this.start_time)
                        this.append([datestr(now),'\t',...
                            this.protocol_name,' was cancelled by user after ',...
                            this.time_taken,'.\n---------------------\n']);
                    end
                else
                    if ~isempty(this.start_time)
                        this.append([datestr(now),'\t',...
                            this.protocol_name,' terminated with an error after ',...
                            this.time_taken,'.\n---------------------\n']);
                    end
                    
                    if this.use_gui
                        e = errordlg(sprintf('Oh no! There was an error when %s.',...
                            this.protocol_name));
                        if this.stop_on_error
                            uiwait(e);
                        end
                    end
                end
                
                this.reset;
            end
        end
        
        function reset(this)
            % Reset the state of the logger
            if this.shouldLog
                % Pop bars off the progress bar until it is closed
                if this.use_gui
                    while this.progress_bar.bar_count > 0
                        this.progress_bar.pop_bar;
                    end
                    % Note that the above automatically removes the ticker 
                    % when the bar count reaches 0.
                
                    % Just in case the user closed the window manually in
                    % the previous run:
                    assignin('base', ['prog_terminate',...
                        num2str(this.progress_bar.window_number)], false);
                    this.window_closed = false;
                end
                
                % Reset the cancel flag in case it got set
                this.cancel = false;
                
                % Close the log file for now
                this.close_logfile;
                
                this.position = [];
                this.npos = [];
                this.timepoint = [];
                this.ntimes = [];
                this.trap = [];
                this.ntraps = [];
                this.protocol_args = '';
                
                % Reset start_time last in case there are any errors with the above
                this.start_time = [];
            end
        end
        
        function cancel_protocol(this,~,~)
            if this.shouldLog
                this.cancel = true;
            end
        end
        
        function update_progress(this,val)
            % Write a wrapper function to update the progress bar and
            % handle the case when the progress window has been closed.
            % Unlike the default behaviour of the Progress class, here we
            % just want window closure to be silent.
            
            if this.shouldLog && this.use_gui
                % Only call the following if the window has not been closed
                % already:
                if ~this.window_closed
                    %Let the swing thread assign the terminate flag if necessary
                    drawnow;
                    %Check for terminate flag
                    if (evalin('base', ['prog_terminate' num2str(this.progress_bar.window_number)]))
                        % The listener and ticker are guaranteed to be cleaned
                        % up by this.reset before the next protocol is run.
                        this.window_closed = true;
                    else
                        this.progress_bar.set_val(val);
                    end
                end
            end
        end
        
        function update_pos(this,~,posUpdateEvent)
            % Callback function for the PositionChanged event
            
            if this.shouldLog
                % Add this babyLogging instance to the cTimelapse:
                posUpdateEvent.cTimelapse.logger = this;
                
                % Start listening for events on the babyTimelapse object:
                this.listenPosLogMsg = ...
                    addlistener(posUpdateEvent.cTimelapse,'LogMsg',@this.log_message);
                this.listenTimepointChanged = ...
                    addlistener(posUpdateEvent.cTimelapse,'TimepointChanged',@this.update_timepoint);
                this.listenTrapChanged = ...
                    addlistener(posUpdateEvent.cTimelapse,'TrapChanged',@this.update_trap);
                
                % Only do anything else if a protocol is running
                if ~isempty(this.start_time)
                    
                    % Update the number of timepoints that we might expect to
                    % process:
                    this.ntimes = length(posUpdateEvent.cTimelapse.timepointsToProcess);
                    % Update the number of traps that we expect to process
                    this.ntraps = 0;
                    if ~isempty(posUpdateEvent.cTimelapse.cTimepoint) ...
                            && isfield(posUpdateEvent.cTimelapse.cTimepoint,'trapInfo')
                        tp = 1;
                        if ~isempty(posUpdateEvent.cTimelapse.timepointsToProcess)
                            tp = posUpdateEvent.cTimelapse.timepointsToProcess(1);
                        end
                        this.ntraps = numel(posUpdateEvent.cTimelapse.cTimepoint(tp).trapInfo);
                    end

                    % Update progress bar:
                    if isempty(this.position)
                        this.position = 1;
                        if this.use_gui
                            this.progress_bar.push_bar('Position',1,this.npos);
                        end
                    else
                        if this.position > 0 && this.position < this.npos
                            % If a timepoint was set, then we should start 
                            % a new line in the command window and pop
                            % a bar off the stack:
                            if ~isempty(this.timepoint)
                                if this.use_gui && this.progress_bar.bar_count > 1
                                    this.progress_bar.pop_bar;
                                end
                                fprintf('\n'); % New line for command window dot tracking
                                this.timepoint = []; % Reset the timepoint tracker
                                this.trap = []; % Reset the trap tracker
                            end
                            this.position = this.position+1;
                            
                            this.update_progress(this.position);
                            
                        else
                            % If the number of positions supplied to
                            % start_protocol was correct, then this should
                            % never be called, but to bullet proof:
                            
                            % Looping of the positions should have finished, so
                            % pop the necessary progress bars:
                            if ~isempty(this.timepoint)
                                if this.use_gui && this.progress_bar.bar_count > 1
                                    this.progress_bar.pop_bar;
                                end
                                fprintf('\n'); % New line for command window dot tracking
                                this.timepoint = []; % Reset the timepoint tracker
                                this.trap = []; % Reset the trap tracker
                            end
                            
                            if this.use_gui && this.progress_bar.bar_count > 1
                                this.progress_bar.pop_bar;
                            end
                        end
                    end
                    
                    % Append a message to the log file
                    this.append(sprintf('%s\tProcessing position %i (%s)\n',...
                        datestr(now),posUpdateEvent.index,posUpdateEvent.label));
                end
            end
        end
        
        function update_timepoint(this,~,~)
            % Callback function for the TimepointChanged event
            
            if this.shouldLog
                % Only do anything if a protocol is running
                if ~isempty(this.start_time)
                    
                    % Update progress bar:
                    if isempty(this.timepoint)
                        this.timepoint = 1;
                        if this.use_gui
                            this.progress_bar.push_bar('Time point',1,this.ntimes);
                        end
                        fprintf('.');
                    else
                        if this.timepoint > 0 && this.timepoint < this.ntimes
                            this.timepoint = this.timepoint+1;
                            this.update_progress(this.timepoint);
                            fprintf('.');
                            % Break line every 60 timepoints
                            if mod(this.timepoint,60)==0
                                fprintf('\n');
                            end
                        else
                            % If this.ntimes was correct, then this should never
                            % be called, but just in case...
                            if this.use_gui
                                this.progress_bar.pop_bar;
                            end
                            this.timepoint = [];
                        end
                    end
                    
                    % Do not log time points to the log file.
                end
            end
        end
        
        function update_trap(this,~,~)
            % Callback function for the TrapChanged event
            
            if this.shouldLog
                % Only do anything if a protocol is running
                if ~isempty(this.start_time)
                    
                    % Update progress bar:
                    if isempty(this.trap)
                        this.trap = 1;
                        if this.use_gui
                            this.progress_bar.push_bar('Trap',1,this.ntraps);
                        end
                        fprintf('.');
                    elseif this.trap > 0 && this.trap < this.ntraps
                        this.trap = this.trap+1;
                        this.update_progress(this.trap);
                        fprintf('.');
                        % Break line every 60 traps
                        if mod(this.trap,60)==0
                            fprintf('\n');
                        end
                    else
                        % If this.ntraps was correct, then this should never
                        % be called, but just in case...
                        if this.use_gui
                            this.progress_bar.pop_bar;
                        end
                        this.trap = [];
                    end
                    
                    % Do not log traps to the log file.
                end
            end
        end
        
        function log_message(this,~,msgEvent)
            % Callback function for the LogMsg event. Just add the message
            % to the log file/command window.
            if this.shouldLog
                this.append([datestr(now),'\t',strrep(msgEvent.message,'\','\\'),'\n']);
            end
        end
        
        function append(this,msg)
            % Ensure that the file is open for appending
            if this.shouldLog
                this.open_logfile;
                fprintf(this.file_handle,msg);
                fprintf(msg); % Also output to command window for backward compatibility
            end
        end
        
        function open_logfile(this)
            if this.shouldLog
                if isempty(this.file_handle)
                    % Open the file and record the handle dir
                    this.file_handle_dir = this.file_dir;
                    this.file_handle = ...
                        fopen(fullfile(this.file_handle_dir,this.file_name),'at');
                elseif ~strcmp(this.file_handle_dir,this.file_dir)
                    % Close the old file handle
                    this.close_logfile;
                    % Update file_handle_dir since it has changed (and also
                    % ensure that we never enter an infinite loop...)
                    this.file_handle_dir = this.file_dir;
                    % Try opening the log file again
                    this.open_logfile;
                else
                    % Otherwise do nothing, file should already be open
                    return
                end
            end
        end
        
        function close_logfile(this)
            % Close the log file handle only if it is not already closed:
            if ~isempty(this.file_handle)
                try
                    fclose(this.file_handle);
                catch
                    warning('Could not close log file in babyLogging class.');
                end
            end
            this.file_handle = [];
        end
        
        function close_gui(this)
            % Close the GUI only if it is not already closed:
            if ~isempty(this.progress_bar_obj)
                try
                    assignin('base', ['prog_terminate',...
                        num2str(this.progress_bar_obj.window_number)], false);
                    while this.progress_bar_obj.bar_count > 0
                        this.progress_bar_obj.pop_bar;
                    end
                    this.progress_bar_obj.frame.dispose;
                    this.window_closed = true;
                catch
                    warning('The babyLogging progress_bar could not be closed properly.');
                end
            end
            this.progress_bar_obj = [];
        end
        
        function delete(this)
            % Close the GUI and file handle.
            % Note that any errors are suppressed in delete functions.
            if ~isempty(this.progress_bar_obj)
                this.progress_bar_obj.frame.dispose;
            end
            if ~isempty(this.file_handle)
                fclose(this.file_handle);
            end
        end
        
        function obj = saveobj(~)
            if this.shouldLog
                % The contents of this class should not be saved
                warning('The babyLogging class should not be saved since it cannot be reloaded properly. Please only instantiate using the constructor.');
                obj = {};
            end
        end
    end
    
    methods (Static)
        function changePos(cExperiment,posIndex,cTimelapse)
            % First check (if we can) whether we should cancel this protocol
            if isa(cExperiment.logger,'babyLogging')
                if cExperiment.logger.cancel
                    error('protocol cancelled by user');
                end
            end
            notify(cExperiment,'PositionChanged',...
                loggingEvents.PosUpdate(posIndex,cExperiment.dirs{posIndex},cTimelapse));
        end
        
        function changeTimepoint(cTimelapse,timepoint)
            % First check (if we can) whether we should cancel this protocol
            if isa(cTimelapse.logger,'babyLogging')
                if cTimelapse.logger.cancel
                    error('protocol cancelled by user');
                end
            end
            notify(cTimelapse,'TimepointChanged',loggingEvents.TimepointUpdate(timepoint));
        end
        
        function changeTrap(cTimelapse,trap)
            % First check (if we can) whether we should cancel this protocol
            if isa(cTimelapse.logger,'babyLogging')
                if cTimelapse.logger.cancel
                    error('protocol cancelled by user');
                end
            end
            notify(cTimelapse,'TrapChanged',loggingEvents.TrapUpdate(trap));
        end
        
        function textout = flatten_struct(textin,value,depth)
            % Recursive function that converts almost any basic type to a
            % formatted string:
            
            if nargin<3
                depth = 0;
            end
            
            if ischar(value)
                textout = [textin,value];
                return
            end
            
            if isnumeric(value)
                if any(size(value)>5)
                    textout = [textin,'<numeric array>'];
                else
                    textout = [textin,num2str(value(:)')];
                end
                return
            end
            
            if iscell(value)
                textout = textin; % just in case the cell has zero length
                if length(value)>0
                    textout = babyLogging.flatten_struct(textout,value{1},depth);
                end
                
                if length(value)>1
                    for i=2:length(value)
                        textout = babyLogging.flatten_struct(...
                            [textout,', '],value{i},depth);
                    end
                end
                return
            end
            
            if isstruct(value)
                % Initialise textout in case there are no fields:
                textout = [textin,'{\n'];
                vnames = fieldnames(value);
                for i=1:length(vnames)
                    textout = [babyLogging.flatten_struct(...
                        [textout,repmat(' ',1,2*(depth+1)),vnames{i},': '],...
                        value.(vnames{i}),depth+1),'\n'];
                end
                textout = [textout,repmat(' ',1,2*depth),'}'];
                return
            end
            
            if isa(value,'function_handle')
                textout = [textin,func2str(value)];
                return
            end
            
            if islogical(value)
                if length(value)==1
                    if value
                        textout = [textin,'yes'];
                    else
                        textout = [textin,'no'];
                    end
                else
                    textout = [textin,num2str(value)];
                end
                return
            end
            
            % Default operation:
            textout = [textin,'<object of class "',class(value),'">'];
        end
        
        function this = loadobj(obj)
            if obj.shouldLog
                % This object cannot be loaded. Throw a warning.
                warning('The babyLogging class cannot be loaded. Please instantiate using the constructor instead.');
                this = obj;
            end
        end
    end
    
end
