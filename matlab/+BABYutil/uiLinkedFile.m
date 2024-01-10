classdef uiLinkedFile < handle
    %uiLinkedFile An auto-wrapping file selection and display control
    
    properties (Dependent)
        link
        filter
        basedir
        Enable
        maxChar
        emptyMessage
    end
    
    properties (Constant)
        % Error message types that get raised by this class
        Errors = struct(...
            'BadParams','BABY:uiLinkedFile:badParams')
        notFoundMessage = '<File No Longer Exists>'
    end
    
    properties (Access=private)
        link_val
        uiLabel
        uiBtn
        uiBox
        props = struct('filter','*.*','basedir','',...
            'Enable','on','maxChar',55,...
            'emptyMessage','<Not yet specified...>');
        prop_names_abbr
    end
    
    events
        UserUpdate % Notify parent that the user has updated the link
    end
    
    methods
        function this = uiLinkedFile(parent,link,varargin)
            %uiLinkedFile GUI control for display and editing of file links
            %   parent: parent object in GUI
            %   link: LinkedFile object to start with (optional)
            %   Follow with named arguments:
            %   - 'filter': wildcard expression to filter files; specify the 
            %       empty string '' to disable editing or '/' if you want a
            %       folder selection dialog
            %   - 'basedir': a directory to determine relative paths from
            %       and to start the file selection dialog in
            %   - 'Enable': whether or not the button is active
            %   - 'maxChar': maximum width of label in characters (it will
            %       auto-wrap at the folder separators)
            %   - 'emptyMessage': message that gets displayed when link is
            %       empty
            
            custom_args = {};
            this.link_val = [];
            if nargin>1
                if isa(link,'LinkedFile') || isempty(link)
                    this.link_val = link;
                    custom_args = varargin;
                else
                    custom_args = [{link},varargin];
                end
            end
            
            % Parse the argument list and overwrite defaults
            nargs = length(custom_args)/2;
            assert(floor(nargs)==nargs,this.Errors.BadParams,...
                'Named arguments must be in pairs.');
            arg_names = custom_args(1:2:(2*nargs));
            assert(iscellstr(arg_names),this.Errors.BadParams,...
                'Some arguments are not named.');
            custom_args = struct(custom_args{:});
            
            % Pre-calculate abbreviated property names:
            this.prop_names_abbr = cellfun(@(f) lower(f(1:3)),...
                fieldnames(this.props),'UniformOutput',false);
            arg_names_abbr = cellfun(@(f) lower(f(1:3)),arg_names,...
                'UniformOutput',false);
            
            % Specify some sensible defaults if not already specified:
            if ~isempty(this.link) && ~ismember('bas',arg_names_abbr) && ...
                    this.link.exists
                % Use the specified LinkedFile to set default basedir:
                this.props.basedir = this.link.basedir;
            end
            if ~ismember('fil',arg_names_abbr)
                % This must be an argument to "this.set" to create the UI
                custom_args.filter = this.props.filter;
            end
            
            % Create a container for this UI
            this.uiBox = uix.HBox('Parent',parent,'Padding',0,'Spacing',3);
            
            % Update all of the properties and create controls
            this.set(custom_args);
        end
        
        function set(this,varargin)
            % Parse the argument list
            if nargin==2 && isstruct(varargin{1})
                new_props = varargin{1};
                arg_names = fieldnames(new_props);
            else                
                nargs = length(varargin)/2;
                assert(floor(nargs)==nargs,this.Errors.BadParams,...
                    'Named arguments must be in pairs.');
                arg_names = varargin(1:2:(2*nargs));
                assert(iscellstr(arg_names),this.Errors.BadParams,...
                    'Some arguments are not named.');
                new_props = struct(varargin{:});
            end
            
            % Simplify to abbreviated names
            arg_names_abbr = cellfun(@(f) lower(f(1:3)),arg_names,...
                'UniformOutput',false);
            
            % Check that all arguments have valid names
            valid_args = ismember(arg_names_abbr,this.prop_names_abbr);
            if ~all(valid_args)
                error(this.Errors.BadParams,...
                    '%s are not valid property specifiers',...
                    strjoin(strcat('"',arg_names(~valid_args),'"'),', '));
            end
            
            % Check that all arguments have valid values and update properties
            for a=1:length(arg_names)
                arg_val = new_props.(arg_names{a});
                switch arg_names_abbr{a}
                    case 'fil' % 'filter'
                        assert(ischar(arg_val),this.Errors.BadParams,...
                            '"filter" must be a character vector');
                        this.props.filter = arg_val;
                    case 'bas' % 'basedir'
                        assert(ischar(arg_val),this.Errors.BadParams,...
                            '"basedir" must be a character vector');
                        assert(isdir(arg_val),this.Errors.BadParams,...
                            '"basedir" must be a valid directory');
                        this.props.basedir = arg_val;
                    case 'ena' % 'Enable'
                        assert(ischar(arg_val) && ismember(arg_val,{'on','off'}),...
                            this.Errors.BadParams,...
                            '"Enable" must be either "on" or "off"');
                        this.props.Enable = arg_val;
                    case 'max' % 'maxChar'
                        assert(isnumeric(arg_val) && floor(arg_val)==arg_val,...
                            this.Errors.BadParams,...
                            '"maxChar" must be an integer');
                        this.props.maxChar = arg_val;
                    case 'emp' % 'emptyMessage'
                        assert(ischar(arg_val),this.Errors.BadParams,...
                            '"emptyMessage" must be a character vector');
                        this.props.emptyMessage = arg_val;
                end
            end
            
            % If a filter has been specified, recreate the UI elements
            if ismember('fil',arg_names_abbr)
                % First delete any existing UI elements
                delete(this.uiBox.Children);
                % Then create a label
                this.uiLabel = uicontrol('Parent',this.uiBox,...
                    'Style','text','HorizontalAlignment','left');
                % If the filter is non-empty, create a "Set" button
                if ~isempty(this.filter)
                    this.uiBtn = uicontrol('Parent',this.uiBox,...
                        'Style','pushbutton','String','Set',...
                        'Callback',@this.changeFile);
                    set(this.uiBox,'Widths',[-1,40]);
                else
                    this.uiBtn = [];
                    set(this.uiBox,'Widths',-1);
                end
            end
            
            % If a basedir has been specified, update the link
            if ismember('bas',arg_names_abbr) && ~isempty(this.link)
                if this.link.exists
                    rel = relativepath(this.link.filename,val);
                    this.link.updateLink(rel,val);
                else
                    try
                        this.link.updateBaseDir(this.basedir);
                    catch err
                        % Fail silently if no file could be found for this
                        % link, otherwise rethrow the error
                        if ~strcmp(err.identifier,BABYutil.LinkedFile.errMissingFile)
                            rethrow(err);
                        end
                    end
                end
            end
            
            % If Enable or filter have been specified, set button status
            if any(ismember({'ena','fil'},arg_names_abbr)) && ~isempty(this.uiBtn)
                this.uiBtn.Enable = this.Enable;
            end
            
            % Always update the label
            this.updateLabel;
        end
        
        function changeFile(this,~,~)
            if this.filter=='/'
                newDir = uigetdir(this.basedir);
                if isequal(newDir,0)
                    % User cancelled, so abort
                    return
                end
                rel = relativepath(newDir,this.basedir);
                if isempty(this.link)
                    this.link = BABYutil.LinkedFile(rel,this.basedir);
                else
                    this.link.updateLink(rel,this.basedir);
                end
            else
                curDir = pwd;
                cd(this.basedir);
                [newFile,newPath] = uigetfile(this.filter);
                cd(curDir);
                if isequal(newFile,0) || isequal(newPath,0)
                    % User cancelled, so abort
                    return
                end
                rel = relativepath(fullfile(newPath,newFile),this.basedir);
                if isempty(this.link)
                    this.link = BABYutil.LinkedFile(rel,this.basedir);
                else
                    this.link.updateLink(rel,this.basedir);
                end
            end
            this.updateLabel;
            notify(this,'UserUpdate');
        end
        
        function set.link(this,val)
            if isempty(val)
                this.link_val = [];
                this.updateLabel;
            elseif isa(val,'BABYutil.LinkedFile')
                this.link_val = val;
                this.updateLabel;
            else
                error(this.Errors.BadParams,...
                    'The link can only be set to empty or a LinkedFile object');
            end
        end
        function val = get.link(this)
            val = this.link_val;
        end
        
        function set.filter(this,val)
            this.set('filter',val);
        end
        function val = get.filter(this)
            val = this.props.filter;
        end
        
        function set.basedir(this,val)
            this.set('basedir',val);
        end
        function val = get.basedir(this)
            if isempty(this.props.basedir)
                val = cd;
            else
                val = this.props.basedir;
            end
        end
        
        function set.Enable(this,val)
            this.set('Enable',val);
        end
        function val = get.Enable(this)
            val = this.props.Enable;
        end
        
        function set.maxChar(this,val)
            this.set('maxChar',val);
        end
        function val = get.maxChar(this)
            val = this.props.maxChar;
        end
        
        function set.emptyMessage(this,val)
            this.set('emptyMessage',val);
        end
        function val = get.emptyMessage(this)
            val = this.props.emptyMessage;
        end
    end
    
    methods (Access=private)
        function updateLabel(this,~,~)
            if isempty(this.link)
                this.uiLabel.String = this.emptyMessage;
            elseif ~this.link.exists
                this.uiLabel.String = this.notFoundMessage;
            else
                this.uiLabel.String = this.wrapFileName(this.link.filename);
            end
        end
        
        function out = wrapFileName(this,in)
            % First check if the file is within the base folder
            if ~isempty(in)
                rel = relativepath(in,this.basedir);
                if isempty(regexp(rel,'\.\.','once'))
                    in = rel;
                end
            end
            % Then, iteratively split the file name into chunks no
            % longer than maxChar:
            out = '';
            while length(in)>this.maxChar
                lastsep = regexp(in(1:this.maxChar),'/[^/]*$');
                if isempty(lastsep)
                    out = [out,in(1:this.maxChar),sprintf('\n')];
                    in = in(this.maxChar+1:end);
                else
                    out = [out,in(1:lastsep),sprintf('\n')];
                    in = in(lastsep+1:end);
                end
            end
            out = [out,in];
        end
    end
end
