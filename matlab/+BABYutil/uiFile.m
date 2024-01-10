classdef uiFile < handle
    %uiFile An auto-wrapping file selection and display control
    
    properties
        relDir
        maxChar
        filter
    end
    
    properties (Dependent)
        fileName
        string
    end
    
    properties (Transient)
        label
        btn
    end
    
    properties (Access=private)
        file
    end
    
    events
        FileChanged
    end
    
    methods
        function this = uiFile(parent,filter,filename,relDir,edit,maxChar)
            %uiFile(parent,filter,filename,edit,relDir,maxChar)
            %Create a file/folder selection control.
            %   filter: wildcard expression to filter files, or '/' if you want
            %   the user to select a folder
            if nargin<4
                relDir = '';
            end
            if nargin<5
                edit = true;
            end
            if nargin<6
                if edit
                    maxChar = 50;
                else
                    maxChar = 60;
                end
            end
            
            this.relDir = relDir;
            this.maxChar = maxChar;
            this.filter = filter;
            
            hbox = uix.HBox('Parent',parent,'Padding',0,'Spacing',3);
            this.label = uicontrol('Parent',hbox,'Style','text',...
                'HorizontalAlignment','left');
            
            % If the filter is an empty string, do not perform
            % auto-wrapping:
            if isempty(filter)
                this.string = filename;
            else
                this.fileName = filename;
            end
            
            if edit
                this.btn = uicontrol('Parent',hbox,'Style','pushbutton',...
                    'String','Change...','Callback',@this.changeFile);
                set(hbox,'Widths',[-1,60]);
            else
                this.btn = [];
            end
        end
        
        function out = wrapFileName(this,in)
            % First check if the file is within the base folder
            if ~isempty(in)
                rel = relativepath(in,this.relDir);
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
        
        function changeFile(this,~,~)
            if this.filter=='/'
                this.fileName = uigetdir(this.relDir);
            else
                curDir = pwd;
                cd(this.relDir);
                [newFile,newPath] = uigetfile(this.filter);
                cd(curDir);
                if isequal(newFile,0) || isequal(newPath,0)
                    return
                end
                this.fileName = fullfile(newPath,newFile);
            end
            this.label.String = this.wrapFileName(this.fileName);
        end
        
        function set.fileName(this,newFile)
            this.file = newFile;
            this.label.String = this.wrapFileName(this.fileName);
            notify(this,'FileChanged');
        end
        
        function field = get.fileName(this)
            field = this.file;
        end
        
        function set.string(this,newString)
            this.file = '';
            this.label.String = newString;
        end
        
        function field = get.string(this)
            field = this.label.String;
        end
    end
end