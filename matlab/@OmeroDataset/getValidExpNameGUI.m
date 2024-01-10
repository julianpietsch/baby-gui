function expName = getValidExpNameGUI(this)
%GETVALIDEXPNAME Prompt user with GUI to enter a cExperiment name
%   expName = obj.getValidExpName(dataset,existing) returns false if the
%   user cancelled the prompts, otherwise it will be a valid string
%   containing a unique name for a cExperiment that is also valid for use
%   in file names. The 'existing' argument is optional, but if provided,
%   should be a cellstr of experiment names that are already present for
%   this dataset (i.e., obtained using LISTCEXPERIMENTS).

existing = this.cExperiments;

% Determine an appropriate message string for the dialog box
msg = 'Enter a name for your new cExperiment. ';
if isempty(existing)
    msg = [msg,'No names have been taken yet.'];
else
    msg = [msg,sprintf('The existing names are:\n%s',...
        strjoin(existing,'\n'))];
end
digit_names = ~cellfun(@isempty,regexp(existing,'^[0-9]{3}$'));
if any(digit_names)
    dflt = sprintf('%03u',max(str2double(existing(digit_names)))+1);
else
    dflt = '001';
end

expName = false;

% Get the user to choose a name for this cExperiment
while true
    inputName = inputdlg({msg},'cExperiment name',1,{dflt});
    if isempty(inputName)
        % User cancelled, so abort creation of experiment
        return
    end
    inputName = inputName{1};
    if isempty(inputName)
        inputName = '001';
    end
    
    % Validate the inputName (it should look like a file name):
    if isempty(regexp(inputName,'^[_a-zA-Z0-9]+$','once'))
        dlg = errordlg('Invalid name. Valid names can only contain letters, numbers and underscores.');
        uiwait(dlg);
    elseif ismember(inputName,existing)
        response = questdlg(['An experiment with suffix "',inputName,...
            '" already exists. Do you want to replace it?'],...
            'Replace existing?','No');
        switch response
            case 'Yes'
                expName = inputName;
                return
            case 'No'
                % Give user the opportunity to specify an alternative
                continue
            otherwise
                % User cancelled, so abort creation of experiment
                return
        end

    else
        % We have a valid inputName
        expName = inputName;
        return
    end
end

end