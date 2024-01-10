function [channel,OK_response] = selectChannelGUI(cTimelapse,name,prompt_string,multiple_selection)
% [channel,Ok_response] = selectChannelGUI(cTimelapse,name,prompt_string,multiple_selection)
% simple channel select GUI.
% multiple selection is a logical of whether to allow selection of multiple
% channels.
% OK_response is 1 if Ok was pressed and zero otherwise


if nargin<4
    multiple_selection = false;
end

if multiple_selection
    select_string = 'multiple';
else
    select_string = 'single';
end

[channel,OK_response] = listdlg('ListString',cTimelapse.channelNames,'Name',name,'PromptString',prompt_string,'SelectionMode',select_string);


