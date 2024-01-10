function keyPress_cb(cDisplay,src,event)
% keyPress_cb(cDisplay,src,event)
%
% Call back for the key press function:
% 
% h             :   (help) displays cTimelapseDisplay.gui_help in a
%                   helpdlg.
%
% in all other cases the key press is just stored like the general key
% press GUI


if strcmp(event.Key,'h')
    helpdlg(cDisplay.gui_help);
end


end

