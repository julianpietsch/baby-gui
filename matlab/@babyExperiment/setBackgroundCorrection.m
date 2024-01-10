function setBackgroundCorrection(cExperiment,BackgroundCorrection,channel,positionsToSet)
%setBackgroundCorrection(cExperiment,BackgroundCorrection,channel,positionsToSet)
%
% set the flat field correction for each cTimelapse specified by
% positionsToSet (defaults to all positions)
%
% BackgroundCorrection  :   Flat field correction: an image of the size of
%                           the images extracted.Image are dot multipled by
%                           this image before being returned by
%                           returnTimepoint.
% channel               :   array of channel indices which should have this
%                           multiplication applied.
%
% See also TIMELAPSETRAPS.RETURNSINGLETIMEPOINT

   
if nargin<4 || isempty(positionsToSet)
    positionsToSet=1:length(cExperiment.dirs);
end


for i=1:length(positionsToSet)
    currentPos=positionsToSet(i);
    cTimelapse = cExperiment.loadCurrentTimelapse(currentPos);
    
    % if channel was not provided, request by gui.
    if i==1 && (nargin<3 || isempty(channel))
        [channel,ok] = listdlg('ListString',cTimelapse.channelNames,...
            'SelectionMode','single',...
            'Name','channel to correct',...
            'PromptString','Please select the channel to which to apply the background correction');

        if ~ok
            return
        end
    end
    
    cTimelapse.BackgroundCorrection{channel} = BackgroundCorrection;
    cExperiment.saveTimelapseExperiment;
end
