function setChannelOffset(cExperiment,positionsToSet,offset)
% setChannelOffset(cExperiment,positionsToSet,offset)
%
% sets the offset field for all positions speciified by the
% positionsToSet.(defaults to all positions). If offset is blank it is set
% by GUI.
%
% This GUI could really do with improvement if this is to be used
% regularly.
%
% See also TIMELAPSETRAPS.RETURNSINGLETIMEPOINT
if nargin<2
    positionsToSet=1:length(cExperiment.dirs);
end

if nargin<3
    cTimelapse = cExperiment.loadCurrentTimelapse(positionsToSet(1));

    num_lines=1;
    dlg_title = 'ChannelOffsets?';
    def=[];prompt=[];
    for i=1:size(cTimelapse.offset,1)
        def{i} = num2str(cTimelapse.offset(i,:));
        prompt{i} = ['offset channel ' num2str(i) ' : ' cTimelapse.channelNames{i}];
    end
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    
    for i=1:size(cTimelapse.offset,1)
        offset(i,:)=str2num(answer{i});
    end
end

cExperiment.setTimelapsesProperty(positionsToSet,'offset',offset);

fprintf('finished setting offsets\n')