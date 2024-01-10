function trapTimepoint=returnSingleTrapTimepoint(cTimelapse,trap_num_to_show,timepoint,channel,type)
%trapTimepoint=returnSingleTrapTimepoint(cTimelapse,trap_num_to_show,timepoint,channel)

if nargin<4 || isempty(channel)
    channel=1;
end

if nargin<5 || isempty(type)
    type=[];
end


trapTimepoint=cTimelapse.returnTrapsTimepoint(trap_num_to_show,timepoint,channel,type);

