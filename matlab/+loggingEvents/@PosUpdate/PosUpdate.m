classdef (ConstructOnLoad) PosUpdate < event.EventData
    %loggingEvents.index Carries data for babyLogging when changing 
    %position
    %   This class stores an index and label for events triggered when
    %   loading positions or loading timepoints.
    
    properties
        index = 0;
        label = '';
        cTimelapse = {};
    end
    
    methods
        function eventData = PosUpdate(index,label,cTimelapse)
            eventData.index = index;
            eventData.label = label;
            eventData.cTimelapse = cTimelapse;
        end
    end
    
end

