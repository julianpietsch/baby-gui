classdef (ConstructOnLoad) TimepointUpdate < event.EventData
    %loggingEvents.TimepointUpdate Carries data for babyLogging when 
    %changing timepoints
    %   This class stores the number of the timepoint for events triggered
    %   when changing timepoints.
    
    properties
        timepoint = 0;
    end
    
    methods
        function eventData = TimepointUpdate(timepoint)
            eventData.timepoint = timepoint;
        end
    end
    
end

