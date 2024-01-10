classdef (ConstructOnLoad) TrapUpdate < event.EventData
    %loggingEvents.TrapUpdate Carries data for babyLogging when 
    %changing traps
    %   This class stores the number of the trap for events triggered when
    %   changing traps.
    
    properties
        trap = 0;
    end
    
    methods
        function eventData = TrapUpdate(trap)
            eventData.trap = trap;
        end
    end
    
end
