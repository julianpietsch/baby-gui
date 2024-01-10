classdef (ConstructOnLoad) msg < event.EventData
    %loggingEvents.msg Carries data to inform babyLogging of a 
    %message to record to its log file.
    
    properties
        message = '';
    end
    
    methods
        function eventData = msg(message)
            eventData.message = message;
        end
    end
    
end

