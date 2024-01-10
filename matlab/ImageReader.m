classdef (Abstract) ImageReader < handle
    properties (Abstract, Dependent, SetAccess=private)
        npos
        posNames
        nchannels
        channels
        imageSize
        pixelSize
        timeInterval
        times
        meta
    end
    
    properties (Abstract, Dependent)
        pos
        posName
    end
    
    methods (Abstract)
        getTimepoint(this,T,Z,C)
    end
end