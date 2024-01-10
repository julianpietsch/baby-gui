classdef MemAware < handle
    methods
        function [memusage,membytes] = getMemUsage(this)
            props = getCopyableProperties(this,class(this));
            totSize = 0;
            
            for ii=1:length(props)
                currentProperty = this.(props{ii});
                s = whos('currentProperty');
                totSize = totSize + s.bytes;
            end
            
            factors = struct('GB',1024^3,'MB',1024^2,'KB',1024,'bytes',1);
            fn = fieldnames(factors);
            f = find(cellfun(@(x) totSize/x>1,struct2cell(factors)),1);
            m = sprintf('%.3f %s', totSize/factors.(fn{f}), fn{f});
            if nargout<1, fprintf('%s\n',m);
            else, memusage = m; membytes = totSize; end
        end
    end
end