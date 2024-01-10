function removeCell(cTimelapse,timepoint,trap_index,cell_index)
% removeCell(cTimelapse,timepoint,trap_index,cell_index)
%
% cTimelapse        :   object of the babyTimelapse class (which will be
%                       edited by the removal of a cell).
% timepoint         :   the timepoint at which a cell should be removed
% trap_index        :   index of the trap from which a cell should be removed
% cell_index        :   indes of the cell to be removed (note, this is not
%                       the label - hich is used for tracking - but the
%                       index.
%
% removes a cell from the provided cTimelapse at the specified timepoint,
% trap_index,cell_index. If the cell is not there nothing will happen.

if ~isempty(cell_index) && cell_index <= length(cTimelapse.cTimepoint(timepoint).trapInfo(trap_index).cell)
    if length(cTimelapse.cTimepoint(timepoint).trapInfo(trap_index).cell)>1
        cTimelapse.cTimepoint(timepoint).trapInfo(trap_index).cell(cell_index)=[];
        if cell_index<=length(cTimelapse.cTimepoint(timepoint).trapInfo(trap_index).cellLabel)
            % if a cell was added by GUI, it will be added to the
            % end and not have a cellLabel
            % left for legacy reasons. cells should always have a
            % label now.
            cTimelapse.cTimepoint(timepoint).trapInfo(trap_index).cellLabel(cell_index)=[];
            
        end
    elseif length(cTimelapse.cTimepoint(timepoint).trapInfo(trap_index).cell)==1
        
        try
            cTimelapse.cTimepoint(timepoint).trapInfo(trap_index) = cTimelapse.trapInfoTemplate;
        catch
            % this slight strange call is because sometimes there are other
            % fields in trapInfo (mostly refinedTrapPixel), and if these
            % are present the above code will fail. This code below
            % acheives the same, preserving unknown fields in trapInfo and
            % setting any known ones to the default value.
            cTimelapse.cTimepoint(timepoint).trapInfo(trap_index) =...
                parse_struct(cTimelapse.trapInfoTemplate,...
                cTimelapse.cTimepoint(timepoint).trapInfo(trap_index));
        end
    end
    
end
end
