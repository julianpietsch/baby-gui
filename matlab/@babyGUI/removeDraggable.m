function removeDraggable(this)
if ~isempty(this.dragPoints)
    for ch=1:length(this.dragPoints)
        if isempty(this.dragPoints{ch}), continue; end
        for p=1:length(this.dragPoints{ch})
            delete(this.dragPoints{ch}{p});
        end
    end
    this.dragPoints = {};
end
if ~isempty(this.dragCentre)
    for ch=1:length(this.dragCentre)
        delete(this.dragCentre{ch});
    end
    this.dragCentre = {};
end
if ~isempty(this.currentOutline)
    for ch=1:length(this.currentOutline)
        delete(this.currentOutline{ch});
    end
    this.currentOutline = {};
end
end