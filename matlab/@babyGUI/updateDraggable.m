function updateDraggable(this)

modenames = fieldnames(this.editModeMap);
outlinemodenames = fieldnames(this.outlineModeMap);
outlinemode = outlinemodenames{this.ctrls.outlinemode.Value};
if ~strcmp('editoutline',modenames{this.ctrls.editmode.Value}) || ...
        strcmp('clickmatch',outlinemode)
    % Ensure all draggable controls are deleted and exit
    this.removeDraggable;
    this.dragUpdate = false;
    return
end

this.dragAsRect = strcmp('dragrect',outlinemode);

%% Assign dimension variables
pad = this.dimensions.impad;
tW = this.tileWidth;
tH = this.tileHeight;
nT = this.tilesContext;
this.tp_pad = pad+(nT+this.dragtp)*(tW+pad);
nC = length(this.currentChannels);
channelInds = 1:nC;
this.yoffs = pad+(channelInds-1)*(tH+pad);

[zoffinit,zscinit] = this.getZoomActiveParams;

boxScale = 2;

[this.centre,this.radii,this.angles] = this.getContourParams;

%% Ensure that placeholder contour plot lines always exist
ax = this.frame.imAx;
nOutlines = length(this.currentOutline);
for ch=1:max(nC,nOutlines)
    if ch>nC
        delete(this.currentOutline{ch});
    else
        if invalid_ind(this.currentOutline,ch)
            this.currentOutline{ch} = plot(ax,[0,1],[1,1],'r','Visible','off');
        end
    end
end
this.currentOutline = this.currentOutline(cellfun(...
    @(x) ~isempty(x) & isvalid(x),this.currentOutline));

%% If there is no cell, then hide the draggable contour
if ~this.isDragVisible || ~this.ctrls.showoutlines.Value || this.dofocusannot
    if ~isempty(this.dragPoints)
        for ch=1:length(this.dragPoints)
            if isempty(this.dragPoints{ch}), continue; end
            for p=1:length(this.dragPoints{ch})
                delete(this.dragPoints{ch}{p});
            end
        end
    end
    if ~isempty(this.dragCentre)
        for ch=1:length(this.dragCentre)
            delete(this.dragCentre{ch});
        end
    end
    if ~isempty(this.currentOutline)
        for ch=1:length(this.currentOutline)
            this.currentOutline{ch}.Visible = 'off';
        end
    end
    return
end

% % Hide centres if adding buds
% if this.ctrls.addasbud.Value
%     if ~isempty(this.dragCentre)
%         for ch=1:length(this.dragCentre)
%             delete(this.dragCentre{ch});
%         end
%     end
% end

%% First update the contours

if this.ctrls.cartesianspline.Value
    [this.px,this.py] = BABYutil.cartesian_spline_from_radii(...
        double(this.radii),double(this.angles),double(this.centre),...
        [this.srcHeight,this.srcWidth]);
else
    [this.px,this.py] = BABYutil.get_full_points_from_radii(...
        double(this.radii),double(this.angles),double(this.centre),...
      [this.srcHeight,this.srcWidth]);
end
update_contours();
for ch=1:length(this.currentOutline)
    this.currentOutline{ch}.Visible = 'on';
end

%% Add and position draggable points if necessary
axconstraints = arrayfun(@(x) makeConstrainToRectFcn(...
    'impoint',this.tp_pad+[-tW,2*tW],x+[-tH,2*tH]),this.yoffs,'uni',0);
if this.dragAsRect
    [lc,lp] = rect_from_radii(this.radii,this.angles);
    this.cpoint = (this.centre + lc)*zscinit-zoffinit;
    this.points = (this.centre(ones(4,1),:) + lp)*zscinit-zoffinit(ones(4,1),:);
else
    this.cpoint = this.centre*zscinit-zoffinit;
    this.points = [...
        (this.centre(1)+this.radii.*cos(this.angles))*zscinit-zoffinit(1),...
        (this.centre(2)+this.radii.*sin(this.angles))*zscinit-zoffinit(2)];
end
this.dragUpdate = false;
for ch=channelInds
    if invalid_ind(this.dragCentre,ch)
        this.dragCentre{ch} = impoint(ax,...
            this.tp_pad+this.cpoint(1),this.yoffs(ch)+this.cpoint(2));
        setColor(this.dragCentre{ch},'g');
        addNewPositionCallback(this.dragCentre{ch},@(x) centre_move(ch,x));
    else
        this.dragCentre{ch}.setPosition([this.tp_pad,this.yoffs(ch)]+this.cpoint);
    end
    setPositionConstraintFcn(this.dragCentre{ch},axconstraints{ch});
    
    if ch>length(this.dragPoints) || isempty(this.dragPoints{ch})
        this.dragPoints{ch} = cell(size(this.points,1),1);
    end
    if length(this.dragPoints{ch})>size(this.points,1)
        for p=(size(this.points,1)+1):length(this.dragPoints{ch})
            delete(this.dragPoints{ch}{p});
        end
    end
    for p=1:size(this.points,1)
        if invalid_ind(this.dragPoints{ch},p)
            this.dragPoints{ch}{p} = impoint(ax,...
                this.tp_pad+this.points(p,1),this.yoffs(ch)+this.points(p,2));
            setColor(this.dragPoints{ch}{p},'b');
            addNewPositionCallback(this.dragPoints{ch}{p},@(x) point_move(ch,p,x));
        else
            this.dragPoints{ch}{p}.setPosition([this.tp_pad,this.yoffs(ch)]+this.points(p,:));
        end
        setPositionConstraintFcn(this.dragPoints{ch}{p},axconstraints{ch});
    end
end
this.dragUpdate = true;

%% Helper functions

    function bad = invalid_ind(h,ind)
        bad = ind>length(h) || isempty(h{ind}) || ~isvalid(h{ind});
    end

    function point_move(cInd,pInd,pos)
        if ~this.dragUpdate, return; end
        this.dragUpdate = false;
        for cch=1:length(this.currentChannels)
            if cch~=cInd
                this.dragPoints{cch}{pInd}.setPosition(...
                    pos+[0,this.yoffs(cch)-this.yoffs(cInd)]);
            end
        end
        this.dragUpdate = true;
        [zoff,zsc] = getZoomActiveParams(this);
        pos = pos-[this.tp_pad,this.yoffs(cInd)];
        realpos = (pos+zoff)/zsc;
        rp = ones(size(this.points,1),1);
        realpoints = (this.points+zoff(rp,:))/zsc;
        if this.dragAsRect
            rc = ones(size(this.radii,1),1);
            coords = this.centre(rc,:)+this.radii(:,[1,1]).*[...
                cos(this.angles(:)),sin(this.angles(:))];
            orthcentre = mean(coords);
            % Transform into PCA space
            [pcaRot,pcaCoords] = pca(coords);
            invPcaRot = inv(pcaRot)';
            pcaBox = (realpoints-orthcentre(rp,:))*invPcaRot;
            pcaCentre = (this.centre-orthcentre)*invPcaRot;
            % Calculate scaling
            fixInd = mod(2+pInd-1,4)+1; % hold the opposite ind fixed
            pcaOffset = pcaBox(fixInd,:);
            pcaPos = (realpos-orthcentre)*invPcaRot;
            scaling = (pcaPos-pcaOffset)./(pcaBox(pInd,:)-pcaOffset);
            this.centre = orthcentre + (pcaOffset + ...
                scaling.*(pcaCentre-pcaOffset))*pcaRot';
            coords = orthcentre(rc,:) + (pcaOffset(rc,:) ...
                + scaling(rc,:).*(pcaCoords-pcaOffset(rc,:)))*pcaRot';
            this.points = (orthcentre(rp,:) + (pcaOffset(rp,:) + ...
                scaling(rp,:).*(pcaBox-pcaOffset(rp,:)))*pcaRot')*zsc-zoff(rp,:);
            this.cpoint = (orthcentre+babyGUI.cpointScale ...
                *max(pcaCoords(:,2))*pcaRot(:,2)')*zsc-zoff;
            this.dragUpdate = false;
            for cch=1:length(this.currentChannels)
                for pJ=setdiff(1:size(this.points,1),[pInd,fixInd])
                    this.dragPoints{cch}{pJ}.setPosition(...
                        this.points(pJ,:)+[this.tp_pad,this.yoffs(cch)]);
                end
                this.dragCentre{cch}.setPosition(...
                    this.cpoint + [this.tp_pad,this.yoffs(cch)]);
            end
            this.dragUpdate = true;
            [r,a] = radii_from_points(coords,this.centre);
            this.radii = r; this.angles = a;
        else
            this.points(pInd,:) = pos;
            realpoints(pInd,:) = realpos;
            [r,a] = radii_from_points(realpoints,this.centre);
        end
        [this.px,this.py] = this.setContourParams(this.centre,r,a);
        update_contours();
    end

    function centre_move(cInd,pos)
        if ~this.dragUpdate, return; end
        this.dragUpdate = false;
        for cch=1:length(this.currentChannels)
            if cch~=cInd
                this.dragCentre{cch}.setPosition(...
                    pos+[0,this.yoffs(cch)-this.yoffs(cInd)]);
            end
        end
        this.dragUpdate = true;
        [zoff,zsc] = getZoomActiveParams(this);
        pos = pos-[this.tp_pad,this.yoffs(cInd)];
        realpos = (pos+zoff)/zsc;
        if this.dragAsRect
            if ismac
                doTranslation = ismember('alt',...
                    get(this.frame.fig,'CurrentModifier'));
            else
                doTranslation = ismember('control',...
                    get(this.frame.fig,'CurrentModifier'));
            end
            if doTranslation
                vmove = pos - this.cpoint;
                this.cpoint = pos;
                % Apply translation to rect
                this.points = this.points + vmove(ones(size(this.points,1),1),:);
                this.dragUpdate = false;
                for cch=1:length(this.currentChannels)
                    for pInd=1:size(this.points,1)
                        this.dragPoints{cch}{pInd}.setPosition(...
                            this.points(pInd,:)+[this.tp_pad,this.yoffs(cch)]);
                    end
                end
                this.dragUpdate = true;
                % Apply translation to contour
                this.centre = this.centre + vmove/zsc;
                r = this.radii;
                a = this.angles;
            else
                coords = this.centre(ones(length(this.radii),1),:) + ...
                    this.radii(:,[1,1]).*[...
                    cos(this.angles(:)),sin(this.angles(:))];
                
                % Determine a 2D rotation matrix using old and new pos:
                orthcentre = mean(coords)*zsc-zoff; % in zoomed coords
                vprev = this.cpoint-orthcentre;
                lprev = norm(vprev);
                vprev = vprev/lprev;
                this.cpoint = pos;
                vnew = this.cpoint-orthcentre;
                lnew = norm(vnew);
                vnew = vnew/lnew;
                scaling = lnew/lprev;
                costheta = dot(vprev,vnew);
                sintheta = vprev(1)*vnew(2)-vprev(2)*vnew(1); % from normed cross product with z=0
                rotMat = scaling*[costheta -sintheta; sintheta costheta]';
                % Apply rotation to rect
                orthcentre_rep = orthcentre(ones(size(this.points,1),1),:);
                this.points = orthcentre_rep+(this.points-orthcentre_rep)*rotMat;
                this.dragUpdate = false;
                for cch=1:length(this.currentChannels)
                    for pInd=1:size(this.points,1)
                        this.dragPoints{cch}{pInd}.setPosition(...
                            this.points(pInd,:)+[this.tp_pad,this.yoffs(cch)]);
                    end
                end
                this.dragUpdate = true;
                
                % Apply rotation to contour
                orthcentre = mean(coords); % in real coords
                this.centre = orthcentre+(this.centre-orthcentre)*rotMat;
                orthcentre_rep = orthcentre(ones(size(coords,1),1),:);
                coords = orthcentre_rep+(coords-orthcentre_rep)*rotMat;
                [r,a] = radii_from_points(coords,this.centre);
                this.radii = r; this.angles = a;
            end
        else
            this.cpoint = pos;
            this.centre = realpos;
            realpoints = (this.points+zoff(ones(size(this.points,1),1),:))/zsc;
            [r,a] = radii_from_points(realpoints,this.centre);
        end
        [this.px,this.py] = this.setContourParams(this.centre,r,a);
        update_contours();
    end

    function update_contours()
        [zoff,zsc] = getZoomActiveParams(this);
        cpx = [this.px(end);this.px]*zsc-zoff(1);
        cpy = [this.py(end);this.py]*zsc-zoff(2);
        for cch=1:length(this.currentChannels)
            set(this.currentOutline{cch},'XData',cpx+this.tp_pad,...
                'YData',cpy+this.yoffs(cch));
        end
    end
end

function [radii,angles] = radii_from_points(points,centre)
dev = (points - centre(ones(size(points,1),1),:))*[-1;1i];
radii = abs(dev); angles = pi - angle(dev);
end

function [cpoint,points] = rect_from_radii(radii,angles)
radii = radii(:);
coords = radii(:,[1,1]).*[cos(angles(:)),sin(angles(:))];
orthcentre = mean(coords);
[pcaRot,pcaCoord] = pca(coords);
pcaBoxMin = min(pcaCoord);
pcaBoxMax = max(pcaCoord);
pcaBox = [pcaBoxMin(1) pcaBoxMin(2); pcaBoxMin(1) pcaBoxMax(2); ...
    pcaBoxMax(1) pcaBoxMax(2); pcaBoxMax(1) pcaBoxMin(2)];
points = orthcentre(ones(size(pcaBox,1),1),:) + ...
    babyGUI.boxScale*pcaBox*pcaRot';
cpoint = orthcentre+babyGUI.cpointScale*pcaBoxMax(2)*pcaRot(:,2)';

% rotmat = [cos(-angles(1)) sin(-angles(1)); -sin(-angles(1)) cos(-angles(1))];
% orthangles = angles - angles(1);
% orthpoints = radii(:,ones(2,1)).*[cos(orthangles),sin(orthangles)];
% pmin = min(orthpoints);
% pmax = max(orthpoints);
% orthcentre = mean([pmin;pmax]);
% orthcentre = orthcentre(ones(4,1),:);
% points = orthcentre + ...
%     ([pmin(1) pmin(2); pmin(1) pmax(2); ...
%     pmax(1) pmax(2); pmax(1) pmin(2)]-orthcentre)*rotmat;
end
