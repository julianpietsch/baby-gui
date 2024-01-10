function refreshTrackPlot(this,~,~)

if ~this.ctrls.refreshtracks.Value
    % Only refresh track plot if option is toggled
    return
end

cla(this.frame.track);
axes(this.frame.track);
hold on;

trap = this.currentTrap;
cLabels = this.cellLabels{trap};
cellLabelMap = zeros(max(cLabels),1);
cellLabelMap(cLabels) = 1:numel(cLabels);
cellMothers = full(this.cTimelapse.cellMothers(trap,cLabels));

ctracks = this.cellTracks;
track_start = NaN(size(ctracks,1),1);
for ct=1:size(ctracks,1)
    starttp = find(ctracks(ct,:),1);
    if ~isempty(starttp), track_start(ct) = starttp; end
end

% Default axis limits
xlim = [min([this.times(:);0]),max([this.times(:);1])];
ylim = [0.5,1.5];

showtracks = 1:size(ctracks,1);
trackmode = fieldnames(this.trackDisplayMap);
trackmode = trackmode{this.ctrls.trackdisp.Value};
if (any(strcmp({'tracks','yloc'},trackmode)) && isempty(showtracks)) ...
        || (strcmp('mdarea',trackmode) && isempty(this.currentCell))
    trackmode = 'empty';
end

% Start by pre-calculating plot values to revise y limits
switch trackmode
    case 'tracks'
        [track_pos,orphans] = organiseTracks();
        if ~this.ctrls.orphantracks.Value
            showtracks = showtracks(~orphans);
        end
        ylim = [0.5,max([track_pos(showtracks);1])+0.5];
    case 'mdarea'
        [motherVol,daughterVols,dLabels] = getMDstats();
        ylim = [0,max([motherVol(:);daughterVols(:);1])];
    case 'yloc'
        if ~this.ctrls.orphantracks.Value
            orphans = findOrphans();
            showtracks = showtracks(~orphans);
        end
        ylocs = this.cellLocs(:,showtracks,2);
        ylim = [min(ylocs(:)),max(ylocs(:))];
        if numel(unique(ylim))==2
            ylim = ylim+0.05*[-1,1]*diff(ylim);
        else
            ylim = max([ylim(:);0])+[0,1];
        end
    case 'empty'
        % do nothing
    otherwise
        error('trackmode not yet implemented!!');
end

% Draw marks before everything else
mtps = this.isCurated(trap);
for mtp=mtps(:)'
    plot(this.times(repmat(mtp,2,1)),ylim,...
        'Color',[0,1,0],'LineWidth',1.5);
end
mtps = this.forCuration(trap);
for mtp=mtps(:)'
    plot(this.times(repmat(mtp,2,1)),ylim,...
        'Color',[0.8,0.5,0],'LineWidth',1.5);
end

switch trackmode
    case 'tracks'
        % Draw lineage lines...
        for c=showtracks
            if cellMothers(c)>0
                motherId = find(cLabels==cellMothers(c),1);
                plot(this.times(max(track_start(c)+[-1,0],1)),...
                    track_pos([motherId,c]),'Color',0.4*ones(1,3),...
                    'LineStyle',':','LineWidth',1);
            end
        end
        % ...then track lines
        for c=showtracks
            endpoints = [...
                find([this.cellTracks(c,1),diff(this.cellTracks(c,:))>0]);...
                find([diff(this.cellTracks(c,:))<0,this.cellTracks(c,end)])]';
            plot(this.times(endpoints),repmat(track_pos(c),size(endpoints)),...
                'Color',this.currentColours(c,:),'LineStyle','-',...
                'Marker','o','MarkerFaceColor',this.currentColours(c,:));
        end
    case 'mdarea'
        % First plot birth events...
        m = this.currentCell;
        btps = find(this.births(m,:));
        for btp=btps(:)'
            plot(this.times(btp(ones(2,1))),ylim,...
                'Color',0.4*ones(1,3),...
                'LineStyle',':','LineWidth',1);
        end
        % ...then mother volumes...
        plot(this.times,motherVol,'.','Color',this.currentColours(m,:));
        % ...and daughter volumes
        daughters = cellLabelMap(dLabels);
        [~,daughterOrder] = sort(track_start(daughters));
        for i=1:numel(daughters)
            d = daughterOrder(i);
            plot(this.times,daughterVols(d,:),'.','Color',...
                this.currentColours(daughters(d),:));
        end
    case 'yloc'
        ylocs = this.cellLocs(:,:,2);
        % Draw lineage lines...
        for c=showtracks
            if cellMothers(c)>0
                motherId = find(cLabels==cellMothers(c),1);
                lintps = max(track_start(c)+[-1,0],1);
                plot(this.times(lintps),[ylocs(lintps(1),motherId),...
                    ylocs(lintps(2),c)],'Color',0.4*ones(1,3),...
                    'LineStyle','-','LineWidth',1);
            end
        end
        % ...then track lines
        for c=showtracks
            plot(this.times,ylocs(:,c),...
                'Color',this.currentColours(c,:),'LineStyle','-');
        end
    case 'empty'
        % do nothing
    otherwise
        error('trackmode not yet implemented!!');
end

% Set up the marker line
curtime = this.times(this.currentTimepoint);
this.trackMarkerLine = plot(repmat(curtime,2,1),ylim,...
    'Color',[0.5,0.5,0.5],'LineStyle','-');

xlabel(this.timelabel);
axis([xlim,ylim]);

    function orphans = findOrphans()
        mothers = unique(cellMothers(cellMothers>0));
        mInds = cellLabelMap(mothers);
        orphans = cellMothers==0; orphans(mInds) = false;
    end

    function [track_order,orphans] = organiseTracks()
        % Choose order for track display by grouping lineage trees
        % First track is the first occuring mother cell without mothers
        mothers = unique(cellMothers(cellMothers>0));
        mInds = cellLabelMap(mothers);
        orphans = cellMothers==0; orphans(mInds) = false;
        seen = orphans;
        % default state is to order last
        linBirthOrder = numel(cLabels)+ones(numel(cLabels));
        % First root_node should be mother 0
        root_nodes = 0;
        for l=1:numel(cLabels)
            for n=1:numel(root_nodes)
                mLbl = root_nodes(n);
                dInds = find(cellMothers==mLbl & ~seen);
                if mLbl>0
                    mInd = cellLabelMap(mLbl);
                    seen(mInd) = true;
                    % Copy lineage tree of mother
                    linBirthOrder(dInds,:) = linBirthOrder(mInd(ones(numel(dInds),1)),:);
                    linBirthOrder(mInd,l) = 0; % order the parent node first
                end
                [~,linBirthOrder(dInds,l)] = sort(track_start(dInds));
            end
            root_nodes = mothers(~seen(mInds) ...
                & ismember(cellMothers(mInds),root_nodes));
            if isempty(root_nodes), break; end
        end
        [~,linBirthOrder(orphans,end)] = sort(track_start(orphans));
        
        % Determine track order from lineage tree
        [~,track_order] = sortrows(linBirthOrder);
        [~,track_order] = sort(track_order);
    end

    function [mVol,dVol,dLabels] = getMDstats()
        % First extract mother-daughter stats for the active cell
        mother = this.currentCell;
        mLabel = cLabels(mother);
        daughterInds = find(cellMothers==mLabel);
        dLabels = cLabels(daughterInds);
        
        motherArea = NaN(1,this.ntimepoints);
        daughterAreas = NaN(numel(dLabels),this.ntimepoints);
        
        for t=1:numel(this.tpInds)
            tp = this.tpInds(t);
            clabs = this.cTimelapse.cTimepoint(tp).trapInfo(trap).cellLabel;
            if ismember(mLabel,clabs)
                motherArea(1,t) = this.getCellArea(trap,t,mother);
            end
            dis = find(ismember(dLabels,clabs));
            for di=dis
                daughterAreas(di,t) = this.getCellArea(trap,t,daughterInds(di));
            end
        end
        
        mVol = 4/3*pi*sqrt(motherArea/pi).^3;
        dVol = 4/3*pi*sqrt(daughterAreas/pi).^3;
    end
end