function refreshTrapIms(this)
% Retrieve and scale intensity of trap images

nT = this.tilesContext;
nC = length(this.currentChannels);

ctp = this.currentTimepoint;
tps = max(1,ctp-nT):min(this.ntimepoints,ctp+nT);

this.trapIms = cell(nC,1);
for c=1:nC
    channel = this.currentChannels{c};
    tempIm = this.imcache.getTrapTimepoints(...
        this.currentPos,this.currentTrap,this.currentChannels{c},tps);
    channel = regexprep(channel,'\W','_');
    if isfield(this.channelNormalisation,channel)
        imrange = this.channelNormalisation.(channel);
        tempIm = (tempIm - imrange(1))/diff(imrange);
        tempIm(tempIm(:)<0) = 0;
        tempIm(tempIm(:)>1) = 1;
    end
    this.trapIms{c} = tempIm*this.alphascale;
end

end
