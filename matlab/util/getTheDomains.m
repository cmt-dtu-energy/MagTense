function [iIn,TheKs] = getTheDomains(xC,vorvx,ExtBorder)

inAgrain = false(size(xC,1),1) ;
COL = hsv(numel(vorvx)) ;
TheKs = zeros(size(xC,1),1) ;
for k=1:numel(vorvx)
    inThis = inhull(xC,vorvx{k}) ;
    iIn{k} = find(inThis) ;
    inAgrain = inAgrain | inThis ;
    TheKs(iIn{k}) = k ;
end

if exist('ExtBorder','var')
iIn{k+1} = find(~inAgrain | ~inhull(xC,ExtBorder.Points)) ; % elements that are not in any of the regions
else
    iIn{k+1} = find(~inAgrain) ;
end
TheKs(iIn{k+1}) = k+1 ;