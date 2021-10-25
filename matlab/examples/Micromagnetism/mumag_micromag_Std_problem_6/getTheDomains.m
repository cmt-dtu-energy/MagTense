function [iIn,TheKs] = getTheDomains(xC,vorvx)

inAgrain = false(size(xC,1),1) ;
COL = hsv(numel(vorvx)) ;
TheKs = zeros(size(xC,1),1) ;
for k=1:numel(vorvx)
    inThis = inhull(xC,vorvx{k}) ;
    iIn{k} = find(inThis) ;
    inAgrain = inAgrain | inThis ;
    TheKs(iIn{k}) = k ;
end

iIn{k+1} = find(~inAgrain) ; % elements that are not in any of the regions
TheKs(iIn{k+1}) = k+1 ;