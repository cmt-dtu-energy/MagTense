function [KU,vorvxU] = CleanUp(K,vorvx,tol)
'';
if ~exist('tol','var')
[vorvxU,iA,iC] = unique(vorvx, 'first' ,'rows') ;
KU = iC(K) ;
else
    Sigma = 10^round(log10(tol)) ;
    [vorvxU,iA,iC] = unique(round(vorvx./Sigma), 'first' ,'rows') ;
    vorvxU = vorvx(iA,:) ;
    KU = iC(K) ;
end