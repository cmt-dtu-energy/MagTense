

function [fig,ax] = plotSphHarmDecomp( sph, fig )

if ~exist('fig','var') || isempty(fig)
    [fig,ax] = getFigure();
end

[srt,ind] = sort( sph.orders(3,:));
sph.C = sph.C(ind);
sph.labels = sph.labels(ind);

plot(1:1:length(sph.C),log10(abs(sph.C)),'o-','linewidth',2,'markersize',10);

xticks(1:1:length(sph.C));
xticklabels(sph.labels);
xtickangle(45);

end