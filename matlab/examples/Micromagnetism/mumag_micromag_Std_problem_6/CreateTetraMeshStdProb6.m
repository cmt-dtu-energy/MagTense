function model = CreateTetraMeshStdProb6(L,Hmax)
%{
[xLim,yLim,zLim] = ndgrid([0,+.5,1].*L(1),[0,1].*L(2),[0,1].*L(3));

xLim = xLim(:);
yLim = yLim(:);
zLim = zLim(:);
left_region=xLim<=L(1)/2;
right_region=xLim>=L(1)/2;
% ConvHullL = convhull(xLim(left_region),yLim(left_region),zLim(left_region));
% ConvHullR = convhull(xLim(right_region),yLim(right_region),zLim(right_region));
ConvHullL = convhull(xLim.*left_region,yLim.*left_region,zLim.*left_region);
ConvHullR = convhull(xLim.*right_region,yLim.*right_region,zLim.*right_region);

nodes = [xLim';yLim';zLim'];
elements = [ConvHullL',ConvHullR'];
%}
% Load a multidomain pre-made mesh with boundary in z-dimension.
load MultidomainMesh3D
% Change boundary to x-dimension
tmp=nodes(1,:);nodes(1,:)=nodes(3,:);nodes(3,:)=nodes(2,:);nodes(2,:)=tmp;
%% Scale dimension. Boundary is in the middle as we want.
loaded_edges=[min(nodes(1,:)),max(nodes(1,:));min(nodes(2,:)),max(nodes(2,:));min(nodes(3,:)),max(nodes(3,:))];
wanted_edges=[zeros(3,1) L'];
% Shift to start from 0
nodes=nodes-loaded_edges(:,1);
% Scale
nodes=nodes.*(wanted_edges(:,2)-wanted_edges(:,1))./(loaded_edges(:,2)-loaded_edges(:,1));

model = createpde();
geometryFromMesh(model,nodes,elements,ElementIdToRegionId);
figure; pdegplot(model,'CellLabels','on')
if ~exist('Hmax','var')
 Hmax = 100e-9 ; 
end
generateMesh(model,'GeometricOrder','linear','Hmax',Hmax);