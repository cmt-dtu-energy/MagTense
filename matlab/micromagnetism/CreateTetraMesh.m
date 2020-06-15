function model = CreateTetraMesh(L,Hmax)

[xLim,yLim,zLim] = ndgrid([-.5,+.5].*L(1),[-.5,+.5].*L(2),[-.5,+.5].*L(3));

xLim = xLim(:);
yLim = yLim(:);
zLim = zLim(:);
ConvHull = convhull(xLim,yLim,zLim);

nodes = [xLim';yLim';zLim'];
elements = ConvHull';

model = createpde();
geometryFromMesh(model,nodes,elements);
if ~exist('Hmax','var')
 Hmax = 100e-9 ; 
end
generateMesh(model,'GeometricOrder','linear','Hmax',Hmax);