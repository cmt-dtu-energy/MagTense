function VoronoiStruct = GenerateVoronoiCells(X,Y,Z,x,y,z,offSetD,fnameSave)
%
% GenerateVoronoiCells   Calculates the shape of the regions of the Voronoi
% tessellation and generates two output files. 
% 
% GenerateVoronoiCells(X,Y,Z,x,y,z,offSetD,fnameSave)
% x,y, and z are the positions of the N generators
% X,Y,Z are either 3-D arrays corresponding to the grid,
% or 2-elements arrays corresponding to the boundaries of the box
% offSetD is the half of the thickness of the intergrain region
% fnameSave is the name of the two output files: 
% a .mat file and an .stl file. Both contain the geometry of the Voronoi regions
% The .stl file is then passed to COMSOL
% The .mat file is used to assign each point of the mesh to its Vornoi region
% (or to the intergrain region)

%%
if ~exist('fnameSave','var')
    fnameSave = 'ThisVoronoiMesh05' ;
end
pos(:,1) = x ; pos(:,2) = y ; pos(:,3) = z ; % pack x,y,z



LLs = [X(end)-X(1),Y(end)-Z(1),Z(end)-Z(1)] ; % Size of Box
LLscale = 10^round(log10(max(LLs))) ; % Global Order of Magnitude


%% all the 8 vertices of box

bnd0(:,1) = [0;0;0;0;1;1;1;1].*LLs(1) + X(1) ; 
bnd0(:,2) = [0;0;1;1;0;0;1;1].*LLs(2) + Y(1) ;
bnd0(:,3) = [0;1;0;1;0;1;0;1].*LLs(3) + Z(1) ;
bnd_pnts = bnd0 ;

%% Compute convex hulls
[vornb,vorvx] = polybnd_voronoi_OffSet(pos./LLscale,bnd_pnts./LLscale,offSetD./LLscale);
for ih = 1:numel(vornb)
    vorvx{ih} =  vorvx{ih}.*LLscale ;
end


%% Remove duplicate points

for i = 1:size(pos,1)
    
    K = convhulln(vorvx{i});
    [~,vorvx{i}] = CleanUp(K,vorvx{i}) ;
    
end

%% Collect all the Voronoi cells

NewVertXall= cell2mat(vorvx.') ;
n = 0 ;

for i = 1:size(pos,1)
    ThatK{i} = convhulln(vorvx{i});
    KK{i} = n+ThatK{i};
    
    n = n + max(ThatK{i}(:)) ;
    
end

KKall = cell2mat(KK.') ;
%% Triangulate the tiles and save the output files
TR = triangulation(KKall,NewVertXall);


stlwrite(TR,[fnameSave,'.stl']) ;
save([fnameSave,'.mat'],'KKall','NewVertXall','ThatK','vorvx','pos','offSetD') ;
VoronoiStruct.KKall = KKall;
VoronoiStruct.NewVertXall = NewVertXall;
VoronoiStruct.ThatK = ThatK;
VoronoiStruct.vorvx = vorvx;
VoronoiStruct.pos = pos;
VoronoiStruct.offSetD = offSetD;

