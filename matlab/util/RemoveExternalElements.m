function [mesh,GridInfo] = RemoveExternalElements(mesh,GridInfo,ExtBorder,VoronoiStruct) 
IIn = inhull([GridInfo.Xel,GridInfo.Yel,GridInfo.Zel],ExtBorder.Points) ;
NotIIn = ~IIn ;

GridInfo.Xel(NotIIn) = [] ;
GridInfo.Yel(NotIIn) = [] ;
GridInfo.Zel(NotIIn) = [] ;
GridInfo.Volumes(NotIIn) = [] ;
GridInfo.TheSigns(NotIIn,:) = [] ;
OrphanFaces = sum(abs(GridInfo.TheSigns),1) == 0  ;
GridInfo.fNormX(OrphanFaces) = [] ;
GridInfo.fNormY(OrphanFaces) = [] ;
GridInfo.fNormZ(OrphanFaces) = [] ;
GridInfo.Xf(OrphanFaces) = [] ;
GridInfo.Yf(OrphanFaces) = [] ;
GridInfo.Zf(OrphanFaces) = [] ;
GridInfo.AreaFaces(OrphanFaces) = [] ;
GridInfo.TheSigns(:,OrphanFaces) = [] ;
GridInfo.DimsF(OrphanFaces,:) = [] ; 

GridInfo.TheTs(:,NotIIn) = [] ;
GridInfo.TheTs(OrphanFaces,:) = [] ;
GridInfo.TheDs(:,NotIIn) = [] ;
GridInfo.TheDs(OrphanFaces,:) = [] ;

mesh.pos_out(NotIIn,:) = [] ;
mesh.dims_out(NotIIn,:) = [] ;
mesh.children_out(NotIIn,:) = [] ;
mesh.nn(NotIIn,:) = [] ;
mesh.xC(NotIIn,:) = [] ;
if exist('VoronoiStruct','var') ;
[iIn,GrainIndex] = getTheDomains(mesh.xC,VoronoiStruct.vorvx,VoronoiStruct.ExtBorder) ;
mesh.iIn = iIn;
mesh.GrainIndex = GrainIndex;
end
% hold on 
% plot3(GridInfo.Xf,GridInfo.Yf,GridInfo.Zf,'.')
'' ;
