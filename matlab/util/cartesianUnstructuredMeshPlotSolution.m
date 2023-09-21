function [hs,hq] = cartesianUnstructuredMeshPlotSolution(GeomTRI,mesh,IntPlane,Mx_arr,My_arr,Mz_arr,hF)
%
% CartesianUnstructuredMeshPlotSolution(GeomTRI,mesh,IntPlane,Mx_arr,My_arr,Mz_arr,hF)
%
% Intersects the geometry with a plane, and shows a cutaway view of the
% solution.
% 
% GeomTri is the triangulation of the geometry (can be obtained using
% convhulln. For this function to work, the geometry is assumed to be
% convex anyway.
% mesh is the Cartesian Unstructured Mesh corresponding to the solution
% IntPlane is a structure with two fields (each of which is a 1 by 3
% array): IntPlane.Normal is the normal to the plane
% IntPlane.Point is a point passing by the plane
% Mx_arr, My_arr, and Mz_arr are the arrays of the solution
%
% Requires the PDE toolbox (used for creating a tetrahedral mesh -
% and the correponding triangular surface mesh on which
% the solution is plotted)

%%
X = [min(GeomTRI.Points(:,1)),max(GeomTRI.Points(:,1))] ;
Y = [min(GeomTRI.Points(:,2)),max(GeomTRI.Points(:,2))] ;
Z = [min(GeomTRI.Points(:,3)),max(GeomTRI.Points(:,3))] ;
IntPlaneN = norm(IntPlane.Normal) ;
IntPlane.Normal = IntPlane.Normal./IntPlaneN ;

%%
ik = [1,2,3,1] ;

% keep only the points on the negative side of the plane
iIn = sum((GeomTRI.Points-repmat(IntPlane.Point,size(GeomTRI.Points,1),1)).*repmat(IntPlane.Normal,size(GeomTRI.Points,1),1),2)<0 ; 
% Cycle over each segment
jN = 0 ; clear xNew
for j=1:size(GeomTRI.ConnectivityList,1)
    for k=1:3
        % find the intersection point between the segment and the
        % intersecting plane
        xA = GeomTRI.Points(GeomTRI.ConnectivityList(j,ik(k)),:) ;
        xB = GeomTRI.Points(GeomTRI.ConnectivityList(j,ik(k+1)),:) ;
        t = -sum((xA-IntPlane.Point).*IntPlane.Normal)/sum((xB-xA).*IntPlane.Normal) ;
        xInt = xA + t.*(xB-xA) ;
        if (t>=0) & (t<=1) ; plot3(xInt(1),xInt(2),xInt(3),'.') ; 
            jN = jN + 1 ;
            xNew(jN,:)  = xInt ;


        end
    end
end

BisectGeomTRI.Points = GeomTRI.Points ; 
BisectGeomTRI.Points(find(~iIn),:) = [] ;
BisectGeomTRI.Points = [BisectGeomTRI.Points;xNew] ;

BisectGeomTRI.ConnectivityList = convhulln(BisectGeomTRI.Points) ;


BiSectVorStruct.pos = sum(BisectGeomTRI.Points,1) ;
BiSectVorStruct.vorvx{1} = BisectGeomTRI.Points ;
BiSectVorStruct.ThatK{1} = BisectGeomTRI.ConnectivityList ;

VoronoiStruct.pos = sum(GeomTRI.Points,1) ;
VoronoiStruct.vorvx{1} = GeomTRI.Points ;
VoronoiStruct.ThatK{1} = GeomTRI.ConnectivityList ;


% PlotVoronoiGeometry(BiSectVorStruct,X,Y,Z,'BisectedHexagonalPlatelet01',hF,[.5,.5,.5])
BisectGeomTRI = triangulation(BisectGeomTRI.ConnectivityList,BisectGeomTRI.Points) ;
stlwrite(BisectGeomTRI,'BisectGeomTRI.stl')
model = createpde(1);
importGeometry(model,'BisectGeomTRI.stl');
delete('BisectGeomTRI.stl') ;
generateMesh(model,'GeometricOrder','linear','Hmax',diff(X)/50) ;
hF0 = figure ;
h = pdeplot3D(model) ;
SurfMesh =  triangulation(h.Faces,h.Vertices) ;
delete(hF0)

%% Mesh for quiver plots
generateMesh(model,'GeometricOrder','linear','Hmin',diff(X)/10) ;
hF0 = figure ;
h = pdeplot3D(model) ;
QuiverMesh =  triangulation(h.Faces,h.Vertices) ;
inQuiv = sum((QuiverMesh.Points-repmat(IntPlane.Point,size(QuiverMesh.Points,1),1)).*repmat(IntPlane.Normal,size(QuiverMesh.Points,1),1),2)> (-.01*diff(X)) ; 
delete(hF0)
%%

Fmx = scatteredInterpolant(mesh.pos_out(:,1),mesh.pos_out(:,2),mesh.pos_out(:,3),Mx_arr) ;
Fmy = scatteredInterpolant(mesh.pos_out(:,1),mesh.pos_out(:,2),mesh.pos_out(:,3),My_arr) ;
Fmz = scatteredInterpolant(mesh.pos_out(:,1),mesh.pos_out(:,2),mesh.pos_out(:,3),Mz_arr) ;

Mx = Fmx(SurfMesh.Points(:,1),SurfMesh.Points(:,2),SurfMesh.Points(:,3)) ;
My = Fmy(SurfMesh.Points(:,1),SurfMesh.Points(:,2),SurfMesh.Points(:,3)) ;
Mz = Fmz(SurfMesh.Points(:,1),SurfMesh.Points(:,2),SurfMesh.Points(:,3)) ;


MxQ = Fmx(QuiverMesh.Points(:,1),QuiverMesh.Points(:,2),QuiverMesh.Points(:,3)) ;
MyQ = Fmy(QuiverMesh.Points(:,1),QuiverMesh.Points(:,2),QuiverMesh.Points(:,3)) ;
MzQ = Fmz(QuiverMesh.Points(:,1),QuiverMesh.Points(:,2),QuiverMesh.Points(:,3)) ;
MxQ(~inQuiv) = 0 ; MyQ(~inQuiv) = 0 ; MzQ(~inQuiv) = 0 ; 
QuiverMesh2.Points = QuiverMesh.Points + (.01*diff(X)).*repmat(IntPlane.Normal,size(QuiverMesh.Points,1),1) ;
ThisColTriangle = ColorFromHorPsiTheta01(Mx,My,Mz) ;
%%


if ~exist('hF','var')
    hF = figure('position',[0 0 600 600],'Color',[1 1 1]);
    ppsz = 1.*[20,19] ;
    ppps = 1.*[0,0,20,19] ;
    set(gcf,'PaperUnits','centimeters','PaperSize',ppsz,'PaperPosition',ppps) ;
else
        ppsz=get(gcf,'PaperSize');
    ppps=get(gcf,'PaperPosition');
    if isequal(get(hF,'type'),'figure')
    set(hF,'PaperUnits','centimeters')

    figure(hF)
%     clf(hF)
    else
       axes(hF) ; 
    end
end

hold on

axis('equal')

%%


% PlotVoronoiGeometry(VoronoiStruct,X,Y,Z,'HexagonalPlatelet01',hF,[.5,.5,.5])
PlotOnlyRealEdges(VoronoiStruct.pos,VoronoiStruct.vorvx,VoronoiStruct.ThatK,1)
PlotOnlyRealEdges(BiSectVorStruct.pos,BiSectVorStruct.vorvx,BiSectVorStruct.ThatK,2)


hs = trisurf(SurfMesh.ConnectivityList,SurfMesh.Points(:,1),SurfMesh.Points(:,2),SurfMesh.Points(:,3),'FaceVertexCData',ThisColTriangle,'linestyle','none','facelighting','flat','SpecularStrength',0,'facealpha',1) ;
hq = quiver3(QuiverMesh2.Points(:,1),QuiverMesh2.Points(:,2),QuiverMesh2.Points(:,3),MxQ,MyQ,MzQ,'color','k','LineWidth',1.5) ;

view(90,30) ;
light('position',[0,0,2].*[X(end),Y(end),Z(end)])
light('position',[0,-2,0].*[X(end),Y(end),Z(end)],'color',.3.*[1,1,1])

set(gca,'xlim',[X(1),X(end)],...
    'ylim',[Y(1),Y(end)],...
    'zlim',[Z(1),Z(end)],'visible','off')

% xlabel('x') ; ylabel('y') ; zlabel('z') ;

end

function PlotOnlyRealEdges(pos,vorvx,ThatK,LW)

TotFaces = 0 ;
if ~exist('col','var')
col = [hsv(size(pos,1));[1,1,1].*0.4];
end
for i = 1:size(pos,1)
    K = ThatK{i} ;
    %             K = convhulln(vorvx{i});
    %             [K,vorvx{i}] = CleanUp(K,vorvx{i}) ;
    TheSharedEdges{i} = zeros(size(K,1),size(K,1)) ;
%     trisurf(K,vorvx{i}(:,1),vorvx{i}(:,2),vorvx{i}(:,3),'FaceColor',col(i,:),'FaceAlpha',.5,'EdgeAlpha',0,'facelighting','flat','SpecularStrength',0)
    hold on;
    %             stlwrite(TR,['testSTLwrite',num2str(i),'.stl']) ;
    for k=1:size(K,1)  % cycle over triangles
        nn{i}{k} = cross((vorvx{i}(K(k,2),:)-vorvx{i}(K(k,1),:)),(vorvx{i}(K(k,3),:)-vorvx{i}(K(k,1),:))) ;
        nn{i}{k} =  nn{i}{k}./norm( nn{i}{k}) ;
        %                line(sigma.*(vorvx{i}(K(k,:),1)-pos(i,1))+pos(i,1),sigma.*(vorvx{i}(K(k,:),2)-pos(i,2)) + pos(i,2),sigma.*(vorvx{i}(K(k,:),3)-pos(i,3)) + pos(i,3))
    end
    ThatDot = zeros(size(K,1),size(K,1)) ;
    for k1=1:size(K,1)
        for k2=(k1+1):size(K,1)
            ThatDot(k1,k2) = dot( nn{i}{k1}, nn{i}{k2} ) ;
        end
    end
    ThatDot = ThatDot+ThatDot.' ; % triangles belonging to the same planes
    ijk = [1,2,3,1] ;
    CommonEdges = zeros(size(K,1),3) ;
    for k1=1:size(K,1)
        for k2=1:size(K,1)
            [C,ia,ib] = intersect( K(k1,:),K(k2,:)) ;
            if numel(C)==2 & abs(abs(ThatDot(k1,k2))-1)<1e-9
                ia = sort(ia) ;
                ib = sort(ib) ;
                TheSide1 = isequal(ia,[1;2]).*1 +isequal(ia,[2;3]).*2 + isequal(ia,[1;3]).*3 ;
                TheSide2 = isequal(ib,[1;2]).*1 +isequal(ib,[2;3]).*2 + isequal(ib,[1;3]).*3 ;
                CommonEdges(k1,TheSide1) = 1 ;
                CommonEdges(k2,TheSide2) = 1 ;
                TheSharedEdges{i}(k1,k2) = 1 ;   TheSharedEdges{i}(k2,k1) = 1 ;
                %                        TheSharedEdges{i}(k2,TheSide2) = k1 ;
                '' ;
            end
        end
    end
    
    '' ;
    %              FindPolygons01(TheSharedEdges{i})
    CCC = zeros(max(K(:)),max(K(:))) ;
    for k=1:size(K,1)
        for j= 1:3
            if ~CommonEdges(k,j)
                CCC(K(k,ijk(j)),K(k,ijk(j+1))) = 1 ;
                CCC(K(k,ijk(j+1)),K(k,ijk(j))) = 1 ;
                %                   line(vorvx{i}(K(k,ijk(j:j+1)),1),vorvx{i}(K(k,ijk(j:j+1)),2),vorvx{i}(K(k,ijk(j:j+1)),3),'color','k','linewidth',2)
            end
        end
    end
    bins = conncomp(graph(TheSharedEdges{i})) ;
    for kb = 1:max(bins)
        BB = bins==kb ; kB = K(BB,:) ;
        [ukB,ikB,ckB] = unique(kB(:)) ;
        
        GG = graph(CCC(ukB,ukB)) ;
        nOrd = dfsearch(GG,1) ;
        iiord = kB(ikB(nOrd)) ;
        iiord = [iiord(:);iiord(1)];
        AllFaces{i,kb} = iiord ;
        TotFaces = TotFaces + 1 ;
        AllSingleFaces{TotFaces} = iiord ;
        AllSingleFacesZones(TotFaces) = i ;
        line(vorvx{i}(iiord,1),vorvx{i}(iiord,2),vorvx{i}(iiord,3),'color','k','linewidth',LW)
    end
    '' ;
    
    
end
end