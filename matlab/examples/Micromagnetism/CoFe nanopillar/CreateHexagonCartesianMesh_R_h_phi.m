function [mesh,GridInfo,mesh_params] = CreateHexagonCartesianMesh_R_h_phi(dimensions, NNs)

Radius = dimensions(1);
Height = dimensions(2);
angle = dimensions(3);

%% Parameters
L = 1 ;                 %--- Number of refinment iterations from base_res for the cartesian grid


%% Define the platelet
PlateletTRI = CreateHexagonalPlatelet_define_R_h_phi(Radius,Height,angle) ;
X = [min(PlateletTRI.Points(:,1)),max(PlateletTRI.Points(:,1))] ;
Y = [min(PlateletTRI.Points(:,2)),max(PlateletTRI.Points(:,2))] ;
Z = [min(PlateletTRI.Points(:,3)),max(PlateletTRI.Points(:,3))] ;

 
%% Generate Cartesian Unstructured Mesh
[VoronoiStruct,x_gen,y_gen,z_gen] = SingleGrainVoronoi(X,Y,Z,PlateletTRI) ; % VoronoiStruct has only one grain, but it's used by refineRectGrid
mesh = refineRectGrid_v6( X,Y,Z,x_gen,y_gen,z_gen,NNs,L,VoronoiStruct,3) ;

%--- Analyze and plot Cartesian Unstructured Mesh
mesh = GetDomainsFromMeshHull(mesh,VoronoiStruct) ;
GridInfo = CartesianUnstructuredMeshAnalysis(mesh.pos_out,mesh.dims_out) ;
[mesh,GridInfo] = RemoveExternalElements(mesh,GridInfo,VoronoiStruct.ExtBorder,VoronoiStruct) ;

%% Save other mesh parameters
mesh.GeomTRI = PlateletTRI;

mesh_params.NumRefinments = L;
mesh_params.NNs = NNs;
mesh_params.resolution = [length(mesh.pos_out) 1 1];
mesh_params.nGrains = 1;
mesh_params.thisGridL = [Radius, Height, angle];

end

function TR = CreateHexagonalPlatelet_define_R_h_phi(R,h, Phi)
phi=Phi*pi/180; %convert angle in degrees to radians

points = [0 0 0 ;   0 0 h   ;   R*sin(phi) R*cos(phi) h ;   R*sin(phi) R*cos(phi) 0; ...
          R+R*sin(phi) R*cos(phi) h   ;   R+R*sin(phi) R*cos(phi) 0   ;   R+2*R*sin(phi) 0 h  ;   R+2*R*sin(phi) 0 0; ...
          R+R*sin(phi) -R*cos(phi) h    ;   R+R*sin(phi) -R*cos(phi) 0   ;  R*sin(phi) -R*cos(phi) h    ;   R*sin(phi) -R*cos(phi) 0;]; %the way I defined it in COMSOL.

points = points + [-(R*sin(phi) + R/2) 0 -h/2]; %shift the coordinates, so middle of geometry is in origin so it looks like Andrea's for phi=30.

rotationmatrix = [cos(pi/2) -sin(pi/2) 0 ; sin(pi/2) cos(pi/2) 0; 0 0 1];

points = points*rotationmatrix; %rotate geometry so it looks excactly like Andrea's for phi=30.

ConnectivityList=convhulln(points);

TR = triangulation(ConnectivityList,points) ;

%plot the triangulated geometry
%extraInputs = {'interpreter','latex','fontsize',13};
%figure; trisurf(TR); axis equal; %view(0,90);
%xlabel('$x$', extraInputs{:}); ylabel('$y$', extraInputs{:}); zlabel('z', extraInputs{:}); 

end

function [VoronoiStruct,x_gen,y_gen,z_gen] = SingleGrainVoronoi(X,Y,Z,PlateletTRI)

    x_gen = mean(X) ;
    y_gen = mean(Y) ;
    z_gen = mean(Z) ;
    
    VoronoiStruct.KKall = PlateletTRI.ConnectivityList ;
    VoronoiStruct.NewVertXall = PlateletTRI.Points ;
    VoronoiStruct.ThatK{1} = PlateletTRI.ConnectivityList ;
    VoronoiStruct.vorvx{1} = PlateletTRI.Points ;
    VoronoiStruct.pos = [x_gen,y_gen,z_gen] ;
    VoronoiStruct.offSetD = 0 ; 
    
    
    VoronoiStruct.ExtBorder = PlateletTRI ;
end
