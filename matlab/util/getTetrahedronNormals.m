
%By Kaspar K. Nielsen, June 2020 (kasparkn@gmail.com)
%finds the outer surface normals of a the tetrahedra in DT tetrahedron 

function NM = getTetrahedronNormals( DT )


%get the free boundary facets
F = DT.freeBoundary;
TR = triangulation(F,DT.Points);
FN = TR.faceNormal;
%gives the normal vectors to each outer face of the mesh
normals = zeros(length(F(:,1)),3);
areas = zeros(length(F(:,1)),1);
centers = zeros(length(F(:,1)),3);
%gives the indices of the outer faces to their respective tetrahedra
normals_list = zeros(length(F(:,1)),1);

for i=1:length(F(:,1))
    for j=1:length(DT.ConnectivityList(:,1))
       f1 = find( F(i,1) == DT.ConnectivityList(j,:), 1 ); 
       f2 = find( F(i,2) == DT.ConnectivityList(j,:), 1 ); 
       f3 = find( F(i,3) == DT.ConnectivityList(j,:), 1 ); 
       if ~isempty(f1) && ~isempty(f2) && ~isempty(f3)
           %found the face!
          normals_list(i) = j;
          normals(i,:) = FN(i,:);
          %vertices of triangle          
          P = DT.Points(F(i,:),:);
          L1 = sqrt(sum( ( P(1,:)-P(2,:) ).^2 ) );
          L2 = sqrt(sum( ( P(1,:)-P(3,:) ).^2 ) );
          L3 = sqrt(sum( ( P(2,:)-P(3,:) ).^2 ) );
          h = (L1+L2+L3)/2;
          %area (Heron's formula) should be checked and verified
          areas(i) = sqrt( h*(h-L1)*(h-L2)*(h-L3) );
          centers(i,:) = mean(P);
          break;
       end
    end
end

NM = struct('normals',normals,'normals_list',normals_list, 'areas', areas, 'centers',centers);
end