function mesh = GetDomainsFromMeshHull(mesh,VoronoiStruct)
%
% GetDomainsFromMeshHull    assigns each point of the mesh to its Voronoi region (or to the intergrain region)
%
% iIn = GetDomainsFromMeshHull(fnameMesh,fnameGeom)
% fnameMesh is the .mat file with the mesh ([fname,'MeshData.mat'])
% fnameGeom is the .mat file with the geometry of the tessellation ([fname,'.mat'])
% iIn is a cell array of K+1 elements (where K is the number of generators)
% Each element is a N-elements array of booleans (where N is the number of
% elements in the mesh). iIn{k}(n) is equal to true if the n-th
% mesh element belongs to the k-th Voronoi region.
% The intergrain region corresponds to iIn{K+1}
% The results are also saved in [fnameMesh(1:end-4),'Analyzed.mat']

%% Compute centers of the mesh-elements
% the mesh may be tetrahedral or Cartesian unstructured 
if strcmp(mesh.type,'tetrahedral') % tetrahedral
    mesh.A = mesh.A+1 ; % vertices-elements connectivity matrix (starts from 1)
    for k=1:size(mesh.A,1)
         mesh.xC(k,:) = mean(mesh.Vert(mesh.A(k,:),:),1) ;
    end
else % Cartesian unstructured 
     mesh.xC = mesh.pos_out ;
end

%% Assign each mesh-element to a Voronoi region
if exist('VoronoiStruct','var')
[iIn,GrainIndex] = getTheDomains(mesh.xC,VoronoiStruct.vorvx) ;
else
    iIn{1} = [1:size(mesh.A,1)].' ;
    iIn{2} = [].' ;
    GrainIndex = ones(size(mesh.A,1),1) ;
    
end
%% Save Results
mesh.iIn = iIn;
mesh.GrainIndex = GrainIndex;

% if strcmp(mesh.type,'tetrahedral')
%     
% %     save([fnameMesh(1:end-4),'Analyzed.mat'],'Vert','A','xC','iIn','GrainIndex') ;
% else
% %     save([fnameMesh(1:end-4),'Analyzed.mat'],'pos_out','iIn','GrainIndex') ;
%     mesh.iIn = iIn;
%     mesh.GrainIndex = GrainIndex;
% end