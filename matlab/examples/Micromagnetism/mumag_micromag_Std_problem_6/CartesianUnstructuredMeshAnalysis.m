function GridInfo = CartesianUnstructuredMeshAnalysis(pos, dims)
%
% CartesianUnstructuredMeshAnalysis analyze a Cartesian Unstructured Mesh
%
% GridInfo = CartesianUnstructuredMeshAnalysis(pos, dims)

% pos is an array with positions of all the mesh-elements
% dims is an array with dimensions of all the mesh-elements
% GridInfo is a structure with different information about the mesh 
% (e.g. normal to the faces, faces-elements connectivity matrix, etc.).
% It is used to compute the differential operators (and to plot the mesh)
%
% start with all the elements centers ("pos") and cell sizes ("dims")
% divide centers and sizes by minimum element-size and round.
% This way we can work with integers.
% The number of elements is already correct, while some of the faces are
% shared or "contained" by other faces. 
% First, we generate 6 faces for each mesh element an connect each of the faces 
% to the corresponding mesh elements (TheSigns matrix)
% There are now too many faces
% So we cycle over each pair of faces.
% For each pair of faces (A,B) these are these possible relations:
%       A equals B
%       otherwise
%          A contains B
%          B contains A
%          otherwise
%               A shares an edge with B
%               A shares no edge with B
% If the face A contains the faces {B_i}, (meaning it is subdivided into smaller faces)
% then A needs to be removed.
% The boundaries of the corresponding element now include all the {B_i}
% TheSigns is equal to -1 since the normals point inwards

%%
Xel = pos(:,1) ;
Yel = pos(:,2) ;
Zel = pos(:,3) ;
Nel = size(pos,1) ;
Volumes = dims(:,1).*dims(:,2).*dims(:,3) ; % not rescaled !
%% Re-scaling and rounding (i.e. converting to integers)
DimsScales = min(dims,[],1)./2 ;
Xel = round(Xel./DimsScales(1)) ;
Yel = round(Yel./DimsScales(2)) ;
Zel = round(Zel./DimsScales(3)) ;
XXel = [Xel,Yel,Zel] ;
dims = dims./repmat(DimsScales,Nel,1) ;

'' ;
%% Constructing all the (unique or non-unique) faces (i.e. 6 for each element)
K = 6*Nel ;
fNormX = zeros(K,1) ;
fNormY = zeros(K,1) ;
fNormZ = zeros(K,1) ;
AreaFaces = zeros(K,1) ;
DimsF = zeros(K,3) ;
TheSigns = repmat(speye(Nel),1,6) ;
ThePM = [-1,+1] ;
TheEE = eye(3) ;
TheNot = [2,3;3,1;1,2] ; % 1st row (x) gives 2 and 3 (i.e. y and z). etc
j=0 ;
for idim = 1:3 % face along x, y, or z direction
    for ipm= 1:2 % face one side of the cube or the opposite side
        j = j + 1 ;
        TheJ(idim,ipm) = j ; % face index (from 1 to 6)
        % the "global" face index is 1+(j-1)*Nel
        fNormX((1+(j-1)*Nel):(j*Nel),1) = TheEE(idim,1).*ThePM(ipm) ;
        fNormY((1+(j-1)*Nel):(j*Nel),1) = TheEE(idim,2).*ThePM(ipm) ;
        fNormZ((1+(j-1)*Nel):(j*Nel),1) = TheEE(idim,3).*ThePM(ipm) ;
        Xf((1+(j-1)*Nel):(j*Nel),1) = Xel+TheEE(idim,1).*ThePM(ipm).*dims(:,1)./2 ;
        Yf((1+(j-1)*Nel):(j*Nel),1) = Yel+TheEE(idim,2).*ThePM(ipm).*dims(:,2)./2 ;
        Zf((1+(j-1)*Nel):(j*Nel),1) = Zel+TheEE(idim,3).*ThePM(ipm).*dims(:,3)./2 ;
        DimsF((1+(j-1)*Nel):(j*Nel),:) = [dims(:,1).*(~TheEE(idim,1)),dims(:,2).*(~TheEE(idim,2)),dims(:,3).*(~TheEE(idim,3))] ; % rescaled
        AreaFaces((1+(j-1)*Nel):(j*Nel),1) = (DimsScales(TheNot(idim,1)).*dims(:,TheNot(idim,1))).*(DimsScales(TheNot(idim,2)).*dims(:,TheNot(idim,2))) ; % not rescaled !
    end
end
XXf = [Xf,Yf,Zf] ;
%% Check which faces are contained by other faces
ItContains = sparse(K,K) ; % ItContains(k1,k2) is one if k1 contains k2, zero otherwise
ItContains = logical(ItContains) ;
for idim = 1:3
    % along a given dimension (x, y, or z)
    % compare all the A faces with ipm = -1 with all the B faces with ipm = +1
    % (meaning that, e.g. one is a "left facing" face and the other one is a
    % "right facing" face.
    Aindex = (1+(TheJ(idim,1)-1)*Nel):(TheJ(idim,1)*Nel) ; % ipm=-1
    Bindex = (1+(TheJ(idim,2)-1)*Nel):(TheJ(idim,2)*Nel) ; % ipm=+1
    % U and V correspond to the other two dimensions, 
    % i.e. y and z for idim=1 (x), and so on.
    % here we get the (integer) coordinates of the vertexes,
    % i.e. coordinate of the face center +/- half of the face side-length
    % along the U or V dimensions.
   
    UminA = XXf(Aindex,TheNot(idim,1)) - .5.*DimsF(Aindex,TheNot(idim,1)) ;
    UmaxA = XXf(Aindex,TheNot(idim,1)) + .5.*DimsF(Aindex,TheNot(idim,1)) ;
    VminA = XXf(Aindex,TheNot(idim,2)) - .5.*DimsF(Aindex,TheNot(idim,2)) ;
    VmaxA = XXf(Aindex,TheNot(idim,2)) + .5.*DimsF(Aindex,TheNot(idim,2)) ;
    
    UminB = XXf(Bindex,TheNot(idim,1)) - .5.*DimsF(Bindex,TheNot(idim,1)) ;
    UmaxB = XXf(Bindex,TheNot(idim,1)) + .5.*DimsF(Bindex,TheNot(idim,1)) ;
    VminB = XXf(Bindex,TheNot(idim,2)) - .5.*DimsF(Bindex,TheNot(idim,2)) ;
    VmaxB = XXf(Bindex,TheNot(idim,2)) + .5.*DimsF(Bindex,TheNot(idim,2)) ;
    
    for kb=1:Nel
        SamePosAlongDim = (XXf(Aindex,idim) == XXf(Bindex(kb),idim)) ; 
        % SamePosAlongDim is true if A and B are aligned along the current (idim) dimension
        UAcontainsUB = ((UminA<= UminB(kb)) & (UmaxA>= UmaxB(kb))) ;
        VAcontainsVB = ((VminA<= VminB(kb)) & (VmaxA>= VmaxB(kb))) ;
        
        UBcontainsUA = ((UminB(kb)<= UminA) & (UmaxB(kb)>= UmaxA)) ;
        VBcontainsVA = ((VminB(kb)<= VminA) & (VmaxB(kb)>= VmaxA)) ;
        
        AcontainsB = SamePosAlongDim & UAcontainsUB & VAcontainsVB ; 
        BcontainsA = SamePosAlongDim & UBcontainsUA & VBcontainsVA ;
        
        ItContains(sub2ind([K,K],Aindex,repmat(Bindex(kb),1,Nel))) = AcontainsB ;
        ItContains(sub2ind([K,K],repmat(Bindex(kb),1,Nel),Aindex)) = BcontainsA ;
        % AcontainsB and BcontainsA can both be true at the same time
        % if the two faces are coincident.
    end
end
ItContains = sparse(ItContains) ;
%% Remove faces containing other faces

ItsContained = (ItContains.') ;
thisBoolean = ItContains&ItsContained ;
MutualContainsInd = find(thisBoolean) ; % MutualContainsInd is true if the two faces are coincident
[k1Mut,k2Mut] = ind2sub([K,K],MutualContainsInd) ; % isequal(sort(k1Mut),sort(k2Mut)) : TRUE
iYes = k1Mut> k2Mut ;
k1Mut = k1Mut(iYes) ; k2Mut = k2Mut(iYes) ; % rermove only one of the two from the list of all faces
if 0 % Old Way
    thisBoolean = not(ItsContained)   ;
    thisBoolean = ItContains&thisBoolean ;
else % New Way
    thisBoolean = ItContains ;
    jkC = find(ItContains) ;
    thisBoolean(jkC) =  not(ItsContained(jkC)) ;
end
NonMutualContainsInd = find(thisBoolean) ; % number of contained faces
[k1NonMut,k2NonMut] = ind2sub([K,K],NonMutualContainsInd) ; % isempty(intersect(k1Mut,k1NonMut)) : TRUE


% each of the faces being removed has:
%   kSurv: one or more surviving contained (or equal) faces
%   one element (nRmv) having the face as one of its 6 original boundaries
%   TheSigns(nRmv,kSurv) = -1 ;
kRmv = [k1Mut;k1NonMut] ;
kSurv = [k2Mut;k2NonMut] ;
nRmv = rem(kRmv-1,Nel)+1 ;

TheSigns(sub2ind([Nel,K],nRmv,kSurv)) = -1 ;


XXf(kRmv,:) = [] ;
DimsF(kRmv,:) = [] ;

fNormX(kRmv) = [] ;
fNormY(kRmv) = [] ;
fNormZ(kRmv) = [] ;
Xf(kRmv) = [] ;
Yf(kRmv) = [] ;
Zf(kRmv) = [] ;
AreaFaces(kRmv) = [] ;
TheSigns(:,unique(kRmv)) = [] ;
K = numel(Xf) ;
% figure ; spy(TheSigns) ;
'' ;
%% Construct T and D matrix
% TheTs is a K times Nel sparse matrix
% The entry TheTs(k,n) is 1 if the n-th element shares at least one edge
% (i.e. two points) with the k-th face
% The entry TheDs(k,n) is 1 if the n-th element shares at least one vertex
% with the k-th face
% each face has 4 (A,B,C,D) points with coordinates(xp,yp,zp)
% a point is on the boundary of an element if
%    xp in [xel-dim(1)/2,xel+dim(1)/2]
%    yp in [yel-dim(2)/2,yel+dim(2)/2]
%    zp in [zel-dim(3)/2,zel+dim(3)/2]
% TheTs = zeros(K,Nel) ;
% TheDs = zeros(K,Nel) ;
TheTs = spalloc(K,Nel,50*Nel) ;
TheDs = spalloc(K,Nel,50*Nel) ;

xVertA = XXf+DimsF./2 ;
xVertC = XXf-DimsF./2 ;
DimsF2 = DimsF ;
iZero = sum(repmat([1,2,3],size(DimsF,1),1).*(DimsF==0),2) ;
theseii = sub2ind([K,3],(1:K).',rem(iZero+1,3)+1) ;
DimsF2(theseii) = -DimsF2(theseii) ;
xVertB = XXf+DimsF2./2 ;
xVertD = XXf-DimsF2./2 ;

xxMinEl = XXel - dims./2 ;
xxMaxEl = XXel + dims./2 ;

for n=1:Nel
    TheseMinXXel = repmat(xxMinEl(n,:),K,1) ;
    TheseMaxXXel = repmat(xxMaxEl(n,:),K,1) ;
    theseA = all((TheseMinXXel <= xVertA ) & (TheseMaxXXel >= xVertA ),2) ;
    theseB = all((TheseMinXXel <= xVertB ) & (TheseMaxXXel >= xVertB ),2) ;
    theseC = all((TheseMinXXel <= xVertC ) & (TheseMaxXXel >= xVertC ),2) ;
    theseD = all((TheseMinXXel <= xVertD ) & (TheseMaxXXel >= xVertD ),2) ;
    theseNum = sum([theseA,theseB,theseC,theseD],2) ;
    theseKs = find(theseNum>=2) ;
    TheTs(sub2ind([K,Nel],theseKs,repmat(n,size(theseKs,1),1))) = 1 ;
    theseLs = find(theseNum>=1) ;
    TheDs(sub2ind([K,Nel],theseLs,repmat(n,size(theseLs,1),1))) = 1 ;
end
% figure ; spy(TheTs) ;

%% re-scaling (inverse)
% elements and faces coordinates are multiplied by the minimum size
% to return to the original scale.
Xel = (Xel.*DimsScales(1)) ;
Yel = (Yel.*DimsScales(2)) ;
Zel = (Zel.*DimsScales(3)) ;
Xf = (Xf.*DimsScales(1)) ;
Yf = (Yf.*DimsScales(2)) ;
Zf = (Zf.*DimsScales(3)) ;
% areas and volumes are already in the original scale
% FaceNormals are already unit-normals
%% Collect every info on the output structure

GridInfo.fNormX = fNormX ;
GridInfo.fNormY = fNormY ;
GridInfo.fNormZ = fNormZ ;
GridInfo.AreaFaces = AreaFaces ;
GridInfo.Volumes = Volumes ;
GridInfo.TheSigns = TheSigns ;
GridInfo.TheTs = TheTs ;
GridInfo.TheDs = TheDs ;
GridInfo.Xel = Xel ;
GridInfo.Yel = Yel ;
GridInfo.Zel = Zel ;
GridInfo.Xf = Xf ;
GridInfo.Yf = Yf ;
GridInfo.Zf = Zf ;