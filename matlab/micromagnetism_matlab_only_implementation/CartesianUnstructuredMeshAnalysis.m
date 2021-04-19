function GridInfo = CartesianUnstructuredMeshAnalysis(pos, dims)
%
% start with all the elements centers ("pos") and cell sizes ("dims")
% round centers and sizes to minimum
% the elements are already the right number
% generate 6 faces for each of them an connect them (TheSigns matrix)
% there are too many faces
% cycle over each pair of faces.
% For each pair (A,B) these are the possible relations
%       A equals B
%       otherwise
%          A contains B
%          B contains A
%          otherwise
%               A shares an edge with B
%               A shares no edge with B
% if A contains {B_i}, then A needs to be removed
% the boundaries of the corresponding element now include all the {B_i}
% TheSigns is equal to -1 since the normals point inwards

%%
Xel = pos(:,1) ;
Yel = pos(:,2) ;
Zel = pos(:,3) ;
Nel = size(pos,1) ;
Volumes = dims(:,1).*dims(:,2).*dims(:,3) ; % not rescaled !
%% re-scaling
DimsScales = min(dims,[],1)./2 ;
Xel = round(Xel./DimsScales(1)) ;
Yel = round(Yel./DimsScales(2)) ;
Zel = round(Zel./DimsScales(3)) ;
XXel = [Xel,Yel,Zel] ;
dims = dims./repmat(DimsScales,Nel,1) ;

'' ;
%% construct all faces (6 for each element)
K = 6*Nel ;
fNormX = zeros(K,1) ;
fNormY = zeros(K,1) ;
fNormZ = zeros(K,1) ;
AreaFaces = zeros(K,1) ;
DimsF = zeros(K,3) ;
TheSigns = repmat(eye(Nel),1,6) ;
ThePM = [-1,+1] ;
TheEE = eye(3) ;
TheNot = [2,3;3,1;1,2] ;
j=0 ;
for idim = 1:3
    for ipm= 1:2
        j = j + 1 ;
        TheJ(idim,ipm) = j ;
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
ItContains = false(K,K) ; % ItContains(k1,k2) is one if k1 contains k2, zero otherwise
for idim = 1:3
    % along a given dimension
    % compare all the A faces with ipm = -1 with all the B faces with ipm=+1
    Aindex = (1+(TheJ(idim,1)-1)*Nel):(TheJ(idim,1)*Nel) ; % ipm=-1
    Bindex = (1+(TheJ(idim,2)-1)*Nel):(TheJ(idim,2)*Nel) ; % ipm=+1
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
        UAcontainsUB = ((UminA<= UminB(kb)) & (UmaxA>= UmaxB(kb))) ;
        VAcontainsVB = ((VminA<= VminB(kb)) & (VmaxA>= VmaxB(kb))) ;
        
        UBcontainsUA = ((UminB(kb)<= UminA) & (UmaxB(kb)>= UmaxA)) ;
        VBcontainsVA = ((VminB(kb)<= VminA) & (VmaxB(kb)>= VmaxA)) ;
        
        AcontainsB = SamePosAlongDim & UAcontainsUB & VAcontainsVB ;
        BcontainsA = SamePosAlongDim & UBcontainsUA & VBcontainsVA ;
        
        ItContains(sub2ind([K,K],Aindex,repmat(Bindex(kb),1,Nel))) = AcontainsB ;
        ItContains(sub2ind([K,K],repmat(Bindex(kb),1,Nel),Aindex)) = BcontainsA ;
    end
end
ItContains = sparse(ItContains) ;
%% Remove faces containing other faces

% figure ; spy(ItContains,'.') ; hold on ;
% spy( ItContains&(ItContains.'),'or') ;
%  spy( ItContains&~(ItContains.'),'og')
ItsContained = (ItContains.') ;
thisBoolean = ItContains&ItsContained ;
MutualContainsInd = find(thisBoolean) ;
[k1Mut,k2Mut] = ind2sub([K,K],MutualContainsInd) ; % isequal(sort(k1Mut),sort(k2Mut)) : TRUE
iYes = k1Mut> k2Mut ;
k1Mut = k1Mut(iYes) ; k2Mut = k2Mut(iYes) ; % rermove only one of the two

thisBoolean = not(ItsContained)   ;
thisBoolean = ItContains&thisBoolean ;
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
%% Construct T matrix
% TheTs is a K times Nel sparse matrix
% The entry TheTs(k,n) is 1 if the n-th element shares at least one edge
% (i.e. two points) with the k-th face
% each face has 4 (A,B,C,D) points with coordinates(xp,yp,zp)
% a point is on the boundary of an element if
%    xp in [xel-dim(1)/2,xel+dim(1)/2]
%    yp in [yel-dim(2)/2,yel+dim(2)/2]
%    zp in [zel-dim(3)/2,zel+dim(3)/2]
TheTs = zeros(K,Nel) ;

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
end
% figure ; spy(TheTs) ;



%% re-scaling (inverse)
Xel = (Xel.*DimsScales(1)) ;
Yel = (Yel.*DimsScales(2)) ;
Zel = (Zel.*DimsScales(3)) ;
Xf = (Xf.*DimsScales(1)) ;
Yf = (Yf.*DimsScales(2)) ;
Zf = (Zf.*DimsScales(3)) ;
% areas and volumes are already in the original scale
% FaceNormals are already unit-normals
%%

GridInfo.fNormX = fNormX ;
GridInfo.fNormY = fNormY ;
GridInfo.fNormZ = fNormZ ;
GridInfo.AreaFaces = AreaFaces ;
GridInfo.Volumes = Volumes ;
GridInfo.TheSigns = TheSigns ;
GridInfo.TheTs = TheTs ;
GridInfo.Xel = Xel ;
GridInfo.Yel = Yel ;
GridInfo.Zel = Zel ;
GridInfo.Xf = Xf ;
GridInfo.Yf = Yf ;
GridInfo.Zf = Zf ;