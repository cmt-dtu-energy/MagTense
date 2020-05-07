function GridInfo = TetrahedralMeshAnalysis(model)


VolumeOfTetrahedron01 = @(a,b,c,d) abs(dot(a-d,cross(b-d,c-d)))/6 ;



%% Dimensions


N = size(model.Mesh.Elements,2) ; % Number of elements (i.e. cells, i.e tetrahedrons)

AllElements = model.Mesh.Elements ; % (10or4) x N indexes to the corresponding nodes
AllNodes = model.Mesh.Nodes ; % 3 x M coordinates
sizeElements = size(AllElements) ;


%% Centers & volumes of elements
AllVolumes3 = zeros(N,1) ;
xC = zeros(N,1) ;
yC = zeros(N,1) ;
zC = zeros(N,1) ;

for j=1:N % runs over the elements
    TheseIj = AllElements(1:4,j) ;
    xC(j,1) = sum(AllNodes(1,TheseIj))/4 ; % centroid of the cell j
    yC(j,1) = sum(AllNodes(2,TheseIj))/4 ; % centroid of the cell j
    zC(j,1) = sum(AllNodes(3,TheseIj))/4 ; % centroid of the cell j
    AllVolumes3(j,1) = VolumeOfTetrahedron01(AllNodes(:,TheseIj(1)),AllNodes(:,TheseIj(2)),AllNodes(:,TheseIj(3)),AllNodes(:,TheseIj(4))) ;
end



%% Initialize arrays
MaxK = N*4 ; % Preliminary over-estimation of the total number of unique faces (for initialization purpose)

totf = 0 ;
lastLink = 0 ;

numKs = 0 ;
kNotJs = zeros(2*MaxK,1) ; % n1 : element index
kNotIs = zeros(2*MaxK,1) ; % n2 : element index
kNotVs = zeros(2*MaxK,1) ; % local face index (1to4) of the face of n1 that is connected with n2

signjs = zeros(1,MaxK) ;
signks = zeros(1,MaxK) ;
signvs = zeros(1,MaxK) ;
fijs = [] ;
fiis = [] ;
fivs = [] ;
js = [] ;
is = [] ;
vs = [] ;

AllFaces = zeros(MaxK,3) ;
axCf = zeros(MaxK,1) ;
ayCf = zeros(MaxK,1) ;
azCf = zeros(MaxK,1) ;

aNormX = zeros(MaxK,1) ;
aNormY = zeros(MaxK,1) ;
aNormZ = zeros(MaxK,1) ;

aFaceArea = zeros(MaxK,1) ;


LocalToGlobalFace = zeros(N,4) ;
%% Build list of unique faces and corresponding information + Sign matrix
for j=1:N % runs over the elements
    
    theseElements = AllElements(:,j) ;
    theseAA = zeros(1,N) ;
    ThisAK = false(N,4)  ;
    for k=1:4
        TheseOnes = find(AllElements == AllElements(k,j) ) ;
        [~,i2] = ind2sub(sizeElements,TheseOnes) ;
        theseAA(i2) = theseAA(i2) + 1 ; % N times 1 array of integers: number of common vertices between element j and i2
        ThisAK(i2,k) = true ;
    end
    iiOK = find(theseAA== 3 ) ; % !!!  neihbor elements : 3 vertices (2 edges) in common
    z = zeros(1,numel(iiOK)) ;
    js = [js,j+z] ;
    is = [is,iiOK] ;
    vs = [vs,1+z] ;
    
    i2 = iiOK ;
    for ii2 = 1:numel(i2) % cycle over found neighbors
        numKs = numKs + 1 ;
        kNotJs(numKs) = j ;
        kNotIs(numKs) = i2(ii2) ;
        kNotVs(numKs) =  find(~ThisAK(i2(ii2),:))  ; % Local Face index (i.e. from 1 to 4: the vertex opposite to the face)
        
        if i2(ii2)>j % New face (can't be the same element, since the number of common vertices is 3 !!!)
            totf = totf + 1 ; % global face index (one for each of the K unique faces)
            
            LocalToGlobalFace(j,kNotVs((kNotJs==j)&(kNotIs==i2(ii2))) ) = totf ;
            AllFaces(totf,:) = theseElements(ThisAK(i2(ii2),:)) ;% indices of the vertices of the new face
            
            [axCf(totf,1),ayCf(totf,1),azCf(totf,1),aNormX(totf,1),aNormY(totf,1),aNormZ(totf,1),aFaceArea(totf,1)] = ComputeFaceInfoTriangle01(AllNodes,AllFaces(totf,:),xC(j),yC(j),zC(j)) ;
            
            lastLink = lastLink + 1 ;
            signjs(lastLink) = j ;
            signks(lastLink) = totf ;
            signvs(lastLink) = 1 ;
            
            fijs = [fijs,j] ;
            fiis = [fiis,i2(ii2)] ;
            fivs = [fivs,totf] ;
            
            '' ;
        end
        
        if i2(ii2)<j % Old face (can't be the same element, since the number of common vertices is 3 !!!)
            thatf = fivs((fijs==i2(ii2))&(fiis==j)) ; % find the global index of the face
            LocalToGlobalFace(j,kNotVs((kNotJs==j)&(kNotIs==i2(ii2))) ) = thatf ;
            lastLink = lastLink + 1 ;
            signjs(lastLink) = j ;
            signks(lastLink) = thatf ;
            signvs(lastLink) = -1 ;
            
        end
        
        
    end
    
    FacesWithNoNeighbor = setdiff(1:4,kNotVs((kNotJs==j))) ;
    if ~isempty(FacesWithNoNeighbor)
        for iik = 1:numel(FacesWithNoNeighbor)
            totf = totf + 1 ;
            LocalToGlobalFace(j,FacesWithNoNeighbor(iik)) = totf ;
            lastLink = lastLink + 1 ;
            signjs(lastLink) = j ;
            signks(lastLink) = totf ;
            signvs(lastLink) = +1 ;
            AllFaces(totf,:) = AllElements(setdiff(1:4,FacesWithNoNeighbor(iik)),j) ;
            [axCf(totf,1),ayCf(totf,1),azCf(totf,1),aNormX(totf,1),aNormY(totf,1),aNormZ(totf,1),aFaceArea(totf,1)] = ComputeFaceInfoTriangle01(AllNodes,AllFaces(totf,:),xC(j),yC(j),zC(j)) ;
            
        end
        '' ;
    end
    
end
%% Remove extra entries
K = totf ; % Actual number of unique faces.
AllFaces(K+1:end,:) = [] ;
axCf(K+1:end) = [] ;
ayCf(K+1:end) = [] ;
azCf(K+1:end) = [] ;

aNormX(K+1:end) = [] ;
aNormY(K+1:end) = [] ;
aNormZ(K+1:end) = [] ;

aFaceArea(K+1:end) = [] ;

signjs(lastLink+1:end) = [] ;
signks(lastLink+1:end) = [] ;
signvs(lastLink+1:end) = [] ;
TheSigns = sparse(signjs,signks,signvs,N,K) ;
A = sparse(is,js,vs,N,N) ;
% A = A + speye(N,N) ; % include self-links
'' ;
%% Generate T matrix
% K times N adjacency matrix
% for each face find all the elements having 2 vertices in common.

ks = [] ;
is = [] ;
vs = logical([]) ;
for k =1:K
    
    
    theseAA = zeros(1,N) ;
    for i=1:3
        TheseOnes = find(AllElements == AllFaces(k,i) ) ;
        [~,i2] = ind2sub(sizeElements,TheseOnes) ;
        theseAA(i2) = theseAA(i2) + 1 ;
    end
    iiOK = find(theseAA> 1 ) ; % !!!  neihbor elements
    z = zeros(1,numel(iiOK)) ;
    zLog = true(1,numel(iiOK)) ;
    ks = [ks,k+z] ;
    is = [is,iiOK] ;
    vs = [vs,zLog] ;
end


T = sparse(ks,is,vs,K,N) ; %figure; spy(T)
'' ;
%%

GridInfo.fNormX = aNormX ;
GridInfo.fNormY = aNormY ;
GridInfo.fNormZ =  aNormZ ;
GridInfo.AreaFaces =  aFaceArea ;
GridInfo.Volumes =  AllVolumes3 ;
GridInfo.TheSigns =  TheSigns ;
GridInfo.Xel =  xC ;
GridInfo.Yel =  yC ;
GridInfo.Zel =  zC ;
GridInfo.Xf =  axCf ;
GridInfo.Yf =  ayCf ;
GridInfo.Zf =  azCf ;
GridInfo.TheTs =  T ;
GridInfo.A =  A ;

end

function [xCf,yCf,zCf,NormX,NormY,NormZ,FaceArea] = ComputeFaceInfoTriangle01(AllNodes,ThoseI,xC,yC,zC)
        xCf = sum(AllNodes(1,ThoseI))/3 ; % centroid of the face k
        yCf = sum(AllNodes(2,ThoseI))/3 ;
        zCf = sum(AllNodes(3,ThoseI))/3 ;
        
        
        DeltaX1 = AllNodes(1,ThoseI(2)) - AllNodes(1,ThoseI(1)) ;
        DeltaY1 = AllNodes(2,ThoseI(2)) - AllNodes(2,ThoseI(1)) ;
        DeltaZ1 = AllNodes(3,ThoseI(2)) - AllNodes(3,ThoseI(1)) ;
        
        DeltaX2 = AllNodes(1,ThoseI(3)) - AllNodes(1,ThoseI(2)) ;
        DeltaY2 = AllNodes(2,ThoseI(3)) - AllNodes(2,ThoseI(2)) ;
        DeltaZ2 = AllNodes(3,ThoseI(3)) - AllNodes(3,ThoseI(2)) ;
        
        NormX = DeltaY1*DeltaZ2 - DeltaZ1*DeltaY2 ; % Normal to the face
        NormY = DeltaZ1*DeltaX2 - DeltaX1*DeltaZ2 ;
        NormZ = DeltaX1*DeltaY2 - DeltaY1*DeltaX2 ;
        
        NormN = sqrt(NormX^2+NormY^2+NormZ^2) ; % Also twice the area !
        FaceArea = NormN./2 ;
                ThatProd = ((xCf-xC)*NormX + (yCf-yC)*NormY + (zCf-zC)*NormZ) ;
%         AllVolumes(j) = AllVolumes(j) + abs(ThatProd)/6 ;
        
        ThatSign = sign(ThatProd) ;
        
        
        NormX = NormX.*ThatSign./NormN ;
        NormY = NormY.*ThatSign./NormN ;
        NormZ = NormZ.*ThatSign./NormN ;
        
end