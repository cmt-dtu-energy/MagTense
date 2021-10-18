function D2 = ComputeDifferentialOperatorsFromMesh_Neumann_FiniteDifference(GridInfo,Aexch)
%
% [DX,DY,DZ,W] = ComputeDifferentialOperatorsFromMesh04_GGDirLap(NX,NY,NZ,Areas,Volumes,Signs,Xel,Yel,Zel,Xf,Yf,Zf,T,interpn,weight,method)
%
% References: "Gradient Calculation Methods on Arbitrary Polyhedral Unstructured Meshes for Cell-Centered CFD Solvers"
%             "DifferentialOperatorMeshes" in the documentation folder in MagTense git repository
%
% Creates differential operators for an unstructured mesh
% The mesh is composed by
% N elements (or tiles, or cells)
% K faces 
% 
% There are two methods in this function: Green Gauss and Direct Laplacian.
% For Green Gauss (Direct Laplacian), the goal is to evaluate the (second) derivative of a function 
% Phi in the center of each element n from the values that Phi assumes in 
% the center of all the elements.
% The goal corresponds to the creation of an N times N sparse matrix DX
% 
% This task is decomposed in two steps:
%
% 1) The value that (the gradient of) Phi assumes on the center of each face k
% is evaluated from the value that Phi assumes in the center of each element
% "close" the face k. This could mean e.g. the simple average of the two
% elements having the face k on their boundaries (see Eq. 15) 
% or the Inverse Distance Weighted (IDW) Face Interpolation among all the elements 
% sharing an edge with the face k (see Eq. 16)
% 
% 2) The (second) derivative of Phi in the center of a given elment n is
% evaluated from the value that (the gradient of) Phi assumes on the center
% of each face that is on the boundary of the element n (see Eq. 14)
% 
% The first step is performed by a K times N sparse matrix W:    W*phi(elements) =  phi(faces)
% The second step is performed by N times K sparse matrix DDX: DDX*phi(faces)    = dphi(elements)
% The matrix DX is thus given by DDX*W: DX*phi(elements) = (DDX*W)*phi(elements) = dphi(elements)
% 
% INPUT ARGUMENTS:
%
% Volumes,Xel,Yel,Zel, are N times 1 arrays with the corresponding properties 
% of the N elements (volumes, and centers)
%
% Areas,NX,NY,Nz,Xf,Yf,Zf are K times 1 arrays with the corresponding properties
% of the K faces (areas, unitg normals, and centers). 
%
% T is a K times N sparse matrix which is used to create the matrix W.
% The entry T(k,n) is 1 if the n-th element shares at least one vertex
% with the k-th face 
%
% Signs is a N times K sparse matrix.
% The entry Signs(n,k) is equal to 1 if the normal to the k-th faces points outwards 
% with respect to the n-th element, and -1 if the normal points inwards.
% If the k-th face is not on the boundary of the n-th element Signs(n,k) is zero.
% 
% interpn is the interpolation scheme used to calculate (the gradient of)
% Phi.
%
% weight is the weighting scheme used in the interpolation. Any number is
% treated as the exponent of reverse distance weighting.
%
% method can be either Green Gauss (GG) or DirectLaplacian.
%
% 
%
% ExtW is an optional argument which can be a function-handle or a K times N sparse matrix.
%    if ExtW is not passed then IDW is used.
%    if ExtW is a function-handle, then it is used to calculate the weights from Xel,Yel,Zel, and Xf,Yf,Zf
%    if ExtW is a matrix, then it is directly used to compute the arithmetic mean of Phi on each face (T is not used) 
%
% OUTPUT ARGUMENTS:
% 
% DX,DY,DZ are N times N sparse matrices performing the derivatives
% W is the K times N sparse matrix performing the averages over the k-faces
%
%% Unpack the variables
NX = GridInfo.fNormX;
NY = GridInfo.fNormY;
NZ = GridInfo.fNormZ;
Areas = GridInfo.AreaFaces;
Volumes = GridInfo.Volumes;
Signs = GridInfo.TheSigns;
Xel = GridInfo.Xel;
Yel = GridInfo.Yel;
Zel = GridInfo.Zel;
Xf = GridInfo.Xf;
Yf = GridInfo.Yf;
Zf = GridInfo.Zf;
T = GridInfo.TheDs;
%if ~exist('Aexch','var')
%    Aexch=ones(size(Xel));
%end

%% Determine dimensionality
if numel(unique(GridInfo.Zel))-1
    dims = 3;
else
    dims = 2;
end
%% Dimensions
N = size(Signs,1) ;
K = size(Signs,2) ;
%% DDX, DDY, DDZ
VolCoeff = 1./Volumes ;
AX = NX.*Areas ;
AY = NY.*Areas ;
AZ = NZ.*Areas ;
s = nonzeros(Signs) ;
[n,k] = find(Signs) ;
ddx = s.*AX(k).*VolCoeff(n) ;
DDX = sparse(n,k,ddx,N,K) ;
ddy = s.*AY(k).*VolCoeff(n) ;
DDY = sparse(n,k,ddy,N,K) ;
if (dims > 2)
    ddz = s.*AZ(k).*VolCoeff(n) ;
    DDZ = sparse(n,k,ddz,N,K) ;
end

el2el=abs(Signs)*abs(Signs'); % compact (nn)

%% W
[n,k] = find(el2el) ;
if ~exist('Aexch','var') %normal method
    if (dims > 2)
        % either 1/dx^2, 1/dy^2 or 1/dz^2 since at all times only one of these
        % are nonzero.
        w =1./((Xel(n)-Xel(k))+(Yel(n)-Yel(k))+(Zel(n)-Zel(k))).^2;
    else
        % either 1/dx^2 or 1/dy^2 since at all times only one of these
        % are nonzero.
        w =1./((Xel(n)-Xel(k))+(Yel(n)-Yel(k))).^2;
    end
else % varying Aexch method
    if (dims > 2)
        % either 1/dx^2, 1/dy^2 or 1/dz^2 since at all times only one of these
        % are nonzero.
        w =Aexch(n)*2./(((Xel(n)-Xel(k))+(Yel(n)-Yel(k))+(Zel(n)-Zel(k))).^2.*(Aexch(n)+Aexch(k)));
    else
        % either 1/dx^2 or 1/dy^2 since at all times only one of these
        % are nonzero.
        w =Aexch(n)*2./(((Xel(n)-Xel(k))+(Yel(n)-Yel(k))).^2.*(Aexch(n)+Aexch(k)));
    end
end

inds2=[find(diff(k));length(n)];
inds1=[1;inds2+1];
vw=zeros(length(w),1);

%% Main loop
for kk=1:N
    ind=inds1(kk):inds2(kk);
    ii=find(n(ind)==k(ind)); % find central element (element kk)
    for i=1:length(ind)
        if i==ii
            continue
        else
            vw(ind(ii))= vw(ind(ii))- w(ind(i)); % central element
            vw(ind(i)) = vw(ind(i)) + w(ind(i)); % neighbour element
        end
    end
end

D2=sparse(k,n,vw);
end
