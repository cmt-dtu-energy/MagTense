function [DX,DY,DZ,W] = ComputeDifferentialOperatorsFromMesh03_GGDirLap(NX,NY,NZ,Areas,Volumes,Signs,Xel,Yel,Zel,Xf,Yf,Zf,T,interpn,weight,method)
%
% [DX,DY,DZ,W] = ComputeDifferentialOperatorsFromMesh03_LSQR(NX,NY,NZ,Areas,Volumes,Signs,X,Y,Z,Xf,Yf,Zf,T,WeightFunct)
%
% Reference: "Gradient Calculation Methods on Arbitrary Polyhedral Unstructured Meshes for Cell-Centered CFD Solvers"
%
% Creates differential operators for an unstructured mesh
% The mesh is composed by
% N elements (or tiles, or cells)
% K faces 
%
% The goal is to evaluate the x-derivative of a function Phi in the center of each element n
% from the values that Phi assumes in the center of all the elements.
% The goal corresponds to the creation of an N times N sparse matrix DX
% 
% This task is decomposed in two steps:
%
% 1) The value that Phi assumes on the center of each face k
% is evaluated from the value that Phi assumes in the center of each element
% is "close" the face k. This could mean e.g. the simple average of the two
% elements having the face k on their boundaries (see Eq. 15) 
% or the Inverse Distance Weighted (IDW) Face Interpolation among all the elements 
% sharing an edge with the face k (see Eq. 16)
% 
% 2) The x-derivative of Phi in the center of a given elment n 
% is evaluated from the value that Phi assumes on the center of each face
% that is on the boundary of the element n (see Eq. 14)
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
% of the K faces (areas, unit normals, and centers). 
%
% T is a K times N sparse matrix which is used to create the matrix W.
% The entry T(k,n) is 1 if the n-th element shares at least one edge
% with the k-th face 
%
% Signs is a N times K sparse matrix.
% The entry Signs(n,k) is equal to 1 if the normal to the k-th faces points outwards 
% with respect to the n-th element, and -1 if the normal points inwards.
% If the k-th face is not on the boundary of the n-th element Signs(n,k) is zero.
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
ddy = s.*AY(k).*VolCoeff(n) ;
ddz = s.*AZ(k).*VolCoeff(n) ;

DDX = sparse(n,k,ddx,N,K) ;
DDY = sparse(n,k,ddy,N,K) ;
DDZ = sparse(n,k,ddz,N,K) ;

if ~exist('interpn','var')
    interpn='extended';
end
if ~exist('weight','var') 
    weight=1;
elseif ~isnumeric(weight)
    warning('Unrecognized weight "%s". Using inverse distance.',weight)
    weight=1;
end
if ~exist('method','var')
    method="DirectLaplacian";
end
switch interpn
    case 'extended'
        el2fa=T'; % extended
    case 'compact'
        el2fa=Signs; % compact. Probably get's stuck in voronoi.
        warning('untested method "compact"')
    otherwise
        warning('unrecognized interpolation scheme "%s". Using extended scheme.',interpn)
        el2fa=T'; % extended
end

%% W

if ~exist('ExtW','var')
%     [k,n] = find(T) ;
    [n,k] = find(el2fa) ;
%     w = (abs(Xel(n)-Xf(k))+abs(Yel(n)-Yf(k))+abs(Zel(n)-Zf(k))).^(-1) ; % 1 entry for each "link"
    w = Volumes(n).*((Xel(n)-Xf(k)).^2+(Yel(n)-Yf(k)).^2+(Zel(n)-Zf(k)).^2).^(-weight/2) ; % 1 entry for each "link"
    W = sparse(k,n,w,K,N) ;
    % Normalization of W
%     ww = sum(W,2) ; % K entries array : Sum over all the elements contributing to the k-th face
%     w = w./ww(k) ;
%     W = sparse(k,n,w,K,N) ;
    % normalize distances
%     Xel=Xel/mean(abs(Xel));Yel=Yel/mean(abs(Yel));Zel=Zel/mean(abs(Zel));
%     Xf=Xf/mean(abs(Xf));Yf=Yf/mean(abs(Yf));Zf=Zf/mean(abs(Zf));
    
%     Xel=Xel/mean(abs(Xel));Yel=Yel/mean(abs(Yel));Zel=Zel/mean(abs(Zel));
%     Xf=Xf/mean(abs(Xf));Yf=Yf/mean(abs(Yf));Zf=Zf/mean(abs(Zf));
%     XEl=repmat(Xel',K,1);
%     YEl=repmat(Yel',K,1);
%     ZEl=repmat(Zel',K,1);
%     XF=repmat(Xf,1,N);
%     YF=repmat(Yf,1,N);
%     ZF=repmat(Zf,1,N);
%     W = (abs(XEl-XF)+abs(YEl-YF)+abs(ZEl-ZF)).^(-1) ;
%     Wnorm=sum(W,'all');
%     W = W/Wnorm;
%     Xel=Xel*Wnorm;Yel=Yel*Wnorm;Zel=Zel*Wnorm;
%     Xf=Xf*Wnorm;Yf=Yf*Wnorm;Zf=Zf*Wnorm;
    % w = ((Y(n)==Yf(k)) & (Z(n)==Zf(k))) | ((Z(n)==Zf(k)) & (X(n)==Xf(k))) | ((X(n)==Xf(k)) & (Y(n)==Yf(k))) ;
else
    if isequal(class(ExtW),'function_handle')
        [k,n] = find(T) ;
        w = ExtW(Xel(n),Yel(n),Zel(n),Xf(k),Yf(k),Zf(k)) ;
        W = sparse(k,n,w,K,N) ;
    else
        [k,n] = find(ExtW) ;
        w =  nonzeros(ExtW) ;
        W =  ExtW ;
    end
end


%% Normalization of W
% 
% ww = sum(W,2) ; % K entries array : Sum over all the elements contributing to the k-th face
% w = w./ww(k) ;
% 
% W = sparse(k,n,w,K,N) ;
vw=zeros(length(w),1);
vx=zeros(length(w),1);vy=zeros(length(w),1);vz=zeros(length(w),1);
counter=0;
warning('off','MATLAB:nearlySingularMatrix') % 
for kk=1:K
    flt=~(k-kk);
%     Gk=[sum(W(kk,n(flt)))                , W(kk,n(flt))*(Xel(n(flt))-Xf(kk))                       , W(kk,n(flt))*(Yel(n(flt))-Yf(kk))                        , W(kk,n(flt))*(Zel(n(flt))-Zf(kk));
%         W(kk,n(flt))*(Xel(n(flt))-Xf(kk)), W(kk,n(flt))*(((Xel(n(flt))-Xf(kk))).^2)                , W(kk,n(flt))*((Xel(n(flt))-Xf(kk)).*(Yel(n(flt))-Yf(kk))), W(kk,n(flt))*((Xel(n(flt))-Xf(kk)).*(Zel(n(flt))-Zf(kk)));
%         W(kk,n(flt))*(Yel(n(flt))-Yf(kk)), W(kk,n(flt))*((Xel(n(flt))-Xf(kk)).*(Yel(n(flt))-Yf(kk))), W(kk,n(flt))*((Yel(n(flt))-Yf(kk)).^2)                  , W(kk,n(flt))*((Yel(n(flt))-Yf(kk)).*(Zel(n(flt))-Zf(kk)));
%         W(kk,n(flt))*(Zel(n(flt))-Zf(kk)), W(kk,n(flt))*((Xel(n(flt))-Xf(kk)).*(Zel(n(flt))-Zf(kk))), W(kk,n(flt))*((Yel(n(flt))-Yf(kk)).*(Zel(n(flt))-Zf(kk))), W(kk,n(flt))*((Zel(n(flt))-Zf(kk)).^2)];
%     Hk=[W(kk,n(flt));
%         W(kk,n(flt)).*((Xel(n(flt))-Xf(kk))');
%         W(kk,n(flt)).*((Yel(n(flt))-Yf(kk))');
%         W(kk,n(flt)).*((Zel(n(flt))-Zf(kk))')];
%     Wkfull=Gk\Hk;
    nk=n(flt);
    dxk=Xel(nk)-Xf(kk);
    dyk=Yel(nk)-Yf(kk);
    dzk=Zel(nk)-Zf(kk);
    nl=length(nk);
    e=ones(nl,1);
    Wk=W(kk,nk);
    Gkl1=[e, dxk, dyk, dzk];
    Gk=[Wk*Gkl1;
        Wk*(dxk.*Gkl1);
        Wk*(dyk.*Gkl1);
        Wk*(dzk.*Gkl1)];
    Hk=Wk.*Gkl1';
    Wkfull=Gk\Hk; 
%     W(kk,n(flt))=Wkfull(1,:);
    vw(counter+(1:nl))=Wkfull(1,:);
    vx(counter+(1:nl))=Wkfull(2,:);
    vy(counter+(1:nl))=Wkfull(3,:);
    vz(counter+(1:nl))=Wkfull(4,:);
    counter=counter+nl;
end
warning('on','MATLAB:nearlySingularMatrix') % 
% ww = sum(W,2) ; % K entries array : Sum over all the elements contributing to the k-th face
% [k,n] = find(W) ;   
% w = nonzeros(W) ;
% w = w./ww(k) ;

if method=="GG"
    W = sparse(k,n,vw,K,N) ;
    %% DX, DY, DZ
    DX = DDX*W ;
    DY = DDY*W ;
    DZ = DDZ*W ;
elseif method=="DirectLaplacian"
    FX = sparse(k,n,vx,K,N) ;
    FY = sparse(k,n,vy,K,N) ;
    FZ = sparse(k,n,vz,K,N) ;
    %% DX, DY, DZ
    DX = DDX*FX ;
    DY = DDY*FY ;
    DZ = DDZ*FZ ;
else
    error('Unrecognized method "%s".', method)
end
