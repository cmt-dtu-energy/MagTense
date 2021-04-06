function [DX,DY,DZ,W] = computeDifferentialOperatorsFromMesh_DirectLap(GridInfo,interpn,weight,method)
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

%% Determine dimensionality
if (numel(unique(GridInfo.Zel)) == 1)
    dims = 2;
else
    dims = 3;
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

if ~exist('interpn','var')
    interpn='extended';
end
if ~exist('weight','var') 
    weight=1;
elseif ~isnumeric(weight)
    try
        weight=str2double(weight);
    catch ME
        warning('Unrecognized weight "%g".',weight)
        rethrow(ME)
    end
end
if ~exist('method','var')
    method="DirectLaplacian";
end
switch interpn
    case 'extended'
        el2fa=T'; % extended
    case 'compact'
        el2fa=Signs; % compact. Probably gets stuck in voronoi.
        warning('untested method "compact"')
    otherwise
        warning('unrecognized interpolation scheme "%s". Using extended scheme.',interpn)
        el2fa=T'; % extended
end

%% W
edgeElWeight=2;
if ~exist('ExtW','var')
%     [k,n] = find(T) ;
    [n,k] = find(el2fa) ;
    % Edge cases. Use adjacent cell as center
    ke = find(sum(abs(Signs),1)==1);
% %     rem=ismember(k,ke);
% %     n(rem)=[];k(rem)=[];
%     w = (abs(Xel(n)-Xf(k))+abs(Yel(n)-Yf(k))+abs(Zel(n)-Zf(k))).^(-1) ; % 1 entry for each "link"
%     w = Volumes(n).*((Xel(n)-Xf(k)).^2+(Yel(n)-Yf(k)).^2+(Zel(n)-Zf(k)).^2).^(-weight/2) ; % 1 entry for each "link"
if (dims == 2)
    w = ((Xel(n)-Xf(k)).^2+(Yel(n)-Yf(k)).^2).^(-weight/2) ; % 1 entry for each "link"
    else
    w = ((Xel(n)-Xf(k)).^2+(Yel(n)-Yf(k)).^2+(Zel(n)-Zf(k)).^2).^(-weight/2) ; % 1 entry for each "link"
end
% %     [ne,~] = find(abs(Signs(:,ke)));
% %     el2el = abs(Signs)*abs(el2fa');
% %     [nee,tmp]=find(el2el(ne,:)');
% %     kee=ke(tmp);
% %     nee2=ne(tmp);
% %     we = ((Xel(nee2)-Xel(nee)).^2+(Yel(nee2)-Yel(nee)).^2+(Zel(nee2)-Zel(nee)).^2).^(-weight/2) ; % 1 entry for each "link"
% %     total=sortrows([k,n,w;
% %            kee',nee,we]);
% %     k=total(:,1);n=total(:,2);w=total(:,3);
    W = sparse(k,n,w,K,N) ;
    '';

    % Normalization of W    !! Does not change anything.
%     ww = sum(W,2) ; % K entries array : Sum over all the elements contributing to the k-th face
%     w = w./ww(k) ;
%     W = sparse(k,n,w,K,N) ;
    % normalize distances   !! Does not work.
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
inds2=[find(diff(k));length(n)];
inds1=[1;inds2+1];
dx=Xel(n)-Xf(k);
dy=Yel(n)-Yf(k);
vw=zeros(length(w),1);
vx=zeros(length(w),1);
vy=zeros(length(w),1);
if (dims > 2)
    dz=Zel(n)-Zf(k);
    vz=zeros(length(w),1);
end

warning('off','MATLAB:nearlySingularMatrix') % 
for kk=1:K
%     flt=~(k-kk);
%     nk=n(flt);
%     dxk=Xel(nk)-Xf(kk);
%     dyk=Yel(nk)-Yf(kk);
%     dzk=Zel(nk)-Zf(kk);
%     nl=length(nk);
%     Wk=W(kk,nk);
%         wadje = edgeElWeight*((Xel(nek)-Xf(kk)).^2+(Yel(nek)-Yf(kk)).^2+(Zel(nek)-Zf(kk)).^2).^(-weight/2); % weight of the face-adjoining element
    ind=inds1(kk):inds2(kk);
    Wk=w(ind)';
% % 	Wk(isinf(Wk))=edgeElWeight*max(Wk(~isinf(Wk))); % weight of the face-adjoining element
    dxk=dx(ind);
    dyk=dy(ind);
    if (dims > 2)
        dzk=dz(ind);
    end
    if sum(abs(Signs(:,kk)))==1 % edge face. Mirror trick to enforce Neumann b.c.
        lind=length(ind);
        e=ones(2*lind,1);
        if (dims == 2)
            nns=[NX(kk),NY(kk)];
            extra=[dxk,dyk]-2*nns.*(sum([dxk,dyk],2)/sum(nns)); % mirror reflection on face plane
        else
            nns=[NX(kk),NY(kk),NZ(kk)];
            extra=[dxk,dyk,dzk]-2*nns.*(sum([dxk,dyk,dzk],2)/sum(nns)); % mirror reflection on face plane
        end
        dxk=[dxk;extra(:,1)];
        dyk=[dyk;extra(:,2)];
        if (dims > 2)
            dzk=[dzk;extra(:,3)];
        end
        Wk=[Wk Wk];
        if (dims == 2)
            Gkl1=[e, dxk, dyk];
            Gk=[Wk*Gkl1; Wk*(dxk.*Gkl1); Wk*(dyk.*Gkl1)];
        else
            Gkl1=[e, dxk, dyk, dzk];
            Gk=[Wk*Gkl1; Wk*(dxk.*Gkl1); Wk*(dyk.*Gkl1); Wk*(dzk.*Gkl1)];
        end
        Hk=Wk.*Gkl1';
    %     Wkfull=Gk\Hk; 
    try
        R=chol(Gk);% Speed upgrade, Gk is positive definite and symmetric
        Wkfull=R\(R'\Hk);
    catch errmsg
        Wkfull=Gk\Hk; % Speed downgrade, most likely a 2D mesh/tile
    end
    %     W(kk,n(flt))=Wkfull(1,:);
        vw(ind)=Wkfull(1,1:lind)+Wkfull(1,lind+1:end);
        vx(ind)=Wkfull(2,1:lind)+Wkfull(2,lind+1:end);
        vy(ind)=Wkfull(3,1:lind)+Wkfull(3,lind+1:end);
        if (dims > 2)
            vz(ind)=Wkfull(4,1:lind)+Wkfull(4,lind+1:end);
        end
    else
        e=ones(length(ind),1);
        if (dims == 2)
            Gkl1=[e, dxk, dyk];
            Gk=[Wk*Gkl1; Wk*(dxk.*Gkl1); Wk*(dyk.*Gkl1)];
        else
            Gkl1=[e, dxk, dyk, dzk];
            Gk=[Wk*Gkl1; Wk*(dxk.*Gkl1); Wk*(dyk.*Gkl1); Wk*(dzk.*Gkl1)];
        end
        Hk=Wk.*Gkl1';
    try
        R=chol(Gk);% Speed upgrade, Gk is positive definite and symmetric
        Wkfull=R\(R'\Hk);
    catch errmsg
        Wkfull=Gk\Hk; % Speed downgrade, most likely a 2D mesh/tile
    end
    %     W(kk,n(flt))=Wkfull(1,:);
        vw(ind)=Wkfull(1,:);
        vx(ind)=Wkfull(2,:);
        vy(ind)=Wkfull(3,:);
        if (dims > 2)
            vz(ind)=Wkfull(4,:);
        end
    end
end
warning('on','MATLAB:nearlySingularMatrix') % 
% ww = sum(W,2) ; % K entries array : Sum over all the elements contributing to the k-th face
% [k,n] = find(W) ;   
% w = nonzeros(W) ;
% w = w./ww(k) ;

if method=="GGNeumann"
    W = sparse(k,n,vw,K,N) ;
    %% DX, DY, DZ
    DX = DDX*W ;
    DY = DDY*W ;
    if (dims > 2)
        DZ = DDZ*W ;
    end
elseif method=="DirectLaplacianNeumann"
    %% DX, DY, DZ
    FX = sparse(k,n,vx,K,N) ;
    DX = DDX*FX ;
    FY = sparse(k,n,vy,K,N) ;
    DY = DDY*FY ;
    if (dims > 2)
        FZ = sparse(k,n,vz,K,N) ;
        DZ = DDZ*FZ ;
    end
else
    error('Unrecognized method "%s".', method)
end
