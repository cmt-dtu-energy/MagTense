function [DX,DY,DZ,W] = ComputeDifferentialOperatorsFromMesh_Neumann04_GGDirLap(GridInfo,interpn,weight,method,Aexch)
%
% [DX,DY,DZ,W] = ComputeDifferentialOperatorsFromMesh04_Neumann04_GGDirLap(GridInfo,interpn,weight,method,Aexch)
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
if ~exist('Aexch','var')
    Aexch=ones(size(Xel));
end

%% Determine dimensionality
if numel(unique(GridInfo.Zel))-1
    dims = 3;
elseif numel(unique(GridInfo.Yel))-1
    dims = 2;
else
    dims = 1;
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
Amat=zeros(length(n),1);
for kk=1:length(k)
    tmp=find(Signs(:,k(kk)));
    tmp2=tmp(tmp~=n(kk));
    % Can be this or 0. Should be 0, honestly.
    if isempty(tmp2)
        Amat(kk)=Aexch(n(kk));
        continue
    elseif length(tmp2)>1
        warning('more than two cells share a face!')
        disp(tmp2')
        disp(' share face ')
        disp(k(kk))
    else
        Amat(kk)=2*Aexch(n(kk))*Aexch(tmp2)/(Aexch(n(kk))+Aexch(tmp2));
    end
end
if method~="Multimaterial"
    ddx = s.*AX(k).*VolCoeff(n) ;
    DDX = sparse(n,k,ddx,N,K) ;
    if (dims > 1)
        ddy = s.*AY(k).*VolCoeff(n) ;
        DDY = sparse(n,k,ddy,N,K) ;
        if (dims > 2)
            ddz = s.*AZ(k).*VolCoeff(n) ;
            DDZ = sparse(n,k,ddz,N,K) ;
        end
    end
else
    ddxA = s.*AX(k).*Amat.*VolCoeff(n) ;
    DDXA = sparse(n,k,ddxA,N,K) ;
    if (dims > 1)
        ddyA = s.*AY(k).*Amat.*VolCoeff(n) ;
        DDYA = sparse(n,k,ddyA,N,K) ;
        if (dims > 2)
            ddzA = s.*AZ(k).*Amat.*VolCoeff(n) ;
            DDZA = sparse(n,k,ddzA,N,K) ;
        end
    end
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
if ~exist('ExtW','var')
%     [k,n] = find(T) ;
    [n,k] = find(el2fa) ;
%     w = (abs(Xel(n)-Xf(k))+abs(Yel(n)-Yf(k))+abs(Zel(n)-Zf(k))).^(-1) ; % 1 entry for each "link"
%     w = Volumes(n).*((Xel(n)-Xf(k)).^2+(Yel(n)-Yf(k)).^2+(Zel(n)-Zf(k)).^2).^(-weight/2) ; % 1 entry for each "link"
if (dims == 1)
    w = ((Xel(n)-Xf(k)).^2).^(-weight/2) ; % 1 entry for each "link"
elseif (dims == 2)
    w = ((Xel(n)-Xf(k)).^2+(Yel(n)-Yf(k)).^2).^(-weight/2) ; % 1 entry for each "link"
else
    w = ((Xel(n)-Xf(k)).^2+(Yel(n)-Yf(k)).^2+(Zel(n)-Zf(k)).^2).^(-weight/2) ; % 1 entry for each "link"
end
    W = sparse(k,n,w,K,N) ;
    '';
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

inds2=[find(diff(k));length(n)];
inds1=[1;inds2+1];
dx=Xel(n)-Xf(k);

vw=zeros(length(w),1);
vx=zeros(length(w),1);

if (dims > 1)
    dy=Yel(n)-Yf(k);
    vy=zeros(length(w),1);
    if (dims > 2)
        dz=Zel(n)-Zf(k);
        vz=zeros(length(w),1);
    end
end
% Scale weights to avoid ill conditioning
wm=max(w(~isinf(w)));
w=w*(wm^(1/weight))/wm;
% dist_scale=mean(abs([dx;dy;dz]));
% dx=dx/dist_scale;
% dy=dy/dist_scale;
% dz=dz/dist_scale;

% warning('off','MATLAB:nearlySingularMatrix') % 
counter=0;
for kk=1:K
    ind=inds1(kk):inds2(kk);
    Wk=w(ind)';
    dxk=dx(ind);
    scale=10^(round(log10(mean(abs(dxk)))));
    if (dims > 1)
        dyk=dy(ind);
        dks=[dxk;dyk];
        scale=10^(round(log10(mean(abs(dks)))));
        if (dims > 2)
            dzk=dz(ind);
            dks=[dks;dzk];
            scale=10^(round(log10(mean(abs(dks)))));
            dzk=dzk/scale;% Scale distances to avoid ill conditioning
        end
        dyk=dyk/scale;% Scale distances to avoid ill conditioning
    end
    dxk=dxk/scale;% Scale distances to avoid ill conditioning
    if sum(abs(Signs(:,kk)))==1 % edge face. Mirror trick to enforce Neumann b.c.
        counter=counter+1;
        lind=length(ind);
        e=ones(2*lind,1);
        if (dims == 1)
            nns=NX(kk);
            if all(nns==0) % z- or y-face. Nothing to see here.
                continue
            end
            extra=dxk-2*nns.*(dxk*nns'); % mirror reflection on face plane
        elseif (dims == 2)
            nns=[NX(kk),NY(kk)];
            if all(nns==[0,0]) % z-face. Nothing to see here.
                continue
            end
            extra=[dxk,dyk]-2*nns.*([dxk,dyk]*nns'); % mirror reflection on face plane
        else
            nns=[NX(kk),NY(kk),NZ(kk)];
            extra=[dxk,dyk,dzk]-2*nns.*([dxk,dyk,dzk]*nns'); % mirror reflection on face plane
        end
        dxk=[dxk;extra(:,1)];
        if (dims > 1)
            dyk=[dyk;extra(:,2)];
            if (dims > 2)
                dzk=[dzk;extra(:,3)];
            end
        end
        Wk=[Wk Wk];
        if (dims == 1 || abs(NX(kk))==1) % hack. What about NY, NZ cases?)
            Gkl1=[e, dxk];
            Gk=[Wk*Gkl1; Wk*(dxk.*Gkl1)];
        elseif (dims == 2)
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
        Wkfull=[Wkfull(1,:);Wkfull(2:end,:)/scale]; % Scale distances to avoid ill conditioning
        vw(ind)=Wkfull(1,1:lind)+Wkfull(1,lind+1:end);
        vx(ind)=Wkfull(2,1:lind)+Wkfull(2,lind+1:end);
        if (dims > 1 && abs(NX(kk))~=1) % hack. What about NY, NZ cases?
            vy(ind)=Wkfull(3,1:lind)+Wkfull(3,lind+1:end);
            if (dims > 2)
                vz(ind)=Wkfull(4,1:lind)+Wkfull(4,lind+1:end);
            end
        end
    else
        e=ones(length(ind),1);
        if (dims == 1 || abs(NX(kk))==1) % hack. What about NY, NZ cases?
            Gkl1=[e, dxk];
            Gk=[Wk*Gkl1; Wk*(dxk.*Gkl1)];
        elseif (dims == 2)
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
        Wkfull=[Wkfull(1,:);Wkfull(2:end,:)/scale];% Scale distances to avoid ill conditioning
        vw(ind)=Wkfull(1,:);
        vx(ind)=Wkfull(2,:);
        if (dims > 1 && abs(NX(kk))~=1) % hack. What about NY, NZ cases?
            vy(ind)=Wkfull(3,:);
            if (dims > 2)
                vz(ind)=Wkfull(4,:);
            end
        end
    end
end
% warning('on','MATLAB:nearlySingularMatrix') % 

if method=="GGNeumann"
    W = sparse(k,n,vw,K,N) ;
    %% DX, DY, DZ
    DX = DDX*W ;
    if (dims > 1)
        DY = DDY*W ;
        if (dims > 2)
            DZ = DDZ*W ;
        end
    end
elseif method=="DirectLaplacianNeumann"
    %% DX, DY, DZ
    FX = sparse(k,n,vx,K,N) ;
    DX = DDX*FX ;
    if (dims > 1)
        FY = sparse(k,n,vy,K,N) ;
        DY = DDY*FY ;
        if (dims > 2)
            FZ = sparse(k,n,vz,K,N) ;
            DZ = DDZ*FZ ;
        end
    end
elseif method=="Multimaterial"
    %% DX, DY, DZ
    FX = sparse(k,n,vx,K,N) ;
    DX = DDXA*FX ;
    if (dims > 1)
        FY = sparse(k,n,vy,K,N) ;
        DY = DDYA*FY ;
        if (dims > 2)
            FZ = sparse(k,n,vz,K,N) ;
            DZ = DDZA*FZ ;
        end
    end
else
    error('Unrecognized method "%s".', method)
end
