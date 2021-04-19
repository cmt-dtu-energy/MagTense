function [DX,DY,DZ,W] = ComputeDifferentialOperatorsFromMesh03_LSQRgrad6(NX,NY,NZ,Areas,Volumes,Signs,Xel,Yel,Zel,Xf,Yf,Zf,T,interpn,weight,method)
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
    weight="1";
end
if ~exist('method','var')
    method="LSQR";
end
switch interpn
    case 'extended'
        el2el=abs(Signs)*abs(T); % extended
    case 'compact'
        el2el=abs(Signs)*abs(Signs'); % compact (nn)
    case 'fn2'
        el2el1=abs(Signs)*abs(Signs'); % compact (nn)
        el2el=el2el1*el2el1; % fn2 (nn + nnn)
    otherwise
        warning('unrecognized interpolation scheme "%s". Using extended scheme.',interpn)
        el2el=abs(Signs)*abs(T); % extended
end

%% W

if ~exist('ExtW','var')
    [n,k] = find(el2el) ;
%     w = (abs(Xel(n)-Xf(k))+abs(Yel(n)-Yf(k))+abs(Zel(n)-Zf(k))).^(-1) ; % 1 entry for each "link"
    sames=n==k;
    n(sames)=[];
    k(sames)=[];
    if ~isempty(str2num(weight))
%     w =(Volumes(k)-Xel(n)*0).*((Xel(n)-Xel(k)).^2+(Yel(n)-Yel(k)).^2+(Zel(n)-Zel(k)).^2).^(-3/4); % Syrakos (2017) magic weight
        w =(Volumes(k)-Xel(n)*0).*((Xel(n)-Xel(k)).^2+(Yel(n)-Yel(k)).^2+(Zel(n)-Zel(k)).^2).^(-str2num(weight)/2);
    elseif weight=="syrakos"
        w =(Volumes(k)-Xel(n)*0).*((Xel(n)-Xel(k)).^2+(Yel(n)-Yel(k)).^2+(Zel(n)-Zel(k)).^2).^(-3/4); % Syrakos (2017) magic weight
    elseif weight=="shima"
% %     % Find faces between elements of n and elements of k
    tmp=Signs(n,:)&Signs(k,:);
    [JJ,II]=find(tmp');
%     (Volumes(k)-Xel(n)*0).*
    if norm(NX.*NY.*NZ)
        disp('tet')
    else
        disp('vor')
    end
    if length(JJ)~=length(n) % doing f2n. Use both faces.
        nII=setdiff(1:length(n),II)';
        
        JJ2=[];
        % Find common neighbour(s) between next-nearest neighbours.
        tmp1=el2el1(n(nII),:)&el2el1(k(nII),:);
        [JJ1,II1]=find(tmp1');
        filt=[~~diff(II1);true];
        JJ1x(1:sum(filt))=JJ1(filt); % Sometimes there are more common neighbours. If so, pick the last one.
        tmp=Signs(JJ1x,:)&Signs(k(nII),:); % Faces between 1st neighbour pairs
        [JJ2(:,1),~]=find(tmp');
        % b: Only use 1st neighbour pairs. a:
        tmp=Signs(JJ1x,:)&Signs(k(nII),:); % Faces between 2nd neighbour pairs
        [JJ2(:,2),~]=find(tmp');
        
%         w=[];
        w=zeros(length(n),1);
        % Weights for nearest neighbours
        w(II) = 4.*(...
        (((Xf(JJ)-Xel(k(II))).*NX(JJ)).^2+((Yf(JJ)-Yel(k(II))).*NY(JJ)).^2+((Zf(JJ)-Zel(k(II))).*NZ(JJ)).^2)./...%l'
        (((Xel(n(II))-Xel(k(II))).*NX(JJ)).^2+((Yel(n(II))-Yel(k(II))).*NY(JJ)).^2+((Zel(n(II))-Zel(k(II))).*NZ(JJ)).^2))...%L'
        .*(Areas(JJ))...% s
        .*((Xel(n(II))-Xel(k(II))).^2+(Yel(n(II))-Yel(k(II))).^2+(Zel(n(II))-Zel(k(II))).^2).^(-1/2) ;% L. % Shima (2013) magic weight
        % Weights for next-nearest neighbours
        w(nII) = sum(4.*(...
        (((Xf(JJ2)-Xel(k(nII))).*NX(JJ2)).^2+((Yf(JJ2)-Yel(k(nII))).*NY(JJ2)).^2+((Zf(JJ2)-Zel(k(nII))).*NZ(JJ2)).^2)./...
        (((Xel(n(nII))-Xel(k(nII))).*NX(JJ2)).^2+((Yel(n(nII))-Yel(k(nII))).*NY(JJ2)).^2+((Zel(n(nII))-Zel(k(nII))).*NZ(JJ2)).^2))...
        .*(Areas(JJ2)/size(JJ2,2))...                                                                                                  % s
        .*((Xel(n(nII))-Xel(k(nII))).^2+(Yel(n(nII))-Yel(k(nII))).^2+(Zel(n(nII))-Zel(k(nII))).^2).^(-1/2)  ,  2) ; % Shima (2013) magic weight
    else
    w = 4.*(...
        (((Xf(JJ)-Xel(k)).*NX(JJ)).^2+((Yf(JJ)-Yel(k)).*NY(JJ)).^2+((Zf(JJ)-Zel(k)).*NZ(JJ)).^2)./...
        (((Xel(n)-Xel(k)).*NX(JJ)).^2+((Yel(n)-Yel(k)).*NY(JJ)).^2+((Zel(n)-Zel(k)).*NZ(JJ)).^2))...
        .*(Areas(JJ))...                                                                                % s
        .*((Xel(n)-Xel(k)).^2+(Yel(n)-Yel(k)).^2+(Zel(n)-Zel(k)).^2).^(-1/2) ; % Shima (2013) magic weight
    end
    else
        warning('unrecognized weighting scheme "%s". Using inverse distance.',weight)
        w =(Volumes(k)-Xel(n)*0).*((Xel(n)-Xel(k)).^2+(Yel(n)-Yel(k)).^2+(Zel(n)-Zel(k)).^2).^(-1/2);
    end
%     w = Volumes(k).*4.*(...
%         (((Xf(Signs(n,:)&Signs(k,:))-Xel(k)).*NX(Signs(n,:)&Signs(k,:))).^2+((Yf(Signs(n,:)&Signs(k,:))-Yel(k)).*NY(Signs(n,:)&Signs(k,:))).^2+((Zf(Signs(n,:)&Signs(k,:))-Zel(k)).*NZ(Signs(n,:)&Signs(k,:))).^2)./...
%         (((Xel(n)-Xel(k)).*NX(Signs(n,:)&Signs(k,:))).^2+((Yel(n)-Yel(k)).*NY(Signs(n,:)&Signs(k,:))).^2+((Zel(n)-Zel(k)).*NZ(Signs(n,:)&Signs(k,:))).^2))...
%         .*((Xel(n)-Xel(k)).^2+(Yel(n)-Yel(k)).^2+(Zel(n)-Zel(k)).^2).^(-1/2) ; % Shima (2013) magic weight
%     w = (abs(Xel(n)-Xf(k)).^10+abs(Yel(n)-Yf(k)).^10+abs(Zel(n)-Zf(k)).^10).^(-1/10) ; % 1 entry for each "link"
else
    if isequal(class(ExtW),'function_handle')
        [n,k] = find(T') ;
        w = ExtW(Xel(n),Yel(n),Zel(n),Xf(k),Yf(k),Zf(k)) ;
        W = sparse(k,n,w,K,N) ;
    else
        [n,k] = find(ExtW') ;
        w =  nonzeros(ExtW') ;
        W =  ExtW ;
    end
end
    dx = (Xel(n)-Xel(k));
    dy = (Yel(n)-Yel(k));
    dz = (Zel(n)-Zel(k));
    inds2=[find(diff(k));length(n)];
    inds1=[1;inds2+1];
% %     if length(JJ)~=length(n)
% %         JJall(II,:)=repmat(JJ,1,size(JJ2,2));
% %         JJall(nII,:)=JJ2;
% %     else
% %         JJall=JJ;
% %     end
%% Normalization of W
% 
% ww = sum(W,2) ; % K entries array : Sum over all the elements contributing to the k-th face
% w = w./ww(k) ;
% 
% W = sparse(k,n,w,K,N) ;
% FDX=sparse(K,N);
% FDY=sparse(K,N);
% FDZ=sparse(K,N);
% DX = spalloc(N,N,length(w));
% DY = spalloc(N,N,length(w));
% DZ = spalloc(N,N,length(w));
vx=zeros(length(w),1);vy=zeros(length(w),1);vz=zeros(length(w),1);
warning('off','MATLAB:nearlySingularMatrix') % 
% diffnorm=sqrt(sqrt(mean(dx.^2+dy.^2+dz.^2)));
diffnorm=1;
% w=w/sum(w);
low_rcond=1;
gg=0;
Cimaxmax=0;
for kk=1:N
    ind=inds1(kk):inds2(kk);
    wk=w(ind)';
%     wk(isinf(wk))=sum(wk(~isinf(wk))); % make 0-cell value count half
    wk(isinf(wk))=max(wk(~isinf(wk))); % b: make 0-cell value count at most half
%     wk(isinf(wk))=0; % c: make 0-cell value not count
%     wk=wk/sum(wk);
%     wk=0*wk+1; % No weights
    dxk=dx(ind)/diffnorm;dyk=dy(ind)/diffnorm;dzk=dz(ind)/diffnorm;
    e=ones(length(ind),1);
%     Gkl1=[e, dxk, dyk, dzk,(0.5*dxk.^2),(dxk.*dyk),(dxk.*dzk),(0.5*dyk.^2),(dyk.*dzk), (0.5*dzk.^2)];
%     Gk=[wk*Gkl1;
%         wk*(dxk.*Gkl1);
%         wk*(dyk.*Gkl1);
%         wk*(dzk.*Gkl1);
%         wk*(0.5*dxk.^2.*Gkl1);
%         wk*(dxk.*dyk.*Gkl1);
%         wk*(dxk.*dzk.*Gkl1);
%         wk*(0.5*dyk.^2.*Gkl1);
%         wk*(dyk.*dzk.*Gkl1);
%         wk*(0.5*dzk.^2.*Gkl1)];
%     Gkl1=[e, dxk, dyk, dzk];
%     Gk=[wk*Gkl1;
%         wk*(dxk.*Gkl1);
%         wk*(dyk.*Gkl1);
%         wk*(dzk.*Gkl1)];
%     xind=2;yind=3;zind=4;
    Gkl1=[dxk, dyk, dzk];
    Gk=[wk*(dxk.*Gkl1);
        wk*(dyk.*Gkl1);
        wk*(dzk.*Gkl1)];
    xind=1;yind=2;zind=3;
    if isnan(rcond(Gk))
        disp('fuck')
    end
    if rcond(Gk) == 0
        if false
        Gkl1=[e, dyk, dzk];
        Gk=[wk*Gkl1;
        wk*(dyk.*Gkl1);
        wk*(dzk.*Gkl1)];
        xind=4;yind=2;zind=3;
        if rcond(Gk) == 0
            Gkl1=[e, dxk, dzk];
            Gk=[wk*Gkl1;
            wk*(dxk.*Gkl1);
            wk*(dzk.*Gkl1)];
            xind=2;yind=4;zind=3;
            if rcond(Gk) == 0
                Gkl1=[e, dxk, dyk];
                Gk=[wk*Gkl1;
                    wk*(dxk.*Gkl1);
                    wk*(dyk.*Gkl1)];
                xind=2;yind=3;zind=4;
                if rcond(Gk) == 0
                    disp('I give up')
                    fprintf('Minimum neighbours: %d\nMaximum neighbours: %d\n',full(min(sum(el2el)-diag(el2el))),full(max(sum(el2el)-diag(el2el))))
                    Gkl1=[e, dxk, dyk, dzk];
                    Gk=[wk*Gkl1;
                        wk*(dxk.*Gkl1);
                        wk*(dyk.*Gkl1);
                        wk*(dzk.*Gkl1)];
                    xind=2;yind=3;zind=4;
                end
            end
        end
        end
    Hk=wk.*Gkl1';
    [Q,R]=qr(Gk);
    Wkfull=R\Q'*Hk;
    else
    Hk=wk.*Gkl1';
    
    Wkfull=Gk\Hk;
    end
    %}

% %     % Check monotonicity condition (Shima, 2013)
% %     faceind=JJall(ind);
% %     Xfk=Xf(faceind);Yfk=Yf(faceind);Zfk=Zf(faceind);
% %     Cimax=max(abs(sum((Xel(kk)-Xfk)'.*Wkfull(xind,:)+(Yel(kk)-Yfk)'.*Wkfull(yind,:)+(Zel(kk)-Zfk)'.*Wkfull(zind,:),1)/size(JJall,2)));
% %     Cimaxmax=max(Cimaxmax,Cimax);
% %     if Cimax>1
% % %         disp(Cimax)
% %         gg=gg+1;
% %         Cimaxmax=max(Cimaxmax,Cimax);
% %     end
    %
%     [wkc,wkr]=size(Wkfull);
%     if wkc == 3
%         Wkfull(4,:)=zeros(1,wkr);
%     end
%     if rcond(Gk) < low_rcond
%         low_rcond=rcond(Gk);
%     end
    % Stop-gap measure for when dxk, dyk or dzk are all-zero:
%     Wkfull(isnan(Wkfull))=0;
    %}
%     if any(isnan(Wkfull))
%         Wkfull(:,:)=0;
%     end
%     W(kk,n(flt))=Wkfull(1,:);
    vx(ind)=Wkfull(xind,:);
    vy(ind)=Wkfull(yind,:);
    vz(ind)=Wkfull(zind,:);
%     DX(kk,n(ind))=Wkfull(xind,:);
%     DY(kk,n(ind))=Wkfull(yind,:);
%     DZ(kk,n(ind))=Wkfull(zind,:);
end
DX=sparse(k,n,vx);
DY=sparse(k,n,vy);
DZ=sparse(k,n,vz);
warning('on','MATLAB:nearlySingularMatrix') % 
% % fprintf('%d GG out of %d. Cij_max: %g\r\n',gg,N,Cimaxmax)
% ww = sum(W,2) ; % K entries array : Sum over all the elements contributing to the k-th face
% [k,n] = find(W) ;   
% w = nonzeros(W) ;
% w = w./ww(k) ;

if method~="LSQR"
    error('Unrecognized method %s.',method)
end

% W = sparse(k,n,w,K,N) ;
% DX = DDX*FDX ;
% DY = DDY*FDY ;
% DZ = DDZ*FDZ ;
%% DX, DY, DZ
% DX = DDX*W ;
% DY = DDY*W ;
% DZ = DDZ*W ;

