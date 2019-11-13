function [AvrgMatrix,CopyMatrix,DemagTensor] = CreateManyRangesDemagTensors3D(N,K,nTimes,Nprod,dx,dy,dz)
%
% THIS IS WITH Z !!
% K : number of grids (finest to coarsest)
% N : number of tiles for each space dimension (finest grid)
% a
% dy = dx ; % FIX THIS
% dz = dx ; % FIX THIS
for k=1:K % iterate over the refinements
    disp([num2str(k),'/',num2str(K)]) ;
    Nfx = nTimes(1,K+1-k) ;
    Nfy = nTimes(2,K+1-k) ;
    Nfz = nTimes(3,K+1-k) ;
    xf = dx.*[-(Nfx-1)/2:1:+(Nfx-1)/2] ;
    yf = dy.*[-(Nfy-1)/2:1:+(Nfy-1)/2] ;
    zf = dz.*[-(Nfz-1)/2:1:+(Nfz-1)/2] ;
    [Xf,Yf,Zf] = ndgrid(xf,yf,zf) ;
    if k < K % all the grids except the coarsest
        % create demag tensor connecting adjacent "super-tiles"
        for jx = 1:3  
            for jy = 1:3
                for jz = 1:3
                DeltaX = (jx-2).*dx*(Nfx) ;
                DeltaY = (jy-2).*dy*(Nfy) ;
                DeltaZ = (jz-2).*dz*(Nfz) ;                
                N_tensorF{k}{jx,jy,jz} = ComputeDemagTensorDelta3D(Xf,Yf,Zf,dx,dy,dz,DeltaX,DeltaY,DeltaZ) ; % fine
                
                end
            end
        end
        [KglobXX{k},KglobXY{k},KglobXZ{k},KglobYX{k},KglobYY{k},KglobYZ{k},KglobZX{k},KglobZY{k},KglobZZ{k}] = DistributeMatrixBlockN3D(N_tensorF{k},[Nfx;Nfy;Nfz],Nprod(:,K-k)) ; % Only short range demag field
    else % coarsest grid

        N_tensorC = ComputeDemagTensor3D(Xf,Yf,Zf,dx,dy,dz) ; % coarse
        [KglobXX{k},KglobXY{k},KglobXZ{k},KglobYX{k},KglobYY{k},KglobYZ{k},KglobZX{k},KglobZY{k},KglobZZ{k}] = PointWiseToComponents(N_tensorC) ;
    end
    Nc  = prod(nTimes(:,1:k-1),2) ; % we are here !!
%     Nc = [Nfx;Nfy;Nfz] ;
% Nc = Nprod(:,K-k) ;
% Nc = prod(nTimes(:,1:K-k+1),2).^2 ;
Nc = prod(nTimes(:,1:K-k+1),2) ;
    if k>1 % all the grids except the finest
        % erase from the demag tensor the elements connecting tiles
        % belonging to the same "super-tile"
        % or to adjacent "super-tiles"
        disp(num2str(size(KglobXX{k})))
        [KglobXX{k},KglobXY{k},KglobXZ{k},KglobYX{k},KglobYY{k},KglobYZ{k},KglobZX{k},KglobZY{k},KglobZZ{k}] = EraseAdjacentBlocksN3D(Nc,KglobXX{k},KglobXY{k},KglobXZ{k},KglobYX{k},KglobYY{k},KglobYZ{k},KglobZX{k},KglobZY{k},KglobZZ{k}) ;
        disp(num2str(size(KglobXX{k}))) 
        AvrgMatrix{k} = ComputeAvrgMatrix3D(N(:)./Nprod(:,K+1-k),Nprod(:,K+1-k),1) ; % Matrices that do the average ( of Sigma ) 
        CopyMatrix{k} = ComputeAvrgMatrix3D(N(:)./Nprod(:,K+1-k),Nprod(:,K+1-k),0).' ; % Matrices that do the copy ( of Heff )

    end
    dx = Nfx*dx ;
    dy = Nfy*dy ;
    dz = Nfz*dz ;
end

DemagTensor.KglobXX = KglobXX;
DemagTensor.KglobXY = KglobXY;
DemagTensor.KglobXZ = KglobXZ;
DemagTensor.KglobYX = KglobYX;
DemagTensor.KglobYY = KglobYY;
DemagTensor.KglobYZ = KglobYZ;
DemagTensor.KglobZX = KglobZX;
DemagTensor.KglobZY = KglobZY;
DemagTensor.KglobZZ = KglobZZ;

if ~exist('AvrgMatrix') 
    AvrgMatrix = [] ;
    CopyMatrix = [] ;
end