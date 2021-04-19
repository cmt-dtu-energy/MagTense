function InteractionMatrices = CalculateInteractionMatrices(ProblemSetupStruct)

%--- Evaluate all variables in the ProblemSetupStruct
names = fieldnames(ProblemSetupStruct);
for i=1:length(names)
    eval([names{i} '=ProblemSetupStruct.' names{i} ';']);
end
grid_n = double(grid_n');

%% All the one-time calculations:
if ~ExternalMesh % Voronoi & Tetra
    Nprod = cumprod(grid_n,2) ;
    K = size(grid_n,2) ;
    
    % N = Nprod(:,end) ;    %--- Remark: Had to change this code
    N = Nprod ;
    
    %% Space Grid
    % NN = prod(N) ; % total number of tiles
    [x,y,z,dx,dy,dz] = MakeTheGrid(N,grid_L(1),grid_L(2),grid_L(3)) ;
    [InteractionMatrices.X,InteractionMatrices.Y,InteractionMatrices.Z] = ndgrid(x,y,z) ; % always ndgrid, never meshgrid. Meshgrid varies y before x.
    
    %% Demag Tensors (approximate method: "Many Ranges")
    [InteractionMatrices.AvrgMatrix,InteractionMatrices.CopyMatrix,InteractionMatrices.DemagTensor] = CreateManyRangesDemagTensors3D(N,K,grid_n,Nprod,dx,dy,dz) ;
else % Voronoi & Tetra
    
    
    
    if isequal(MeshType,'Voronoi')
%         load('CompareInteractionMatrices.mat') ; InteractionMatricesCompare =  InteractionMatrices ;
        load(ExternalMeshFileName,'voronoi_map') ; % Voronoi & Tetra
            N = [size(voronoi_map.pos,1),1,1] ;
        InteractionMatrices.X = voronoi_map.pos(:,1) ; % Voronoi & Tetra
        InteractionMatrices.Y = voronoi_map.pos(:,2) ; % Voronoi & Tetra
        InteractionMatrices.Z = voronoi_map.pos(:,3) ; % Voronoi & Tetra
    
        
        dx = voronoi_map.dims(:,1) ; % Voronoi & Tetra
        dy = voronoi_map.dims(:,2) ; % Voronoi & Tetra
        dz = voronoi_map.dims(:,3) ; % Voronoi & Tetra

% InteractionMatrices.X = InteractionMatricesCompare.X(:) ;
% InteractionMatrices.Y = InteractionMatricesCompare.Y(:) ;
% InteractionMatrices.Z = InteractionMatricesCompare.Z(:) ;
% 
%             dx = repmat(InteractionMatricesCompare.X(2,1,1)-InteractionMatricesCompare.X(1,1,1),size(voronoi_map.dims,1)) ;
%             dy = repmat(InteractionMatricesCompare.Y(1,2,1)-InteractionMatricesCompare.Y(1,1,1),size(voronoi_map.dims,1)) ;
%             dz = repmat(InteractionMatricesCompare.Z(1,1,2)-InteractionMatricesCompare.Z(1,1,1),size(voronoi_map.dims,1)) ;
%         N_tensor = ComputeDemagTensorCartesianMesh(InteractionMatricesCompare.X(:),InteractionMatricesCompare.Y(:),InteractionMatricesCompare.Z(:),dx,dy,dz) ; % Voronoi & Tetra        
        N_tensor = ComputeDemagTensorCartesianMesh(InteractionMatrices.X,InteractionMatrices.Y,InteractionMatrices.Z,dx,dy,dz) ; % Voronoi & Tetra
        
        
        GridInfo = CartesianUnstructuredMeshAnalysis(voronoi_map.pos,voronoi_map.dims) ; % Voronoi & Tetra
        %         [DX0,DY0,DZ0,~] = ComputeDifferentialOperatorsFromMesh03_IDW(GridInfo.fNormX,GridInfo.fNormY,GridInfo.fNormZ,GridInfo.AreaFaces,GridInfo.Volumes,GridInfo.TheSigns,GridInfo.Xel,GridInfo.Yel,GridInfo.Zel,GridInfo.Xf,GridInfo.Yf,GridInfo.Zf,GridInfo.TheTs,abs(GridInfo.TheSigns.')) ;
        %         InteractionMatrices.A2 = (DX0*DX0 + DY0*DY0 + DZ0*DZ0) ; % Voronoi & Tetra
    elseif isequal(MeshType,'Tetra')
        load(ExternalMeshFileName,'model') ; % Voronoi & Tetra
        GridInfo = TetrahedralMeshAnalysis(model) ;
        InteractionMatrices.X = GridInfo.Xel ;
        InteractionMatrices.Y = GridInfo.Yel ;
        InteractionMatrices.Z = GridInfo.Zel ;

        InteractionMatrices.model = model ; %  REMOVE THIS ??
        N = numel(GridInfo.Xel) ;
        
        dx = (GridInfo.Volumes).^(1/3) ;
        dy = dx ;
        dz = dx ;
        N_tensor = zeros(3*N,3*N) ;
        
        tile = getDefaultMagTile();
        tile.magnetType = getMagnetType('hard');
        tile.tileType = getMagTileType('tetrahedron');
        for n1=1:N
            TheseIj = model.Mesh.Elements(1:4,n1) ;
            v = model.Mesh.Nodes(:,TheseIj) ;
%             for n2 =1:N
%                 r = [InteractionMatrices.X(n2);InteractionMatrices.Y(n2);InteractionMatrices.Z(n2)] ;
%                 [N_tensor(((n2-1)*3+1):(n2*3),((n1-1)*3+1):(n1*3))] = getNTetrahedron_Matlab( r, v ) ;
%                 '' ;
%             end
            
            %Calculate the demag tensor directly from Fortran
            tile.vertices = v;
            sjask = [InteractionMatrices.X(:),InteractionMatrices.Y(:),InteractionMatrices.Z(:)];
            N_fortran_temp = getNFromTile_mex( tile, sjask, int32( length( sjask(:,1) )) );
            N_tensor(:,((n1-1)*3+1):(n1*3)) = reshape(shiftdim(N_fortran_temp,2),[3*N,3]);
        end
        
    end

    [InteractionMatrices.DemagTensor.KglobXX{1},...
        InteractionMatrices.DemagTensor.KglobXY{1},...
        InteractionMatrices.DemagTensor.KglobXZ{1},...
        ~,...
        InteractionMatrices.DemagTensor.KglobYY{1},...
        InteractionMatrices.DemagTensor.KglobYZ{1},...
        ~,...
        ~,...
        InteractionMatrices.DemagTensor.KglobZZ{1}] = PointWiseToComponents(N_tensor) ; % Voronoi & Tetra
    InteractionMatrices.AvrgMatrix = [] ; % Voronoi & Tetra
    InteractionMatrices.CopyMatrix = [] ; % Voronoi & Tetra

    '' ;
    InteractionMatrices.GridInfo = GridInfo ;
end % Voronoi & Tetra
disp('Demag Done') ;

disp('Demag Done') ;
%% FFT
if ~all(FFTdims==0)
    [InteractionMatrices] = ConvertTensorToFourier(InteractionMatrices,N,FFTdims) ;
end
%% Anisotropy Matrix

%     if Nmax >1 % Remark: more than one "grain" (Voronoi) -- not yet implemented
%         [Kxx,Kxy,Kxz,Kyx,Kyy,Kyz,Kzx,Kzy,Kzz,xPoint,yPoint,Hx,Hy,Hz] = CreateAnisotropyMatrixGrains(X,Y,Nmax,.05,InPlane) ;
%     else % one "grain"
%         [InteractionMatrices.Kxx,InteractionMatrices.Kxy,InteractionMatrices.Kxz,InteractionMatrices.Kyx, ...
%          InteractionMatrices.Kyy,InteractionMatrices.Kyz,InteractionMatrices.Kzx,InteractionMatrices.Kzy, ...
%          InteractionMatrices.Kzz] = CreateAnisotropyMatrixBase(Kx,Ky,Kz) ;
[InteractionMatrices.Kxx,InteractionMatrices.Kxy,InteractionMatrices.Kxz,InteractionMatrices.Kyx, ...
    InteractionMatrices.Kyy,InteractionMatrices.Kyz,InteractionMatrices.Kzx,InteractionMatrices.Kzy, ...
    InteractionMatrices.Kzz] = CreateAnisotropyMatrixBase(u_ea(:,1),u_ea(:,2),u_ea(:,3)) ;

%     end

%% Exchange Interaction Matrix

% TheTopology = 'Flat' ; % Torus yet to be implemented
% if 0 %  exist(['Ising2DAdjMatrix',TheTopology,num2str(N),'.mat'])
%     load(['Ising2DAdjMatrix',TheTopology,num2str(N),'.mat'],'A2') ;
% else
%     ComputeExchangeTerm
if ~ExternalMesh % Voronoi & Tetra
    InteractionMatrices.A2 = ComputeExchangeTerm3D(N,dx,dy,dz) ;
else % Voronoi & Tetra
    % thS = 1e15 ;
    %     InteractionMatrices.A2 = ComputeDifferentialOperatorsFromMesh03_IDW_2nd(GridInfo.fNormX,GridInfo.fNormY,GridInfo.fNormZ,(thS^2).*GridInfo.AreaFaces,(thS^3).*GridInfo.Volumes,GridInfo.TheSigns,(thS).*GridInfo.Xel,(thS).*GridInfo.Yel,(thS).*GridInfo.Zel,(thS).*GridInfo.Xf,(thS).*GridInfo.Yf,(thS).*GridInfo.Zf,GridInfo.TheTs,abs(GridInfo.TheSigns.')) ;
    
    if  isequal(MeshType,'Voronoi')
        InteractionMatrices.A2 = ComputeDifferentialOperatorsFromMesh03_IDW_2nd(GridInfo.fNormX,GridInfo.fNormY,GridInfo.fNormZ,GridInfo.AreaFaces,GridInfo.Volumes,GridInfo.TheSigns,GridInfo.Xel,GridInfo.Yel,GridInfo.Zel,GridInfo.Xf,GridInfo.Yf,GridInfo.Zf,GridInfo.TheTs,abs(GridInfo.TheSigns.')) ;
    else
        [DX,DY,DZ] = ComputeDifferentialOperatorsFromMesh03_IDW(GridInfo.fNormX,GridInfo.fNormY,GridInfo.fNormZ,GridInfo.AreaFaces,GridInfo.Volumes,GridInfo.TheSigns,GridInfo.Xel,GridInfo.Yel,GridInfo.Zel,GridInfo.Xf,GridInfo.Yf,GridInfo.Zf,GridInfo.TheTs,abs(GridInfo.TheSigns.')) ;
        InteractionMatrices.A2 = DX*DX + DY*DY+ DZ*DZ ;
        % InteractionMatrices.A2 = ComputeDifferentialOperatorsFromMesh03_IDW_2nd(GridInfo.fNormX,GridInfo.fNormY,GridInfo.fNormZ,GridInfo.AreaFaces,GridInfo.Volumes,GridInfo.TheSigns,GridInfo.Xel,GridInfo.Yel,GridInfo.Zel,GridInfo.Xf,GridInfo.Yf,GridInfo.Zf,GridInfo.TheTs,GridInfo.TheTs) ;
    end
    InteractionMatrices.A2 = (InteractionMatrices.A2)-diag(diag((InteractionMatrices.A2))) ;
    
end % Voronoi & Tetra
%     save(['Ising2DAdjMatrix',TheTopology,num2str(N),'.mat'],'A2') ;
% end

InteractionMatrices.N = N;
%     InteractionMatrices.X = X;
%     InteractionMatrices.Y = Y;
%     InteractionMatrices.Z = Z;
InteractionMatrices.dx = dx;
InteractionMatrices.dy = dy;
InteractionMatrices.dz = dz;
%     InteractionMatrices.A2 = A2;
%     InteractionMatrices.Kxx = Kxx;
%     InteractionMatrices.Kxy = Kxy;
%     InteractionMatrices.Kxz = Kxz;
%     InteractionMatrices.Kyx = Kyx;
%     InteractionMatrices.Kyy = Kyy;
%     InteractionMatrices.Kyz = Kyz;
%     InteractionMatrices.Kzx = Kzx;
%     InteractionMatrices.Kzy = Kzy;
%     InteractionMatrices.Kzz = Kzz;
%     InteractionMatrices.AvrgMatrix = AvrgMatrix;
%     InteractionMatrices.CopyMatrix = CopyMatrix;

if (ProblemSetupStruct.dem_appr == 4 )  %DemagApproximationThresholdFraction
    if ~all(FFTdims==0)
        [ProblemSetupStruct.threshold] = GetThreshold(InteractionMatrices.FFT.DemagTensor,ProblemSetupStruct.thresholdFract ) ;
    else
        [ProblemSetupStruct.threshold] = GetThreshold(InteractionMatrices.DemagTensor,ProblemSetupStruct.thresholdFract ) ;
    end
end

if (ProblemSetupStruct.dem_appr == 2) %DemagApproximationThreshold
    if ~all(FFTdims==0)
        [InteractionMatrices.FFT.DemagTensor] = FilterTensorThreshold(InteractionMatrices.FFT.DemagTensor,ProblemSetupStruct.threshold) ;
    else
        [InteractionMatrices.DemagTensor] = FilterTensorThreshold(InteractionMatrices.DemagTensor,ProblemSetupStruct.threshold) ;
        %             InteractionMatrices.DemagTensor.KglobXX{1,1}(abs(InteractionMatrices.DemagTensor.KglobXX{1,1}) < ProblemSetupStruct.threshold) = 0;
        %         InteractionMatrices.DemagTensor.KglobXY{1,1}(abs(InteractionMatrices.DemagTensor.KglobXY{1,1}) < ProblemSetupStruct.threshold) = 0;
        %         InteractionMatrices.DemagTensor.KglobXZ{1,1}(abs(InteractionMatrices.DemagTensor.KglobXZ{1,1}) < ProblemSetupStruct.threshold) = 0;
        %         InteractionMatrices.DemagTensor.KglobYY{1,1}(abs(InteractionMatrices.DemagTensor.KglobYY{1,1}) < ProblemSetupStruct.threshold) = 0;
        %         InteractionMatrices.DemagTensor.KglobYZ{1,1}(abs(InteractionMatrices.DemagTensor.KglobYZ{1,1}) < ProblemSetupStruct.threshold) = 0;
        %         InteractionMatrices.DemagTensor.KglobZZ{1,1}(abs(InteractionMatrices.DemagTensor.KglobZZ{1,1}) < ProblemSetupStruct.threshold) = 0;
    end
    %--- Test if any matrices are zero, and if they are remove them
    if ~any(any(InteractionMatrices.DemagTensor.KglobXX{1,1}))
        InteractionMatrices.DemagTensor.KglobXX{1,1} = 0;
    end
    if ~any(any(InteractionMatrices.DemagTensor.KglobXY{1,1}))
        InteractionMatrices.DemagTensor.KglobXY{1,1} = 0;
    end
    if ~any(any(InteractionMatrices.DemagTensor.KglobXZ{1,1}))
        InteractionMatrices.DemagTensor.KglobXZ{1,1} = 0;
    end
    if ~any(any(InteractionMatrices.DemagTensor.KglobYY{1,1}))
        InteractionMatrices.DemagTensor.KglobYY{1,1} = 0;
    end
    if ~any(any(InteractionMatrices.DemagTensor.KglobYZ{1,1}))
        InteractionMatrices.DemagTensor.KglobYZ{1,1} = 0;
    end
    if ~any(any(InteractionMatrices.DemagTensor.KglobZZ{1,1}))
        InteractionMatrices.DemagTensor.KglobZZ{1,1} = 0;
    end
end

if (ProblemSetupStruct.use_single)
    InteractionMatrices.DemagTensor.KglobXX{1,1} = single(InteractionMatrices.DemagTensor.KglobXX{1,1});
    InteractionMatrices.DemagTensor.KglobXY{1,1} = single(InteractionMatrices.DemagTensor.KglobXY{1,1});
    InteractionMatrices.DemagTensor.KglobXZ{1,1} = single(InteractionMatrices.DemagTensor.KglobXZ{1,1});
    InteractionMatrices.DemagTensor.KglobYY{1,1} = single(InteractionMatrices.DemagTensor.KglobYY{1,1});
    InteractionMatrices.DemagTensor.KglobYZ{1,1} = single(InteractionMatrices.DemagTensor.KglobYZ{1,1});
    InteractionMatrices.DemagTensor.KglobZZ{1,1} = single(InteractionMatrices.DemagTensor.KglobZZ{1,1});
    InteractionMatrices.A2 = single(full(InteractionMatrices.A2));
end

if (ProblemSetupStruct.use_sparse)
    if (ProblemSetupStruct.use_sparse && ProblemSetupStruct.use_single)
        error('A matrix in Matlab cannot be both sparse and single')
    end
    if ~all(FFTdims==0)
        [InteractionMatrices.FFT.DemagTensor] = ConvertTensorToSparse(InteractionMatrices.FFT.DemagTensor) ;
    else
        [InteractionMatrices.DemagTensor] = ConvertTensorToSparse(InteractionMatrices.DemagTensor) ;
        
        %         InteractionMatrices.DemagTensor.KglobXX{1,1} = sparse(InteractionMatrices.DemagTensor.KglobXX{1,1});
        %         InteractionMatrices.DemagTensor.KglobXY{0,1} = sparse(InteractionMatrices.DemagTensor.KglobXY{1,1});
        %         InteractionMatrices.DemagTensor.KglobXZ{1,1} = sparse(InteractionMatrices.DemagTensor.KglobXZ{1,1});
        %         InteractionMatrices.DemagTensor.KglobYY{1,1} = sparse(InteractionMatrices.DemagTensor.KglobYY{1,1});
        %         InteractionMatrices.DemagTensor.KglobYZ{1,1} = sparse(InteractionMatrices.DemagTensor.KglobYZ{1,1});
        %         InteractionMatrices.DemagTensor.KglobZZ{1,1} = sparse(InteractionMatrices.DemagTensor.KglobZZ{1,1});
    end
    InteractionMatrices.A2 = sparse(InteractionMatrices.A2);
end

if (ProblemSetupStruct.useCuda)
    InteractionMatrices.DemagTensor.KglobXX{1,1} = gpuArray(InteractionMatrices.DemagTensor.KglobXX{1,1});
    InteractionMatrices.DemagTensor.KglobXY{1,1} = gpuArray(InteractionMatrices.DemagTensor.KglobXY{1,1});
    InteractionMatrices.DemagTensor.KglobXZ{1,1} = gpuArray(InteractionMatrices.DemagTensor.KglobXZ{1,1});
    InteractionMatrices.DemagTensor.KglobYY{1,1} = gpuArray(InteractionMatrices.DemagTensor.KglobYY{1,1});
    InteractionMatrices.DemagTensor.KglobYZ{1,1} = gpuArray(InteractionMatrices.DemagTensor.KglobYZ{1,1});
    InteractionMatrices.DemagTensor.KglobZZ{1,1} = gpuArray(InteractionMatrices.DemagTensor.KglobZZ{1,1});
    InteractionMatrices.A2 = gpuArray(InteractionMatrices.A2);
end
% InteractionMatrices.DemagTensor = InteractionMatricesCompare.DemagTensor ;

% InteractionMatrices.DemagTensor.KglobXX = InteractionMatricesCompare.DemagTensor.KglobXX ;
% InteractionMatrices.DemagTensor.KglobYY = InteractionMatricesCompare.DemagTensor.KglobYY ;
% InteractionMatrices.DemagTensor.KglobZZ = InteractionMatricesCompare.DemagTensor.KglobZZ ;
% InteractionMatrices.DemagTensor.KglobXY = InteractionMatricesCompare.DemagTensor.KglobXY ;
% InteractionMatrices.DemagTensor.KglobXZ = InteractionMatricesCompare.DemagTensor.KglobXZ ;
% InteractionMatrices.DemagTensor.KglobYZ = InteractionMatricesCompare.DemagTensor.KglobYZ ;
% InteractionMatrices.A2 = InteractionMatricesCompare.A2 ;
%     DoCompare01  


if ischar(ProblemSetupStruct.DemagTensorFileName)
    save(DemagTensorFileName,'InteractionMatrices') ;
end
end