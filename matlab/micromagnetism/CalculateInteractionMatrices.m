function InteractionMatrices = CalculateInteractionMatrices(ProblemSetupStruct)

    %--- Evaluate all variables in the ProblemSetupStruct
    names = fieldnames(ProblemSetupStruct);
    for i=1:length(names)
        eval([names{i} '=ProblemSetupStruct.' names{i} ';']);
    end

    %% All the one-time calculations:
    Nprod = cumprod(nGrid,2) ;
    K = size(nGrid,2) ;

    % N = Nprod(:,end) ;    %--- Remark: Had to change this code
    N = Nprod ;  

    %% Space Grid
    % NN = prod(N) ; % total number of tiles
    [x,y,z,dx,dy,dz] = MakeTheGrid(N,Lx,Ly,Lz) ;
    [InteractionMatrices.X,InteractionMatrices.Y,InteractionMatrices.Z] = ndgrid(x,y,z) ; % always ndgrid, never meshgrid. Meshgrid varies y before x.

    %% Demag Tensors (approximate method: "Many Ranges")
    [InteractionMatrices.AvrgMatrix,InteractionMatrices.CopyMatrix,InteractionMatrices.DemagTensor] = CreateManyRangesDemagTensors3D(N,K,nGrid,Nprod,dx,dy,dz) ;
    disp('Demag Done') ;

    %% Anisotropy Matrix

    if Nmax >1 % Remark: more than one "grain" (Voronoi) -- not yet implemented
        [Kxx,Kxy,Kxz,Kyx,Kyy,Kyz,Kzx,Kzy,Kzz,xPoint,yPoint,Hx,Hy,Hz] = CreateAnisotropyMatrixGrains(X,Y,Nmax,.05,InPlane) ;
    else % one "grain"
        [InteractionMatrices.Kxx,InteractionMatrices.Kxy,InteractionMatrices.Kxz,InteractionMatrices.Kyx, ...
         InteractionMatrices.Kyy,InteractionMatrices.Kyz,InteractionMatrices.Kzx,InteractionMatrices.Kzy, ...
         InteractionMatrices.Kzz] = CreateAnisotropyMatrixBase(Kx,Ky,Kz) ;
    end

    %% Exchange Interaction Matrix

    % TheTopology = 'Flat' ; % Torus yet to be implemented
    % if 0 %  exist(['Ising2DAdjMatrix',TheTopology,num2str(N),'.mat'])
    %     load(['Ising2DAdjMatrix',TheTopology,num2str(N),'.mat'],'A2') ;
    % else
    %     ComputeExchangeTerm
    InteractionMatrices.A2 = ComputeExchangeTerm3D(N,dx,dy,dz) ;
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
if ischar(ProblemSetupStruct.DemagTensorFileName) 
    save(DemagTensorFileName,'InteractionMatrices') ;
end
end