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
    %% FFT
    if ~all(FFTdims==0)
    [InteractionMatrices] = ConvertTensorToFourier(InteractionMatrices,N,FFTdims) ; 
    end
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

    if (ProblemSetupStruct.thresholdFract > 0)
        if ~all(FFTdims==0)
        [ProblemSetupStruct.threshold] = GetThreshold(InteractionMatrices.FFT.DemagTensor,ProblemSetupStruct.thresholdFract ) ;        
        else
        [ProblemSetupStruct.threshold] = GetThreshold(InteractionMatrices.DemagTensor,ProblemSetupStruct.thresholdFract ) ;
        end
    end

    if (ProblemSetupStruct.threshold > 0) 
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
%         InteractionMatrices.DemagTensor.KglobXY{1,1} = sparse(InteractionMatrices.DemagTensor.KglobXY{1,1});
%         InteractionMatrices.DemagTensor.KglobXZ{1,1} = sparse(InteractionMatrices.DemagTensor.KglobXZ{1,1});
%         InteractionMatrices.DemagTensor.KglobYY{1,1} = sparse(InteractionMatrices.DemagTensor.KglobYY{1,1});
%         InteractionMatrices.DemagTensor.KglobYZ{1,1} = sparse(InteractionMatrices.DemagTensor.KglobYZ{1,1});
%         InteractionMatrices.DemagTensor.KglobZZ{1,1} = sparse(InteractionMatrices.DemagTensor.KglobZZ{1,1});
        end
        InteractionMatrices.A2 = sparse(InteractionMatrices.A2);
    end
    
    if (ProblemSetupStruct.use_gpuArray) 
        InteractionMatrices.DemagTensor.KglobXX{1,1} = gpuArray(InteractionMatrices.DemagTensor.KglobXX{1,1});
        InteractionMatrices.DemagTensor.KglobXY{1,1} = gpuArray(InteractionMatrices.DemagTensor.KglobXY{1,1});
        InteractionMatrices.DemagTensor.KglobXZ{1,1} = gpuArray(InteractionMatrices.DemagTensor.KglobXZ{1,1});
        InteractionMatrices.DemagTensor.KglobYY{1,1} = gpuArray(InteractionMatrices.DemagTensor.KglobYY{1,1});
        InteractionMatrices.DemagTensor.KglobYZ{1,1} = gpuArray(InteractionMatrices.DemagTensor.KglobYZ{1,1});
        InteractionMatrices.DemagTensor.KglobZZ{1,1} = gpuArray(InteractionMatrices.DemagTensor.KglobZZ{1,1});
        InteractionMatrices.A2 = gpuArray(InteractionMatrices.A2);
    end

    
    if ischar(ProblemSetupStruct.DemagTensorFileName) 
        save(DemagTensorFileName,'InteractionMatrices') ;
    end
end