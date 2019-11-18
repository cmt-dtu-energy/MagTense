function Sigma = InitialSigma(ProblemSetupStruct,InteractionMatrices)

    %--- Evaluate all variables in the ProblemSetupStruct
    names = fieldnames(ProblemSetupStruct);
    for i=1:length(names)
        eval([names{i} '=ProblemSetupStruct.' names{i} ';']);
    end
    
    %--- Evaluate some variables in the ProblemSetupStruct.
    N = InteractionMatrices.N;
    X = InteractionMatrices.X;
    Y = InteractionMatrices.Y;
    Z = InteractionMatrices.Z;

    
%% Initial Sigma
if ischar(ProblemSetupStruct.SigmaInitialFileName) 
    if exist(ProblemSetupStruct.SigmaInitialFileName,'file')
    %--- Remark: We should setting on the structure of how it is saved...
    load(ProblemSetupStruct.SigmaInitialFileName,'SigmaInitial') ;
    NNold = numel(SigmaInitial)/3;
    SigmaX = SigmaInitial(0*NNold+[1:NNold]) ;
    SigmaY = SigmaInitial(1*NNold+[1:NNold]) ;
    SigmaZ = SigmaInitial(2*NNold+[1:NNold]) ;
    end
else
    switch InitialState
        case 'Input' % Function specified by user
            SigmaX = SigmaXfun(X,Y,Z) ;
            SigmaY = SigmaYfun(X,Y,Z) ;
            SigmaZ = SigmaZfun(X,Y,Z) ;
            
        case 'rand' % random state
            SigmaPh = (2*pi).*(rand(N(1),N(2),N(3))-0.5) ;
            ThetaPh = (pi  ).*(rand(N(1),N(2),N(3))    ) ;
            SigmaX = sin(ThetaPh).*cos(SigmaPh) ;
            SigmaY = sin(ThetaPh).*sin(SigmaPh) ;
            SigmaZ = cos(ThetaPh) ;
            
        case 'OldSol' % Configuration specified by user
            NNold = numel(SigmaIN)/3 ;
            SigmaXIN = SigmaIN(0*NNold+[1:NNold]) ;
            SigmaYIN = SigmaIN(1*NNold+[1:NNold]) ;
            SigmaZIN = SigmaIN(2*NNold+[1:NNold]) ;
    
            if ~isequal(N,Nold) % Interpolation must be done. Doesn't work if N or Nold has only 1 layer in x- or y-direction.
                if min(Nold) == 1 && min(N) ~= 1 % The hacky and bad way           
                    disp('Warning!! New grid is 3D and old grid is (maximally) 2D. Assuming Nold(3) == 1')
                    [xold,yold,zold,~,~,~] = MakeTheGrid([Nold(1);Nold(2);2],Lx,Ly,Lz) ;
                    [Xold,Yold,Zold] = ndgrid(xold,yold,zold) ;
                    SigmaX = interpn(Xold,Yold,Zold,reshape([SigmaXIN;SigmaXIN],Nold(1),Nold(2),2),X,Y,Z) ;
                    SigmaY = interpn(Xold,Yold,Zold,reshape([SigmaYIN;SigmaYIN],Nold(1),Nold(2),2),X,Y,Z) ;
                    SigmaZ = interpn(Xold,Yold,Zold,reshape([SigmaZIN;SigmaZIN],Nold(1),Nold(2),2),X,Y,Z) ;
                else % The almost good way
                    [xold,yold,zold,~,~,~] = MakeTheGrid(Nold,Lx,Ly,Lz) ;
                    [Xold,Yold,Zold] = ndgrid(xold,yold,zold) ;
                    SigmaX = interpn(Xold,Yold,Zold,reshape(SigmaXIN,Nold(1),Nold(2),Nold(3)),X,Y,Z) ;
                    SigmaY = interpn(Xold,Yold,Zold,reshape(SigmaYIN,Nold(1),Nold(2),Nold(3)),X,Y,Z) ;
                    SigmaZ = interpn(Xold,Yold,Zold,reshape(SigmaZIN,Nold(1),Nold(2),Nold(3)),X,Y,Z) ;
                end
            else
                SigmaX = SigmaXIN ;
                SigmaY = SigmaYIN ;
                SigmaZ = SigmaZIN ;
            end
    end
end
Sigma = NormalizeSigmaCAT([SigmaX(:);SigmaY(:);SigmaZ(:)]) ;
% 
% SigmaX = SigmaX(:) ; 
% SigmaY = SigmaY(:) ;
% SigmaZ = SigmaZ(:) ;
% 
% [SigmaX,SigmaY,SigmaZ] = NormalizeSigma(SigmaX,SigmaY,SigmaZ) ;
% 
% Sigma(1,1,:) = SigmaX;
% Sigma(1,2,:) = SigmaY;
% Sigma(1,3,:) = SigmaZ;
end