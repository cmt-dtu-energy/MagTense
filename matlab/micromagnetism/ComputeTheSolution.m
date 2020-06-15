function [SigmaSol,VV,InteractionMatrices] = ComputeTheSolution(MySim)
VV = [];
tic

%% Setup problem variables

[ProblemSetupStruct, PreSetFileName] = SetupProblem(MySim);


%% Calculate or load the matrices that only needs to be calculated once

if ischar(ProblemSetupStruct.DemagTensorFileName) & (~ProblemSetupStruct.RecomputeInteractionMatrices)
    if exist([ProblemSetupStruct.DemagTensorFileName,'.mat'],'file') || exist(ProblemSetupStruct.DemagTensorFileName,'file')
        %--- Remark: We should settle on the structure of how it is saved...
        load(ProblemSetupStruct.DemagTensorFileName,'InteractionMatrices') ;
        MatricesHaveBeenLoaded = 1 ;
    else
        MatricesHaveBeenLoaded = 0 ;
    end
else
   MatricesHaveBeenLoaded = 0 ; 
end

if ~MatricesHaveBeenLoaded % ComputeThem
    InteractionMatrices = CalculateInteractionMatrices(ProblemSetupStruct);
end
N = InteractionMatrices.N;


%% The starting guess for the magnetization
% This is randomized or loaded from a file
Sigma = MySim.m0(:);
% Sigma = InitialSigma(ProblemSetupStruct,InteractionMatrices);

%% Static, "implicit" equilibrium solution

%UseImplicitSolver
if ProblemSetupStruct.SolverType == 1
    
    [SigmaSol]  = LandauLifshitzEvolveDirectHaODEinit(ProblemSetupStruct,InteractionMatrices,Sigma);
    if ProblemSetupStruct.SaveTheResult
        elapsedTime = toc;
        save([PreSetFileName,'solution'],'N','SigmaSol','elapsedTime');
    end
    
%UseImplicitStepsSolver
elseif ProblemSetupStruct.SolverType == 2
    
    [SigmaSol,ImplicitFailed,VV] = LandauLifshitzEvolveImplicitSteps(ProblemSetupStruct,InteractionMatrices,Sigma);
    if ProblemSetupStruct.SaveTheResult
        elapsedTime = toc;
        save([PreSetFileName,'solution'],'N','SigmaSol','ImplicitFailed','VV','elapsedTime');
    end
    
%UseExplicitSolver
elseif ProblemSetupStruct.SolverType == 3
    
    [SigmaSol, VV] = LandauLifshitzEvolveStepByStepInit(ProblemSetupStruct,InteractionMatrices,Sigma);
    if ProblemSetupStruct.SaveTheResult
        elapsedTime = toc;
        save([PreSetFileName,'solution.mat'],'N','SigmaSol','VV','elapsedTime');
    end
    
%UseDynamicSolver
elseif ProblemSetupStruct.SolverType == 4
    
    [SigmaSol] = LandauLifshitzEvolveCombined(ProblemSetupStruct,InteractionMatrices,Sigma) ;
    if ProblemSetupStruct.SaveTheResult
        elapsedTime = toc;
        save([PreSetFileName,'solution'],'N','SigmaSol','elapsedTime');
    end
    
end
end