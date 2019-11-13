function [SigmaSol,AppliedField,VV,InteractionMatrices] = ComputeTheSolution(MySim)

%% Setup problem variables

[ProblemSetupStruct, PreSetFileName] = SetupProblem(MySim);


%% Calculate or load the matrices that only needs to be calculated once

if ischar(ProblemSetupStruct.DemagTensorFileName)
    if exist([ProblemSetupStruct.DemagTensorFileName,'.mat'],'file') || exist(ProblemSetupStruct.DemagTensorFileName,'file')
        %--- Remark: We should settle on the structure of how it is saved...
        load(ProblemSetupStruct.DemagTensorFileName,'InteractionMatrices') ;
    else
        InteractionMatrices = CalculateInteractionMatrices(ProblemSetupStruct);
    end
else
    InteractionMatrices = CalculateInteractionMatrices(ProblemSetupStruct);
end
N = InteractionMatrices.N;

%% The starting guess for the magnetization
% This is randomized or loaded from a file

Sigma = InitialSigma(ProblemSetupStruct,InteractionMatrices);


%% Static, "implicit" equilibrium solution

if ProblemSetupStruct.UseImplicitSolver == 1
    [SigmaSol, AppliedField]  = LandauLifshitzEvolveDirectHaODEinit(ProblemSetupStruct,InteractionMatrices,Sigma);
    if ProblemSetupStruct.SaveTheResult
        save([PreSetFileName,'solution'],'N','SigmaSol','AppliedField');
    end
    
elseif ProblemSetupStruct.UseImplicitStepsSolver == 1
    [SigmaSol, AppliedField,ImplicitFailed,VV] = LandauLifshitzEvolveImplicitSteps(ProblemSetupStruct,InteractionMatrices,Sigma);
    if ProblemSetupStruct.SaveTheResult
        save([PreSetFileName,'solution'],'N','SigmaSol','ImplicitFailed','VV');
    end
    
elseif ProblemSetupStruct.UseExplicitSolver == 1
    
    [SigmaSol, AppliedField, VV] = LandauLifshitzEvolveStepByStepInit(ProblemSetupStruct,InteractionMatrices,Sigma);
    if ProblemSetupStruct.SaveTheResult
        save([PreSetFileName,'solution.mat'],'N','SigmaSol','AppliedField','VV');
    end
    
elseif ProblemSetupStruct.UseDynamicSolver == 1
    
    [SigmaSol] = LandauLifshitzEvolveCombined(ProblemSetupStruct,InteractionMatrices,Sigma) ;
    if ProblemSetupStruct.SaveTheResult
        save([PreSetFileName,'solution'],'N','SigmaSol');
    end
    
end
end