function [SigmaSol,AppliedField,VV,prehapsN, prehapsmatrics, other stuff] = ComputeTheSolution(MySim)

%% Setup problem variables

[ProblemSetupStruct, PreSetFileName] = SetupProblem(MySim);


%% Calculate or load the matrices that only needs to be calculated once

if ischar(ProblemSetupStruct.DemagTensorFileName)
    if exist([ProblemSetupStruct.DemagTensorFileName,'.mat'],'file') || exist(ProblemSetupStruct.DemagTensorFileName,'file')
        %--- Remark: We should setting on the structure of how it is saved...
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
%     N = InteractionMatrices.N;
    if ProblemSetupStruct.SaveTheResult
        save([PreSetFileName,'solution'],'N','SigmaSol','AppliedField');
    end
elseif ProblemSetupStruct.UseExplicitSolver == 1
    

    if ProblemSetupStruct.CalcEigenvalue
        [SigmaSol, AppliedField,VV] = LandauLifshitzEvolveStepByStepInit(ProblemSetupStruct,InteractionMatrices,Sigma);
        if ProblemSetupStruct.SaveTheResult
            save([PreSetFileName,'solution.mat'],'N','SigmaSol','AppliedField','VV');
        end
    else
        [SigmaSol, AppliedField] = LandauLifshitzEvolveStepByStepInit(ProblemSetupStruct,InteractionMatrices,Sigma);
        if ProblemSetupStruct.SaveTheResult
            save([PreSetFileName,'solution.mat'],'N','SigmaSol','AppliedField');
        end
    end
    %--- other remark: Should we save the below calculation?
    %     NetMagnetization = CalculateMagnetization(SigmaSol.') ;
    %     save([SimulationName(1:end-4),'HystLoop.mat'],'tt','AppliedField','NetMagnetization') ;
    
elseif ProblemSetupStruct.UseDynamicSolver == 1
    
    [SigmaSol] = LandauLifshitzEvolveCombined(ProblemSetupStruct,InteractionMatrices,Sigma) ;
    if ProblemSetupStruct.SaveTheResult
        save([PreSetFileName,'solution'],'N','SigmaSol');
    end
elseif ProblemSetupStruct.UseImplicitStepsSolver == 1
    if ProblemSetupStruct.CalcEigenvalue
        [SigmaSol, AppliedField,ImplicitFailed,VV] = LandauLifshitzEvolveImplicitSteps(ProblemSetupStruct,InteractionMatrices,Sigma);
        if ProblemSetupStruct.SaveTheResult
            save([PreSetFileName,'solution'],'N','SigmaSol','ImplicitFailed','VV');
        end
    else
        [SigmaSol, AppliedField,ImplicitFailed] = LandauLifshitzEvolveImplicitSteps(ProblemSetupStruct,InteractionMatrices,Sigma);
        if ProblemSetupStruct.SaveTheResult
            save([PreSetFileName,'solution'],'N','SigmaSol','ImplicitFailed');
        end
    end
end
end