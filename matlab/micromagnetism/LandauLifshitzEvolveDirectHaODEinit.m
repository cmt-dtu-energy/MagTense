function [AllSigmas] = LandauLifshitzEvolveDirectHaODEinit(ProblemSetupStruct,InteractionMatrices,Sigma)

%--- Evaluate all variables in the ProblemSetupStruct
names = fieldnames(ProblemSetupStruct);
for i=1:length(names)
    eval([names{i} '=ProblemSetupStruct.' names{i} ';']);
end
    
ProblemSetupStructSingleStep = ProblemSetupStruct;
ProblemSetupStructSingleStep.Hext = ProblemSetupStruct.Hext(:,:).*0 + ProblemSetupStruct.Hext(1,:);
ProblemSetupStructSingleStep.Hext(:,1) = linspace(ProblemSetupStructSingleStep.t(1), ProblemSetupStructSingleStep.t(2), numel(ProblemSetupStructSingleStep.Hext(:,1)));

%% Initial equilibrium configuration found by explicit solver

ProblemSetupStructSingleStep.SolverType = 3;    %'UseExplicitSolver'

[Sigma] = LandauLifshitzEvolveCombined(ProblemSetupStructSingleStep,InteractionMatrices,Sigma);
    
StartingSigma = Sigma(end,:)';  % We only want the converged value

ProblemSetupStructImplicit = ProblemSetupStruct;
ProblemSetupStructImplicit.SolverType = 1;
ProblemSetupStructImplicit.t = ProblemSetupStruct.Hext(:,1)';

%% Remaining steps taken by implicit solver
[AllSigmas] = LandauLifshitzEvolveCombined(ProblemSetupStructImplicit,InteractionMatrices,StartingSigma);
end