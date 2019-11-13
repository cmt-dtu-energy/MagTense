function [AllSigmas, AllAppliedFields] = LandauLifshitzEvolveDirectHaODEinit(ProblemSetupStruct,InteractionMatrices,Sigma)

%--- Evaluate all variables in the ProblemSetupStruct
names = fieldnames(ProblemSetupStruct);
for i=1:length(names)
    eval([names{i} '=ProblemSetupStruct.' names{i} ';']);
end
    
%--- Remark: What is this code?
if exist('t','var')
    minT = t(1) ;
else
    minT = -1 ;
end
originalHsX = HsX ;
originalHsY = HsY ;
originalHsZ = HsZ ;
originalnT = nT;

ProblemSetupStruct.HsX = @(t) originalHsX(minT) + 0.*t;
ProblemSetupStruct.HsY = @(t) originalHsY(minT) + 0.*t;
ProblemSetupStruct.HsZ = @(t) originalHsZ(minT) + 0.*t;
ProblemSetupStruct.nT = 2;

%% Initial equilibrium configuration found by explicit solver
ProblemSetupStruct.UseExplicitSolver = 1;
ProblemSetupStruct.UseImplicitSolver = 0;
if ~( ProblemSetupStruct.AlreadyEquilibrium & ProblemSetupStruct.InitialState)
    [Sigma] = LandauLifshitzEvolveCombined(ProblemSetupStruct,InteractionMatrices,Sigma);
    
    StartingSigma = Sigma(end,:)';  % We only want the converged value
else
   StartingSigma =  ProblemSetupStruct.SigmaIN ;
end

ProblemSetupStruct.HsX = originalHsX;
ProblemSetupStruct.HsY = originalHsY;
ProblemSetupStruct.HsZ = originalHsZ;
ProblemSetupStruct.nT = originalnT;
ProblemSetupStruct.t = t;

%% Remaining steps taken by implicit solver
ProblemSetupStruct.UseExplicitSolver = 0;
ProblemSetupStruct.UseImplicitSolver = 1;
[AllSigmas,AllAppliedFields] = LandauLifshitzEvolveCombined(ProblemSetupStruct,InteractionMatrices,StartingSigma);
end