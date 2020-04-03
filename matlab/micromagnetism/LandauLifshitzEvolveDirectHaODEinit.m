function [AllSigmas] = LandauLifshitzEvolveDirectHaODEinit(ProblemSetupStruct,InteractionMatrices,Sigma)

%--- Evaluate all variables in the ProblemSetupStruct
names = fieldnames(ProblemSetupStruct);
for i=1:length(names)
    eval([names{i} '=ProblemSetupStruct.' names{i} ';']);
end
    
ProblemSetupStructSingleStep = ProblemSetupStruct;
% ProblemSetupStructSingleStep.t = linspace(ProblemSetupStruct.t(1),10*ProblemSetupStruct.t(1),2);
% ProblemSetupStructSingleStep.nt = 2;
ProblemSetupStructSingleStep.Hext = ProblemSetupStruct.Hext(:,:).*0 + ProblemSetupStruct.Hext(1,:);
ProblemSetupStructSingleStep.Hext(:,1) = linspace(ProblemSetupStructSingleStep.t(1), ProblemSetupStructSingleStep.t(2), numel(ProblemSetupStructSingleStep.Hext(:,1)));
%time-dependent applied field

% originalHsX = HsX ;
% originalHsY = HsY ;
% originalHsZ = HsZ ;
% originalnT = nT;
% 
% ProblemSetupStruct.HsX = @(t) originalHsX(minT) + 0.*t;
% ProblemSetupStruct.HsY = @(t) originalHsY(minT) + 0.*t;
% ProblemSetupStruct.HsZ = @(t) originalHsZ(minT) + 0.*t;
% ProblemSetupStruct.nT = 2;

%% Initial equilibrium configuration found by explicit solver
% ProblemSetupStructSingleStep = ProblemSetupStructSingleStep.setSolverType( 'UseExplicitSolver' );
ProblemSetupStructSingleStep.SolverType = 3;
% if ~( ProblemSetupStruct.AlreadyEquilibrium & ProblemSetupStruct.InitialState)

    [Sigma] = LandauLifshitzEvolveCombined(ProblemSetupStructSingleStep,InteractionMatrices,Sigma);
    
    StartingSigma = Sigma(end,:)';  % We only want the converged value
% else
%    StartingSigma =  ProblemSetupStruct.SigmaIN ;
% end

ProblemSetupStructImplicit = ProblemSetupStruct;
ProblemSetupStructImplicit.SolverType = 1;
ProblemSetupStructImplicit.t = ProblemSetupStruct.Hext(:,1)';

% ProblemSetupStruct.HsX = originalHsX;
% ProblemSetupStruct.HsY = originalHsY;
% ProblemSetupStruct.HsZ = originalHsZ;
% ProblemSetupStruct.nT = originalnT;
% ProblemSetupStruct.t = t;

%% Remaining steps taken by implicit solver
% ProblemSetupStruct.UseExplicitSolver = 0;
% ProblemSetupStruct.UseImplicitSolver = 1;
[AllSigmas] = LandauLifshitzEvolveCombined(ProblemSetupStructImplicit,InteractionMatrices,StartingSigma);
end