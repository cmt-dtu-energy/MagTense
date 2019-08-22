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
% originalt = t;
% % originalMaxT = MaxT;
originalnT = nT;

HsX = @(t) originalHsX(minT) + 0.*t ;
HsY = @(t) originalHsY(minT) + 0.*t ;
HsZ = @(t) originalHsZ(minT) + 0.*t ;

ProblemSetupStruct.HsX = HsX;
ProblemSetupStruct.HsY = HsY;
ProblemSetupStruct.HsZ = HsZ;
% ProblemSetupStruct.t = t(1);
% % ProblemSetupStruct.MaxT = 0;
ProblemSetupStruct.nT = 2;

%% Initial equilibrium configuration found by explicit solver
ProblemSetupStruct.UseExplicitSolver = 1;
ProblemSetupStruct.UseImplicitSolver = 0;
if ~( ProblemSetupStruct.AlreadyEquilibrium & ProblemSetupStruct.InitialState)
[Sigma] = LandauLifshitzEvolveCombined(ProblemSetupStruct,InteractionMatrices,Sigma) ;
close(gcf)
%%
% We only want the converged value
StartingSigma = Sigma(end,:)';
else
   StartingSigma =  ProblemSetupStruct.SigmaIN ;
end
HsX = originalHsX ;
HsY = originalHsY ;
HsZ = originalHsZ ;
% t = originalt;
% MaxT = originalMaxT;
nT = originalnT;

ProblemSetupStruct.HsX = HsX;
ProblemSetupStruct.HsY = HsY;
ProblemSetupStruct.HsZ = HsZ;
ProblemSetupStruct.t = t;
% ProblemSetupStruct.MaxT = 0;
ProblemSetupStruct.nT = nT;

% ThisSolutionFileName = ['InputSolutionStepPreSim.mat'] ;
% save(ThisSolutionFileName,'Sigma') ;

% ThisPreSetFileName2 = [ThisPreSetFileName(1:end-4),'B.mat'] ;
%     save(ThisPreSetFileName2) ;

%% Remaining steps taken by implicit solver
ProblemSetupStruct.UseExplicitSolver = 0;
ProblemSetupStruct.UseImplicitSolver = 1;
[AllSigmas,AllAppliedFields] = LandauLifshitzEvolveCombined(ProblemSetupStruct,InteractionMatrices,StartingSigma);

% if ProblemSetupStruct.SaveTheResult
% save([SimulationName(1:end-4),'solution.mat'],'AllSigmas','AllAppliedFields') ;
% end
end