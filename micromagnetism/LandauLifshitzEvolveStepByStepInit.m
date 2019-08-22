function [SigmaSol, AppliedField, AllVV] = LandauLifshitzEvolveStepByStepInit(ProblemSetupStruct,InteractionMatrices,Sigma)

%--- Evaluate all variables in the ProblemSetupStruct
names = fieldnames(ProblemSetupStruct);
for i=1:length(names)
    eval([names{i} '=ProblemSetupStruct.' names{i} ';']);
end

% load(ThisPreSetFileName) ;

% if ~exist('MaxT','var')
%     MaxT = 3 ;
% end
% if ~exist('MaxTstep','var')
%     MaxTstep = 5 ;
% end
% originalMaxT = MaxT ;
% if ~exist('tt','var')
%     tt = linspace(0,originalMaxT,4) ;
% end

originalHsX = HsX ;
originalHsY = HsY ;
originalHsZ = HsZ ;

% DemagTensorFileName = [ThisPreSetFileName(1:end-4),'DemagTens.mat'] ;
for k=1:numel(tt)
    disp(['  ',num2str(k),'/',num2str(numel(tt))]) ;
    ProblemSetupStruct.HsX = @(t) 0.*t + originalHsX(tt(k)) ;
    ProblemSetupStruct.HsY = @(t) 0.*t + originalHsY(tt(k)) ;
    ProblemSetupStruct.HsZ = @(t) 0.*t + originalHsZ(tt(k)) ;
    %     save('ThisStuff.mat')
    %     ThisPresetFileNameK = [ThisPreSetFileName(1:end-4),'Step',num2str(k),'.mat'] ;
    %     HeffLimMagn = -21 ; % was - 20 ;
    %     save(ThisPresetFileNameK) ;
    if CalcEigenvalue
        if k==1
            %         LandauLifshitzEvolve3D(ThisPresetFileNameK) ;
            [SigmaSol,~,VV] = LandauLifshitzEvolveCombined(ProblemSetupStruct,InteractionMatrices,Sigma) ;
        else
            %         LandauLifshitzEvolve3D(ThisPresetFileNameK,ThisSolutionFileName) ;
            [SigmaSol,~,VV] = LandauLifshitzEvolveCombined(ProblemSetupStruct,InteractionMatrices,AllSigmas(k-1,:)) ;
        end
    else
        if k==1
            %         LandauLifshitzEvolve3D(ThisPresetFileNameK) ;
            [SigmaSol] = LandauLifshitzEvolveCombined(ProblemSetupStruct,InteractionMatrices,Sigma) ;
        else
            %         LandauLifshitzEvolve3D(ThisPresetFileNameK,ThisSolutionFileName) ;
            [SigmaSol] = LandauLifshitzEvolveCombined(ProblemSetupStruct,InteractionMatrices,AllSigmas(k-1,:)) ;
        end
        
    end
    close all ;
    %     load([ThisPresetFileNameK(1:end-4),'solution']) ;
    AllSigmas(k,:) = SigmaSol(end,:) ;
    if CalcEigenvalue
        AllVV(k,:) = VV(:) ;
    end
    %     ThisSolutionFileName = ['InputSolutionStep',num2str(k),'.mat'] ;
    %     Nold = N ;
    %     save(ThisSolutionFileName,'SigmaIN','Nold') ;
    %     if k==3
    %        '' ;
    %     end
end

%% Erase Files and save Solution (Deprecated. What to to about CalcEigs?)
% for k=1:numel(tt)
%     ThisSolutionFileName = ['InputSolutionStep',num2str(k),'.mat'] ;
%     ThisPresetFileNameK = [ThisPreSetFileName(1:end-4),'Step',num2str(k),'.mat'] ;
%     load([ThisPresetFileNameK(1:end-4),'solution']) ;
%     AllSigmas(k,:) = SigmaSol(end,:) ;
%     delete([ThisPresetFileNameK(1:end-4),'solution.mat']) ;
%     load([ThisPresetFileNameK(1:end-4),'HystLoop.mat']) ;
%     AllApplied_fields(k,:) = Applied_field(end,:) ;
%     NetMagnetization(k,:) = CalculateMagnetization(AllSigmas(k,:).') ;
%     delete([ThisPresetFileNameK(1:end-4),'HystLoop.mat']) ;
%     delete(ThisSolutionFileName) ;
%     delete(ThisPresetFileNameK) ;
%     if exist('CalcEigs','var')
%         if CalcEigs
%             load([ThisPresetFileNameK(1:end-4),'Eigs.mat']) ;
%             AllVV(:,k) = VV ;
%             delete([ThisPresetFileNameK(1:end-4),'Eigs.mat']) ;
%         end
%     end
% end
Hfact = 1/MU0 ; % AppliedField is provided in T
AppliedField = [Hfact*originalHsX(tt), Hfact*originalHsY(tt),Hfact*originalHsZ(tt)+0.*tt] ;
SigmaSol = AllSigmas ;
% NetMagnetization = CalculateMagnetization(AllSigmas.') ;
% Applied_field = AllApplied_fields ;

%% What to do about Eigs?
% DoCalcEigs = 0 ;
% if exist('CalcEigs','var')
%     if CalcEigs
%         DoCalcEigs = 1 ;
%     end
% end

% if DoCalcEigs
%     save([ThisPreSetFileName(1:end-4),'solution.mat'],'tt','SigmaSol','AllVV')
% else
    %     save([ThisPreSetFileName(1:end-4),'solution.mat'],'tt','SigmaSol')
    %     save([SimulationName(1:end-4),'solution.mat'],'AllSigmas','AllAppliedFields') ;
end
% save([SimulationName(1:end-4),'HystLoop.mat'],'tt','Applied_field','NetMagnetization') ;
%% Launch simulation
%
%[SigmaIN,Nold] = LandauLifshitzEvolve3D(ThisPreSetFileName) ;
%outFileName = 'OutRes01.mat' ;
%save(outFileName,'SigmaIN','Nold') ;
% LandauLifshitzEvolve3D(ThisPreSetFileName) ;

% load()