function [SigmaSol, AppliedField, AllVV] = LandauLifshitzEvolveStepByStepInit(ProblemSetupStruct,InteractionMatrices,Sigma)

%--- Evaluate all variables in the ProblemSetupStruct
names = fieldnames(ProblemSetupStruct);
for i=1:length(names)
    eval([names{i} '=ProblemSetupStruct.' names{i} ';']);
end

originalHsX = HsX ;
originalHsY = HsY ;
originalHsZ = HsZ ;

for k=1:numel(tt)
    disp(['  ',num2str(k),'/',num2str(numel(tt))]) ;
    ProblemSetupStruct.HsX = @(t) 0.*t + originalHsX(tt(k)) ;
    ProblemSetupStruct.HsY = @(t) 0.*t + originalHsY(tt(k)) ;
    ProblemSetupStruct.HsZ = @(t) 0.*t + originalHsZ(tt(k)) ;

    if k==1
        [SigmaSol,~,VV] = LandauLifshitzEvolveCombined(ProblemSetupStruct,InteractionMatrices,Sigma) ;
    else
        [SigmaSol,~,VV] = LandauLifshitzEvolveCombined(ProblemSetupStruct,InteractionMatrices,AllSigmas(k-1,:)) ;
    end

    AllSigmas(k,:) = SigmaSol(end,:) ;
    if CalcEigenvalue
        AllVV(k,:) = VV(:) ;
    else
        AllVV = [];
    end
end

Hfact = 1/MU0 ; % AppliedField is provided in T
AppliedField = [Hfact*originalHsX(tt), Hfact*originalHsY(tt),Hfact*originalHsZ(tt)+0.*tt] ;
SigmaSol = AllSigmas ;
end
