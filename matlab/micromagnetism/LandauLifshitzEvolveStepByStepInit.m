function [SigmaSol, AllVV] = LandauLifshitzEvolveStepByStepInit(ProblemSetupStruct,InteractionMatrices,Sigma)
mu0 = 4*pi*1e-7;

%--- Evaluate all variables in the ProblemSetupStruct
names = fieldnames(ProblemSetupStruct);
for i=1:length(names)
    eval([names{i} '=ProblemSetupStruct.' names{i} ';']);
end

for k=1:numel(t_explicit)
    disp(['  ',num2str(k),'/',num2str(numel(t_explicit))]) ;
    
    ProblemSetupStructExplicit = ProblemSetupStruct ;
    ProblemSetupStructExplicit.Hext(:,1)   = linspace(ProblemSetupStructExplicit.t(1),ProblemSetupStructExplicit.t(end),length(ProblemSetupStruct.Hext(:,2)));
    ProblemSetupStructExplicit.Hext(:,2:4) = ProblemSetupStruct.Hext(k,2:4).*ones(length(ProblemSetupStruct.Hext(:,2)),1);
    
    if k==1
        [SigmaSol,VV] = LandauLifshitzEvolveCombined(ProblemSetupStructExplicit,InteractionMatrices,Sigma) ;
    else
        [SigmaSol,VV] = LandauLifshitzEvolveCombined(ProblemSetupStructExplicit,InteractionMatrices,AllSigmas(k-1,:)) ;
    end

    AllSigmas(k,:) = SigmaSol(end,:) ;
    if CalcEigenvalue
        AllVV(k,:) = VV(:) ;
    else
        AllVV = [];
    end
end

SigmaSol = AllSigmas ;
end
