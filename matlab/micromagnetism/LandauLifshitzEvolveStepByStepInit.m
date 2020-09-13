function [SigmaSol, AllVV] = LandauLifshitzEvolveStepByStepInit(ProblemSetupStruct,InteractionMatrices,Sigma)
mu0 = 4*pi*1e-7;

%--- Evaluate all variables in the ProblemSetupStruct
names = fieldnames(ProblemSetupStruct);
for i=1:length(names)
    eval([names{i} '=ProblemSetupStruct.' names{i} ';']);
end

for k=1:length(ProblemSetupStruct.Hext(:,1))
    disp(['  ',num2str(k),'/',num2str(length(ProblemSetupStruct.Hext(:,1)))]) ;
    
    ProblemSetupStructExplicit = ProblemSetupStruct ;
    ProblemSetupStructExplicit.Hext(:,1)   = linspace(ProblemSetupStructExplicit.t(1),ProblemSetupStructExplicit.t(end),length(ProblemSetupStruct.Hext(:,2)));
    ProblemSetupStructExplicit.Hext(:,2:4) = ProblemSetupStruct.Hext(k,2:4).*ones(length(ProblemSetupStruct.Hext(:,2)),1);
    
    if k==1
        [SigmaSol,VV] = LandauLifshitzEvolveCombined(ProblemSetupStructExplicit,InteractionMatrices,Sigma) ;
    else
        Sigmastart = AllSigmas(k-1,:);        
        Sigmastart = Sigmastart+randn(size(Sigmastart))*1;
        
        NN = round(numel(Sigmastart)/3) ;
        SigmaXorg = Sigmastart(0*NN+[1:NN]) ;
        SigmaYorg = Sigmastart(1*NN+[1:NN]) ;
        SigmaZorg = Sigmastart(2*NN+[1:NN]) ;
        
        SigmaX = SigmaXorg./sqrt(SigmaXorg.^2+SigmaYorg.^2+SigmaZorg.^2);
        SigmaY = SigmaYorg./sqrt(SigmaXorg.^2+SigmaYorg.^2+SigmaZorg.^2);
        SigmaZ = SigmaZorg./sqrt(SigmaXorg.^2+SigmaYorg.^2+SigmaZorg.^2);
        
        Sigmastart = [SigmaX SigmaY SigmaZ];
        
        [SigmaSol,VV] = LandauLifshitzEvolveCombined(ProblemSetupStructExplicit,InteractionMatrices,Sigmastart) ;
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
