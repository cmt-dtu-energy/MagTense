function [SigmaSol, ImplicitFailed, AllVV] = LandauLifshitzEvolveImplicitSteps(ProblemSetupStruct,InteractionMatrices,Sigma)

%--- Evaluate all variables in the ProblemSetupStruct
names = fieldnames(ProblemSetupStruct);
for i=1:length(names)
    eval([names{i} '=ProblemSetupStruct.' names{i} ';']);
end

tt = t ;
originalHsX = HsX ;
originalHsY = HsY ;
originalHsZ = HsZ ;

for k=1:numel(tt)
    disp(['  ',num2str(k),'/',num2str(numel(tt))]) ;
    
    if k==1
        ProblemSetupStruct.UseExplicitSolver = 1;
        ProblemSetupStruct.UseImplicitSolver = 0;
        ProblemSetupStruct.HsX = @(t) 0.*t + originalHsX(tt(k)) ;
        ProblemSetupStruct.HsY = @(t) 0.*t + originalHsY(tt(k)) ;
        ProblemSetupStruct.HsZ = @(t) 0.*t + originalHsZ(tt(k)) ;
        ProblemSetupStruct.MaxT = MaxTstep ;
        ProblemSetupStruct.nT = 2;
        
        [SigmaSol,~,VV] = LandauLifshitzEvolveCombined(ProblemSetupStruct,InteractionMatrices,Sigma) ;
    else
        
        try
            ProblemSetupStruct.UseExplicitSolver = 0;
            ProblemSetupStruct.UseImplicitSolver = 1;
            ProblemSetupStruct.HsX = originalHsX ;
            ProblemSetupStruct.HsY = originalHsY ;
            ProblemSetupStruct.HsZ = originalHsZ ;
            ProblemSetupStruct.t = [tt(k-1),tt(k)] ;
            
            [SigmaSol, ~,VV]  = LandauLifshitzEvolveCombined(ProblemSetupStruct,InteractionMatrices,AllSigmas(k-1,:));
            
            ImplicitFailed(k-1) = 0 ;
        catch % Implicit Solver cannot go past this point
            ProblemSetupStruct.UseExplicitSolver = 1;
            ProblemSetupStruct.UseImplicitSolver = 0;
            ProblemSetupStruct.HsX = @(t) 0.*t + originalHsX(tt(k)) ;
            ProblemSetupStruct.HsY = @(t) 0.*t + originalHsY(tt(k)) ;
            ProblemSetupStruct.HsZ = @(t) 0.*t + originalHsZ(tt(k)) ;
            ProblemSetupStruct.MaxT = MaxTstep ;
            ProblemSetupStruct.nT = 2;
            
            [SigmaSol,~,VV] = LandauLifshitzEvolveCombined(ProblemSetupStruct,InteractionMatrices,AllSigmas(k-1,:)) ;
            
            ImplicitFailed(k-1) = 1 ;
        end
    end
    
    AllSigmas(k,:) = SigmaSol(end,:) ;
    
    if ProblemSetupStruct.CalcEigenvalue
        AllVV(k,:) = VV(:) ;
    else
        AllVV = [];
    end
end

SigmaSol = AllSigmas ;
end