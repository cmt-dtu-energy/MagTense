function status = MyOutputFunction(t,SigmaSol,flag,UseImplicitSolver,dSigmaMagTol)
if isempty(flag) 
    if UseImplicitSolver
        global TheData
        LastP = TheData.p ;
%         disp([num2str(t(end)),'  Last P : ',num2str(LastP)])

        if (LastP > 0)
            status = 0 ; % was 1
        else
            status = 0 ; 
        end
    else
        try
            global TheData
            LastRMS = TheData.dSigmaRMS ;
        catch
            LastRMS = inf ;
        end

%         disp(['Log10dS ',num2str(log10(LastRMS))]) ;
        if (log10(LastRMS) < dSigmaMagTol)
            status = 1 ; % was 1
        else
            status = 0 ; 
        end
    end
end
end