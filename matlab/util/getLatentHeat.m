function [L_PM,ind_PM,L_FM,ind_FM] = getLatentHeat( TT, HH, pp, S_pm, S_fm )
    
    L_PM = zeros( length(HH), length(pp) );
    ind_PM = L_PM;
    L_FM = L_PM;
    ind_FM = L_PM;
    for i=1:length(HH)
        for j=1:length(pp)
            [L_PM(i,j),ind_PM(i,j)] = hasLatentHeat( TT, S_pm(i,:,j) );
            [L_FM(i,j),ind_FM(i,j)] = hasLatentHeat( TT, S_fm(i,:,j) );
            %convert from J/kgK to J/kg
            L_PM(i,j) = L_PM(i,j) .* TT(ind_PM(i,j))';
            L_FM(i,j) = L_FM(i,j) .* TT(ind_FM(i,j))';
    
        end
    end
    
end