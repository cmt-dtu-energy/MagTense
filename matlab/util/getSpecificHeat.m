
%%Returns the specific heat as a function of T at constant field
%%T is an input array with dimensions [n,1] containing temperature
%%S is an input array with dimensions [n,m]
function [cp] = getSpecificHeat( T, S )
    
    dT = 0.001;

	T_low = T - dT;
    T_high = T + dT;
    
    S_low = interp1( T, S, T_low, 'linear','extrap');
    S_high = interp1( T, S, T_high, 'linear', 'extrap' );
    
    cp = zeros( size(S) );
    
    for i=1:length(cp(1,:))
       cp(:,i) = T .* ( S_high(:,i) - S_low(:,i) ) ./ (T_high - T_low); 
    end
end