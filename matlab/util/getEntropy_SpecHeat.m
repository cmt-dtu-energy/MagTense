
function [S,cp] = getEntropy_SpecHeat( T, H, HystMap, stateFcn )

n = length( stateFcn );

cp = zeros( n, 1 );
S = zeros( n, 1 );


for i=1:n

    S_cool = interp2( stateFcn(i).T, stateFcn(i).H, stateFcn(i).SCool, T(i), H(i) );

    S_heat = interp2( stateFcn.T(i), stateFcn.H(i), stateFcn(i).SHeat, T(i), H(i) );


    cp_cool = interp2( stateFcn(i).T, stateFcn(i).H, stateFcn(i).CCool, T(i), H(i) );

    cp_heat = interp2( stateFcn.T(i), stateFcn.H(i), stateFcn(i).CHeat, T(i), H(i) );

    
    %means use the cooling curve (paramagnetic)
    if HystMap(i) == 0
        cp( i ) = cp_cool;
        S( I ) = S_cool;
    elseif HystMap(i) == 1
     %means use the heating curve (ferromagnetic)
        cp( i ) = cp_heat;
        S( i ) = S_heat;
    end


end