function [L,ind1] = hasLatentHeat( T, S )

    err = 0.99;
    
    T1 = (( T(2:end) + T(1:end-1) ) .* 0.5);
    c1 = T1 .* ( S(2:end) - S(1:end-1) )./ ( T(2:end) - T(1:end-1) );

    cMax1 = max(c1);
    
    %Downsample T to make a coarser grid
    T_ds = T(1:2:end);
    S_ds = S(1:2:end);
    
    T1_ds = (( T_ds(2:end) + T_ds(1:end-1) ) .* 0.5);
    c1_ds = T1_ds .* ( S_ds(2:end) - S_ds(1:end-1) )./ ( T_ds(2:end) - T_ds(1:end-1) );
    
    cMax_ds = max(c1_ds);
    
    if cMax_ds / cMax1 <= err
        %Latent heat, find it and return it
        ind1 = find( c1 == cMax1 );
        ind2 = ind1 + 1;
        L = S(ind2) - S(ind1);
    else
        %Zero latent heat
        L = 0;
        ind1 = 1;
    end
    
end