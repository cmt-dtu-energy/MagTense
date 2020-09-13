function [Smag] = getSmagBR(T,H,M,eta,T0,material)
    
    [J,g,Ns,muB,mu0,kB,rho,Na,mol,Tdeb,gammaSommerfeld,kappa] = getMFTParams( material );

    indGTZ = find(M>0);
    indEQZ = find(M==0);
    %Get Smag from Magnetization
    A = (2*J+1)/(2*J);
    B = 1/(2*J);
    NV = Ns * rho;
    gamma = g * muB;
    beta = sqrt( 1/40. * eta * 1 / ( NV * kB * T0 * kappa ) * ( (2*J+1)^4-1 )/(J*(J+1))^2 );
    lambda0 = T0 * 3 * kB / ( J * ( J + 1 ) * gamma^2 ) * 1/NV;
    
    x = zeros(size(M));
    
    x(indGTZ) = gamma * J/(kB*T) .* ( H(indGTZ) + lambda0.*M(indGTZ) + 1/2*lambda0^2*beta^2*kappa.*M(indGTZ).^3 );
    
    x(indEQZ) = 0;
    Br = zeros(size(M));
    
    Br(indGTZ) = A .* coth( A.*x(indGTZ) ) - B.*coth(B.*x(indGTZ));
    
    Smag = zeros( size(M) );
        
    Smag(indGTZ) = kB*NV .* ( log( sinh(A.*x(indGTZ))./sinh(B.*x(indGTZ)) )  - x(indGTZ) .* Br(indGTZ) );
    
    Smag(indEQZ) = kB*NV .* log( 2 * J + 1 );

end