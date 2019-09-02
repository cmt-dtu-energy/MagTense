function [M_cl,S_cl, M_ht, S_ht] = getBR_MF_properties( T, H, p, Tc, eta, material )

[J,g,Ns,muB,mu0,kB,rhos,Na,mol,Tdeb,gamma,kappa,Na_c, T0_] = getMFTParams( material );

BR_params = struct();
BR_params.kappa = kappa;
BR_params.eta = eta;
BR_params.J = J;
BR_params.g = g;
BR_params.Ns = Ns;
BR_params.rho = rhos;
BR_params.signBeta = int32(1);


Tin = zeros( length(T), length(H) );
Hin = zeros( length(T), length(H) );

for i=1:length(H)
   Tin(:,i) = T; 
end
for i=1:length(T)
    Hin(i,:) = H;
end

Tin = reshape( Tin, [numel(Tin), 1] );
Hin = reshape( Hin, [numel(Hin), 1] );

T0 = zeros(numel(Tin),1);
T0(:) = Tc;


%    BR_params.eta = eta/T0_c * (T0_c+dTc(j));
n = numel(Tin);        
[M_,S_,deltaF] = BR_LIB_Matlab( BR_params, int32(n), Tin, Hin, p, T0 );

M_ = reshape(M_,[length(T),length(H),2]);
S_ = reshape(S_,[length(T),length(H),2]);

M_cl = M_(:,:,2);
S_cl = S_(:,:,2);

M_ht = M_(:,:,1);
S_ht = S_(:,:,1);
    
end

