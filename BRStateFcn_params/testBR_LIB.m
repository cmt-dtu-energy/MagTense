clear vars
%close all

mu0 = 4*pi*1e-7;

[J,g,Ns,muB,mu0,kB,rhos,Na,mol,Tdeb,gamma,kappa,Na_c, T0] = getMFTParams( 'LaFeSi' );

eta = 1.4;

BR_params = struct();
BR_params.kappa = kappa;
BR_params.eta = eta;
BR_params.J = J;
BR_params.g = g;
BR_params.Ns = Ns;
BR_params.rho = rhos;
BR_params.signBeta = int32(1);

T = 100:1:320;
H = [0.1 1.0];%0.1:0.1:1;

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

%T0, the Tc if no volume change
T0 = zeros(numel(Tin));

%pressure in Pa
p = 0;
n_Tc = 1000;
sigma = 1;%K

n = numel(Tin);

dTc = random( 'normal', 0, sigma, [n_Tc,1] );

M1 = zeros( length(T), length(H) );
M2 = zeros( length(T), length(H) );

S1 = zeros( length(T), length(H) );
S2 = zeros( length(T), length(H) );

for i=1:n_Tc
    T0(:) = 220 + dTc(i);

    [M_,S_,deltaF] = BR_LIB_Matlab( BR_params, int32(n), Tin, Hin, p, T0 );

    M_ = reshape(M_,[length(T),2,length(H)]);
    S_ = reshape(S_,[length(T),2,length(H)]);
    
    M1 = M1 + squeeze(M_(:,1,:));
    M2 = M2 + squeeze(M_(:,2,:));
    
    S1 = S1 + squeeze(S_(:,1,:));
    S2 = S2 + squeeze(S_(:,2,:));
    
end

M1 = M1 ./ n_Tc;
M2 = M2 ./ n_Tc;

S1 = S1 ./ n_Tc;
S2 = S2 ./ n_Tc;

plot(T,-(S1(:,1)-S1(:,2)))
plot(T,-(S2(:,1)-S2(:,2)),'--')

