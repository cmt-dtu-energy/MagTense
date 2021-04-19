

function [x,Hmin,Kn] = iterate_H_const_perm_iso(mur)
    %param struct
    params = struct();
    %relative permeability
    params.mur = mur;

    %applied field vector
    params.Happ = [1,0,0];
    
    %demag tensor
    params.N = -1.*[[1./3.,0,0];[0,1./3.,0];[0,0,1./3.]];

   [Hmin,Kn] = Knorm_min_analy( params);
    
   K_fct = @(H) Knorm( H, params );
    
   [x,fval] = fminbnd( K_fct, 0, 2 );
   
end

%%returns the norm of K-vector, which is defined as
%% K = (Happ - H) * Happ_unit + H * (mur-1) * N dot Happ_unit
%% with Happ and H denoting the norms of the applied and local fields, 
%% respectively, and mur is the relative permeability and N the demag tensor
function Knorm = Knorm( H, params )

    %Applied field
    Happ = params.Happ;
    %relative permeability
    mur = params.mur;
    %demag tensor
    N = params.N;
    
    %norm of applied field
    Happ_nrm = sqrt(sum(Happ.^2));
    %unit vector of applied field
    Happ_un = Happ ./ Happ_nrm;
    
    %K-vector 
    K = (Happ_nrm - H) .* Happ_un + (H .* (mur-1) .* ( N * Happ_un' ) )';
    
    %norm of K-vector
    Knorm = sqrt(sum(K.^2));
end

function [Hmin,Knorm] = Knorm_min_analy( params)
     %Applied field
    Happ = params.Happ;
    %relative permeability
    mur = params.mur;
    %demag tensor
    N = params.N;
    
    %norm of applied field
    Happ_nrm = sqrt(sum(Happ.^2));
    %unit vector of applied field
    Happ_un = Happ ./ Happ_nrm;
    
    NHapp = N * Happ_un';
    
    v1 = ((mur-1) .* NHapp)' - Happ_un;
    
    v2 = Happ_un;
    Hmin = - Happ_nrm * dot(v1,v2)/sum(v1.^2);
    Knorm = Hmin^2 * sum(v1.^2) + Happ_nrm^2 * sum(v2.^2) + 2*Hmin*Happ_nrm * dot(v1,v2);
end