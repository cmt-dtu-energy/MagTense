function [E_dem_red, E_exc_red, E_ani_red, E_ext_red] = computeMagneticEnergy(solution,Vols,problem,Ms)

    mu0 = 4*pi*1e-7;
    
    %--- Compute the magnetic moment for a general mesh
    m = zeros(size(solution.M));
    for i = 1:length(solution.M(:,1,1,1))
        for j = 1:length(solution.M(1,1,:,1))
            for k = 1:length(solution.M(1,1,1,:))
                m(i,:,j,k) = mu0*Ms.*(Vols'.*solution.M(i,:,j,k));
            end
        end
    end
    
    
    %--- Calculate the energy terms normalized to the volume elements      
    E_exc = sum((1/2)*(m(:,:,:,1).*solution.H_exc(:,:,:,1) + m(:,:,:,2).*solution.H_exc(:,:,:,2) + m(:,:,:,3).*solution.H_exc(:,:,:,3)),2) ; % [J]       
    E_ext = sum(      (m(:,:,:,1).*solution.H_ext(:,:,:,1) + m(:,:,:,2).*solution.H_ext(:,:,:,2) + m(:,:,:,3).*solution.H_ext(:,:,:,3)),2) ; % [J]
    E_dem = sum((1/2)*(m(:,:,:,1).*solution.H_dem(:,:,:,1) + m(:,:,:,2).*solution.H_dem(:,:,:,2) + m(:,:,:,3).*solution.H_dem(:,:,:,3)),2) ; % [J]
    E_ani = sum((1/2)*(m(:,:,:,1).*solution.H_ani(:,:,:,1) + m(:,:,:,2).*solution.H_ani(:,:,:,2) + m(:,:,:,3).*solution.H_ani(:,:,:,3)),2) ; % [J]
    
    %--- The anistropy energy must be added the constant volume term of the integral of K0 over the volume 
    %--- See "Nonlinear Magnetization Dynamics in Thin-films and Nanoparticles" by Massimiliano d'Aquino, Eq. (1.43)
    E_ani = E_ani + sum(problem.K0.*Vols); 
    
    %--- Calcute the energy densities
    E_exc_dens = E_exc/prod(problem.grid_L);  % [J/m^3]
    E_ext_dens = E_ext/prod(problem.grid_L);  % [J/m^3]
    E_dem_dens = E_dem/prod(problem.grid_L);  % [J/m^3]
    E_ani_dens = E_ani/prod(problem.grid_L);  % [J/m^3]
    
    Km = (1/2)*(Ms^2)*mu0;  % [J/m^3] 
    %--- See Eq. (29) from "Standard Problems in Micromagnetics", Donald Porter and Michael Donahue (2018)
    
    %--- Reduced energies, i.e. normalized by Km
    E_exc_red = E_exc_dens./Km;
    E_ext_red = E_ext_dens./Km;
    E_dem_red = E_dem_dens./Km;
    E_ani_red = E_ani_dens./Km;

end