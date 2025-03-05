function Cylinder_and_cylindrical_slice_test_limits

    clearvars
    close all
    
    info.use_matlab = 1;
    
    %--- Magnetization directions
    M(1) = (-1 + 2.*rand(1,1))*10;
    M(2) = (-1 + 2.*rand(1,1))*10;
    M(3) = (-1 + 2.*rand(1,1))*10;
    
    [r_t, phi_t, z_t, Ri_t, Ro_t, phi1_t, phi2_t, z1_t, z2_t, slice_t] = test_cases();
    
    for i = 1:length(r_t)
        r    = r_t(i);
        phi  = phi_t(i);
        z    = z_t(i);
        Ri   = Ri_t(i); 
        Ro   = Ro_t(i);
        phi1 = phi1_t(i); 
        phi2 = phi2_t(i);
        z1   = z1_t(i);
        z2   = z2_t(i);
        slice = slice_t(i);
        
        [Psi_tot, Psi_z_arr, Psi_r_arr, Psi_phi_arr] = MagTense_Validation_cylinder_and_clyndrical_slice_potential(Ro, Ri, phi1, phi2, z1, z2, M, r, phi, z, slice, info);
            
        disp([Psi_tot(1,1) Psi_z_arr(1,1) Psi_r_arr(1,1) Psi_phi_arr(1,1)]);
        disp([Psi_tot(1,2) Psi_z_arr(1,2) Psi_r_arr(1,2) Psi_phi_arr(1,2)]);
        disp([Psi_tot(1,3) Psi_z_arr(1,3) Psi_r_arr(1,3) Psi_phi_arr(1,3)]);
        
        accuracy = abs((Psi_tot(1,1)-Psi_tot(1,3))/Psi_tot(1,1)*100);
        if (accuracy < 1e-4)
            str = ['<a href="">PASS - difference: ' num2str(accuracy) '% </a>'];
            disp(str)
        else
            str = ['FAIL - difference: ' num2str(accuracy) '%\n'];
            fprintf(2,str)
            fprintf(2,'\n')
        end
    end
end


%% Test limits and singularities
function [r, phi, z, Ri, Ro, phi1, phi2, z1, z2, slice] = test_cases()

    %--- r = 0, slice
    r(1)      = 0;
    phi(1)    = 2*pi*rand(1,1);
    z(1)      = rand(1,1);
    Ri(1)     = rand(1,1);
    Ro(1)     = Ri(1) + (1-Ri(1))*rand(1,1);
    phi1(1)   = 2*pi*rand(1,1);
    phi2(1)   = phi1(1) + (2*pi-phi1(1))*rand(1,1);
    z1(1)     = -rand(1,1);
    z2(1)     = +rand(1,1);
    slice(1) = 1;

    %--- r = 0, full
    r(2)      = 0;
    phi(2)    = 2*pi*rand(1,1);
    z(2)      = rand(1,1);
    Ri(2)     = 0;
    Ro(2)     = rand(1,1);
    phi1(2)   = 0;
    phi2(2)   = 2*pi;
    z1(2)     = -rand(1,1);
    z2(2)     = +rand(1,1);
    slice(2) = 0;

    %--- z-z' = 0
    r(3)      = rand(1,1);
    phi(3)    = 2*pi*rand(1,1);
    z(3)      = -rand(1,1);
    Ri(3)     = rand(1,1);
    Ro(3)     = Ri(3) + (1-Ri(3))*rand(1,1);
    phi1(3)   = 2*pi*rand(1,1);
    phi2(3)   = phi1(3) + (2*pi-phi1(3))*rand(1,1);
    z1(3)     = z(3);
    z2(3)     = +rand(1,1);
    slice(3) = 1;

    %--- z-z' = 0, full
    r(4)      = rand(1,1);
    phi(4)    = 2*pi*rand(1,1);
    z(4)      = -rand(1,1);
    Ri(4)     = 0;
    Ro(4)     = rand(1,1);
    phi1(4)   = 0;
    phi2(4)   = 2*pi;
    z1(4)     = z(4);
    z2(4)     = +rand(1,1);
    slice(4) = 0;

    %--- r = 0 && z-z' = 0
    r(5)      = 0;
    phi(5)    = 2*pi*rand(1,1);
    z(5)      = -rand(1,1);
    Ri(5)     = rand(1,1);
    Ro(5)     = Ri(5) + (1-Ri(5))*rand(1,1);
    phi1(5)   = 2*pi*rand(1,1);
    phi2(5)   = phi1(5) + (2*pi-phi1(5))*rand(1,1);
    z1(5)     = z(5);
    z2(5)     = +rand(1,1);
    slice(5) = 1;

    %--- r = 0 && z-z' = 0, full
    r(6)      = 0;
    phi(6)    = 2*pi*rand(1,1);
    z(6)      = -rand(1,1);
    Ri(6)     = 0;
    Ro(6)     = rand(1,1);
    phi1(6)   = 0;
    phi2(6)   = 2*pi;
    z1(6)     = z(6);
    z2(6)     = +rand(1,1);
    slice(6) = 0;
end