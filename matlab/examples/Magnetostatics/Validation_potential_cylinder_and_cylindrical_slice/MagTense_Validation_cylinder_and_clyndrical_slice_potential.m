function [Psi_tot, Psi_z_arr, Psi_r_arr, Psi_phi_arr] = MagTense_Validation_cylinder_and_clyndrical_slice_potential(Ro, Ri, phi1, phi2, z1, z2, M, r_arr, phi_arr, z_arr, slice, info)

% Input parameters
% Ro      The outer radius of the cylinder slice or full cylinder
% Ri      The innter radius of the cylinder slice. For the full cylinder this is zero
% phi1    The start angle of the cylinder slice. For the full cylinder this is not a parameter
% phi2    The end angle of the cylinder slice. For the full cylinder this is not a parameter
% z1      The lower z coordinate of the cylinder slice or cylinder
% z2      The upper z coordinate of the cylinder slice or cylinder
% M       The magnetization in cartesian coordinates, specified as a vector of [Mx My Mz]
% r_arr   The radial coordinate of the points at which the scalar potential should be evaluated
% phi_arr The angular coordinate of the points at which the scalar potential should be evaluated
% z_arr   The height coordinate of the points at which the scalar potential should be evaluated
% slice   The parameters specifying if the tile is a cylinder slice or a full cylinder
% info    The struct used to specify options on how the calculations must be done

if ~isfield(info,'evaluate_old_expressions')
    info.evaluate_old_expressions = false;              %--- Evaluate the original expressions for the integral coefficients
    info.evaluate_small_q_expressions = false;          %--- Evaluate the small q expressions for the integral coefficients
end

%--- Get the cylindrical coordinates
[phi_M,theta_M,M_mag] = cart2sph(M(1),M(2),M(3));
theta_M = pi/2 - theta_M;                               %--- This is necessary because Matlab defines the elevation angle from the xy-place, which is not the ISO-standard

%--- Preallocating the arrays for storing the magnetostatic scalar potential
Psi_z_arr   = zeros(length(phi_arr),3);
Psi_r_arr   = zeros(length(phi_arr),3);
Psi_phi_arr = zeros(length(phi_arr),3);
Psi_tot     = zeros(length(phi_arr),3);

%--- Loop over the coordinates and calculate the magnetic scalar potential
for i = 1:length(phi_arr)
     phi = phi_arr(i);
     r   = r_arr(i);
     z   = z_arr(i);

    %%-------------------------------------------------------------------------
    %%------------- The end face integral in r',phi' (constant z') ------------
    %%-------------------------------------------------------------------------   
    
    %--- Double integral
    Psi_z1 = -1/(4*pi)*integral2(@(rp,phip) z0_full_func(rp, phip, z1, r, phi, z, M_mag, theta_M, phi_M), Ri, Ro, phi1, phi2);
    Psi_z2 =  1/(4*pi)*integral2(@(rp,phip) z0_full_func(rp, phip, z2, r, phi, z, M_mag, theta_M, phi_M), Ri, Ro, phi1, phi2);

    Psi_z_arr(i,1) = Psi_z1+Psi_z2;

    if (slice == 1)
        if (info.evaluate_small_q_expressions)
            Psi_z1_ana   = -((analytical_z0(Ro, phi2, z1, r, phi, z, M_mag, theta_M, phi_M, slice, info) - analytical_z0(Ro, phi1, z1, r, phi, z, M_mag, theta_M, phi_M, slice, info)) ...
                           - (analytical_z0(Ri, phi2, z1, r, phi, z, M_mag, theta_M, phi_M, slice, info) - analytical_z0(Ri, phi1, z1, r, phi, z, M_mag, theta_M, phi_M, slice, info)));
        
            Psi_z2_ana   =   (analytical_z0(Ro, phi2, z2, r, phi, z, M_mag, theta_M, phi_M, slice, info) - analytical_z0(Ro, phi1, z2, r, phi, z, M_mag, theta_M, phi_M, slice, info)) ...
                           - (analytical_z0(Ri, phi2, z2, r, phi, z, M_mag, theta_M, phi_M, slice, info) - analytical_z0(Ri, phi1, z2, r, phi, z, M_mag, theta_M, phi_M, slice, info));
        
            Psi_z_arr(i,2) = Psi_z1_ana + Psi_z2_ana;
        end

        [~,Q_1] = analytical_z0(Ro, phi2, z1, r, phi, z, M_mag, theta_M, phi_M, slice, info);
        [~,Q_2] = analytical_z0(Ro, phi1, z1, r, phi, z, M_mag, theta_M, phi_M, slice, info);
        [~,Q_3] = analytical_z0(Ri, phi2, z1, r, phi, z, M_mag, theta_M, phi_M, slice, info);
        [~,Q_4] = analytical_z0(Ri, phi1, z1, r, phi, z, M_mag, theta_M, phi_M, slice, info);

        [~,Q_5] = analytical_z0(Ro, phi2, z2, r, phi, z, M_mag, theta_M, phi_M, slice, info);
        [~,Q_6] = analytical_z0(Ro, phi1, z2, r, phi, z, M_mag, theta_M, phi_M, slice, info);
        [~,Q_7] = analytical_z0(Ri, phi2, z2, r, phi, z, M_mag, theta_M, phi_M, slice, info);
        [~,Q_8] = analytical_z0(Ri, phi1, z2, r, phi, z, M_mag, theta_M, phi_M, slice, info);

        N_z1 = -((Q_1 - Q_2) - (Q_3 - Q_4));
        N_z2 = +((Q_5 - Q_6) - (Q_7 - Q_8));

        Psi_z_arr(i,3) = sum(N_z1.*M) + sum(N_z2.*M);
    end

    if (slice == 0)
        if (info.evaluate_small_q_expressions)
            Psi_z1_ana_full_cyl = -(- 2*analytical_z0(Ro, [], z1, r, phi, z, M_mag, theta_M, phi_M, slice, info) ...
                                    + 2*analytical_z0(Ri, [], z1, r, phi, z, M_mag, theta_M, phi_M, slice, info));
    
            Psi_z2_ana_full_cyl =   - 2*analytical_z0(Ro, [], z2, r, phi, z, M_mag, theta_M, phi_M, slice, info) ...
                                    + 2*analytical_z0(Ri, [], z2, r, phi, z, M_mag, theta_M, phi_M, slice, info);
    
            Psi_z_arr(i,2) = Psi_z1_ana_full_cyl + Psi_z2_ana_full_cyl;
        end

        [~,Q_1] = analytical_z0(Ro, [], z1, r, phi, z, M_mag, theta_M, phi_M, slice, info);
        [~,Q_2] = analytical_z0(Ri, [], z1, r, phi, z, M_mag, theta_M, phi_M, slice, info);
        
        [~,Q_3] = analytical_z0(Ro, [], z2, r, phi, z, M_mag, theta_M, phi_M, slice, info);
        [~,Q_4] = analytical_z0(Ri, [], z2, r, phi, z, M_mag, theta_M, phi_M, slice, info);

        N_z = -(-2*Q_1 + 2*Q_2) + (-2*Q_3 + 2*Q_4);
        Psi_z_arr(i,3) = sum(N_z.*M);
    end

    
    %%-------------------------------------------------------------------------
    %%------------- The tube integral in z',phi' (constant r') ----------------
    %%-------------------------------------------------------------------------
    
    %--- Double integral
    Psi_r_1  =  -1/(4*pi)*integral2(@(phip,zp) r0_full_func(phip, zp, Ri, r, phi, z, M_mag, theta_M, phi_M), phi1, phi2, z1, z2);
    Psi_r_2  =   1/(4*pi)*integral2(@(phip,zp) r0_full_func(phip, zp, Ro, r, phi, z, M_mag, theta_M, phi_M), phi1, phi2, z1, z2);

    Psi_r_arr(i,1) = Psi_r_1 + Psi_r_2;


    %--- Full analytical
    if (slice == 1)
        if (info.evaluate_small_q_expressions)
            Psi_r1_ana = -((analytical_r0(Ri, phi2, z2, r, phi, z, M_mag, theta_M, phi_M, slice, info) - analytical_r0(Ri, phi1, z2, r, phi, z, M_mag, theta_M, phi_M, slice, info)) ...
                          -(analytical_r0(Ri, phi2, z1, r, phi, z, M_mag, theta_M, phi_M, slice, info) - analytical_r0(Ri, phi1, z1, r, phi, z, M_mag, theta_M, phi_M, slice, info)));
    
            Psi_r2_ana = (analytical_r0(Ro, phi2, z2, r, phi, z, M_mag, theta_M, phi_M, slice, info) - analytical_r0(Ro, phi1, z2, r, phi, z, M_mag, theta_M, phi_M, slice, info)) ...
                        -(analytical_r0(Ro, phi2, z1, r, phi, z, M_mag, theta_M, phi_M, slice, info) - analytical_r0(Ro, phi1, z1, r, phi, z, M_mag, theta_M, phi_M, slice, info));
    
            Psi_r_arr(i,2) = Psi_r1_ana + Psi_r2_ana;
        end
        
        [~, Q_1] = analytical_r0(Ri, phi2, z2, r, phi, z, M_mag, theta_M, phi_M, slice, info);
        [~, Q_2] = analytical_r0(Ri, phi1, z2, r, phi, z, M_mag, theta_M, phi_M, slice, info);
        [~, Q_3] = analytical_r0(Ri, phi2, z1, r, phi, z, M_mag, theta_M, phi_M, slice, info);
        [~, Q_4] = analytical_r0(Ri, phi1, z1, r, phi, z, M_mag, theta_M, phi_M, slice, info);
        
        [~, Q_5] = analytical_r0(Ro, phi2, z2, r, phi, z, M_mag, theta_M, phi_M, slice, info);
        [~, Q_6] = analytical_r0(Ro, phi1, z2, r, phi, z, M_mag, theta_M, phi_M, slice, info);
        [~, Q_7] = analytical_r0(Ro, phi2, z1, r, phi, z, M_mag, theta_M, phi_M, slice, info);
        [~, Q_8] = analytical_r0(Ro, phi1, z1, r, phi, z, M_mag, theta_M, phi_M, slice, info);

        N_r1 = -((Q_1 - Q_2) - (Q_3 - Q_4));
        N_r2 = +((Q_5 - Q_6) - (Q_7 - Q_8));

        Psi_r_arr(i,3) = sum(N_r1.*M) + sum(N_r2.*M);
    end

    if (slice == 0)
        if (info.evaluate_small_q_expressions)
            Psi_r2_ana_full_cyl = - 2*analytical_r0(Ro, [], z2, r, phi, z, M_mag, theta_M, phi_M, slice, info) ...
                                  + 2*analytical_r0(Ro, [], z1, r, phi, z, M_mag, theta_M, phi_M, slice, info);
    
            Psi_r_arr(i,2) = Psi_r2_ana_full_cyl;
        end

        [~,Q_1] = analytical_r0(Ro, [], z2, r, phi, z, M_mag, theta_M, phi_M, slice, info);
        [~,Q_2] = analytical_r0(Ro, [], z1, r, phi, z, M_mag, theta_M, phi_M, slice, info);

        N_r = -2*Q_1 + 2*Q_2;
        Psi_r_arr(i,3) = sum(N_r.*M);
    end
   
    %%-------------------------------------------------------------------------
    %%------- The vertical face integral in z',r' (constant phi') -------------
    %%-------------------------------------------------------------------------
    
    if (slice == 1)
        %--- Double integral
        Psi_phi_1  =  -1/(4*pi)*integral2(@(rp,zp) phi0_full_func(rp, zp, phi1, r, phi, z, M_mag, theta_M, phi_M), Ri, Ro, z1, z2);
        Psi_phi_2  =   1/(4*pi)*integral2(@(rp,zp) phi0_full_func(rp, zp, phi2, r, phi, z, M_mag, theta_M, phi_M), Ri, Ro, z1, z2);

        Psi_phi_arr(i,1) = Psi_phi_1 + Psi_phi_2;

        if (info.evaluate_small_q_expressions)
            %--- Full analytical
            Psi_phi1_ana = -((analytical_phi0(Ro, phi1, z2, r, phi, z, M_mag, theta_M, phi_M) - analytical_phi0(Ro, phi1, z1, r, phi, z, M_mag, theta_M, phi_M)) ...
                            -(analytical_phi0(Ri, phi1, z2, r, phi, z, M_mag, theta_M, phi_M) - analytical_phi0(Ri, phi1, z1, r, phi, z, M_mag, theta_M, phi_M)));
        
            Psi_phi2_ana =   (analytical_phi0(Ro, phi2, z2, r, phi, z, M_mag, theta_M, phi_M) - analytical_phi0(Ro, phi2, z1, r, phi, z, M_mag, theta_M, phi_M)) ...
                            -(analytical_phi0(Ri, phi2, z2, r, phi, z, M_mag, theta_M, phi_M) - analytical_phi0(Ri, phi2, z1, r, phi, z, M_mag, theta_M, phi_M));
    
            Psi_phi_arr(i,2) = Psi_phi1_ana + Psi_phi2_ana;
        end

        [~, Q_1] = analytical_phi0(Ro, phi1, z2, r, phi, z, M_mag, theta_M, phi_M); 
        [~, Q_2] = analytical_phi0(Ro, phi1, z1, r, phi, z, M_mag, theta_M, phi_M);
        [~, Q_3] = analytical_phi0(Ri, phi1, z2, r, phi, z, M_mag, theta_M, phi_M);
        [~, Q_4] = analytical_phi0(Ri, phi1, z1, r, phi, z, M_mag, theta_M, phi_M);
    
        [~, Q_5] = analytical_phi0(Ro, phi2, z2, r, phi, z, M_mag, theta_M, phi_M);
        [~, Q_6] = analytical_phi0(Ro, phi2, z1, r, phi, z, M_mag, theta_M, phi_M);
        [~, Q_7] = analytical_phi0(Ri, phi2, z2, r, phi, z, M_mag, theta_M, phi_M);
        [~, Q_8] = analytical_phi0(Ri, phi2, z1, r, phi, z, M_mag, theta_M, phi_M);

        N_phi1 = -((Q_1 - Q_2) - (Q_3 - Q_4));
        N_phi2 = +((Q_5 - Q_6) - (Q_7 - Q_8));

        Psi_phi_arr(i,3) = sum(N_phi1.*M) + sum(N_phi2.*M);
    
    else
        Psi_phi_arr(i,1) = 0;
        Psi_phi_arr(i,2) = 0;
        Psi_phi_arr(i,3) = 0;
    end
    
    %%-------------------------------------------------------------------------
    %%------------------ Add the potential components -------------------------
    %%-------------------------------------------------------------------------
    
    %--- Double integral
    Psi_tot(i,1) = Psi_z_arr(i,1) + Psi_r_arr(i,1) + Psi_phi_arr(i,1);

    %--- Full analytical
    Psi_tot(i,2) = Psi_z_arr(i,2) + Psi_r_arr(i,2) + Psi_phi_arr(i,2);

    %--- Full analytical, demagnetization vector expression
    Psi_tot(i,3) = Psi_z_arr(i,3) + Psi_r_arr(i,3) + Psi_phi_arr(i,3);
end
end

%% Auxiliary functions
%%-------------------------------------------------------------------------
%%---------------- Functions for numerical integrals ----------------------
%%-------------------------------------------------------------------------
function Psi = z0_full_func(rp, phip, zp, r, phi, z, M_mag, theta_M, phi_M) 
    Psi = M_mag*cos(theta_M)*rp./(sqrt(r.^2+rp.^2-2.*r.*rp.*cos(phi-phip)+(z-zp).^2));
end

function Psi = r0_full_func(phip, zp, rp, r, phi, z, M_mag, theta_M, phi_M) 
    Psi = M_mag*sin(theta_M).*cos(phi_M-phip)*rp./(sqrt(r.^2+rp.^2-2.*r.*rp.*cos(phi-phip)+(z-zp).^2));
end

function Psi = phi0_full_func(rp, zp, phip, r, phi, z, M_mag, theta_M, phi_M) 
    %--- Notice there is no rp in the dr*drz in this integral
    Psi = M_mag*sin(theta_M).*sin(phi_M-phip)*1./(sqrt(r.^2+rp.^2-2.*r.*rp.*cos(phi-phip)+(z-zp).^2));
end

%%-------------------------------------------------------------------------
%%------------- The end face integral in r',phi' (constant z') ------------
%%-------------------------------------------------------------------------
function [totalint, totalint_V] = analytical_z0(rp, phip, zp, r, phi, z, M_mag, theta_M, phi_M, slice, info)
do_terms = 1;

a = -2*r*rp;
b = r^2+rp^2+(z-zp)^2;
c = r^2 - 3*rp^2 + 3*b;
d = sqrt(4*r^2*b - a^2);
e = -2*(rp^2 - b)^2 + 3/2*a^2 - 6*b*r^2;

[E_int, F_int, Pi1p_int, Pi1n_int, Pi2p_int, Pi2n_int, ~] = calculate_elliptic_integrals(slice, rp, phip, zp, r, phi, z, a, b, c, d, e, 'z', info.use_matlab);

if (r == 0)
    if (slice)
        part2 = -(phi-phip)*sqrt(rp^2 + (z-zp)^2);
    else
        part2 = -(E_int + F_int)*sqrt(rp^2 + (z-zp)^2);
    end
    do_terms = 0;
end

if (abs(z-zp) < 1e-3)
    if (slice)
        part2 = -abs(r - rp)*E_int + sign(r - rp)*(r + rp)*F_int + r*atanh((-rp + r*cos(phi - phip))/sqrt(b + a*cos(phi - phip)))*sin(phi - phip);
    else
        part2 = -abs(r - rp)*E_int + sign(r - rp)*(r + rp)*F_int;
    end
    do_terms = 0;
end

if ((abs(z-zp) < 1e-3) && r == 0)
    if (slice)
        part2 = -(phi - phip)*rp;
    else
        part2 = -pi*rp;
    end
    do_terms = 0;
end

if ((rp == 0) && (r ~= 0) && abs(z-zp) > 1e-3)
    % part2 = -sqrt(b)*pi - 2*r*sqrt(-1)*(Pi2p_int-Pi2n_int);
 
    np = (a + b)/(a - b)*(z-zp)^2/(+d - 2*r^2 - (z-zp)^2);
    nn = (a + b)/(a - b)*(z-zp)^2/(-d - 2*r^2 - (z-zp)^2);
    part2 = -pi*(sqrt(b) + r*(-(sqrt(np)-1)/(1-np)+(sqrt(nn)-1)/(1-nn)));

    do_terms = 0;
end

%--- All limit terms must have a factor of sign(phi-phip) in front if we do
%--- not use the trick of multiplying by two in the total integral
if (do_terms == 1)
    term0 = -2*sqrt(a + b)/r*E_int;
    
    if (slice)
        term1 = -atanh((rp - r*cos(phi - phip))/sqrt(b + a*cos(phi - phip)))*sin(phi - phip);
    else
        term1 = 0;
    end
    
    term2_X = A_factor(+1,rp, phip, zp, r, phi, z, a, b, c, d, e)*E_int;
    
    term3_X = A_factor(-1,rp, phip, zp, r, phi, z, a, b, c, d, e)*F_int;
    
    term4_X = B_factor(+1,rp, phip, zp, r, phi, z, a, b, c, d, e)*E_int;
    
    term5_X = B_factor(-1,rp, phip, zp, r, phi, z, a, b, c, d, e)*E_int;
    
    term6_X = C_factor(+1,rp, phip, zp, r, phi, z, a, b, c, d, e)*F_int;
    
    term7_X = C_factor(-1,rp, phip, zp, r, phi, z, a, b, c, d, e)*F_int;
    
    term8_X = D_factor(+1,rp, phip, zp, r, phi, z, a, b, c, d, e)*Pi1p_int;
    
    term9_X = D_factor(-1,rp, phip, zp, r, phi, z, a, b, c, d, e)*Pi1n_int;
    
    term10_X = E_factor(+1,rp, phip, zp, r, phi, z, a, b, c, d, e)*Pi2p_int;
    
    term11_X = -E_factor(+1,rp, phip, zp, r, phi, z, a, b, c, d, e)*Pi2n_int;
	
    part2_X = r*sum([term0 term1 term2_X term3_X term4_X term5_X term6_X term7_X term8_X term9_X term10_X term11_X]);
    
    if (info.evaluate_old_expressions)
        term2 = + (2*r*cos((phi - phip)/2)^2*ellipticE((phi - phip)/2, -((2*(-a))/(b - (-a))))*(b - (-a) + (b - a)*tan((phi - phip)/2)^2))/((r^2 + 2*r*rp + rp^2 - b - (-a))*sqrt(cos((phi - phip)/2)^(2)*(b - (-a) + (b + (-a))*tan((phi - phip)/2)^2))*sqrt((cos((phi - phip)/2)^(2)*(b - (-a) + (b + (-a))*tan((phi - phip)/2)^2))/(b - (-a))));
    
        term3 = - (2*r*(b + (-a))*cos((phi - phip)/2)^(2)*ellipticF((phi - phip)/2, -((2*(-a))/(b - (-a))))*(b - (-a) + (b + (-a))*tan((phi - phip)/2)^2))/((r^2 + 2*r*rp + rp^2 - b - (-a))*(b - (-a))*sqrt(cos((phi - phip)/2)^(2)*(b - (-a) + (b + (-a))*tan((phi - phip)/2)^2))*sqrt((cos((phi - phip)/2)^(2)*(b - (-a) + (b + (-a))*tan((phi - phip)/2)^2))/(b - (-a))));
        
        term4   = - (r*(r^2 - 2*r*rp - 3*rp^2 + 3*b - a - (2*r^(3)*rp - 2*rp^4 - 2*b^2 - 2*r*rp*(rp^2 - b - 2*(-a)) + r^(2)*(2*rp^2 - 6*b - (-a)) - b*(-a) - (-a)^2 + rp^(2)*(4*b + (-a)))/sqrt(4*r^(2)*b - 4*r*rp*(-a) + (-a)^2))*cos((phi - phip)/2)^(2)*ellipticE((phi - phip)/2, -((2*(-a))/(b - (-a))))*(b - (-a) + (b + (-a))*tan((phi - phip)/2)^2))/((r^2 + 2*r*rp + rp^2 - b - (-a))*(2*r^2 + 2*r*rp - (-a) + sqrt(4*r^(2)*b - 4*r*rp*(-a) + (-a)^2))*sqrt(cos((phi - phip)/2)^(2)*(b - (-a) + (b + (-a))*tan((phi - phip)/2)^2))*sqrt((cos((phi - phip)/2)^(2)*(b - (-a) + (b + (-a))*tan((phi - phip)/2)^2))/(b - (-a))));
        
        term5   = - (r*(r^2 - 2*r*rp - 3*rp^2 + 3*b - a + (2*r^(3)*rp - 2*rp^4 - 2*b^2 - 2*r*rp*(rp^2 - b - 2*(-a)) + r^(2)*(2*rp^2 - 6*b - (-a)) - b*(-a) - (-a)^2 + rp^(2)*(4*b + (-a)))/sqrt(4*r^(2)*b - 4*r*rp*(-a) + (-a)^2))*cos((phi - phip)/2)^(2)*ellipticE((phi - phip)/2, -((2*(-a))/(b - (-a))))*(b - (-a) + (b + (-a))*tan((phi - phip)/2)^2))/((r^2 + 2*r*rp + rp^2 - b - (-a))*(2*r^2 + 2*r*rp - (-a) - sqrt(4*r^(2)*b - 4*r*rp*(-a) + (-a)^2))*sqrt(cos((phi - phip)/2)^(2)*(b - (-a) + (b + (-a))*tan((phi - phip)/2)^2))*sqrt((cos((phi - phip)/2)^(2)*(b - (-a) + (b + (-a))*tan((phi - phip)/2)^2))/(b - (-a))));
        
        term6 = + (r*(r^2 - 2*r*rp - 3*rp^2 + 3*b - a - (2*r^(3)*rp - 2*rp^4 - 2*b^2 - 2*r*rp*(rp^2 - b - 2*(-a)) + r^(2)*(2*rp^2 - 6*b - (-a)) - b*(-a) - (-a)^2 + rp^(2)*(4*b + (-a)))/sqrt(4*r^(2)*b - 4*r*rp*(-a) + (-a)^2))*(2*rp^(2)*(-a) + 2*r^(2)*(b + 2*(-a)) + 2*r*rp*(b + 3*(-a)) - (b + (-a))*(3*(-a) - sqrt(4*r^(2)*b - 4*r*rp*(-a) + (-a)^2)))*cos((phi - phip)/2)^(2)*ellipticF((phi - phip)/2, -((2*(-a))/(b - (-a))))*(b - (-a) + (b + (-a))*tan((phi - phip)/2)^2))/((r^2 + 2*r*rp + rp^2 - b - (-a))*(b - (-a))*(2*r^2 + 2*r*rp - (-a) + sqrt(4*r^(2)*b - 4*r*rp*(-a) + (-a)^2))^(2)*sqrt(cos((phi - phip)/2)^(2)*(b - (-a) + (b + (-a))*tan((phi - phip)/2)^2))*sqrt((cos((phi - phip)/2)^(2)*(b - (-a) + (b + (-a))*tan((phi - phip)/2)^2))/(b - (-a))));
        
        term7 = + (r*(r^2 - 2*r*rp - 3*rp^2 + 3*b - a + (2*r^(3)*rp - 2*rp^4 - 2*b^2 - 2*r*rp*(rp^2 - b - 2*(-a)) + r^(2)*(2*rp^2 - 6*b - (-a)) - b*(-a) - (-a)^2 + rp^(2)*(4*b + (-a)))/sqrt(4*r^(2)*b - 4*r*rp*(-a) + (-a)^2))*(2*rp^(2)*(-a) + 2*r^(2)*(b + 2*(-a)) + 2*r*rp*(b + 3*(-a)) - (b + (-a))*(3*(-a) + sqrt(4*r^(2)*b - 4*r*rp*(-a) + (-a)^2)))*cos((phi - phip)/2)^(2)*ellipticF((phi - phip)/2, -((2*(-a))/(b - (-a))))*(b - (-a) + (b + (-a))*tan((phi - phip)/2)^2))/((r^2 + 2*r*rp + rp^2 - b - (-a))*(b - (-a))*(2*r^2 + 2*r*rp - (-a) - sqrt(4*r^(2)*b - 4*r*rp*(-a) + (-a)^2))^(2)*sqrt(cos((phi - phip)/2)^(2)*(b - (-a) + (b + (-a))*tan((phi - phip)/2)^2))*sqrt((cos((phi - phip)/2)^(2)*(b - (-a) + (b + (-a))*tan((phi - phip)/2)^2))/(b - (-a))));
        
        term8 = + (2*r*(r^2 + 2*r*rp + rp^2 - b + a)*(-a)*(r^2 - 2*r*rp - 3*rp^2 + 3*b + (-a) + (2*r^(3)*rp - 2*rp^4 - 2*b^2 - 2*r*rp*(rp^2 - b - 2*(-a)) + r^(2)*(2*rp^2 - 6*b - (-a)) - b*(-a) - (-a)^2 + rp^(2)*(4*b + (-a)))/sqrt(4*r^(2)*b - 4*r*rp*(-a) + (-a)^2))*cos((phi - phip)/2)^(2)*ellipticPi(1 + (r^2 + 2*r*rp + rp^2 - b - (-a))/(r^2 - rp^2 + b - sqrt(4*r^(2)*b - 4*r*rp*(-a) + (-a)^2)), (phi - phip)/2, -((2*(-a))/(b - (-a))))*(b - (-a) + (b + (-a))*tan((phi - phip)/2)^2))/((b - (-a))*(r^2 - rp^2 + b - sqrt(4*r^(2)*b - 4*r*rp*(-a) + (-a)^2))*(2*r^2 + 2*r*rp - (-a) - sqrt(4*r^(2)*b - 4*r*rp*(-a) + (-a)^2))^(2)*sqrt(cos((phi - phip)/2)^(2)*(b - (-a) + (b + (-a))*tan((phi - phip)/2)^2))*sqrt((cos((phi - phip)/2)^(2)*(b - (-a) + (b + (-a))*tan((phi - phip)/2)^2))/(b - (-a))));
        
        term9 = + (2*r*(r^2 + 2*r*rp + rp^2 - b + a)*(-a)*(r^2 - 2*r*rp - 3*rp^2 + 3*b + (-a) - (2*r^(3)*rp - 2*rp^4 - 2*b^2 - 2*r*rp*(rp^2 - b - 2*(-a)) + r^(2)*(2*rp^2 - 6*b - (-a)) - b*(-a) - (-a)^2 + rp^(2)*(4*b + (-a)))/sqrt(4*r^(2)*b - 4*r*rp*(-a) + (-a)^2))*cos((phi - phip)/2)^(2)*ellipticPi(1 + (r^2 + 2*r*rp + rp^2 - b - (-a))/(r^2 - rp^2 + b + sqrt(4*r^(2)*b - 4*r*rp*(-a) + (-a)^2)), (phi - phip)/2, -((2*(-a))/(b - (-a))))*(b - (-a) + (b + (-a))*tan((phi - phip)/2)^2))/((b - (-a))*(r^2 - rp^2 + b + sqrt(4*r^(2)*b - 4*r*rp*(-a) + (-a)^2))*(2*r^2 + 2*r*rp - (-a) + sqrt(4*r^(2)*b - 4*r*rp*(-a) + (-a)^2))^(2)*sqrt(cos((phi - phip)/2)^(2)*(b - (-a) + (b + (-a))*tan((phi - phip)/2)^2))*sqrt((cos((phi - phip)/2)^(2)*(b - (-a) + (b + (-a))*tan((phi - phip)/2)^2))/(b - (-a))));
        
        term10 = - (2*sqrt(-b + (-a))*(2*r*b - rp*(-a))*ellipticPi(-(((r^2 + 2*r*rp + rp^2 - b - (-a))*(b - (-a)))/((b + (-a))*(r^2 - rp^2 + b - sqrt(4*r^(2)*b - 4*r*rp*(-a) + (-a)^2)))), asin((sqrt(b + (-a))*tan((phi - phip)/2))/sqrt(-b + (-a))), (b - (-a))/(b + (-a)))*sqrt(1 + ((b + (-a))*tan((phi - phip)/2)^2)/(b - (-a))))/(sqrt(b + (-a))*sqrt(4*r^(2)*b - 4*r*rp*(-a) + (-a)^2)*sqrt(sec((phi - phip)/2)^2)*sqrt(cos((phi - phip)/2)^(2)*(b - (-a) + (b + (-a))*tan((phi - phip)/2)^2)));
        
        term11 = + (2*sqrt(-b + (-a))*(2*r*b - rp*(-a))*ellipticPi(-(((r^2 + 2*r*rp + rp^2 - b - (-a))*(b - (-a)))/((b + (-a))*(r^2 - rp^2 + b + sqrt(4*r^(2)*b - 4*r*rp*(-a) + (-a)^2)))), asin((sqrt(b + (-a))*tan((phi - phip)/2))/sqrt(-b + (-a))), (b - (-a))/(b + (-a)))*sqrt(1 + ((b + (-a))*tan((phi - phip)/2)^2)/(b - (-a))))/(sqrt(b + (-a))*sqrt(4*r^(2)*b - 4*r*rp*(-a) + (-a)^2)*sqrt(sec((phi - phip)/2)^2)*sqrt(cos((phi - phip)/2)^(2)*(b - (-a) + (b + (-a))*tan((phi - phip)/2)^2)));
	    
        part2_o = r*sum([term0 term1 term2 term3 term4 term5 term6 term7 term8 term9 term10 term11]);
    end

    part2 = part2_X;
end

totalint = 1/(4*pi)*M_mag*cos(theta_M)*part2;

totalint_V(3) = 1/(4*pi)*part2;

end

function A_factor = A_factor(s, rp, phip, zp, r, phi, z, a, b, c, d, e)
    A_factor = -(s*sqrt(b*(a^2+d^2)) + 2*r*a)/((z-zp)^2*sqrt(a + b));
end

function B_factor = B_factor(s, rp, phip, zp, r, phi, z, a, b, c, d, e)
    B_factor = r*sqrt(a + b)*(s*d*c + 2*((z-zp)^2 + 2*r^2)^2 + 2*r^2*(z-zp)^2)/(s*d*(s*d + 2*r^2)*(z-zp)^2);
end

function C_factor = C_factor(s, rp, phip, zp, r, phi, z, a, b, c, d, e)
    C_factor = -(r*(s*d*c - e)*(s*d*(b - a) - 2*a*rp^2 + 2*b*r^2 - 4*a*r^2 + 2*a*b))/(s*d*(s*d + 2*r^2)^2*(z-zp)^2*sqrt(a + b));
end

function D_factor = D_factor(s, rp, phip, zp, r, phi, z, a, b, c, d, e)
    D_factor = -(2*a*r*(s*d*c + e)*(z-zp)^2)/(s*d*(s*d - 2*r^2)^2*(s*d - b + rp^2 - r^2)*sqrt(a + b));
end

function E_factor = E_factor(s, rp, phip, zp, r, phi, z, a, b, c, d, e)
    E_factor = - (2*sqrt(-1)*(a*rp + 2*b*r))/(d*sqrt(b - a));
end


%%-------------------------------------------------------------------------
%%------------- The tube integral in z',phi' (constant r') ----------------
%%-------------------------------------------------------------------------
function [totalint, totalint_V] = analytical_r0(rp, phip, zp, r, phi, z, M_mag, theta_M, phi_M, slice, info)
do_terms = 1;

part2_V(1:3) = 0;

a2 = -2*r*rp./(z-zp)^2;
b2 = (r^2+rp^2+(z-zp)^2)./(z-zp)^2;

a = -2*r*rp;
b = r^2+rp^2+(z-zp)^2;

[E_int, F_int, ~, ~, ~, ~, Pi3_int] = calculate_elliptic_integrals(slice, rp, phip, zp, r, phi, z, a, b, [], [], [], 'r', info.use_matlab);

if (r == 0)
    if (slice)
        part2 = -log((z - zp) + sqrt(rp^2 + (z - zp)^2))*sin(phi_M-phip);

        factor_1 = log((z - zp) + sqrt(rp^2 + (z - zp)^2));
        part2_V(1) = +factor_1*sin(phip);
        part2_V(2) = -factor_1*cos(phip);
    else
        part2 = 0;
    end
    do_terms = 0;
end

if (abs(z-zp) < 1e-9)
    if (slice)
        part2 = sqrt(-1)*b/a*pi/2*sin(phi_M-phi);

        factor_1 = sqrt(-1)*b/a*pi/2;
        part2_V(1) = -factor_1*sin(phi);
        part2_V(2) = +factor_1*cos(phi);
    else
        %--- When doing the full cylinder, the imaginary part should change
        %--- sign and dissappear. Thus we take the real value which is zero
        % part2 = real(b/a*sqrt(-1)*pi/2*sin(phi_M-phi));
        part2 = 0;
    end
    do_terms = 0;
end

if ((abs(z-zp) < 1e-9) && r == 0)
    if (slice)
        part2 = -log(rp)*sin(phi_M - phip);

        part2_V(1) = +log(rp)*sin(phip);
        part2_V(2) = -log(rp)*cos(phip);
    else
        part2 = 0;
    end
    do_terms = 0;
end


if (do_terms == 1)
    %--- Integral with cos
    F_factor = -(z-zp)*sqrt(a + b)/a;
    term1_X = F_factor*E_int;
    
    G_factor = -(z-zp)*((z-zp)^2 - 2*b)/(a*sqrt(a + b));
    term2_X = G_factor*F_int;

    H_factor = (z-zp)*((z-zp)^2 + a - b)/(a*sqrt(a + b));
    term3_X = H_factor*Pi3_int;
    
    if info.evaluate_old_expressions
        term1 = (sqrt(b2 + a2*cos((phi - phip)))*ellipticE((phi - phip)/2, (2*a2)/(a2 + b2)))/(a2*sqrt((b2 + a2*cos((phi - phip)))/(a2 + b2)));
    
        term2 = + ((1 - 2*b2)*sqrt((b2 + a2*cos((phi - phip)))/(a2 + b2))*ellipticF((phi - phip)/2, (2*a2)/(a2 + b2)))/(a2*sqrt(b2 + a2*cos((phi - phip))));
        
        term3 = - ((1 + a2 - b2)*sqrt((b2 + a2*cos((phi - phip)))/(a2 + b2))*ellipticPi(-((2*a2)/(1 - a2 - b2)), (phi - phip)/2, (2*a2)/(a2 + b2)))/(a2*sqrt(b2 + a2*cos((phi - phip))));
    end

    if (slice)
        I_factor = - atanh((z-zp)/sqrt(b + a*cos(phi - phip)))*sin(phi - phip);	 
	         
        %--- Integral with sin
        J_factor = (b - (z-zp)^2)/a*atanh(sqrt(b + a*cos(phi - phip))/(z-zp));
	        
        K_factor = atanh((z-zp)/sqrt(b + a*cos(phi - phip)))*cos(phi - phip);
        
        L_factor = (z-zp)*sqrt(b + a*cos(phi - phip))/a;
    
        part2_X = cos(phi_M-phi).*sum([term1_X term2_X term3_X I_factor]) - sin(phi_M-phi).*sum([J_factor K_factor L_factor]);

        part2_V(1) = sum([term1_X term2_X term3_X I_factor])*cos(phi) + sum([J_factor K_factor L_factor])*sin(phi);
        part2_V(2) = sum([term1_X term2_X term3_X I_factor])*sin(phi) - sum([J_factor K_factor L_factor])*cos(phi);
        part2_V(3) = 0;

        if info.evaluate_old_expressions
            term4 = + atanh(1/sqrt(b2 + a2*cos((phi - phip))))*sin((phi - phip));
        
            term5 = ((1 - b2)*atanh(sqrt(b2 + a2*cos((phi - phip)))))/a2; 
        
            term6 = - atanh(1/sqrt(b2 + a2*cos((phi - phip))))*cos((phi - phip));
       
            term7 = - sqrt(b2 + a2*cos((phi - phip)))/a2;

            part2_o = -sign(z-zp)*(cos(phi_M-phi).*sum([term1 term2 term3 term4]) - sin(phi_M-phi).*sum([term5 term6 term7]));
        end

        % disp([term1 term2 term3 term4 term5 term6 term7])
    else
        part2_X = cos(phi_M-phi).*sum([term1_X term2_X term3_X]);
        
        part2_V(1) = sum([term1_X term2_X term3_X])*cos(phi);
        part2_V(2) = sum([term1_X term2_X term3_X])*sin(phi);
    end
    
    part2 = part2_X;
end

totalint   = -1/(4*pi)*rp*part2*M_mag.*sin(theta_M);

totalint_V = -1/(4*pi)*rp*part2_V;
end


%%-------------------------------------------------------------------------
%%------- The vertical face integral in z',r' (constant phi') -------------
%%-------------------------------------------------------------------------
function [totalint, totalint_V] = analytical_phi0(rp, phip, zp, r, phi, z, M_mag, theta_M, phi_M)
do_terms = 1;

a = -2*r*rp;
b = r^2+rp^2+(z-zp)^2;

if (r == 0)
    part2 = -((z-zp)*atanh(rp/sqrt(rp^2 + (z - zp)^2)) + rp*log(z - zp + sqrt(rp^2 + (z - zp)^2)));
    do_terms = 0;
end

if ((abs(z-zp) < 1e-9) && r == 0)
    part2 = -rp*log(rp); 
    do_terms = 0;
end

if (do_terms == 1)
    M_factor = sqrt(-1)*r*sin(phi-phip)*atanh(((z-zp)*(rp - r*cos(phi-phip)))/(r*sin(phi-phip)*sqrt(-b - a*cos(phi-phip))));
    N_factor =  (z-zp)*atanh((r*cos(phi-phip) - rp)/sqrt(b + a*cos(phi-phip)));
    O_factor = -rp*atanh((z-zp)/sqrt(b + a*cos(phi-phip)));
    P_factor =  r*cos(phi-phip)*atanh(sqrt(b + a*cos(phi-phip))/(z-zp));

    part2 = sum([M_factor N_factor O_factor P_factor]);
end

totalint = 1/(4*pi)*M_mag.*sin(theta_M).*sin(phi_M - phip)*part2;

totalint_V(1) = -1/(4*pi)*part2*sin(phip);
totalint_V(2) = +1/(4*pi)*part2*cos(phip);
totalint_V(3) = 0;

end


%%-------------------------------------------------------------------------
%%-------------- The evaluation of the elliptic integrals -----------------
%%-------------------------------------------------------------------------

function [E_int, F_int, Pi1p_int, Pi1n_int, Pi2p_int, Pi2n_int, Pi3_int] = calculate_elliptic_integrals(slice, rp, phip, zp, r, phi, z, a, b, c, d, e, surf_type, use_matlab)
    Pi1p_int = 0;
    Pi1n_int = 0;
    Pi2p_int = 0;
    Pi2n_int = 0;
    Pi3_int  = 0;

    if (slice)
        int_argument_1 = (phi - phip)/2;
        int_argument_2 = 2*a/(a + b);

        if (use_matlab)
            E_int = ellipticE(int_argument_1, int_argument_2);
            F_int = ellipticF(int_argument_1, int_argument_2);
        else
            [F_int,E_int] = elliptic123(int_argument_1,int_argument_2);
        end
        if (surf_type == 'z')
            Pi1p_int = Pi1_int(+1, rp, phip, zp, r, phi, z, a, b, c, d, e, use_matlab);
            Pi1n_int = Pi1_int(-1, rp, phip, zp, r, phi, z, a, b, c, d, e, use_matlab);
            Pi2p_int = Pi2_int(+1, rp, phip, zp, r, phi, z, a, b, c, d, e, use_matlab);
            Pi2n_int = Pi2_int(-1, rp, phip, zp, r, phi, z, a, b, c, d, e, use_matlab);
        else
            if (((z-zp)^2 - a - b) ~= 0)
                int_argument_1 = -2*a/((z-zp)^2 - a - b);
                int_argument_2 = (phi - phip)/2;
                int_argument_3 = 2*a/(a + b);

                if (use_matlab)
                    Pi3_int = ellipticPi(int_argument_1, int_argument_2, int_argument_3);
                else
                    [~,~,Pi3_int] = elliptic123(int_argument_2,int_argument_3,int_argument_1);
                end
            else
                Pi3_int = 0;
            end
        end
    else
        int_argument_1 = 2*a/(a + b);
        if (use_matlab)
            E_int = ellipticE(int_argument_1);
            F_int = ellipticK(int_argument_1);
        else
            [F_int,E_int] = elliptic123(int_argument_1);
        end
        if (surf_type == 'z')
            int_argument_1 = 1 + (z-zp)^2/(+d - 2*r^2 - (z-zp)^2);
            int_argument_2 = 1 + (z-zp)^2/(-d - 2*r^2 - (z-zp)^2);
            int_argument_3 = 2*a/(a + b);
            if (int_argument_1 == 0 && int_argument_2 == 0 && int_argument_3 == 0)
                Pi1p_int = 0;
                Pi1n_int = 0;
            elseif (isnan(int_argument_1) || isnan(int_argument_2) || isnan(int_argument_3))
                Pi1p_int = 0;
                Pi1n_int = 0;
            else
                if (use_matlab)
                    Pi1p_int = ellipticPi(int_argument_1, int_argument_3);
                    Pi1n_int = ellipticPi(int_argument_2, int_argument_3);
                else
                    [~,~,Pi1p_int] = elliptic123(int_argument_3,int_argument_1);
                    [~,~,Pi1n_int] = elliptic123(int_argument_3,int_argument_2);
                end 
            end
            Pi2p_int = Pi2_int_full(+1, rp, phip, zp, r, phi, z, a, b, c, d, e, use_matlab);
            Pi2n_int = Pi2_int_full(-1, rp, phip, zp, r, phi, z, a, b, c, d, e, use_matlab);
        else
            if (((z-zp)^2 - a - b) ~= 0)
                int_argument_1 = -2*a/((z-zp)^2 - a - b);
                int_argument_2 = 2*a/(a + b);
                if (use_matlab)
                    Pi3_int = ellipticPi(int_argument_1, int_argument_2);
                else
                    [~,~,Pi3_int] = elliptic123(int_argument_2, int_argument_1);
                end
            else
                Pi3_int = 0;
            end
        end
    end
end

%%-------------------------------------------------------------------------
%%--------------------------- Evaluation of P1 ----------------------------
%%-------------------------------------------------------------------------
function Pi1_int = Pi1_int(s, rp, phip, zp, r, phi, z, a, b, c, d, e, use_matlab)
    int_argument_1 = 1 + (z-zp)^2/(s*d - 2*r^2 - (z-zp)^2);
    int_argument_2 = (phi - phip)/2;
    int_argument_3 = 2*a/(a + b);

    if (isnan(int_argument_1))
        Pi1_int = 0;
        return
    end

    if (int_argument_1 == 1)
        % This in fact is not relevant as there is a seperate considertion
        % for z-zp = 0 in the face-integral, and z-zp = 0 is the only
        % condition that can cause this.
        Pi1_int = (sqrt(1-int_argument_3*sin(int_argument_2)^2)*tan(int_argument_2)-ellipticE(int_argument_2,int_argument_3))/(1-int_argument_3)+ellipticF(int_argument_2,int_argument_3);
        return
    end

    if (use_matlab)
        Pi1_int = ellipticPi(int_argument_1, int_argument_2, int_argument_3);
    else
        [~,~,Pi1_int] = elliptic123(int_argument_2,int_argument_3,int_argument_1);
    end
end

%%-------------------------------------------------------------------------
%%--------------------------- Evaluation of P2 ----------------------------
%%-------------------------------------------------------------------------
function Pi2_int = Pi2_int(s, rp, phip, zp, r, phi, z, a, b, c, d, e, use_matlab)
    int_argument_1 = (a + b)*(z-zp)^2/((a - b)*(s*d - 2*r^2 - (z-zp)^2));
    int_argument_2 = -asin(sqrt((a - b)/(a + b))*tan((phi - phip)/2));
    int_argument_3 = -(a + b)/(a - b);

    if (isnan(int_argument_1))
        Pi2_int = 0;
        return
    end

    if (use_matlab)
        Pi2_int  = ellipticPi(int_argument_1, int_argument_2, int_argument_3);
    else
        %--- As the first argument is complex, the numerical implementation
        %--- in elliptic123 does not work and we have to use Matlab's
        %--- build-in expression to evaluate the integral
        %-------
        %--- Using https://dlmf.nist.gov/19.7#E5 the argument can be
        %--- changed from imaginary to non-imaginary. However, the expression
        %--- given in https://dlmf.nist.gov/19.7#E5 is not correct!
        %--- Instead, using the formula in
        %--- https://personal.math.ubc.ca/~cbm/aands/abramowitz_and_stegun.pdf
        %--- Eq. 17.4.8, we can see that we need to modify the expression such
        %--- that the last argument is changed. We get
        psi = atan(sinh(imag(int_argument_2))); 
        alpha = asin(sqrt(int_argument_3));
        Pi2_int = sqrt(-1)*((ellipticF(psi, sin(1/2*pi-alpha)^2)) - int_argument_1*ellipticPi(1-int_argument_1, psi, sin(1/2*pi-alpha)^2))/(1-int_argument_1);
    end

    %--- Note: The Pi2p_int and Pi2n_int functions are not continous across
    % the boundary phip-phi > pi because of the asin(tan()) factor. 
    % These function cannot be modified to be this.
    % The problem can be fixed by instead first integrating up to phip-phi = pi 
    % and then from phi-pi to the limit, which will be at 2*(phi-pi)+(angle_arr(i)-pi). 
    % In that way we avoid the singularity.
    % Thus we can define the following angle if (phi1-phi < -pi) then
    % phi1_new = ((phi1-phi)-pi)+(phi-pi) and phi2_new = phi2
    % and if (phi2-phi > pi) then phi1_new=phi1 and phi2_new = ((phi2-phi)-pi)+(phi-pi); 
    % and then integrate as follows
    %     q_z1_ana_part_1 = -((analytical_z0(Ro, phi+pi  , -Height/2, r, phi, z, M_mag, theta_M, phi_M, slice) - analytical_z0(Ro, phi1_new, -Height/2, r, phi, z, M_mag, theta_M, phi_M, slice)) ...
    %                       - (analytical_z0(Ri, phi+pi  , -Height/2, r, phi, z, M_mag, theta_M, phi_M, slice) - analytical_z0(Ri, phi1_new, -Height/2, r, phi, z, M_mag, theta_M, phi_M, slice)));
    % 
    %     q_z1_ana_part_2 = -((analytical_z0(Ro, phi2_new, -Height/2, r, phi, z, M_mag, theta_M, phi_M, slice) - analytical_z0(Ro, phi-pi  , -Height/2, r, phi, z, M_mag, theta_M, phi_M, slice)) ...
    %                       - (analytical_z0(Ri, phi2_new, -Height/2, r, phi, z, M_mag, theta_M, phi_M, slice) - analytical_z0(Ri, phi-pi  , -Height/2, r, phi, z, M_mag, theta_M, phi_M, slice)));
    % 
    %     q_z1_ana = q_z1_ana_part_1 + q_z1_ana_part_2;
    % However, this corresponds to evaluating analytical_z0 at exactly
    % phi+pi, which is what we do for the full cylinder integral anyways!
    % As it is only the Pi2p_int and Pi2n_int inegrals that cause problems,
    % we can just consider the normal phi1 and phi2 integral limits, if we
    % remember to add in the additional factor that results from the split
    % integral above evaluated at phi+pi. But that additional factor is
    % then only the Pi2p_int and Pi2n_int terms as the other terms are included in the phi1 and phi2 limits!
    % As we need both the phi+pi and phi-pi, we must add Pi2p_int and Pi2n_int
    % twice. And they must be substrated because of the minus sign above.
    if (phip-phi < -pi) || (phip-phi > pi)
        Pi2_int = Pi2_int - 2*Pi2_int_full(s, rp, phip, zp, r, phi, z, a, b, c, d, e, use_matlab);
    end  
end

function Pi2_int_full = Pi2_int_full(s, rp, phip, zp, r, phi, z, a, b, c, d, e, use_matlab)
    n = (a + b)/(a - b)*(z-zp)^2/(s*d - 2*r^2 - (z-zp)^2);
    int_argument_1 = 2*a/(a - b);

    if (n > (1-1e-9)) && (n < (1+1e-9))
        if (use_matlab)
            Pi2_int_full = sqrt(-1)*(-ellipticK(int_argument_1) + (a - b)/(2*a)*(ellipticK(int_argument_1) - ellipticE(int_argument_1)));
        else
            [ellipticK_int,ellipticE_int] = elliptic123(int_argument_1);
            Pi2_int_full = sqrt(-1)*(-ellipticK_int + (a - b)/(2*a)*(ellipticK_int - ellipticE_int));
        end
    else
        if (isnan(int_argument_1) || isnan(n))
            Pi2_int_full = 0;
            return
        end
        if (use_matlab)
            Pi2_int_full = sqrt(-1)/(1-n)*(n*ellipticPi(1-n,int_argument_1)-ellipticK(int_argument_1));
        else
            [~,~,ellipticPi_int] = elliptic123(int_argument_1,1-n);
            [ellipticK_int,~] = elliptic123(int_argument_1);
            Pi2_int_full = sqrt(-1)/(1-n)*(n*ellipticPi_int-ellipticK_int);
         end
    end
end

%% Notes
%--- Some special values for the elliptic integrals can be found from here:
% https://functions.wolfram.com/EllipticIntegrals/EllipticPi/03/01/02/
% https://functions.wolfram.com/EllipticIntegrals/EllipticPi3/03/01/03/
% https://functions.wolfram.com/HypergeometricFunctions/Hypergeometric2F1/03/
% https://reference.wolfram.com/language/ref/EllipticE.html
% https://reference.wolfram.com/language/ref/EllipticF.html
% https://reference.wolfram.com/language/ref/EllipticPi.html

%--- The numerical fast evaluation of the elliptic integrals is taken from
% https://github.com/moiseevigor/elliptic?tab=readme-ov-file#elliptic3-incomplete-elliptic-integral-of-the-third-kind