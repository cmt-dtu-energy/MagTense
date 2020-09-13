function [J,g,Ns,muB,mu0,kB,rhos,Na,mol,Tdeb,gamma,kappa,Na_c, T0, beta] = getMFTParams( material )
%default values

%Bohr magneton
muB = 9.274e-24;

%Avogadro's number
Na = 6.02e23;

%defaul molar mass of the material (kg/mol)
mol = 0.157;
%Boltzmann's constant
kB = 1.38e-23;

%vacuum permeability
mu0 = 4.*pi*1e-7;
%Default angular momentum
J = 3.5;
%Default Landé factor
g = 2.;
%default spin density
Ns = 3.83e24;
%default mass density
rhos = 6000.;
%default Debye temperature
Tdeb = 169.;
%default Sommerfeld constant
gamma = 6.93e-2;
%default compressibility
kappa = 7.3e-3 *1e-9; %1/Pa

%Atoms per unit cell
Na_c = 1;

switch material
    case'Gd'
        J = 3.5;
        g = 2.;  
        Ns = 3.83e24;
        rhos = 7900.;
        mol = 0.157;
        Tdeb = 169.;
        gamma = 6.93e-2;
        Na_c = 1;
        T0 = 293;
    case'LCM'
        nMn = 1;
        %La_0.67Ca0.23MnO3
        J = 1.83;
        g = 2;
        kappa = 1e-13;%Pa^-1
        rhos = 6034;%kg/m^3
        mol = 0.209224;%kg/mol
        Ns = 1/mol * Na * nMn;
        Tdeb = 353;%K
        gamma = 0.025;
        Na_c = 5;
        T0 = 266.6;
	case'LCSM_1'
        nMn = 1;
        %La_0.67Ca0.23MnO3
        J = 1.83;
        g = 2;
        kappa = 1e-13;%Pa^-1
        rhos = 6034;%kg/m^3
        mol = 0.209224;%kg/mol
        Ns = 1/mol * Na * nMn;
        Tdeb = 353;%K
        gamma = 0.025;
        Na_c = 5;
        T0 = 278.6;
	case'LCSM_2'
        nMn = 1;
        %La_0.67Ca0.23MnO3
        J = 1.83;
        g = 2;
        kappa = 1e-13;%Pa^-1
        rhos = 6034;%kg/m^3
        mol = 0.209224;%kg/mol
        Ns = 1/mol * Na * nMn;
        Tdeb = 353;%K
        gamma = 0.025;
        Na_c = 5;
        T0 = 282.6;
	case'LCSM_3'
        nMn = 1;
        %La_0.67Ca0.23MnO3
        J = 1.83;
        g = 2;
        kappa = 1e-13;%Pa^-1
        rhos = 6034;%kg/m^3
        mol = 0.209224;%kg/mol
        Ns = 1/mol * Na * nMn;
        Tdeb = 353;%K
        gamma = 0.025;
        Na_c = 5;
        T0 = 286.5;
	case'LCSM_4'
        nMn = 1;
        %La_0.67Ca0.23MnO3
        J = 1.83;
        g = 2;
        kappa = 1e-13;%Pa^-1
        rhos = 6034;%kg/m^3
        mol = 0.209224;%kg/mol
        Ns = 1/mol * Na * nMn;
        Tdeb = 353;%K
        gamma = 0.025;
        Na_c = 5;
        T0 = 289.4;
	case'LCSM_5'
        nMn = 1;
        %La_0.67Ca0.23MnO3
        J = 1.83;
        g = 2;
        kappa = 1e-13;%Pa^-1
        rhos = 6034;%kg/m^3
        mol = 0.209224;%kg/mol
        Ns = 1/mol * Na * nMn;
        Tdeb = 353;%K
        gamma = 0.025;
        Na_c = 5;
        T0 = 292.5;
	case'LCSM_6'
        nMn = 1;
        %La_0.67Ca0.23MnO3
        J = 1.83;
        g = 2;
        kappa = 1e-13;%Pa^-1
        rhos = 6034;%kg/m^3
        mol = 0.209224;%kg/mol
        Ns = 1/mol * Na * nMn;
        Tdeb = 353;%K
        gamma = 0.025;
        Na_c = 5;
        T0 = 298.6;
end

if ( strcmp(material,'LaFeSi') )
    %La(Fe_0.88Si_0.12)_13H_y
    J = 1;
    g = 2;
    kappa = 7.3e-3 *1e-9;
    
    x = 0.88;
    nFe = 13;
    
    rhos = 7200;
    mol = 0.821;
    
    Ns = 1 / mol * Na * nFe * x;
    
    Tdeb = 350;
    gamma = 0.235;
    Na_c = 14;
    T0 = 300;
end
if ( strcmp(material,'Gd_like_1st_order_Anders'))
    J = 3.5;
    g = 2;
    kappa = 1/(37.91e9);
    rhos = 7900;
    Ns = 3.83e28/rhos;
    Tdeb = 169;
    gamma = 0.0693;
    T0 = 293;
    beta = 1;
end
if ( strcmp(material,'LaFeSi_dehyd') )
    %VAC70 sample that has been dehydrogenated and
    %has had M and XRD as a function of field measured
    %as well as M as a function of applied pressure
   J = 0.92;
   g = 2;
   kappa = 7.5e-11;%Pa^-1
   
   rhos = 7200;
   Ns = 8.45e24;% kg^-1
   Tdeb = 350;%K
   gamma = 0.235;
   T0 = 151;
   beta = 14.1;
end

if ( strcmp(material,'Gd5Si2Ge2') )
    J = 3.5;%????
    
    gamma = 0.025;
    
	
	%Gd5Si2Ge2
	nGd = 5;
    g = 2;    
    mol = 0.988;    
    Na_c = 9;
	kappa = 1e-13;%Pa^-1
	Ns = 1/mol * Na * nGd;
	rhos = 5700;%kg/m^3   
    Tdeb = 250;%K, Svitelskiy, PRB 74, 184105, 2006
    
end

if ( strcmp(material,'Fe') ) 
   %J = 3.5;
   %g = 2.;  
   %Ns = 3.83e24;
   rhos = 7100.;
   %mol = 0.157;
   %Tdeb = 169.;
   %gamma = 6.93e-2;
end



end