

%%This functions shows how to use MagTense to find the field from a finite
%%sized magnet and calculate the force on a particle with a given
%%magnetization assumed to be parallel to the local field and with a given
%%magnitude. The equations of motions are then integrated as a function of
%%time including Stokes drag on the particle assuming it is suspended in a
%%viscous fluid.
function [] = MagTense_Particle()
close all
addpath('../../util/');
addpath('../../Mex_files/');

%Define magnet geometry using MagTense. The magnet considered is a
%cube with side lengths 10 mm

%side length of the magnet
d = 0.01;%m

%remanence of the magnet
Mrem = 1.2;%R

%vacuum permeability
mu0 = 4*pi*1e-7;

%get the default tile
tile1 = getDefaultMagTile();
%set the tile type to cylinder
tile1.tileType = getMagTileType('prism');
%set the remanence 
tile1.Mrem = Mrem/mu0;%A/m
%the easy axis permeability
tile1.mu_r_ea = 1.05;
%the hard axes permeabilities
tile1.mu_r_oa = 1.18;
%ensure the tile is a permanent magnet
tile1.magnetType = getMagnetType('hard');
%the size of the tile
tile1.abc = [d,d,d];
%set the easy axis and the hard axes
tile1.u_ea = [1,0,0];
tile1.u_oa1 = [0,1,0];
tile1.u_oa2 = [0,0,1];


%boundle the tiles into an array
tiles = [tile1];

%find the self-consistent magnetization (in this case only caused by the
%demagnetizing field)
tiles = IterateMagnetization( tiles, [], [], 1e-6, 100 );

%%Now find the field in a set of points
%define a range of points spanning the xy plane at z=0
x = linspace( d/2+1e-4, d, 51 );
y = linspace( 0, d, 50 );
z = 0.0001;


%use meshgrid to fill out the span
[X,Y,Z] = meshgrid(x,y,z);


%get the field
[H,Hnorm] = getHMagTense( tiles, X, Y, Z );
H = H .* mu0;
Hnorm = Hnorm .* mu0;

%plot the field
surf_and_con(squeeze(X),squeeze(Y),Hnorm);

%now get on with integration of the equations of motion:
%ms * d^2x/dt^2 = F_mag_x - C * R * |dx/dt| * (dx/dt / |dx/dt| )
%ms * d^2y/dt^2 = F_mag_y - C * R * |dy/dt| * (dy/dt / |dy/dt| )
%with the mass m, F_max_x (y) the magnetic force in the x (y) direction, C
%the coefficient for Stokes drag on a sphere (=6*pi*mu with mu being the
%viscosity of the fluid) and the drag force defined such that it points
%away from the direction of motion. The magnetic force can be found on a
%point particle by F_mag = grad( m dot B ) where m is the magnetization of
%the particle (m=M * V with V = volume of the particle and M = the
%magnetization in A/m) and B the local magnetic flux density in T. Finally,
%ms is the mass of the particle
%We now get

%d^2x/dt^2 = grad( m dot B )|_x / ms - 6*pi*mu*R/ms * dx/dt
%d^2y/dt^2 = grad( m dot B )|_y / ms - 6*pi*mu*R/ms * dy/dt

%define constants
mu = 0.001;%Pa * s (water)
%particle diameter
R = 2.5e-6;%m
%particle mass

%Volume of magnetite
V_Fe = 7.24e-19;%m^3
%magnetite density
rho_Fe = 5200;%kg/m^3
%polymer volume
V_p = 4/3*pi*R^3 - V_Fe;
%polymer density
rho_p = 1000;%kg/m^3
%mass of the particle
ms = V_Fe * rho_Fe + V_p * rho_p;
%magnetization for the sphere
m = 2.83e-13;%Am^2
%Drag coefficient
C = 6 * pi * mu;


%bundle into struct
dat = struct( 'ms', ms, 'C', C, 'R', R, 'm', m, 'd', d );

%note that y_init(1) and y_init(3) are the initial positions in x,y while
%y_init(2) and y_init(4) are initial velocities
y_init = [d/2 + 0.002, 0, 0.001, 0];

%Make a grid in x and y (really z and y) for Hnorm and use this to 
%interpolate in       
nx = 100;
ny = 101;

x_tmp = linspace( d/2, y_init(1) * 1.1, nx );
y_tmp = linspace( 0.00001, y_init(3)*1.1, ny );
[dat.X,dat.Y] = meshgrid(x_tmp,y_tmp);
pts = zeros(numel(dat.Y),3);
pts(:,1) = dat.X(:);
pts(:,2) = dat.Y(:);
pts(:,3) = 0;
H = getHFromTiles_mex( tiles, pts, int32( length(tiles) ), int32( length(pts(:,1)) ) );
H = H .* mu0;
dat.Hnorm = sqrt(sum(H.^2,2));%note that H and B are identical outside the magnet
dat.Hnorm = reshape( dat.Hnorm, [ny,nx] );    

%anonymous function to pass arguments
fct = @(t,y) ode_func( t, y, dat );
%time grid for the solution
t_in = linspace(0,10,300);


opts = odeset('OutputFcn',@plotODE,'RelTol',1e-3);
figure
subplot(3,1,1);
xlabel('x pos [mm]');
ylabel('y pos [mm]');
hold on
subplot(3,1,2);
xlabel('Time [s]');
ylabel('x velocity [mm/s]');
hold on
subplot(3,1,3);
xlabel('Time [s]');
ylabel('y velocity [mm/s]');
hold on

[t,y] = ode45( fct, t_in, y_init, opts );

end

function status = plotODE( t, y, flag )
    if isempty(flag)
        cols = lines(2);
        subplot(3,1,1);
        plot(y(1)*1e3,y(3)*1e3,'.','color',cols(1,:));
        
        
        subplot(3,1,2);        
        plot(t(1),y(2)*1e3,'.','color',cols(1,:));
        subplot(3,1,3);        
        plot(t(1),y(4)*1e3,'o','color',cols(2,:));
        
        drawnow
    end
    status = 0;
end

%returns the derivatives of the system of ODE's to be solved
%t and y are the standard inputs from the ode solver (time and current y
%value) while dat is a struct with the extra stuff needed for the equations
%The four elements in y denote the following:
%y(1,t) is x(t)
%y(2,t) is gamma(t)
%y(3,t) is y(t)
%y(4,t) is beta(t)
%gamma and beta are included to transform the second order ODE system to a
%first order system such that gamma = dx/dt and beta = dy/dt
function [yp] = ode_func(t,y,dat)

    ms = dat.ms;
    C = dat.C;
    R = dat.R;
    m = dat.m;
    
    
    %get the magnetic flux through numerical differentiation of 
    %Fmag = grad( m dot B ) assuming m = |m| B/|B| =>
    %Fmag = grad( |m| * |B| ) = |m| * grad( |B| ) as |m| is assumed
    %constant
    %distance to differentiate over
    delta = 1e-6;%m
    xx = [y(1)-delta, y(1)+delta];
    %y-direction
    yy = [y(3)-delta, y(3)+delta];
    %x-direction (modeled as z due to the nature of the cylindrical tile
    %and the laziness of the programmer not to include rotation)
    zz = 0;
       
    
    Hnorm_x = interp2( dat.X, dat.Y, dat.Hnorm, xx, [y(3), y(3)] );
    Hnorm_y = interp2( dat.X, dat.Y, dat.Hnorm, [y(1), y(1)], yy );
   
    Fmag = zeros(2,1);    
    Fmag(1) = m * ( Hnorm_x(2) - Hnorm_x(1) ) / ( xx(2) - xx(1) );
    Fmag(2) = m * ( Hnorm_y(2) - Hnorm_y(1) ) / ( yy(2) - yy(1) );
    
    if ~isfinite(Fmag) | ~isfinite(Hnorm_x) | ~isfinite(Hnorm_y)
        disp('NAN!!!');
    end
    %re-write d^2x/dt^2 = F/ms - C/ms*R*dx/dt to two first order equations
    %dx/dt = gamma (yp(1))
    %dgamma/dt = d^2x/dt^2 (yp(2))
    %and similarly for the y-component
    %dy/dt = beta (yp(3))
    %dbeta/dt = d^2y/dt^2 (yp(4))   
    
    yp = zeros(4,1);
    %x-direction
    %check if the particle has reached the boundary
    if y(1) <= min(dat.X(:))
        yp(1) = 0;
        yp(2) = 0;
    else
        yp(1) = y(2);
        yp(2) = Fmag(1) / ms - C/ms * R * y(2); %note that dx/dt = gamma = y(2)
    end
    %y-direction
    if y(3) <= min(dat.Y(:))
        yp(3) = 0;
        yp(4) = 0;
    else                
        yp(3) = y(4);
        yp(4) = Fmag(2) / ms - C/ms * R * y(4);%note that dy/dt = beta = y(4)
    end
end
