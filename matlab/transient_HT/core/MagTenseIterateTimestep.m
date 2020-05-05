
%By Kaspar K. Nielsen, kasparkn@gmail.com
%Latest update 20200502
%This function progresses the temporal solution of the coupled heat
%transfer and magnetic field problem iteratively, which means that it both
%solves for the temperature and magnetic field (T and H) at the new
%timestep (t+dt) but also includes the various properties (conductivity,
%density, specific heat, entropy etc.) as averages over the timestep, i.e.
%as weighted functions of the current (t) and the future (t+dt) timesteps
%The idea is that given dt and t the iteration works like this:
%while err>maxErr && nIte<maxIte
%c = 0.5 * ( c(T(t),H(t)) + c(T(t+dt),H(t+dt) )
%etc
%then set up the equations and find T* = T(t+t) and check the relative err
%err = max( abs( ( T(t+dt) - T(t) ) / (T(t)) )

%end
%argument solution is an object of type MagTenseTransietSolution that will
%also contain the resulting temperature field (new solution) and is thus
%modified and returned.
%argument geom is an object of type MagTenseTransientGeometry
%argument setts is an object of type MagTenseTransientSettings
function [solution, outFlag] = MagTenseIterateTimestep( solution, geom, setts, debug )

%iteration settings
maxErr = setts.maxErr;
maxIte = setts.maxIte;
dt = setts.dt;

%the geometrical part of the thermal resistance.
R_geom = geom.R_geom;
dV = geom.dV;

%no. of cells or finite volumes
n = length(dV);

%the solution at time t
T_old = solution.T;
c_old = solution.c;
k_old = solution.k;
rho_old = solution.rho;

%initialization of parameters
err = 10 * maxErr;
nIte = 0;

%the temperature at the new tims, t+dt, is denoted T_new
%the temperature at the current time, t, is denoted T_old
%initially the new temperature is set to the initial temperature
c_new = c_old;
k_new = k_old;
rho_new = rho_old;


%while loop that exits after a max. no. of iteration or when the solution
%has converged
while err>maxErr && nIte<maxIte
   
    %find resulting thermal properties (average over the timestep)
    %these all have to be single-column vectors (length n)
    c = 0.5 * ( c_old + c_new );
    k = 0.5 * ( k_old + k_new );
    rho = 0.5 * ( rho_old + rho_new );
    
    %thermal resistance matrix (each element R_ij = dl_i/(A_ij * k_i ) for
    %the i'th tile and it's surface area, A_ij, with the j'th cell)
    R = R_geom ./ repmat( k, 1, n );
    
    %setup matrix for inversion (A in A*T = b )
    %first find all thermal resistance terms that apply to all but the
    %diagonal. Note the minus sign as this is the part of the expression
    %multiplied onto the off-diagonal tile
    A = -1.*( 1./R + 1./R');
    A(~isfinite(A)) = 0.;
    %find the diagonal terms as minus the sum of the other terms
    A(1:n+1:end) = -sum(A,2);
    %find the capacity term
    cap = dV .* rho .* c ./ dt;
    
    A(1:n+1:end) = A(1:n+1:end) + cap';
    
    %setup result-vector b (right-hand side in A*T = b)
    b = cap .* T_old;
    %add the boundary conditions here
    [lhs_bdry,rhs_bdry] = geom.getBoundaryConditions(  k );
    
    A(1:n+1:end) = A(1:n+1:end) + lhs_bdry';
    b = b + rhs_bdry;
    
    %solve the equation system using \
    T_new = A \ b;
   
    %find the magnetic field at t+dt
    H_new = zeros(size(T_new));
    
    %find the pressure at t+dt
    p_new = zeros(size(T_new));
    
    %find thermal properties at the new temperature
    c_new = setts.c(T_new,H_new,p_new);
    k_new = setts.k(T_new,H_new,p_new);
    rho_new = setts.rho(T_new,H_new,p_new);
    
    %find the relative error
    err = max ( abs ( (T_new - T_old)./T_old ) );
    
    if debug
        disp( ['err = ' num2str(err,'%7.5f')] );
        disp( ['Iteration no. ' num2str(nIte)] );
    end
    
    nIte = nIte + 1;
end

%set output flag
if err<=maxErr
    outFlag = 1;
    %update the solution
    solution.T = T_new;
    solution.H = H_new;
    solution.p = p_new;
else
    outFlag = -1;
end

end