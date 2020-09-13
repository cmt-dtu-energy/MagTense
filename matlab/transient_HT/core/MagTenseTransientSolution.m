

classdef MagTenseTransientSolution
   
    properties 
        T %single column vector of length n containing the current temperature
        H %(n,3) matrix containing column vectors for the Hx, Hy and Hz components of the magnetic field
        p %single column vector of length n containing the current pressure.
        M %(n,3) matrix containing the column vectors for the Mx, My and Mz components of the magnetization
        hyst %(n,1) array with either FM or PM indicated by a constant as defined in MagTenseStateFunction.FM and .PM
        
        
        
        Happ_old %(n,3) applied field at time t
        Happ_new %(n,3) applied field at time t+dt
        
        %thermal properties
        c %single column length n vector containing the current specific heat
        k %single column length n vector containing the current thermal conductivity
        rho %single column length n vector containing the current mass density
        
        
    end
    
    properties (SetAccess=private,GetAccess=public)
       nSnap %the current snapshot ID
       T_sol %cell array with the solution at given points during the temporal evolution
       t_sol %double array with the time at which the solution was saved (indices correspond to T_sol)
       H_sol %(n,3) array containing the field
       M_sol %(n,3) array containing the magnetization
       hyst_sol % (n,1) array containing the hysteresis at the current solution
       k_sol
       c_sol
       rho_sol
       p_sol
    end
    
    methods
        function obj = MagTenseTransientSolution( n )
            obj.t_sol = 0;
            obj.nSnap = 1;
            hyst = zeros(n,1);
        end
        
        %called when the solution is to be stored at this given time, t
        function obj = StoreSolution( obj, t )
            obj.t_sol(obj.nSnap) = t;
            obj.T_sol{obj.nSnap} = obj.T;
            obj.H_sol{obj.nSnap} = obj.H;
            obj.M_sol{obj.nSnap} = obj.M;
            obj.hyst_sol{obj.nSnap} = obj.hyst;
            obj.c_sol{obj.nSnap} = obj.c;
            obj.k_sol{obj.nSnap} = obj.k;
            obj.rho_sol{obj.nSnap} = obj.rho;
            obj.p_sol{obj.nSnap} = obj.p;
            obj.nSnap = obj.nSnap + 1;
        end
        
    end
    
end