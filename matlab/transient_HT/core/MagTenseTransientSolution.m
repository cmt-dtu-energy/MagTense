

classdef MagTenseTransientSolution
   
    properties 
        T %single column vector of length n containing the current temperature
        H %(n,3) matrix containing column vectors for the Hx, Hy and Hz components of the magnetic field
        p %single column vector of length n containing the current pressure.
        
        %thermal properties
        c %single column length n vector containing the current specific heat
        k %single column length n vector containing the current thermal conductivity
        rho %single column length n vector containing the current mass density
        
        
    end
    
    properties (SetAccess=private,GetAccess=public)
       nSnap %the current snapshot ID
       T_sol %cell array with the solution at given points during the temporal evolution
       t_sol %double array with the time at which the solution was saved (indices correspond to T_sol)
    end
    
    methods
        function obj = MagTenseTransientSolution(  )
            obj.t_sol = 0;
            obj.nSnap = 1;
        end
        
        %called when the solution is to be stored at this given time, t
        function obj = StoreSolution( obj, t )
            obj.t_sol(obj.nSnap) = t;
            obj.T_sol{obj.nSnap} = obj.T;
            obj.nSnap = obj.nSnap + 1;
        end
        
    end
    
end