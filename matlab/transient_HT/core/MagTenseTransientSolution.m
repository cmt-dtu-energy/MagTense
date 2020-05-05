

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
    
    methods
        function obj = MagTenseTransientSolution(  )
           
        end
    end
    
end