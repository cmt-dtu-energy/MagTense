

classdef MagTenseTransientSettings
   
    properties
        maxIte %maximum no. of iterations
        maxErr %maximul relative err for iteration step
        dt %current length of the timestep in seconds
        t_tot %the total time / physical duration of the simulation
        t %current time (s)
        %thermal property functions
        c;
        k;
        rho;
    end
    
    methods
        
        
        %constructor takes in function handles for the thermal properties
        function obj = MagTenseTransientSettings( c_in, k_in, rho_in )
           obj.maxIte = 10;
           obj.maxErr = 0.001;
           obj.dt = 0.001;%seconds
           
           if isa( c_in, 'function_handle' )
            obj.c = c_in;
           end
           if isa( k_in, 'function_handle' )
            obj.k = k_in;
           end
           if isa( rho_in, 'function_handle' )
            obj.rho = rho_in;
           end
        end
    end
    
end