

classdef MagTenseTransientSettings
   
    properties
        maxIte %maximum no. of iterations
        maxErr %maximul relative err for iteration step
        dt %current length of the timestep in seconds
        t_tot %the total time / physical duration of the simulation
        t %current time (s)
        %thermal property functions
        c
        k
        rho
        %magnetization function
        M
        %applied field function
        Happ
        %state function function
        stateFunction MagTenseStateFunction
    end
    
    methods
        
        
        %constructor takes in function handles for the thermal properties
        function obj = MagTenseTransientSettings( c_in, k_in, rho_in, M_in, Happ_in, stFct )
           %default values, can be changed directly on the object when
           %running
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
           if exist('M_in','var' ) && isa( M_in, 'function_handle' )
            obj.M = M_in;
           end
           if exist('Happ_in','var' ) && isa( Happ_in, 'function_handle' )
            obj.Happ = Happ_in;
           end
           if exist('stFct','var' ) && isa( stFct, 'MagTenseStateFunction' )
            obj.stateFunction = stFct;
           end
            
        end
        
        %finds the hysteretic state based on H, T and p and the state
        %function given in each tile
        function hyst = Hyst( obj, solution )
            %if clearly in FM then set hyst = 0 (constant in
            %MagTenseStateFunction)
            %if clearly in PM then set hyst = 1 (constant in
            %MagTenseStateFunction)
            %if not clearly in either state then don't change hyst as the
            %state is then given by history
            Hnorm = sqrt(sum( solution.H.^2, 2 ) );
            
            %initialize hysteresis to current state
            hyst = solution.hyst;
            
            [Tcrit_heat, Tcrit_cool] = obj.stateFunction.getTcrit( Hnorm );
            
            %set FM
            hyst( solution.T < Tcrit_cool ) = MagTenseStateFunction.FM;
            %set PM
            hyst( solution.T > Tcrit_heat ) = MagTenseStateFunction.PM;
            
        end
        
    end
    
end