
%Kaspar Kirstein Nielsen, kasparkn@gmail.com, June 2020
%This class defines the class that abstracts the state function for the
%generalized mean field model by Bean & Rodbell (1962)
classdef MagTenseStateFunction
    
    properties
        %The following are the independent variables for the state
        %functions. Each element is a cell array of size (n,2) where n is
        %the no. of state functions and e.g. H0{i,1} represents the field
        %values for the i'th state function in ferromagnetic mode while
        %H0{i,2} is for the paramagnetic mode.
        H0 %field values
        T0 %Temperature
        p0 %pressure
        M %magnetization (note the format for each element {i,j}, i=1..n and j=1,2, is (nH,nT,np) with
        %nH = numel(H0{i,j}) etc.
        S %Entropy, same format as for M
        
        Tcrit %cell array with critical temperature for a given state function and heating or cooling
        Hcrit %similiar to Tcrit but for the field
        
        stateFctInd %array containing the indices of the state functions. Size (n,1) where n is the no. of tiles
    end
    
    properties(Constant)
       FM = 0
       PM = 1
    end
    
    methods(Static)

        %By Kaspar Kirstein Nielsen, kasparkn@gmail.com, 12 August 2020
        %Finds the critical temperature at each value of the field. If no
        %transition is found the critical temperature is set to NaN
        %Takes vector T of size (n,1) and vector H of size (m,1) and array M of
        %size (n,m)
        function [Tcrit,Hcrit] = getCriticalPoint( T, H, M )

            %critical field
            Hcrit = H;

            Tcrit = zeros( length(Hcrit), 1 );
            Tcrit(:) = NaN;

            for i=1:length(Hcrit)
                %find the second derivative
                dm2d2T = 1./( 0.5*(T(3:end)+T(2:end-1)) - 0.5*(T(2:end-1)+T(1:end-2)) ) .* ( M(3:end, i) - M(2:end-1,i) ) ./ ( T(3:end) - T(2:end-1) ) - ( M(2:end-1, i) - M(1:end-2,i) ) ./ ( T(2:end-1) - T(1:end-2) );
                
                %add one since the temperatures corresponding to dm2d2T are
                %T(2:end-1)
                Tcrit(i) = T(1 + find(dm2d2T==min(dm2d2T),1));
                
            end

        end 
    end
    
    methods
        
        %method for returning critical temperature at given field for
        %heating and cooling
        function [Tcrit_heat, Tcrit_cool] = getTcrit( obj, Hnorm )
            Tcrit_heat = zeros(length(Hnorm),1);
            Tcrit_cool = Tcrit_heat;
            for i=1:length(Hnorm)
               Tcrit_heat(i) = interp1( obj.Hcrit{obj.stateFctInd(i),1}(:,1), obj.Tcrit{obj.stateFctInd(i),1}(:,1), Hnorm(i) );
               Tcrit_cool(i) = interp1( obj.Hcrit{obj.stateFctInd(i),2}(:,1), obj.Tcrit{obj.stateFctInd(i),2}(:,1), Hnorm(i) );
            end
        end
        
        %constructor
        function obj = MagTenseStateFunction( ind )
            
           obj.stateFctInd = ind;
            
           %simply constructor with a single BR state function later to be expanded into n such functions with individual parameters
           nH = 50;
           nT = 1000;
           np = 2;
           
           %units of tesla
           %ferro
           obj.H0{1,1} = linspace(0,2,nH)';
           %para
           obj.H0{1,2} = linspace(0,2,nH)';
           
           %Kelvin
           obj.T0{1,1} = linspace(290,305,nT)';
           obj.T0{1,2} = linspace(290,305,nT)';
           
           %Pa
           obj.p0{1,1} = [1e5 1e6];
           obj.p0{1,2} = [1e5 1e6];
           
           
           obj.M{1,1} = zeros( length(obj.T0{1,1}), length(obj.H0{1,1}), length(obj.p0{1,1}) );
           obj.M{1,2} = zeros( length(obj.T0{1,2}), length(obj.H0{1,2}), length(obj.p0{1,2}) );
           
           obj.Tcrit{1,1} = zeros( length(obj.H0{1,1}), length(obj.p0{1,1}) );
           obj.Hcrit{1,1} = obj.Tcrit{1,1};
           
           obj.Tcrit{1,2} = zeros( length(obj.H0{1,2}), length(obj.p0{1,2}) );
           obj.Hcrit{1,2} = obj.Tcrit{1,2};
           
           %for each pressure find the magnetization for heating and
           %cooling as a function of temperature and field
           for i=1:length(obj.p0{1,1})
               [M_cl,S_cl, M_ht, S_ht] = getBR_MF_properties( obj.T0{1,1}, obj.H0{1,1}, obj.p0{1,1}(i), 295, 1.3, 'LaFeSi' );

               obj.M{1,1}(:,:,i) = M_ht; 
               obj.M{1,2}(:,:,i) = M_cl;
               
               %find the critical temperature for each field and for
               %heating and cooling
               [obj.Tcrit{1,1}(:,i),obj.Hcrit{1,1}(:,i)] = MagTenseStateFunction.getCriticalPoint( obj.T0{1,1}, obj.H0{1,1}, obj.M{1,1}(:,:,i) );
               [obj.Tcrit{1,2}(:,i),obj.Hcrit{1,2}(:,i)] = MagTenseStateFunction.getCriticalPoint( obj.T0{1,2}, obj.H0{1,2}, obj.M{1,2}(:,:,i) );
               
           end
        end
        
        %interpolates and finds M as a function of H, T, p, hysteresis and
        %local state function as given by the indices in ind.
        function Mnorm = getMnorm( obj, Hnorm, T, p, hyst )
            
            %we treat hyst=0 as FM (ferromagnetic) and hyst=1 as PM
            %(paramagnetic). FM is understood as the state which the local
            %lump of the sample is reset at a low temperature, i.e. FM is
            %also known as the heating curve. PM is understood as the state
            %in which the sample is reset at a high temperature, i.e. the
            %cooling curve.
            
            %allocate Mnorm
            Mnorm = zeros(length(Hnorm),1);
            
            %indices into the various state functions
            ind = obj.stateFctInd;
            
            %find unique indices into the state function array
            un_st = unique( ind );
            %loop over the unique state functions applied
            for i=1:length(un_st)
                %find the positions of a primitive is FM and uses the
                %current state function
                FM_IND = find( hyst==MagTenseStateFunction.FM & ind==un_st(i) );
                %similar to above but for paramagnetic
                PM_IND = find( hyst==MagTenseStateFunction.PM & ind==un_st(i) );

                if ~isempty( FM_IND )
                    Mnorm(FM_IND) = interp3( obj.H0{un_st(i),1}, obj.T0{un_st(i),1}, obj.p0{un_st(i),1},...
                                             obj.M{un_st(i),1}, Hnorm(FM_IND), T(FM_IND), p(FM_IND) );
                end
                if ~isempty( PM_IND )
                    Mnorm(PM_IND) = interp3( obj.H0{un_st(i),2}, obj.T0{un_st(i),2}, obj.p0{un_st(i),2},...
                                             obj.M{un_st(i),2}, Hnorm(PM_IND), T(PM_IND), p(PM_IND) );
                end
            end
        end
        
    end
end