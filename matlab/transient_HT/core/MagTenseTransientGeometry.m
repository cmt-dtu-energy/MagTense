

classdef MagTenseTransientGeometry
   
    properties
       dV %volume array containing the volume of each cell. Should be a single column vector of length n
       
       R_geom %sparse matrix of dimensions (n,n) that contains the geometrical part of the thermal resistance from the center of i'th cell to the center of the boundary between the i'th cell and the j'th cell
        
       FaceConditions %is an array with n rows and m columns where m is the largest no. of face-conditions (boundary conditions) a single primitive in the given problem has
       FaceConditionsBoundaryStatic % static part of boundary condition values, size n x m (sparse)
       FaceConditionsBoundaryDynamic %dynamic (changes potentially at each time step) part of the boundary conditions, sparse n x m
    end
    
    properties (Constant)
       %these indices are used for identifying what the entries in FaceConditions mean
       FC_DIRICHLET = 1; %Dirichlet boundary condition (constant value at boundary)
       FC_NEUMANN = 2; %von Neumann boundary condition (constant flux)
       FC_INTERNAL = 11; %interval boundary condition, i.e. contact with another primitive
    end
    
    
    methods
        function obj = MagTenseTransientGeometry( tiles )
        
            obj = initializeBoundaryConditions(obj, tiles);
            
        end
        
        %Takes an array of tiles as input and expects each tile to carry a
        %description of its boundary conditions (to the
        %ambient/environment) as well as connectivity with adjacent tiles
        function obj = initializeBoundaryConditions( obj, tiles )
            
            %no. of tiles
            n = length(tiles);
            
            %assuming no more than 100 interactions through the boundaries
            %of a single tile we make space for the sparse matrices
            %containing the internal and external boundary conditions
            
            obj.FaceConditions = sparse(n,100);
            obj.FaceConditionsBoundaryStatic = obj.FaceConditions;
            obj.FaceConditionsBoundaryDynamic = obj.FaceConditions;
            
            %volume of each tile, should be n x 1
            obj.dV = MagTenseTilesUtil.getVolume( tiles );
            
            %thermal resistance, i.e. internal boundary condition due to
            %finite heat transfer
            obj.R_geom = sparse(n,n);
            
           %setup the static part of the boundary conditions
           for i=1:n
               
               for j=1:length(tiles(i).bdryCdts)
                  bdry = tiles(i).bdryCdts(j);
                  switch bdry.Type
                      case MagTenseTransientGeometry.FC_DIRICHLET
                          %add a Dirichlet condition
                          obj.FaceConditions(i,j) = MagTenseTransientGeometry.FC_DIRICHLET;
                          %represents the thermal resistance to the face of
                          %the given boundary condition
                          obj.FaceConditionsBoundaryStatic(i,j) = bdry.A/bdry.l;
                      case MagTenseTransientGeometry.FC_NEUMANN
                          %add a von Neumann condition
                          obj.FaceConditions(i,j) = MagTenseTransientGeometry.FC_NEUMANN;
                          %only a dynamic part, i.e. the (heat) flux across
                          %the boundary
                      case MagTenseTransientGeometry.FC_INTERNAL
                          obj.FaceConditions(i,j) = MagTenseTransientGeometry.FC_INTERNAL;
                          obj.R_geom(i,bdry.n_ind) = bdry.l/bdry.A;
                  end
                  
               end
           end
        end
        
        %should be called once per time step (at the beginning) in order to
        %update / evolve the boundary conditions in time
        function obj = updateBoundaryConditions( obj, t )
            %find all Dirichlet conditions
            [indx,indy] = find ( obj.FaceConditions == MagTenseTransientGeometry.FC_DIRICHLET );
            obj.FaceConditionsBoundaryDynamic(indx,indy) = 290;
        end
        
        %returns the boundary conditions in a part that goes to the
        %left-handside of the equation (implicit part) and a part that goes
        %to right-hand side (explicit part)
        %takes the thermal conductivity, k (n x 1) as
        %input
        function [lhs_bdry,rhs_bdry] = getBoundaryConditions( obj, k )
            n = length( obj.dV );
            lhs_bdry = zeros(n,1);
            rhs_bdry = zeros(n,1);
            
            %find all Dirichlet conditions
            [indx,indy] = find ( obj.FaceConditions == MagTenseTransientGeometry.FC_DIRICHLET );
            
            %implementation of the Dirichlet condition where a fixed
            %temperature is imposed on a boundary. This becomes a flux on
            %the right-hand side of the unsteady equation: dV*rho*c*dT/dt =
            %...-(T-T_bdry) * A * k / l
            %FaceConditionsBoundaryDynamic is supposed to contain the
            %boundary temperature in this case
            lhs_bdry(indx) = k(indx) .* sum( obj.FaceConditionsBoundaryStatic(indx,indy), 2 );
            rhs_bdry(indx) = lhs_bdry(indx) .* sum( obj.FaceConditionsBoundaryDynamic(indx,indy), 2 ) ;
        end
    end
    
end