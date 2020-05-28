%Kaspar K. Nielsen, kasparkn@gmail.com, 
%Latest modified 3 May 2020
%Function that performs the main time loop in the transient
%HT-field-structural mechanics solver.
function [solution, geom, setts] = MagTenseTransientSolver( solution, geom, setts, debug )

cnt = 0;
while setts.t < setts.t_tot
    
    %update dynamic boundary conditions
    geom = geom.updateBoundaryConditions( setts.t );
    
    [solution, outFlag] = MagTenseIterateTimestep( solution, geom, setts, debug );
    if outFlag == 1
        setts.t = setts.t + setts.dt;
        if cnt==100
            solution = solution.StoreSolution( setts.t );
            cnt = 0;
            disp(['Storing... t = ' num2str(setts.t) ', t_tot = ' num2str(setts.t_tot)]);
        end
        cnt = cnt + 1;
    else
        disp('Iteration not converging, exiting');
        break;
    end
end


end