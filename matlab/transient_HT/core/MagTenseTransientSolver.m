%Kaspar K. Nielsen, kasparkn@gmail.com, 
%Latest modified 3 May 2020
%Function that performs the main time loop in the transient
%HT-field-structural mechanics solver.
function [solution, geom, setts] = MagTenseTransientSolver( solution, geom, setts, debug )

cnt = 0;
%initial applied field
solution.Happ_old = setts.Happ(setts.t);
while setts.t < setts.t_tot
    
    %update dynamic boundary conditions
    geom = geom.updateBoundaryConditions( setts.t );
    %update the state (hysteresis)
    solution.hyst = setts.Hyst( solution );
    %set the applied field at time t+dt
    solution.Happ_new = setts.Happ(setts.t+setts.dt);
    %run the time step
    [solution, outFlag] = MagTenseIterateTimestep( solution, geom, setts, debug );
    
    
    %save data from the time step
    if outFlag == 1
        %save the new applied field as the old one
        solution.Happ_old = solution.Happ_new;
        %update timestep
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