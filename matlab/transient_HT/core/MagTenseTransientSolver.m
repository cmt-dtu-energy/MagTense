%Kaspar K. Nielsen, kasparkn@gmail.com, 
%Latest modified 3 May 2020
%Function that performs the main time loop in the transient
%HT-field-structural mechanics solver.
function [solution, geom, setts] = MagTenseTransientSolver( solution, geom, setts, debug )

setts.dt = 0.001;
%figure;
%hold on;
cnt = 0;
while setts.t < setts.t_tot
    
    %update dynamic boundary conditions
    geom = geom.updateBoundaryConditions( setts.t );
    
    [solution, outFlag] = MagTenseIterateTimestep( solution, geom, setts, debug );
    if outFlag == 1
        setts.t = setts.t + setts.dt;
        if cnt==100
            plot( setts.t, solution.T(1), '.k');
            plot( setts.t, solution.T(2), 'vr');
            plot( setts.t, solution.T(2), 'squareb');
            drawnow;
            cnt = 0;
        end
        cnt = cnt + 1;
    else
        disp('Iteration not converging, exiting');
        break;
    end
end


end