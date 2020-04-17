%% Benchmark MagTense, both Matlab and Fortran versions, in terms of computation time using Standard problem #4

addpath('../NIST_micromag_Std_problem_4');

%--- Use field X from the NIST example
NIST_field = 2;
res_arr = 1:9;

%--- Loop over using the code in both CPU and CUDA mode
for j = 1:2
    if j == 1
        use_CUDA = false;
        j_s = 'no_CUDA';
        dir_f = 'CPU';
        dir_m = 'CPU';
        disp('Running no CUDA')
    else
        use_CUDA = true;
        j_s = 'with_CUDA';
        dir_f = 'CUDA';
        dir_m = 'GPU';
         disp('Running with CUDA')
    end
    
    %--- Loop over the resolution of the problem
    for i = res_arr

        clearvars -except i use_CUDA j_s dem_thres dem_string dir_f dir_m NIST_field res_arr
        
        %--- Vary the resolution
        resolution = [i*20,i*5,1];
    
        ShowTheResult = 0;
        SaveTheResult = 1;
    
        [elapsedTime_part1,elapsedTime_part2,problem_ini,solution_ini,problem_dym,solution_dym] = Standard_problem_4( NIST_field, resolution, use_CUDA, ShowTheResult, SaveTheResult );
    
        mkdir(['Field_' num2str(NIST_field) '\' dir_f ])
        save(['Field_' num2str(NIST_field) '\' dir_f '\Solution_' num2str(resolution(1)) '_' num2str(resolution(2)) '_' num2str(resolution(3)) '_' j_s '.mat'],'problem_ini','solution_ini','problem_dym','solution_dym','elapsedTime_part1','elapsedTime_part2')

    end
end


