%% Benchmark MagTense, both Matlab and Fortran versions, in terms of computation time using Standard problem #4

addpath('../mumag_micromag_Std_problem_2');

%--- Use field X from the mumag example
res_arr = 1:9;

%--- Loop over using the code in both CPU and CUDA mode
for j = 2%1:2
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
        SaveTheResult = 0;
        run_single_curve = 0;
        run_MrHc = 1;
    
        [elapsedTime,problem,solution,results] = Standard_problem_2(resolution, use_CUDA, ShowTheResult, SaveTheResult, run_single_curve, run_MrHc );

        mkdir([ dir_f ])
        save([ dir_f '\Solution_' num2str(resolution(1)) '_' num2str(resolution(2)) '_' num2str(resolution(3)) '_' j_s '.mat'],'results','elapsedTime')

    end
end


