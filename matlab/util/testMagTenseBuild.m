util_dir = pwd;

%% Test magnetostatic part by running validation scripts and comparing with FEM simulations
USE_CUDA_arr    = [true,false];
USE_CVODE_arr   = [true,false];
USE_RELEASE_arr = [true,false];

for i = 1:length(USE_CVODE_arr)
    for j = 1:length(USE_CUDA_arr)
        for k = 1:length(USE_RELEASE_arr)

        clear functions
    
        USE_CUDA_var = USE_CUDA_arr(i);
        USE_CVODE_var = USE_CVODE_arr(j);
        USE_RELEASE_var = USE_RELEASE_arr(k);
        
        cd(util_dir)
        cd '..'
        buildMagTenseMEX('USE_CUDA',USE_CUDA_var,'USE_CVODE',USE_CVODE_var,'VS_STUDIO',true,'USE_RELEASE',USE_RELEASE_var)

        %% Test micromagnetic part by running the mumag standard problems and comparing with published solutions
        cd(util_dir)
        cd('..\examples\Micromagnetism\mumag_micromag_Std_problem_4')
        [~,~,~,~,~,~,rel_int_error] = Standard_problem_4( 1, [36,9,1], 'USE_CUDA', USE_CUDA_var, 'USE_CVODE', USE_CVODE_var );
        disp(rel_int_error)
        assert(all(rel_int_error < 100))
        end
    end
end


