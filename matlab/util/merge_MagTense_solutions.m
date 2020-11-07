function merge_MagTense_solutions(org_file,con_file)

    load(org_file)
    elapsedTime2 = elapsedTime;
    problem2 = problem;
    solution2 = solution;
    results2 = results;
    
    load(con_file)
    
    elapsedTime = elapsedTime + elapsedTime2;
    
    results.Mxr = results2.Mxr;
    results.Myr = results2.Myr;

    mu0 = 4*pi*1e-7;
    HystDir = 1/mu0*[1,1,1]/sqrt(3) ;
    HextFct = @(t) HystDir .* t';
    problem = problem.setHext( HextFct, sort(unique([0.5:-0.0005:-0.09]),'descend') );

    size_org = size(solution2.M);
    size_con = size(solution.M);

    comb_arr = zeros([size_org(1) size_org(2) size_org(3)+size_con(3) size_org(4)]);
    comb_arr(:,:,1:size_org(3),:)       = solution2.M;
    comb_arr(:,:,(size_org(3)+1):end,:) = solution.M;
    solution.M = comb_arr;
    
    comb_arr(:,:,1:size_org(3),:)       = solution2.H_exc;
    comb_arr(:,:,(size_org(3)+1):end,:) = solution.H_exc;
    solution.H_exc = comb_arr;
    
    comb_arr(:,:,1:size_org(3),:)       = solution2.H_ext;
    comb_arr(:,:,(size_org(3)+1):end,:) = solution.H_ext;
    solution.H_ext = comb_arr;
    
    comb_arr(:,:,1:size_org(3),:)       = solution2.H_dem;
    comb_arr(:,:,(size_org(3)+1):end,:) = solution.H_dem;
    solution.H_dem = comb_arr;
    
    comb_arr(:,:,1:size_org(3),:)       = solution2.H_ani;
    comb_arr(:,:,(size_org(3)+1):end,:) = solution.H_ani;
    solution.H_ani = comb_arr;

    save('Combined.mat','elapsedTime','problem','problem2','results','solution')
end