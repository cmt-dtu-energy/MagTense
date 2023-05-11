function Find_energy_minimum
    clearvars
    options = optimset('TolFun',1e-3,'TolX',1e-3);
    k = 1;
    for i = 5:20
        tic
        [x,~,~,~] = fminsearchbnd(@Find_energy_cross,[8.5 i],[8 i],[9 i],options);
        E_time = toc;

        L_arr(k,:) = [i x E_time];
        k = k+1;
        save('Standard_problem_3_E_cross.mat','L_arr')
    end
end

function E_diff = Find_energy_cross(arr)
    L = arr(1); i = arr(2);
    options.ShowTheResult = false;
    [~,~,~,E_arr,~] = Standard_problem_3( [i i i], L, options ); 
    E_diff = abs(sum(E_arr(:,:,1),1)-sum(E_arr(:,:,2),1));
    
    if (numel(E_diff) > 1)
        1
    end
    if (isnan(E_diff))
        1
    end
end
