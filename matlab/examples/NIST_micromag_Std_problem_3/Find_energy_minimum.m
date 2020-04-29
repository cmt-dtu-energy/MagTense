function Find_energy_minimum
    clearvars
    
    k = 1;
    for i = 5:15
        try
            tic
            [x,~,~,~] = fminsearchbnd(@Find_energy_cross,[8.5 i],[8 i],[9 i]);
            E_time = toc;
            
            L_arr(k,:) = [x toc];
            k = k+1;
        catch
            1
        end
    end
end

function E_diff = Find_energy_cross(arr)
    L = arr(1); i = arr(2);
    [~,~,~,E_arr,~] = Standard_problem_3( [i i i], 1, 0, 0, 0, 0, L ); 
    E_diff = abs(sum(E_arr(:,:,1),1)-sum(E_arr(:,:,2),1));
    
    if (numel(E_diff) > 1)
        1
    end
    if (isnan(E_diff))
        1
    end
end
