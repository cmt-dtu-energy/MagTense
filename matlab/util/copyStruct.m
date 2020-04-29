function str = copyStruct( str_in )

    str = struct();    
    for fn = fieldnames(str_in)'
        str.(fn{1}) = str_in.(fn{1});
    end
end