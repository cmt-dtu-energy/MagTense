function solver = getMicroMagSolver( nm )
%default value
solver = int32(-1);
%the following maps from naming to internal (fortran) representation of the
%solver type
    switch nm
        case 'Explicit'
            solver = int32(1);
        case 'Implicit'
            solver = int32(3);
        case 'Dynamic'
            solver = int32(2);
    end
    
end
