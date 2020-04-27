

function mode = getMicroMagProblemMode( nm )

mode = int32(-1);
%converts from name to internal (fortran) representation
switch( nm )
    case 'new'
        mode = int32(1);
    case 'old'
        mode = int32(2);
end


end