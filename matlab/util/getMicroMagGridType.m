

function type = getMicroMagGridType( nm )
%% maps the grid type from name to internal int value

type = int32(-1);
switch ( nm )
    case 'uniform' 
        type = int32(1);
    case 'tetrahedron' 
        type = int32(2);
    case 'unstructuredPrisms' 
        type = int32(3);
end
end