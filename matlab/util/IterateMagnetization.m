%%Calls the fortran code and returns the iterated tiles.
%%will convert any colors not represented as rgb triplets to rgb triplets
%%for compatibility
%%Tiles is an array of tiles to be passed to the solver
%%stateFunction is a struct encapsulating the state function (used for soft
%%magnets)
%%T is the temperature of the tiles. Needs to be updated at a later time to
%%be an array
%%max_err is the maximum relative error to be acceptable before terminating
%%the iteration
%%nIte is the max.no. of iterations to run before returning.
function tiles = IterateMagnetization( tiles, stateFunction, T, max_err, nIte )

    if isempty( stateFunction )
        %%Make a default stateFunction
        stateFunction = getStateFunctionSoftIron();
    end
    
    if isempty(T)
        T = 300;
    end

    for i=1:length(tiles)
        if ~isa(tiles(i).color,'double')
            %assume it is a char and convert
            [clr,tiles(i).color] = colornames('Matlab',tiles(i).color);
        end
    end

    tiles = IterateMagnetization_mex( tiles, int32(length(tiles)), stateFunction, int32(length(stateFunction)), T, max_err, int32(nIte));
end