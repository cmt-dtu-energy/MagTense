%%Calls the fortran code and returns the iterated tiles.
%%will convert any colors not represented as rgb triplets to rgb triplets
%%for compatibility
%%resumeIteration is optional. If present, it will be passed (0 means no
%otherwise yes, resume the iteration (assuming the tiles' current
%magnetization vectors).
function tiles = IterateMagnetization( tiles, stateFunction, T, max_err, nIte, resumeIteration )

    for i=1:length(tiles)
        if ~isa(tiles(i).color,'double')
            %assume it is a char and convert
            [clr,tiles(i).color] = colornames('Matlab',tiles(i).color);
        end
    end
    
    if  ~exist('resumeIteration','var')
        resumeIteration = 0.;
    end

    tiles = IterateMagnetization_mex( tiles, int32(length(tiles)), stateFunction, int32(length(stateFunction)), T, max_err, int32(nIte), resumeIteration);
end