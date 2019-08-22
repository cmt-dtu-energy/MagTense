%%This function is meant to iterate the tiles (input) for finding the
%%self-consistent solution in terms of the tiles' magnetization vectors and
%%their mutual interaction. If the second argument, N, is not provided the
%%demagnetization tensor field will be calculated first (and returned for
%%subsequent use). If N is provided the iteration will be done directly
%%using N thus saving significant computational time.
%%This method is only verified for cylindrical piece tiles at the moment
%%(25-2-2019, kaki).
function [tiles,N] = iterateFixedGeometry( tiles, N, verb )

    %maximum allowed relative error
    err_max = 1e-10;
    
    %no. of tiles
    n = length(tiles);
    if ~exist( 'N', 'var' ) || isempty(N)
        %center points of the tiles
        pts = zeros(n,3);
        r0 = [tiles.r0];
        theta0 = [tiles.theta0];

        pts(:,1) = r0 .* cos(theta0);
        pts(:,2) = r0 .* sin(theta0);
        pts(:,3) = [tiles.z0];

        %add the offset
        pts = pts + reshape( [tiles.offset], [3,n] )';
        %get the N-tensor
        N = getNTensor_experiment( tiles, pts );
    end
    
    %Set the initial field
    H = zeros(n,3);
    M_norm_old = zeros(n,1);
    err = 1;
    %iteration loop
    while err>err_max
        %update the magnetization of the tiles
        for i=1:n
            tiles(i) = getM_HardMagnet_1( H(i,:), tiles(i) );
        end

        %get the magnetization        
        M = reshape( [tiles.M], [3*n,1] );
        
        %get B from the tensor
        B = N * M;
        %get H by subtracting M
        H = B - M;
        %reshape H to be (n,3)
        H = reshape( H, [3,n] )';
        Mnorm = sqrt( sum( reshape([tiles.M],[3,n])'.^2,2));
        
        k = abs( (Mnorm-M_norm_old)./M_norm_old );
        if isempty(find(isfinite(k)))
            err = 1;
        else
            err = max(k(isfinite(k)));
        end
        M_norm_old = Mnorm;
        if exist('verb','var')
            disp(['err = ' num2str(err)]);
        end
    end

end