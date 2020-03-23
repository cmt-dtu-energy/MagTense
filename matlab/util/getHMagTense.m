

function [Hout,Hnorm] = getHMagTense( tiles, X,Y,Z )

    %it is expected that X,Y and Z have the same dimensions

    %convert to (n,3) array acceptable to MagTense
    pts = zeros( numel(X), 3 );
    pts(:,1) = reshape(X,[numel(X),1]);
    pts(:,2) = reshape(Y,[numel(Y),1]);
    pts(:,3) = reshape(Z,[numel(Z),1]);

    H = getHFromTiles_mex( tiles, pts, int32( length(tiles) ), int32( length(pts(:,1)) ) );
    
    %convert H back to the same format as the input X,Y and Z
    sx = size(X);
    
    nx = sx(1);
    ny = 1;
    nz = 1;
    
    if length(sx)==3
        nz = sx(3);
    end
    if length(sx)>=2
        ny = sx(2);
    end
    
        
    
    Hout = zeros( nx, ny, nz, 3 );
    
    Hout(:,:,:,1) = reshape( H(:,1), nx, ny, nz );
    Hout(:,:,:,2) = reshape( H(:,2), nx, ny, nz );
    Hout(:,:,:,3) = reshape( H(:,3), nx, ny, nz );
    
    %Find the norm of the field
    Hnorm = squeeze( sqrt( sum(Hout.^2,4) ) );
    
    Hout = squeeze(Hout);
end