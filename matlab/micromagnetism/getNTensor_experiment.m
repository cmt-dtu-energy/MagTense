%% tiles is an array of length m consisting of the tiles that produce the
%% field at the n pts ( (n,3) array )
%% returns N which is a matrix size (3m,3n) that works like this:
%% N dot M = H, where M is a column vector size 3m with the magnetization
%% components of each tile and H is a column vector with the field
%% components at each measured point.
function [N] = getNTensor_experiment( tiles, pts )

%get the no. of tiles
m = length( tiles );

%get the no. of measured points
n = length( pts(:,1) );

%Make the output matrix
N = zeros( 3*n, 3*m );
% N2 = zeros( 3*n, 3*m );

%loop over each tile
for i=1:m
    %remember to multiply by -1 as the demag tensor for prisms is
    %defined as H = - N*M in the code
    N_ = 1.*getNFromTile_mex( tiles(i), pts, int32(n) );
%     for j=1:n
%         N( (j-1)*3+1:j*3, (i-1)*3+1:i*3 ) = squeeze( N_(j,:,:) );
%     end
    N( :, (i-1)*3+1:i*3 ) = reshape(permute(N_,[3 1 2]),n*3,3);
end

end