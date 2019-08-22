%given the magnetic field H (n,1) vector and a tile this method returns and
%updated tile where the magnetization (tile.M) has been calculated from the
%permeability tensor. Not norm conserving.
function tile = getM_HardMagnet_1( H, tile )
    %Define the unit vectors of the global coordinate system
    un_x = [1,0,0];
    un_y = [0,1,0];    
    un_z = [0,0,1];    
    
    %Find the M vector in the i'th tile given the total internal field
    %at the i'th tile's center
    %Magnetization in the orthonormal basis of the easy axis
    M = zeros(3,1);
    %Magnetic field along the easy axis
    H_par = dot( H, tile.u_ea );
           
    %Magnetization along the easy axis       
    M(1) = tile.Mrem + (tile.mu_r_ea - 1) * H_par;

    %Magnetic field orthogonal to the easy axis
    H_trans_1 = dot( H, tile.u_oa1 );
    %Magnetization orthogonal to the easy axis       
    M(2) = (tile.mu_r_oa - 1) * H_trans_1;

    %the other magnetic field component orthogonal to the easy axis
    H_trans_2 = dot( H, tile.u_oa2 );
    %Magnetization orthogonal to the easy axis       
    M(3) = (tile.mu_r_oa - 1) * H_trans_2;
    %magnetization of the i'th tile (defined with respect to the global
    %coordinate system)
    tile.M(1) = dot( (M(1) * tile.u_ea + M(2) * tile.u_oa1 + M(3) * tile.u_oa2), un_x );
    tile.M(2) = dot( (M(1) * tile.u_ea + M(2) * tile.u_oa1 + M(3) * tile.u_oa2), un_y );
    tile.M(3) = dot( (M(1) * tile.u_ea + M(2) * tile.u_oa1 + M(3) * tile.u_oa2), un_z );
    
end