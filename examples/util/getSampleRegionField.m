function [H,N,pts] = getSampleRegionField( tiles, N )
% 
%     r_disc = 0.002;
%     theta=linspace(0.001,pi-0.001,20);
%     phi=linspace(0.001,2*pi-0.001,32);
%     [PHI,THETA]=meshgrid(phi,theta);
%     x = cos(PHI) .* sin(THETA) .* r_disc;
%     y = sin(PHI) .* sin(THETA) .* r_disc;
%     z = cos(THETA) .* r_disc;
%     pts(:,1) = reshape( x,[640,1] );
%     pts(:,2) = reshape( y,[640,1] );
%     pts(:,3) = reshape( z,[640,1] );


    r_disc = 0.002;
    phi = linspace(0.001,2*pi-0.001,32);
    z = linspace(-r_disc,r_disc,20);
    
    [PHI,Z] = meshgrid( phi, z );
    x = r_disc .* cos( PHI );
    y = r_disc .* sin( PHI );
    z = Z;
    
    pts = zeros(640,1);
    pts(:,1) = reshape( x, [640,1] );
    pts(:,2) = reshape( y, [640,1] );
    pts(:,3) = reshape( z, [640,1] );

%     nx = 10;
%     ny = 10;
%     nz = 20;
%     x = linspace( -2.5e-3,2.5e-3,nx);
%     y = linspace( -2.6e-3,2.6e-3,ny);
%     z = linspace( -3.5e-3,3.5e-3,nz);
%     [X,Y,Z] = meshgrid(x,y,z);
% 
%     pts = zeros( nx*ny*nz, 3 );
% 
%     pts(:,1) = reshape( X, [nx*ny*nz,1] );
%     pts(:,2) = reshape( Y, [nx*ny*nz,1] );
%     pts(:,3) = reshape( Z, [nx*ny*nz,1] );
    if ~exist( 'N', 'var' ) || isempty(N)
        N = getNTensor_experiment( tiles, pts );
    end
   
    H = reshape( N * [tiles.M]',[3,length(pts(:,1))])';
    Hnorm = sqrt(sum( H.^2,2) ); %for testing
end