function A = ComputeExchangeTerm3D(N,dx,dy,dz)
% Returns exchange matrix divided by finite difference factors.
% Calculations use central-difference and Neumann b.c. (equivalent to no
% b.c. when the diagonal is not present, see below).
% N is the number of tiles along each space dimension. dx, dy and dz are
% distances between grid points. The diagonal is empty since the cross
% products of the LL equation make it vanish anyway.
x1=repmat([ones(N(1)-1,1)/dx^2;0],N(2)*N(3),1); % x-dimension neighbors
x2=repmat([0;ones(N(1)-1,1)/dx^2],N(2)*N(3),1);
y1=repmat([ones((N(2)-1)*N(1),1)/dy^2;zeros(N(1),1)],N(3),1); % y-dimension neighbors
y2=repmat([zeros(N(1),1);ones((N(2)-1)*N(1),1)/dy^2],N(3),1);
z1=ones(prod(N),1)/dz^2; % z-dimension neighbors
z2=ones(prod(N),1)/dz^2;
A = spdiags([z1,y1,x1,x2,y2,z2],[-N(1)*N(2),-N(1),-1,+1,+N(1),N(1)*N(2)],prod(N),prod(N));

'' ;
end