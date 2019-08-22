function [x,y,z,dx,dy,dz] = MakeTheGrid(N,Lx,Ly,Lz)
if N(1)>1
    x = linspace(-Lx/2,+Lx/2,N(1)) ;
    dx = x(2)-x(1) ; % dx = 2/(N-1)
else
    x = 0 ;
    dx = Lx ;
end

if N(2)>1
    y = linspace(-Ly/2,+Ly/2,N(2)) ;
    dy = y(2)-y(1) ; % dx = 2/(N-1)
else
    y = 0 ;
    dy = Ly ;
end

if N(3)>1
    z = linspace(-Lz/2,+Lz/2,N(3)) ;
    dz = z(2)-z(1) ; % dx = 2/(N-1)
else
    z = 0 ;
    dz = Lz ;
end