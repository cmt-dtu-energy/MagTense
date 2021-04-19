function [x,y,z,dx,dy,dz] = MakeTheGrid(N,Lx,Ly,Lz)
if N(1)>1
    dx = Lx/N(1);
    x = linspace(-Lx/2+dx/2,+Lx/2-dx/2,N(1)) ;
    %dx = x(2)-x(1) ; % dx = 2/(N-1)
else
    x = 0 ;
    dx = Lx ;
end

if N(2)>1
    dy = Ly/N(2);
    y = linspace(-Ly/2+dy/2,+Ly/2-dy/2,N(2)) ;
    %dy = y(2)-y(1) ; % dx = 2/(N-1)
else
    y = 0 ;
    dy = Ly ;
end

if N(3)>1
    dz = Lz/N(3);
    z = linspace(-Lz/2+dz/2,+Lz/2-dz/2,N(3)) ;
    %dz = z(2)-z(1) ; % dx = 2/(N-1)
else
    z = 0 ;
    dz = Lz ;
end