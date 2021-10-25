function [x,y,z,InThis] = DoLloydIteration(X,Y,Z,x,y,z,IJK)
% 
% DoLloydIteration    Performs Lloyd Iteration on a set of generators in 3D
%
% [x,y,z,InThis] = DoLloydIteration(X,Y,Z,x,y,z) 
% x,y, and z are the positions of the N generators
% X,Y,Z are either 3-D arrays corresponding to the grid over which the
% iteration is performed, or 2-elements arrays corresponding to the
% boundaries of the box over which the iteration is performed
% InThis is a N-element cell array of booleans, indicating whether a point of the grid 
% belongs to a given cell or not
%
% [x,y,z,InThis] = DoLloydIteration(X,Y,Z,x,y,z,IJK)
% additionally specify the number of steps IJK (default value is 50)

%%
if ~exist('IJK','var')
IJK = 50 ;
end
N = numel(x) ;

GenPos = [x(:),y(:),z(:)] ; % pack x,y,z 
%% Check if {X,Y,Z} is a grid or a box
if (numel(X)==2) & (numel(Y)==2) & (numel(Z)==2) % if X,Y,Z are the boundaries of the box (i.e. 2-elemtns arrays)
    NN = 1000 ; % Total Number of Points in the grid
    Lx = (X(2)-X(1)) ;
    Ly = (Y(2)-Y(1)) ;
    Lz = (Z(2)-Z(1)) ;
    Vtot = Lx*Ly*Lz ;
    dV = Vtot/NN ;
    dx = dV^(1/3) ;
    Nx = ceil(Lx/dx) ;
    Ny = ceil(Ly/dx) ;
    Nz = ceil(Lz/dx) ;
    xx = linspace(X(1),X(2),Nx) ;
    yy = linspace(Y(1),Y(2),Ny) ;
    zz = linspace(Z(1),Z(2),Nz) ;
    [X,Y,Z] = ndgrid(xx,yy,zz) ;
    
end
%% Perform Lloyd Iteration
for ijk = 1:IJK
    
    AllTheDists = zeros(size(X,1),size(X,2),size(X,3),size(GenPos,1)) ;
    for i=1:N
        AllTheDists(:,:,:,i) = sqrt((X-GenPos(i,1)).^2+(Y-GenPos(i,2)).^2+(Z-GenPos(i,3)).^2) ;
    end
    [~,Kmax] = min(AllTheDists,[],4) ;
    for i=1:N
        InThis{i} = Kmax==i ;
        Nin(i,ijk) = sum(InThis{i}(:)) ;
        NewPos(i,1) = sum(X(:).*InThis{i}(:))/Nin(i,ijk) ;
        NewPos(i,2) = sum(Y(:).*InThis{i}(:))/Nin(i,ijk) ;
        NewPos(i,3) = sum(Z(:).*InThis{i}(:))/Nin(i,ijk) ;
        
    end

    GenPos = NewPos ;

end
% unpack x,y,z 
x = GenPos(:,1) ;
y = GenPos(:,2) ;
z = GenPos(:,3) ;
