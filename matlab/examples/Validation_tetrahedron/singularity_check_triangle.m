%close all
clear all


v = [[1,0,0];[0,1,0];[-1,0,0];[0,0,-1]]';

M = [0,0,1]';

%x = linspace(-.205,-.195,1000);
x = linspace(-1.2,1.5,1000);
y = 1.2;%linspace(-1.2,1.5,100);

[X,Y]=meshgrid(x,y);

H = zeros(length(x),length(y),3);
r = [0,0,0]';

for i=1:length(x)
    for j=1:length(y)
        r(1)=x(i);
        r(2)=y(j);
        
        [N,P] = getN_Triangle( v, r );
        
        H(i,j,:) = N * M;
        
    end
end
length(find( imag(H)~=0 | ~isfinite(H)))
H(imag(H)~=0) = NaN;
Hnorm = sqrt(sum(H.^2,3));
[fg,ax] = getFigure(true);
surf_and_con(X,Y,H(:,:,1)',ax);
[fg,ax] = getFigure(true);
surf_and_con(X,Y,H(:,:,2)',ax);
[fg,ax] = getFigure(true);
surf_and_con(X,Y,H(:,:,3)',ax);
[fg,ax] = getFigure(true);
surf_and_con(X,Y,Hnorm',ax);

[fg,ax] = getFigure(true);
surf_and_con(X,Y,Hnorm',ax);
quiver(X,Y,H(:,:,1)',H(:,:,2)',10);

axis equal

