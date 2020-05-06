function [x,y,z,InThis] = DoLloydIteration(X,Y,Z,x,y,z)
IJK = 50 ;
GenPos = [x(:),y(:),z(:)] ;
for ijk = 1:IJK
    
%     [vornb,vorvx] = polybnd_voronoi(pos,bnd_pnts);
    
    AllTheDists = zeros(size(X,1),size(X,2),size(X,3),size(GenPos,1)) ;
    for i=1:size(GenPos,1)
        AllTheDists(:,:,:,i) = sqrt((X-GenPos(i,1)).^2+(Y-GenPos(i,2)).^2+(Z-GenPos(i,3)).^2) ;
    end
    [~,Kmax] = min(AllTheDists,[],4) ;
    for i=1:size(GenPos,1)
        InThis{i} = Kmax==i ;
        Nin(i,ijk) = sum(InThis{i}(:)) ;
        NewPos(i,1) = sum(X(:).*InThis{i}(:))/Nin(i,ijk) ;
        NewPos(i,2) = sum(Y(:).*InThis{i}(:))/Nin(i,ijk) ;
        NewPos(i,3) = sum(Z(:).*InThis{i}(:))/Nin(i,ijk) ;
        
    end
%     disp([num2str(ijk),') ',num2str(sqrt(sum( (NewPos(:)-GenPos(:)).^2)/size(GenPos,1)))]) ;
    GenPos = NewPos ;

end

x = GenPos(:,1) ;
y = GenPos(:,2) ;
z = GenPos(:,3) ;
