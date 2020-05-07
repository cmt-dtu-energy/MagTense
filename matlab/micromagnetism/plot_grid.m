function plot_grid(pos_out, dims_out, ndim, Tc_out, vor_c, dTc)
deltaDraw = min(dims_out(:)) ;
if ndim == 2
    
    x = vor_c(:,1);
    y = vor_c(:,2);
    %         z = vor_c(:,3);
    
    [v,c] = voronoin(vor_c);
    
    border = 0:0.01:1;
    for j_inner = 1:4
        clear vor;
        for i_inner = 1:length(x)
            switch j_inner
                case 1
                    vor(i_inner,:) = sqrt((x(i_inner)-border).^2+(y(i_inner)-0).^2);
                case 2
                    vor(i_inner,:) = sqrt((x(i_inner)-border).^2+(y(i_inner)-1).^2);
                case 3
                    vor(i_inner,:) = sqrt((x(i_inner)-0).^2+(y(i_inner)-border).^2);
                case 4
                    vor(i_inner,:) = sqrt((x(i_inner)-1).^2+(y(i_inner)-border).^2);
            end
        end
        for i_inner = 1:length(vor(1,:))
            [~, be(i_inner)] = min(vor(:,i_inner));
        end;
        indx = find(logical(diff(be)));
        indx = [1 indx length(be)-1];
        
        for i_inner = 1:length(indx)
            switch j_inner
                case 1
                    v(end+1,1) = (border(indx(i_inner)+1)-border(indx(i_inner)))/2+border(indx(i_inner));
                    v(end,2) = 0;
                case 2
                    v(end+1,1) = (border(indx(i_inner)+1)-border(indx(i_inner)))/2+border(indx(i_inner));
                    v(end,2) = 1;
                case 3
                    v(end+1,1) = 0;
                    v(end,2) = (border(indx(i_inner)+1)-border(indx(i_inner)))/2+border(indx(i_inner));
                case 4
                    v(end+1,1) = 1;
                    v(end,2) = (border(indx(i_inner)+1)-border(indx(i_inner)))/2+border(indx(i_inner));
            end
            c{be(indx(i_inner))}   = [c{be(indx(i_inner))} length(v(:,1))];
            c{be(indx(i_inner)+1)} = [c{be(indx(i_inner)+1)} length(v(:,1))];
        end
    end
    
    %%make a plot of the Voronoi polyhedrons
%     col=jet(length(dTc));
    col = hsv(length(dTc)) ;
    figure;
    hold on;
    
    for i = 1:length(c)
        %         if all(c{i}~=1)   % If at least one of the indices is 1,
        % then it is an open region and we can't
        % patch that.
        if ndim == 2
            temp = c{i};
            temp = temp(temp > 1);
            c{i} = temp;
            
            phi = atan2(v(c{i},2)-y(i),v(c{i},1)-x(i));
            phi(phi < 0) = phi(phi < 0)+2*pi;
            [~,indx] = sort(phi);
            temp = c{i};
            c{i} = temp(indx);
            patch(v(c{i},1),v(c{i},2),col(i,:));
            plot(x,y,'k.')
        else
            patch(v(c{i},1),v(c{i},2),v(c{i},3),col(i,:));
            plot3(x,y,z,'k.')
        end
    end
    xlim([0,1]);
    ylim([0,1]);
    if ndim==3
        zlim([0,1]);
    end
    
    
    hFgrid = figure ;
    hold on
    
    for i=1:length(pos_out(:,1))
        colInd = find(Tc_out(i) == dTc);
        rectangle('position',[pos_out(i,1)-dims_out(i,1)/2,pos_out(i,2)-dims_out(i,2)/2,dims_out(i,1),dims_out(i,2)],'facecolor',col(colInd,:));
    end
    
    figure
    hold on
    
    for i=1:length(pos_out(:,1))
        colInd = find(Tc_out(i) == dTc);
        rectangle('position',[pos_out(i,1)-dims_out(i,1)/2,pos_out(i,2)-dims_out(i,2)/2,dims_out(i,1),dims_out(i,2)],'facecolor',col(colInd,:),'linestyle','none');
    end
else
    if ndim == 3
        hFgrid = figure('color','w') ;
        hold on
        if 0 
        n=100;
        ind = round(rand([n,1]) * ( length(pos_out(:,1)) - 1 ) + 1);
        
        col = hsv(length(unique(Tc_out)));
        
        colIndex = round( ( Tc_out(ind)-min(Tc_out) )./ ( max(Tc_out) - min(Tc_out) ).*( length(col(:,1)) - 1 ) + 1 );
        
        for i=1:size(pos_out,1)%length(ind)
            %             drawCube(pos_out(ind(i),:),dims_out(ind(i),1),col(colIndex(i),:));
            drawCube(pos_out((i),:),dims_out((i),:),col(colIndex(i),:));
            
        end
        
        else % draw all of them
%             col=jet(length(dTc));
            col=hsv(length(dTc));
            col = 0.*col + (.5 + col)./2
             for i=1:length(pos_out(:,1))
        colInd = find(Tc_out(i) == dTc);
        
         drawCube(pos_out((i),:),dims_out((i),:)-deltaDraw/8,col(colInd,:),max(dims_out(:)));

%         rectangle('position',[pos_out(i,1)-dims_out(i,1)/2,pos_out(i,2)-dims_out(i,2)/2,dims_out(i,1),dims_out(i,2)],'facecolor',col(colInd,:),'linestyle','none');
             end
             
             
                     hold on
%         line(pos_out(:,1),pos_out(:,2),pos_out(:,3),'marker','o','linestyle','none') ;

    axis equal
    view(30,30)
            light('position',2.*[0,1,1])
        %         light('position',[1,1,0])
        light('position',2.*[1,0,1])
        set(gca,'visible','off') ;
%     xlabel('x') ;
%     ylabel('y') ;
%     zlabel('z') ;
'' ;
        end
        
        
    end
end
try
    figure(hFgrid)
end
end

function drawCube(pos,dim,col,TheScale)
TheAlpha = 1 ;
TheLS = 'none' ;
TheLW = 1 ;
% figure ;
% hold on
u = [0,1,1,0,0]-1/2 ;
v = [0,0,1,1,0]-1/2 ;
w = [1,1,1,1,1]-1/2 ;
onetwothree =[1,2,3] ;
twodims = @(j) onetwothree(onetwothree~=j) ;
for idim=1:3
    clear TheXs
    TheDims = twodims(idim) ;
    TheXs(:,TheDims(1)) = u ;
    TheXs(:,TheDims(2)) = v ;
    TheXs(:,idim) = w ;
    if 1 % max(dim)>0
        thisLW = 2*max(dim)/TheScale ;
%         disp(num2str(thisLW))
    else
        thisLW = TheLW ;
    end
    for iii=1:2
        patch('xdata',pos(1)+dim(1).*TheXs(:,1),'ydata',pos(2)+dim(2).*TheXs(:,2),'zdata',pos(3)+dim(3).*TheXs(:,3),'facecolor',col,'facelighting','flat','SpecularStrength',0,'facealpha',TheAlpha,'linestyle',TheLS,'linewidth',thisLW) ;
        
        TheXs(:,idim) = -w ;
    end
    %     patch('xdata',pos(1)+dim(1).*TheXs(:,1),'ydata',pos(2)+dim(2).*TheXs(:,2),'zdata',pos(3)+dim(3).*TheXs(:,3),'facecolor',col,'facelighting','flat','SpecularStrength',0,'facealpha',TheAlpha,'linestyle',TheLS) ;
    
    
    %     hold on
end
'' ;

end
