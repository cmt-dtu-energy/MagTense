function CartesianUnstructuredMeshPlot(pos, dims, GridInfo, iIn, fnamesave)
%
% CartesianUnstructuredMeshPlot   plots Cartesian unstructured mesh


hF = figure('position',[0 0 600 600],'Color',[1 1 1]);
% ppsz = .2.*[20,19] ;
% ppps = .2.*[0,0,20,19] ;

% set(gcf,'PaperUnits','centimeters','PaperSize',ppsz,'PaperPosition',ppps) ;

view(30,30) ;
hold on
%%
cols = [hsv(numel(iIn)-1);[1,1,1].*0.4];
iOne = sum(abs(GridInfo.TheSigns),1) == 1 ;
iOneInd = find(iOne) ;
for ik=1:numel(iOneInd)
    k = iOneInd(ik) ;
    n= find(GridInfo.TheSigns(:,k)) ;
    if numel(iIn)>1
    for j=1:numel(iIn)
        if ismember(n,iIn{j})
            thisCol = cols(j,:) ;
        end
    end
    TheLS = '-' ;
    ColProperty = 'facecolor' ;
    else
        if size(iIn{1},2)== 3
            thisCol = iIn{1}(n,:,1) ;
            ColProperty = 'facecolor' ;
            TheLS = 'none' ;
        else
        thisCol = (iIn{1}(n)/max(iIn{1}(:))) ;      
            ColProperty = 'cdata' ;
            TheLS = '-' ;
        end
    end
    if abs(GridInfo.fNormX(k))
      hP(ik) = patch(pos(n,1)+[1,1,1,1].*GridInfo.fNormX(k).*dims(n,1)./2,...
            pos(n,2)+([0,0,1,1]-1/2).*dims(n,2),...
            pos(n,3)+([0,1,1,0]-1/2).*dims(n,3),ik,ColProperty,thisCol,'linestyle',TheLS) ;
    end
    if abs(GridInfo.fNormY(k))
        hP(ik) = patch(pos(n,1)+([0,0,1,1]-1/2).*dims(n,1),...
            pos(n,2)+[1,1,1,1].*GridInfo.fNormY(k).*dims(n,2)./2,...
            pos(n,3)+([0,1,1,0]-1/2).*dims(n,3),ik,ColProperty,thisCol,'linestyle',TheLS) ;
    end
    
    if abs(GridInfo.fNormZ(k))
        hP(ik) = patch(pos(n,1)+([0,1,1,0]-1/2).*dims(n,1),...
            pos(n,2)+([0,0,1,1]-1/2).*dims(n,2),...
            pos(n,3)+[1,1,1,1].*GridInfo.fNormZ(k).*dims(n,3)./2,ik,ColProperty,thisCol,'linestyle',TheLS) ;
    end
    '' ;
%     drawnow ;
    '' ;
end
%%
dd = 0 ; % max(dims(:))*100 ;
set(gca,'xlim',[min(GridInfo.Xf(:))-dd/20,max(GridInfo.Xf(:))+dd/20],...
    'ylim',[min(GridInfo.Yf(:))-dd/20,max(GridInfo.Yf(:))+dd/20],...
    'zlim',[min(GridInfo.Zf(:))-dd/20,max(GridInfo.Zf(:))+dd/20])
view(30,30) ;
axis equal
hold on

light('position',[0,0,2].*[max(GridInfo.Xf(:)),max(GridInfo.Yf(:)),max(GridInfo.Zf(:))])
light('position',[0,-2,0].*[max(GridInfo.Xf(:)),max(GridInfo.Yf(:)),max(GridInfo.Zf(:))],'color',.3.*[1,1,1])
xlabel('x') ; ylabel('y') ; zlabel('z') ;
%     plot3(GridInfo.Xf(:),GridInfo.Yf(:),GridInfo.Zf(:),'.')
%%
if numel(iIn)==1
    if size(iIn{1},3)> 1
        if exist('fnamesave','var')
            SaveGifFrames(hF,[fnamesave,'Animation'],1) ;
        end
        for ii = 2:size(iIn{1},3)
            for ik=1:numel(iOneInd)
                k = iOneInd(ik) ;
                n= find(GridInfo.TheSigns(:,k)) ;
                
                set(hP(ik),'facecolor',iIn{1}(n,:,ii)) ;
            end
            
            if exist('fnamesave','var')
                SaveGifFrames(hF,[fnamesave,'Animation'],ii) ;
            end
        end
    end
end

%%
if exist('fnamesave','var')
%     set(hF,'PaperUnits','centimeters','PaperSize',ppsz,'PaperPosition',ppps) ;
%     set(hF,'PaperUnits','inches') ;
    figure(hF) ;
    eval(['print -dtiff ',fnamesave,'Cartesian.tiff']) ;
    savefig(hF,[fnamesave,'Cartesian.fig'])
end
end

function SaveGifFrames(hF,filename,n)

filename = [filename,'.gif'] ;
DelayTime = .1 ;
frame = getframe(hF);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
if n == 1
    if exist(filename,'file')
        delete(filename) ;
    end
    imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',DelayTime);
else
    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',DelayTime);
end
end