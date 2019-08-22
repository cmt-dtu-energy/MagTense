function PlotStruct = CreatePlotStructure3D(X,Y,Z,PlotStruct)


%% Figure and axes
TheLW = 2 ;
set(gcf,'units','pixel') ;
SqSz = 600 ;
PlotStruct.hF = gcf ;

if PlotStruct.MaxHn~=0 | isfield(PlotStruct,'HystDir')
    set(gcf,'position',[100,100,2*SqSz,SqSz]) ;
    PlotStruct.hA2 = axes('units','pixel','position',[SqSz,0,SqSz,SqSz],'xtick',[-1,0,+1],'ytick',[-1,0,+1],'xticklabel',[],'yticklabel',[],'xgrid','on','ygrid','on') ;
    xlabel(PlotStruct.hA2,'H') ;
    ylabel(PlotStruct.hA2,'M') ;
    
else
    set(gcf,'position',[100,100,SqSz,SqSz]) ;
end
PlotStruct.hA1 = axes('units','pixel','position',[0,0,SqSz,SqSz]) ;

%% Color-surface and arrows

% PlotStruct.hS = surface(X,Y,-.01+0.*X,'linestyle','none','cdatamapping','direct','parent',PlotStruct.hA1) ;
hold on
%PlotStruct.hQ = quiver3(X(:),Y(:),Z(:),1+0.*X(:),0.*X(:),0.*X(:),.5,'color','k','parent',PlotStruct.hA1) ;
PlotStruct.hQ = quiver3(X,Y,Z,1+0.*X,0.*X,0.*X,.5,'color','k','parent',PlotStruct.hA1) ;

set(PlotStruct.hA1,'xlim',[-1,+1],'ylim',[-1,+1],'zlim',[-1,+1],'visible','off') ;

TheData.LastT = 0 ;
TheData.LastPlottedT = 0 ;
set(gcf,'userdata',TheData) ;
axis equal
set(PlotStruct.hA1,'color',.5.*[1,1,1],'xtick',[],'ytick',[])

%% Create the hysteresis loop plot structure
if isfield(PlotStruct,'HystDir')
    PlotStruct.hL6x = line(nan,nan,'linewidth',TheLW,'color','k','parent',PlotStruct.hA2) ;
    PlotStruct.hL7x = line(nan,nan,'linewidth',TheLW,'color','k','marker','o','markerfacecolor','w','parent',PlotStruct.hA2) ;
else
    if PlotStruct.MaxHn~=0
        
        PlotStruct.hL6x = line(nan,nan,'linewidth',TheLW,'color','r','parent',PlotStruct.hA2) ;
        PlotStruct.hL7x = line(nan,nan,'linewidth',TheLW,'color','r','marker','o','markerfacecolor','w','parent',PlotStruct.hA2) ;
        PlotStruct.hL6y = line(nan,nan,'linewidth',TheLW,'color','b','parent',PlotStruct.hA2) ;
        PlotStruct.hL7y = line(nan,nan,'linewidth',TheLW,'color','b','marker','o','markerfacecolor','w','parent',PlotStruct.hA2) ;
        PlotStruct.hL6z = line(nan,nan,'linewidth',TheLW,'color',[0,.5,0],'parent',PlotStruct.hA2) ;
        PlotStruct.hL7z = line(nan,nan,'linewidth',TheLW,'color',[0,.5,0],'marker','o','markerfacecolor','w','parent',PlotStruct.hA2) ;
        
        legend([PlotStruct.hL6x,PlotStruct.hL6y,PlotStruct.hL6z],{'x','y','z'},'location','southeast') ;
    end
end