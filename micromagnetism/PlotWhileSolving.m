function PlotStruct = PlotWhileSolving(N,X,Y,Z,Lx,Ly,Lz,Sigma,MaxHx,MaxHy,MaxHz,K0,Kx,Ky,Kz,MU0,DrawIt,DeltaT,SaveGif,GifFilename)

%% Initialize Plots
PlotStruct.DrawIt = DrawIt ;
PlotStruct.MaxHn = sqrt(MaxHx^2+MaxHy^2+MaxHz^2)/MU0 ; % TODO: Fix inelegance - Use MU0 only one place.
% if exist('HystDir', 'var')
%     PlotStruct.HystDir = HystDir ;
% end
if DrawIt
    hF = figure('units','normalized','toolbar','none') ;
    set(hF,'color',[.5,.5,.5]) ;
    PlotStruct = CreatePlotStructure3D(X,Y,Z,PlotStruct) ;
end
PlotStruct.Lx = Lx ;
PlotStruct.Ly = Ly ;
PlotStruct.Lz = Lz ;

PlotStruct.DeltaT = DeltaT ;
% PlotStruct.DeltaT = inf ; % shows nothing
PlotStruct.SaveGif = SaveGif ;
PlotStruct.GifFilename = GifFilename ; % Specify the output file name

%% Initialize gif
ThatS = 1.5 ;
if PlotStruct.SaveGif
%     ThisCData =  reshape(ColorFromHorPsiTheta01(SigmaX,SigmaY,SigmaZ),N,N,N,3) ;
%     set(PlotStruct.hS,'cdata',ThisCData) ;
    NN = round(numel(Sigma)/3) ;
    SigmaX = Sigma(0*NN+[1:NN]) ;
    SigmaY = Sigma(1*NN+[1:NN]) ;
    SigmaZ = Sigma(2*NN+[1:NN]) ;
    set(PlotStruct.hQ,'udata',reshape(SigmaX,N(1),N(2),N(3)),'vdata',reshape(SigmaY,N(1),N(2),N(3)),'wdata',reshape(SigmaZ,N(1),N(2),N(3))) ;
    %set(PlotStruct.hQ,'udata',SigmaX,'vdata',SigmaY,'wdata',SigmaZ) ;
    set(gca,'xlim',ThatS.*(PlotStruct.Lx/2).*[-1,+1],'ylim',ThatS.*(PlotStruct.Ly/2).*[-1,+1],'zlim',ThatS.*(PlotStruct.Lz/2).*[-1,+1]) ;
    drawnow
% disp(num2str(mean(sqrt(
    frame = getframe(gcf);
    im = frame2im(frame);
    [Agif,map] = rgb2ind(im,256);
    imwrite(Agif,map,PlotStruct.GifFilename,'gif','LoopCount',Inf,'DelayTime',1/30);
    
end
%% Plot Grains
Nmax = 0; % Multigrain plotting not supported
if Nmax > 1
    TheQLW = 2 ;
    [vx,vy] = voronoi(xPoint,yPoint) ;
    for jk = 1:size(vx,2)
        line(vx(:,jk),vy(:,jk),.2.*[1,1],'color','k','linewidth',2,'linestyle','--','parent',PlotStruct.hA1) ;
    end
    hold on
    
    quiver3(xPoint,yPoint,.2+0.*yPoint,+Hx,+Hy,+Hz,.2,'color','k','linewidth',TheQLW,'parent',PlotStruct.hA1) ;
    quiver3(xPoint,yPoint,.2+0.*yPoint,-Hx,-Hy,-Hz,.2,'color','k','linewidth',TheQLW,'parent',PlotStruct.hA1) ;
    set(gca,'xlim',[-1,+1],'ylim',[-1,+1],'zlim',[-1,+1]) ;
else
    TheQLW = 2 ;
    if abs(K0)>0
        TheSc = .2 ;
        TheSc0 = 1 ;
        quiver3(0,0,.2,+TheSc*Kx,+TheSc*Ky,0,TheSc0,'color','k','linewidth',TheQLW,'parent',PlotStruct.hA1) ;
        quiver3(0,0,.2,-TheSc*Kx,-TheSc*Ky,0,TheSc0,'color','k','linewidth',TheQLW,'parent',PlotStruct.hA1) ;
        
        set(gca,'xlim',[-1,+1],'ylim',[-1,+1],'zlim',[-1,+1]) ;
    end
end
end