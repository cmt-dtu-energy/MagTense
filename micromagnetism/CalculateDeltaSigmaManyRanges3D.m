function dSigma = CalculateDeltaSigmaManyRanges3D(t,Sigma,N,AvrgMatrix,CopyMatrix,Mfact,DemagTensor,alpha,gamma,PlotStruct,AA,dx,dy,dz)

NN = round(numel(Sigma)/3) ;
% N = round(NN^(1/3)) ;
SigmaX = Sigma(0*NN+[1:NN]) ;
SigmaY = Sigma(1*NN+[1:NN]) ;
SigmaZ = Sigma(2*NN+[1:NN]) ;

%% Demagnetization: Long range (coarse grids)

HmXcTOT = zeros(NN,1) ;
HmYcTOT = zeros(NN,1) ;
HmZcTOT = zeros(NN,1) ;
for k=2:numel(AvrgMatrix)
    % get sigma over coarse grid (average)
    SigmaXC = AvrgMatrix{k}*SigmaX ;
    SigmaYC = AvrgMatrix{k}*SigmaY ;
    SigmaZC = AvrgMatrix{k}*SigmaZ ;
    % calculate demag. field over coarse grid   

    HmXc = DemagTensor.KglobXX{k}*SigmaXC+DemagTensor.KglobXY{k}*SigmaYC+DemagTensor.KglobXZ{k}*SigmaZC ;
    HmYc = DemagTensor.KglobYX{k}*SigmaXC+DemagTensor.KglobYY{k}*SigmaYC+DemagTensor.KglobYZ{k}*SigmaZC ;
    HmZc = DemagTensor.KglobZX{k}*SigmaXC+DemagTensor.KglobZY{k}*SigmaYC+DemagTensor.KglobZZ{k}*SigmaZC ;
    '' ; 

    % get demag. field over fine grid (copy)   
    HmXcTOT = HmXcTOT - Mfact*CopyMatrix{k}*HmXc ;  % Coarser
    HmYcTOT = HmYcTOT - Mfact*CopyMatrix{k}*HmYc ;  % Coarser
    HmZcTOT = HmZcTOT - Mfact*CopyMatrix{k}*HmZc ;  % Coarser
end

%%  Exchange 
ThisHjX = AA.HjX(SigmaX,SigmaY,SigmaZ,t) ;
ThisHjY = AA.HjY(SigmaX,SigmaY,SigmaZ,t) ;
ThisHjZ = AA.HjZ(SigmaX,SigmaY,SigmaZ,t) ;

%% External field
ThisHhX = AA.HhX(SigmaX,SigmaY,SigmaZ,t) ;
ThisHhY = AA.HhY(SigmaX,SigmaY,SigmaZ,t) ;
ThisHhZ = AA.HhZ(SigmaX,SigmaY,SigmaZ,t) ;
%% Demagnetization
if numel(CopyMatrix)>0
    ThisHmX = AA.HmX(SigmaX,SigmaY,SigmaZ,t) + HmXcTOT ; % fine + coarser
    ThisHmY = AA.HmY(SigmaX,SigmaY,SigmaZ,t) + HmYcTOT ; % fine + coarser
    ThisHmZ = AA.HmZ(SigmaX,SigmaY,SigmaZ,t) + HmZcTOT ; % fine + coarser
else
    ThisHmX = AA.HmX(SigmaX,SigmaY,SigmaZ,t) ; % fine 
    ThisHmY = AA.HmY(SigmaX,SigmaY,SigmaZ,t) ; % fine
    ThisHmZ = AA.HmZ(SigmaX,SigmaY,SigmaZ,t) ; % fine
end
%% Anisotropy
ThisHkX = AA.HkX(SigmaX,SigmaY,SigmaZ,t) ;
ThisHkY = AA.HkY(SigmaX,SigmaY,SigmaZ,t) ;
ThisHkZ = AA.HkZ(SigmaX,SigmaY,SigmaZ,t) ;

ThisHeffX = ThisHjX + ThisHhX + ThisHmX + ThisHkX ;
ThisHeffY = ThisHjY + ThisHhY + ThisHmY + ThisHkY ;
ThisHeffZ = ThisHjZ + ThisHhZ + ThisHmZ + ThisHkZ ;

%% Calculate Precession and Damping terms from Heff

% Compute m x heff (Precession term)
TheCrossX = -(SigmaY.*ThisHeffZ - SigmaZ.*ThisHeffY) ;
TheCrossY = -(SigmaZ.*ThisHeffX - SigmaX.*ThisHeffZ) ;
TheCrossZ = -(SigmaX.*ThisHeffY - SigmaY.*ThisHeffX) ;
% Compute m x m x heff (Damping term)
ThisHeffX2 = +SigmaY.*TheCrossZ-SigmaZ.*TheCrossY ;
ThisHeffY2 = +SigmaZ.*TheCrossX-SigmaX.*TheCrossZ ;
ThisHeffZ2 = +SigmaX.*TheCrossY-SigmaY.*TheCrossX ;

%% Calculate time-derivative of Sigma

if isequal(class(alpha),'double') % alpha is time-independent
    dSigmaX = alpha.*ThisHeffX2 + gamma.*TheCrossX ;
    dSigmaY = alpha.*ThisHeffY2 + gamma.*TheCrossY ;
    dSigmaZ = alpha.*ThisHeffZ2 + gamma.*TheCrossZ ;
    dSigma = [dSigmaX;dSigmaY;dSigmaZ] ;
    dSigmaRMS = sqrt(sum((dSigma./alpha).^2)/NN) ;
else % alpha is time-dependent
    dSigmaX = alpha(t).*ThisHeffX2 + gamma.*TheCrossX ;
    dSigmaY = alpha(t).*ThisHeffY2 + gamma.*TheCrossY ;
    dSigmaZ = alpha(t).*ThisHeffZ2 + gamma.*TheCrossZ ;
    dSigma = [dSigmaX;dSigmaY;dSigmaZ] ;
    dSigmaRMS = sqrt(sum((dSigma./alpha(t)).^2)/NN) ;
end

%% 

global TheData
% TheData = get(PlotStruct.hF,'userdata') ;
LastT = TheData.LastT ;
LastPlottedT = TheData.LastPlottedT ;
TheData.LastT = t ;
if isfield(TheData,'LastDSigma')
LastDSigma = TheData.LastDSigma ;
UsePrevDSigma = 1 ;
else
UsePrevDSigma = 0 ;    
end
TheData.LastDSigma = dSigma ;
% if UsePrevDSigma
% dSigma = dSigma + .9.*LastDSigma ;
% % disp('y')
% end
TheData.dSigmaRMS = dSigmaRMS ;
% set(PlotStruct.hF,'userdata',TheData) ;

%% Exit ??
if t>.01
   ''   ;
end
if ~((t-LastT)> 0) 
    return
end

if ~isfinite(PlotStruct.DeltaT)
    return
end

%% Update surface and arrows plot

% ThisCData =  reshape(ColorFromHorPsiTheta01(SigmaX,SigmaY,SigmaZ),N,N,N,3) ;
% set(PlotStruct.hS,'cdata',ThisCData) ;
% TheUdata = reshape(SigmaX,N(1),N(2),N(3)) ;
% TheVdata = reshape(SigmaY,N(1),N(2),N(3)) ;
% TheWdata = reshape(SigmaZ,N(1),N(2),N(3)) ;
% disp(num2str(mean(sqrt(TheUdata(:).^2+TheVdata(:).^2+TheWdata(:).^2)))) ;


%% update  hysteresis loop 
if PlotStruct.DrawIt
    if isfield(PlotStruct,'HystDir')
        SsX = -sum(SigmaX)./prod(N) ;
        SsY = -sum(SigmaY)./prod(N) ;
        SsZ = -sum(SigmaZ)./prod(N) ;
        
        SsK = SsX*PlotStruct.HystDir(1) + SsY*PlotStruct.HystDir(2) + SsZ*PlotStruct.HystDir(3) ;
        
        SsK = [get(PlotStruct.hL6x,'ydata'),SsK] ;
        %         SsY = [get(PlotStruct.hL6y,'ydata'),SsY] ;
        %         SsZ = [get(PlotStruct.hL6z,'ydata'),SsZ] ;
        %
        ThisHhK(1) = ThisHhX(1)*PlotStruct.HystDir(1) + ThisHhX(2)*PlotStruct.HystDir(2) + ThisHhX(3)*PlotStruct.HystDir(3) ;
        HsK = [get(PlotStruct.hL6x,'xdata'),ThisHhK(1)./PlotStruct.MaxHn] ;
        set(PlotStruct.hL6x,'xdata',HsK,'ydata',SsK) ;
        set(PlotStruct.hL7x,'xdata',HsK(end),'ydata',SsK(end)) ;
        set(PlotStruct.hA2,'xlim',1.1.*[-1,+1],'ylim',1.1.*[-1,+1]) ;
    else
    if PlotStruct.MaxHn~=0
        SsX = -sum(SigmaX)./prod(N) ;
        SsY = -sum(SigmaY)./prod(N) ;
        SsZ = -sum(SigmaZ)./prod(N) ;
        
        SsX = [get(PlotStruct.hL6x,'ydata'),SsX] ;
        SsY = [get(PlotStruct.hL6y,'ydata'),SsY] ;
        SsZ = [get(PlotStruct.hL6z,'ydata'),SsZ] ;
        
        HsX = [get(PlotStruct.hL6x,'xdata'),ThisHhX(1)./PlotStruct.MaxHn] ;
        HsY = [get(PlotStruct.hL6y,'xdata'),ThisHhY(1)./PlotStruct.MaxHn] ;
        HsZ = [get(PlotStruct.hL6z,'xdata'),ThisHhZ(1)./PlotStruct.MaxHn] ;
    end
    
    if PlotStruct.MaxHn~=0
        set(PlotStruct.hL6x,'xdata',HsX,'ydata',SsX) ;
        set(PlotStruct.hL7x,'xdata',HsX(end),'ydata',SsX(end)) ;
        set(PlotStruct.hL6y,'xdata',HsY,'ydata',SsY) ;
        set(PlotStruct.hL7y,'xdata',HsY(end),'ydata',SsY(end)) ;
        set(PlotStruct.hL6z,'xdata',HsZ,'ydata',SsZ) ;
        set(PlotStruct.hL7z,'xdata',HsZ(end),'ydata',SsZ(end)) ;
        set(PlotStruct.hA2,'xlim',1.1.*[-1,+1],'ylim',1.1.*[-1,+1]) ;
    end
    end
end
%% Exit ?

if ~((t-LastPlottedT)> PlotStruct.DeltaT)
    return
end

disp(['Here ',num2str(t)]) ;
if PlotStruct.DrawIt
    set(PlotStruct.hQ,'udata',reshape(SigmaX,N(1),N(2),N(3)),'vdata',reshape(SigmaY,N(1),N(2),N(3)),'wdata',reshape(SigmaZ,N(1),N(2),N(3))) ;
    %set(PlotStruct.hQ,'udata',SigmaX,'vdata',SigmaY,'wdata',SigmaZ) ;
    set(gcf,'name',[' t = ',num2str(t)]) ;

    % view(t*300,45) ;
    drawnow ;
    TheData.LastPlottedT = t ;
    % set(PlotStruct.hF,'userdata',TheData) ;
end

%% Save gif frame
if PlotStruct.SaveGif

    frame = getframe(gcf);
    im = frame2im(frame);
    [Agif,map] = rgb2ind(im,256);

    if t ~= 0
        imwrite(Agif,map,PlotStruct.GifFilename,'gif','WriteMode','append','DelayTime',1/30);
    end
end

