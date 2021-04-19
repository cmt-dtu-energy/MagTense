function ddSigma = CalculateDeltaSigmaEquilibriumExplicitOut(t,Sigma,HHess,AA,MaxComputationalTimePerStep)
ThisHHess5ttBlock = -HHess(t,Sigma) ;
try 
    ThisTime = toc ;
catch
    ThisTime = 0 ;
end
if  ThisTime>MaxComputationalTimePerStep
    error('Maximum Computational Time Exceeded') ; '' ;
    
%    break ; 
end
%%
Sigma = NormalizeSigmaCAT(Sigma) ;
NN = round(numel(Sigma)/3) ;
% N = round(NN^(1/3)) ;
% Sigma = Sigma./NormalizeSigmaCAT(Sigma) ;
SigmaX = Sigma(0*NN+[1:NN]) ;
SigmaY = Sigma(1*NN+[1:NN]) ;
SigmaZ = Sigma(2*NN+[1:NN]) ;
% ddSigma = -MtimesMtimes(Sigma,ddSigma) ;

% ddSigmaX = ddSigma(0*NN+[1:NN]) ;
% ddSigmaY = ddSigma(1*NN+[1:NN]) ;
% ddSigmaZ = ddSigma(2*NN+[1:NN]) ;
%%

eZ = [zeros(2*NN,1);ones(NN,1)] ;
e1 = Sigma ;
e2 = CrossProdCAT(e1,eZ) ; 
e2 = NormalizeSigmaCAT(e2) ;
e3 = CrossProdCAT(e1,e2) ;

% figure ; plot(e2,'.') ; hold on ; plot(Sigma,'.r') ; hold on ; plot(e3,'.k') ;
% DotProdCAT(e2,e3)
% TT = [[diag(e1(0*NN+[1:NN]));diag(e1(1*NN+[1:NN]));diag(e1(2*NN+[1:NN]))],...
%       [diag(e2(0*NN+[1:NN]));diag(e2(1*NN+[1:NN]));diag(e2(2*NN+[1:NN]))],...
%       [diag(e3(0*NN+[1:NN]));diag(e3(1*NN+[1:NN]));diag(e3(2*NN+[1:NN]))]];
% TT = sparse(TT) ;
TT = spdiags([[e1(2*NN+(1:NN));zeros(2*NN,1)],...
              [e1(1*NN+(1:NN));e2(2*NN+(1:NN));zeros(NN,1)],...
              [e1(0*NN+(1:NN));e2(1*NN+(1:NN));e3(2*NN+(1:NN))],...
              [zeros(NN,1);e2(0*NN+(1:NN));e3(1*NN+(1:NN))],...
              [zeros(2*NN,1);e3(0*NN+(1:NN))]],[-2*NN,-NN,0,NN,2*NN],3*NN,3*NN);
  '' ;
%%
%Am
% [ThisHeffX,ThisHeffY,ThisHeffZ] = CalculateHeffTerm(t,SigmaX,SigmaY,SigmaZ,AvrgMatrix,CopyMatrix,Mfact,KcoarXX,KcoarXY,KcoarXZ,KcoarYX,KcoarYY,KcoarYZ,KcoarZX,KcoarZY,KcoarZZ,AA) ;
% [ThisHeffX,ThisHeffY,ThisHeffZ] = CalculateHeffTerm(0,SigmaX,SigmaY,SigmaZ,AvrgMatrix,CopyMatrix,Mfact,KcoarXX,KcoarXY,KcoarXZ,KcoarYX,KcoarYY,KcoarYZ,KcoarZX,KcoarZY,KcoarZZ,AA) ;

%Adm
% [ddThisHeffX,ddThisHeffY,ddThisHeffZ] = CalculateHeffTermWithoutHappl(t,ddSigmaX,ddSigmaY,ddSigmaZ,AvrgMatrix,CopyMatrix,Mfact,KcoarXX,KcoarXY,KcoarXZ,KcoarYX,KcoarYY,KcoarYZ,KcoarZX,KcoarZY,KcoarZZ,AA) ;

%dAm
ddThisHhX = AA.ddHhX(SigmaX,SigmaY,SigmaZ,t) ; 
ddThisHhY = AA.ddHhY(SigmaX,SigmaY,SigmaZ,t) ;
ddThisHhZ = AA.ddHhZ(SigmaX,SigmaY,SigmaZ,t) ;

% Hf
% [HfX,HfY,HfZ] = CalculateHeffTermWithoutHappl(t,SigmaX,SigmaY,SigmaZ,AvrgMatrix,CopyMatrix,Mfact,KcoarXX,KcoarXY,KcoarXZ,KcoarYX,KcoarYY,KcoarYZ,KcoarZX,KcoarZY,KcoarZZ,AA) ;
% Ha
% HaX = AA.HhX(SigmaX,SigmaY,SigmaZ,t) ;
% HaY = AA.HhY(SigmaX,SigmaY,SigmaZ,t) ;
% HaZ = AA.HhZ(SigmaX,SigmaY,SigmaZ,t) ;


mx = SigmaX ;
my = SigmaY ;
mz = SigmaZ ;

[bx,by,bz] = FixedTermFstarComponents(mx,my,mz,ddThisHhX,ddThisHhY,ddThisHhZ) ;
b = -[bx;by;bz] ;
% ddRes = ThisHHess5*ddSigma + [bx;by;bz] ;
%%
% ThisHHess5tt = TT.'*ThisHHess5*TT ;
ThisHHess5ttBlock = BasisChange(ThisHHess5ttBlock,TT) ;
btt = TT.'*b ;
% ThisHHess5ttBlock = ThisHHess5tt(NN+1:end,NN+1:end) ;
ThisHHess5ttBlock = ThisHHess5ttBlock(NN+1:end,NN+1:end) ;
% theeigs = eigs(ThisHHess5ttBlock,1,'smallestreal') ;
% theeigs = eigs(ThisHHess5ttBlock,1,'largestreal') ;

xtt = ThisHHess5ttBlock\btt(NN+1:end) ;
x = TT(:,NN+1:end)*xtt ; 
% sum(abs(ThisHHess5*x - b))
% ddSigma = ThisHHess5\b ;
ddSigma = x ;
%%

% ddSigma = -MtimesMtimes(Sigma,ddSigma) ;
% ddSigma = linsolve(ThisHHess5,b) ;

% 
% dd1 = CrossCrossFunction00(ddSigmaX,ddSigmaY,ddSigmaZ,mx,my,mz,ThisHeffX,ThisHeffY,ThisHeffZ) ;
% dd2 = CrossCrossFunction00(mx,my,mz,ddSigmaX,ddSigmaY,ddSigmaZ,ThisHeffX,ThisHeffY,ThisHeffZ) ;
% dd3 = CrossCrossFunction00(mx,my,mz,mx,my,mz,ddThisHeffX,ddThisHeffY,ddThisHeffZ) ;
% dd4 = CrossCrossFunction00(mx,my,mz,mx,my,mz,ddThisHhX,ddThisHhY,ddThisHhZ) ;
% 
% [ddx,ddy,ddz] = TestImplicitDerivativeSymbolicComponents(SigmaX,SigmaY,SigmaZ,...
%     HaX,HaY,HaZ,...
%     HfX,HfY,HfZ,...
%     ddThisHeffX,ddThisHeffY,ddThisHeffZ,...
%     ddThisHhX,ddThisHhY,ddThisHhZ,...
%     ddSigmaX,ddSigmaY,ddSigmaZ) ;


% [ddx,ddy,ddz] = ThisTestDD2(ddThisHeffX,ddThisHeffY,ddThisHeffZ,...
%     ThisHeffX,ThisHeffY,ThisHeffZ,...
%     ddThisHhX,ddThisHhY,ddThisHhZ,...
%     ddSigmaX,ddSigmaY,ddSigmaZ,...
%     SigmaX,SigmaY,SigmaZ) ;

% ddRes = [ddx;ddy;ddz] ;
% '' ; 
%     figure
%     plot(ddRes,'ok')
%     hold on
%     plot(dd1,'.r')
%     plot(dd2,'.m')
%     plot(dd3,'.b')
%     plot(dd4,'.g')
%     plot(dd1+dd2+dd3+dd4,'^k')    
'' ;
% if theeigs > 0
%     '' ;
% end
[~,p] = chol(-ThisHHess5ttBlock) ;
pp = p ;
LastT = t ;

global TheData
%     TheData = get(PlotStruct.hF,'userdata') ;
if isfield(TheData,'p')
    if  TheData.p > 0
        LastT = TheData.LastT ;
        p = TheData.p ;
    end
end

TheData.p = p ;
TheData.LastT = LastT ;
% if t>LastT
%     ddSigma = 0.*ddSigma ;
% end
% set(PlotStruct.hF,'userdata',TheData)
% disp([num2str(t),'  ',num2str(log10(-theeigs)),' ',num2str(p)])
% plot(t,theeigs,'.k') ; drawnow ;
disp([num2str(t),'  ',num2str(p),'  ',num2str(pp)])

return ;

