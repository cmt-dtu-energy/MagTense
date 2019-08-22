function status = MyOutputFunction(t,SigmaSol,flag,dSigmaMagTol)
if isempty(flag) 
% ThisSigma = SigmaSol(:,end) ;
% NN = numel(ThisSigma)/3 ;
% N = sqrt(NN) ;
% SigmaX = ThisSigma(0*NN+[1:NN]) ;
% SigmaY = ThisSigma(1*NN+[1:NN]) ;
% SigmaZ = ThisSigma(2*NN+[1:NN]) ;
% 
% 
% ThisCData =  reshape(ColorFromHorPsiTheta01(SigmaX,SigmaY,SigmaZ),N,N,3) ;
% set(hS,'cdata',ThisCData) ;
% drawnow ;
% 
try
    global TheData
%     TheData = get(PlotStruct.hF,'userdata') ;
LastRMS = TheData.dSigmaRMS ;
catch
    LastRMS = inf ;
end
% LastRMS = LastRMS(end) ;
disp(['Log10dS ',num2str(log10(LastRMS))]) ;
if (log10(LastRMS) < dSigmaMagTol)
%     disp('Tolerance reached') ;
    status = 1 ; % was 1
else
    status = 0 ; 
end
end