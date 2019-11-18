function status = MyOutputFunctionDirect(t,SigmaSol,flag)
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
% try
global TheData
%     TheData = get(PlotStruct.hF,'userdata') ;
LastP = TheData.p ;
% catch
%     LastP = 0 ;
% end
% LastRMS = LastRMS(end) ;
disp([num2str(t(end)),'  Last P : ',num2str(LastP)])
% disp(['Log10dS ',num2str(log10(LastRMS))]) ;
if (LastP > 0)
%     disp('Tolerance reached') ;
    status = 0 ; % was 1
else
    status = 0 ; 
end
end