function makePDFResults(filename)

load([filename 'Problem.mat']);
load([filename 'Solution.mat']);


Ms = median(problem.Ms);
A0 = median(problem.A0);
K0 = median(problem.K0);
Msig = info.Msig;
A0igFactor = info.A0igFactor;
coneAngleDegree = info.coneAngleDegree;
nGrains = mesh_params.nGrains;
offSetD = mesh_params.offSetD;
NumRefinments = mesh_params.NumRefinments;
thisGridL = mesh_params.thisGridL;
NNs = mesh_params.NNs;
H_N = 2*K0/(Ms) ; % [T]
Hvalues = problem.Hext(:,1);
mu0 = 4*pi*1e-7;

fname = [filename];
fnameSave = [filename];

%% Plot
hF = figure('color','w') ;
 ppsz = [2*20,19] ;
    ppps = [0,0,2*20,19] ;
    set(hF,'PaperUnits','centimeters','PaperSize',ppsz,'PaperPosition',ppps) ;

hA1 = subplot(1,2,1) ;

cartesianUnstructuredMeshPlot(problem.grid_pts,problem.grid_abc,GridInfo,mesh_cart.iIn,fname,hA1)  ;
set(hA1,'visible','off') ;

hA2 = subplot(2,2,2) ;



plot(H,M/Ms,'.-k','linewidth',1.5,'markersize',12,'parent',hA2) ;
hold on ; 
hL2 = plot(Hc,0,'hw','markerfacecolor','r','markersize',11,'parent',hA2) ;
text(Hc,0,'   H_C','horizontalalignment','left','verticalalignment','bottom','color','r','parent',hA2,'fontsize',13)

B = (H+4*pi*1e-7*M);
[BH_max,indx] = max(-H/(4*pi*1e-7).*B);
hL3 = plot(H(indx),M(indx)/Ms,'hw','markerfacecolor','b','markersize',11,'parent',hA2) ;
text(H(indx),M(indx)/Ms-0.35,'   (BH)_{max}','horizontalalignment','left','verticalalignment','bottom','color','b','parent',hA2,'fontsize',13)

% hL3 = plot(-H_N,0,'hw','markerfacecolor','b','markersize',9,'parent',hA2) ;
% text(-H_N,0,'H_N   ','horizontalalignment','right','verticalalignment','bottom','color','b','parent',hA2,'fontsize',13) 
set(gca,'ylim',[-1,+1].*1.2,'xlim',[min(Hvalues),max(Hvalues)]) ;
% legend([hL2,hL3],{'Coercive Field','Nucleation Field'},'fontsize',13,'location','northwest')
legend([hL2,hL3],{'Coercive Field','Maximum energy product'},'fontsize',13,'location','SouthEast')
grid on
xlabel('H [T]','fontsize',13)
ylabel('M / M_s [-]','fontsize',13)
set(hA2,'linewidth',1.5,'fontsize',13)
set(gcf,'units','normalized','position',[0.05,0.05,.9,.9]) ;

annotation(hF,'textbox','position',[.6-.02,.39+.05,.3,.1],'string',['Hc = ' num2str(Hc),' [T]'],'LineStyle','none','fontsize',13,'FontWeight','bold','FontAngle','italic') ; 
annotation(hF,'textbox','position',[.573-.02,.35+.05,.3,.1],'string',['(BH)_{max} = ' num2str(BH_max/1000),' [kJm^{-3}]'],'LineStyle','none','fontsize',13,'FontWeight','bold','FontAngle','italic') ; 

annotation(hF,'textbox','position',[.82-.02,.39+.05,.3,.1],'string',[num2str(nGrains),' Grains'],'LineStyle','none','fontsize',13,'FontWeight','bold') ; 

theTextBox1 = annotation(hF,'textbox','position',[.6,.3+.05,.3,.1],'string',[],'LineStyle','none','fontsize',13) ; 
set(theTextBox1,'string',[...
    'M_s^{ } = ',DoTheUnitThing(Ms*mu0,'T'),newline,...
    'A_0^{ } = ',DoTheUnitThing(A0,'J/m^3'),newline,...
    'K_0^{ } = ',DoTheUnitThing(K0,'J/m^3'),newline,...
    '\theta_{cone} = ',[num2str(coneAngleDegree) '^\circ'],newline,...
    ]) ;
if offSetD > 0 
theTextBoxB = annotation(hF,'textbox','position',[.77-.02,.33+.05,.3,.1],'string','Inter-Grain properties','LineStyle','none','fontsize',13,'FontWeight','bold') ; 
theTextBox2 = annotation(hF,'textbox','position',[.75,.3+.05,.3,.1],'string',[],'LineStyle','none','fontsize',13) ; 
set(theTextBox2,'string',[...
    'M_s^* = ',DoTheUnitThing(Msig*mu0,'T'),newline,...
    'A_0^* = ',DoTheUnitThing(A0*A0igFactor,'J/m^3'),newline,...
    'K_0^* = ','0 J/m^3',newline,...
    '\delta_{} = ',DoTheUnitThing(2*offSetD,'m'),newline,...
     ]) ;
end
 theTextBoxC = annotation(hF,'textbox','position',[.6-.02,.11+.05,.3,.1],'string',['Unstructured Cartesian Mesh (',num2str(NumRefinments),' refinements)'],'LineStyle','none','fontsize',13,'FontWeight','bold') ; 
 
 theTextBox3 = annotation(hF,'textbox','position',[.6,.07+.05,.3,.1],'string',[],'LineStyle','none','fontsize',13) ; 
set(theTextBox3,'string',[...
    'L_x = ',DoTheUnitThing(thisGridL(1),'m'),' ',newline,...
    'L_y = ',DoTheUnitThing(thisGridL(2),'m'),' ',newline,...
    'L_z = ',DoTheUnitThing(thisGridL(3),'m'),' ',newline,...
    ]) ;

 theTextBox4 = annotation(hF,'textbox','position',[.75,.07+.05,.3,.1],'string',[],'LineStyle','none','fontsize',13) ; 
set(theTextBox4,'string',[...
    'N_x base = ',num2str(NNs(1)),' ',newline,...
    'N_y base = ',num2str(NNs(2)),' ',newline,...
    'N_z base = ',num2str(NNs(3)),' ',newline,...
    ]) ;

    set(gcf,'PaperUnits','centimeters','PaperSize',ppsz,'PaperPosition',ppps) ;
    set(gcf,'PaperUnits','inches','Renderer','painters') ;
    figure(gcf) ;
    eval(['print -dpdf ',fnameSave,'Result','.pdf']) ;

    drawnow ;

end