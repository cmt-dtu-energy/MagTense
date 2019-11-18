close all
clear vars

figure1= figure('PaperPositionMode', 'auto'); fig1 = axes('Parent',figure1,'Layer','top','FontSize',12); grid on; hold all; box on
figure2= figure('PaperPositionMode', 'auto'); fig2 = axes('Parent',figure2,'Layer','top','FontSize',12); grid on; hold all; box on

col_map = lines(3);

%legend for figure 2
for kk = 1:3
    plot(fig2,-10,-10,'-','linewidth',2,'color',col_map(kk,:));
end
plot(fig2,-10,-10,'ko','linewidth',1,'markersize',6);

mu0 = 4*pi*1e-7;
%%test the tetrahedron implementation
%define four vertices that are not co-planar or colinear
v = [[2.5,3,1];[2,1,4];[1.5,4,3];[4.5,5,2]]';

%define center where the field is plotted along the axis
c =  [3,3,2.5];

%define the point of interest (where the field is evaluated)
x = linspace(-10,10,10000);

%define the remanent magnetization
M = [0.324264068, 0.734846928, 0.891545179]';

nms = {'x','y','z'};
for kk=1:3
    tic
    H = zeros(3,length(x));
    for i=1:length(x)
        N = zeros(3,3);
        r = c';
        r(kk) = x(i);
        for j=1:4
            k = circshift( v, j-1, 2 );
            [tmp,P] = getN_Triangle( k, r );
            N = N + tmp;        
        end
        H(:,i) = N * M;
    end
    toc
    %load comsol solution
    cms_path = '..\..\..\documentation\examples_FEM_validation\Validation_tetrahedron\';
    H_comsol = load( [cms_path 'Validation_tetrahedron_normH_' char(nms(kk)) '.txt'] );

    cla(fig1);
    plot(fig1,H_comsol(:,1),H_comsol(:,2),'o','displayname','comsol','linewidth',2,'markersize',10);
    plot(fig1,x,sqrt(sum(H.^2,1)),'displayname','Magtense','linewidth',2,'color','k');
    xlabel(fig1,[ char(nms(kk)) ' position  []']);
    ylabel(fig1,'Field norm [T]');
    outfile = ['comsol_magtense_tetrahedron_comp_dir_ ' char(nms(kk))];
%     print('-dpng',[outfile '.png'] );
%     saveas(fig1,[outfile '.fig']);
    
    plot(fig2,H_comsol(:,1),H_comsol(:,2),'ko','linewidth',1,'markersize',6);
    plot(fig2,x,sqrt(sum(H.^2,1)),'-','linewidth',2,'color',col_map(kk,:));
    xlabel(fig2, 'x, y or z [mm]');
    ylabel(fig2,' |H| [A m^{-1}]');
end
xlim(fig2,[0 6])
legend(fig2,'MagTense along x','MagTense along y','MagTense along z','Comsol','Location','NorthWest')
figure(figure2)
print('-depsc','Tetrahedron_Comsol.eps' );
