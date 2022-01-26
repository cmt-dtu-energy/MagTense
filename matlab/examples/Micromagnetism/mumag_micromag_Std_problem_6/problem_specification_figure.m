figure
numps=20;xs=linspace(0,80,numps);ys=linspace(0,1,1);
[X,Y]=meshgrid(xs,ys);
U=2*((X<40)-0.5);V=Y*0+0.3;
h=fill([0,40,40,0],[0,0,2.3,2.3],[0.3,0.3,1],[40,80,80,40],[0,0,2.3,2.3],[1,0.3,0.3]);
set(h,'facealpha',0.5);
hold on
quiver(X,Y,U,V,0.1,'k')
axis equal
text(10,6,'$$A^L$$, $$K^L$$, $$|{\bf{M}}|^L$$','color',[0.3,0.3,1],'fontsize',14,'Interpreter','latex');
text(50,6,'$$A^R$$, $$K^R$$, $$|{\bf{M}}|^R$$','color',[1,0.3,0.3],'fontsize',14,'Interpreter','latex');
an=annotation('doublearrow');an2=annotation('doublearrow');
an.X=[0.13,0.775/2+0.13];an.Y=[0.48,0.48];
an2.X=[0.775/2+0.13,0.775+0.13];an2.Y=[0.48,0.48];
text(15,-6,'$$40$$ nm','fontsize',14,'Interpreter','latex');
text(55,-6,'$$40$$ nm','fontsize',14,'Interpreter','latex');
ylim([-8,9])
xlim([0,80])
pbaspect([40,8.5,1])
an.X=[0.13,0.775/2+0.13];an.Y=[0.412,0.412];
an2.X=[0.775/2+0.13,0.775+0.13];an2.Y=[0.412,0.412];
axis off