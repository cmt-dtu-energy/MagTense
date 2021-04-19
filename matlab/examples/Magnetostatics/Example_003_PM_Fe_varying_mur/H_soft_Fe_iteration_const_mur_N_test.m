clear all
close all

%demag const
N = -1./3.;
%const perm.
mur = 2;
%Applied field (in a particular direction)
Happ = 1;%A/m

%no. of iterations
n_ite = 100;
%output field
Hy_a = Happ ./ (1. - N .* ( mur - 1 ) );

n = 100;
Hy = zeros(n,1);
Hy(1) = 0;

for i=1:n-1

    Hy(i+1) = Happ + N * (mur-1) * Hy(i);
end

getFigure();
plot(Hy,'-o','linewidth',2);
plot([0,n],[Hy_a Hy_a],'--','linewidth',2)
%ylim([-1,4])