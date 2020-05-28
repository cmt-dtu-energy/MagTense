clear all
close all

%load the results from the optimization of ideal remanence distribution
k = load('bulk_optim_nr_3_nz_6_r1_0.02.mat');

%relative standard deviation on the magnitude of Mrem
Mrem_std = 0.01;

%for each configuration find the optimal distribution of Gaussianly
%perturbed magnets
x_pert = cell(length(k.L),length(k.r2));
p2p_pert = zeros(size(x_pert));
Hmean_pert=p2p_pert;
MremDistr = x_pert;
for i=4:4%1:length(k.L)
    for j=4:4%1:length(k.r2)
        %run the optimizer
        [x_pert{i,j},p2p_pert(i,j),Hmean_pert(i,j),MremDistr{i,j}] = ga_optimHB_perturbed_remanence( k.conf{i,j}, Mrem_std );
     
    end
end
save('bulk_optim_nr3_nz_6_gaussian_perturbed_Mrem_std_0.03.mat');