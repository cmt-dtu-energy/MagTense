
clear all
close all
addpath('util/');
mu0 = 4*pi*1e-7;

%standard deviation on remanence
Mrem_std = 0.05;%0.01;%relative
Mdir_std = 0;
%no. of samples to do
n = 50;
%no. of radial rings
nr = 3;
%no. of axial rings
nz = 6;
%load model with optimized remanence distribution

load('bulk_optim_perturbed_and_ideal.mat');
for ii=1:length(k.conf(:,1))
    for jj=1:length(k.conf(1,:))
        conf = k.conf{ii,jj};
        Nseg = conf(1).n_seg;
        tiles = [];
        for i=1:length(conf)    
            t = getSegmentedHalbach( conf(i), false );  
            tiles = [tiles refineTiles(  moveHalbachTiles( t, conf(i) ), conf(i).res ) ];
        end

        %find the self-consisten solution and the corresponding N tensor field
        [tiles,N] = iterateFixedGeometry( tiles );

        %find the N tensor field for the points of interest
        [H,N_sample] = getSampleRegionField(tiles);
        Hnorm = sqrt(sum( H.^2,2) );
        p2p = zeros(n,1);

        p2p(1) = ( max(Hnorm) - min(Hnorm) ) / max(Hnorm);

        %initial distribution, empty, will be generated on the first call to
        %modifyMremHB
        [tl,cf,Mrem_distr] = modifyMremHB( conf, nr, nz, Mrem_std, Mdir_std, [] );

        %loop over each sample
        cf = cell(n,1);
        Mrem_var = cf;
        x_opt = cf;
        p2p_opt = cf;
        Hmean_opt = cf;
        Mrem_distr_opt = cf;
        parfor i=2:n
           %tl = tiles;
           %set the remanence distribution of the current tiles 

           %make a random permutation of the initial distribution
           %x_init is the random distribution index to be used the first time

           x_init = randi(factorial(4),1,nr*nz/2*Nseg/2);
           %also permute the mother-distribution for the overall problem
           k = randperm(numel(Mrem_distr));
           Mrem_distr_perm = reshape( Mrem_distr(k), [4,Nseg*nr*nz/4] )';
           [tl,cf{i},Mrem_var{i}] = modifyMremHB( conf, nr, nz, Mrem_std, Mdir_std, Mrem_distr_perm, x_init );

           %get the iterated solution 
           tl = iterateFixedGeometry( tl, N );

           %get the field
           H = getSampleRegionField(tl,N_sample);
           Hnorm = sqrt(sum( H.^2,2) );
           %and the p2p
           p2p(i) = ( max(Hnorm) - min(Hnorm) ) / max(Hnorm); 

           %get the optimized solution
           [x_opt{i},p2p_opt{i},Hmean_opt{i},Mrem_distr_opt{i}] = ga_optimHB_perturbed_remanence( conf, Mrem_std, Mrem_distr, x_init );


        end
        save(['bulk_HB_optim_perturbed_re_optim_i_j_' num2str(ii) '_' num2str(jj) '_0.05.mat']);
    end
end