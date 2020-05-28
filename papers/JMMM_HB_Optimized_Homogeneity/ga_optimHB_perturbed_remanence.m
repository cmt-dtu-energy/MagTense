function [x,p2p,Hmean,Mrem_distr] = ga_optimHB_perturbed_remanence( conf, Mrem_std, Mrem_distr, x_init )
mu0 = 4*pi*1e-7;

%standard deviation on remanence
%Mrem_std = 0.01;%relative
Mdir_std = 0;

%no. of radial rings
nr = 3;
%no. of axial rings
nz = 6;

%load model with optimized remanence distribution
%model = load('results/optim_sol/25-Feb-2019_8_29_10n_18_rnd_0087992_input.mat');
%tiles = model.res.tiles;
%conf = model.conf;
%no. of segments
%create the tiles for this particular configuration
tiles = [];
for i=1:length(conf)    
    
    t = getSegmentedHalbach( conf(i), false );  
    tiles = [tiles refineTiles(  moveHalbachTiles( t, conf(i) ), conf(i).res ) ];
end
        
Nseg = conf(1).n_seg;
%find the self-consisten solution and the corresponding N tensor field
[tiles,N] = iterateFixedGeometry( tiles );

%find the N tensor field for the points of interest
[H,N_sample] = getSampleRegionField(tiles);

%initial distribution, empty, will be generated on the first call to
%modifyMremHB
%[tl,cf,Mrem_distr] = modifyMremHB( conf, nr, nz, Mrem_std, Mdir_std, [] );



%no. of variables to optimize over
ntot = nr * Nseg/2 * nz/2;
%boundaries on each are from 1 to ntot!
lb = ones(ntot,1);
ub = ones(ntot,1) * factorial(4);

%objective function
fct = @(x) optim_Mrem_distr(x, conf, N, N_sample,nr, nz, Mrem_distr);

opts = optimoptions(@ga, 'UseParallel',false,...
                    'PopulationSize',200,...
                    'MaxGenerations',100,...                                                      
                    'PlotFcn', {@gaplotselection,@gaplotbestf,@gaplotbestindiv},...                    
                    'Display','iter',...
                    'InitialPopulation',x_init);

IntCon = 1:1:ntot;
                
[x,fval,exitflag,output,population,scores] = ga( fct, ntot, [], [], [],[], lb, ub, [],IntCon, opts );
[p2p,Hmean] = fct( x );
end
%x is expected to have the size ntot,1 and contains the indices of the
%permutations of the remanences of each group (which there are ntot of) of
%nominally identical magnets
function [p2p,Hmean] = optim_Mrem_distr( x, conf, N, N_sample, nr, nz, Mrem_distr )
   
   %make a random permutation of the initial distribution
   tl = modifyMremHB( conf, nr, nz, 0, 0, Mrem_distr, x );
   
   %get the iterated solution 
   tl = iterateFixedGeometry( tl, N );
   
   %get the field
   H = getSampleRegionField(tl,N_sample);
   Hnorm = sqrt(sum( H.^2,2) );
   Hmean = mean(Hnorm);
   %return peak-2-peak value
   p2p = ( max(Hnorm) - min(Hnorm) ) / max(Hnorm);
end