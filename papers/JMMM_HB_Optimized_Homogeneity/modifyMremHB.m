%%modifies the Mrem for each magnet segment in a Halbach configuration. It
%%is assumed that conf is a struct with n HB ring configurations where n =
%%nr * nz. Each configuration is assumed to consist of the same no. of
%%segments as specified in conf(i).Nseg (should be equal for all
%%i=1..length(conf) ).
%%Mnorm_std is the desired standard deviation on the norm of Mrem (relative) while
%%Mdir_std is the corresponding std dev on the direction of Mrem (in deg.)
%%distrIndex, if present, is the index from 1 to n! where n=nz*2 that
%%defines the index of the permutation to use for a group of nominally
%%identical magnets. The permutation determines how these identical magnets
%%have their individual remanence perturbed by the distribution.
%%Example: if nz = 6 where are 12! = 480e6 permutations of the numbers
%%1,2..12 and each of these numbers are interpreted as an index into the
%%pre-defined distribution. So, given a value of the distrIndex the order
%%of the numbers is found as, e.g., 1,2,4,3,5,6,7,8,9,10,11,12 etc. 
%%distrIndex is expected to have the size (nr*nseg/2,1) indicating that
%%there are nr*nseg/2 groups of identical magnets each with 2nz members
function [tiles,conf_out,Mrem_var] = modifyMremHB( conf_in, nr, nz, Mnorm_std, Mdir_std, Mrem_var, distrIndex )

    %%it is assumed that the no. of segments is even and the HB is
    %%configured so that the segments are pairwise identical following that
    %%segment i is identical to segment Nseg, i+1 to Nseg-i etc. In the
    %%case of 16 segments one can see that segment 1 and seg. 16 are
    %%identical if they are rotated about the radial direction, thus they
    %%are nominally identical magnets

    %no. of segments
    Nseg = conf_in(1).n_seg;

    conf_out = conf_in;
    
    %three modes possible:
    %1st is random drawing for each magnet segment
    %2nd is a pre-defined distribution where a random permutation is chosen
    %from this
    %3rd is a predetermined list of remanence values for each group and
    %then a permutation is chosen. This mode is not implemented yet
    
    if ~exist('Mrem_var','var') || isempty(Mrem_var)
        %draw a random distribution
        Mrem_var = randn(Nseg*nr*nz/4,4) * Mnorm_std + 1;
        
    elseif ~exist('distrIndex','var')
        %permute an existing distribution
        k = randperm(numel(Mrem_var));
        Mrem_var = reshape( Mrem_var(k), [4,Nseg*nr*nz/4] )';
        
    else        
        x = distrIndex;
        %draw a specific permutation from the pre-made distribution based
        %on the index values given in distrIndex (size(nr*nseg/2,1)).            
        for i=1:nz/2
            for j=1:nr
                for k=1:Nseg/2
                    %get the current permutation (there are four possible sites
                    %for each magnet piece)
                    perm = getPermutation( distrIndex( (i-1)*nr*Nseg/2 + (j-1)*Nseg/2+k), 4 );
                    %take out a vector of the current 4 nominally identical
                    %magnets
                    Mrem_nom_ident = Mrem_var( (i-1)*nr*Nseg/2 + (j-1)*Nseg/2+k, : )';

                    Mrem_var( (i-1)*nr*Nseg/2 + (j-1)*Nseg/2+k , : ) = Mrem_nom_ident( perm ) ;                    
                end            
            end
        end
    end
    
   
    
    %go through each layer axially
    
    for i=1:nz/2
        %through each radial ring
        for j=1:nr
            %index into conf array, which is assumed to be organized so
            %that conf(1..nr) provide the first nr radial HB rings and
            %conf(nr+1..2nr) the next etc.
            indC = 2*(i-1)*nr + 2*(j-1)+1;
            %through each segment pair
            for k=1:Nseg/2
                conf_out(indC).M0(k) = conf_in(indC).M0(k) * Mrem_var( (i-1)*nr*Nseg/2 + (j-1)*Nseg/2+k ,1 );
                conf_out(indC).M0(Nseg-k+1) = conf_in(indC).M0(Nseg-k+1) * Mrem_var( (i-1)*nr*Nseg/2 + (j-1)*Nseg/2+k , 2 );
                
                conf_out(indC+1).M0(k) = conf_in(indC+1).M0(k) * Mrem_var( (i-1)*nr*Nseg/2 + (j-1)*Nseg/2+k ,3 );
                conf_out(indC+1).M0(Nseg-k+1) = conf_in(indC+1).M0(Nseg-k+1) * Mrem_var( (i-1)*nr*Nseg/2 + (j-1)*Nseg/2+k , 4 );
                
            end
        end
    end
    %create the Halbach rings
    tiles = [];
    for i=1:length(conf_out)    
        t = getSegmentedHalbach( conf_out(i), false );  
        tiles = [tiles refineTiles(  moveHalbachTiles( t, conf_out(i) ), conf_out(i).res ) ];
    end
end