
%mur is the set relative permeability
%Ms is the saturation magnetization in T
function stateFcn = MakeMH_Fe_const_mur( mur, Ms )

    mu0 = 4 * pi * 1e-7;
    
    nH = 100;
    
    H = 10.^(linspace( 0, log10(Ms*2/mu0), nH ) );
    H(1) = 0;
    %H = linspace( 0, Ms*2/mu0, nH );

    M = (mur-1) * H;
    
    M( M>Ms/mu0 ) = Ms/mu0;
    
    plot(H*mu0,M*mu0,'-o');
    
    out = zeros( nH, 2 );
    out(:,1) = H;
    out(:,2) = M;
    
    dlmwrite( ['Fe_mur_' num2str(mur) '_Ms_' num2str(Ms) '.txt'], out, 'delimiter','\t','precision','%15.7f');
    
    
    stateFcn.nT = int32(3);
    stateFcn.nH = int32(nH);
    stateFcn.T = [200,300,400];
    stateFcn.H = out(:,1);
    stateFcn.M = zeros( stateFcn.nT, stateFcn.nH );
    stateFcn.M(1,:) = out(:,2);
    stateFcn.M(2,:) = out(:,2);
    stateFcn.M(3,:) = out(:,2);

    
end