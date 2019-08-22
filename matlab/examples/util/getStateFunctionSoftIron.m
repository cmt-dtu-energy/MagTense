%%Name is not used currently, but enables the possibility for having
%%different state functions
function stateFcn = getStateFunctionSoftIron(name)

    stateFcn = getDefaultMagStatStateFunction();
    %%setup the state function
    
    k = load('MH_Comsol_low_carb_annealed_extrap.txt');            
    
    stateFcn.nT = int32(3);
    stateFcn.nH = int32(length(k(:,1)));
    stateFcn.T = [200,300,400];
    stateFcn.H = k(:,1);
    stateFcn.M = zeros( stateFcn.nT, stateFcn.nH );
    stateFcn.M(1,:) = k(:,2);
    stateFcn.M(2,:) = k(:,2);
    stateFcn.M(3,:) = k(:,2);

end