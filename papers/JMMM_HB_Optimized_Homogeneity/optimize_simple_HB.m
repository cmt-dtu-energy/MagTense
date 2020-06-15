function res = optimize_simple_HB()
addpath('util/');

%bounds on possible values of Mrem
Mrem_bnds = [1.0, 1.5];

%optimization orderInd. 0 for p2p, else a vector of the orders to optimize
%for
orderInd = 0;

%setup the default settings
setts = struct();
%inner radius, always fixed
setts.r1 = 0.0175;
setts.doIteration = true;
setts.Mnorm_std = 0;
setts.Mdir_std = 0;
setts.includePoleShoes = false;
setts.dir = ['results/' getDateTimeString() '_HB_simple_optim_2.0T/']; 
gr = getMagnetGrade('N40UH');
%loop over number of segments and length to optimize the outer radius for
%specific field requirement
Nseg = 8:1:16;
L = linspace( 0.08, 0.3, 20 );
ns = length(Nseg);
nl = length(L);
res = cell( ns, nl );%

parfor i=1:ns
    
    for j=1:nl
        st = setts;
        st.HBLength = L(j);
        st.Nseg = Nseg(i);
        st.r2 = 0.15;
        st.MagGrade = cell( Nseg(i), 1 );
        for k=1:length(st.MagGrade)
            st.MagGrade{k} = gr;
        end
        fct = @(x) runModel( x, st, orderInd );
        options = optimset('Display','iter','PlotFcns', {@optimplotx,@optimplotfval},'FunctionTolerance',1e-4);
        %minimize the mass of the magnet given a minimum average field of
        %2.0 T
        %[x,fval,flag] = fminbnd( fct, st.r1*1.1, 0.5, options );
        
        %lower bounds
        lb = [st.r1*1.1];%,ones(1,Nseg(i))*Mrem_bnds(1)];
        %upper bounds
        ub = [0.5];%,ones(1,Nseg(i))*Mrem_bnds(2)];
        [x,fval,flag] = fmincon( fct, st.r1*2, [], [], [], [], lb, ub, [], options );
        %retrieve the essential parameters
        st.r2 = x;
        res(i,j) = run_singleHB_OI( st );
    end    
end
end
function val = runModel( x, setts, orderInd )
    if length(x) > 1
        %setup the problem for also optimizing on orders and allowing for
        %mag grade variation
        %outer radius
        setts.r2 = x(1);
        for i=1:length(x(2:end))
            setts.MagGrade{i}.Mrem = x(i);
        end
        res = run_singleHB_OI( setts );
        
        if orderInd(1) == 0
            %optimize for p2p only
            B = res.Bnorm;
            p2p = (max(B)-min(B))/max(B);
            goal = 1e-4;
            
            if p2p > goal
                w = p2p / goal;
            else
                w = 1;
            end
        else
            %optimize for the given orders
        end
        
        val = w * abs( mean(res.Bnorm) - 2.0 );
    else
        setts.r2 = x;
        res = run_singleHB_OI( setts );
        val = abs( mean(res.Bnorm) - 2.0 );
    end
end
