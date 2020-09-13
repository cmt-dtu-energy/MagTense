function [S_fm,S_pm,M_fm,M_pm,Tarr,Harr,parr] = getBR_StateFunction(Tc,eta)
    Tarr = linspace(304,310,1000);
    Harr = linspace(0,1,10);
    parr = linspace(0,1,4);
    
    material = 'LaFeSi';
    
    [cNonField,sNonField] = getNonFieldSC( material, Tarr );
    
    sNonField = repmat(sNonField,length(Harr),1);
    %%note that the dimensions of the data table should be [length(Y),
    %%length(X), length(Z)] since Matlab has its funny ideas
    M_pm = zeros(length(Harr),length(Tarr),length(parr));
    M_fm = M_pm;
    S_pm = M_pm;
    S_fm = M_pm;
    
    for i=1:length(parr)
        [M_p,S_p, M_f, S_f] = getBR_MF_properties( Tarr, Harr, parr(i), Tc, eta, -1, material );
        M_pm(:,:,i) = M_p';
        M_fm(:,:,i) = M_f';
                
        S_pm(:,:,i) = S_p' + sNonField;
        S_fm(:,:,i) = S_f' + sNonField;
    end
    %%get latent heat and modify the entropy function
     [L_PM,ind_PM,L_FM,ind_FM] = getLatentHeat( Tarr, Harr, parr, S_pm, S_fm );
end


function S = smooth_S_fit( T, S, ind )
%%Makes S continuous and differentiable for all T given the 1st order phase
%%transition at T(ind)

%define entropy on either side of the phase transition
S_low = S(ind);
S_high = S(ind+1);

%define the entropy at the middle of the transition
g = 0.5 * ( S_low + S_high );

T0 = T(ind);

%define steepness scale
d = 10;

%define temperature interval outside which the S curve is assumed linear
dT = 0.5;

%fit linear expressions to the lower and upper parts of the S curve
ft_low = fit( T(T<T(ind)-dT)'-T0, S(T<T(ind)-dT)', 'poly2' );

ft_high = fit( T(T>T(ind)+dT)'-T0, S(T>T(ind)+dT)', 'poly1' );
getFigure();plot(T,S,'--');
plot(T,feval(ft_low,T-T0));
plot(T,feval(ft_high,T-T0));

%assume the expression f(T) = C * asinh( d * (T-T0) ) + g for the curve on
%the interval T-dT_L <= T <= T + dT_H
%Enforcing continuity and differentiability this gives 
%%


fct = @( d2, d4, x ) S_fit_fct( [d2,d4], x, T0, ft_low, ft_high );

fit( T', S', fct, 'StartPoint', [10, 10 ], 'Lower', [0.1,0.1], 'Upper', [1000,1000] );


end

function f = S_fit_fct( x, T, T0, ft_low, ft_high )
%% x is the parameter vector that is being fitted to
%% for T<=T0: f(T) = A(T-T0)^2 + B(T-T0) + C + k1( tanh(k2(T-T0)) - tanh(k2Tmin) )
%% for T>=T0: f(T) = D(T-T0) + E + k3( tanh(k4(T-T0)) - tanh(k4Tmax) )
%% with Tmin = T(1) and Tmax = T(end)
    f = zeros( length(T), 1 );

    A = ft_low.p1;
    B = ft_low.p2;
    C = ft_low.p3;
    D = ft_high.p1;
    E = ft_high.p2;
    k2 = x(1);
    k4 = x(2);
    Tmin = T(1);
    Tmax = T(end);
    
    %% find k3 through enforcing continuity and differentiability at T=T0
    k3 = ( C - E - (D-B)/k2 * tanh(k2*Tmin) ) / ( k4/k2 * tanh(k2*Tmin) - tanh(k4*Tmax) );
    k1 = (D-B+k3*k4)/k2;
    
    %below T0    
    f( T<=T0 ) = A * (T(T<=T0)-T0).^2 + B * (T(T<=T0)-T0) + C + k1 * ( tanh(k2*(T(T<=T0)-T0)) - tanh(k2*Tmin) );
    
    %above T0
    f( T>T0 ) = D * (T(T>T0)-T0) + E + k3 * ( tanh(k4*(T(T>T0)-T0)) - tanh(k4*Tmin) );
     
    
end


function S = smooth_S_fct(T,S,ind)
%%T is the temperature array
%%S the entropy 
%%ind the index in T at which the discrete transition occurs
    %temperature to move on either side of the discontinuity
    dT = 1;
    if ind>1
        %scale of tanh in order to get it close to -1 and +1 (the steepness of
        %tanh)
        sclLow = -1;
        sclHigh = 2;
        %low side fix point
        Tlow = T(ind)-dT;
        indLow = find(T-Tlow>=0,1);
        Tlow = T(indLow);
        Slow = S(indLow);
        %high side fix point
        Thigh = T(ind)+dT;
        indHigh = find(T-Thigh>=0,1);
        Thigh = T(indHigh);
        Shigh = S(indHigh);

        fct = @(x) asinh(x);
        figure; hold on;grid;
        plot(T,S);
        S(indLow:indHigh) = ( fct( (T(indLow:indHigh)-Tlow)./(Thigh-Tlow) * (sclHigh-sclLow)+sclLow ) - fct(sclLow) ) ./(fct(sclHigh)-fct(sclLow)) .* ( Shigh - Slow ) + Slow;
        plot(T,S,'-o');
    end
end