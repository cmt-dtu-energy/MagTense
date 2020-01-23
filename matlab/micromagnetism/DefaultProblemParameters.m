function DefaultSetupStruct = DefaultProblemParameters()
    MU0 = pi*4e-7 ; % Vacuum permeability

    %% Grid Parameters
%     nHalfTimes = [0,1,2] ;  % from coarsest (1) to finest (K)

    %% "H" : External field parameters
    MaxHx = 0 ; % max x component of external field vector [T]
    MaxHy = 0 ; % max y component of external field vector [T]
    MaxHz = 0 ; % max z component of external field vector [T]
    FreqH = 1*(.25+1/2).*pi ;  % oscillation frequency of external field
    Ms = 800e3 ; % Saturation magnetism % Vacuum permeability

    %% "J" : Exchange parameters
    A0 = 13e-12 ; % strength of exchange term
    
    %% "M" : Demagnetization parameters
%     dz = .1 ;

    %% "K" : Anisotropy parameters
    Nmax = 1 ; % number of "grains"     %--- Remark: This needs to be determined in the program 
    InPlane = 1 ;
    K0 = 0 ; % strength of anisotropy term
    Kdir = [0 0 0];
    Kdir = Kdir/norm(Kdir); % anisotropy direction specified by 3D vector
    Kx = Kdir(1);
    Ky = Kdir(2);
    Kz = Kdir(3);
    Kx = 0 ; % x component of easy axis vector
    Ky = 0 ; % y component of easy axis vector
    Kz = 0 ; % z component of easy axis vector
    
    %% Time-evolution parameters
    MaxT = 5 ;
    nT = 1000 ;
    gamma = 0 ; % strength of precession term
    MaxT0 = 1 ; 
    alpha = @(t) -10000*(10.^(5*min(t,MaxT0)/MaxT0));  % strength of damping term
    
    %% Sigma initial state
    InitialState = 'rand' ;
    % InitialState = 'magY' ;

    %% Figure Options
    DrawIt = 0 ;
    DeltaT = 1/500  ;
    SaveGif = 0 ;

    %% Applied Field time-dependence (function handles)
    HsX = @(t) MaxHx.*sin(FreqH.*t) ;
    HsY = @(t) MaxHy.*cos(FreqH.*t) ;
    HsZ = @(t) MaxHz ;
    
    %% Dimensions
    Lx = 2 ;
    Ly = 2 ;
    Lz = .5 ;
    
    %--- Define all variables in a struct    
    DefaultSetupStruct.MU0 = MU0;
    DefaultSetupStruct.FreqH = FreqH;
    DefaultSetupStruct.Nmax = Nmax;
    DefaultSetupStruct.InPlane = InPlane;
    DefaultSetupStruct.K0 = K0;
    DefaultSetupStruct.Kdir = Kdir;
    DefaultSetupStruct.nT = nT;
    DefaultSetupStruct.gamma = gamma;
    DefaultSetupStruct.alpha = alpha;
    DefaultSetupStruct.DrawIt = DrawIt ;
    DefaultSetupStruct.DeltaT = DeltaT  ;
    DefaultSetupStruct.SaveGif = SaveGif ;
    DefaultSetupStruct.MaxT = MaxT;
    DefaultSetupStruct.K0 = K0;
    DefaultSetupStruct.Kx = Kx;
    DefaultSetupStruct.Ky = Ky;
    DefaultSetupStruct.Kz = Kz;
    DefaultSetupStruct.A0 = A0;
    DefaultSetupStruct.Ms = Ms;
    DefaultSetupStruct.MaxHx = MaxHx;
    DefaultSetupStruct.MaxHy = MaxHy;
    DefaultSetupStruct.MaxHz = MaxHz;
    DefaultSetupStruct.HsX = HsX;
    DefaultSetupStruct.HsY = HsY;
    DefaultSetupStruct.HsZ = HsZ;
%     DefaultSetupStruct.nHalfTimes = nHalfTimes;
    DefaultSetupStruct.InitialState = InitialState;
    DefaultSetupStruct.Lx = Lx;
    DefaultSetupStruct.Ly = Ly;
    DefaultSetupStruct.Lz = Lz;
    DefaultSetupStruct.UseImplicitSolver = 0;
    DefaultSetupStruct.UseExplicitSolver = 0;
    DefaultSetupStruct.UseDynamicSolver = 0 ;
    DefaultSetupStruct.UseImplicitStepsSolver = 0 ;
    DefaultSetupStruct.DemagTensorFileName = 0;
    DefaultSetupStruct.SigmaInitialFileName = 0;
    DefaultSetupStruct.CalcEigenvalue = 0;
    DefaultSetupStruct.AlreadyEquilibrium = 0;
    DefaultSetupStruct.SaveTheResult = 1 ;
    DefaultSetupStruct.MaxComputationalTimePerStep = inf ;
    DefaultSetupStruct.use_sparse = 0;
    DefaultSetupStruct.use_gpuArray = 0;
    DefaultSetupStruct.threshold = 0; 
    DefaultSetupStruct.thresholdFract = 0;
    DefaultSetupStruct.use_single = 0;
    % FFT
    DefaultSetupStruct.FFTdims = [0,0,0] ;

    
end