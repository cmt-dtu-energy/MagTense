! ODELibtester.f90
    !
    !
    !
    !
    !
    !
    
    
    
    program ODELibTester
    use ODE_Solvers
    
    implicit none
    
    procedure(dydt_fct), pointer :: fct      
    
    fct => dydt_ML
    
    !kaki: should be "mex'ed"
    call MagTense_ODE( fct, t, y0, t_out, y_out )
    
    
    end program ODELibTester