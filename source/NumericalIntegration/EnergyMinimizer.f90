!>-----------------------------------------------------------------------------------------------
!> @author Rasmus Bjoerk, rabj@dtu.dk, DTU, 2026
!> @brief
!> Energy minimizer for a set of unit vectors: steepest descent on the unit sphere with
!> Barzilai-Borwein step lengths, following Exl et al., J. Appl. Phys. 115, 17D118 (2014).
!>
!> The module knows nothing about micromagnetics. The caller provides a routine that returns the
!> effective field H for a given state m, and one that returns the energy terms of the state the
!> field was last evaluated at. The state has the layout used everywhere in MagTense:
!> mx = m(1:n), my = m(n+1:2n), mz = m(2n+1:3n), and H is laid out the same way.
!>
!> The descent direction of vector i is m_i x (m_i x H_i), the same vector that drives the damping
!> term of the Landau-Lifshitz equation, and the update
!>     m1 = ( (1 - tau^2 |t|^2/4) m0 - tau m0 x t ) / (1 + tau^2 |t|^2/4),   t = m0 x H,
!> is the exact solution of the midpoint rule m1 = m0 - (tau/2) (m0 + m1) x t, so that |m| = 1 is
!> preserved to rounding. One field evaluation per iteration.
!>
!> Convergence is declared when max_i |m_i x H_i| / H_scale < tol. Two safeguards keep the
!> non-monotone Barzilai-Borwein iteration in check: no vector rotates by more than maxrot in one
!> iteration, and an energy watchdog halves the step when the energy would rise above the largest
!> of the last n_hist iterates (Grippo-Lampariello-Lucidi style). If the iteration cap is hit, or the
!> watchdog stalls, the caller's ODE integration is run once from the current state (when fallback is
!> set) and the minimizer is restarted from its result - this is what carries a hysteresis loop
!> through a switching event, where the state sits near a saddle.
!>
!> A vanishing torque is also what a saddle point looks like, and a symmetric starting state sits on
!> one. When saddle_check is set, a converged state is therefore nudged by a small random rotation of
!> every vector and relaxed again; a minimum takes the nudge back and is returned unperturbed, a
!> saddle lets the state fall to a lower energy, which is then checked the same way.
!>-----------------------------------------------------------------------------------------------
module EnergyMinimizer

    use ODE_Solvers
    use integrationDataTypes
    implicit none

    !> Interface of the field routine: H is the effective field at the state m, both of size 3n
    abstract interface
        subroutine heff_fct( m, H )
            real,dimension(:),intent(in) :: m
            real,dimension(:),intent(out) :: H
        end subroutine heff_fct
    end interface

    !> Interface of the energy routine: the energy terms [J] of the state the field was last
    !> evaluated at. Only the sum is used here; four terms are carried so the caller can reuse it.
    abstract interface
        subroutine energy_fct( E )
            real,dimension(4),intent(out) :: E
        end subroutine energy_fct
    end interface

    !> Settings of the minimizer
    type MinimizerSettings
        real :: tol = 1.0e-5          !> Convergence criterion on max_i |m_i x H_i| / H_scale
        integer :: maxiter = 10000    !> Iteration cap per call (per round when the fallback is used)
        real :: maxrot = 0.3          !> Largest rotation of any vector in one iteration [rad]
        integer :: fallback = 1       !> 1: fall back to the ODE integration when stalled, 0: give up
        integer :: saddle_check = 1   !> 0: accept a converged state as it is, 1: nudge it and relax again, 2: lowest Hessian eigenvalue by Lanczos, descending along the eigenvector when it is negative
        real :: saddle_pert = 1.0e-2  !> Size of that nudge, the rotation of each vector [rad]
        real :: H_scale = 1.0         !> Field scale that makes the torque criterion dimensionless, e.g. max(Ms)
        real :: E_scale = 1.0         !> Energy scale for the watchdog tolerance, e.g. 1/2 mu0 Ms^2 V
        integer :: display_every = 100 !> Progress message every this many iterations
        real :: tau_init = 0.0        !> Step length of the first iteration [rad per field unit]; 0 rotates the most-torqued vector by rot_init instead
        real :: eig_tol = 1.0e-2      !> saddle_check 2: relative accuracy of the lowest Hessian eigenvalue (absolute floor eig_tol * H_scale / 100)
        integer :: eig_maxiter = 60   !> saddle_check 2: cap on the Lanczos steps per check (one field evaluation each)
        real,dimension(:),allocatable :: weights !> saddle_check 2: weight of every vector in the energy, Ms_i V_i; uniform when not allocated
        real,dimension(:),allocatable :: eig_start !> saddle_check 2: Lanczos start vector (3n), e.g. the eigenvector of a similar state; random when not allocated
    end type MinimizerSettings

    !> What the minimizer reports back
    type MinimizerResult
        integer :: n_iter = 0         !> Iterations over all rounds
        real :: torque = 0.0          !> Final max_i |m_i x H_i| / H_scale
        integer :: status = 2         !> 0 converged, 1 converged after an ODE fallback, 2 not converged
        real :: tau_last = 0.0        !> Step length of the last iteration, a curvature estimate the next call can start from (tau_init)
        real :: eig_min = 0.0         !> saddle_check 2: lowest Hessian eigenvalue of the returned state, relative to H_scale
        integer :: n_eig = 0          !> saddle_check 2: Lanczos steps spent over all checks
        real,dimension(:),allocatable :: eig_vec !> saddle_check 2: the eigenvector of eig_min (3n), tangent to the returned state
    end type MinimizerResult

    integer,parameter,private :: n_hist = 10         !> Window of the non-monotone energy watchdog
    integer,parameter,private :: n_shrink_max = 30   !> Consecutive step halvings before the iteration is declared stalled
    integer,parameter,private :: n_round_max = 3     !> Minimizer rounds: the first plus up to two restarts after a fallback
    integer,parameter,private :: n_saddle_max = 3    !> Saddle checks per round before the state is accepted as it is
    integer,parameter,private :: n_saddle_iter = 100 !> Iterations a perturbed state gets to fall below the unperturbed energy
    real,parameter,private :: rot_init = 1.0e-3      !> Rotation of the most-torqued vector in the first step [rad]
    real,parameter,private :: E_rise_rel = 1.0e-6    !> Allowed energy rise per step, relative to E_scale
    real,parameter,private :: eig_step = 1.0e-4      !> Rotation of the most displaced vector in the finite-difference Hessian product [rad]
    logical,parameter,private :: eig_debug = .false. !> Report every Lanczos step (development aid)
    integer,parameter,private :: eig_minit = 3       !> Fewest Lanczos steps before the convergence test is trusted, warm start
    integer,parameter,private :: eig_minit_moved = 10 !> The same after a push off a saddle, when the start vector is mostly random

    contains

    !>-----------------------------------------
    !> Minimizes the energy starting from m0 and returns the final state in m.
    !> @param[in] heff routine returning the effective field at a state
    !> @param[in] energy routine returning the energy terms of the last evaluated state
    !> @param[in] m0 initial state (3n)
    !> @param[inout] m final state (3n)
    !> @param[in] n number of vectors
    !> @param[in] settings tolerances and scales, see MinimizerSettings
    !> @param[inout] result iterations, final torque and status
    !> @param[in] callback progress messages
    !> The remaining arguments are only used by the fallback and are passed straight to MagTense_ODE:
    !> the derivative routine fct, the thermal routine fct_thermal (never called, the fallback runs
    !> without a thermal field), the output times t, the ODE tolerances and the convergence settings.
    !>-----------------------------------------
    subroutine MagTense_Minimize( heff, energy, m0, m, n, settings, result, callback, &
                                  fct, fct_thermal, t, tol, thres_value, useCVODE, t_conv, conv_tol, callback_display )
    procedure(heff_fct) :: heff
    procedure(energy_fct) :: energy
    real,dimension(:),intent(in) :: m0
    real,dimension(:),intent(inout) :: m
    integer,intent(in) :: n
    type(MinimizerSettings),intent(in) :: settings
    type(MinimizerResult),intent(inout) :: result
    procedure(callback_fct), pointer :: callback
    procedure(dydt_fct), pointer :: fct
    procedure(no_argument_fct), pointer :: fct_thermal
    real,dimension(:),intent(in) :: t, t_conv
    real,intent(in) :: tol, thres_value, conv_tol
    integer,intent(in) :: useCVODE, callback_display

    real,dimension(:),allocatable :: H, m_prev, g, g_prev, tq, tq_prev, m_best, u, w
    real,dimension(:),allocatable :: t_out_fb
    real,dimension(:,:),allocatable :: y_fb
    real,dimension(n_hist) :: E_hist
    real,dimension(4) :: E4
    real :: E, E_best, tq_max, tq_rel, tq_best, tau, E_stop, lambda, umax
    integer :: n_eig
    logical :: eig_ok
    integer :: k, i_round, i_sc, n_iter_total, status, nt, n_display, k_cap
    logical :: converged, stalled
    character*(256) :: prog_str

    allocate( H(3*n), m_prev(3*n), g(3*n), g_prev(3*n), tq(3*n), tq_prev(3*n), m_best(3*n) )
    n_display = max( 1, settings%display_every )

    m = m0
    status = 2
    n_iter_total = 0
    tq_rel = 0.0

    do i_round = 1, n_round_max

        k_cap = settings%maxiter
        E_stop = -huge(1.0)
        call startState()
        call descend()
        n_iter_total = n_iter_total + k

        if ( converged .and. settings%saddle_check .eq. 1 ) then
            !A vanishing torque also marks a saddle point, and a symmetric starting state - the
            !canonical vortex of standard problem 3, or a magnetization exactly antiparallel to the
            !field in a hysteresis loop - sits on one. Steepest descent then converges onto it
            !instead of falling off, which the time integration escapes only through rounding
            !noise. So the converged state is nudged by a small random rotation and relaxed again:
            !a minimum takes the nudge back and is kept as it was, a saddle lets the state fall to
            !a lower energy, which is then checked in the same way.
            !The perturbed descent is not run to convergence: a minimum is recognised by the energy
            !staying above the unperturbed value for n_saddle_iter iterations, a saddle by the energy
            !falling below it, which ends the trial descent at once.
            do i_sc = 1, n_saddle_max
                m_best = m
                E_best = E
                tq_best = tq_rel
                call perturb( m, n, settings%saddle_pert )
                k_cap = n_saddle_iter
                E_stop = E_best - E_rise_rel * settings%E_scale
                call startState()
                call descend()
                n_iter_total = n_iter_total + k
                if ( E .lt. E_stop ) then
                    write(prog_str,'(A,ES10.3,A,I4,A)') 'Minimizer: saddle check lowered the energy by ', E_best - E, ' J after ', k, ' iterations, continuing'
                    call callback( trim(prog_str), -1 )
                    !Relax the escaped state to its minimum, then check that one as well
                    k_cap = settings%maxiter
                    E_stop = -huge(1.0)
                    call startState()
                    call descend()
                    n_iter_total = n_iter_total + k
                    if ( .not. converged ) exit
                else
                    !Back to the unperturbed minimum, whichever way the trial descent ended
                    m = m_best
                    E = E_best
                    tq_rel = tq_best
                    converged = .true.
                    exit
                endif
            end do
        endif

        if ( converged .and. settings%saddle_check .eq. 2 ) then
            !The rigorous version of the check above: the lowest eigenvalue of the energy Hessian in
            !the tangent space, from Lanczos on finite-difference Hessian products (one field
            !evaluation each). Positive means minimum. Negative means saddle, and the eigenvector is
            !the direction of steepest descent out of it, so the state is pushed along it and
            !relaxed again, then checked again.
            allocate( u(3*n), w(n) )
            if ( allocated(settings%weights) ) then
                w = settings%weights * ( real(n) / sum(settings%weights) )
            else
                w = 1.0
            endif
            do i_sc = 1, n_saddle_max
                if ( i_sc .eq. 1 .and. allocated(settings%eig_start) ) then
                    u = settings%eig_start
                else if ( i_sc .eq. 1 ) then
                    u = 0.0
                endif
                !The finite differences need the gradient at exactly this state, so it is evaluated
                !afresh rather than taken from the descent (one field evaluation).
                call heff( m, H )
                call torqueAndGradient( m, H, n, tq, g, tq_max )
                !u enters as the start vector (the eigenvector of the previous check or field, or zero
                !for a random one) and leaves as the eigenvector
                !Between neighbouring fields the handed-over eigenvector is close and a small random
                !admixture with a few steps suffices; after a push off a saddle the state has moved and
                !the search starts from an equal mix with more steps.
                if ( i_sc .eq. 1 ) then
                    call lowestEigenpair( heff, m, g, n, w, settings%H_scale, settings%eig_tol, settings%eig_maxiter, &
                                          eig_minit, 0.1, lambda, u, n_eig, eig_ok, callback )
                else
                    call lowestEigenpair( heff, m, g, n, w, settings%H_scale, settings%eig_tol, settings%eig_maxiter, &
                                          eig_minit_moved, 1.0, lambda, u, n_eig, eig_ok, callback )
                endif
                result%eig_vec = u
                result%n_eig = result%n_eig + n_eig
                result%eig_min = lambda / settings%H_scale
                if ( .not. eig_ok ) then
                    write(prog_str,'(A,I4,A,ES10.3)') 'Minimizer: Hessian eigenvalue not converged in ', n_eig, ' steps, lambda/H ', lambda / settings%H_scale
                    call callback( trim(prog_str), -1 )
                endif
                if ( lambda .ge. -1.0e-2 * settings%eig_tol * settings%H_scale ) exit
                write(prog_str,'(A,ES10.3,A)') 'Minimizer: saddle, lowest Hessian eigenvalue ', lambda / settings%H_scale, ' H_scale, descending'
                call callback( trim(prog_str), -1 )
                !Push along the eigenvector so that the most displaced vector rotates by saddle_pert
                umax = sqrt( maxval( u(1:n)**2 + u(n+1:2*n)**2 + u(2*n+1:3*n)**2 ) )
                m = m + ( settings%saddle_pert / max( umax, tiny(1.0) ) ) * u
                call normalizeVectors( m, n )
                k_cap = settings%maxiter
                E_stop = -huge(1.0)
                call startState()
                call descend()
                n_iter_total = n_iter_total + k
                if ( .not. converged ) exit
            end do
            deallocate( u, w )
        endif

        if ( converged ) then
            if ( i_round .eq. 1 ) then
                status = 0
            else
                status = 1
            endif
            exit
        endif

        if ( stalled ) then
            write(prog_str,'(A,I7,A,ES9.2)') 'Minimizer stalled after iter ', k, ', torque ', tq_rel
        else
            write(prog_str,'(A,I7,A,ES9.2)') 'Minimizer hit the iteration cap ', k, ', torque ', tq_rel
        endif
        call callback( trim(prog_str), -1 )

        if ( settings%fallback .eq. 0 .or. i_round .eq. n_round_max ) exit

        !Fall back to the ODE integration over the requested time window from the current state,
        !then try the minimizer again from where it ends up
        call callback( 'Falling back to LL time integration, then restarting the minimizer', -1 )
        nt = size( t )
        allocate( t_out_fb(nt), y_fb(3*n, nt) )
        call MagTense_ODE( fct, t, m, t_out_fb, y_fb, .false., fct_thermal, callback, &
                callback_display, tol, thres_value, useCVODE, t_conv, conv_tol )
        m = y_fb(:,nt)
        deallocate( t_out_fb, y_fb )
    end do

    result%n_iter = result%n_iter + n_iter_total
    result%torque = tq_rel
    result%status = status
    result%tau_last = tau

    deallocate( H, m_prev, g, g_prev, tq, tq_prev, m_best )

    contains

    !> Field, torque, energy and first step length at the current state m
    subroutine startState()
        call heff( m, H )
        call torqueAndGradient( m, H, n, tq, g, tq_max )
        call energy( E4 )
        E = sum( E4 )
        E_hist = E
        tq_rel = tq_max / settings%H_scale
        converged = tq_rel .lt. settings%tol
        stalled = .false.
        !The first step: rot_init on the most-torqued vector, unless the caller hands over the step
        !length of a previous, similar minimization - a start close to the minimum (the secant
        !predictor) has a torque far smaller than the curvature, and a fixed rotation overshoots it
        if ( settings%tau_init .gt. 0.0 ) then
            tau = min( settings%tau_init, settings%maxrot / max( tq_max, tiny(1.0) ) )
        else
            tau = rot_init / max( tq_max, tiny(1.0) )
        endif
        k = 0
    end subroutine startState

    !> The Barzilai-Borwein descent from the current state until converged, stalled or capped
    subroutine descend()
        real :: E_ref, tau_cap, ss, sy, yy
        integer :: n_shrink

        do while ( .not. converged .and. k .lt. k_cap )
            k = k + 1
            m_prev = m
            g_prev = g
            tq_prev = tq

            !Trial step, halved until the energy watchdog is satisfied
            n_shrink = 0
            do
                call cayleyUpdate( m_prev, tq_prev, g_prev, tau, n, m )
                call heff( m, H )
                call energy( E4 )
                E = sum( E4 )
                E_ref = maxval( E_hist )
                if ( E .le. E_ref + E_rise_rel * settings%E_scale ) exit
                if ( n_shrink .ge. n_shrink_max ) exit
                tau = 0.5 * tau
                n_shrink = n_shrink + 1
            end do
            if ( n_shrink .ge. n_shrink_max ) then
                !No step length lowers the energy: the gradient is not a descent direction to within
                !the noise of the field evaluation, so the state is at a saddle or the minimum is
                !resolved as well as the field precision allows
                stalled = .true.
                exit
            endif
            E_hist(1:n_hist-1) = E_hist(2:n_hist)
            E_hist(n_hist) = E
            if ( E .lt. E_stop ) exit

            call torqueAndGradient( m, H, n, tq, g, tq_max )
            tq_rel = tq_max / settings%H_scale
            if ( tq_rel .lt. settings%tol ) then
                converged = .true.
                exit
            endif

            !Barzilai-Borwein step for the next iteration, alternating between the two rules
            ss = sum( (m - m_prev)**2 )
            sy = sum( (m - m_prev) * (g - g_prev) )
            yy = sum( (g - g_prev)**2 )
            if ( sy .gt. 0.0 .and. yy .gt. 0.0 ) then
                if ( mod(k,2) .eq. 0 ) then
                    tau = ss / sy
                else
                    tau = sy / yy
                endif
            else
                !Negative curvature along the step: grow the step and let the cap and the watchdog bound it
                tau = 2.0 * tau
            endif
            tau_cap = settings%maxrot / max( tq_max, tiny(1.0) )
            tau = min( tau, tau_cap )

            !if ( mod(k, n_display) .eq. 0 ) then
            !    write(prog_str,'(A,I7,A,ES9.2,A,ES13.6)') 'Minimizer iter ', k, ' torque ', tq_rel, ' E [J] ', E
            !    call callback( trim(prog_str), -1 )
            !endif
        end do
    end subroutine descend

    end subroutine MagTense_Minimize

    !>-----------------------------------------
    !> Renormalizes every vector of m to unit length
    !>-----------------------------------------
    subroutine normalizeVectors( m, n )
    real,dimension(:),intent(inout) :: m
    integer,intent(in) :: n
    real,dimension(:),allocatable :: nrm
    allocate( nrm(n) )
    nrm = sqrt( m(1:n)**2 + m(n+1:2*n)**2 + m(2*n+1:3*n)**2 )
    m(1:n)       = m(1:n) / nrm
    m(n+1:2*n)   = m(n+1:2*n) / nrm
    m(2*n+1:3*n) = m(2*n+1:3*n) / nrm
    deallocate( nrm )
    end subroutine normalizeVectors

    !>-----------------------------------------
    !> Removes from every v_i its component along the unit vector m_i
    !>-----------------------------------------
    subroutine projectTangent( m, v, n )
    real,dimension(:),intent(in) :: m
    real,dimension(:),intent(inout) :: v
    integer,intent(in) :: n
    real,dimension(:),allocatable :: c
    allocate( c(n) )
    c = m(1:n) * v(1:n) + m(n+1:2*n) * v(n+1:2*n) + m(2*n+1:3*n) * v(2*n+1:3*n)
    v(1:n)       = v(1:n)       - c * m(1:n)
    v(n+1:2*n)   = v(n+1:2*n)   - c * m(n+1:2*n)
    v(2*n+1:3*n) = v(2*n+1:3*n) - c * m(2*n+1:3*n)
    deallocate( c )
    end subroutine projectTangent

    !>-----------------------------------------
    !> Hessian-vector product in field units at the state m, for a tangent displacement v,
    !>     (K v)_i = -P_i (A v)_i + (m_i . H_i) v_i,
    !> i.e. the linearised negative field plus the curvature of the sphere, obtained from one field
    !> evaluation as the one-sided difference of the gradient g = m x (m x H),
    !>     K v = P [ g(m + eps v) - g(m) ] / eps,
    !> with eps chosen so that the most displaced vector rotates by eig_step. The difference also
    !> covers a nonlinear (cubic) anisotropy, which a product with the field operator would not.
    !> g0 is the gradient at m, which the caller already has.
    !>-----------------------------------------
    subroutine hessianProduct( heff, m, g0, v, n, Kv )
    procedure(heff_fct) :: heff
    real,dimension(:),intent(in) :: m, g0, v
    integer,intent(in) :: n
    real,dimension(:),intent(out) :: Kv
    real,dimension(:),allocatable :: m1, H1, tq1, g1
    real :: eps, vmax, tqmax

    vmax = sqrt( maxval( v(1:n)**2 + v(n+1:2*n)**2 + v(2*n+1:3*n)**2 ) )
    if ( vmax .le. tiny(1.0) ) then
        Kv = 0.0
        return
    endif
    allocate( m1(3*n), H1(3*n), tq1(3*n), g1(3*n) )
    eps = eig_step / vmax
    m1 = m + eps * v
    call normalizeVectors( m1, n )
    call heff( m1, H1 )
    call torqueAndGradient( m1, H1, n, tq1, g1, tqmax )
    Kv = ( g1 - g0 ) / eps
    call projectTangent( m, Kv, n )
    deallocate( m1, H1, tq1, g1 )
    end subroutine hessianProduct

    !>-----------------------------------------
    !> Lowest eigenpair of the energy Hessian at the state m, restricted to the tangent space, by
    !> the Lanczos iteration with the matrix-free products of hessianProduct. K is self-adjoint in
    !> the inner product weighted with w_i = Ms_i V_i (mean 1), so the Lanczos vectors are
    !> orthonormalised in that product, with full reorthogonalisation (V is 3n x maxit). The
    !> eigenvalues are those of the energy Hessian divided by mu0 Ms_i V_i, in field units, and
    !> have the same signs: lambda > 0 means minimum, lambda < 0 saddle, with u the direction of
    !> steepest descent out of it. Converged when the residual beta_j |y_j| falls below
    !> eig_tol * H_scale or the Krylov space is exhausted; ok is false when maxit steps did not
    !> get there (lambda and u are then the best available Ritz pair).
    !>-----------------------------------------
    subroutine lowestEigenpair( heff, m, g0, n, w, H_scale, eig_tol, maxit, minit, admix, lambda, u, n_iter, ok, callback )
    procedure(heff_fct) :: heff
    real,dimension(:),intent(in) :: m, g0, w
    integer,intent(in) :: n, maxit, minit   !> minit: fewest steps before the convergence test is trusted
    real,intent(in) :: H_scale, eig_tol, admix   !> admix: weight of a random vector added to the (unit) start vector
    real,intent(out) :: lambda
    real,dimension(:),intent(inout) :: u   !> in: start vector (zero for random), out: eigenvector
    integer,intent(out) :: n_iter
    logical,intent(out) :: ok
    procedure(callback_fct), pointer :: callback
    character*(256) :: prog_str

    real,dimension(:,:),allocatable :: V, T, Z
    real,dimension(:),allocatable :: r, Kv, ev, alpha, beta
    real :: nrm, resid, c, lambda_prev
    integer :: j, i, jmin

    allocate( V(3*n, maxit), r(3*n), Kv(3*n), alpha(maxit), beta(maxit) )
    alpha = 0.0
    beta = 0.0
    ok = .false.
    lambda = 0.0
    n_iter = 0

    !Start vector: the one handed over, projected onto the tangent space and normalised, plus a
    !random tangent vector of weight admix so that every mode is represented - Lanczos only finds
    !the lowest eigenvalue if the start vector has a component along its eigenvector, and the
    !eigenvector of a neighbouring state need not have one when the state has changed much.
    r = u
    call projectTangent( m, r, n )
    nrm = sqrt( wdot( r, r ) )
    if ( nrm .gt. 1.0e-3 ) then
        r = r / nrm
    else
        r = 0.0
    endif
    call random_number( Kv )
    Kv = 2.0 * Kv - 1.0
    call projectTangent( m, Kv, n )
    nrm = sqrt( wdot( Kv, Kv ) )
    if ( nrm .gt. tiny(1.0) ) r = r + max( admix, merge( 1.0, 0.0, sum(r**2) .eq. 0.0 ) ) * Kv / nrm
    nrm = sqrt( wdot( r, r ) )
    u = 0.0
    lambda_prev = huge(1.0)
    if ( nrm .le. tiny(1.0) ) then
        deallocate( V, r, Kv, alpha, beta )
        return
    endif
    V(:,1) = r / nrm

    do j = 1, maxit
        call hessianProduct( heff, m, g0, V(:,j), n, Kv )
        alpha(j) = wdot( V(:,j), Kv )
        r = Kv - alpha(j) * V(:,j)
        if ( j .gt. 1 ) r = r - beta(j-1) * V(:,j-1)
        do i = 1, j
            c = wdot( V(:,i), r )
            r = r - c * V(:,i)
        end do
        call projectTangent( m, r, n )
        beta(j) = sqrt( wdot( r, r ) )

        !Lowest eigenpair of the tridiagonal T_j
        allocate( T(j,j), Z(j,j), ev(j) )
        T = 0.0
        do i = 1, j
            T(i,i) = alpha(i)
            if ( i .lt. j ) then
                T(i,i+1) = beta(i)
                T(i+1,i) = beta(i)
            endif
        end do
        call jacobiEigen( T, j, ev, Z )
        jmin = minloc( ev, 1 )
        lambda = ev(jmin)
        resid = beta(j) * abs( Z(j,jmin) )
        n_iter = j
        if ( eig_debug ) then
            write(prog_str,'(A,I4,A,ES12.4,A,ES12.4,A,ES12.4)') 'Lanczos ', j, ' lambda/H ', lambda / H_scale, ' resid/H ', resid / H_scale, ' beta/H ', beta(j) / H_scale
            call callback( trim(prog_str), -1 )
        endif
        !Bauer-Fike: the true eigenvalue lies within resid of the Ritz value, so this bounds the
        !relative error of lambda by eig_tol, with an absolute floor near lambda = 0. The Krylov
        !space is exhausted when beta vanishes (a single vector has only two tangent directions).
        !At least minit steps, and the lowest Ritz value must have stopped falling (it decreases
        !monotonically with j), so that a start vector that happens to be an eigenvector of a
        !higher mode is not accepted on its own residual.
        if ( ( j .ge. minit .and. resid .lt. eig_tol * max( abs(lambda), 1.0e-2 * H_scale ) .and. &
               lambda_prev - lambda .lt. eig_tol * max( abs(lambda), 1.0e-2 * H_scale ) ) .or. &
             beta(j) .le. 1.0e-12 * max( H_scale, abs(lambda) ) ) ok = .true.
        lambda_prev = lambda
        if ( ok .or. j .eq. maxit ) then
            !Ritz vector of the lowest eigenvalue
            u = 0.0
            do i = 1, j
                u = u + Z(i,jmin) * V(:,i)
            end do
            call projectTangent( m, u, n )
            deallocate( T, Z, ev )
            exit
        endif
        deallocate( T, Z, ev )
        V(:,j+1) = r / beta(j)
    end do
    deallocate( V, r, Kv, alpha, beta )

    contains

        !> Inner product weighted per vector with w
        real function wdot( a, b )
        real,dimension(:),intent(in) :: a, b
        wdot = sum( w * ( a(1:n) * b(1:n) + a(n+1:2*n) * b(n+1:2*n) + a(2*n+1:3*n) * b(2*n+1:3*n) ) )
        end function wdot

    end subroutine lowestEigenpair

    !>-----------------------------------------
    !> Eigenvalues ev and eigenvectors (the columns of Z) of the symmetric matrix A by cyclic Jacobi
    !> rotations. A is destroyed. Meant for the small Lanczos tridiagonal, not for large matrices.
    !>-----------------------------------------
    subroutine jacobiEigen( A, n, ev, Z )
    integer,intent(in) :: n
    real,dimension(n,n),intent(inout) :: A
    real,dimension(n),intent(out) :: ev
    real,dimension(n,n),intent(out) :: Z
    integer :: sweep, p, q, k
    real :: off, diag, theta, t, c, s, tau, apq, app, aqq, akp, akq, zkp, zkq

    Z = 0.0
    do k = 1, n
        Z(k,k) = 1.0
    end do
    do sweep = 1, 100
        off = 0.0
        diag = 0.0
        do p = 1, n
            diag = diag + A(p,p)**2
            do q = p+1, n
                off = off + A(p,q)**2
            end do
        end do
        if ( off .le. 1.0e-30 * max( diag, tiny(1.0) ) ) exit
        do p = 1, n-1
            do q = p+1, n
                apq = A(p,q)
                if ( abs(apq) .le. tiny(1.0) ) cycle
                theta = ( A(q,q) - A(p,p) ) / ( 2.0 * apq )
                t = sign( 1.0, theta ) / ( abs(theta) + sqrt( theta**2 + 1.0 ) )
                c = 1.0 / sqrt( t**2 + 1.0 )
                s = t * c
                tau = s / ( 1.0 + c )
                app = A(p,p)
                aqq = A(q,q)
                A(p,p) = app - t * apq
                A(q,q) = aqq + t * apq
                A(p,q) = 0.0
                A(q,p) = 0.0
                do k = 1, n
                    if ( k .ne. p .and. k .ne. q ) then
                        akp = A(k,p)
                        akq = A(k,q)
                        A(k,p) = akp - s * ( akq + tau * akp )
                        A(k,q) = akq + s * ( akp - tau * akq )
                        A(p,k) = A(k,p)
                        A(q,k) = A(k,q)
                    endif
                    zkp = Z(k,p)
                    zkq = Z(k,q)
                    Z(k,p) = zkp - s * ( zkq + tau * zkp )
                    Z(k,q) = zkq + s * ( zkp - tau * zkq )
                end do
            end do
        end do
    end do
    do k = 1, n
        ev(k) = A(k,k)
    end do
    end subroutine jacobiEigen

    !>-----------------------------------------
    !> Rotates every vector by a random angle of the order of pert [rad] and renormalizes. Uses the
    !> intrinsic generator, so the sequence follows whatever seeding the caller has done.
    !>-----------------------------------------
    subroutine perturb( m, n, pert )
    real,dimension(:),intent(inout) :: m
    integer,intent(in) :: n
    real,intent(in) :: pert
    real,dimension(:),allocatable :: r, nrm
    real,dimension(3) :: r0

    !A common tilt of all vectors, which is smooth and cheap to relax back, plus a ten times
    !smaller independent tilt of every vector, which reaches unstable modes that are orthogonal to
    !a uniform one by symmetry
    allocate( r(3*n), nrm(n) )
    call random_number( r0 )
    call random_number( r )
    m(1:n)       = m(1:n)       + pert * ( 2.0 * r0(1) - 1.0 )
    m(n+1:2*n)   = m(n+1:2*n)   + pert * ( 2.0 * r0(2) - 1.0 )
    m(2*n+1:3*n) = m(2*n+1:3*n) + pert * ( 2.0 * r0(3) - 1.0 )
    m = m + 0.1 * pert * ( 2.0 * r - 1.0 )
    nrm = sqrt( m(1:n)**2 + m(n+1:2*n)**2 + m(2*n+1:3*n)**2 )
    m(1:n)       = m(1:n) / nrm
    m(n+1:2*n)   = m(n+1:2*n) / nrm
    m(2*n+1:3*n) = m(2*n+1:3*n) / nrm
    deallocate( r, nrm )
    end subroutine perturb

    !>-----------------------------------------
    !> Torque tq = m x H and descent direction g = m x tq = m x (m x H) of every vector, both in the
    !> (3n) layout of m. tq_max is the largest |tq_i|.
    !>-----------------------------------------
    subroutine torqueAndGradient( m, H, n, tq, g, tq_max )
    real,dimension(:),intent(in) :: m, H
    integer,intent(in) :: n
    real,dimension(:),intent(out) :: tq, g
    real,intent(out) :: tq_max

    tq(1:n)       = m(n+1:2*n) * H(2*n+1:3*n) - m(2*n+1:3*n) * H(n+1:2*n)
    tq(n+1:2*n)   = m(2*n+1:3*n) * H(1:n) - m(1:n) * H(2*n+1:3*n)
    tq(2*n+1:3*n) = m(1:n) * H(n+1:2*n) - m(n+1:2*n) * H(1:n)

    g(1:n)       = m(n+1:2*n) * tq(2*n+1:3*n) - m(2*n+1:3*n) * tq(n+1:2*n)
    g(n+1:2*n)   = m(2*n+1:3*n) * tq(1:n) - m(1:n) * tq(2*n+1:3*n)
    g(2*n+1:3*n) = m(1:n) * tq(n+1:2*n) - m(n+1:2*n) * tq(1:n)

    tq_max = sqrt( maxval( tq(1:n)**2 + tq(n+1:2*n)**2 + tq(2*n+1:3*n)**2 ) )
    end subroutine torqueAndGradient

    !>-----------------------------------------
    !> The sphere-preserving descent step
    !>     m1 = ( (1 - tau^2 |tq|^2 / 4) m0 - tau g ) / (1 + tau^2 |tq|^2 / 4),   tq = m0 x H,  g = m0 x tq,
    !> which is the exact solution of the midpoint rule m1 = m0 - (tau/2) (m0 + m1) x tq. Every vector is
    !> rotated about its own torque axis by the angle 2 atan(tau |tq|/2), so |m1| = |m0| to rounding.
    !>-----------------------------------------
    subroutine cayleyUpdate( m0, tq, g, tau, n, m1 )
    real,dimension(:),intent(in) :: m0, tq, g
    real,intent(in) :: tau
    integer,intent(in) :: n
    real,dimension(:),intent(out) :: m1
    real,dimension(:),allocatable :: q, denom

    allocate( q(n), denom(n) )
    q = 0.25 * tau**2 * ( tq(1:n)**2 + tq(n+1:2*n)**2 + tq(2*n+1:3*n)**2 )
    denom = 1.0 / ( 1.0 + q )
    m1(1:n)       = ( (1.0 - q) * m0(1:n)       - tau * g(1:n) )       * denom
    m1(n+1:2*n)   = ( (1.0 - q) * m0(n+1:2*n)   - tau * g(n+1:2*n) )   * denom
    m1(2*n+1:3*n) = ( (1.0 - q) * m0(2*n+1:3*n) - tau * g(2*n+1:3*n) ) * denom
    deallocate( q, denom )
    end subroutine cayleyUpdate

end module EnergyMinimizer
