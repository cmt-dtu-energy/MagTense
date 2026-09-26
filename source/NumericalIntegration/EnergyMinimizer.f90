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
        integer :: saddle_check = 1   !> 1: nudge a converged state and relax again to make sure it is a minimum, 0: skip
        real :: saddle_pert = 1.0e-2  !> Size of that nudge, the rotation of each vector [rad]
        real :: H_scale = 1.0         !> Field scale that makes the torque criterion dimensionless, e.g. max(Ms)
        real :: E_scale = 1.0         !> Energy scale for the watchdog tolerance, e.g. 1/2 mu0 Ms^2 V
        integer :: display_every = 100 !> Progress message every this many iterations
        real :: tau_init = 0.0        !> Step length of the first iteration [rad per field unit]; 0 rotates the most-torqued vector by rot_init instead
    end type MinimizerSettings

    !> What the minimizer reports back
    type MinimizerResult
        integer :: n_iter = 0         !> Iterations over all rounds
        real :: torque = 0.0          !> Final max_i |m_i x H_i| / H_scale
        integer :: status = 2         !> 0 converged, 1 converged after an ODE fallback, 2 not converged
        real :: tau_last = 0.0        !> Step length of the last iteration, a curvature estimate the next call can start from (tau_init)
    end type MinimizerResult

    integer,parameter,private :: n_hist = 10         !> Window of the non-monotone energy watchdog
    integer,parameter,private :: n_shrink_max = 30   !> Consecutive step halvings before the iteration is declared stalled
    integer,parameter,private :: n_round_max = 3     !> Minimizer rounds: the first plus up to two restarts after a fallback
    integer,parameter,private :: n_saddle_max = 3    !> Saddle checks per round before the state is accepted as it is
    integer,parameter,private :: n_saddle_iter = 100 !> Iterations a perturbed state gets to fall below the unperturbed energy
    real,parameter,private :: rot_init = 1.0e-3      !> Rotation of the most-torqued vector in the first step [rad]
    real,parameter,private :: E_rise_rel = 1.0e-6    !> Allowed energy rise per step, relative to E_scale

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

    real,dimension(:),allocatable :: H, m_prev, g, g_prev, tq, tq_prev, m_best
    real,dimension(:),allocatable :: t_out_fb
    real,dimension(:,:),allocatable :: y_fb
    real,dimension(n_hist) :: E_hist
    real,dimension(4) :: E4
    real :: E, E_best, tq_max, tq_rel, tq_best, tau, E_stop
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

        if ( converged .and. settings%saddle_check .ne. 0 ) then
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
