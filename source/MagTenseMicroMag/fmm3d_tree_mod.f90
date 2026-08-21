module fmm3d_tree_mod
    use omp_lib
    use omp_mod
    use trace_mod
        implicit none

    ! On Windows the python extension is compiled with
    ! /assume:underscore /names:lowercase (EXTRA_FFLAGS in the root Makefile)
    ! while this library is not, so the extension looks for these procedures
    ! under Linux-style names while ifx would emit FMM3D_TREE_MOD_mp_BUILD1
    ! and friends. Each procedure below therefore carries an explicit ALIAS,
    ! the same convention already used for the other f2py entry points (see
    ! LandauLifshitzEquationSolver.f90). The alias strings are identical to
    ! ifx's default mangling on Linux, so they change nothing there.

      type :: FMM3DTree
            logical :: is_built = .false.
            logical :: keep_tree = .true.

            integer :: nterms_in = -1

            integer :: nlmax = 51
            integer :: nlevels = 0
            integer :: nboxes = 0
            integer(8) :: ltree = 0
            integer :: nlmin = 0
            integer :: iper = 0
            integer :: ifunif = 1

            integer :: ifnear = 1
            integer :: nexpc = 0
            integer :: nadd = 0
            integer :: ntj = 0

            integer :: mnlist1 = 0
            integer :: mnlist2 = 0
            integer :: mnlist3 = 0
            integer :: mnlist4 = 0
            integer :: mnbors = 27
            integer :: isep = 1

            integer :: ifprint = 0

            integer :: nthd = 1

            !------------------- pointer to input and output 
            double precision :: targ(3,1)
            double precision, contiguous, pointer ::  source(:,:)
            double precision, contiguous, pointer ::  dipvec(:,:,:)
            double precision, contiguous, pointer ::  grad(:,:,:)
            double precision, contiguous, pointer :: sourcesort(:,:)
            double precision, contiguous, pointer :: dipvecsort(:,:,:)
            double precision, contiguous, pointer :: gradsort(:,:,:)
            !-----------------------------------------
            integer nsource,ntarg 
            integer nd,ier
            double precision eps
            !-----------------------------------------
            integer :: ndiv = 0 
            integer :: idivflag
            integer(8), dimension(8) :: ipointer
            double precision expc(3)
            integer :: iexpc
            integer, contiguous, pointer :: itree(:)
            double precision, contiguous, pointer :: boxsize(:), treecenters(:,:), centers(:,:)
            integer, contiguous, pointer :: isrcse(:,:),iexpcse(:,:)
            integer, contiguous, pointer :: isrc(:),itarg(:)
            !-----------------------------------------
            double precision b0,b0inv,b0inv2,b0inv3

            integer :: nmax
            integer, contiguous, pointer :: nterms(:)
            integer :: lmptemp
            integer(8) :: lmptot
            integer(8), contiguous, pointer :: iaddr(:,:)
            double precision, contiguous, pointer :: mptemp(:),mptemp2(:)
            double precision, contiguous, pointer :: rmlexp(:)
            double precision, contiguous, pointer :: scales(:)



            double precision expcsort(3)
            double precision scjsort(1),radexp
            double precision timeinfo(6)


            !-----
            double precision thresh
            integer, contiguous, pointer :: list1(:,:),nlist1(:)
            integer, contiguous, pointer :: list2(:,:),nlist2(:)
            integer, contiguous, pointer :: list3(:,:),nlist3(:)
            integer, contiguous, pointer :: list4(:,:),nlist4(:)

            integer :: nlams

            double precision, contiguous, pointer :: rlams(:), whts(:)
            integer, contiguous, pointer :: nphysical(:), nfourier(:)



        double precision, contiguous, pointer  :: carray(:,:), dc(:,:)
        double precision, contiguous, pointer  :: cs(:,:),fact(:),rdplus(:,:,:)
        double precision, contiguous, pointer  :: rdminus(:,:,:), rdsq3(:,:,:)
        double precision, contiguous, pointer  :: rdmsq3(:,:,:)
        double precision, contiguous, pointer  :: rscpow(:)
        double precision, contiguous, pointer  :: rlsc(:,:,:)


        double precision, contiguous, pointer :: scarray(:,:)



        integer nexptot, nexptotp, nthmax, nphmax
        double complex, contiguous, pointer  :: xshift(:,:)
        double complex, contiguous, pointer  :: yshift(:,:)
        double precision, contiguous, pointer :: zshift(:,:)

        double complex, contiguous, pointer  :: fexpe(:),fexpo(:),fexpback(:)
        double complex, contiguous, pointer  :: mexp(:,:,:,:)
        double complex, contiguous, pointer  :: mexpf1(:,:,:),mexpf2(:,:,:)
        double complex, contiguous, pointer  :: &
        &    mexpp1(:,:,:),mexpp2(:,:,:),mexppall(:,:,:,:)

        double complex, contiguous, pointer  :: tmp(:,:,:,:)
        double precision, contiguous, pointer :: mptmp(:,:)


        integer, contiguous, pointer :: list4ct(:)
        integer, contiguous, pointer :: ilist4(:)
        integer :: cntlist4



        integer, contiguous, pointer ::  laddr(:,:)


        integer :: nlege = 100
        integer :: lw7 = 40000
        double precision, contiguous, pointer ::  wlege(:)
        integer :: lused7, lca



        double complex, contiguous, pointer :: gboxmexp(:,:,:)
        double complex, contiguous, pointer :: gboxwexp(:,:,:,:,:)
        double complex, contiguous, pointer :: pgboxwexp(:,:,:,:)
        double precision, contiguous, pointer :: gboxsubcenters(:,:,:)
        double precision, contiguous, pointer :: gboxsort(:,:,:)
        integer, contiguous, pointer :: gboxind(:,:)
        integer, contiguous, pointer :: gboxfl(:,:,:)
        double precision, contiguous, pointer :: gboxcgsort(:,:,:)
        double precision, contiguous, pointer :: gboxdpsort(:,:,:,:)

        integer, contiguous, pointer :: uall(:,:),dall(:,:),nall(:,:)
        integer, contiguous, pointer :: sall(:,:),eall(:,:),wall(:,:)
        integer, contiguous, pointer :: u1234(:,:),d5678(:,:)
        integer, contiguous, pointer :: n1256(:,:),s3478(:,:)
        integer, contiguous, pointer :: e1357(:,:),w2468(:,:)
        integer, contiguous, pointer :: n12(:,:),n56(:,:),s34(:,:),s78(:,:)
        integer, contiguous, pointer :: e13(:,:),e57(:,:),w24(:,:),w68(:,:)
        integer, contiguous, pointer :: e1(:,:),e3(:,:),e5(:,:),e7(:,:)
        integer, contiguous, pointer :: w2(:,:),w4(:,:),w6(:,:),w8(:,:)

        !----------- list 3 variables 
        double complex, contiguous, pointer :: iboxlexp(:,:,:)
        double precision, contiguous, pointer :: iboxsubcenters(:,:,:)
        double precision, contiguous, pointer :: iboxgrad(:,:,:,:)
        double precision, contiguous, pointer :: iboxsrc(:,:,:)
        integer, contiguous, pointer :: iboxsrcind(:,:)
        integer, contiguous, pointer :: iboxfl(:,:,:)
        !-------- 

        contains
          procedure :: build1
          procedure :: reset_sort_arg
          procedure :: dealloc
          procedure :: build2
          procedure :: lfmm3dmain_tree
          procedure :: reset_expansion_coeff

          procedure :: eval_local
          procedure :: eval_direct

          procedure :: build_tree
          procedure :: make_and_eval
          procedure :: reorder_dipvec
      end type FMM3DTree
    contains 


      subroutine build_tree(self, source, eps, ndiv, ier, ifunif, nlmin, nlmax)
        !DEC$ ATTRIBUTES ALIAS:"fmm3d_tree_mod_mp_build_tree_" :: build_tree
        class(FMM3DTree), intent(inout) :: self
        double precision, contiguous, pointer :: source(:,:)
        double precision eps
        integer ndiv
        integer ier
        integer :: ifunif, nlmin, nlmax
        !------------------------------------------------
        integer :: i
        integer :: ibox, ilev
        integer, save :: itimer = 0
        !-----------------

        call trace%begin( "FMM3DTree_build_tree", itimer=itimer, verbose=2 )

        if (self%is_built) then
            call trace%end( "FMM3DTree_build_tree", itimer=itimer, verbose=2 )
            return
        end if


        if (ndiv .gt.0) then
          self%ndiv = ndiv 
        end if

        self%source => source
        self%nsource = size(source,2)
        self%eps = eps
        self%nd = 1
        self%ier = ier
        self%ntarg = 0
        self%ifunif = ifunif
        self%nlmin = nlmin
        self%nlmax = nlmax

       
        !----------- use omp_mod to get number of threads -----------------
        self%nthd = omp%numthreads()
        !-----------------------------------------------------------------
        

        call self%build1()
        call self%build2()


        !Tree diagnostics. These are only of interest when debugging the FMM setup, and the tree is
        !rebuilt for every problem, so they are gated on the trace flag rather than printed
        !unconditionally to stdout - which, under Matlab, is not even reliably visible.
        if ( trace%enabled .and. trace%verbose .ge. 2 ) then
            print *, " built tree with ", self%nlevels, " levels and ", self%nboxes, " boxes "
            print *, " Number of multipole expansion terms:", self%nterms(0)
            print *, "Number of boxes per level:"
            do i = 0, self%nlevels
                print *, "Level ", i, ": ", self%laddr(2, i) - self%laddr(1, i) + 1, " boxes"
            enddo
            print *, " scaling factor b0 = ", self%b0, " b0inv = ", self%b0inv
        endif


        self%is_built = .true.
        call trace%end( "FMM3DTree_build_tree", itimer=itimer, verbose=2 )
      end subroutine 

      subroutine make_and_eval(self, dipvec, grad)
        !DEC$ ATTRIBUTES ALIAS:"fmm3d_tree_mod_mp_make_and_eval_" :: make_and_eval
        class(FMM3DTree), intent(inout) :: self
        double precision, contiguous, pointer :: dipvec(:,:,:)
        double precision, contiguous, pointer :: grad(:,:,:)
        integer, save :: itimer = 0
        !------------------------------------------------
        call trace%begin( "FMM3DTree_make_and_eval", itimer=itimer, verbose=2 )


        self%grad => grad
        self%dipvec => dipvec

        call self%lfmm3dmain_tree()
        call self%eval_local()


        !======== commented out - all near field is handled by MagTense analytical tensor ==========
        !call self%eval_direct()
        !===========================================================================================

        call dreorderi(3*self%nd,self%nsource,self%gradsort,self%grad,self%isrc)
        call drescale(self%nd*3*self%nsource,self%grad,self%b0inv)

      call trace%end( "FMM3DTree_make_and_eval", itimer=itimer, verbose=2 )


      end subroutine make_and_eval


    subroutine build1(self)
        !DEC$ ATTRIBUTES ALIAS:"fmm3d_tree_mod_mp_build1_" :: build1
        class(FMM3DTree), intent(inout) :: self
        !------------------------------------------------
        !----------------
        integer, contiguous, pointer :: itree(:)
        double precision, contiguous, pointer :: boxsize(:), treecenters(:,:)
        integer, contiguous, pointer :: isrcse(:,:),iexpcse(:,:)
        integer, contiguous, pointer :: isrc(:),itarg(:)
        !--------------------------------
        double precision, contiguous, pointer :: sourcesort(:,:)
        double precision, contiguous, pointer :: dipvecsort(:,:,:)
        double precision, contiguous, pointer :: gradsort(:,:,:) 
        !----------------
        integer, contiguous, pointer :: nterms(:)
        integer(8), contiguous, pointer :: iaddr(:,:)
        double precision, contiguous, pointer :: mptemp(:),mptemp2(:)

        double precision, contiguous, pointer :: rmlexp(:)
         integer :: iert
        double precision, contiguous, pointer :: scales(:)
 
        !----------------
        integer :: i, ilev
        !------------------------------------------------

    
        if (self%ndiv .eq. 0) then
          call lndiv(self%eps,self%nsource,0,0,1,2, &
                  &    0,self%ndiv,self%idivflag) 
        end if


        call pts_tree_mem(self%source,self%nsource,self%targ,0,self%idivflag,self%ndiv,self%nlmin, & 
            &    self%nlmax,self%iper,self%ifunif,self%nlevels,self%nboxes,self%ltree)



        allocate(itree(self%ltree))
        allocate(boxsize(0:self%nlevels))
        allocate(treecenters(3,self%nboxes))

        call pts_tree_build(self%source,self%nsource,self%targ,self%ntarg,self%idivflag,self%ndiv, &
        &    self%nlmin,self%nlmax,self%iper,self%ifunif,self%nlevels,self%nboxes,self%ltree,itree,self%ipointer, &
        &    treecenters,boxsize)
    

        allocate(isrcse(2,self%nboxes),iexpcse(2,self%nboxes))
        allocate(isrc(self%nsource),itarg(self%ntarg))

        call pts_tree_sort(self%nsource,self%source,itree,self%ltree,self%nboxes,self%nlevels, &
        &    self%ipointer,treecenters,isrc,isrcse)
        
        
        call pts_tree_sort(self%nexpc,self%expc,itree,self%ltree,self%nboxes,self%nlevels, &
        &    self%ipointer,treecenters,self%iexpc,iexpcse)

        self%itree => itree
        self%boxsize => boxsize
        self%treecenters => treecenters
        self%centers => treecenters ! same 
        self%isrcse => isrcse
        self%iexpcse => iexpcse
        self%isrc => isrc
        self%itarg => itarg

        !-------- end of tree build --------------------
        !------ set scaling 

        self%b0 = self%boxsize(0)
        self%b0inv = 1.0d0  /self%b0
        self%b0inv2 = self%b0inv**2
        !self%b0inv3 = self%b0inv2*self%b0inv !not used 





        !-------------- allocate sorted source and targ arrays------
        allocate(sourcesort(3,self%nsource))
        allocate(dipvecsort(self%nd,3,self%nsource))
        allocate(gradsort(self%nd,3,self%nsource))
            !------------------------------------------------------
        self%sourcesort => sourcesort
        self%dipvecsort => dipvecsort
        self%gradsort => gradsort
        !--------------------------------------------------------

        allocate(nterms(0:self%nlevels))
        self%nmax = 0
        do i=0,self%nlevels
            !-----------if nterms_in is explicitly set, use that, otherwise compute nterms based on eps ------------------------------
            if (self%nterms_in .gt. 0) then
                nterms(i) = self%nterms_in
            else
              call l3dterms(self%eps,nterms(i))
            end if 
            !-------------------------------------------------------------------------------------------------------------------------
            if(nterms(i).gt.self%nmax) self%nmax = nterms(i)
        enddo
        self%nterms => nterms
!       
      call dreorderf(3,self%nsource,self%source,self%sourcesort,self%isrc)
      call drescale(3*self%nsource,self%sourcesort,self%b0inv)
      call drescale(3*self%nboxes,self%treecenters,self%b0inv)
      call drescale(self%nlevels+1,self%boxsize,self%b0inv)

      allocate(iaddr(2,self%nboxes))
      self%lmptemp = (self%nmax+1)*(2*self%nmax+1)*2*self%nd
      allocate(mptemp(self%lmptemp),mptemp2(self%lmptemp))
      self%mptemp => mptemp
      self%mptemp2 => mptemp2 

      call mpalloc(self%nd,self%itree(self%ipointer(1)),iaddr,self%nlevels,self%lmptot,self%nterms)

        self%iaddr => iaddr

        allocate(rmlexp(self%lmptot),stat=iert)
        self%rmlexp => rmlexp

        allocate(scales(0:self%nlevels))
        do ilev = 0,self%nlevels
            scales(ilev) = boxsize(ilev)
        enddo
        self%scales => scales


        self%laddr(1:2,0:self%nlevels) => self%itree(self%ipointer(1) : self%ipointer(1)+(self%nlevels + 1)*2-1)
    end subroutine build1

    subroutine reorder_dipvec(self)
        !DEC$ ATTRIBUTES ALIAS:"fmm3d_tree_mod_mp_reorder_dipvec_" :: reorder_dipvec
        class(FMM3DTree), intent(inout) :: self
        !------------------------------------------------   
        call dreorderf(3*self%nd,self%nsource,self%dipvec,self%dipvecsort, self%isrc)
        call drescale(3*self%nd*self%nsource,self%dipvecsort,self%b0inv2)
    end subroutine reorder_dipvec


    subroutine build2(self)
        !DEC$ ATTRIBUTES ALIAS:"fmm3d_tree_mod_mp_build2_" :: build2
        class(FMM3DTree), intent(inout) :: self
        !------------------------------------------------
        integer :: i, nn, ilev
        integer :: iert 
        integer :: ibox, nchild, istart, iend 
        integer(kind=8) :: bigint
        integer :: nmaxt, npts
        !----------------

        
        self%thresh = 2.0d0**(-51)*self%boxsize(0)


        call computemnlists(self%nlevels,self%nboxes,self%itree(self%ipointer(1)),self%boxsize, &
            &    self%centers,self%itree(self%ipointer(3)),self%itree(self%ipointer(4)), &
            &    self%itree(self%ipointer(5)),self%isep,self%itree(self%ipointer(6)),self%mnbors, &
            &    self%itree(self%ipointer(7)),self%iper,self%mnlist1,self%mnlist2,self%mnlist3,self%mnlist4)



        allocate(self%list1(self%mnlist1,self%nboxes),self%nlist1(self%nboxes))
        allocate(self%list2(self%mnlist2,self%nboxes),self%nlist2(self%nboxes))
        allocate(self%list3(self%mnlist3,self%nboxes),self%nlist3(self%nboxes))
        allocate(self%list4(self%mnlist4,self%nboxes),self%nlist4(self%nboxes))


        call computelists(self%nlevels,self%nboxes,self%itree(self%ipointer(1)),self%boxsize, &
            &    self%centers,self%itree(self%ipointer(3)),self%itree(self%ipointer(4)), &
            &    self%itree(self%ipointer(5)),self%isep,self%itree(self%ipointer(6)),self%mnbors, &
            &    self%itree(self%ipointer(7)),self%iper,self%nlist1,self%mnlist1,self%list1,self%nlist2, &
            &    self%mnlist2,self%list2,self%nlist3,self%mnlist3,self%list3, &
            &    self%nlist4,self%mnlist4,self%list4)
            

            
        if(self%isep.eq.1) then
            if(self%eps.ge.0.5d-2) self%nlams = 12
            if(self%eps.lt.0.5d-2.and.self%eps.ge.0.5d-3) self%nlams = 12
            if(self%eps.lt.0.5d-3.and.self%eps.ge.0.5d-6) self%nlams = 20
            if(self%eps.lt.0.5d-6.and.self%eps.ge.0.5d-9) self%nlams = 29
            if(self%eps.lt.0.5d-9) self%nlams = 37
        endif
        if(self%isep.eq.2) then
            if(self%eps.ge.0.5d-3) self%nlams = 9
            if(self%eps.lt.0.5d-3.and.self%eps.ge.0.5d-6) self%nlams = 15
            if(self%eps.lt.0.5d-6.and.self%eps.ge.0.5d-9) self%nlams = 22
            if(self%eps.lt.0.5d-9) self%nlams = 29
        endif

        allocate(self%rlams(self%nlams),self%whts(self%nlams))
        allocate(self%nphysical(self%nlams),self%nfourier(self%nlams))




     self%nmax = 0
      do i=0,self%nlevels
         if(self%nmax.lt.self%nterms(i)) self%nmax = self%nterms(i)
      enddo

       allocate(self%rscpow(0:self%nmax))
       allocate(self%carray(4*self%nmax+1,4*self%nmax+1))
       allocate(self%dc(0:4*self%nmax,0:4*self%nmax))
       allocate(self%rdplus(0:self%nmax,0:self%nmax,-self%nmax:self%nmax))
       allocate(self%rdminus(0:self%nmax,0:self%nmax,-self%nmax:self%nmax))
       allocate(self%rdsq3(0:self%nmax,0:self%nmax,-self%nmax:self%nmax))
       allocate(self%rdmsq3(0:self%nmax,0:self%nmax,-self%nmax:self%nmax))
       allocate(self%rlsc(0:self%nmax,0:self%nmax,self%nlams))


      call getpwrotmat(self%nmax,self%carray,self%rdplus,self%rdminus,self%rdsq3,self%rdmsq3,self%dc)
      call vwts(self%rlams,self%whts,self%nlams)
      call numthetahalf(self%nfourier,self%nlams)
      call numthetafour(self%nphysical,self%nlams)
      call rlscini(self%rlsc,self%nlams,self%rlams,self%nmax)


      nn = 10*(self%nmax+2)**2
      allocate(self%scarray(nn,0:self%nlevels))
      do ilev=0,self%nlevels
        call l3dtaevalhessdini(self%nterms(ilev),self%scarray(1,ilev))
      enddo

      self%nexptotp = 0
      self%nexptot = 0
      self%nthmax = 0
      self%nphmax = 0
      nn = 0
      do i=1,self%nlams
         self%nexptot = self%nexptot + self%nfourier(i)
         self%nexptotp = self%nexptotp + self%nphysical(i)
         if(self%nfourier(i).gt.self%nthmax) self%nthmax = self%nfourier(i)
         if(self%nphysical(i).gt.self%nphmax) self%nphmax = self%nphysical(i)
         nn = nn + self%nphysical(i)*self%nfourier(i)
      enddo
      



      allocate(self%fexpe(nn),self%fexpo(nn),self%fexpback(nn))
      allocate(self%tmp(self%nd,0:self%nmax,-self%nmax:self%nmax,self%nthd))
      allocate(self%mptmp(self%lmptemp,self%nthd))

      allocate(self%xshift(-5:5,self%nexptotp))
      allocate(self%yshift(-5:5,self%nexptotp))
      allocate(self%zshift(5,self%nexptotp))

      allocate(self%mexpf1(self%nd,self%nexptot,self%nthd),self%mexpf2(self%nd,self%nexptot,self%nthd), &
     &    self%mexpp1(self%nd,self%nexptotp,self%nthd))
      allocate(self%mexpp2(self%nd,self%nexptotp,self%nthd),self%mexppall(self%nd,self%nexptotp,16,self%nthd))
      !
      bigint = 0
      bigint = self%nboxes
      bigint = bigint*6
      bigint = bigint*self%nexptotp*self%nd

      if(self%ifprint.ge.1) print *, "mexp memory=",bigint/1.0d9, " GB "

    
      allocate(self%mexp(self%nd,self%nexptotp,self%nboxes,6),stat=iert)
      if(iert.ne.0) then
        print *, "Cannot allocate pw expansion workspace"
        print *, "bigint=", bigint
        self%ier = 8
        return
      endif

      allocate(self%list4ct(self%nboxes))
      allocate(self%ilist4(self%nboxes))
      do i=1,self%nboxes
        self%list4ct(i)=0
        self%ilist4(i)=0
      enddo
      self%cntlist4=0
      
      call mkexps(self%rlams,self%nlams,self%nphysical,self%nexptotp,self%xshift,self%yshift,self%zshift)
      call mkfexp(self%nlams,self%nfourier,self%nphysical,self%fexpe,self%fexpo,self%fexpback)


!
!      set scjsort
!
      do ilev=0,self%nlevels
        !$OMP TASKLOOP DEFAULT(SHARED) &
        !$OMP PRIVATE(ibox,nchild,istart,iend,i)
         do ibox=self%laddr(1,ilev),self%laddr(2,ilev)
            nchild = self%itree(self%ipointer(4)+ibox-1)
            if(nchild.gt.0) then
               istart = self%iexpcse(1,ibox)
               iend = self%iexpcse(2,ibox) 
               do i=istart,iend
                  self%scjsort(i) = self%scales(ilev)
               enddo
            endif
         enddo
        !$OMP END TASKLOOP
      enddo


      allocate( self%wlege(self%lw7) )
      call ylgndrfwini(self%nlege,self%wlege,self%lw7,self%lused7)

      do ilev=1,self%nlevels-1
         do ibox=self%laddr(1,ilev),self%laddr(2,ilev)
            if(self%nlist3(ibox).gt.0) then
              self%cntlist4=self%cntlist4+1
              self%list4ct(ibox)=self%cntlist4
              self%ilist4(self%cntlist4)=ibox
            endif
         enddo
      enddo

      self%lca = 4 * self%nmax


      allocate(self%pgboxwexp(self%nd,self%nexptotp,self%cntlist4,6))
      allocate(self%gboxmexp(self%nd*(self%nterms(self%nlevels)+1)*(2*self%nterms(self%nlevels)+1) &
                        ,8,self%cntlist4))
      allocate(self%gboxsubcenters(3,8,self%nthd)) 
      allocate(self%gboxfl(2,8,self%nthd))

      self%pgboxwexp=0d0
      self%gboxmexp=0d0

      nmaxt = 0 

      !!$OMP TASKLOOP DEFAULT(SHARED) PRIVATE(ibox,istart,iend,npts) &
      !!$OMP REDUCTION(max:nmaxt)
            do ibox=1,self%nboxes
              if(self%list4ct(ibox).gt.0) then
                istart = self%isrcse(1,ibox)
                iend = self%isrcse(2,ibox)
                npts = iend-istart+1
                if(npts.gt.nmaxt) nmaxt = npts
              endif
            enddo
      !!$OMP END TASKLOOP

      allocate(self%gboxind(nmaxt,self%nthd))
      allocate(self%gboxsort(3,nmaxt,self%nthd))
      allocate(self%gboxwexp(self%nd,self%nexptotp,6,8,self%nthd))
      allocate(self%gboxcgsort(self%nd,nmaxt,self%nthd))
      allocate(self%gboxdpsort(self%nd,3,nmaxt,self%nthd))




      allocate(self%uall(200,self%nthd), self%dall(200,self%nthd), self%nall(120,self%nthd))
      allocate(self%sall(120,self%nthd), self%eall(72,self%nthd), self%wall(72,self%nthd))
      allocate(self%u1234(36,self%nthd), self%d5678(36,self%nthd), self%n1256(24,self%nthd))
      allocate(self%s3478(24,self%nthd))
      allocate(self%e1357(16,self%nthd), self%w2468(16,self%nthd), self%n12(20,self%nthd))
      allocate(self%n56(20,self%nthd), self%s34(20,self%nthd), self%s78(20,self%nthd))
      allocate(self%e13(20,self%nthd), self%e57(20,self%nthd), self%w24(20,self%nthd), self%w68(20,self%nthd))
      allocate(self%e1(20,self%nthd), self%e3(5,self%nthd), self%e5(5,self%nthd), self%e7(5,self%nthd))
      allocate(self%w2(5,self%nthd), self%w4(5,self%nthd), self%w6(5,self%nthd), self%w8(5,self%nthd))
      allocate(self%iboxsubcenters(3,8,self%nthd))
      allocate(self%iboxfl(2,8,self%nthd))

    !  figure out allocations needed for iboxsrc,iboxsrcind,iboxpot

      nmaxt = 0
    !$OMP TASKLOOP DEFAULT(SHARED) PRIVATE(ibox,istart,iend,npts) &
    !$OMP REDUCTION(max:nmaxt)
          do ibox=1,self%nboxes
            if(self%nlist3(ibox).gt.0) then
              istart = self%isrcse(1,ibox)
              iend = self%isrcse(2,ibox)
              npts = iend-istart+1
              if(npts.gt.nmaxt) nmaxt = npts
            endif
          enddo
    !$OMP END TASKLOOP

      allocate(self%iboxsrcind(nmaxt,self%nthd))
      allocate(self%iboxsrc(3,nmaxt,self%nthd))
      allocate(self%iboxgrad(self%nd,3,nmaxt,self%nthd))
    end subroutine build2



    subroutine reset_expansion_coeff(self)
        !DEC$ ATTRIBUTES ALIAS:"fmm3d_tree_mod_mp_reset_expansion_coeff_" :: reset_expansion_coeff
        class(FMM3DTree), intent(inout) :: self
        !------------------------------------------------
        integer :: ilev, ibox
        integer :: i,j,k,idim
        !------------------------------------------------

      do ilev = 0,self%nlevels
        !$OMP TASKLOOP DEFAULT(SHARED) &
        !$OMP PRIVATE(ibox)
        do ibox=self%laddr(1,ilev),self%laddr(2,ilev)
          call mpzero(self%nd,self%rmlexp(self%iaddr(1,ibox)),self%nterms(ilev))
          call mpzero(self%nd,self%rmlexp(self%iaddr(2,ibox)),self%nterms(ilev))
        enddo
        !$OMP END TASKLOOP
      enddo
      !-------------------------------------------
      !$OMP TASKLOOP collapse(4) DEFAULT(SHARED) &
      !$OMP PRIVATE(i,j,k,idim)
            do k=1,6
              do i=1,self%nboxes
                do j=1,self%nexptotp
                  do idim=1,self%nd
                    self%mexp(idim,j,i,k) = 0.0d0
                  enddo
                enddo
              enddo
            enddo
      !$OMP END TASKLOOP

    end subroutine reset_expansion_coeff




    subroutine reset_sort_arg(self)
        !DEC$ ATTRIBUTES ALIAS:"fmm3d_tree_mod_mp_reset_sort_arg_" :: reset_sort_arg
        class(FMM3DTree), intent(inout) :: self
        !------------------------------------------------
        integer i,idim
        !------------------------------------------------
    !$OMP TASKLOOP collapse(2) DEFAULT(SHARED) PRIVATE(i,idim)
            do i=1,self%nsource
            do idim=1,self%nd
                self%gradsort(idim,1,i) = 0.0
                self%gradsort(idim,2,i) = 0.0
                self%gradsort(idim,3,i) = 0.0
            enddo
            enddo
    !$OMP END TASKLOOP
    end subroutine reset_sort_arg


    subroutine dealloc(self)
        !DEC$ ATTRIBUTES ALIAS:"fmm3d_tree_mod_mp_dealloc_" :: dealloc
        class(FMM3DTree), intent(inout) :: self
        !------------------------------------------------


        !The source positions are handed over by build_tree and referenced by the tree for as long
        !as it lives, so the tree releases them here rather than the caller releasing them while
        !self%source still points at them.
        if (associated(self%source)) deallocate(self%source)
        if (associated(self%itree)) deallocate(self%itree)
        if (associated(self%boxsize)) deallocate(self%boxsize)
        if (associated(self%treecenters)) deallocate(self%treecenters)
        if (associated(self%isrcse)) deallocate(self%isrcse)
        if (associated(self%iexpcse)) deallocate(self%iexpcse)
        if (associated(self%isrc)) deallocate(self%isrc)
        if (associated(self%sourcesort)) deallocate(self%sourcesort)
        if (associated(self%dipvecsort)) deallocate(self%dipvecsort)
        if (associated(self%gradsort)) deallocate(self%gradsort)
        if (associated(self%nterms)) deallocate(self%nterms)
        if (associated(self%iaddr)) deallocate(self%iaddr)
        if (associated(self%mptemp)) deallocate(self%mptemp)
        if (associated(self%mptemp2)) deallocate(self%mptemp2)
        if (associated(self%rmlexp)) deallocate(self%rmlexp)    
    end subroutine dealloc



     subroutine lfmm3dmain_tree(self) 
      implicit none
        !DEC$ ATTRIBUTES ALIAS:"fmm3d_tree_mod_mp_lfmm3dmain_tree_" :: lfmm3dmain_tree
        class(FMM3DTree), intent(inout) :: self
      integer nd
      integer ier
      double precision eps
      integer nsource,ntarg,nexpc
      integer ndiv,nlevels

      double precision, contiguous, pointer :: sourcesort(:,:)
      double precision, contiguous, pointer :: dipvecsort(:,:,:)


      double precision, contiguous, pointer ::  pot(:,:),grad(:,:,:),hess(:,:,:)
      double precision, contiguous, pointer ::  pottarg(:,:),gradtarg(:,:,:),hesstarg(:,:,:)

      integer ntj
      integer ifnear
      double precision expcsort(3,self%nexpc)
      double complex, contiguous, pointer ::  tsort(:,:,:,:)
      double precision scjsort(self%nexpc)

      integer nboxes
      integer(kind=8) lmptot
      integer(kind=8), contiguous, pointer ::  iaddr(:,:)



      integer lmptemp
      double precision, contiguous, pointer ::  rmlexp(:)
      double precision, contiguous, pointer ::  mptemp(:)
      double precision, contiguous, pointer ::  mptemp2(:)

      double precision thresh
       
      double precision timeinfo(6)
      double precision, contiguous, pointer :: centers(:,:)

      integer isep,iper
      integer, contiguous, pointer ::  laddr(:,:)
      integer, contiguous, pointer ::  nterms(:)
      integer(kind=8) ipointer(8),ltree
      integer, contiguous, pointer ::  itree(:)
      double precision, contiguous, pointer :: rscales(:)
      double precision, contiguous, pointer :: boxsize(:)
      integer, contiguous, pointer :: isrcse(:,:),iexpcse(:,:)

      integer, contiguous, pointer :: nlist1(:),list1(:,:)
      integer, contiguous, pointer :: nlist2(:),list2(:,:)
      integer, contiguous, pointer :: nlist3(:),list3(:,:)
      integer, contiguous, pointer :: nlist4(:),list4(:,:)

      integer nuall,ndall,nnall,nsall,neall,nwall
      integer nu1234,nd5678,nn1256,ns3478,ne1357,nw2468
      integer nn12,nn56,ns34,ns78,ne13,ne57,nw24,nw68
      integer ne1,ne3,ne5,ne7,nw2,nw4,nw6,nw8

      integer, contiguous, pointer :: uall(:,:),dall(:,:),nall(:,:)
      integer, contiguous, pointer :: sall(:,:),eall(:,:),wall(:,:)
      integer, contiguous, pointer :: u1234(:,:),d5678(:,:)
      integer, contiguous, pointer :: n1256(:,:),s3478(:,:)
      integer, contiguous, pointer :: e1357(:,:),w2468(:,:)
      integer, contiguous, pointer :: n12(:,:),n56(:,:),s34(:,:),s78(:,:)
      integer, contiguous, pointer :: e13(:,:),e57(:,:),w24(:,:),w68(:,:)
      integer, contiguous, pointer :: e1(:,:),e3(:,:),e5(:,:),e7(:,:)
      integer, contiguous, pointer :: w2(:,:),w4(:,:),w6(:,:),w8(:,:)

!     temp variables
      integer i,j,k,l,ii,jj,kk,ll,m,idim,igbox
      integer ibox,jbox,ilev,npts,npts0,kbox,dir
      integer nchild

      integer istart,iend,istarts,iends
      integer istartt,iendt,istarte,iende
      integer isstart,isend,jsstart,jsend
      integer jstart,jend

      integer ifprint

      double precision d,time1,time2,second,omp_get_wtime
      double precision pottmp,fldtmp(3),hesstmp(3)

!     PW variables
      integer nexpmax, nlams, nmax, nthmax, nphmax,nmax2,nmaxt
      integer lca
      double precision, contiguous, pointer :: carray(:,:), dc(:,:)
      double precision, contiguous, pointer :: cs(:,:),fact(:),rdplus(:,:,:)
      double precision, contiguous, pointer :: rdminus(:,:,:), rdsq3(:,:,:)
      double precision, contiguous, pointer :: rdmsq3(:,:,:)
  
        double precision, contiguous, pointer :: rlams(:),whts(:)

      double precision, contiguous, pointer :: rlsc(:,:,:)
      integer, contiguous, pointer :: nfourier(:), nphysical(:)
      integer nexptot, nexptotp
      double complex, contiguous, pointer:: xshift(:,:)
      double complex, contiguous, pointer :: yshift(:,:)
      double precision, contiguous, pointer :: zshift(:,:)

      double complex, contiguous, pointer :: fexpe(:),fexpo(:),fexpback(:)
      double complex, contiguous, pointer :: mexp(:,:,:,:)
      double complex, contiguous, pointer :: mexpf1(:,:,:),mexpf2(:,:,:)
      double complex, contiguous, pointer :: &
     &    mexpp1(:,:,:),mexpp2(:,:,:),mexppall(:,:,:,:)

      double complex, contiguous, pointer:: tmp(:,:,:,:)
      double precision, contiguous, pointer :: mptmp(:,:)

      double precision sourcetmp(3)
      double complex chargetmp

      integer ix,iy,iz,ictr
      double precision rtmp
      double complex zmul

      integer nlege, lw7, lused7, itype

      double precision, contiguous, pointer :: wlege(:)
      integer nterms_eval(4,0:self%nlevels)

      integer mnlist1, mnlist2,mnlist3,mnlist4,mnbors
      double complex eye, ztmp
      double precision alphaj
      integer ctr,nn,iptr1,iptr2
      double precision, contiguous, pointer :: rscpow(:)
      double precision pi,errtmp
      double complex ima

      double precision ctmp(3)

!     list 3 variables
      double complex, contiguous, pointer :: iboxlexp(:,:,:)
      double precision, contiguous, pointer :: iboxsubcenters(:,:,:)
      double precision, contiguous, pointer :: iboxpot(:,:,:)
      double precision, contiguous, pointer :: iboxgrad(:,:,:,:)
      double precision, contiguous, pointer :: iboxhess(:,:,:,:)
      double precision, contiguous, pointer :: iboxsrc(:,:,:)
      integer, contiguous, pointer :: iboxsrcind(:,:)
      integer, contiguous, pointer :: iboxfl(:,:,:)
!     end of list 3 variables
!     list 4 variables
      integer cntlist4
      integer, contiguous, pointer :: list4ct(:),ilist4(:)
      double complex, contiguous, pointer :: gboxmexp(:,:,:)
      double complex, contiguous, pointer :: gboxwexp(:,:,:,:,:)
      double complex, contiguous, pointer :: pgboxwexp(:,:,:,:)
      double precision, contiguous, pointer :: gboxsubcenters(:,:,:)
      double precision, contiguous, pointer :: gboxsort(:,:,:)
      integer, contiguous, pointer :: gboxind(:,:)
      integer, contiguous, pointer :: gboxfl(:,:,:)
      double precision, contiguous, pointer :: gboxcgsort(:,:,:)
      double precision, contiguous, pointer :: gboxdpsort(:,:,:,:)
!     end of list 4 variables

!
!   hessian variables
!
      double precision, contiguous, pointer :: scarray(:,:)
      integer *8 bigint
      integer iert
      data ima/(0.0d0,1.0d0)/
      integer nthd,ithd
      integer, save :: itimer = 0

      call trace%begin("FMM3DTree:lfmm3dmain_tree", itimer=itimer, verbose=3)
      !--------
      nthd = self%nthd
      ifprint=0

      !--------- load self variables... replacement of argument list -----
      nd = self%nd
      eps = self%eps
      nsource = self%nsource
      sourcesort => self%sourcesort
      dipvecsort => self%dipvecsort

      ntarg = self%ntarg
      nexpc = self%nexpc
      expcsort(:,1) = self%expcsort(:)



      iaddr => self%iaddr
      rmlexp => self%rmlexp
      lmptot = self%lmptot
      mptemp => self%mptemp
      mptemp2 => self%mptemp2
      lmptemp = self%lmptemp

      itree => self%itree
      ltree = self%ltree
      ipointer = self%ipointer
      ndiv = self%ndiv
      nlevels = self%nlevels

      nboxes = self%nboxes
      iper = self%iper
      boxsize => self%boxsize
      centers => self%treecenters
      isrcse => self%isrcse
      iexpcse => self%iexpcse

      rscales => self%scales 
      laddr => self%laddr
      nterms => self%nterms

      grad => self%gradsort
      ntj = self%ntj

      scjsort = self%scjsort
      ifnear = self%ifnear
      timeinfo = self%timeinfo
      ier = self%ier

      !--------------------------------------------------------------------

      !--------------- load self variables into local variables -------
        thresh = self%thresh
        mnlist1 = self%mnlist1
        mnlist2 = self%mnlist2
        mnlist3 = self%mnlist3
        mnlist4 = self%mnlist4
        mnbors  = self%mnbors

        isep = self%isep
      
        list1 => self%list1
        nlist1 => self%nlist1
        list2 => self%list2
        nlist2 => self%nlist2
        list3 => self%list3
        nlist3 => self%nlist3
        list4 => self%list4
        nlist4 => self%nlist4
        nlams = self%nlams
        rlams => self%rlams
        whts => self%whts
        nphysical => self%nphysical
        nfourier => self%nfourier



        rscpow => self%rscpow
        carray => self%carray
        dc => self%dc
        rdplus => self%rdplus
        rdminus => self%rdminus
        rdsq3 => self%rdsq3
        rdmsq3 => self%rdmsq3
        rlsc => self%rlsc
         nmax = self%nmax
         nmaxt = nmax


         scarray => self%scarray
         
        fexpe => self%fexpe
        fexpo => self%fexpo
        fexpback => self%fexpback
        tmp => self%tmp
        mptmp => self%mptmp
        xshift => self%xshift
        yshift => self%yshift
        zshift => self%zshift
        mexpf1 => self%mexpf1
        mexpf2 => self%mexpf2
        mexpp1 => self%mexpp1
        mexpp2 => self%mexpp2
        mexppall => self%mexppall

        nexptotp = self%nexptotp
        nexptot = self%nexptot
        nthmax = self%nthmax
        nphmax = self%nphmax


        mexp => self%mexp
        list4ct => self%list4ct
        ilist4 => self%ilist4
        cntlist4 = self%cntlist4


        wlege => self%wlege

        lca = self%lca

        pgboxwexp => self%pgboxwexp
        gboxmexp => self%gboxmexp
        pgboxwexp=0d0
        gboxmexp=0d0
        
        gboxsubcenters => self%gboxsubcenters
        gboxfl => self%gboxfl
        gboxind => self%gboxind
        gboxsort => self%gboxsort
        gboxwexp => self%gboxwexp
        gboxcgsort => self%gboxcgsort
        gboxdpsort => self%gboxdpsort

        uall => self%uall
        dall => self%dall
        nall => self%nall
        sall => self%sall
        eall => self%eall
        wall => self%wall
        u1234 => self%u1234
        d5678 => self%d5678
        n1256 => self%n1256
        s3478 => self%s3478
        e1357 => self%e1357
        w2468 => self%w2468
        n12 => self%n12
        n56 => self%n56
        s34 => self%s34
        s78 => self%s78
        e13 => self%e13
        e57 => self%e57
        w24 => self%w24
        w68 => self%w68
        e1 => self%e1
        e3 => self%e3
        e5 => self%e5
        e7 => self%e7
        w2 => self%w2
        w4 => self%w4
        w6 => self%w6
        w8 => self%w8
        iboxsubcenters => self%iboxsubcenters
        iboxfl => self%iboxfl


        iboxsrcind => self%iboxsrcind
        iboxsrc => self%iboxsrc
        iboxgrad => self%iboxgrad

        nlege = self%nlege

        self%gradsort = 0.0
        
        self%rmlexp = 0.0
        self%mexp = 0.0


    !-----------------------------------------------

      call self%reorder_dipvec()

!
!====== not tested and only for non-complete trees
!

!     form mexp for all list4 type box at first ghost box center
!      do ilev=1,nlevels-1
!
!         rscpow(0) = 1.0d0/boxsize(ilev+1)
!         rtmp = rscales(ilev+1)/boxsize(ilev+1)
!         do i=1,nterms(ilev+1)
!            rscpow(i) = rscpow(i-1)*rtmp
!         enddo
!
!!!$OMP TASKLOOP DEFAULT(SHARED) &
!!!$OMP PRIVATE(ibox,istart,iend,jbox,jstart,jend,npts,npts0,i) &
!!!$OMP PRIVATE(ithd)
!         do ibox=laddr(1,ilev),laddr(2,ilev)
!            ithd = 0
!           !ithd=omp_get_thread_num()
!           ! ithd = omp%thread_id()
!            ithd = ithd + 1
!            if(list4ct(ibox).gt.0) then
!              istart=isrcse(1,ibox)
!              iend=isrcse(2,ibox)
!              npts = iend-istart+1
!
!              if(npts.gt.0) then
!                call subdividebox(sourcesort(1,istart),npts, &
!     &    centers(1,ibox),boxsize(ilev+1), &
!     &    gboxind(1,ithd),gboxfl(1,1,ithd), &
!     &    gboxsubcenters(1,1,ithd))
!                call dreorderf(3,npts,sourcesort(1,istart), &
!     &    gboxsort(1,1,ithd),gboxind(1,ithd))
!                  call dreorderf(3*nd,npts,dipvecsort(1,1,istart), &
!     &    gboxdpsort(1,1,1,ithd),gboxind(1,ithd))
!                do i=1,8
!                  if(gboxfl(1,i,ithd).gt.0) then
!                    jstart=gboxfl(1,i,ithd)
!                    jend=gboxfl(2,i,ithd)
!                    npts0=jend-jstart+1
!                    jbox=list4ct(ibox)
!
!                      call l3dformmpd(nd,rscales(ilev+1), &
!     &    gboxsort(1,jstart,ithd), &
!     &    gboxdpsort(1,1,jstart,ithd), &
!     &    npts0,gboxsubcenters(1,i,ithd),nterms(ilev+1), &
!     &    gboxmexp(1,i,jbox),wlege,nlege)          
!
!                    call l3dmpmp(nd,rscales(ilev+1), &
!     &    gboxsubcenters(1,i,ithd),gboxmexp(1,i,jbox), &
!     &    nterms(ilev+1),rscales(ilev),centers(1,ibox), &
!     &    rmlexp(iaddr(1,ibox)),nterms(ilev),dc,lca)
!     
!                    call mpscale(nd,nterms(ilev+1),gboxmexp(1,i,jbox), &
!     &    rscpow,tmp(1,0,-nmax,ithd))
!!
!!c                process up down for current box
!!
!                    call mpoletoexp(nd,tmp(1,0,-nmax,ithd), &
!     &    nterms(ilev+1),nlams, &
!     &    nfourier,nexptot,mexpf1(1,1,ithd), &
!     &    mexpf2(1,1,ithd),rlsc)
!
!                    call ftophys(nd,mexpf1(1,1,ithd), &
!     &    nlams,rlams,nfourier, &
!     &    nphysical,nthmax,gboxwexp(1,1,1,i,ithd), &
!     &    fexpe,fexpo)
!
!                    call ftophys(nd,mexpf2(1,1,ithd), &
!     &    nlams,rlams,nfourier, &
!     &    nphysical,nthmax,gboxwexp(1,1,2,i,ithd), &
!     &    fexpe,fexpo)
!
!                    call processgboxudexp(nd,gboxwexp(1,1,1,i,ithd), &
!     &    gboxwexp(1,1,2,i,ithd),i,nexptotp, &
!     &    pgboxwexp(1,1,jbox,1),pgboxwexp(1,1,jbox,2), &
!     &    xshift,yshift,zshift)
!!
!!c                process north-south for current box
!!
!                    call rotztoy(nd,nterms(ilev+1),tmp(1,0,-nmax,ithd), &
!     &    mptmp(1,ithd),rdminus)
!                    call mpoletoexp(nd,mptmp(1,ithd), &
!     &    nterms(ilev+1),nlams, &
!     &    nfourier,nexptot,mexpf1(1,1,ithd), &
!     &    mexpf2(1,1,ithd),rlsc)
!
!                    call ftophys(nd,mexpf1(1,1,ithd), &
!     &    nlams,rlams,nfourier, &
!     &    nphysical,nthmax,gboxwexp(1,1,3,i,ithd), &
!     &    fexpe,fexpo)
!
!                    call ftophys(nd,mexpf2(1,1,ithd), &
!     &    nlams,rlams,nfourier, &
!     &    nphysical,nthmax,gboxwexp(1,1,4,i,ithd), &
!     &    fexpe,fexpo)
!
!                    call processgboxnsexp(nd,gboxwexp(1,1,3,i,ithd), &
!     &    gboxwexp(1,1,4,i,ithd),i,nexptotp, &
!     &    pgboxwexp(1,1,jbox,3),pgboxwexp(1,1,jbox,4), &
!     &    xshift,yshift,zshift)
!!
!!c                process east-west for current box
!!
!                    call rotztox(nd,nterms(ilev+1),tmp(1,0,-nmax,ithd), &
!     &    mptmp(1,ithd),rdplus)
!                    call mpoletoexp(nd,mptmp(1,ithd), &
!     &    nterms(ilev+1),nlams, &
!     &    nfourier,nexptot,mexpf1(1,1,ithd), &
!     &    mexpf2(1,1,ithd),rlsc)
!
!                    call ftophys(nd,mexpf1(1,1,ithd), &
!     &    nlams,rlams,nfourier, &
!     &    nphysical,nthmax,gboxwexp(1,1,5,i,ithd), &
!     &    fexpe,fexpo)
!
!                    call ftophys(nd,mexpf2(1,1,ithd), &
!     &    nlams,rlams,nfourier, &
!     &    nphysical,nthmax,gboxwexp(1,1,6,i,ithd), &
!     &    fexpe,fexpo)
!                
!                    call processgboxewexp(nd,gboxwexp(1,1,5,i,ithd), &
!     &    gboxwexp(1,1,6,i,ithd),i,nexptotp, &
!     &    pgboxwexp(1,1,jbox,5),pgboxwexp(1,1,jbox,6), &
!     &    xshift,yshift,zshift)
!                  endif
!                enddo
!              endif
!            endif
!         enddo
!!!$OMP END TASKLOOP
!      enddo


!------------------ step 1 ??? -----------------------------------------------------------------
!       ... step 1, locate all charges, assign them to boxes, and
!       form multipole expansions
      do ilev=2,nlevels
!$OMP TASKLOOP DEFAULT(NONE) &
!$OMP PRIVATE(ibox,npts,istart,iend,nchild) &
!$OMP SHARED(ilev, laddr, isrcse, itree, ipointer, list4ct, nd, rscales, sourcesort, dipvecsort, centers, nterms, rmlexp, iaddr, wlege, nlege)
            do ibox=laddr(1,ilev),laddr(2,ilev)

               istart = isrcse(1,ibox) 
               iend = isrcse(2,ibox) 
               npts = iend-istart+1

               nchild = itree(ipointer(4)+ibox-1)

               if(npts.gt.0.and.nchild.eq.0.and.list4ct(ibox).eq.0) then      
                  call l3dformmpd_new(nd,rscales(ilev), &
     &    sourcesort(1,istart), &
     &    dipvecsort(1,1,istart),npts, &
     &    centers(1,ibox),nterms(ilev), &
     &    rmlexp(iaddr(1,ibox)),wlege,nlege)    
               endif
            enddo
!$OMP END TASKLOOP
      enddo

      !----------------------------------------------------------------------------------------------------

      do ilev=nlevels-1,0,-1
!$OMP TASKLOOP DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,i,jbox,istart,iend,npts)
         do ibox = laddr(1,ilev),laddr(2,ilev)
            do i=1,8
               jbox = itree(ipointer(5)+8*(ibox-1)+i-1)
               if(jbox.gt.0) then
                  istart = isrcse(1,jbox)
                  iend = isrcse(2,jbox)
                  npts = iend-istart+1
                  if(npts.gt.0) then
                     call l3dmpmp(nd,rscales(ilev+1), &
     &    centers(1,jbox),rmlexp(iaddr(1,jbox)), &
     &    nterms(ilev+1),rscales(ilev),centers(1,ibox), &
     &    rmlexp(iaddr(1,ibox)),nterms(ilev),dc,lca)
                  endif
               endif
            enddo
         enddo
!$OMP END TASKLOOP
      enddo


!-----------

      do ilev=2,nlevels
        allocate(iboxlexp(nd*(nterms(ilev)+1)* &
     &    (2*nterms(ilev)+1),8,nthd))

         rscpow(0) = 1.0d0/boxsize(ilev)
         rtmp = rscales(ilev)/boxsize(ilev)
         do i=1,nterms(ilev)
            rscpow(i) = rscpow(i-1)*rtmp
         enddo

!$OMP TASKLOOP DEFAULT (SHARED) &
!$OMP PRIVATE(ibox,istart,iend,npts) &
!$OMP PRIVATE(ithd)
         do ibox=laddr(1,ilev),laddr(2,ilev)
            !ithd = 0
            ithd=omp%thread_id()
            ithd = ithd + 1
            istart = isrcse(1,ibox) 
            iend = isrcse(2,ibox)

            npts = iend-istart+1

            if(npts.gt.0) then
!            rescale the multipole expansion

                call mpscale(nd,nterms(ilev),rmlexp(iaddr(1,ibox)), &
     &    rscpow,tmp(1,0,-nmax,ithd))
!
!c                process up down for current box
!
                call mpoletoexp(nd,tmp(1,0,-nmax,ithd),nterms(ilev), &
     &    nlams,nfourier, &
     &    nexptot,mexpf1(1,1,ithd),mexpf2(1,1,ithd),rlsc)

                call ftophys(nd,mexpf1(1,1,ithd),nlams,rlams,nfourier, &
     &    nphysical,nthmax,mexp(1,1,ibox,1),fexpe,fexpo)

                call ftophys(nd,mexpf2(1,1,ithd),nlams,rlams,nfourier, &
     &    nphysical,nthmax,mexp(1,1,ibox,2),fexpe,fexpo)


!
!c                process north-south for current box
!
                call rotztoy(nd,nterms(ilev),tmp(1,0,-nmax,ithd), &
     &    mptmp(1,ithd),rdminus)
                call mpoletoexp(nd,mptmp(1,ithd),nterms(ilev), &
     &    nlams,nfourier, &
     &    nexptot,mexpf1(1,1,ithd),mexpf2(1,1,ithd),rlsc)

                call ftophys(nd,mexpf1(1,1,ithd),nlams,rlams,nfourier, &
     &    nphysical,nthmax,mexp(1,1,ibox,3),fexpe,fexpo)

                call ftophys(nd,mexpf2(1,1,ithd),nlams,rlams,nfourier, &
     &    nphysical,nthmax,mexp(1,1,ibox,4),fexpe,fexpo)

!
!c                process east-west for current box

                call rotztox(nd,nterms(ilev),tmp(1,0,-nmax,ithd), &
     &    mptmp(1,ithd),rdplus)
                call mpoletoexp(nd,mptmp(1,ithd), &
     &    nterms(ilev),nlams,nfourier, &
     &    nexptot,mexpf1(1,1,ithd), &
     &    mexpf2(1,1,ithd),rlsc)

                call ftophys(nd,mexpf1(1,1,ithd),nlams,rlams,nfourier, &
     &    nphysical,nthmax,mexp(1,1,ibox,5),fexpe,fexpo)


                call ftophys(nd,mexpf2(1,1,ithd),nlams,rlams,nfourier, &
     &    nphysical,nthmax,mexp(1,1,ibox,6),fexpe,fexpo)

            endif

         enddo
!$OMP END TASKLOOP
!
!
!c         loop over parent boxes and ship plane wave
!          expansions to the first child of parent 
!          boxes. 
!          The codes are now written from a gathering perspective
!
!          so the first child of the parent is the one
!          recieving all the local expansions
!          coming from all the lists
!
!          
!
         rscpow(0) = 1.0d0
         rtmp = rscales(ilev)/boxsize(ilev)
         do i=1,nterms(ilev)
            rscpow(i) = rscpow(i-1)*rtmp
         enddo
!$OMP TASKLOOP DEFAULT (SHARED) &
!$OMP PRIVATE(ibox,istart,iend,npts,nchild) &
!$OMP PRIVATE(nuall,ndall,nnall,nsall) &
!$OMP PRIVATE(neall,nwall,nu1234,nd5678) &
!$OMP PRIVATE(nn1256,ns3478,ne1357,nw2468) &
!$OMP PRIVATE(nn12,nn56,ns34,ns78,ne13,ne57) &
!$OMP PRIVATE(nw24,nw68,ne1,ne3,ne5,ne7) &
!$OMP PRIVATE(nw2,nw4,nw6,nw8) &
!$OMP PRIVATE(npts0,ctmp,jstart,jend,i) &
!$OMP PRIVATE(ithd)
         do ibox = laddr(1,ilev-1),laddr(2,ilev-1)
           !ithd = 0
           ithd=omp%thread_id()
           ithd = ithd + 1
           npts = 0

           istart = iexpcse(1,ibox) 
           iend = iexpcse(2,ibox) 
           npts = npts + iend-istart+1

           nchild = itree(ipointer(4)+ibox-1)

          istart = isrcse(1,ibox) 
          iend = isrcse(2,ibox) 
          npts = npts + iend-istart+1
           


           if(npts.gt.0.and.nchild.gt.0) then

               call getpwlistall(ibox,boxsize(ilev),nboxes, &
     &    itree(ipointer(6)+ibox-1),itree(ipointer(7)+ &
     &    mnbors*(ibox-1)),nchild,itree(ipointer(5)),centers, &
     &    isep,nuall,uall(1,ithd),ndall,dall(1,ithd), &
     &    nnall,nall(1,ithd),nsall,sall(1,ithd), &
     &    neall,eall(1,ithd),nwall,wall(1,ithd), &
     &    nu1234,u1234(1,ithd),nd5678,d5678(1,ithd), &
     &    nn1256,n1256(1,ithd),ns3478,s3478(1,ithd), &
     &    ne1357,e1357(1,ithd),nw2468,w2468(1,ithd), &
     &    nn12,n12(1,ithd),nn56,n56(1,ithd),ns34,s34(1,ithd), &
     &    ns78,s78(1,ithd),ne13,e13(1,ithd),ne57,e57(1,ithd), &
     &    nw24,w24(1,ithd),nw68,w68(1,ithd),ne1,e1(1,ithd), &
     &    ne3,e3(1,ithd),ne5,e5(1,ithd),ne7,e7(1,ithd), &
     &    nw2,w2(1,ithd),nw4,w4(1,ithd),nw6,w6(1,ithd), &
     &    nw8,w8(1,ithd))

               call processudexp(nd,ibox,ilev,nboxes,centers, &
     &    itree(ipointer(5)),rscales(ilev),boxsize(ilev), &
     &    nterms(ilev), &
     &    iaddr,rmlexp,rlams,whts, &
     &    nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp, &
     &    nuall,uall(1,ithd),nu1234,u1234(1,ithd), &
     &    ndall,dall(1,ithd),nd5678,d5678(1,ithd), &
     &    mexpf1(1,1,ithd),mexpf2(1,1,ithd), &
     &    mexpp1(1,1,ithd),mexpp2(1,1,ithd),mexppall(1,1,1,ithd), &
     &    mexppall(1,1,2,ithd),mexppall(1,1,3,ithd), &
     &    mexppall(1,1,4,ithd),xshift, &
     &    yshift,zshift,fexpback,rlsc,rscpow, &
     &    pgboxwexp,cntlist4,list4ct, &
     &    nlist4,list4,mnlist4)
               
               call processnsexp(nd,ibox,ilev,nboxes,centers, &
     &    itree(ipointer(5)),rscales(ilev),boxsize(ilev), &
     &    nterms(ilev), &
     &    iaddr,rmlexp,rlams,whts, &
     &    nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp, &
     &    nnall,nall(1,ithd),nn1256,n1256(1,ithd), &
     &    nn12,n12(1,ithd),nn56,n56(1,ithd),nsall,sall(1,ithd), &
     &    ns3478,s3478(1,ithd),ns34,s34(1,ithd),ns78,s78(1,ithd), &
     &    mexpf1(1,1,ithd),mexpf2(1,1,ithd), &
     &    mexpp1(1,1,ithd),mexpp2(1,1,ithd),mexppall(1,1,1,ithd), &
     &    mexppall(1,1,2,ithd),mexppall(1,1,3,ithd), &
     &    mexppall(1,1,4,ithd), &
     &    mexppall(1,1,5,ithd),mexppall(1,1,6,ithd), &
     &    mexppall(1,1,7,ithd), &
     &    mexppall(1,1,8,ithd),rdplus,xshift,yshift,zshift, &
     &    fexpback,rlsc,rscpow, &
     &    pgboxwexp,cntlist4,list4ct, &
     &    nlist4,list4,mnlist4)

               
               call processewexp(nd,ibox,ilev,nboxes,centers, &
     &    itree(ipointer(5)),rscales(ilev),boxsize(ilev), &
     &    nterms(ilev), &
     &    iaddr,rmlexp,rlams,whts, &
     &    nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp, &
     &    neall,eall(1,ithd),ne1357,e1357(1,ithd), &
     &    ne13,e13(1,ithd),ne57,e57(1,ithd),ne1,e1(1,ithd), &
     &    ne3,e3(1,ithd),ne5,e5(1,ithd), &
     &    ne7,e7(1,ithd),nwall,wall(1,ithd), &
     &    nw2468,w2468(1,ithd), &
     &    nw24,w24(1,ithd),nw68,w68(1,ithd), &
     &    nw2,w2(1,ithd),nw4,w4(1,ithd),nw6,w6(1,ithd), &
     &    nw8,w8(1,ithd), &
     &    mexpf1(1,1,ithd),mexpf2(1,1,ithd), &
     &    mexpp1(1,1,ithd),mexpp2(1,1,ithd),mexppall(1,1,1,ithd), &
     &    mexppall(1,1,2,ithd),mexppall(1,1,3,ithd), &
     &    mexppall(1,1,4,ithd), &
     &    mexppall(1,1,5,ithd),mexppall(1,1,6,ithd), &
     &    mexppall(1,1,7,ithd),mexppall(1,1,8,ithd), &
     &    mexppall(1,1,9,ithd), &
     &    mexppall(1,1,10,ithd),mexppall(1,1,11,ithd), &
     &    mexppall(1,1,12,ithd), &
     &    mexppall(1,1,13,ithd),mexppall(1,1,14,ithd), &
     &    mexppall(1,1,15,ithd), &
     &    mexppall(1,1,16,ithd),rdminus,xshift,yshift,zshift, &
     &    fexpback,rlsc,rscpow, &
     &    pgboxwexp,cntlist4,list4ct,nlist4,list4,mnlist4)


            endif

            if(nlist3(ibox).gt.0.and.npts.gt.0) then

              call getlist3pwlistall(ibox,boxsize(ilev),nboxes, &
     &    nlist3(ibox),list3(1,ibox),isep, &
     &    centers,nuall,uall(1,ithd),ndall,dall(1,ithd), &
     &    nnall,nall(1,ithd), &
     &    nsall,sall(1,ithd),neall,eall(1,ithd), &
     &    nwall,wall(1,ithd))
              iboxlexp(:,:,ithd)=0.0d0

              call processlist3udexplong(nd,ibox,nboxes,centers, &
     &    boxsize(ilev),nterms(ilev),iboxlexp(1,1,ithd),rlams, &
     &    whts,nlams,nfourier,nphysical,nthmax,nexptot, &
     &    nexptotp,mexp,nuall,uall(1,ithd),ndall,dall(1,ithd), &
     &    mexpf1(1,1,ithd),mexpf2(1,1,ithd), &
     &    mexpp1(1,1,ithd),mexpp2(1,1,ithd), &
     &    mexppall(1,1,1,ithd),mexppall(1,1,2,ithd), &
     &    xshift,yshift,zshift,fexpback,rlsc,rscpow)

              call processlist3nsexplong(nd,ibox,nboxes,centers, &
     &    boxsize(ilev),nterms(ilev),iboxlexp(1,1,ithd),rlams, &
     &    whts,nlams,nfourier,nphysical,nthmax,nexptot, &
     &    nexptotp,mexp,nnall,nall(1,ithd),nsall,sall(1,ithd), &
     &    mexpf1(1,1,ithd),mexpf2(1,1,ithd), &
     &    mexpp1(1,1,ithd),mexpp2(1,1,ithd), &
     &    mexppall(1,1,1,ithd),mexppall(1,1,2,ithd),rdplus, &
     &    xshift,yshift,zshift,fexpback,rlsc,rscpow)

              call processlist3ewexplong(nd,ibox,nboxes,centers, &
     &    boxsize(ilev),nterms(ilev),iboxlexp(1,1,ithd),rlams, &
     &    whts,nlams,nfourier,nphysical,nthmax,nexptot, &
     &    nexptotp,mexp,neall,eall(1,ithd),nwall,wall(1,ithd), &
     &    mexpf1(1,1,ithd),mexpf2(1,1,ithd), &
     &    mexpp1(1,1,ithd),mexpp2(1,1,ithd), &
     &    mexppall(1,1,1,ithd),mexppall(1,1,2,ithd),rdminus, &
     &    xshift,yshift,zshift,fexpback,rlsc,rscpow)


                istart = isrcse(1,ibox)
                iend = isrcse(2,ibox) 
                npts = iend-istart+1
                if(npts.gt.0) then
                  call subdividebox(sourcesort(1,istart),npts, &
     &    centers(1,ibox),boxsize(ilev), &
     &    iboxsrcind(1,ithd),iboxfl(1,1,ithd), &
     &    iboxsubcenters(1,1,ithd))
                  call dreorderf(3,npts,sourcesort(1,istart), &
     &    iboxsrc(1,1,ithd),iboxsrcind(1,ithd))
                  call dreorderf(3*nd,npts,grad(1,1,istart), &
     &    iboxgrad(1,1,1,ithd),iboxsrcind(1,ithd))
                  do i=1,8
                    if(iboxfl(1,i,ithd).gt.0) then
                      jstart=iboxfl(1,i,ithd)
                      jend=iboxfl(2,i,ithd)
                      npts0=jend-jstart+1
                      if(npts0.gt.0) then

                        call l3dtaevalg_grad(nd,rscales(ilev), &
     &    iboxsubcenters(1,i,ithd),iboxlexp(1,i,ithd), &
     &    nterms(ilev),iboxsrc(1,jstart,ithd),npts0, &
     &    iboxgrad(1,1,jstart,ithd),wlege,nlege)
                      endif
                    endif
                  enddo
                  call dreorderi(3*nd,npts,iboxgrad(1,1,1,ithd), &
     &    grad(1,1,istart),iboxsrcind(1,ithd))
                endif
!
!  continue from here
            endif
         enddo
!$OMP END TASKLOOP
        deallocate(iboxlexp)  
      enddo

      !----------------------------------------------------


      !------------- local to local translations ---------

      do ilev = 2,nlevels-1

!$OMP TASKLOOP DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,i,jbox,istart,iend,npts)
         do ibox = laddr(1,ilev),laddr(2,ilev)

            npts = 0

            istart = iexpcse(1,ibox) 
            iend = iexpcse(2,ibox) 
            npts = npts + iend-istart+1

            istart = isrcse(1,ibox)
            iend = isrcse(2,ibox) 
            npts = npts + iend-istart+1

            if(npts.gt.0) then
               do i=1,8
                  jbox = itree(ipointer(5)+8*(ibox-1)+i-1)
                  if(jbox.gt.0) then
                     call l3dlocloc(nd,rscales(ilev), &
     &    centers(1,ibox),rmlexp(iaddr(2,ibox)), &
     &    nterms(ilev),rscales(ilev+1),centers(1,jbox), &
     &    rmlexp(iaddr(2,jbox)),nterms(ilev+1),dc,lca)
                  endif
               enddo
            endif
         enddo
!$OMP END TASKLOOP
      enddo

      !--------------------------------------------------------------------

      call trace%end('FMM3DTree:lfmm3dmain_tree', itimer=itimer, verbose=3)
      
    end subroutine lfmm3dmain_tree


      subroutine eval_local(self)
        !DEC$ ATTRIBUTES ALIAS:"fmm3d_tree_mod_mp_eval_local_" :: eval_local
        class(FMM3DTree), intent(inout) :: self
        !--------------------------------------------
        integer :: ilev,ibox,istart,iend,i,npts
        integer :: nchild
        integer, save :: itimer = 0
        !--------------------------------------------

        call trace%begin('FMM3DTree:eval_local', itimer=itimer, verbose=3)

        do ilev = 0,self%nlevels
          !$OMP TASKLOOP DEFAULT(SHARED) &
          !$OMP PRIVATE(ibox,nchild,istart,iend,i,npts) 
          do ibox = self%laddr(1,ilev),self%laddr(2,ilev)
            nchild=self%itree(self%ipointer(4)+ibox-1)
            if(nchild.eq.0) then 
              istart = self%isrcse(1,ibox) 
              iend = self%isrcse(2,ibox)
              npts = iend-istart+1

              call l3dtaevalg_grad(self%nd,self%scales(ilev),self%centers(1,ibox), &
     &    self%rmlexp(self%iaddr(2,ibox)),self%nterms(ilev),self%sourcesort(1,istart), &
     &    npts,self%gradsort(1,1,istart),self%wlege,self%nlege)
            endif
          enddo
          !$OMP END TASKLOOP
      enddo

      call trace%end('FMM3DTree:eval_local', itimer=itimer, verbose=3)

      end subroutine eval_local


      subroutine eval_direct(self)
        !DEC$ ATTRIBUTES ALIAS:"fmm3d_tree_mod_mp_eval_direct_" :: eval_direct
        class(FMM3DTree), intent(inout) :: self
        !--------------------------------------------
        integer :: ilev,ibox,istarts,iends,npts0,i
        integer :: jbox,jstart,jend,npts
        !--------------------------------------------
        !TODO - Needs testing in parallel!!!
        do ilev=0,self%nlevels
            !!$OMP TASKLOOP DEFAULT(SHARED) &
            !!$OMP PRIVATE(ibox,istarts,iends,npts0,i,jbox,jstart,jend,npts) 
            do ibox = self%laddr(1,ilev),self%laddr(2,ilev)
              istarts = self%isrcse(1,ibox)
              iends = self%isrcse(2,ibox)
              npts0 = iends-istarts+1
              
              do i=1,self%nlist1(ibox)
                jbox = self%list1(i,ibox)
                jstart = self%isrcse(1,jbox)
                jend = self%isrcse(2,jbox)
                npts = jend-jstart+1
                     call l3ddirectdg_grad(self%nd,self%sourcesort(1,jstart), &
     &    self%dipvecsort(1,1,jstart),npts,self%sourcesort(1,istarts), &
     &    npts0,self%gradsort(1,1,istarts),self%thresh, self%isrc, istarts, jstart)     
              enddo
            enddo
           ! !$OMP END TASKLOOP
      enddo
      if ( trace%enabled .and. trace%verbose .ge. 2 ) print *, " finished evaluating direct interactions"
      end subroutine eval_direct

end module fmm3d_tree_mod



!------------_ IMPORTANT NOTE ----------------
! the below l3ddirectdg_grad_vec function calls a c++ version to do direct evaluation of dipole gradients
! but, direct evaluation should be done with MagTense CUDA
! therefore it is commented out so the fast FMM kernels does not need to be compiled
! This should make linking easier. 
!------------------------------------------------------------------------------------------


! !***********************************************************************
!       subroutine l3ddirectdg_grad_vec(nd,sources, &
!                  dipvec,ns,ztarg,nt,grad,thresh)
! !**********************************************************************
! !
! !     This subroutine evaluates the potential and gradient due to a 
! !     collection of sources and adds to existing quantities.
! !   
! !     grad(x) = grad(x) + Gradient( sum  
! !                                    j
! !
! !                            \nabla 1|/|x-x_{j}| \cdot v_{j}
! !                            )
! !                                   
! !      where v_{j} is the dipole orientation vector, 
! !      \nabla denotes the gradient is with respect to the x_{j} 
! !      variable, and Gradient denotes the gradient with respect to
! !      the x variable
! !      If |r| < thresh 
! !          then the subroutine does not update the potential
! !          (recommended value = |boxsize(0)|*machine precision
! !           for boxsize(0) is the size of the computational domain) 
! !
! !
! !-----------------------------------------------------------------------
! !     INPUT:
! !
! !     nd     :    number of charge and dipole densities
! !     sources:    source locations
! !     dipvec :    dipole orientation vector
! !     ns     :    number of sources
! !     ztarg  :    target locations
! !     ntarg  :    number of targets
! !     thresh :    threshold for updating potential,
! !                 potential at target won't be updated if
! !                 |t - s| <= thresh, where t is the target
! !                 location and, and s is the source location 
! !                 
! !-----------------------------------------------------------------------
! !     OUTPUT:
! !
! !     pot    :    updated potential at ztarg 
! !     grad   :    updated gradient at ztarg 
! !
! !-----------------------------------------------------------------------
!       implicit none
! !f2py intent(in) nd,sources,dipvec,ns,ztarg,nt,thresh
! !f2py intent(out) pot,grad
! !
! !c      calling sequence variables
! !  
!       integer ns,nt,nd
!       real *8 sources(3,ns),ztarg(3,nt),dipvec(nd,3,ns)
!       real *8 grad(nd,3,nt)
!       real *8 thresh
!       call l3ddirectdg_cpp_grad(nd,sources, &
!                  dipvec,ns,ztarg,nt,grad,thresh)

!       return
!       end



!**********************************************************************
      subroutine l3dtaevalg_grad(nd,rscale,center,mpole,nterms, &
     		ztarg,ntarg,grad,wlege,nlege)
!**********************************************************************
!
!
!     this subroutine evaluates the potentials and gradients due to  
!     an incoming local expansion and increments inputs accordingly:
!
!     grad =  grad + 
!              Gradient( sum sum mpole(n,m)r^{n}Y_nm(theta,phi)/sqrt(2n+1))
!                         n   m
!-----------------------------------------------------------------------
!     INPUT:
!
!     nd     :    number of multipole expansions
!     rscale :    scaling parameter 
!     center :    expansion center
!     mpole  :    local expansion 
!     nterms :    order of the multipole expansion
!     ztarg  :    target location
!     ntarg  :    number of target locations
!     wlege  :    precomputed array of scaling coeffs for Pnm
!     nlege  :    dimension parameter for wlege
!-----------------------------------------------------------------------
!     OUTPUT:
!
!     pot    :   updated potentials at targets
!     grad   :   updated gradients at targets 
!----------------------------------------------------------------------
      implicit none
!
!c     calling sequence variables
!
      integer nterms,nlege,ntarg,nd
      real *8 rscale,center(3),ztarg(3,ntarg)
      real *8 grad(nd,3,ntarg)
      complex *16 mpole(nd,0:nterms,-nterms:nterms)
      real *8 wlege(0:nlege,0:nlege)
!
!c     temporary variables
!
      integer idim
      real *8, allocatable :: ynm(:,:),ynmd(:,:),fr(:),frder(:)
      complex *16, allocatable :: ephi(:)
      integer i,j,k,l,m,n,itarg
      real *8 done,r,theta,phi,zdiff(3)
      real *8 ctheta,stheta,cphi,sphi
      real *8 d,rx,ry,rz,thetax,thetay,thetaz,phix,phiy,phiz,rs
      real *8 rtmp1,rtmp2,rtmp3,rtmp4,rtmp5,rtmp6
      complex *16 ephi1
      real *8 ur(nd),utheta(nd),uphi(nd)
!
      complex *16 eye
      complex *16 ztmp1,ztmp2,ztmp3,ztmpsum,z
      real *8 rscaleinv
!
      data eye/(0.0d0,1.0d0)/
!
      done=1.0d0
      allocate(ephi(0:nterms+1))
      allocate(fr(0:nterms+1),frder(0:nterms))
      allocate(ynm(0:nterms,0:nterms))
      allocate(ynmd(0:nterms,0:nterms))
!
      do itarg=1,ntarg
        zdiff(1)=ztarg(1,itarg)-center(1)
        zdiff(2)=ztarg(2,itarg)-center(2)
        zdiff(3)=ztarg(3,itarg)-center(3)
!
        call cart2polar(zdiff,r,theta,phi)

        ctheta = dcos(theta)
        stheta = dsin(theta)
        cphi = dcos(phi)
        sphi = dsin(phi)
        ephi1 = dcmplx(cphi,sphi)
!
!     compute exp(eye*m*phi) array
!
        ephi(0)=done
        ephi(1)=ephi1
        d = r/rscale
        fr(0) = 1.0d0
        fr(1) = fr(0)*d
        do i=2,nterms+1
          fr(i) = fr(i-1)*d
          ephi(i)=ephi(i-1)*ephi1
        enddo
        frder(0) = 0
        do i=1,nterms
          frder(i) = i*fr(i-1)/rscale
        enddo
!
!    get the associated Legendre functions:
!
        call ylgndr2sfw(nterms,ctheta,ynm,ynmd,wlege,nlege)
        do l = 0,nterms
          rs = sqrt(1.0d0/(2*l+1))
          do m=0,l
            ynm(l,m) = ynm(l,m)*rs
            ynmd(l,m) = ynmd(l,m)*rs
          enddo
        enddo
!
!     compute coefficients in change of variables from spherical
!     to Cartesian gradients. In phix, phiy, we leave out the 
!     1/sin(theta) contribution, since we use values of Ynm (which
!     multiplies phix and phiy) that are scaled by 
!     1/sin(theta).
!
!
!     NOTE: sphereical derivative needs to be fixed for r=0
!

        rscaleinv = 1.0d0/rscale
        rx = stheta*cphi
        thetax = ctheta*cphi*rscaleinv
        phix = -sphi*rscaleinv
        ry = stheta*sphi
        thetay = ctheta*sphi*rscaleinv
        phiy = cphi*rscaleinv
        rz = ctheta
        thetaz = -stheta*rscaleinv
        phiz = 0.0d0
!
        do idim=1,nd
          ur(idim) = real(mpole(idim,0,0))*frder(0)
          utheta(idim) = 0.0d0
          uphi(idim) = 0.0d0
        enddo
!
        do n=1,nterms
          rtmp1 = fr(n)*ynm(n,0)
          rtmp2 = frder(n)*ynm(n,0)
          rtmp3 = -fr(n-1)*ynmd(n,0)*stheta
          do idim=1,nd
            ur(idim)=ur(idim)+real(mpole(idim,n,0))*rtmp2
            utheta(idim)=utheta(idim)+real(mpole(idim,n,0))*rtmp3
          enddo
!
	      do m=1,n
            rtmp1 = fr(n)*ynm(n,m)*stheta
            rtmp4 = frder(n)*ynm(n,m)*stheta
            rtmp5 = -fr(n-1)*ynmd(n,m)
            rtmp6 = -m*fr(n-1)*ynm(n,m)

            do idim=1,nd
              rtmp2 = 2*real(mpole(idim,n,m)*ephi(m)) 

              ur(idim) = ur(idim) + rtmp4*rtmp2
              utheta(idim) = utheta(idim)+rtmp5*rtmp2
              rtmp2 = 2*dimag(mpole(idim,n,m)*ephi(m))
              uphi(idim) = uphi(idim) + rtmp6*rtmp2
            enddo
          enddo
        enddo
!
        do idim=1,nd
          grad(idim,1,itarg)=grad(idim,1,itarg)+ur(idim)*rx+ &
               utheta(idim)*thetax+uphi(idim)*phix
          grad(idim,2,itarg)=grad(idim,2,itarg)+ur(idim)*ry+ &
               utheta(idim)*thetay+uphi(idim)*phiy
          grad(idim,3,itarg)=grad(idim,3,itarg)+ur(idim)*rz+ &
               utheta(idim)*thetaz+uphi(idim)*phiz
        enddo
 1000 continue
      enddo
      return
      end




!***********************************************************************
      subroutine l3ddirectdg_grad(nd,sources, &
                 dipvec,ns,ztarg,nt,grad,thresh, isrc, istart, jstart)
!**********************************************************************
!
!     This subroutine evaluates the potential and gradient due to a 
!     collection of sources and adds to existing quantities.
!
!   
!     grad(x) = grad(x) + Gradient( sum  
!                                    j
!
!                            \nabla 1|/|x-x_{j}| \cdot v_{j}
!                            )
!                                   
!      where v_{j} is the dipole orientation vector, 
!      \nabla denotes the gradient is with respect to the x_{j} 
!      variable, and Gradient denotes the gradient with respect to
!      the x variable
!      If |r| < thresh 
!          then the subroutine does not update the potential
!          (recommended value = |boxsize(0)|*machine precision
!           for boxsize(0) is the size of the computational domain) 
!
!
!-----------------------------------------------------------------------
!     INPUT:
!
!     nd     :    number of charge and dipole densities
!     sources:    source locations
!     dipvec :    dipole orientation vector
!     ns     :    number of sources
!     ztarg  :    target locations
!     ntarg  :    number of targets
!     thresh :    threshold for updating potential,
!                 potential at target won't be updated if
!                 |t - s| <= thresh, where t is the target
!                 location and, and s is the source location 
!                 
!-----------------------------------------------------------------------
!     OUTPUT:
!
!     pot    :    updated potential at ztarg 
!     grad   :    updated gradient at ztarg 
!
!-----------------------------------------------------------------------
      implicit none
!cf2py intent(in) nd,sources,dipvec,ns,ztarg,nt,thresh
!cf2py intent(out) pot,grad
      intent(in) nd,sources,dipvec,ns,ztarg,nt,thresh
      integer, intent(in) :: isrc(nt), istart, jstart
      
      intent(out) grad
!c
!cc      calling sequence variables
!c  
      integer ns,nt,nd
      real *8 sources(3,ns),ztarg(3,nt),dipvec(nd,3,ns)
      real *8 grad(nd,3,nt)
      real *8 thresh
      
!c
!cc     temporary variables
!c
      real *8 zdiff(3),dd,d,dinv,dinv2,dotprod
      real *8 cd,cd2,cd3,cd4
      real *8 threshsq
      integer i,j,idim

      threshsq = thresh**2
      do i=1,nt

    !  print *, " ====== iteration ", i, " fmm target ", isrc(i+istart-1), " ====== "
      !!$OMP TASKLOOP default(none) &
      !!$omp private(j,zdiff, dd, dinv,dinv2,dotprod,cd,cd2,cd3,cd4,idim) &
      !!$omp shared(i, nd,ns,sources,ztarg,dipvec,grad,threshsq)
        do j=1,ns
          zdiff(1) = ztarg(1,i)-sources(1,j)
          zdiff(2) = ztarg(2,i)-sources(2,j)
          zdiff(3) = ztarg(3,i)-sources(3,j)

          dd = zdiff(1)**2 + zdiff(2)**2 + zdiff(3)**2
         ! print *, " nbor ", j, " real index", isrc(j+jstart-1), "zdiff = ", zdiff 
          if(dd.lt.threshsq) goto 1000

          dinv2 = 1/dd
          dinv = sqrt(dinv2)
          cd = dinv
          cd2 = -cd*dinv2
          cd3 = -3*cd*dinv2*dinv2




          do idim=1,nd
          
            dotprod = zdiff(1)*dipvec(idim,1,j)+ &
                    zdiff(2)*dipvec(idim,2,j)+ &
                    zdiff(3)*dipvec(idim,3,j)
            cd4 = cd3*dotprod

            grad(idim,1,i) = grad(idim,1,i) + (cd4*zdiff(1) - &
              cd2*dipvec(idim,1,j))
            grad(idim,2,i) = grad(idim,2,i) + (cd4*zdiff(2) - &
              cd2*dipvec(idim,2,j))
            grad(idim,3,i) = grad(idim,3,i) + (cd4*zdiff(3) - &
              cd2*dipvec(idim,3,j))
          enddo
 1000     continue
        enddo
      enddo


      return
      end



subroutine l3dformmpd_new(nd, rscale, sources, dipvec, ns, center, &
                                nterms, mpole, wlege, nlege)
    implicit none

    ! --- Arguments ---
    integer, intent(in) :: nd, ns, nterms, nlege
    real(kind(0.d0)), intent(in) :: rscale, center(3)
    real(kind(0.d0)), intent(in) :: sources(3, ns)
    real(kind(0.d0)), intent(in) :: dipvec(nd, 3, ns)
    real(kind(0.d0)), intent(in) :: wlege(0:nlege, 0:nlege)
    ! mpole is modified (added to)
    complex(kind(0.d0)), intent(inout) :: mpole(nd, 0:nterms, -nterms:nterms)

    ! --- Local Parameters ---
    complex(kind(0.d0)), parameter :: eye = (0.0d0, 1.0d0)
    real(kind(0.d0)), parameter    :: one = 1.0d0

    ! --- Local Variables (Thread-Safe on the Stack) ---
    integer :: i, j, k, l, m, n, isrc, idim
    real(kind(0.d0)) :: zdiff(3), r, theta, phi
    real(kind(0.d0)) :: stheta, ctheta, sphi, cphi
    real(kind(0.d0)) :: rx, ry, rz, thetax, thetay, thetaz, phix, phiy, phiz
    real(kind(0.d0)) :: fruse, d
    complex(kind(0.d0)) :: ephi1, ur, utheta, uphi, ux, uy, uz, zzz

    ! Allocatable scratch arrays (Automatically thread-private)
    real(kind(0.d0)), allocatable :: ynm(:,:), fr(:), rfac(:), frder(:), ynmd(:,:)
    complex(kind(0.d0)), allocatable :: ephi(:)

    ! Allocate local memory
    allocate(ynm(0:nterms, 0:nterms), fr(0:nterms+1))
    allocate(frder(0:nterms), ynmd(0:nterms, 0:nterms))
    allocate(ephi(-nterms-1:nterms+1))
    allocate(rfac(0:nterms))

    ! Initialize scaling factors
    do i = 0, nterms
        rfac(i) = one / sqrt(2.0d0 * i + one)
    end do

    ! Loop over sources
    do isrc = 1, ns
        zdiff(1) = sources(1, isrc) - center(1)
        zdiff(2) = sources(2, isrc) - center(2)
        zdiff(3) = sources(3, isrc) - center(3)

        ! Convert to polar coordinates
        call cart2polar(zdiff, r, theta, phi)
        
        ctheta = cos(theta)
        stheta = sin(theta)
        cphi   = cos(phi)
        sphi   = sin(phi)
        ephi1  = cmplx(cphi, sphi, kind=kind(0.d0))

        ! Compute e^(i*m*phi) and radial scaling powers
        ephi(0)  = one
        ephi(1)  = ephi1
        ephi(-1) = conjg(ephi1)
        
        fr(0) = one
        d = r / rscale
        fr(1) = d
        
        do i = 2, nterms + 1
            fr(i) = fr(i-1) * d
            ephi(i) = ephi(i-1) * ephi1
            ephi(-i) = ephi(-i+1) * ephi(-1)
        end do

        frder(0) = 0.0d0
        do i = 1, nterms
            frder(i) = i * fr(i-1) / rscale
        end do

        ! Spherical to Cartesian gradient transformation components
        rx = stheta * cphi
        thetax = ctheta * cphi
        phix = -sphi
        ry = stheta * sphi
        thetay = ctheta * sphi
        phiy = cphi
        rz = ctheta
        thetaz = -stheta
        phiz = 0.0d0

        ! Legendre functions (Ensure this routine is also thread-safe/recursive)
        call ylgndr2sfw(nterms, ctheta, ynm, ynmd, wlege, nlege)
        
        do i = 0, nterms
            do j = 0, nterms
                ynm(j, i)  = ynm(j, i) * rfac(j)
                ynmd(j, i) = ynmd(j, i) * rfac(j)
            end do
        end do

        ! --- N=0 Contribution ---
        ur = ynm(0, 0) * frder(0)
        ux = ur * rx 
        uy = ur * ry 
        uz = ur * rz
        do idim = 1, nd
            zzz = dipvec(idim, 1, isrc)*ux + dipvec(idim, 2, isrc)*uy + &
                  dipvec(idim, 3, isrc)*uz
            mpole(idim, 0, 0) = mpole(idim, 0, 0) + zzz
        end do

        ! --- N > 0 Contributions ---
        do n = 1, nterms
            fruse = fr(n-1) / rscale
            ur = ynm(n, 0) * frder(n)
            utheta = -fruse * ynmd(n, 0) * stheta
            ux = ur*rx + utheta*thetax 
            uy = ur*ry + utheta*thetay 
            uz = ur*rz + utheta*thetaz
            
            do idim = 1, nd
                zzz = dipvec(idim, 1, isrc)*ux + dipvec(idim, 2, isrc)*uy + &
                      dipvec(idim, 3, isrc)*uz
                mpole(idim, n, 0) = mpole(idim, n, 0) + zzz
            end do

            do m = 1, n
                ! m positive terms (using conjugate ephi for m > 0 logic)
                ur = frder(n) * ynm(n, m) * stheta * ephi(-m)
                utheta = -ephi(-m) * fruse * ynmd(n, m)
                uphi = -eye * m * ephi(-m) * fruse * ynm(n, m)
                ux = ur*rx + utheta*thetax + uphi*phix
                uy = ur*ry + utheta*thetay + uphi*phiy
                uz = ur*rz + utheta*thetaz + uphi*phiz
                do idim = 1, nd
                    zzz = dipvec(idim, 1, isrc)*ux + dipvec(idim, 2, isrc)*uy + &
                          dipvec(idim, 3, isrc)*uz
                    mpole(idim, n, m) = mpole(idim, n, m) + zzz
                end do

                ! m negative terms
                ur = frder(n) * ynm(n, m) * stheta * ephi(m)
                utheta = -ephi(m) * fruse * ynmd(n, m)
                uphi = eye * m * ephi(m) * fruse * ynm(n, m)
                ux = ur*rx + utheta*thetax + uphi*phix
                uy = ur*ry + utheta*thetay + uphi*phiy
                uz = ur*rz + utheta*thetaz + uphi*phiz
                do idim = 1, nd
                    zzz = dipvec(idim, 1, isrc)*ux + dipvec(idim, 2, isrc)*uy + &
                          dipvec(idim, 3, isrc)*uz
                    mpole(idim, n, -m) = mpole(idim, n, -m) + zzz
                end do
            end do
        end do
    end do
    
    ! Clean up memory to avoid leaks in the TASKLOOP
    deallocate(ynm, fr, frder, ynmd, ephi, rfac)

end subroutine l3dformmpd_new