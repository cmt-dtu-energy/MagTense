module fmm3d_tree_mod
        implicit none
      type :: FMM3DTree
            integer :: nlmax = 51
            integer :: nlevels = 0
            integer :: nboxes = 0
            integer(8) :: ltree = 0
            integer :: nlmin = 0
            integer :: iper = 0
            integer :: ifunif = 0

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
            double precision, contiguous, pointer ::  source(:,:),targ(:,:)
            double precision, contiguous, pointer ::  charge(:,:)
            double precision, contiguous, pointer ::  dipvec(:,:,:)
            double precision, contiguous, pointer ::  pot(:,:),grad(:,:,:),hess(:,:,:)
            double precision, contiguous, pointer ::  pottarg(:,:),gradtarg(:,:,:),hesstarg(:,:,:)
            double precision, contiguous, pointer :: sourcesort(:,:),targsort(:,:)
            double precision, contiguous, pointer :: chargesort(:,:)
            double precision, contiguous, pointer :: dipvecsort(:,:,:)
            double precision, contiguous, pointer :: potsort(:,:),gradsort(:,:,:)
            double precision, contiguous, pointer :: hesssort(:,:,:)
            double precision, contiguous, pointer :: pottargsort(:,:)
            double precision, contiguous, pointer :: gradtargsort(:,:,:)
            double precision, contiguous, pointer :: hesstargsort(:,:,:)
            !-----------------------------------------
            integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
            integer nd,ier
            double precision eps
            !-----------------------------------------
            integer ndiv, idivflag
            integer(8), dimension(8) :: ipointer
            double precision expc(3)
            integer :: iexpc
            integer, contiguous, pointer :: itree(:)
            double precision, contiguous, pointer :: boxsize(:), treecenters(:,:), centers(:,:)
            integer, contiguous, pointer :: isrcse(:,:),itargse(:,:),iexpcse(:,:)
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
            double complex texpssort(100)
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
        double precision, contiguous, pointer :: iboxpot(:,:,:)
        double precision, contiguous, pointer :: iboxgrad(:,:,:,:)
        double precision, contiguous, pointer :: iboxhess(:,:,:,:)
        double precision, contiguous, pointer :: iboxsrc(:,:,:)
        integer, contiguous, pointer :: iboxsrcind(:,:)
        integer, contiguous, pointer :: iboxfl(:,:,:)
        !-------- 




        contains
          procedure :: full_fmm
          procedure :: build1
          procedure :: reset_sort_arg
          procedure :: dealloc
          procedure :: build2
          procedure :: lfmm3dmain_tree
          procedure :: reset_expansion_coeff
      end type FMM3DTree
    contains 


    subroutine full_fmm(self,nd,eps,nsource,source, &
            dipvec,pot,grad,ier)
        class(FMM3DTree), intent(inout) :: self
        !------------------------------------------------
        double precision eps

        integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
        integer nd,iper,ier
        
        double precision, target :: source(3,nsource),targ(3,1)
        double precision, target :: charge(nd,1)
        
        double precision, target :: dipvec(nd,3,nsource)

        double precision, target :: pot(nd,nsource),grad(nd,3,nsource)
        double precision, target :: pottarg(nd,1),gradtarg(nd,3,1)
        double precision, target :: hess(nd,6,1),hesstarg(nd,6,1)


        integer, contiguous, pointer ::  laddr(:,:)

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 0

      ntarg = 0

      ier = 0

        !---------- setting self varialbes ----------------
        self%source => source
        self%nsource = nsource
        self%dipvec => dipvec

        self%pot => pot
        self%grad => grad
        self%eps = eps
        self%nd = nd
        self%ier = ier
        self%ifcharge = ifcharge
        self%ifdipole = ifdipole
        self%ntarg = ntarg
        self%targ => targ
        self%ifpgh = ifpgh
        self%ifpghtarg = ifpghtarg
        self%pottarg => pottarg
        self%gradtarg => gradtarg
        self%hess => hess
        self%hesstarg => hesstarg
        !--------------------------------------------------


      print *, " calling fmm tree "




        !self%nthd = omp_get_num_threads()


        call self%build1()
        self%laddr(1:2,0:self%nlevels) => self%itree(self%ipointer(1) : self%ipointer(1)+(self%nlevels + 1)*2-1)
        call self%build2()

        
        
        !print *, " laddr ", laddr
        !allocate (laddr(2,0:self%nlevels))
        !laddr = self%itree(self%ipointer(1):self%ipointer(1)+(self%nlevels + 1)*2-1)    
        !print *, " laddr ", laddr


    !   call lfmm3d_tree(nd,eps,nsource,source,ifcharge,charge, &
    !        ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,targ, &
    !        ifpghtarg,pottarg,gradtarg,hesstarg,ier)




        call self%lfmm3dmain_tree(self%nd,self%eps, &
            &   self%nsource,self%sourcesort, &
            &   self%ifcharge,self%chargesort, &
            &   self%ifdipole,self%dipvecsort, &
            &   self%ntarg,self%targsort,self%nexpc,self%expcsort, &
            &   self%iaddr,self%rmlexp,self%lmptot,self%mptemp,self%mptemp2,self%lmptemp, &
            &   self%itree,self%ltree,self%ipointer,self%ndiv,self%nlevels, & 
            &   self%nboxes,self%iper,self%boxsize,self%treecenters,self%isrcse,self%itargse,self%iexpcse, &
            &   self%scales,self%laddr,self%nterms, &
            &   self%ifpgh,self%potsort,self%gradsort,self%hesssort, &
            &   self%ifpghtarg,self%pottargsort,self%gradtargsort,self%hesstargsort,self%ntj, &
            &   self%texpssort,self%scjsort,self%ifnear,self%timeinfo,self%ier)

        print *, "self%ipointer(1) ", self%ipointer

        ! call lfmm3dmain_tree_old(self%nd,self%eps, &
        !     &   self%nsource,self%sourcesort, &
        !     &   self%ifcharge,self%chargesort, &
        !     &   self%ifdipole,self%dipvecsort, &
        !     &   self%ntarg,self%targsort,self%nexpc,self%expcsort, &
        !     &   self%iaddr,self%rmlexp,self%lmptot,self%mptemp,self%mptemp2,self%lmptemp, &
        !     &   self%itree,self%ltree,self%ipointer,self%ndiv,self%nlevels, &
        !     &   self%nboxes,self%iper,self%boxsize,self%treecenters,self%isrcse,self%itargse,self%iexpcse, &
        !     &   self%scales,self%itree(self%ipointer(1)),self%nterms, &
        !     &   self%ifpgh,self%potsort,self%gradsort,self%hesssort, &
        !     &   self%ifpghtarg,self%pottargsort,self%gradtargsort,self%hesstargsort,self%ntj, &
        !     &   self%texpssort,self%scjsort,self%ifnear,self%timeinfo,self%ier)

        if(self%ifpgh.ge.1) then
            call dreorderi(self%nd,self%nsource,self%potsort,self%pot,self%isrc)
        endif
        if(self%ifpgh.ge.2) then 
            call dreorderi(3*self%nd,self%nsource,self%gradsort,self%grad,self%isrc)
            call drescale(self%nd*3*self%nsource,self%grad,self%b0inv)
        endif

        if(self%ifpgh.ge.3) then 
            call dreorderi(6*self%nd,self%nsource,self%hesssort,self%hess,self%isrc)
            call drescale(self%nd*6*self%nsource,self%hess,self%b0inv2)
        endif


        if(self%ifpghtarg.ge.1) then
            call dreorderi(self%nd,self%ntarg,self%pottargsort,self%pottarg,self%itarg)
        endif

        if(self%ifpghtarg.ge.2) then
            call dreorderi(3*self%nd,self%ntarg,self%gradtargsort,self%gradtarg,self%itarg)
            call drescale(self%nd*3*self%ntarg,self%gradtarg,self%b0inv)
        endif

        if(self%ifpghtarg.ge.3) then
            call dreorderi(6*self%nd,self%ntarg,self%hesstargsort,self%hesstarg,self%itarg)
            call drescale(self%nd*6*self%ntarg,self%hesstarg,self%b0inv2)
        endif


            call self%dealloc()
    end subroutine full_fmm


    subroutine build1(self)
        class(FMM3DTree), intent(inout) :: self
        !------------------------------------------------
        !----------------
        integer, contiguous, pointer :: itree(:)
        double precision, contiguous, pointer :: boxsize(:), treecenters(:,:)
        integer, contiguous, pointer :: isrcse(:,:),itargse(:,:),iexpcse(:,:)
        integer, contiguous, pointer :: isrc(:),itarg(:)
        !--------------------------------
        double precision, contiguous, pointer :: sourcesort(:,:),targsort(:,:)
        double precision, contiguous, pointer :: chargesort(:,:) 
        double precision, contiguous, pointer :: dipvecsort(:,:,:)
        double precision, contiguous, pointer :: potsort(:,:),gradsort(:,:,:)
        double precision, contiguous, pointer :: hesssort(:,:,:)
        double precision, contiguous, pointer :: pottargsort(:,:)
        double precision, contiguous, pointer :: gradtargsort(:,:,:)
        double precision, contiguous, pointer :: hesstargsort(:,:,:)
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
        print *, " starting tree build "
    
        call lndiv(self%eps,self%nsource,self%ntarg,self%ifcharge,self%ifdipole,self%ifpgh, &
                &    self%ifpghtarg,self%ndiv,self%idivflag) 


        call pts_tree_mem(self%source,self%nsource,self%targ,self%ntarg,self%idivflag,self%ndiv,self%nlmin, &
            &    self%nlmax,self%iper,self%ifunif,self%nlevels,self%nboxes,self%ltree)



        allocate(itree(self%ltree))
        allocate(boxsize(0:self%nlevels))
        allocate(treecenters(3,self%nboxes))

        call pts_tree_build(self%source,self%nsource,self%targ,self%ntarg,self%idivflag,self%ndiv, &
        &    self%nlmin,self%nlmax,self%iper,self%ifunif,self%nlevels,self%nboxes,self%ltree,itree,self%ipointer, &
        &    treecenters,boxsize)
    

        allocate(isrcse(2,self%nboxes),itargse(2,self%nboxes),iexpcse(2,self%nboxes))
        allocate(isrc(self%nsource),itarg(self%ntarg))

        call pts_tree_sort(self%nsource,self%source,itree,self%ltree,self%nboxes,self%nlevels, &
        &    self%ipointer,treecenters,isrc,isrcse)
        
        call pts_tree_sort(self%ntarg,self%targ,itree,self%ltree,self%nboxes,self%nlevels, &
        &    self%ipointer,treecenters,itarg,itargse)
        
        call pts_tree_sort(self%nexpc,self%expc,itree,self%ltree,self%nboxes,self%nlevels, &
        &    self%ipointer,treecenters,self%iexpc,iexpcse)

        self%itree => itree
        self%boxsize => boxsize
        self%treecenters => treecenters
        self%centers => treecenters ! same 
        self%isrcse => isrcse
        self%itargse => itargse
        self%iexpcse => iexpcse
        self%isrc => isrc
        self%itarg => itarg

        !-------- end of tree build --------------------

        print *, " finished tree build "
        !------ set scaling 

        self%b0 = self%boxsize(0)
        self%b0inv = 1.0d0/self%b0
        self%b0inv2 = self%b0inv**2
        self%b0inv3 = self%b0inv2*self%b0inv





        !-------------- allocate sorted source and targ arrays------

        print *, " allocating sorted arrays "

        allocate(sourcesort(3,self%nsource))
        allocate(targsort(3,self%ntarg))

        print *, " source and targ"
        if(self%ifcharge.eq.1) allocate(chargesort(self%nd,self%nsource))

        if(self%ifdipole.eq.1) then
            allocate(dipvecsort(self%nd,3,self%nsource))
        endif

        if(self%ifpgh.eq.1) then 
            allocate(potsort(self%nd,self%nsource),gradsort(self%nd,3,1),hesssort(self%nd,6,1))
        else if(self%ifpgh.eq.2) then
            allocate(potsort(self%nd,self%nsource),gradsort(self%nd,3,self%nsource), &
        &    hesssort(self%nd,6,1))
        else if(self%ifpgh.eq.3) then
            allocate(potsort(self%nd,self%nsource),gradsort(self%nd,3,self%nsource), &
        &    hesssort(self%nd,6,self%nsource))
        else
            allocate(potsort(self%nd,1),gradsort(self%nd,3,1),hesssort(self%nd,6,1))
        endif
        print *, " halfway"

        if(self%ifpghtarg.eq.1) then
            allocate(pottargsort(self%nd,self%ntarg),gradtargsort(self%nd,3,1), &
        &    hesstargsort(self%nd,6,1))
        else if(self%ifpghtarg.eq.2) then
            allocate(pottargsort(self%nd,self%ntarg),gradtargsort(self%nd,3,self%ntarg), &
        &    hesstargsort(self%nd,6,1))
        else if(self%ifpghtarg.eq.3) then
            allocate(pottargsort(self%nd,self%ntarg),gradtargsort(self%nd,3,self%ntarg), &
        &    hesstargsort(self%nd,6,self%ntarg))
        else
            allocate(pottargsort(self%nd,1),gradtargsort(self%nd,3,1), &
        &    hesstargsort(self%nd,6,1))
        endif
            !------------------------------------------------------


        self%sourcesort => sourcesort
        self%targsort => targsort
        self%chargesort => chargesort
        self%dipvecsort => dipvecsort
        self%potsort => potsort
        self%gradsort => gradsort
        self%hesssort => hesssort
        self%pottargsort => pottargsort
        self%gradtargsort => gradtargsort
        self%hesstargsort => hesstargsort
        !--------------------------------------------------------
        call self%reset_sort_arg()

        print *, " finished allocating sorted arrays "


        allocate(nterms(0:self%nlevels))

        self%nmax = 0
        do i=0,self%nlevels
            call l3dterms(self%eps,nterms(i))
            if(nterms(i).gt.self%nmax) self%nmax = nterms(i)
        enddo
        self%nterms => nterms
!       



        call dreorderf(3,self%nsource,self%source,self%sourcesort,self%isrc)

      call drescale(3*self%nsource,self%sourcesort,self%b0inv)


        if(self%ifcharge.eq.1) then
            call dreorderf(self%nd,self%nsource,self%charge,self%chargesort, &
        &    self%isrc)
            call drescale(self%nd*self%nsource,self%chargesort,self%b0inv)
        endif


        if(self%ifdipole.eq.1) then
            call dreorderf(3*self%nd,self%nsource,self%dipvec,self%dipvecsort, &
        &    self%isrc)
            call drescale(3*self%nd*self%nsource,self%dipvecsort,self%b0inv2)
        endif

        call dreorderf(3,self%ntarg,self%targ,self%targsort,self%itarg)
        call drescale(3*self%ntarg,self%targsort,self%b0inv)



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

        print *, " done with build1 "

    end subroutine build1


    subroutine build2(self)
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
      print *, " FMM3DTree: max number of terms = ", self%nmax


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

      print *, " done with build2 "


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
        !$OMP PARALLEL DO DEFAULT(SHARED) &
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
        !$OMP END PARALLEL DO
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

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,istart,iend,npts) &
      !$OMP REDUCTION(max:nmaxt)
            do ibox=1,self%nboxes
              if(self%list4ct(ibox).gt.0) then
                istart = self%isrcse(1,ibox)
                iend = self%isrcse(2,ibox)
                npts = iend-istart+1
                if(npts.gt.nmaxt) nmaxt = npts
              endif
            enddo
      !$OMP END PARALLEL DO

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
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,istart,iend,npts) &
    !$OMP REDUCTION(max:nmaxt)
          do ibox=1,self%nboxes
            if(self%nlist3(ibox).gt.0) then
              istart = self%isrcse(1,ibox)
              iend = self%isrcse(2,ibox)
              npts = iend-istart+1
              if(npts.gt.nmaxt) nmaxt = npts

              istart = self%itargse(1,ibox)
              iend = self%itargse(2,ibox)
              npts = iend - istart + 1
              if(npts.gt.nmaxt) nmaxt = npts
            endif
          enddo
    !$OMP END PARALLEL DO

      allocate(self%iboxsrcind(nmaxt,self%nthd))
      allocate(self%iboxsrc(3,nmaxt,self%nthd))
      allocate(self%iboxpot(self%nd,nmaxt,self%nthd))
      allocate(self%iboxgrad(self%nd,3,nmaxt,self%nthd))
      allocate(self%iboxhess(self%nd,6,nmaxt,self%nthd))




    end subroutine build2



    subroutine reset_expansion_coeff(self)
        class(FMM3DTree), intent(inout) :: self
        !------------------------------------------------
        integer :: ilev, ibox
        integer :: i,j,k,idim
        !------------------------------------------------


        !------- reset texpssort to zero --------------
        ! NOTO - tsort in lfmm3dmain_tree
        self%texpssort = 0.0d0
        !-------------------------------------------


      !-------- reset rmlexp to zero --------------
      do ilev = 0,self%nlevels
        !$OMP PARALLEL DO DEFAULT(SHARED) &
        !$OMP PRIVATE(ibox)
        do ibox=self%laddr(1,ilev),self%laddr(2,ilev)
          call mpzero(self%nd,self%rmlexp(self%iaddr(1,ibox)),self%nterms(ilev))
          call mpzero(self%nd,self%rmlexp(self%iaddr(2,ibox)),self%nterms(ilev))
        enddo
        !$OMP END PARALLEL DO
      enddo
      !-------------------------------------------


      !$OMP PARALLEL DO DEFAULT(SHARED) &
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
      !$OMP END PARALLEL DO

            

    end subroutine reset_expansion_coeff




    subroutine reset_sort_arg(self)
        class(FMM3DTree), intent(inout) :: self
        !------------------------------------------------
        integer i,idim
        !------------------------------------------------
        double precision, pointer :: sourcesort(:,:),targsort(:,:)
        double precision, pointer :: chargesort(:,:)
        double precision, pointer :: dipvecsort(:,:,:)
        double precision, pointer :: potsort(:,:),gradsort(:,:,:)
        double precision, pointer :: hesssort(:,:,:)
        double precision, pointer :: pottargsort(:,:)
        double precision, pointer :: gradtargsort(:,:,:)
        double precision, pointer :: hesstargsort(:,:,:)
        integer ifpgh,ifpghtarg, nsource,ntarg,nd
        !------------------------------------------------

        sourcesort => self%sourcesort
        targsort => self%targsort
        chargesort => self%chargesort
        dipvecsort => self%dipvecsort
        potsort => self%potsort
        gradsort => self%gradsort
        hesssort => self%hesssort
        pottargsort => self%pottargsort
        gradtargsort => self%gradtargsort

        ifpgh = self%ifpgh
        ifpghtarg = self%ifpghtarg
        nsource = self%nsource
        ntarg = self%ntarg
        nd = self%nd

        if(ifpgh.eq.1) then
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)
            do i=1,nsource
            do idim=1,nd
                potsort(idim,i) = 0
            enddo
            enddo
    !$OMP END PARALLEL DO
        endif

        if(ifpgh.eq.2) then
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)

            do i=1,nsource
            do idim=1,nd
                potsort(idim,i) = 0
                gradsort(idim,1,i) = 0
                gradsort(idim,2,i) = 0
                gradsort(idim,3,i) = 0
            enddo
            enddo
    !$OMP END PARALLEL DO
        endif


        if(ifpgh.eq.3) then
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)
            do i=1,nsource
            do idim=1,nd
                potsort(idim,i) = 0
                gradsort(idim,1,i) = 0
                gradsort(idim,2,i) = 0
                gradsort(idim,3,i) = 0
                hesssort(idim,1,i) = 0
                hesssort(idim,2,i) = 0
                hesssort(idim,3,i) = 0
                hesssort(idim,4,i) = 0
                hesssort(idim,5,i) = 0
                hesssort(idim,6,i) = 0
            enddo
            enddo
    !$OMP END PARALLEL DO
        endif



    !
    !c       initialize potential and gradient  at targ
    !        locations
    !
        if(ifpghtarg.eq.1) then
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)
            do i=1,ntarg
            do idim=1,nd
                pottargsort(idim,i) = 0
            enddo
            enddo
    !$OMP END PARALLEL DO
        endif

        if(ifpghtarg.eq.2) then
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)
            do i=1,ntarg
            do idim=1,nd
                pottargsort(idim,i) = 0
                gradtargsort(idim,1,i) = 0
                gradtargsort(idim,2,i) = 0
                gradtargsort(idim,3,i) = 0
            enddo
            enddo
    !$OMP END PARALLEL DO
        endif

        if(ifpghtarg.eq.3) then
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)
            do i=1,ntarg
            do idim=1,nd
                pottargsort(idim,i) = 0
                gradtargsort(idim,1,i) = 0
                gradtargsort(idim,2,i) = 0
                gradtargsort(idim,3,i) = 0
                hesstargsort(idim,1,i) = 0
                hesstargsort(idim,2,i) = 0
                hesstargsort(idim,3,i) = 0
                hesstargsort(idim,4,i) = 0
                hesstargsort(idim,5,i) = 0
                hesstargsort(idim,6,i) = 0
            enddo
            enddo
    !$OMP END PARALLEL DO
        endif

    end subroutine reset_sort_arg


    subroutine dealloc(self)
        class(FMM3DTree), intent(inout) :: self
        !------------------------------------------------


        if (associated(self%itree)) deallocate(self%itree)
        if (associated(self%boxsize)) deallocate(self%boxsize)
        if (associated(self%treecenters)) deallocate(self%treecenters)
        if (associated(self%isrcse)) deallocate(self%isrcse)
        if (associated(self%itargse)) deallocate(self%itargse)
        if (associated(self%iexpcse)) deallocate(self%iexpcse)
        if (associated(self%isrc)) deallocate(self%isrc)
        if (associated(self%itarg)) deallocate(self%itarg)
        if (associated(self%sourcesort)) deallocate(self%sourcesort)
        if (associated(self%targsort)) deallocate(self%targsort)
        if (associated(self%chargesort)) deallocate(self%chargesort)
        if (associated(self%dipvecsort)) deallocate(self%dipvecsort)
        if (associated(self%potsort)) deallocate(self%potsort)
        if (associated(self%gradsort)) deallocate(self%gradsort)
        if (associated(self%hesssort)) deallocate(self%hesssort)
        if (associated(self%pottargsort)) deallocate(self%pottargsort)
        if (associated(self%gradtargsort)) deallocate(self%gradtargsort)
        if (associated(self%hesstargsort)) deallocate(self%hesstargsort)
        if (associated(self%nterms)) deallocate(self%nterms)
        if (associated(self%iaddr)) deallocate(self%iaddr)
        if (associated(self%mptemp)) deallocate(self%mptemp)
        if (associated(self%mptemp2)) deallocate(self%mptemp2)
        if (associated(self%rmlexp)) deallocate(self%rmlexp)    
    end subroutine dealloc



     subroutine lfmm3dmain_tree(self,nd,eps, &
     &     nsource,sourcesort, &
     &     ifcharge,chargesort, &
     &     ifdipole,dipvecsort, &
     &     ntarg,targsort,nexpc,expcsort, &
     &     iaddr,rmlexp,lmptot,mptemp,mptemp2,lmptemp, &
     &     itree,ltree,ipointer,ndiv,nlevels,  &
     &     nboxes,iper,boxsize,centers,isrcse,itargse,iexpcse, &
     &     rscales,laddr,nterms, &
     &     ifpgh,pot,grad,hess, &
     &     ifpghtarg,pottarg,gradtarg,hesstarg,ntj, &
     &     tsort,scjsort,ifnear,timeinfo,ier)
      implicit none
        class(FMM3DTree), intent(inout) :: self
      integer nd
      integer ier
      double precision eps
      integer nsource,ntarg,nexpc
      integer ndiv,nlevels

      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg

      double precision sourcesort(3,nsource)

      double precision chargesort(nd,*)
      double precision dipvecsort(nd,3,*)

      double precision targsort(3,ntarg)

      double precision pot(nd,*),grad(nd,3,*),hess(nd,6,*)
      double precision pottarg(nd,*),gradtarg(nd,3,*),hesstarg(nd,6,*)

      integer ntj
      integer ifnear
      double precision expcsort(3,nexpc)
      double complex tsort(nd,0:ntj,-ntj:ntj,nexpc)
      double precision scjsort(nexpc)

      integer nboxes
      integer *8 iaddr(2,nboxes), lmptot
      integer lmptemp
      double precision rmlexp(lmptot)
      double precision mptemp(lmptemp)
      double precision mptemp2(lmptemp)

      double precision thresh
       
      double precision timeinfo(6)
      double precision centers(3,nboxes)

      integer isep,iper
      integer laddr(2,0:nlevels)
      integer nterms(0:nlevels)
      integer *8 ipointer(8),ltree
      integer itree(ltree)
      !double precision rscales(0:nlevels)
      double precision, contiguous, pointer :: rscales(:)
      double precision boxsize(0:nlevels)
      integer isrcse(2,nboxes),itargse(2,nboxes),iexpcse(2,nboxes)
    !   integer, allocatable :: nlist1(:),list1(:,:)
    !   integer, allocatable :: nlist2(:),list2(:,:)
    !   integer, allocatable :: nlist3(:),list3(:,:)
    !   integer, allocatable :: nlist4(:),list4(:,:)

      integer, contiguous, pointer :: nlist1(:),list1(:,:)
      integer, contiguous, pointer :: nlist2(:),list2(:,:)
      integer, contiguous, pointer :: nlist3(:),list3(:,:)
      integer, contiguous, pointer :: nlist4(:),list4(:,:)

      integer nuall,ndall,nnall,nsall,neall,nwall
      integer nu1234,nd5678,nn1256,ns3478,ne1357,nw2468
      integer nn12,nn56,ns34,ns78,ne13,ne57,nw24,nw68
      integer ne1,ne3,ne5,ne7,nw2,nw4,nw6,nw8

      ! integer, allocatable :: uall(:,:),dall(:,:),nall(:,:)
      ! integer, allocatable :: sall(:,:),eall(:,:),wall(:,:)
      ! integer, allocatable :: u1234(:,:),d5678(:,:)
      ! integer, allocatable :: n1256(:,:),s3478(:,:)
      ! integer, allocatable :: e1357(:,:),w2468(:,:)
      ! integer, allocatable :: n12(:,:),n56(:,:),s34(:,:),s78(:,:)
      ! integer, allocatable :: e13(:,:),e57(:,:),w24(:,:),w68(:,:)
      ! integer, allocatable :: e1(:,:),e3(:,:),e5(:,:),e7(:,:)
      ! integer, allocatable :: w2(:,:),w4(:,:),w6(:,:),w8(:,:)

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
    !   double precision, allocatable :: carray(:,:), dc(:,:)
    !   double precision, allocatable :: cs(:,:),fact(:),rdplus(:,:,:)
    !   double precision, allocatable :: rdminus(:,:,:), rdsq3(:,:,:)
    !   double precision, allocatable :: rdmsq3(:,:,:)
      double precision, contiguous, pointer :: carray(:,:), dc(:,:)
      double precision, contiguous, pointer :: cs(:,:),fact(:),rdplus(:,:,:)
      double precision, contiguous, pointer :: rdminus(:,:,:), rdsq3(:,:,:)
      double precision, contiguous, pointer :: rdmsq3(:,:,:)
  
  
      !double precision, allocatable :: rlams(:),whts(:)
      double precision, contiguous, pointer :: rlams(:),whts(:)

      !double precision, allocatable :: rlsc(:,:,:)
      double precision, contiguous, pointer :: rlsc(:,:,:)
      !integer, allocatable :: nfourier(:), nphysical(:)
      integer, contiguous, pointer :: nfourier(:), nphysical(:)
      integer nexptot, nexptotp
    !   double complex, allocatable :: xshift(:,:)
    !   double complex, allocatable :: yshift(:,:)
    !   double precision, allocatable :: zshift(:,:)
      double complex, contiguous, pointer:: xshift(:,:)
      double complex, contiguous, pointer :: yshift(:,:)
      double precision, contiguous, pointer :: zshift(:,:)


    !   double complex, allocatable :: fexpe(:),fexpo(:),fexpback(:)
    !   double complex, allocatable :: mexp(:,:,:,:)
    !   double complex, allocatable :: mexpf1(:,:,:),mexpf2(:,:,:)
    !   double complex, allocatable :: &
    !  &    mexpp1(:,:,:),mexpp2(:,:,:),mexppall(:,:,:,:)

      double complex, contiguous, pointer :: fexpe(:),fexpo(:),fexpback(:)
      double complex, contiguous, pointer :: mexp(:,:,:,:)
      double complex, contiguous, pointer :: mexpf1(:,:,:),mexpf2(:,:,:)
      double complex, contiguous, pointer :: &
     &    mexpp1(:,:,:),mexpp2(:,:,:),mexppall(:,:,:,:)

    !   double complex, allocatable :: tmp(:,:,:,:)
    !   double precision, allocatable :: mptmp(:,:)

      double complex, contiguous, pointer:: tmp(:,:,:,:)
      double precision, contiguous, pointer :: mptmp(:,:)

      double precision sourcetmp(3)
      double complex chargetmp

      integer ix,iy,iz,ictr
      double precision rtmp
      double complex zmul

      integer nlege, lw7, lused7, itype
      !double precision wlege(40000)

      double precision, contiguous, pointer :: wlege(:)
      integer nterms_eval(4,0:nlevels)

      integer mnlist1, mnlist2,mnlist3,mnlist4,mnbors
      double complex eye, ztmp
      double precision alphaj
      integer ctr,nn,iptr1,iptr2
      !ouble precision, allocatable :: rscpow(:)
      double precision, contiguous, pointer :: rscpow(:)
      double precision pi,errtmp
      double complex ima

      double precision ctmp(3)

!     list 3 variables
      ! double complex, allocatable :: iboxlexp(:,:,:)
      ! double precision, allocatable :: iboxsubcenters(:,:,:)
      ! double precision, allocatable :: iboxpot(:,:,:)
      ! double precision, allocatable :: iboxgrad(:,:,:,:)
      ! double precision, allocatable :: iboxhess(:,:,:,:)
      ! double precision, allocatable :: iboxsrc(:,:,:)
      ! integer, allocatable :: iboxsrcind(:,:)
      ! integer, allocatable :: iboxfl(:,:,:)
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
      !integer, allocatable :: list4ct(:),ilist4(:)
      integer, contiguous, pointer :: list4ct(:),ilist4(:)
      ! double complex, allocatable :: gboxmexp(:,:,:)
      ! double complex, allocatable :: gboxwexp(:,:,:,:,:)
      ! double complex, allocatable :: pgboxwexp(:,:,:,:)
      ! double precision, allocatable :: gboxsubcenters(:,:,:)
      ! double precision, allocatable :: gboxsort(:,:,:)
      ! integer, allocatable :: gboxind(:,:)
      ! integer, allocatable :: gboxfl(:,:,:)
      ! double precision, allocatable :: gboxcgsort(:,:,:)
      ! double precision, allocatable :: gboxdpsort(:,:,:,:)


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
      !double precision, allocatable :: scarray(:,:)
      double precision, contiguous, pointer :: scarray(:,:)
      integer *8 bigint
      integer iert
      data ima/(0.0d0,1.0d0)/

      integer nthd,ithd
      integer omp_get_max_threads,omp_get_thread_num
      nthd = 1
!$    nthd=omp_get_max_threads()

!       pi = 4.0d0*atan(1.0d0)

!       thresh = 2.0d0**(-51)*boxsize(0)


! !     ifprint is an internal information printing flag. 
! !     Suppressed if ifprint=0.
! !     Prints timing breakdown and other things if ifprint=1.
! !     Prints timing breakdown, list information, 
! !     and other things if ifprint=2.
! !       
       ifprint=0
      
! !
! !   initialize various tree lists
! !
!       mnlist1 = 0
!       mnlist2 = 0
!       mnlist3 = 0
!       mnlist4 = 0
!       mnbors = 27

!       isep = 1
      
!       call computemnlists(nlevels,nboxes,itree(ipointer(1)),boxsize, &
!      &    centers,itree(ipointer(3)),itree(ipointer(4)), &
!      &    itree(ipointer(5)),isep,itree(ipointer(6)),mnbors, &
!      &    itree(ipointer(7)),iper,mnlist1,mnlist2,mnlist3,mnlist4)
      
!       allocate(list1(mnlist1,nboxes),nlist1(nboxes))
!       allocate(list2(mnlist2,nboxes),nlist2(nboxes))
!       allocate(list3(mnlist3,nboxes),nlist3(nboxes))
!       allocate(list4(mnlist4,nboxes),nlist4(nboxes))

!       call computelists(nlevels,nboxes,itree(ipointer(1)),boxsize, &
!      &    centers,itree(ipointer(3)),itree(ipointer(4)), &
!      &    itree(ipointer(5)),isep,itree(ipointer(6)),mnbors, &
!      &    itree(ipointer(7)),iper,nlist1,mnlist1,list1,nlist2, &
!      &    mnlist2,list2,nlist3,mnlist3,list3, &
!      &    nlist4,mnlist4,list4)
      

! !     Initialize routines for plane wave mp loc translation
 
!       if(isep.eq.1) then
!          if(eps.ge.0.5d-2) nlams = 12
!          if(eps.lt.0.5d-2.and.eps.ge.0.5d-3) nlams = 12
!          if(eps.lt.0.5d-3.and.eps.ge.0.5d-6) nlams = 20
!          if(eps.lt.0.5d-6.and.eps.ge.0.5d-9) nlams = 29
!          if(eps.lt.0.5d-9) nlams = 37
!       endif
!       if(isep.eq.2) then
!          if(eps.ge.0.5d-3) nlams = 9
!          if(eps.lt.0.5d-3.and.eps.ge.0.5d-6) nlams = 15
!          if(eps.lt.0.5d-6.and.eps.ge.0.5d-9) nlams = 22
!          if(eps.lt.0.5d-9) nlams = 29
!       endif
!       allocate(rlams(nlams),whts(nlams))
!       allocate(nphysical(nlams),nfourier(nlams))


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

        rscales => self%scales


        wlege => self%wlege

        lca = self%lca

        call self%reset_expansion_coeff()


        pgboxwexp => self%pgboxwexp
        gboxmexp => self%gboxmexp
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
        iboxpot => self%iboxpot
        iboxgrad => self%iboxgrad
        iboxhess => self%iboxhess

    !------------------------------------------------



      !  print *, " self%nmax=",self%nmax
    !   nmax = 0
    !   do i=0,nlevels
    !      if(nmax.lt.nterms(i)) nmax = nterms(i)
    !   enddo

    !  ! print *, " nmax=",nmax
    !     nmaxt = nmax
    !   allocate(rscpow(0:nmax))
    !   allocate(carray(4*nmax+1,4*nmax+1))
    !   allocate(dc(0:4*nmax,0:4*nmax))
    !   allocate(rdplus(0:nmax,0:nmax,-nmax:nmax))
    !   allocate(rdminus(0:nmax,0:nmax,-nmax:nmax))
    !   allocate(rdsq3(0:nmax,0:nmax,-nmax:nmax))
    !   allocate(rdmsq3(0:nmax,0:nmax,-nmax:nmax))
    !   allocate(rlsc(0:nmax,0:nmax,nlams))


!     generate rotation matrices and carray
     ! call getpwrotmat(nmax,carray,rdplus,rdminus,rdsq3,rdmsq3,dc)


!     generate rlams and weights (these are the nodes
!     and weights for the lambda integral)

     ! call vwts(rlams,whts,nlams)


!     generate the number of fourier modes required to represent the
!     moment function in fourier space

      !call numthetahalf(nfourier,nlams)
 
!     generate the number of fourier modes in physical space
!     required for the exponential representation
     ! call numthetafour(nphysical,nlams)

!     Generate powers of lambda for the exponential basis
     ! call rlscini(rlsc,nlams,rlams,nmax)

!
!
!
    !   nn = 10*(nmax+2)**2
    !   allocate(scarray(nn,0:nlevels))
    !   do ilev=0,nlevels
    !     call l3dtaevalhessdini(nterms(ilev),scarray(1,ilev))
    !   enddo

!     Compute total number of plane waves
    !   nexptotp = 0
    !   nexptot = 0
    !   nthmax = 0
    !   nphmax = 0
    !   nn = 0
    !   do i=1,nlams
    !      nexptot = nexptot + nfourier(i)
    !      nexptotp = nexptotp + nphysical(i)
    !      if(nfourier(i).gt.nthmax) nthmax = nfourier(i)
    !      if(nphysical(i).gt.nphmax) nphmax = nphysical(i)
    !      nn = nn + nphysical(i)*nfourier(i)
    !   enddo


    !   allocate(fexpe(nn),fexpo(nn),fexpback(nn))
    !   allocate(tmp(nd,0:nmax,-nmax:nmax,nthd))
    !   allocate(mptmp(lmptemp,nthd))

    !   allocate(xshift(-5:5,nexptotp))
    !   allocate(yshift(-5:5,nexptotp))
    !   allocate(zshift(5,nexptotp))

    !   allocate(mexpf1(nd,nexptot,nthd),mexpf2(nd,nexptot,nthd), &
    !  &    mexpp1(nd,nexptotp,nthd))
    !   allocate(mexpp2(nd,nexptotp,nthd),mexppall(nd,nexptotp,16,nthd))

!
!c      NOTE: there can be some memory savings here
!
      ! bigint = 0
      ! bigint = nboxes
      ! bigint = bigint*6
      ! bigint = bigint*nexptotp*nd

      ! if(ifprint.ge.1) print *, "mexp memory=",bigint/1.0d9

      ! allocate(mexp(nd,nexptotp,nboxes,6),stat=iert)
      ! if(iert.ne.0) then
      !   print *, "Cannot allocate pw expansion workspace"
      !   print *, "bigint=", bigint
      !   ier = 8
      !   return
      ! endif

      ! allocate(list4ct(nboxes))
      ! allocate(ilist4(nboxes))
      ! do i=1,nboxes
      !   list4ct(i)=0
      !   ilist4(i)=0
      ! enddo
      ! cntlist4=0

!     Precompute table for shifting exponential coefficients in 
!     physical domain
!      call mkexps(rlams,nlams,nphysical,nexptotp,xshift,yshift,zshift)

!     Precompute table of exponentials for mapping from
!     fourier to physical domain
!      call mkfexp(nlams,nfourier,nphysical,fexpe,fexpo,fexpback)
      
!
!c    compute array of factorials

     
      ! nmax2 = 2*nmax
      ! allocate(fact(0:nmax2),cs(0:nmax,-nmax:nmax))
      
      ! d = 1
      ! fact(0) = d
      ! do i=1,nmax2
      !   d=d*sqrt(i+0.0d0)
      !   fact(i) = d
      ! enddo

      ! cs(0,0) = 1.0d0
      ! do l=1,nmax
      !   do m=0,l
      !     cs(l,m) = ((-1)**l)/(fact(l-m)*fact(l+m))
      !     cs(l,-m) = cs(l,m)
      !   enddo
      ! enddo
      
!       if(ifprint.ge.1)  &
!      &    call prin2('end of generating plane wave info*',i,0)
! !
!
!     ... set the expansion coefficients to zero
!
! !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,idim)
!       do i=1,nexpc
!         do k=-ntj,ntj
!           do j = 0,ntj
!             do idim=1,nd
!               tsort(idim,j,k,i)=0
!             enddo
!           enddo
!         enddo
!       enddo
! !$OMP END PARALLEL DO

!       
      ! do i=1,6
      !   timeinfo(i)=0
      ! enddo

!
!       ... set all multipole and local expansions to zero
!

!       do ilev = 0,nlevels
! !$OMP PARALLEL DO DEFAULT(SHARED) &
! !$OMP PRIVATE(ibox)
!         do ibox=laddr(1,ilev),laddr(2,ilev)
!           call mpzero(nd,rmlexp(iaddr(1,ibox)),nterms(ilev))
!           call mpzero(nd,rmlexp(iaddr(2,ibox)),nterms(ilev))
!         enddo
! !$OMP END PARALLEL DO
!       enddo

!
!      set scjsort
!
!       do ilev=0,nlevels
! !$OMP PARALLEL DO DEFAULT(SHARED) &
! !$OMP PRIVATE(ibox,nchild,istart,iend,i)
!          do ibox=laddr(1,ilev),laddr(2,ilev)
!             nchild = itree(ipointer(4)+ibox-1)
!             if(nchild.gt.0) then
!                istart = iexpcse(1,ibox)
!                iend = iexpcse(2,ibox) 
!                do i=istart,iend
!                   scjsort(i) = rscales(ilev)
!                enddo
!             endif
!          enddo
! !$OMP END PARALLEL DO
!       enddo


!    initialize legendre function evaluation routines
      ! nlege = 100
      ! lw7 = 40000
      ! call ylgndrfwini(nlege,wlege,lw7,lused7)

!
!     count number of boxes are in list4
     ! lca = 4*nmax
    !   if(ifprint.ge.1) &
    !  &   call prinf('=== STEP 0 list4===*',i,0)
    !   call cpu_time(time1)
!$    time1=omp_get_wtime()
      ! do ilev=1,nlevels-1
      !    do ibox=laddr(1,ilev),laddr(2,ilev)
      !       if(nlist3(ibox).gt.0) then
      !         cntlist4=cntlist4+1
      !         list4ct(ibox)=cntlist4
      !         ilist4(cntlist4)=ibox
      !       endif
      !    enddo
      ! enddo
      ! if(ifprint.ge.1) print *,"nboxes:",nboxes,"cntlist4:",cntlist4
    !   allocate(pgboxwexp(nd,nexptotp,cntlist4,6))
    !   allocate(gboxmexp(nd*(nterms(nlevels)+1)* &
    !  &    (2*nterms(nlevels)+1),8,cntlist4))



    !   allocate(gboxsubcenters(3,8,nthd))
    !   allocate(gboxfl(2,8,nthd))

!       nmaxt = 0
! !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,istart,iend,npts) &
! !$OMP REDUCTION(max:nmaxt)
!       do ibox=1,nboxes
!         if(list4ct(ibox).gt.0) then
!           istart = isrcse(1,ibox)
!           iend = isrcse(2,ibox)
!           npts = iend-istart+1
!           if(npts.gt.nmaxt) nmaxt = npts
!         endif
!       enddo
! !$OMP END PARALLEL DO

!       allocate(gboxind(nmaxt,nthd))
!       allocate(gboxsort(3,nmaxt,nthd))
!       allocate(gboxwexp(nd,nexptotp,6,8,nthd))
!       allocate(gboxcgsort(nd,nmaxt,nthd))
!       allocate(gboxdpsort(nd,3,nmaxt,nthd))

! !   note gboxmexp is an array not scalar
!       pgboxwexp=0d0
!       gboxmexp=0d0

!     form mexp for all list4 type box at first ghost box center
      do ilev=1,nlevels-1

         rscpow(0) = 1.0d0/boxsize(ilev+1)
         rtmp = rscales(ilev+1)/boxsize(ilev+1)
         do i=1,nterms(ilev+1)
            rscpow(i) = rscpow(i-1)*rtmp
         enddo

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,istart,iend,jbox,jstart,jend,npts,npts0,i) &
!$OMP PRIVATE(ithd)
         do ibox=laddr(1,ilev),laddr(2,ilev)
            ithd = 0
!$          ithd=omp_get_thread_num()
            ithd = ithd + 1
            if(list4ct(ibox).gt.0) then
              istart=isrcse(1,ibox)
              iend=isrcse(2,ibox)
              npts = iend-istart+1

              if(npts.gt.0) then
                call subdividebox(sourcesort(1,istart),npts, &
     &    centers(1,ibox),boxsize(ilev+1), &
     &    gboxind(1,ithd),gboxfl(1,1,ithd), &
     &    gboxsubcenters(1,1,ithd))
                call dreorderf(3,npts,sourcesort(1,istart), &
     &    gboxsort(1,1,ithd),gboxind(1,ithd))
                if(ifcharge.eq.1) then
                  call dreorderf(nd,npts,chargesort(1,istart), &
     &    gboxcgsort(1,1,ithd),gboxind(1,ithd))
                endif
                if(ifdipole.eq.1) then
                  call dreorderf(3*nd,npts,dipvecsort(1,1,istart), &
     &    gboxdpsort(1,1,1,ithd),gboxind(1,ithd))
                endif
                do i=1,8
                  if(gboxfl(1,i,ithd).gt.0) then
                    jstart=gboxfl(1,i,ithd)
                    jend=gboxfl(2,i,ithd)
                    npts0=jend-jstart+1
                    jbox=list4ct(ibox)
                    if(ifcharge.eq.1.and.ifdipole.eq.0) then
                      call l3dformmpc(nd,rscales(ilev+1), &
     &    gboxsort(1,jstart,ithd), &
     &    gboxcgsort(1,jstart,ithd), &
     &    npts0,gboxsubcenters(1,i,ithd),nterms(ilev+1), &
     &    gboxmexp(1,i,jbox),wlege,nlege)          
                    endif
                    if(ifcharge.eq.0.and.ifdipole.eq.1) then
                      call l3dformmpd(nd,rscales(ilev+1), &
     &    gboxsort(1,jstart,ithd), &
     &    gboxdpsort(1,1,jstart,ithd), &
     &    npts0,gboxsubcenters(1,i,ithd),nterms(ilev+1), &
     &    gboxmexp(1,i,jbox),wlege,nlege)          
                    endif
                    if(ifcharge.eq.1.and.ifdipole.eq.1) then
                      call l3dformmpcd(nd,rscales(ilev+1), &
     &    gboxsort(1,jstart,ithd), &
     &    gboxcgsort(1,jstart,ithd), &
     &    gboxdpsort(1,1,jstart,ithd), &
     &    npts0,gboxsubcenters(1,i,ithd),nterms(ilev+1), &
     &    gboxmexp(1,i,jbox),wlege,nlege)          
                    endif
                    call l3dmpmp(nd,rscales(ilev+1), &
     &    gboxsubcenters(1,i,ithd),gboxmexp(1,i,jbox), &
     &    nterms(ilev+1),rscales(ilev),centers(1,ibox), &
     &    rmlexp(iaddr(1,ibox)),nterms(ilev),dc,lca)
     
                    call mpscale(nd,nterms(ilev+1),gboxmexp(1,i,jbox), &
     &    rscpow,tmp(1,0,-nmax,ithd))
!
!c                process up down for current box
!
                    call mpoletoexp(nd,tmp(1,0,-nmax,ithd), &
     &    nterms(ilev+1),nlams, &
     &    nfourier,nexptot,mexpf1(1,1,ithd), &
     &    mexpf2(1,1,ithd),rlsc)

                    call ftophys(nd,mexpf1(1,1,ithd), &
     &    nlams,rlams,nfourier, &
     &    nphysical,nthmax,gboxwexp(1,1,1,i,ithd), &
     &    fexpe,fexpo)

                    call ftophys(nd,mexpf2(1,1,ithd), &
     &    nlams,rlams,nfourier, &
     &    nphysical,nthmax,gboxwexp(1,1,2,i,ithd), &
     &    fexpe,fexpo)

                    call processgboxudexp(nd,gboxwexp(1,1,1,i,ithd), &
     &    gboxwexp(1,1,2,i,ithd),i,nexptotp, &
     &    pgboxwexp(1,1,jbox,1),pgboxwexp(1,1,jbox,2), &
     &    xshift,yshift,zshift)
!
!c                process north-south for current box
!
                    call rotztoy(nd,nterms(ilev+1),tmp(1,0,-nmax,ithd), &
     &    mptmp(1,ithd),rdminus)
                    call mpoletoexp(nd,mptmp(1,ithd), &
     &    nterms(ilev+1),nlams, &
     &    nfourier,nexptot,mexpf1(1,1,ithd), &
     &    mexpf2(1,1,ithd),rlsc)

                    call ftophys(nd,mexpf1(1,1,ithd), &
     &    nlams,rlams,nfourier, &
     &    nphysical,nthmax,gboxwexp(1,1,3,i,ithd), &
     &    fexpe,fexpo)

                    call ftophys(nd,mexpf2(1,1,ithd), &
     &    nlams,rlams,nfourier, &
     &    nphysical,nthmax,gboxwexp(1,1,4,i,ithd), &
     &    fexpe,fexpo)

                    call processgboxnsexp(nd,gboxwexp(1,1,3,i,ithd), &
     &    gboxwexp(1,1,4,i,ithd),i,nexptotp, &
     &    pgboxwexp(1,1,jbox,3),pgboxwexp(1,1,jbox,4), &
     &    xshift,yshift,zshift)
!
!c                process east-west for current box
!
                    call rotztox(nd,nterms(ilev+1),tmp(1,0,-nmax,ithd), &
     &    mptmp(1,ithd),rdplus)
                    call mpoletoexp(nd,mptmp(1,ithd), &
     &    nterms(ilev+1),nlams, &
     &    nfourier,nexptot,mexpf1(1,1,ithd), &
     &    mexpf2(1,1,ithd),rlsc)

                    call ftophys(nd,mexpf1(1,1,ithd), &
     &    nlams,rlams,nfourier, &
     &    nphysical,nthmax,gboxwexp(1,1,5,i,ithd), &
     &    fexpe,fexpo)

                    call ftophys(nd,mexpf2(1,1,ithd), &
     &    nlams,rlams,nfourier, &
     &    nphysical,nthmax,gboxwexp(1,1,6,i,ithd), &
     &    fexpe,fexpo)
                
                    call processgboxewexp(nd,gboxwexp(1,1,5,i,ithd), &
     &    gboxwexp(1,1,6,i,ithd),i,nexptotp, &
     &    pgboxwexp(1,1,jbox,5),pgboxwexp(1,1,jbox,6), &
     &    xshift,yshift,zshift)
                  endif
                enddo
              endif
            endif
         enddo
!$OMP END PARALLEL DO
      enddo
      !deallocate(gboxfl,gboxsubcenters,gboxwexp,gboxcgsort)
      !deallocate(gboxdpsort,gboxind,gboxsort)

      !call cpu_time(time2)
!$    time2=omp_get_wtime()
      !if(ifprint.ge.1) print *,"mexp list4 time:",time2-time1
      !timeinfo(3)=time2-time1
!     end of count number of boxes are in list4
!

!
!
    !   if(ifprint .ge. 1)  &
    !  &   call prinf('=== STEP 1 (form mp) ====*',i,0)
    !     call cpu_time(time1)
!$        time1=omp_get_wtime()
!
!       ... step 1, locate all charges, assign them to boxes, and
!       form multipole expansions


      do ilev=2,nlevels
         if(ifcharge.eq.1.and.ifdipole.eq.0) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,npts,istart,iend,nchild)
            do ibox=laddr(1,ilev),laddr(2,ilev)

               istart = isrcse(1,ibox)
               iend = isrcse(2,ibox) 
               npts = iend-istart+1

               nchild = itree(ipointer(4)+ibox-1)

               if(npts.gt.0.and.nchild.eq.0.and.list4ct(ibox).eq.0) then
                  call l3dformmpc(nd,rscales(ilev), &
     &    sourcesort(1,istart),chargesort(1,istart),npts, &
     &    centers(1,ibox),nterms(ilev), &
     &    rmlexp(iaddr(1,ibox)),wlege,nlege)          
               endif
            enddo
!$OMP END PARALLEL DO
         endif

         if(ifcharge.eq.0.and.ifdipole.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,npts,istart,iend,nchild)
            do ibox=laddr(1,ilev),laddr(2,ilev)

               istart = isrcse(1,ibox) 
               iend = isrcse(2,ibox) 
               npts = iend-istart+1

               nchild = itree(ipointer(4)+ibox-1)

               if(npts.gt.0.and.nchild.eq.0.and.list4ct(ibox).eq.0) then
                  call l3dformmpd(nd,rscales(ilev), &
     &    sourcesort(1,istart), &
     &    dipvecsort(1,1,istart),npts, &
     &    centers(1,ibox),nterms(ilev), &
     &    rmlexp(iaddr(1,ibox)),wlege,nlege)          
               endif
            enddo
!$OMP END PARALLEL DO
         endif

         if(ifdipole.eq.1.and.ifcharge.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,npts,istart,iend,nchild)
            do ibox=laddr(1,ilev),laddr(2,ilev)

               istart = isrcse(1,ibox) 
               iend = isrcse(2,ibox)
               npts = iend-istart+1

               nchild = itree(ipointer(4)+ibox-1)

               if(npts.gt.0.and.nchild.eq.0.and.list4ct(ibox).eq.0) then
                  call l3dformmpcd(nd,rscales(ilev), &
     &    sourcesort(1,istart),chargesort(1,istart), &
     &    dipvecsort(1,1,istart),npts, &
     &    centers(1,ibox),nterms(ilev), &
     &    rmlexp(iaddr(1,ibox)),wlege,nlege)          
               endif
            enddo
!$OMP END PARALLEL DO
         endif
      enddo
      if(ifprint.ge.1) print *,"nboxes:",nboxes,"leaf:",cntlist4



      call cpu_time(time2)
!$    time2=omp_get_wtime()
      timeinfo(1)=time2-time1

      lca = 4*nmax


!       
      if(ifprint .ge. 1) &
     &      call prinf('=== STEP 2 (merge mp) ====*',i,0)
      call cpu_time(time1)
!$    time1=omp_get_wtime()
!
      do ilev=nlevels-1,0,-1
!$OMP PARALLEL DO DEFAULT(SHARED) &
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
!$OMP END PARALLEL DO
      enddo

      call cpu_time(time2)
!$    time2=omp_get_wtime()
      timeinfo(2)=time2-time1

      if(ifprint.ge.1) &
     &    call prinf('=== Step 3 (mp to loc+formta+mpeval) ===*',i,0)
!      ... step 3, convert multipole expansions into local
!       expansions

      call cpu_time(time1)
!$        time1=omp_get_wtime()

!
!c     zero out mexp
! 

! !$OMP PARALLEL DO DEFAULT(SHARED) &
! !$OMP PRIVATE(i,j,k,idim)
!       do k=1,6
!         do i=1,nboxes
!           do j=1,nexptotp
!             do idim=1,nd
!               mexp(idim,j,i,k) = 0.0d0
!             enddo
!           enddo
!         enddo
!       enddo
! !$OMP END PARALLEL DO

!     init uall,dall,...,etc arrays
      ! allocate(uall(200,nthd),dall(200,nthd),nall(120,nthd))
      ! allocate(sall(120,nthd),eall(72,nthd),wall(72,nthd))
      ! allocate(u1234(36,nthd),d5678(36,nthd),n1256(24,nthd))
      ! allocate(s3478(24,nthd))
      ! allocate(e1357(16,nthd),w2468(16,nthd),n12(20,nthd))
      ! allocate(n56(20,nthd),s34(20,nthd),s78(20,nthd))
      ! allocate(e13(20,nthd),e57(20,nthd),w24(20,nthd),w68(20,nthd))
      ! allocate(e1(20,nthd),e3(5,nthd),e5(5,nthd),e7(5,nthd))
      ! allocate(w2(5,nthd),w4(5,nthd),w6(5,nthd),w8(5,nthd))
      ! allocate(iboxsubcenters(3,8,nthd))
      ! allocate(iboxfl(2,8,nthd))
!
!  figure out allocations needed for iboxsrc,iboxsrcind,iboxpot
!  and so on
!
!       nmaxt = 0
! !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,istart,iend,npts) &
! !$OMP REDUCTION(max:nmaxt)
!       do ibox=1,nboxes
!         if(nlist3(ibox).gt.0) then
!           istart = isrcse(1,ibox)
!           iend = isrcse(2,ibox)
!           npts = iend-istart+1
!           if(npts.gt.nmaxt) nmaxt = npts

!           istart = itargse(1,ibox)
!           iend = itargse(2,ibox)
!           npts = iend - istart + 1
!           if(npts.gt.nmaxt) nmaxt = npts
!         endif
!       enddo
! !$OMP END PARALLEL DO

!       allocate(iboxsrcind(nmaxt,nthd))
!       allocate(iboxsrc(3,nmaxt,nthd))
!       allocate(iboxpot(nd,nmaxt,nthd))
!       allocate(iboxgrad(nd,3,nmaxt,nthd))
!       allocate(iboxhess(nd,6,nmaxt,nthd))

      do ilev=2,nlevels
        allocate(iboxlexp(nd*(nterms(ilev)+1)* &
     &    (2*nterms(ilev)+1),8,nthd))

         rscpow(0) = 1.0d0/boxsize(ilev)
         rtmp = rscales(ilev)/boxsize(ilev)
         do i=1,nterms(ilev)
            rscpow(i) = rscpow(i-1)*rtmp
         enddo

!$OMP PARALLEL DO DEFAULT (SHARED) &
!$OMP PRIVATE(ibox,istart,iend,npts) &
!$OMP PRIVATE(ithd)
         do ibox=laddr(1,ilev),laddr(2,ilev)
            ithd = 0
!$          ithd=omp_get_thread_num()
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
!$OMP END PARALLEL DO
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
!$OMP PARALLEL DO DEFAULT (SHARED) &
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
           ithd = 0
!$         ithd=omp_get_thread_num()
           ithd = ithd + 1
           npts = 0
           if(ifpghtarg.gt.0) then
             istart = itargse(1,ibox)
             iend = itargse(2,ibox) 
             npts = npts + iend-istart+1
           endif

           istart = iexpcse(1,ibox) 
           iend = iexpcse(2,ibox) 
           npts = npts + iend-istart+1

           nchild = itree(ipointer(4)+ibox-1)

           if(ifpgh.gt.0) then
             istart = isrcse(1,ibox) 
             iend = isrcse(2,ibox) 
             npts = npts + iend-istart+1
           endif


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
              do i=1,8
                call mpzero(nd,iboxlexp(1,i,ithd),nterms(ilev))
              enddo

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

              if(ifpgh.eq.1) then
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
                  call dreorderf(nd,npts,pot(1,istart), &
     &    iboxpot(1,1,ithd),iboxsrcind(1,ithd))
                  do i=1,8
                    if(iboxfl(1,i,ithd).gt.0) then
                      jstart=iboxfl(1,i,ithd)
                      jend=iboxfl(2,i,ithd)
                      npts0=jend-jstart+1
                      if(npts0.gt.0) then
                        call l3dtaevalp(nd,rscales(ilev), &
     &    iboxsubcenters(1,i,ithd),iboxlexp(1,i,ithd), &
     &    nterms(ilev),iboxsrc(1,jstart,ithd),npts0, &
     &    iboxpot(1,jstart,ithd),wlege,nlege)
                      endif
                    endif
                  enddo
                  call dreorderi(nd,npts,iboxpot(1,1,ithd), &
     &    pot(1,istart),iboxsrcind(1,ithd))
                endif
              endif

              if(ifpgh.eq.2) then
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
                  call dreorderf(nd,npts,pot(1,istart), &
     &    iboxpot(1,1,ithd),iboxsrcind(1,ithd))
                  call dreorderf(3*nd,npts,grad(1,1,istart), &
     &    iboxgrad(1,1,1,ithd),iboxsrcind(1,ithd))
                  do i=1,8
                    if(iboxfl(1,i,ithd).gt.0) then
                      jstart=iboxfl(1,i,ithd)
                      jend=iboxfl(2,i,ithd)
                      npts0=jend-jstart+1
                      if(npts0.gt.0) then
                        call l3dtaevalg(nd,rscales(ilev), &
     &    iboxsubcenters(1,i,ithd),iboxlexp(1,i,ithd), &
     &    nterms(ilev),iboxsrc(1,jstart,ithd),npts0, &
     &    iboxpot(1,jstart,ithd), &
     &    iboxgrad(1,1,jstart,ithd),wlege,nlege)
                      endif
                    endif
                  enddo
                  call dreorderi(nd,npts,iboxpot(1,1,ithd), &
     &    pot(1,istart),iboxsrcind(1,ithd))
                  call dreorderi(3*nd,npts,iboxgrad(1,1,1,ithd), &
     &    grad(1,1,istart),iboxsrcind(1,ithd))
                endif
              endif
!
!  continue from here
!
              

              if(ifpgh.eq.3) then
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
                  call dreorderf(nd,npts,pot(1,istart), &
     &    iboxpot(1,1,ithd),iboxsrcind(1,ithd))
                  call dreorderf(3*nd,npts,grad(1,1,istart), &
     &    iboxgrad(1,1,1,ithd),iboxsrcind(1,ithd))
                  call dreorderf(6*nd,npts,hess(1,1,istart), &
     &    iboxhess(1,1,1,ithd),iboxsrcind(1,ithd))
           
                  do i=1,8
                    if(iboxfl(1,i,ithd).gt.0) then
                      jstart=iboxfl(1,i,ithd)
                      jend=iboxfl(2,i,ithd)
                      npts0=jend-jstart+1
                      if(npts0.gt.0) then
                        call l3dtaevalh(nd,rscales(ilev), &
     &    iboxsubcenters(1,i,ithd),iboxlexp(1,i,ithd), &
     &    nterms(ilev),iboxsrc(1,jstart,ithd),npts0, &
     &    iboxpot(1,jstart,ithd), &
     &    iboxgrad(1,1,jstart,ithd), &
     &    iboxhess(1,1,jstart,ithd),scarray(1,ilev))
                      endif
                    endif
                  enddo
                  call dreorderi(nd,npts,iboxpot(1,1,ithd), &
     &    pot(1,istart),iboxsrcind(1,ithd))
                  call dreorderi(3*nd,npts,iboxgrad(1,1,1,ithd), &
     &    grad(1,1,istart),iboxsrcind(1,ithd))
                  call dreorderi(6*nd,npts,iboxhess(1,1,1,ithd), &
     &    hess(1,1,istart),iboxsrcind(1,ithd))
                endif
              endif


              if(ifpghtarg.eq.1) then
                istart = itargse(1,ibox) 
                iend = itargse(2,ibox) 
                npts = iend-istart+1
                if(npts.gt.0) then
                  call subdividebox(targsort(1,istart),npts, &
     &    centers(1,ibox),boxsize(ilev), &
     &    iboxsrcind(1,ithd),iboxfl(1,1,ithd), &
     &    iboxsubcenters(1,1,ithd))
                  call dreorderf(3,npts,targsort(1,istart), &
     &    iboxsrc(1,1,ithd),iboxsrcind(1,ithd))
                  call dreorderf(nd,npts,pottarg(1,istart), &
     &    iboxpot(1,1,ithd),iboxsrcind(1,ithd))
                  do i=1,8
                    if(iboxfl(1,i,ithd).gt.0) then
                      jstart=iboxfl(1,i,ithd)
                      jend=iboxfl(2,i,ithd)
                      npts0=jend-jstart+1
                      if(npts0.gt.0) then
                        call l3dtaevalp(nd,rscales(ilev), &
     &    iboxsubcenters(1,i,ithd),iboxlexp(1,i,ithd), &
     &    nterms(ilev),iboxsrc(1,jstart,ithd),npts0, &
     &    iboxpot(1,jstart,ithd),wlege,nlege)
                      endif
                    endif
                  enddo
                  call dreorderi(nd,npts,iboxpot(1,1,ithd), &
     &    pottarg(1,istart),iboxsrcind(1,ithd))
                endif
              endif

              if(ifpghtarg.eq.2) then
                istart = itargse(1,ibox) 
                iend = itargse(2,ibox) 
                npts = iend-istart+1
                if(npts.gt.0) then
                  call subdividebox(targsort(1,istart),npts, &
     &    centers(1,ibox),boxsize(ilev), &
     &    iboxsrcind(1,ithd),iboxfl(1,1,ithd), &
     &    iboxsubcenters(1,1,ithd))
                  call dreorderf(3,npts,targsort(1,istart), &
     &    iboxsrc(1,1,ithd),iboxsrcind(1,ithd))
                  call dreorderf(nd,npts,pottarg(1,istart), &
     &    iboxpot(1,1,ithd),iboxsrcind(1,ithd))
                  call dreorderf(3*nd,npts,gradtarg(1,1,istart), &
     &    iboxgrad(1,1,1,ithd),iboxsrcind(1,ithd))
                  do i=1,8
                    if(iboxfl(1,i,ithd).gt.0) then
                      jstart=iboxfl(1,i,ithd)
                      jend=iboxfl(2,i,ithd)
                      npts0=jend-jstart+1
                      if(npts0.gt.0) then
                        call l3dtaevalg(nd,rscales(ilev), &
     &    iboxsubcenters(1,i,ithd),iboxlexp(1,i,ithd), &
     &    nterms(ilev),iboxsrc(1,jstart,ithd),npts0, &
     &    iboxpot(1,jstart,ithd), &
     &    iboxgrad(1,1,jstart,ithd),wlege,nlege)
                      endif
                    endif
                  enddo
                  call dreorderi(nd,npts,iboxpot(1,1,ithd), &
     &    pottarg(1,istart),iboxsrcind(1,ithd))
                  call dreorderi(3*nd,npts,iboxgrad(1,1,1,ithd), &
     &    gradtarg(1,1,istart),iboxsrcind(1,ithd))
                endif
              endif

              if(ifpghtarg.eq.3) then
                istart = itargse(1,ibox) 
                iend = itargse(2,ibox) 
                npts = iend-istart+1
                if(npts.gt.0) then
                  call subdividebox(targsort(1,istart),npts, &
     &    centers(1,ibox),boxsize(ilev), &
     &    iboxsrcind(1,ithd),iboxfl(1,1,ithd), &
     &    iboxsubcenters(1,1,ithd))
                  call dreorderf(3,npts,targsort(1,istart), &
     &    iboxsrc(1,1,ithd),iboxsrcind(1,ithd))
                  call dreorderf(nd,npts,pottarg(1,istart), &
     &    iboxpot(1,1,ithd),iboxsrcind(1,ithd))
                  call dreorderf(3*nd,npts,gradtarg(1,1,istart), &
     &    iboxgrad(1,1,1,ithd),iboxsrcind(1,ithd))
                  call dreorderf(6*nd,npts,hesstarg(1,1,istart), &
     &    iboxhess(1,1,1,ithd),iboxsrcind(1,ithd))
                  do i=1,8
                    if(iboxfl(1,i,ithd).gt.0) then
                      jstart=iboxfl(1,i,ithd)
                      jend=iboxfl(2,i,ithd)
                      npts0=jend-jstart+1
                      if(npts0.gt.0) then
                        call l3dtaevalh(nd,rscales(ilev), &
     &    iboxsubcenters(1,i,ithd),iboxlexp(1,i,ithd), &
     &    nterms(ilev),iboxsrc(1,jstart,ithd),npts0, &
     &    iboxpot(1,jstart,ithd), &
     &    iboxgrad(1,1,jstart,ithd), &
     &    iboxhess(1,1,jstart,ithd),scarray(1,ilev))
                      endif
                    endif
                  enddo
                  call dreorderi(nd,npts,iboxpot(1,1,ithd), &
     &    pottarg(1,istart),iboxsrcind(1,ithd))
                  call dreorderi(3*nd,npts,iboxgrad(1,1,1,ithd), &
     &    gradtarg(1,1,istart),iboxsrcind(1,ithd))
                  call dreorderi(6*nd,npts,iboxhess(1,1,1,ithd), &
     &    hesstarg(1,1,istart),iboxsrcind(1,ithd))
                endif
              endif

            endif
         enddo
!$OMP END PARALLEL DO
        deallocate(iboxlexp)  
      enddo

      ! deallocate(iboxsrcind,iboxsrc,iboxpot,iboxgrad,iboxhess)
      ! deallocate(iboxsubcenters,iboxfl)
      ! deallocate(uall,dall,nall,sall,eall,wall)
      ! deallocate(u1234,d5678,n1256,s3478)
      ! deallocate(e1357,w2468,n12,n56,s34,s78)
      ! deallocate(e13,e57,w24,w68)
      ! deallocate(e1,e3,e5,e7,w2,w4,w6,w8)
      ! deallocate(tmp,mptmp)
      ! call cpu_time(time2)
!$        time2=omp_get_wtime()
      ! timeinfo(3) = timeinfo(3) + time2-time1


    !   if(ifprint.ge.1) &
    !  &    call prinf('=== Step 4 (split loc) ===*',i,0)

    !   call cpu_time(time1)
!$        time1=omp_get_wtime()
      do ilev = 2,nlevels-1

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,i,jbox,istart,iend,npts)
         do ibox = laddr(1,ilev),laddr(2,ilev)

            npts = 0

            if(ifpghtarg.gt.0) then
               istart = itargse(1,ibox)
               iend = itargse(2,ibox) 
               npts = npts + iend-istart+1
            endif

            istart = iexpcse(1,ibox) 
            iend = iexpcse(2,ibox) 
            npts = npts + iend-istart+1

            if(ifpgh.gt.0) then
               istart = isrcse(1,ibox)
               iend = isrcse(2,ibox) 
               npts = npts + iend-istart+1
            endif

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
!$OMP END PARALLEL DO
      enddo
      ! call cpu_time(time2)
!$        time2=omp_get_wtime()
      ! timeinfo(4) = time2-time1


    !   if(ifprint.ge.1) &
    !  &    call prinf('=== step 5 (eval lo) ===*',i,0)

!     ... step 6, evaluate all local expansions
!

      ! call cpu_time(time1)
!$        time1=omp_get_wtime()
!

!
!c       shift local expansion to local epxanion at expansion centers
!        (note: this part is not relevant for particle codes.
!        it is relevant only for qbx codes)

      do ilev = 0,nlevels
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,nchild,istart,iend,i) &
!$OMP SCHEDULE(DYNAMIC)
         do ibox = laddr(1,ilev),laddr(2,ilev)
            nchild=itree(ipointer(4)+ibox-1)
            if(nchild.eq.0) then 
               istart = iexpcse(1,ibox) 
               iend = iexpcse(2,ibox) 
               do i=istart,iend

                  call l3dlocloc(nd,rscales(ilev), &
     &    centers(1,ibox),rmlexp(iaddr(2,ibox)), &
     &    nterms(ilev),rscales(ilev),expcsort(1,i), &
     &    tsort(1,0,-ntj,i),ntj,dc,lca)
               enddo
            endif
         enddo
!$OMP END PARALLEL DO
      enddo






!
!c        evaluate local expansion at source and target
!         locations
!
      do ilev = 0,nlevels
        if(ifpgh.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,nchild,istart,iend,i,npts) &
!$OMP SCHEDULE(DYNAMIC)
          do ibox = laddr(1,ilev),laddr(2,ilev)
            nchild=itree(ipointer(4)+ibox-1)
            if(nchild.eq.0) then 
              istart = isrcse(1,ibox) 
              iend = isrcse(2,ibox)
              npts = iend-istart+1
              call l3dtaevalp(nd,rscales(ilev),centers(1,ibox), &
     &    rmlexp(iaddr(2,ibox)),nterms(ilev),sourcesort(1,istart), &
     &    npts,pot(1,istart),wlege,nlege)
            endif
          enddo
!$OMP END PARALLEL DO
        endif

        if(ifpgh.eq.2) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,nchild,istart,iend,i,npts) &
!$OMP SCHEDULE(DYNAMIC)
          do ibox = laddr(1,ilev),laddr(2,ilev)
            nchild=itree(ipointer(4)+ibox-1)
            if(nchild.eq.0) then 
              istart = isrcse(1,ibox) 
              iend = isrcse(2,ibox)
              npts = iend-istart+1
              call l3dtaevalg(nd,rscales(ilev),centers(1,ibox), &
     &    rmlexp(iaddr(2,ibox)),nterms(ilev),sourcesort(1,istart), &
     &    npts,pot(1,istart),grad(1,1,istart),wlege,nlege)
            endif
          enddo
!$OMP END PARALLEL DO
        endif


        if(ifpgh.eq.3) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,nchild,istart,iend,i,npts) &
!$OMP SCHEDULE(DYNAMIC)
          do ibox = laddr(1,ilev),laddr(2,ilev)
            nchild=itree(ipointer(4)+ibox-1)
            if(nchild.eq.0) then 
              istart = isrcse(1,ibox)
              iend = isrcse(2,ibox)
              npts = iend-istart+1
              call l3dtaevalh(nd,rscales(ilev),centers(1,ibox), &
     &    rmlexp(iaddr(2,ibox)),nterms(ilev),sourcesort(1,istart), &
     &    npts,pot(1,istart),grad(1,1,istart),hess(1,1,istart), &
     &    scarray(1,ilev))
            endif
          enddo
!$OMP END PARALLEL DO
        endif

        if(ifpghtarg.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,nchild,istart,iend,i,npts) &
!$OMP SCHEDULE(DYNAMIC)
          do ibox = laddr(1,ilev),laddr(2,ilev)
            nchild=itree(ipointer(4)+ibox-1)
            if(nchild.eq.0) then 
              istart = itargse(1,ibox)
              iend = itargse(2,ibox)
              npts = iend-istart+1
              call l3dtaevalp(nd,rscales(ilev),centers(1,ibox), &
     &    rmlexp(iaddr(2,ibox)),nterms(ilev),targsort(1,istart), &
     &    npts,pottarg(1,istart),wlege,nlege)
            endif
          enddo
!$OMP END PARALLEL DO
        endif

        if(ifpghtarg.eq.2) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,nchild,istart,iend,i,npts) &
!$OMP SCHEDULE(DYNAMIC)
          do ibox = laddr(1,ilev),laddr(2,ilev)
            nchild=itree(ipointer(4)+ibox-1)
            if(nchild.eq.0) then 
              istart = itargse(1,ibox)
              iend = itargse(2,ibox)
              npts = iend-istart+1

              call l3dtaevalg(nd,rscales(ilev),centers(1,ibox), &
     &    rmlexp(iaddr(2,ibox)),nterms(ilev),targsort(1,istart), &
     &    npts,pottarg(1,istart),gradtarg(1,1,istart),wlege,nlege)
            endif
          enddo
!$OMP END PARALLEL DO
        endif

        if(ifpghtarg.eq.3) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,nchild,istart,iend,i,npts) &
!$OMP SCHEDULE(DYNAMIC)
          do ibox = laddr(1,ilev),laddr(2,ilev)
            nchild=itree(ipointer(4)+ibox-1)
            if(nchild.eq.0) then 
              istart = itargse(1,ibox)
              iend = itargse(2,ibox)
              npts = iend-istart+1

              call l3dtaevalh(nd,rscales(ilev),centers(1,ibox), &
     &    rmlexp(iaddr(2,ibox)),nterms(ilev),targsort(1,istart), &
     &    npts,pottarg(1,istart),gradtarg(1,1,istart), &
     &    hesstarg(1,1,istart),scarray(1,ilev))
            endif
          enddo
!$OMP END PARALLEL DO
        endif
      enddo

    
      ! call cpu_time(time2)
!$        time2=omp_get_wtime()
      ! timeinfo(5) = time2 - time1


    !   if(ifprint .ge. 1) &
    !  &     call prinf('=== STEP 6 (direct) =====*',i,0)
    !   call cpu_time(time1)
!$        time1=omp_get_wtime()

      if(ifnear.eq.0) goto 1000
!
!c       directly form local expansions for list1 sources
!        at expansion centers. 
!        (note: this part is not relevant for particle codes.
!         It is relevant only for qbx codes)


      do ilev=0,nlevels
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,istarte,iende,i,jbox) &
!$OMP PRIVATE(jstart,jend)
         do ibox = laddr(1,ilev),laddr(2,ilev)
            istarte = iexpcse(1,ibox) 
            iende = iexpcse(2,ibox) 

            

            do i =1,nlist1(ibox)
               jbox = list1(i,ibox)
               jstart = isrcse(1,jbox) 
               jend = isrcse(2,jbox) 

               call lfmm3dexpc_direct_tree(nd,jstart,jend,istarte, &
     &    iende,sourcesort,ifcharge,chargesort,ifdipole, &
     &    dipvecsort,expcsort,tsort,scjsort,ntj, &
     &    wlege,nlege)
            enddo
         enddo
!$OMP END PARALLEL DO
      enddo

!
!c        directly evaluate potential at sources and targets 
!         due to sources in list1

      do ilev=0,nlevels
!
!c           evaluate at the sources
!

        if(ifpgh.eq.1) then
          if(ifcharge.eq.1.and.ifdipole.eq.0) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,istarts,iends,npts0,i,jbox,jstart,jend,npts) &
!$OMP SCHEDULE(DYNAMIC)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istarts = isrcse(1,ibox) 
              iends = isrcse(2,ibox)
              npts0 = iends-istarts+1

              do i=1,nlist1(ibox)
                jbox = list1(i,ibox) 
                jstart =  isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectcp(nd,sourcesort(1,jstart), &
     &    chargesort(1,jstart),npts,sourcesort(1,istarts), &
     &    npts0,pot(1,istarts),thresh)          
              enddo
            enddo
!$OMP END PARALLEL DO
          endif

          if(ifcharge.eq.0.and.ifdipole.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,istarts,iends,npts0,i,jbox,jstart,jend,npts) &
!$OMP SCHEDULE(DYNAMIC)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istarts = isrcse(1,ibox)
              iends = isrcse(2,ibox)
              npts0 = iends-istarts+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart =  isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectdp(nd,sourcesort(1,jstart), &
     &    dipvecsort(1,1,jstart),npts,sourcesort(1,istarts), &
     &    npts0,pot(1,istarts),thresh)          
              enddo
            enddo
!$OMP END PARALLEL DO
          endif

          if(ifcharge.eq.1.and.ifdipole.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,istarts,iends,npts0,i,jbox,jstart,jend,npts) &
!$OMP SCHEDULE(DYNAMIC)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istarts = isrcse(1,ibox)
              iends = isrcse(2,ibox)
              npts0 = iends-istarts+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectcdp(nd,sourcesort(1,jstart), &
     &    chargesort(1,jstart), &
     &    dipvecsort(1,1,jstart),npts,sourcesort(1,istarts), &
     &    npts0,pot(1,istarts),thresh)          
              enddo
            enddo
!$OMP END PARALLEL DO
          endif
        endif

        if(ifpgh.eq.2) then
          if(ifcharge.eq.1.and.ifdipole.eq.0) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,istarts,iends,npts0,i,jbox,jstart,jend,npts) &
!$OMP SCHEDULE(DYNAMIC)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istarts = isrcse(1,ibox)
              iends = isrcse(2,ibox)
              npts0 = iends-istarts+1
              

              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectcg(nd,sourcesort(1,jstart), &
     &    chargesort(1,jstart),npts,sourcesort(1,istarts), &
     &    npts0,pot(1,istarts),grad(1,1,istarts),thresh)   
              enddo
            enddo
!$OMP END PARALLEL DO
          endif

          if(ifcharge.eq.0.and.ifdipole.eq.1) then

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,istarts,iends,npts0,i,jbox,jstart,jend,npts) &
!$OMP SCHEDULE(DYNAMIC)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istarts = isrcse(1,ibox)
              iends = isrcse(2,ibox)
              npts0 = iends-istarts+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectdg(nd,sourcesort(1,jstart), &
     &    dipvecsort(1,1,jstart),npts,sourcesort(1,istarts), &
     &    npts0,pot(1,istarts),grad(1,1,istarts),thresh)     
              enddo
            enddo
!$OMP END PARALLEL DO
          endif

          if(ifcharge.eq.1.and.ifdipole.eq.1) then

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,istarts,iends,npts0,i,jbox,jstart,jend,npts) &
!$OMP SCHEDULE(DYNAMIC)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istarts = isrcse(1,ibox)
              iends = isrcse(2,ibox)
              npts0 = iends-istarts+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectcdg(nd,sourcesort(1,jstart), &
     &    chargesort(1,jstart), &
     &    dipvecsort(1,1,jstart),npts,sourcesort(1,istarts), &
     &    npts0,pot(1,istarts),grad(1,1,istarts),thresh)      
              enddo
            enddo
!$OMP END PARALLEL DO
          endif
        endif


        if(ifpgh.eq.3) then
          if(ifcharge.eq.1.and.ifdipole.eq.0) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,istarts,iends,npts0,i,jbox,jstart,jend,npts) &
!$OMP SCHEDULE(DYNAMIC)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istarts = isrcse(1,ibox)
              iends = isrcse(2,ibox)
              npts0 = iends-istarts+1
              

              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectch(nd,sourcesort(1,jstart), &
     &    chargesort(1,jstart),npts,sourcesort(1,istarts), &
     &    npts0,pot(1,istarts),grad(1,1,istarts), &
     &    hess(1,1,istarts),thresh)   
              enddo
            enddo
!$OMP END PARALLEL DO
          endif

          if(ifcharge.eq.0.and.ifdipole.eq.1) then

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,istarts,iends,npts0,i,jbox,jstart,jend,npts) &
!$OMP SCHEDULE(DYNAMIC)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istarts = isrcse(1,ibox)
              iends = isrcse(2,ibox)
              npts0 = iends-istarts+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectdh(nd,sourcesort(1,jstart), &
     &    dipvecsort(1,1,jstart),npts,sourcesort(1,istarts), &
     &    npts0,pot(1,istarts),grad(1,1,istarts), &
     &    hess(1,1,istarts),thresh)          
              enddo
            enddo
!$OMP END PARALLEL DO
          endif

          if(ifcharge.eq.1.and.ifdipole.eq.1) then

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,istarts,iends,npts0,i,jbox,jstart,jend,npts) &
!$OMP SCHEDULE(DYNAMIC)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istarts = isrcse(1,ibox)
              iends = isrcse(2,ibox)
              npts0 = iends-istarts+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectcdh(nd,sourcesort(1,jstart), &
     &    chargesort(1,jstart), &
     &    dipvecsort(1,1,jstart),npts,sourcesort(1,istarts), &
     &    npts0,pot(1,istarts),grad(1,1,istarts), &
     &    hess(1,1,istarts),thresh)          
              enddo
            enddo
!$OMP END PARALLEL DO
          endif
        endif

        if(ifpghtarg.eq.1) then
          if(ifcharge.eq.1.and.ifdipole.eq.0) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,istartt,iendt,npts0,i,jbox,jstart,jend,npts) &
!$OMP SCHEDULE(DYNAMIC)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istartt = itargse(1,ibox)
              iendt = itargse(2,ibox)
              npts0 = iendt-istartt+1
              

              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectcp(nd,sourcesort(1,jstart), &
     &    chargesort(1,jstart),npts,targsort(1,istartt), &
     &    npts0,pottarg(1,istartt),thresh)          
              enddo
            enddo
!$OMP END PARALLEL DO
          endif

          if(ifcharge.eq.0.and.ifdipole.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,istartt,iendt,npts0,i,jbox,jstart,jend,npts) &
!$OMP SCHEDULE(DYNAMIC)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istartt = itargse(1,ibox)
              iendt = itargse(2,ibox)
              npts0 = iendt-istartt+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectdp(nd,sourcesort(1,jstart), &
     &    dipvecsort(1,1,jstart),npts,targsort(1,istartt), &
     &    npts0,pottarg(1,istartt),thresh)          
              enddo
            enddo
!$OMP END PARALLEL DO
          endif

          if(ifcharge.eq.1.and.ifdipole.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,istartt,iendt,npts0,i,jbox,jstart,jend,npts) &
!$OMP SCHEDULE(DYNAMIC)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istartt = itargse(1,ibox)
              iendt = itargse(2,ibox)
              npts0 = iendt-istartt+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectcdp(nd,sourcesort(1,jstart), &
     &    chargesort(1,jstart), &
     &    dipvecsort(1,1,jstart),npts,targsort(1,istartt), &
     &    npts0,pottarg(1,istartt),thresh)          
              enddo
            enddo
!$OMP END PARALLEL DO
          endif
        endif

        if(ifpghtarg.eq.2) then
          if(ifcharge.eq.1.and.ifdipole.eq.0) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,istartt,iendt,npts0,i,jbox,jstart,jend,npts) &
!$OMP SCHEDULE(DYNAMIC)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istartt = itargse(1,ibox)
              iendt = itargse(2,ibox)
              npts0 = iendt-istartt+1
              

              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectcg(nd,sourcesort(1,jstart), &
     &    chargesort(1,jstart),npts,targsort(1,istartt), &
     &    npts0,pottarg(1,istartt),gradtarg(1,1,istartt), &
     &    thresh)   
              enddo
            enddo
!$OMP END PARALLEL DO
          endif

          if(ifcharge.eq.0.and.ifdipole.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,istartt,iendt,npts0,i,jbox,jstart,jend,npts) &
!$OMP SCHEDULE(DYNAMIC)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istartt = itargse(1,ibox)
              iendt = itargse(2,ibox)
              npts0 = iendt-istartt+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectdg(nd,sourcesort(1,jstart), &
     &    dipvecsort(1,1,jstart),npts,targsort(1,istartt), &
     &    npts0,pottarg(1,istartt),gradtarg(1,1,istartt), &
     &    thresh)          
              enddo
            enddo
!$OMP END PARALLEL DO
          endif

          if(ifcharge.eq.1.and.ifdipole.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,istartt,iendt,npts0,i,jbox,jstart,jend,npts) &
!$OMP SCHEDULE(DYNAMIC)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istartt = itargse(1,ibox)
              iendt = itargse(2,ibox)
              npts0 = iendt-istartt+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectcdg(nd,sourcesort(1,jstart), &
     &    chargesort(1,jstart), &
     &    dipvecsort(1,1,jstart),npts,targsort(1,istartt), &
     &    npts0,pottarg(1,istartt),gradtarg(1,1,istartt), &
     &    thresh)          
              enddo
            enddo
!$OMP END PARALLEL DO
          endif
        endif

        if(ifpghtarg.eq.3) then
          if(ifcharge.eq.1.and.ifdipole.eq.0) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,istartt,iendt,npts0,i,jbox,jstart,jend,npts) &
!$OMP SCHEDULE(DYNAMIC)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istartt = itargse(1,ibox)
              iendt = itargse(2,ibox)
              npts0 = iendt-istartt+1
              

              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectch(nd,sourcesort(1,jstart), &
     &    chargesort(1,jstart),npts,targsort(1,istartt), &
     &    npts0,pottarg(1,istartt),gradtarg(1,1,istartt), &
     &    hesstarg(1,1,istartt),thresh)   
              enddo
            enddo
!$OMP END PARALLEL DO
          endif

          if(ifcharge.eq.0.and.ifdipole.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,istartt,iendt,npts0,i,jbox,jstart,jend,npts) &
!$OMP SCHEDULE(DYNAMIC)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istartt = itargse(1,ibox)
              iendt = itargse(2,ibox)
              npts0 = iendt-istartt+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectdh(nd,sourcesort(1,jstart), &
     &    dipvecsort(1,1,jstart),npts,targsort(1,istartt), &
     &    npts0,pottarg(1,istartt),gradtarg(1,1,istartt), &
     &    hesstarg(1,1,istartt),thresh)          
              enddo
            enddo
!$OMP END PARALLEL DO
          endif

          if(ifcharge.eq.1.and.ifdipole.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,istartt,iendt,npts0,i,jbox,jstart,jend,npts) &
!$OMP SCHEDULE(DYNAMIC)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istartt = itargse(1,ibox)
              iendt = itargse(2,ibox)
              npts0 = iendt-istartt+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectcdh(nd,sourcesort(1,jstart), &
     &    chargesort(1,jstart), &
     &    dipvecsort(1,1,jstart),npts,targsort(1,istartt), &
     &    npts0,pottarg(1,istartt),gradtarg(1,1,istartt), &
     &    hesstarg(1,1,istartt),thresh)          
              enddo
            enddo
!$OMP END PARALLEL DO
          endif
        endif
      enddo
 1000 continue      
 
      call cpu_time(time2)
!$        time2=omp_get_wtime()
      timeinfo(6) = time2-time1
      if(ifprint.ge.1) call prin2('timeinfo=*',timeinfo,6)
      d = 0
      do i = 1,6
         d = d + timeinfo(i)
      enddo

      if(ifprint.ge.1) call prin2('sum(timeinfo)=*',d,1)

      return
      end
!------------------------------------------------
      subroutine lfmm3dexpc_direct_tree(nd,istart,iend,jstart,jend, &
     &     source,ifcharge,charge,ifdipole, &
     &     dipvec,expc,texps,scj,ntj,wlege,nlege)
!--------------------------------------------------------------------
!     This subroutine adds the local expansions due to sources
!     istart to iend in the source array at the expansion centers
!     jstart to jend in the expansion center array to the existing 
!     local expansions at the corresponding expansion centers.
!
!     INPUT arguments
!-------------------------------------------------------------------
!     nd           in: integer
!                  number of charge densities
!
!     istart       in:Integer
!                  Starting index in source array whose expansions
!                  we wish to add
!
!     iend         in:Integer
!                  Last index in source array whose expansions
!                  we wish to add
!
!     jstart       in: Integer
!                  First index in the expansion center array at 
!                  which we  wish to compute the expansions
! 
!     jend         in:Integer
!                  Last index in expansion center array at 
!                  which we wish to compute the expansions
! 
!     scjsort      in: double precision(*)
!                  Scale of expansions formed at the expansion centers
!
!     source       in: double precision(3,ns)
!                  Source locations
!
!     ifcharge     in: Integer
!                  flag for including expansions due to charges
!                  The expansion due to charges will be included
!                  if ifcharge == 1
!
!     charge       in: double precision
!                  Charge at the source locations
!
!     ifdipole     in: Integer
!                 flag for including expansions due to dipoles
!                 The expansion due to dipoles will be included
!                 if ifdipole == 1
!
!     dipvec      in: double precision(3,ns)
!                 Dipole orientation vector at the source locations
!
!     expc        in: double precision(3,nexpc)
!                 Expansion center locations
!
!     ntj         in: Integer
!                 Number of terms in expansion
!
!     wlege       in: double precision(0:nlege,0:nlege)
!                 precomputed array of recurrence relation
!                 coeffs for Ynm calculation.
!
!    nlege        in: integer
!                 dimension parameter for wlege
!------------------------------------------------------------
!     OUTPUT
!
!   Updated expansions at the targs
!   texps       out: double complex(0:ntj,-ntj:ntj,expc) 
!                 coeffs for local expansions
!-------------------------------------------------------               
        implicit none
!
        integer istart,iend,jstart,jend,ns,j, nlege
        integer ifcharge,ifdipole,ier,nd
        double precision source(3,*)
        double precision scj(*)
        double precision wlege(*)
        double precision charge(nd,*)
        double precision dipvec(nd,3,*)
        double precision expc(3,*)

        integer nlevels,ntj
!
        double complex texps(nd,0:ntj,-ntj:ntj,*)
        
!
        ns = iend - istart + 1
        if(ifcharge.eq.1.and.ifdipole.eq.0) then
          do j=jstart,jend
            call l3dformtac(nd,scj(j), &
     &    source(1,istart),charge(1,istart),ns, &
     &    expc(1,j),ntj,texps(1,0,-ntj,j),wlege,nlege)
           enddo
         endif

         if(ifcharge.eq.0.and.ifdipole.eq.1) then
          do j=jstart,jend
            call l3dformtad(nd,scj(j), &
     &    source(1,istart), &
     &    dipvec(1,1,istart),ns,expc(1,j),ntj,texps(1,0,-ntj,j), &
     &    wlege,nlege)
           enddo
         endif

         if(ifcharge.eq.1.and.ifdipole.eq.1) then
          do j=jstart,jend
            call l3dformtacd(nd,scj(j), &
     &    source(1,istart),charge(1,istart), &
     &    dipvec(1,1,istart),ns,expc(1,j),ntj,texps(1,0,-ntj,j), &
     &    wlege,nlege)
           enddo
         endif

!
        return
        end


end module fmm3d_tree_mod
!
       subroutine lfmm3d_tree(nd,eps,nsource,source,ifcharge, &
     &    charge,ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg, &
     &    targ,ifpghtarg,pottarg,gradtarg,hesstarg,ier)
!
!        Laplace FMM in R^{3}: evaluate all pairwise particle
!        interactions (ignoring self-interactions) and interactions
!        with targs.
!
!        We use (1/r) for the Green's function, without the
!        1/(4 \pi) scaling.
!
!
!        Input parameters:
!
!   nd:   number of densities
!
!   eps:  requested precision
!
!   nsource in: integer  
!                number of sources
!
!   source  in: double precision (3,nsource)
!                source(k,j) is the kth component of the jth
!                source locations
!
!   ifcharge  in: integer  
!             charge computation flag
!              ifcharge = 1   =>  include charge contribution
!                                     otherwise do not
! 
!   charge    in: double precision (nsource) 
!              charge strengths
!
!   ifdipole   in: integer
!              dipole computation flag
!              ifdipole = 1   =>  include dipole contribution
!                                     otherwise do not
!
!
!   dipvec   in: double precision (3,nsource) 
!              dipole orientation vectors
!   iper    in: integer
!             flag for periodic implmentations. Currently unused
!
!   ifpgh   in: integer
!              flag for evaluating potential/gradient at the sources
!              ifpgh = 1, only potential is evaluated
!              ifpgh = 2, potential and gradients are evaluated
!              ifpgh = 3, potential, gradients, and hessian are
!                 evaluated
!
!   ntarg  in: integer  
!                 number of targs 
!
!   targ  in: double precision (3,ntarg)
!               targ(k,j) is the kth component of the jth
!               targ location
!
!   ifpghtarg   in: integer
!              flag for evaluating potential/gradient at the targs
!              ifpghtarg = 1, only potential is evaluated
!              ifpghtarg = 2, potential and gradient are evaluated
!              ifpghtarg = 3, potential, gradient, and hess are
!                  evaluated
!
!
!     OUTPUT parameters:
!
!
!   pot:    out: double precision(nd,nsource) 
!               potential at the source locations
!
!   grad:   out: double precision(nd,3,nsource)
!               gradient at the source locations
!
!   hess    out: double precision(nd,6,nsource)
!               hessian at the source locations
!
!   pottarg:    out: double precision(nd,ntarg) 
!               potential at the targ locations
!
!   gradtarg:   out: double precision(nd,3,ntarg)
!               gradient at the targ locations
!
!   hesstarg    out: double precision(nd,6,ntarg)
!                hessian at the target locations 
!   ier         out: integer
!                error flag
!                ier = 0, for successful execution
!                ier = 4, if failed to allocate workspace
!                      for multipole and local expansions
!                ier = 8, if failed to allocate workspace
!                      for plane wave expansions
!     
!     
       implicit none

       integer nd,ier,iper

       double precision eps

       integer ifcharge,ifdipole
       integer ifpgh,ifpghtarg

       integer ntarg,nsource
       

       double precision source(3,*),targ(3,*)
       double precision charge(nd,*)
       double precision dipvec(nd,3,*)

       double precision pot(nd,*),grad(nd,3,*),hess(nd,6,*)
       double precision pottarg(nd,*),gradtarg(nd,3,*),hesstarg(nd,6,*)

       double precision timeinfo(6)

!
!c       tree variables
!
       integer idivflag,ndiv,nboxes,nlevels
       integer nlmax,nlmin,ifunif
       integer *8 ipointer(8),ltree
       integer, allocatable :: itree(:)
       integer, allocatable :: isrcse(:,:),itargse(:,:),isrc(:)
       integer, allocatable :: itarg(:)
       integer, allocatable :: iexpcse(:,:)
       integer iexpc
       double precision, allocatable :: treecenters(:,:),boxsize(:)
       double precision b0,b0inv,b0inv2,b0inv3

!
!c       temporary sorted arrays
!
       double precision, allocatable :: sourcesort(:,:),targsort(:,:)
       double precision, allocatable :: radsrc(:)
       double precision, allocatable :: chargesort(:,:)
       double precision, allocatable :: dipvecsort(:,:,:)

       double precision, allocatable :: potsort(:,:),gradsort(:,:,:)
       double precision, allocatable :: hesssort(:,:,:)
       double precision, allocatable :: pottargsort(:,:)
       double precision, allocatable :: gradtargsort(:,:,:)
       double precision, allocatable :: hesstargsort(:,:,:)
!
!c        temporary fmm arrays
!
       integer, allocatable :: nterms(:)
       integer *8, allocatable :: iaddr(:,:)
       double precision, allocatable :: scales(:)
       double precision, allocatable :: rmlexp(:)

       integer lmptemp,nmax
       integer *8 lmptot
       double precision, allocatable :: mptemp(:),mptemp2(:)

!
!c       temporary variables not main fmm routine but
!        not used in particle code
       double precision expc(3),scjsort(1),radexp
       double complex texpssort(100)
       double precision expcsort(3)
       integer ntj,nexpc,nadd,ifnear

!
!c         other temporary variables
!
        integer i,iert,ifprint,ilev,idim
        double precision time1,time2,omp_get_wtime,second

!
!
!     ifprint is an internal information printing flag. 
!     Suppressed if ifprint=0.
!     Prints timing breakdown and other things if ifprint=1.
!      

      call cpu_time(time1)
!$     time1=omp_get_wtime()      
      ifprint=0

!
!c        figure out tree structure
!
!
!c         set criterion for box subdivision
!
      call lndiv(eps,nsource,ntarg,ifcharge,ifdipole,ifpgh, &
     &    ifpghtarg,ndiv,idivflag) 


!
!       turn on computation of list 1
!
      ifnear = 1


       nexpc = 0
       nadd = 0
       ntj = 0


!
!c      set tree flags
! 
       nlmax = 51
       nlevels = 0
       nboxes = 0
       ltree = 0
       nlmin = 0
       iper = 0
       ifunif = 0

!
!c     memory management code for contructing level restricted tree
      call pts_tree_mem(source,nsource,targ,ntarg,idivflag,ndiv,nlmin, &
     &    nlmax,iper,ifunif,nlevels,nboxes,ltree)
      

        if(ifprint.ge.1) print *, ltree/1.0d9

        allocate(itree(ltree))
        allocate(boxsize(0:nlevels))
        allocate(treecenters(3,nboxes))

!       Call tree code
      call pts_tree_build(source,nsource,targ,ntarg,idivflag,ndiv, &
     &    nlmin,nlmax,iper,ifunif,nlevels,nboxes,ltree,itree,ipointer, &
     &    treecenters,boxsize)
      

      allocate(isrcse(2,nboxes),itargse(2,nboxes),iexpcse(2,nboxes))
      allocate(isrc(nsource),itarg(ntarg))

      call pts_tree_sort(nsource,source,itree,ltree,nboxes,nlevels, &
     &    ipointer,treecenters,isrc,isrcse)
      
      call pts_tree_sort(ntarg,targ,itree,ltree,nboxes,nlevels, &
     &    ipointer,treecenters,itarg,itargse)
      
      call pts_tree_sort(nexpc,expc,itree,ltree,nboxes,nlevels, &
     &    ipointer,treecenters,iexpc,iexpcse)

!
!   End of tree build
!

!
!  Set rescaling parameters
!
      b0 = boxsize(0)
      b0inv = 1.0d0/b0
      b0inv2 = b0inv**2
      b0inv3 = b0inv2*b0inv

!     Allocate sorted source and targ arrays      

      allocate(sourcesort(3,nsource))
      allocate(targsort(3,ntarg))
      if(ifcharge.eq.1) allocate(chargesort(nd,nsource))

      if(ifdipole.eq.1) then
         allocate(dipvecsort(nd,3,nsource))
      endif

      if(ifpgh.eq.1) then 
        allocate(potsort(nd,nsource),gradsort(nd,3,1),hesssort(nd,6,1))
      else if(ifpgh.eq.2) then
        allocate(potsort(nd,nsource),gradsort(nd,3,nsource), &
     &    hesssort(nd,6,1))
      else if(ifpgh.eq.3) then
        allocate(potsort(nd,nsource),gradsort(nd,3,nsource), &
     &    hesssort(nd,6,nsource))
      else
        allocate(potsort(nd,1),gradsort(nd,3,1),hesssort(nd,6,1))
      endif

      if(ifpghtarg.eq.1) then
        allocate(pottargsort(nd,ntarg),gradtargsort(nd,3,1), &
     &    hesstargsort(nd,6,1))
      else if(ifpghtarg.eq.2) then
        allocate(pottargsort(nd,ntarg),gradtargsort(nd,3,ntarg), &
     &    hesstargsort(nd,6,1))
      else if(ifpghtarg.eq.3) then
        allocate(pottargsort(nd,ntarg),gradtargsort(nd,3,ntarg), &
     &    hesstargsort(nd,6,ntarg))
      else
        allocate(pottargsort(nd,1),gradtargsort(nd,3,1), &
     &    hesstargsort(nd,6,1))
      endif


!
!c      initialize potential and gradient at source
!       locations
!
      if(ifpgh.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)
        do i=1,nsource
          do idim=1,nd
            potsort(idim,i) = 0
          enddo
        enddo
!$OMP END PARALLEL DO
      endif

      if(ifpgh.eq.2) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)

        do i=1,nsource
          do idim=1,nd
            potsort(idim,i) = 0
            gradsort(idim,1,i) = 0
            gradsort(idim,2,i) = 0
            gradsort(idim,3,i) = 0
          enddo
        enddo
!$OMP END PARALLEL DO
      endif


      if(ifpgh.eq.3) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)
        do i=1,nsource
          do idim=1,nd
            potsort(idim,i) = 0
            gradsort(idim,1,i) = 0
            gradsort(idim,2,i) = 0
            gradsort(idim,3,i) = 0
            hesssort(idim,1,i) = 0
            hesssort(idim,2,i) = 0
            hesssort(idim,3,i) = 0
            hesssort(idim,4,i) = 0
            hesssort(idim,5,i) = 0
            hesssort(idim,6,i) = 0
          enddo
        enddo
!$OMP END PARALLEL DO
      endif



!
!c       initialize potential and gradient  at targ
!        locations
!
      if(ifpghtarg.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)
        do i=1,ntarg
          do idim=1,nd
            pottargsort(idim,i) = 0
          enddo
        enddo
!$OMP END PARALLEL DO
      endif

      if(ifpghtarg.eq.2) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)
        do i=1,ntarg
          do idim=1,nd
            pottargsort(idim,i) = 0
            gradtargsort(idim,1,i) = 0
            gradtargsort(idim,2,i) = 0
            gradtargsort(idim,3,i) = 0
          enddo
        enddo
!$OMP END PARALLEL DO
      endif

      if(ifpghtarg.eq.3) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)
        do i=1,ntarg
          do idim=1,nd
            pottargsort(idim,i) = 0
            gradtargsort(idim,1,i) = 0
            gradtargsort(idim,2,i) = 0
            gradtargsort(idim,3,i) = 0
            hesstargsort(idim,1,i) = 0
            hesstargsort(idim,2,i) = 0
            hesstargsort(idim,3,i) = 0
            hesstargsort(idim,4,i) = 0
            hesstargsort(idim,5,i) = 0
            hesstargsort(idim,6,i) = 0
          enddo
        enddo
!$OMP END PARALLEL DO
      endif

      allocate(nterms(0:nlevels))

!     Compute length of expansions at each level      
      nmax = 0
      do i=0,nlevels
         call l3dterms(eps,nterms(i))
         if(nterms(i).gt.nmax) nmax = nterms(i)
      enddo
!       
!     Multipole and local expansions will be held in workspace
!     in locations pointed to by array iaddr(2,nboxes).
!
!     iaddr is pointer to iaddr array, itself contained in workspace.
!     imptemp is pointer for single expansion (dimensioned by nmax)
!
!       ... allocate iaddr and temporary arrays
!

      allocate(iaddr(2,nboxes))
      lmptemp = (nmax+1)*(2*nmax+1)*2*nd
      allocate(mptemp(lmptemp),mptemp2(lmptemp))

!
!c     reorder sources 
!
      call dreorderf(3,nsource,source,sourcesort,isrc)

!
!       rescale sources to be contained in unit box
!
      call drescale(3*nsource,sourcesort,b0inv)

      if(ifcharge.eq.1) then
        call dreorderf(nd,nsource,charge,chargesort, &
     &    isrc)
        call drescale(nd*nsource,chargesort,b0inv)
      endif


      if(ifdipole.eq.1) then
         call dreorderf(3*nd,nsource,dipvec,dipvecsort, &
     &    isrc)
         call drescale(3*nd*nsource,dipvecsort,b0inv2)
      endif

!
!c      reorder and rescale targs
!
      call dreorderf(3,ntarg,targ,targsort,itarg)
      call drescale(3*ntarg,targsort,b0inv)


!
!        update tree centers and boxsize
!
      call drescale(3*nboxes,treecenters,b0inv)
      call drescale(nlevels+1,boxsize,b0inv)

!
!     allocate memory need by multipole, local expansions at all
!     levels
!     irmlexp is pointer for workspace need by various fmm routines,
!
      call mpalloc(nd,itree(ipointer(1)),iaddr,nlevels,lmptot,nterms)
      if(ifprint.ge. 1) print *, "lmptot =",lmptot/1.0d9


      allocate(rmlexp(lmptot),stat=iert)
      if(iert.ne.0) then
         print *, "Cannot allocate mpole expansion workspace"
         print *, "lmptot=", lmptot
         ier = 4
         return
      endif

!     Memory allocation is complete. 
!     scaling factor for multipole and local expansions at all levels
!
      allocate(scales(0:nlevels))
      do ilev = 0,nlevels
        scales(ilev) = boxsize(ilev)
      enddo

      call cpu_time(time2)
!$     time2=omp_get_wtime()      

      if(ifprint.ge.1)  &
     &    call prin2('time before fmm main=*',time2-time1,1)
!     Call main fmm routine

      call cpu_time(time1)
!$      time1=omp_get_wtime()
    !   call lfmm3dmain_tree(nd,eps, &
    !  &   nsource,sourcesort, &
    !  &   ifcharge,chargesort, &
    !  &   ifdipole,dipvecsort, &
    !  &   ntarg,targsort,nexpc,expcsort, &
    !  &   iaddr,rmlexp,lmptot,mptemp,mptemp2,lmptemp, &
    !  &   itree,ltree,ipointer,ndiv,nlevels, &
    !  &   nboxes,iper,boxsize,treecenters,isrcse,itargse,iexpcse, &
    !  &   scales,itree(ipointer(1)),nterms, &
    !  &   ifpgh,potsort,gradsort,hesssort, &
    !  &   ifpghtarg,pottargsort,gradtargsort,hesstargsort,ntj, &
    !  &   texpssort,scjsort,ifnear,timeinfo,ier)
    !   if(ier.ne.0) return

      call cpu_time(time2)
!$        time2=omp_get_wtime()
      if( ifprint .eq. 1 ) call prin2('time in fmm main=*', &
     &    time2-time1,1)



      if(ifpgh.ge.1) then
        call dreorderi(nd,nsource,potsort,pot,isrc)
      endif
      if(ifpgh.ge.2) then 
        call dreorderi(3*nd,nsource,gradsort,grad,isrc)
        call drescale(nd*3*nsource,grad,b0inv)
      endif

      if(ifpgh.ge.3) then 
        call dreorderi(6*nd,nsource,hesssort,hess,isrc)
        call drescale(nd*6*nsource,hess,b0inv2)
      endif


      if(ifpghtarg.ge.1) then
        call dreorderi(nd,ntarg,pottargsort,pottarg,itarg)
      endif

      if(ifpghtarg.ge.2) then
        call dreorderi(3*nd,ntarg,gradtargsort,gradtarg,itarg)
        call drescale(nd*3*ntarg,gradtarg,b0inv)
      endif

      if(ifpghtarg.ge.3) then
        call dreorderi(6*nd,ntarg,hesstargsort,hesstarg,itarg)
        call drescale(nd*6*ntarg,hesstarg,b0inv2)
      endif

      return
      end

!       
!---------------------------------------------------------------
!
      subroutine lfmm3dmain_tree_old(nd,eps, &
     &     nsource,sourcesort, &
     &     ifcharge,chargesort, &
     &     ifdipole,dipvecsort, &
     &     ntarg,targsort,nexpc,expcsort, &
     &     iaddr,rmlexp,lmptot,mptemp,mptemp2,lmptemp, &
     &     itree,ltree,ipointer,ndiv,nlevels,  &
     &     nboxes,iper,boxsize,centers,isrcse,itargse,iexpcse, &
     &     rscales,laddr,nterms, &
     &     ifpgh,pot,grad,hess, &
     &     ifpghtarg,pottarg,gradtarg,hesstarg,ntj, &
     &     tsort,scjsort,ifnear,timeinfo,ier)
      implicit none
      integer nd
      integer ier
      double precision eps
      integer nsource,ntarg,nexpc
      integer ndiv,nlevels

      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg

      double precision sourcesort(3,nsource)

      double precision chargesort(nd,*)
      double precision dipvecsort(nd,3,*)

      double precision targsort(3,ntarg)

      double precision pot(nd,*),grad(nd,3,*),hess(nd,6,*)
      double precision pottarg(nd,*),gradtarg(nd,3,*),hesstarg(nd,6,*)

      integer ntj
      integer ifnear
      double precision expcsort(3,nexpc)
      double complex tsort(nd,0:ntj,-ntj:ntj,nexpc)
      double precision scjsort(nexpc)

      integer nboxes
      integer *8 iaddr(2,nboxes), lmptot
      integer lmptemp
      double precision rmlexp(lmptot)
      double precision mptemp(lmptemp)
      double precision mptemp2(lmptemp)

      double precision thresh
       
      double precision timeinfo(6)
      double precision centers(3,nboxes)

      integer isep,iper
      integer laddr(2,0:nlevels)
      integer nterms(0:nlevels)
      integer *8 ipointer(8),ltree
      integer itree(ltree)
      double precision rscales(0:nlevels)
      double precision boxsize(0:nlevels)
      integer isrcse(2,nboxes),itargse(2,nboxes),iexpcse(2,nboxes)
      integer, allocatable :: nlist1(:),list1(:,:)
      integer, allocatable :: nlist2(:),list2(:,:)
      integer, allocatable :: nlist3(:),list3(:,:)
      integer, allocatable :: nlist4(:),list4(:,:)

      integer nuall,ndall,nnall,nsall,neall,nwall
      integer nu1234,nd5678,nn1256,ns3478,ne1357,nw2468
      integer nn12,nn56,ns34,ns78,ne13,ne57,nw24,nw68
      integer ne1,ne3,ne5,ne7,nw2,nw4,nw6,nw8

      integer, allocatable :: uall(:,:),dall(:,:),nall(:,:)
      integer, allocatable :: sall(:,:),eall(:,:),wall(:,:)
      integer, allocatable :: u1234(:,:),d5678(:,:)
      integer, allocatable :: n1256(:,:),s3478(:,:)
      integer, allocatable :: e1357(:,:),w2468(:,:)
      integer, allocatable :: n12(:,:),n56(:,:),s34(:,:),s78(:,:)
      integer, allocatable :: e13(:,:),e57(:,:),w24(:,:),w68(:,:)
      integer, allocatable :: e1(:,:),e3(:,:),e5(:,:),e7(:,:)
      integer, allocatable :: w2(:,:),w4(:,:),w6(:,:),w8(:,:)

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
      double precision, allocatable :: carray(:,:), dc(:,:)
      double precision, allocatable :: cs(:,:),fact(:),rdplus(:,:,:)
      double precision, allocatable :: rdminus(:,:,:), rdsq3(:,:,:)
      double precision, allocatable :: rdmsq3(:,:,:)
  
      double precision, allocatable :: rlams(:),whts(:)

      double precision, allocatable :: rlsc(:,:,:)
      integer, allocatable :: nfourier(:), nphysical(:)
      integer nexptot, nexptotp
      double complex, allocatable :: xshift(:,:)
      double complex, allocatable :: yshift(:,:)
      double precision, allocatable :: zshift(:,:)

      double complex, allocatable :: fexpe(:),fexpo(:),fexpback(:)
      double complex, allocatable :: mexp(:,:,:,:)
      double complex, allocatable :: mexpf1(:,:,:),mexpf2(:,:,:)
      double complex, allocatable :: &
     &    mexpp1(:,:,:),mexpp2(:,:,:),mexppall(:,:,:,:)

      double complex, allocatable :: tmp(:,:,:,:)
      double precision, allocatable :: mptmp(:,:)

      double precision sourcetmp(3)
      double complex chargetmp

      integer ix,iy,iz,ictr
      double precision rtmp
      double complex zmul

      integer nlege, lw7, lused7, itype
      double precision wlege(40000)
      integer nterms_eval(4,0:nlevels)

      integer mnlist1, mnlist2,mnlist3,mnlist4,mnbors
      double complex eye, ztmp
      double precision alphaj
      integer ctr,nn,iptr1,iptr2
      double precision, allocatable :: rscpow(:)
      double precision pi,errtmp
      double complex ima

      double precision ctmp(3)

!     list 3 variables
      double complex, allocatable :: iboxlexp(:,:,:)
      double precision, allocatable :: iboxsubcenters(:,:,:)
      double precision, allocatable :: iboxpot(:,:,:)
      double precision, allocatable :: iboxgrad(:,:,:,:)
      double precision, allocatable :: iboxhess(:,:,:,:)
      double precision, allocatable :: iboxsrc(:,:,:)
      integer, allocatable :: iboxsrcind(:,:)
      integer, allocatable :: iboxfl(:,:,:)
!     end of list 3 variables
!     list 4 variables
      integer cntlist4
      integer, allocatable :: list4ct(:),ilist4(:)
      double complex, allocatable :: gboxmexp(:,:,:)
      double complex, allocatable :: gboxwexp(:,:,:,:,:)
      double complex, allocatable :: pgboxwexp(:,:,:,:)
      double precision, allocatable :: gboxsubcenters(:,:,:)
      double precision, allocatable :: gboxsort(:,:,:)
      integer, allocatable :: gboxind(:,:)
      integer, allocatable :: gboxfl(:,:,:)
      double precision, allocatable :: gboxcgsort(:,:,:)
      double precision, allocatable :: gboxdpsort(:,:,:,:)
!     end of list 4 variables

!
!   hessian variables
!
      double precision, allocatable :: scarray(:,:)

      integer *8 bigint
      integer iert
      data ima/(0.0d0,1.0d0)/

      integer nthd,ithd
      integer omp_get_max_threads,omp_get_thread_num
      nthd = 1
!$    nthd=omp_get_max_threads()

      pi = 4.0d0*atan(1.0d0)

      thresh = 2.0d0**(-51)*boxsize(0)


!     ifprint is an internal information printing flag. 
!     Suppressed if ifprint=0.
!     Prints timing breakdown and other things if ifprint=1.
!     Prints timing breakdown, list information, 
!     and other things if ifprint=2.
!       
      ifprint=0
      

      print *, " laddr = ", laddr
!
!   initialize various tree lists
!
      mnlist1 = 0
      mnlist2 = 0
      mnlist3 = 0
      mnlist4 = 0
      mnbors = 27

      isep = 1
      
      call computemnlists(nlevels,nboxes,itree(ipointer(1)),boxsize, &
     &    centers,itree(ipointer(3)),itree(ipointer(4)), &
     &    itree(ipointer(5)),isep,itree(ipointer(6)),mnbors, &
     &    itree(ipointer(7)),iper,mnlist1,mnlist2,mnlist3,mnlist4)
      
      allocate(list1(mnlist1,nboxes),nlist1(nboxes))
      allocate(list2(mnlist2,nboxes),nlist2(nboxes))
      allocate(list3(mnlist3,nboxes),nlist3(nboxes))
      allocate(list4(mnlist4,nboxes),nlist4(nboxes))

      call computelists(nlevels,nboxes,itree(ipointer(1)),boxsize, &
     &    centers,itree(ipointer(3)),itree(ipointer(4)), &
     &    itree(ipointer(5)),isep,itree(ipointer(6)),mnbors, &
     &    itree(ipointer(7)),iper,nlist1,mnlist1,list1,nlist2, &
     &    mnlist2,list2,nlist3,mnlist3,list3, &
     &    nlist4,mnlist4,list4)
      

!     Initialize routines for plane wave mp loc translation
 
      if(isep.eq.1) then
         if(eps.ge.0.5d-2) nlams = 12
         if(eps.lt.0.5d-2.and.eps.ge.0.5d-3) nlams = 12
         if(eps.lt.0.5d-3.and.eps.ge.0.5d-6) nlams = 20
         if(eps.lt.0.5d-6.and.eps.ge.0.5d-9) nlams = 29
         if(eps.lt.0.5d-9) nlams = 37
      endif
      if(isep.eq.2) then
         if(eps.ge.0.5d-3) nlams = 9
         if(eps.lt.0.5d-3.and.eps.ge.0.5d-6) nlams = 15
         if(eps.lt.0.5d-6.and.eps.ge.0.5d-9) nlams = 22
         if(eps.lt.0.5d-9) nlams = 29
      endif

      allocate(rlams(nlams),whts(nlams))
      allocate(nphysical(nlams),nfourier(nlams))

      nmax = 0
      do i=0,nlevels
         if(nmax.lt.nterms(i)) nmax = nterms(i)
      enddo
      allocate(rscpow(0:nmax))
      allocate(carray(4*nmax+1,4*nmax+1))
      allocate(dc(0:4*nmax,0:4*nmax))
      allocate(rdplus(0:nmax,0:nmax,-nmax:nmax))
      allocate(rdminus(0:nmax,0:nmax,-nmax:nmax))
      allocate(rdsq3(0:nmax,0:nmax,-nmax:nmax))
      allocate(rdmsq3(0:nmax,0:nmax,-nmax:nmax))
      allocate(rlsc(0:nmax,0:nmax,nlams))


!     generate rotation matrices and carray
      call getpwrotmat(nmax,carray,rdplus,rdminus,rdsq3,rdmsq3,dc)


!     generate rlams and weights (these are the nodes
!     and weights for the lambda integral)

      call vwts(rlams,whts,nlams)


!     generate the number of fourier modes required to represent the
!     moment function in fourier space

      call numthetahalf(nfourier,nlams)
 
!     generate the number of fourier modes in physical space
!     required for the exponential representation
      call numthetafour(nphysical,nlams)

!     Generate powers of lambda for the exponential basis
      call rlscini(rlsc,nlams,rlams,nmax)

!
!
!
      nn = 10*(nmax+2)**2
      allocate(scarray(nn,0:nlevels))
      do ilev=0,nlevels
        call l3dtaevalhessdini(nterms(ilev),scarray(1,ilev))
      enddo

!     Compute total number of plane waves
      nexptotp = 0
      nexptot = 0
      nthmax = 0
      nphmax = 0
      nn = 0
      do i=1,nlams
         nexptot = nexptot + nfourier(i)
         nexptotp = nexptotp + nphysical(i)
         if(nfourier(i).gt.nthmax) nthmax = nfourier(i)
         if(nphysical(i).gt.nphmax) nphmax = nphysical(i)
         nn = nn + nphysical(i)*nfourier(i)
      enddo

      allocate(fexpe(nn),fexpo(nn),fexpback(nn))
      allocate(tmp(nd,0:nmax,-nmax:nmax,nthd))
      allocate(mptmp(lmptemp,nthd))

      allocate(xshift(-5:5,nexptotp))
      allocate(yshift(-5:5,nexptotp))
      allocate(zshift(5,nexptotp))

      allocate(mexpf1(nd,nexptot,nthd),mexpf2(nd,nexptot,nthd), &
     &    mexpp1(nd,nexptotp,nthd))
      allocate(mexpp2(nd,nexptotp,nthd),mexppall(nd,nexptotp,16,nthd))

!
!c      NOTE: there can be some memory savings here
!
      bigint = 0
      bigint = nboxes
      bigint = bigint*6
      bigint = bigint*nexptotp*nd

      if(ifprint.ge.1) print *, "mexp memory=",bigint/1.0d9

      allocate(mexp(nd,nexptotp,nboxes,6),stat=iert)
      if(iert.ne.0) then
        print *, "Cannot allocate pw expansion workspace"
        print *, "bigint=", bigint
        ier = 8
        return
      endif

      allocate(list4ct(nboxes))
      allocate(ilist4(nboxes))
      do i=1,nboxes
        list4ct(i)=0
        ilist4(i)=0
      enddo
      cntlist4=0

!     Precompute table for shifting exponential coefficients in 
!     physical domain
      call mkexps(rlams,nlams,nphysical,nexptotp,xshift,yshift,zshift)

!     Precompute table of exponentials for mapping from
!     fourier to physical domain
      call mkfexp(nlams,nfourier,nphysical,fexpe,fexpo,fexpback)
      
!
!c    compute array of factorials

     
      nmax2 = 2*nmax
      allocate(fact(0:nmax2),cs(0:nmax,-nmax:nmax))
      
      d = 1
      fact(0) = d
      do i=1,nmax2
        d=d*sqrt(i+0.0d0)
        fact(i) = d
      enddo

      cs(0,0) = 1.0d0
      do l=1,nmax
        do m=0,l
          cs(l,m) = ((-1)**l)/(fact(l-m)*fact(l+m))
          cs(l,-m) = cs(l,m)
        enddo
      enddo


      
      if(ifprint.ge.1)  &
     &    call prin2('end of generating plane wave info*',i,0)
!
!
!     ... set the expansion coefficients to zero
!
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,idim)
      do i=1,nexpc
        do k=-ntj,ntj
          do j = 0,ntj
            do idim=1,nd
              tsort(idim,j,k,i)=0
            enddo
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

!       
      do i=1,6
        timeinfo(i)=0
      enddo

!
!       ... set all multipole and local expansions to zero
!

      do ilev = 0,nlevels
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox)
        do ibox=laddr(1,ilev),laddr(2,ilev)
          call mpzero(nd,rmlexp(iaddr(1,ibox)),nterms(ilev))
          call mpzero(nd,rmlexp(iaddr(2,ibox)),nterms(ilev))
        enddo
!$OMP END PARALLEL DO
      enddo

!
!      set scjsort
!
      do ilev=0,nlevels
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,nchild,istart,iend,i)
         do ibox=laddr(1,ilev),laddr(2,ilev)
            nchild = itree(ipointer(4)+ibox-1)
            if(nchild.gt.0) then
               istart = iexpcse(1,ibox)
               iend = iexpcse(2,ibox) 
               do i=istart,iend
                  scjsort(i) = rscales(ilev)
               enddo
            endif
         enddo
!$OMP END PARALLEL DO
      enddo


!    initialize legendre function evaluation routines
      nlege = 100
      lw7 = 40000
      call ylgndrfwini(nlege,wlege,lw7,lused7)

!
!     count number of boxes are in list4
      lca = 4*nmax
      if(ifprint.ge.1) &
     &   call prinf('=== STEP 0 list4===*',i,0)
      call cpu_time(time1)
!$    time1=omp_get_wtime()
      do ilev=1,nlevels-1
         do ibox=laddr(1,ilev),laddr(2,ilev)
            if(nlist3(ibox).gt.0) then
              cntlist4=cntlist4+1
              list4ct(ibox)=cntlist4
              ilist4(cntlist4)=ibox
            endif
         enddo
      enddo
      if(ifprint.ge.1) print *,"nboxes:",nboxes,"cntlist4:",cntlist4
      allocate(pgboxwexp(nd,nexptotp,cntlist4,6))
      allocate(gboxmexp(nd*(nterms(nlevels)+1)* &
     &    (2*nterms(nlevels)+1),8,cntlist4))



      allocate(gboxsubcenters(3,8,nthd))
      allocate(gboxfl(2,8,nthd))

      nmaxt = 0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,istart,iend,npts) &
!$OMP REDUCTION(max:nmaxt)
      do ibox=1,nboxes
        if(list4ct(ibox).gt.0) then
          istart = isrcse(1,ibox)
          iend = isrcse(2,ibox)
          npts = iend-istart+1
          if(npts.gt.nmaxt) nmaxt = npts
        endif
      enddo
!$OMP END PARALLEL DO

      allocate(gboxind(nmaxt,nthd))
      allocate(gboxsort(3,nmaxt,nthd))
      allocate(gboxwexp(nd,nexptotp,6,8,nthd))
      allocate(gboxcgsort(nd,nmaxt,nthd))
      allocate(gboxdpsort(nd,3,nmaxt,nthd))

!   note gboxmexp is an array not scalar
      pgboxwexp=0d0
      gboxmexp=0d0

!     form mexp for all list4 type box at first ghost box center
      do ilev=1,nlevels-1

         rscpow(0) = 1.0d0/boxsize(ilev+1)
         rtmp = rscales(ilev+1)/boxsize(ilev+1)
         do i=1,nterms(ilev+1)
            rscpow(i) = rscpow(i-1)*rtmp
         enddo

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,istart,iend,jbox,jstart,jend,npts,npts0,i) &
!$OMP PRIVATE(ithd)
         do ibox=laddr(1,ilev),laddr(2,ilev)
            ithd = 0
!$          ithd=omp_get_thread_num()
            ithd = ithd + 1
            if(list4ct(ibox).gt.0) then
              istart=isrcse(1,ibox)
              iend=isrcse(2,ibox)
              npts = iend-istart+1

              if(npts.gt.0) then
                call subdividebox(sourcesort(1,istart),npts, &
     &    centers(1,ibox),boxsize(ilev+1), &
     &    gboxind(1,ithd),gboxfl(1,1,ithd), &
     &    gboxsubcenters(1,1,ithd))
                call dreorderf(3,npts,sourcesort(1,istart), &
     &    gboxsort(1,1,ithd),gboxind(1,ithd))
                if(ifcharge.eq.1) then
                  call dreorderf(nd,npts,chargesort(1,istart), &
     &    gboxcgsort(1,1,ithd),gboxind(1,ithd))
                endif
                if(ifdipole.eq.1) then
                  call dreorderf(3*nd,npts,dipvecsort(1,1,istart), &
     &    gboxdpsort(1,1,1,ithd),gboxind(1,ithd))
                endif
                do i=1,8
                  if(gboxfl(1,i,ithd).gt.0) then
                    jstart=gboxfl(1,i,ithd)
                    jend=gboxfl(2,i,ithd)
                    npts0=jend-jstart+1
                    jbox=list4ct(ibox)
                    if(ifcharge.eq.1.and.ifdipole.eq.0) then
                      call l3dformmpc(nd,rscales(ilev+1), &
     &    gboxsort(1,jstart,ithd), &
     &    gboxcgsort(1,jstart,ithd), &
     &    npts0,gboxsubcenters(1,i,ithd),nterms(ilev+1), &
     &    gboxmexp(1,i,jbox),wlege,nlege)          
                    endif
                    if(ifcharge.eq.0.and.ifdipole.eq.1) then
                      call l3dformmpd(nd,rscales(ilev+1), &
     &    gboxsort(1,jstart,ithd), &
     &    gboxdpsort(1,1,jstart,ithd), &
     &    npts0,gboxsubcenters(1,i,ithd),nterms(ilev+1), &
     &    gboxmexp(1,i,jbox),wlege,nlege)          
                    endif
                    if(ifcharge.eq.1.and.ifdipole.eq.1) then
                      call l3dformmpcd(nd,rscales(ilev+1), &
     &    gboxsort(1,jstart,ithd), &
     &    gboxcgsort(1,jstart,ithd), &
     &    gboxdpsort(1,1,jstart,ithd), &
     &    npts0,gboxsubcenters(1,i,ithd),nterms(ilev+1), &
     &    gboxmexp(1,i,jbox),wlege,nlege)          
                    endif
                    call l3dmpmp(nd,rscales(ilev+1), &
     &    gboxsubcenters(1,i,ithd),gboxmexp(1,i,jbox), &
     &    nterms(ilev+1),rscales(ilev),centers(1,ibox), &
     &    rmlexp(iaddr(1,ibox)),nterms(ilev),dc,lca)
     
                    call mpscale(nd,nterms(ilev+1),gboxmexp(1,i,jbox), &
     &    rscpow,tmp(1,0,-nmax,ithd))
!
!c                process up down for current box
!
                    call mpoletoexp(nd,tmp(1,0,-nmax,ithd), &
     &    nterms(ilev+1),nlams, &
     &    nfourier,nexptot,mexpf1(1,1,ithd), &
     &    mexpf2(1,1,ithd),rlsc)

                    call ftophys(nd,mexpf1(1,1,ithd), &
     &    nlams,rlams,nfourier, &
     &    nphysical,nthmax,gboxwexp(1,1,1,i,ithd), &
     &    fexpe,fexpo)

                    call ftophys(nd,mexpf2(1,1,ithd), &
     &    nlams,rlams,nfourier, &
     &    nphysical,nthmax,gboxwexp(1,1,2,i,ithd), &
     &    fexpe,fexpo)

                    call processgboxudexp(nd,gboxwexp(1,1,1,i,ithd), &
     &    gboxwexp(1,1,2,i,ithd),i,nexptotp, &
     &    pgboxwexp(1,1,jbox,1),pgboxwexp(1,1,jbox,2), &
     &    xshift,yshift,zshift)
!
!c                process north-south for current box
!
                    call rotztoy(nd,nterms(ilev+1),tmp(1,0,-nmax,ithd), &
     &    mptmp(1,ithd),rdminus)
                    call mpoletoexp(nd,mptmp(1,ithd), &
     &    nterms(ilev+1),nlams, &
     &    nfourier,nexptot,mexpf1(1,1,ithd), &
     &    mexpf2(1,1,ithd),rlsc)

                    call ftophys(nd,mexpf1(1,1,ithd), &
     &    nlams,rlams,nfourier, &
     &    nphysical,nthmax,gboxwexp(1,1,3,i,ithd), &
     &    fexpe,fexpo)

                    call ftophys(nd,mexpf2(1,1,ithd), &
     &    nlams,rlams,nfourier, &
     &    nphysical,nthmax,gboxwexp(1,1,4,i,ithd), &
     &    fexpe,fexpo)

                    call processgboxnsexp(nd,gboxwexp(1,1,3,i,ithd), &
     &    gboxwexp(1,1,4,i,ithd),i,nexptotp, &
     &    pgboxwexp(1,1,jbox,3),pgboxwexp(1,1,jbox,4), &
     &    xshift,yshift,zshift)
!
!c                process east-west for current box
!
                    call rotztox(nd,nterms(ilev+1),tmp(1,0,-nmax,ithd), &
     &    mptmp(1,ithd),rdplus)
                    call mpoletoexp(nd,mptmp(1,ithd), &
     &    nterms(ilev+1),nlams, &
     &    nfourier,nexptot,mexpf1(1,1,ithd), &
     &    mexpf2(1,1,ithd),rlsc)

                    call ftophys(nd,mexpf1(1,1,ithd), &
     &    nlams,rlams,nfourier, &
     &    nphysical,nthmax,gboxwexp(1,1,5,i,ithd), &
     &    fexpe,fexpo)

                    call ftophys(nd,mexpf2(1,1,ithd), &
     &    nlams,rlams,nfourier, &
     &    nphysical,nthmax,gboxwexp(1,1,6,i,ithd), &
     &    fexpe,fexpo)
                
                    call processgboxewexp(nd,gboxwexp(1,1,5,i,ithd), &
     &    gboxwexp(1,1,6,i,ithd),i,nexptotp, &
     &    pgboxwexp(1,1,jbox,5),pgboxwexp(1,1,jbox,6), &
     &    xshift,yshift,zshift)
                  endif
                enddo
              endif
            endif
         enddo
!$OMP END PARALLEL DO
      enddo
      deallocate(gboxfl,gboxsubcenters,gboxwexp,gboxcgsort)
      deallocate(gboxdpsort,gboxind,gboxsort)

      call cpu_time(time2)
!$    time2=omp_get_wtime()
      if(ifprint.ge.1) print *,"mexp list4 time:",time2-time1
      timeinfo(3)=time2-time1
!     end of count number of boxes are in list4
!

!
!
      if(ifprint .ge. 1)  &
     &   call prinf('=== STEP 1 (form mp) ====*',i,0)
        call cpu_time(time1)
!$        time1=omp_get_wtime()
!
!       ... step 1, locate all charges, assign them to boxes, and
!       form multipole expansions


      do ilev=2,nlevels
         if(ifcharge.eq.1.and.ifdipole.eq.0) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,npts,istart,iend,nchild)
            do ibox=laddr(1,ilev),laddr(2,ilev)

               istart = isrcse(1,ibox)
               iend = isrcse(2,ibox) 
               npts = iend-istart+1

               nchild = itree(ipointer(4)+ibox-1)

               if(npts.gt.0.and.nchild.eq.0.and.list4ct(ibox).eq.0) then
                  call l3dformmpc(nd,rscales(ilev), &
     &    sourcesort(1,istart),chargesort(1,istart),npts, &
     &    centers(1,ibox),nterms(ilev), &
     &    rmlexp(iaddr(1,ibox)),wlege,nlege)          
               endif
            enddo
!$OMP END PARALLEL DO
         endif

         if(ifcharge.eq.0.and.ifdipole.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,npts,istart,iend,nchild)
            do ibox=laddr(1,ilev),laddr(2,ilev)

               istart = isrcse(1,ibox) 
               iend = isrcse(2,ibox) 
               npts = iend-istart+1

               nchild = itree(ipointer(4)+ibox-1)

               if(npts.gt.0.and.nchild.eq.0.and.list4ct(ibox).eq.0) then
                  call l3dformmpd(nd,rscales(ilev), &
     &    sourcesort(1,istart), &
     &    dipvecsort(1,1,istart),npts, &
     &    centers(1,ibox),nterms(ilev), &
     &    rmlexp(iaddr(1,ibox)),wlege,nlege)          
               endif
            enddo
!$OMP END PARALLEL DO
         endif

         if(ifdipole.eq.1.and.ifcharge.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,npts,istart,iend,nchild)
            do ibox=laddr(1,ilev),laddr(2,ilev)

               istart = isrcse(1,ibox) 
               iend = isrcse(2,ibox)
               npts = iend-istart+1

               nchild = itree(ipointer(4)+ibox-1)

               if(npts.gt.0.and.nchild.eq.0.and.list4ct(ibox).eq.0) then
                  call l3dformmpcd(nd,rscales(ilev), &
     &    sourcesort(1,istart),chargesort(1,istart), &
     &    dipvecsort(1,1,istart),npts, &
     &    centers(1,ibox),nterms(ilev), &
     &    rmlexp(iaddr(1,ibox)),wlege,nlege)          
               endif
            enddo
!$OMP END PARALLEL DO
         endif
      enddo
      if(ifprint.ge.1) print *,"nboxes:",nboxes,"leaf:",cntlist4



      call cpu_time(time2)
!$    time2=omp_get_wtime()
      timeinfo(1)=time2-time1

      lca = 4*nmax


!       
      if(ifprint .ge. 1) &
     &      call prinf('=== STEP 2 (merge mp) ====*',i,0)
      call cpu_time(time1)
!$    time1=omp_get_wtime()
!
      do ilev=nlevels-1,0,-1
!$OMP PARALLEL DO DEFAULT(SHARED) &
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
!$OMP END PARALLEL DO
      enddo

      call cpu_time(time2)
!$    time2=omp_get_wtime()
      timeinfo(2)=time2-time1

      if(ifprint.ge.1) &
     &    call prinf('=== Step 3 (mp to loc+formta+mpeval) ===*',i,0)
!      ... step 3, convert multipole expansions into local
!       expansions

      call cpu_time(time1)
!$        time1=omp_get_wtime()

!
!c     zero out mexp
! 

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(i,j,k,idim)
      do k=1,6
        do i=1,nboxes
          do j=1,nexptotp
            do idim=1,nd
              mexp(idim,j,i,k) = 0.0d0
            enddo
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

!     init uall,dall,...,etc arrays
      allocate(uall(200,nthd),dall(200,nthd),nall(120,nthd))
      allocate(sall(120,nthd),eall(72,nthd),wall(72,nthd))
      allocate(u1234(36,nthd),d5678(36,nthd),n1256(24,nthd))
      allocate(s3478(24,nthd))
      allocate(e1357(16,nthd),w2468(16,nthd),n12(20,nthd))
      allocate(n56(20,nthd),s34(20,nthd),s78(20,nthd))
      allocate(e13(20,nthd),e57(20,nthd),w24(20,nthd),w68(20,nthd))
      allocate(e1(20,nthd),e3(5,nthd),e5(5,nthd),e7(5,nthd))
      allocate(w2(5,nthd),w4(5,nthd),w6(5,nthd),w8(5,nthd))
      allocate(iboxsubcenters(3,8,nthd))
      allocate(iboxfl(2,8,nthd))
!
!  figure out allocations needed for iboxsrc,iboxsrcind,iboxpot
!  and so on
!
      nmaxt = 0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,istart,iend,npts) &
!$OMP REDUCTION(max:nmaxt)
      do ibox=1,nboxes
        if(nlist3(ibox).gt.0) then
          istart = isrcse(1,ibox)
          iend = isrcse(2,ibox)
          npts = iend-istart+1
          if(npts.gt.nmaxt) nmaxt = npts

          istart = itargse(1,ibox)
          iend = itargse(2,ibox)
          npts = iend - istart + 1
          if(npts.gt.nmaxt) nmaxt = npts
        endif
      enddo
!$OMP END PARALLEL DO

      allocate(iboxsrcind(nmaxt,nthd))
      allocate(iboxsrc(3,nmaxt,nthd))
      allocate(iboxpot(nd,nmaxt,nthd))
      allocate(iboxgrad(nd,3,nmaxt,nthd))
      allocate(iboxhess(nd,6,nmaxt,nthd))

      do ilev=2,nlevels
        allocate(iboxlexp(nd*(nterms(ilev)+1)* &
     &    (2*nterms(ilev)+1),8,nthd))

         rscpow(0) = 1.0d0/boxsize(ilev)
         rtmp = rscales(ilev)/boxsize(ilev)
         do i=1,nterms(ilev)
            rscpow(i) = rscpow(i-1)*rtmp
         enddo

!$OMP PARALLEL DO DEFAULT (SHARED) &
!$OMP PRIVATE(ibox,istart,iend,npts) &
!$OMP PRIVATE(ithd)
         do ibox=laddr(1,ilev),laddr(2,ilev)
            ithd = 0
!$          ithd=omp_get_thread_num()
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
!$OMP END PARALLEL DO
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
!$OMP PARALLEL DO DEFAULT (SHARED) &
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
           ithd = 0
!$         ithd=omp_get_thread_num()
           ithd = ithd + 1
           npts = 0
           if(ifpghtarg.gt.0) then
             istart = itargse(1,ibox)
             iend = itargse(2,ibox) 
             npts = npts + iend-istart+1
           endif

           istart = iexpcse(1,ibox) 
           iend = iexpcse(2,ibox) 
           npts = npts + iend-istart+1

           nchild = itree(ipointer(4)+ibox-1)

           if(ifpgh.gt.0) then
             istart = isrcse(1,ibox) 
             iend = isrcse(2,ibox) 
             npts = npts + iend-istart+1
           endif


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
              do i=1,8
                call mpzero(nd,iboxlexp(1,i,ithd),nterms(ilev))
              enddo

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

              if(ifpgh.eq.1) then
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
                  call dreorderf(nd,npts,pot(1,istart), &
     &    iboxpot(1,1,ithd),iboxsrcind(1,ithd))
                  do i=1,8
                    if(iboxfl(1,i,ithd).gt.0) then
                      jstart=iboxfl(1,i,ithd)
                      jend=iboxfl(2,i,ithd)
                      npts0=jend-jstart+1
                      if(npts0.gt.0) then
                        call l3dtaevalp(nd,rscales(ilev), &
     &    iboxsubcenters(1,i,ithd),iboxlexp(1,i,ithd), &
     &    nterms(ilev),iboxsrc(1,jstart,ithd),npts0, &
     &    iboxpot(1,jstart,ithd),wlege,nlege)
                      endif
                    endif
                  enddo
                  call dreorderi(nd,npts,iboxpot(1,1,ithd), &
     &    pot(1,istart),iboxsrcind(1,ithd))
                endif
              endif

              if(ifpgh.eq.2) then
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
                  call dreorderf(nd,npts,pot(1,istart), &
     &    iboxpot(1,1,ithd),iboxsrcind(1,ithd))
                  call dreorderf(3*nd,npts,grad(1,1,istart), &
     &    iboxgrad(1,1,1,ithd),iboxsrcind(1,ithd))
                  do i=1,8
                    if(iboxfl(1,i,ithd).gt.0) then
                      jstart=iboxfl(1,i,ithd)
                      jend=iboxfl(2,i,ithd)
                      npts0=jend-jstart+1
                      if(npts0.gt.0) then
                        call l3dtaevalg(nd,rscales(ilev), &
     &    iboxsubcenters(1,i,ithd),iboxlexp(1,i,ithd), &
     &    nterms(ilev),iboxsrc(1,jstart,ithd),npts0, &
     &    iboxpot(1,jstart,ithd), &
     &    iboxgrad(1,1,jstart,ithd),wlege,nlege)
                      endif
                    endif
                  enddo
                  call dreorderi(nd,npts,iboxpot(1,1,ithd), &
     &    pot(1,istart),iboxsrcind(1,ithd))
                  call dreorderi(3*nd,npts,iboxgrad(1,1,1,ithd), &
     &    grad(1,1,istart),iboxsrcind(1,ithd))
                endif
              endif
!
!  continue from here
!
              

              if(ifpgh.eq.3) then
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
                  call dreorderf(nd,npts,pot(1,istart), &
     &    iboxpot(1,1,ithd),iboxsrcind(1,ithd))
                  call dreorderf(3*nd,npts,grad(1,1,istart), &
     &    iboxgrad(1,1,1,ithd),iboxsrcind(1,ithd))
                  call dreorderf(6*nd,npts,hess(1,1,istart), &
     &    iboxhess(1,1,1,ithd),iboxsrcind(1,ithd))
           
                  do i=1,8
                    if(iboxfl(1,i,ithd).gt.0) then
                      jstart=iboxfl(1,i,ithd)
                      jend=iboxfl(2,i,ithd)
                      npts0=jend-jstart+1
                      if(npts0.gt.0) then
                        call l3dtaevalh(nd,rscales(ilev), &
     &    iboxsubcenters(1,i,ithd),iboxlexp(1,i,ithd), &
     &    nterms(ilev),iboxsrc(1,jstart,ithd),npts0, &
     &    iboxpot(1,jstart,ithd), &
     &    iboxgrad(1,1,jstart,ithd), &
     &    iboxhess(1,1,jstart,ithd),scarray(1,ilev))
                      endif
                    endif
                  enddo
                  call dreorderi(nd,npts,iboxpot(1,1,ithd), &
     &    pot(1,istart),iboxsrcind(1,ithd))
                  call dreorderi(3*nd,npts,iboxgrad(1,1,1,ithd), &
     &    grad(1,1,istart),iboxsrcind(1,ithd))
                  call dreorderi(6*nd,npts,iboxhess(1,1,1,ithd), &
     &    hess(1,1,istart),iboxsrcind(1,ithd))
                endif
              endif


              if(ifpghtarg.eq.1) then
                istart = itargse(1,ibox) 
                iend = itargse(2,ibox) 
                npts = iend-istart+1
                if(npts.gt.0) then
                  call subdividebox(targsort(1,istart),npts, &
     &    centers(1,ibox),boxsize(ilev), &
     &    iboxsrcind(1,ithd),iboxfl(1,1,ithd), &
     &    iboxsubcenters(1,1,ithd))
                  call dreorderf(3,npts,targsort(1,istart), &
     &    iboxsrc(1,1,ithd),iboxsrcind(1,ithd))
                  call dreorderf(nd,npts,pottarg(1,istart), &
     &    iboxpot(1,1,ithd),iboxsrcind(1,ithd))
                  do i=1,8
                    if(iboxfl(1,i,ithd).gt.0) then
                      jstart=iboxfl(1,i,ithd)
                      jend=iboxfl(2,i,ithd)
                      npts0=jend-jstart+1
                      if(npts0.gt.0) then
                        call l3dtaevalp(nd,rscales(ilev), &
     &    iboxsubcenters(1,i,ithd),iboxlexp(1,i,ithd), &
     &    nterms(ilev),iboxsrc(1,jstart,ithd),npts0, &
     &    iboxpot(1,jstart,ithd),wlege,nlege)
                      endif
                    endif
                  enddo
                  call dreorderi(nd,npts,iboxpot(1,1,ithd), &
     &    pottarg(1,istart),iboxsrcind(1,ithd))
                endif
              endif

              if(ifpghtarg.eq.2) then
                istart = itargse(1,ibox) 
                iend = itargse(2,ibox) 
                npts = iend-istart+1
                if(npts.gt.0) then
                  call subdividebox(targsort(1,istart),npts, &
     &    centers(1,ibox),boxsize(ilev), &
     &    iboxsrcind(1,ithd),iboxfl(1,1,ithd), &
     &    iboxsubcenters(1,1,ithd))
                  call dreorderf(3,npts,targsort(1,istart), &
     &    iboxsrc(1,1,ithd),iboxsrcind(1,ithd))
                  call dreorderf(nd,npts,pottarg(1,istart), &
     &    iboxpot(1,1,ithd),iboxsrcind(1,ithd))
                  call dreorderf(3*nd,npts,gradtarg(1,1,istart), &
     &    iboxgrad(1,1,1,ithd),iboxsrcind(1,ithd))
                  do i=1,8
                    if(iboxfl(1,i,ithd).gt.0) then
                      jstart=iboxfl(1,i,ithd)
                      jend=iboxfl(2,i,ithd)
                      npts0=jend-jstart+1
                      if(npts0.gt.0) then
                        call l3dtaevalg(nd,rscales(ilev), &
     &    iboxsubcenters(1,i,ithd),iboxlexp(1,i,ithd), &
     &    nterms(ilev),iboxsrc(1,jstart,ithd),npts0, &
     &    iboxpot(1,jstart,ithd), &
     &    iboxgrad(1,1,jstart,ithd),wlege,nlege)
                      endif
                    endif
                  enddo
                  call dreorderi(nd,npts,iboxpot(1,1,ithd), &
     &    pottarg(1,istart),iboxsrcind(1,ithd))
                  call dreorderi(3*nd,npts,iboxgrad(1,1,1,ithd), &
     &    gradtarg(1,1,istart),iboxsrcind(1,ithd))
                endif
              endif

              if(ifpghtarg.eq.3) then
                istart = itargse(1,ibox) 
                iend = itargse(2,ibox) 
                npts = iend-istart+1
                if(npts.gt.0) then
                  call subdividebox(targsort(1,istart),npts, &
     &    centers(1,ibox),boxsize(ilev), &
     &    iboxsrcind(1,ithd),iboxfl(1,1,ithd), &
     &    iboxsubcenters(1,1,ithd))
                  call dreorderf(3,npts,targsort(1,istart), &
     &    iboxsrc(1,1,ithd),iboxsrcind(1,ithd))
                  call dreorderf(nd,npts,pottarg(1,istart), &
     &    iboxpot(1,1,ithd),iboxsrcind(1,ithd))
                  call dreorderf(3*nd,npts,gradtarg(1,1,istart), &
     &    iboxgrad(1,1,1,ithd),iboxsrcind(1,ithd))
                  call dreorderf(6*nd,npts,hesstarg(1,1,istart), &
     &    iboxhess(1,1,1,ithd),iboxsrcind(1,ithd))
                  do i=1,8
                    if(iboxfl(1,i,ithd).gt.0) then
                      jstart=iboxfl(1,i,ithd)
                      jend=iboxfl(2,i,ithd)
                      npts0=jend-jstart+1
                      if(npts0.gt.0) then
                        call l3dtaevalh(nd,rscales(ilev), &
     &    iboxsubcenters(1,i,ithd),iboxlexp(1,i,ithd), &
     &    nterms(ilev),iboxsrc(1,jstart,ithd),npts0, &
     &    iboxpot(1,jstart,ithd), &
     &    iboxgrad(1,1,jstart,ithd), &
     &    iboxhess(1,1,jstart,ithd),scarray(1,ilev))
                      endif
                    endif
                  enddo
                  call dreorderi(nd,npts,iboxpot(1,1,ithd), &
     &    pottarg(1,istart),iboxsrcind(1,ithd))
                  call dreorderi(3*nd,npts,iboxgrad(1,1,1,ithd), &
     &    gradtarg(1,1,istart),iboxsrcind(1,ithd))
                  call dreorderi(6*nd,npts,iboxhess(1,1,1,ithd), &
     &    hesstarg(1,1,istart),iboxsrcind(1,ithd))
                endif
              endif

            endif
         enddo
!$OMP END PARALLEL DO
        deallocate(iboxlexp)  
      enddo

      deallocate(iboxsrcind,iboxsrc,iboxpot,iboxgrad,iboxhess)
      deallocate(iboxsubcenters,iboxfl)
      deallocate(uall,dall,nall,sall,eall,wall)
      deallocate(u1234,d5678,n1256,s3478)
      deallocate(e1357,w2468,n12,n56,s34,s78)
      deallocate(e13,e57,w24,w68)
      deallocate(e1,e3,e5,e7,w2,w4,w6,w8)
      deallocate(tmp,mptmp)
      call cpu_time(time2)
!$        time2=omp_get_wtime()
      timeinfo(3) = timeinfo(3) + time2-time1


      if(ifprint.ge.1) &
     &    call prinf('=== Step 4 (split loc) ===*',i,0)

      call cpu_time(time1)
!$        time1=omp_get_wtime()
      do ilev = 2,nlevels-1

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,i,jbox,istart,iend,npts)
         do ibox = laddr(1,ilev),laddr(2,ilev)

            npts = 0

            if(ifpghtarg.gt.0) then
               istart = itargse(1,ibox)
               iend = itargse(2,ibox) 
               npts = npts + iend-istart+1
            endif

            istart = iexpcse(1,ibox) 
            iend = iexpcse(2,ibox) 
            npts = npts + iend-istart+1

            if(ifpgh.gt.0) then
               istart = isrcse(1,ibox)
               iend = isrcse(2,ibox) 
               npts = npts + iend-istart+1
            endif

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
!$OMP END PARALLEL DO
      enddo
      call cpu_time(time2)
!$        time2=omp_get_wtime()
      timeinfo(4) = time2-time1


      if(ifprint.ge.1) &
     &    call prinf('=== step 5 (eval lo) ===*',i,0)

!     ... step 6, evaluate all local expansions
!

      call cpu_time(time1)
!$        time1=omp_get_wtime()
!

!
!c       shift local expansion to local epxanion at expansion centers
!        (note: this part is not relevant for particle codes.
!        it is relevant only for qbx codes)

      do ilev = 0,nlevels
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,nchild,istart,iend,i) &
!$OMP SCHEDULE(DYNAMIC)
         do ibox = laddr(1,ilev),laddr(2,ilev)
            nchild=itree(ipointer(4)+ibox-1)
            if(nchild.eq.0) then 
               istart = iexpcse(1,ibox) 
               iend = iexpcse(2,ibox) 
               do i=istart,iend

                  call l3dlocloc(nd,rscales(ilev), &
     &    centers(1,ibox),rmlexp(iaddr(2,ibox)), &
     &    nterms(ilev),rscales(ilev),expcsort(1,i), &
     &    tsort(1,0,-ntj,i),ntj,dc,lca)
               enddo
            endif
         enddo
!$OMP END PARALLEL DO
      enddo

!
!c        evaluate local expansion at source and target
!         locations
!
      do ilev = 0,nlevels
        if(ifpgh.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,nchild,istart,iend,i,npts) &
!$OMP SCHEDULE(DYNAMIC)
          do ibox = laddr(1,ilev),laddr(2,ilev)
            nchild=itree(ipointer(4)+ibox-1)
            if(nchild.eq.0) then 
              istart = isrcse(1,ibox) 
              iend = isrcse(2,ibox)
              npts = iend-istart+1
              call l3dtaevalp(nd,rscales(ilev),centers(1,ibox), &
     &    rmlexp(iaddr(2,ibox)),nterms(ilev),sourcesort(1,istart), &
     &    npts,pot(1,istart),wlege,nlege)
            endif
          enddo
!$OMP END PARALLEL DO
        endif

        if(ifpgh.eq.2) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,nchild,istart,iend,i,npts) &
!$OMP SCHEDULE(DYNAMIC)
          do ibox = laddr(1,ilev),laddr(2,ilev)
            nchild=itree(ipointer(4)+ibox-1)
            if(nchild.eq.0) then 
              istart = isrcse(1,ibox) 
              iend = isrcse(2,ibox)
              npts = iend-istart+1
              call l3dtaevalg(nd,rscales(ilev),centers(1,ibox), &
     &    rmlexp(iaddr(2,ibox)),nterms(ilev),sourcesort(1,istart), &
     &    npts,pot(1,istart),grad(1,1,istart),wlege,nlege)
            endif
          enddo
!$OMP END PARALLEL DO
        endif


        if(ifpgh.eq.3) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,nchild,istart,iend,i,npts) &
!$OMP SCHEDULE(DYNAMIC)
          do ibox = laddr(1,ilev),laddr(2,ilev)
            nchild=itree(ipointer(4)+ibox-1)
            if(nchild.eq.0) then 
              istart = isrcse(1,ibox)
              iend = isrcse(2,ibox)
              npts = iend-istart+1
              call l3dtaevalh(nd,rscales(ilev),centers(1,ibox), &
     &    rmlexp(iaddr(2,ibox)),nterms(ilev),sourcesort(1,istart), &
     &    npts,pot(1,istart),grad(1,1,istart),hess(1,1,istart), &
     &    scarray(1,ilev))
            endif
          enddo
!$OMP END PARALLEL DO
        endif

        if(ifpghtarg.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,nchild,istart,iend,i,npts) &
!$OMP SCHEDULE(DYNAMIC)
          do ibox = laddr(1,ilev),laddr(2,ilev)
            nchild=itree(ipointer(4)+ibox-1)
            if(nchild.eq.0) then 
              istart = itargse(1,ibox)
              iend = itargse(2,ibox)
              npts = iend-istart+1
              call l3dtaevalp(nd,rscales(ilev),centers(1,ibox), &
     &    rmlexp(iaddr(2,ibox)),nterms(ilev),targsort(1,istart), &
     &    npts,pottarg(1,istart),wlege,nlege)
            endif
          enddo
!$OMP END PARALLEL DO
        endif

        if(ifpghtarg.eq.2) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,nchild,istart,iend,i,npts) &
!$OMP SCHEDULE(DYNAMIC)
          do ibox = laddr(1,ilev),laddr(2,ilev)
            nchild=itree(ipointer(4)+ibox-1)
            if(nchild.eq.0) then 
              istart = itargse(1,ibox)
              iend = itargse(2,ibox)
              npts = iend-istart+1

              call l3dtaevalg(nd,rscales(ilev),centers(1,ibox), &
     &    rmlexp(iaddr(2,ibox)),nterms(ilev),targsort(1,istart), &
     &    npts,pottarg(1,istart),gradtarg(1,1,istart),wlege,nlege)
            endif
          enddo
!$OMP END PARALLEL DO
        endif

        if(ifpghtarg.eq.3) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,nchild,istart,iend,i,npts) &
!$OMP SCHEDULE(DYNAMIC)
          do ibox = laddr(1,ilev),laddr(2,ilev)
            nchild=itree(ipointer(4)+ibox-1)
            if(nchild.eq.0) then 
              istart = itargse(1,ibox)
              iend = itargse(2,ibox)
              npts = iend-istart+1

              call l3dtaevalh(nd,rscales(ilev),centers(1,ibox), &
     &    rmlexp(iaddr(2,ibox)),nterms(ilev),targsort(1,istart), &
     &    npts,pottarg(1,istart),gradtarg(1,1,istart), &
     &    hesstarg(1,1,istart),scarray(1,ilev))
            endif
          enddo
!$OMP END PARALLEL DO
        endif
      enddo

    
      call cpu_time(time2)
!$        time2=omp_get_wtime()
      timeinfo(5) = time2 - time1


      if(ifprint .ge. 1) &
     &     call prinf('=== STEP 6 (direct) =====*',i,0)
      call cpu_time(time1)
!$        time1=omp_get_wtime()

      if(ifnear.eq.0) goto 1000
!
!c       directly form local expansions for list1 sources
!        at expansion centers. 
!        (note: this part is not relevant for particle codes.
!         It is relevant only for qbx codes)


      do ilev=0,nlevels
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,istarte,iende,i,jbox) &
!$OMP PRIVATE(jstart,jend)
         do ibox = laddr(1,ilev),laddr(2,ilev)
            istarte = iexpcse(1,ibox) 
            iende = iexpcse(2,ibox) 

            

            do i =1,nlist1(ibox)
               jbox = list1(i,ibox)
               jstart = isrcse(1,jbox) 
               jend = isrcse(2,jbox) 

               call lfmm3dexpc_direct_tree(nd,jstart,jend,istarte, &
     &    iende,sourcesort,ifcharge,chargesort,ifdipole, &
     &    dipvecsort,expcsort,tsort,scjsort,ntj, &
     &    wlege,nlege)
            enddo
         enddo
!$OMP END PARALLEL DO
      enddo

!
!c        directly evaluate potential at sources and targets 
!         due to sources in list1

      do ilev=0,nlevels
!
!c           evaluate at the sources
!

        if(ifpgh.eq.1) then
          if(ifcharge.eq.1.and.ifdipole.eq.0) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,istarts,iends,npts0,i,jbox,jstart,jend,npts) &
!$OMP SCHEDULE(DYNAMIC)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istarts = isrcse(1,ibox) 
              iends = isrcse(2,ibox)
              npts0 = iends-istarts+1

              do i=1,nlist1(ibox)
                jbox = list1(i,ibox) 
                jstart =  isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectcp(nd,sourcesort(1,jstart), &
     &    chargesort(1,jstart),npts,sourcesort(1,istarts), &
     &    npts0,pot(1,istarts),thresh)          
              enddo
            enddo
!$OMP END PARALLEL DO
          endif

          if(ifcharge.eq.0.and.ifdipole.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,istarts,iends,npts0,i,jbox,jstart,jend,npts) &
!$OMP SCHEDULE(DYNAMIC)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istarts = isrcse(1,ibox)
              iends = isrcse(2,ibox)
              npts0 = iends-istarts+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart =  isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectdp(nd,sourcesort(1,jstart), &
     &    dipvecsort(1,1,jstart),npts,sourcesort(1,istarts), &
     &    npts0,pot(1,istarts),thresh)          
              enddo
            enddo
!$OMP END PARALLEL DO
          endif

          if(ifcharge.eq.1.and.ifdipole.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,istarts,iends,npts0,i,jbox,jstart,jend,npts) &
!$OMP SCHEDULE(DYNAMIC)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istarts = isrcse(1,ibox)
              iends = isrcse(2,ibox)
              npts0 = iends-istarts+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectcdp(nd,sourcesort(1,jstart), &
     &    chargesort(1,jstart), &
     &    dipvecsort(1,1,jstart),npts,sourcesort(1,istarts), &
     &    npts0,pot(1,istarts),thresh)          
              enddo
            enddo
!$OMP END PARALLEL DO
          endif
        endif

        if(ifpgh.eq.2) then
          if(ifcharge.eq.1.and.ifdipole.eq.0) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,istarts,iends,npts0,i,jbox,jstart,jend,npts) &
!$OMP SCHEDULE(DYNAMIC)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istarts = isrcse(1,ibox)
              iends = isrcse(2,ibox)
              npts0 = iends-istarts+1
              

              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectcg(nd,sourcesort(1,jstart), &
     &    chargesort(1,jstart),npts,sourcesort(1,istarts), &
     &    npts0,pot(1,istarts),grad(1,1,istarts),thresh)   
              enddo
            enddo
!$OMP END PARALLEL DO
          endif

          if(ifcharge.eq.0.and.ifdipole.eq.1) then

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,istarts,iends,npts0,i,jbox,jstart,jend,npts) &
!$OMP SCHEDULE(DYNAMIC)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istarts = isrcse(1,ibox)
              iends = isrcse(2,ibox)
              npts0 = iends-istarts+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectdg(nd,sourcesort(1,jstart), &
     &    dipvecsort(1,1,jstart),npts,sourcesort(1,istarts), &
     &    npts0,pot(1,istarts),grad(1,1,istarts),thresh)     
              enddo
            enddo
!$OMP END PARALLEL DO
          endif

          if(ifcharge.eq.1.and.ifdipole.eq.1) then

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,istarts,iends,npts0,i,jbox,jstart,jend,npts) &
!$OMP SCHEDULE(DYNAMIC)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istarts = isrcse(1,ibox)
              iends = isrcse(2,ibox)
              npts0 = iends-istarts+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectcdg(nd,sourcesort(1,jstart), &
     &    chargesort(1,jstart), &
     &    dipvecsort(1,1,jstart),npts,sourcesort(1,istarts), &
     &    npts0,pot(1,istarts),grad(1,1,istarts),thresh)      
              enddo
            enddo
!$OMP END PARALLEL DO
          endif
        endif


        if(ifpgh.eq.3) then
          if(ifcharge.eq.1.and.ifdipole.eq.0) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,istarts,iends,npts0,i,jbox,jstart,jend,npts) &
!$OMP SCHEDULE(DYNAMIC)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istarts = isrcse(1,ibox)
              iends = isrcse(2,ibox)
              npts0 = iends-istarts+1
              

              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectch(nd,sourcesort(1,jstart), &
     &    chargesort(1,jstart),npts,sourcesort(1,istarts), &
     &    npts0,pot(1,istarts),grad(1,1,istarts), &
     &    hess(1,1,istarts),thresh)   
              enddo
            enddo
!$OMP END PARALLEL DO
          endif

          if(ifcharge.eq.0.and.ifdipole.eq.1) then

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,istarts,iends,npts0,i,jbox,jstart,jend,npts) &
!$OMP SCHEDULE(DYNAMIC)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istarts = isrcse(1,ibox)
              iends = isrcse(2,ibox)
              npts0 = iends-istarts+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectdh(nd,sourcesort(1,jstart), &
     &    dipvecsort(1,1,jstart),npts,sourcesort(1,istarts), &
     &    npts0,pot(1,istarts),grad(1,1,istarts), &
     &    hess(1,1,istarts),thresh)          
              enddo
            enddo
!$OMP END PARALLEL DO
          endif

          if(ifcharge.eq.1.and.ifdipole.eq.1) then

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,istarts,iends,npts0,i,jbox,jstart,jend,npts) &
!$OMP SCHEDULE(DYNAMIC)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istarts = isrcse(1,ibox)
              iends = isrcse(2,ibox)
              npts0 = iends-istarts+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectcdh(nd,sourcesort(1,jstart), &
     &    chargesort(1,jstart), &
     &    dipvecsort(1,1,jstart),npts,sourcesort(1,istarts), &
     &    npts0,pot(1,istarts),grad(1,1,istarts), &
     &    hess(1,1,istarts),thresh)          
              enddo
            enddo
!$OMP END PARALLEL DO
          endif
        endif

        if(ifpghtarg.eq.1) then
          if(ifcharge.eq.1.and.ifdipole.eq.0) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,istartt,iendt,npts0,i,jbox,jstart,jend,npts) &
!$OMP SCHEDULE(DYNAMIC)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istartt = itargse(1,ibox)
              iendt = itargse(2,ibox)
              npts0 = iendt-istartt+1
              

              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectcp(nd,sourcesort(1,jstart), &
     &    chargesort(1,jstart),npts,targsort(1,istartt), &
     &    npts0,pottarg(1,istartt),thresh)          
              enddo
            enddo
!$OMP END PARALLEL DO
          endif

          if(ifcharge.eq.0.and.ifdipole.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,istartt,iendt,npts0,i,jbox,jstart,jend,npts) &
!$OMP SCHEDULE(DYNAMIC)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istartt = itargse(1,ibox)
              iendt = itargse(2,ibox)
              npts0 = iendt-istartt+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectdp(nd,sourcesort(1,jstart), &
     &    dipvecsort(1,1,jstart),npts,targsort(1,istartt), &
     &    npts0,pottarg(1,istartt),thresh)          
              enddo
            enddo
!$OMP END PARALLEL DO
          endif

          if(ifcharge.eq.1.and.ifdipole.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,istartt,iendt,npts0,i,jbox,jstart,jend,npts) &
!$OMP SCHEDULE(DYNAMIC)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istartt = itargse(1,ibox)
              iendt = itargse(2,ibox)
              npts0 = iendt-istartt+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectcdp(nd,sourcesort(1,jstart), &
     &    chargesort(1,jstart), &
     &    dipvecsort(1,1,jstart),npts,targsort(1,istartt), &
     &    npts0,pottarg(1,istartt),thresh)          
              enddo
            enddo
!$OMP END PARALLEL DO
          endif
        endif

        if(ifpghtarg.eq.2) then
          if(ifcharge.eq.1.and.ifdipole.eq.0) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,istartt,iendt,npts0,i,jbox,jstart,jend,npts) &
!$OMP SCHEDULE(DYNAMIC)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istartt = itargse(1,ibox)
              iendt = itargse(2,ibox)
              npts0 = iendt-istartt+1
              

              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectcg(nd,sourcesort(1,jstart), &
     &    chargesort(1,jstart),npts,targsort(1,istartt), &
     &    npts0,pottarg(1,istartt),gradtarg(1,1,istartt), &
     &    thresh)   
              enddo
            enddo
!$OMP END PARALLEL DO
          endif

          if(ifcharge.eq.0.and.ifdipole.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,istartt,iendt,npts0,i,jbox,jstart,jend,npts) &
!$OMP SCHEDULE(DYNAMIC)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istartt = itargse(1,ibox)
              iendt = itargse(2,ibox)
              npts0 = iendt-istartt+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectdg(nd,sourcesort(1,jstart), &
     &    dipvecsort(1,1,jstart),npts,targsort(1,istartt), &
     &    npts0,pottarg(1,istartt),gradtarg(1,1,istartt), &
     &    thresh)          
              enddo
            enddo
!$OMP END PARALLEL DO
          endif

          if(ifcharge.eq.1.and.ifdipole.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,istartt,iendt,npts0,i,jbox,jstart,jend,npts) &
!$OMP SCHEDULE(DYNAMIC)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istartt = itargse(1,ibox)
              iendt = itargse(2,ibox)
              npts0 = iendt-istartt+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectcdg(nd,sourcesort(1,jstart), &
     &    chargesort(1,jstart), &
     &    dipvecsort(1,1,jstart),npts,targsort(1,istartt), &
     &    npts0,pottarg(1,istartt),gradtarg(1,1,istartt), &
     &    thresh)          
              enddo
            enddo
!$OMP END PARALLEL DO
          endif
        endif

        if(ifpghtarg.eq.3) then
          if(ifcharge.eq.1.and.ifdipole.eq.0) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,istartt,iendt,npts0,i,jbox,jstart,jend,npts) &
!$OMP SCHEDULE(DYNAMIC)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istartt = itargse(1,ibox)
              iendt = itargse(2,ibox)
              npts0 = iendt-istartt+1
              

              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectch(nd,sourcesort(1,jstart), &
     &    chargesort(1,jstart),npts,targsort(1,istartt), &
     &    npts0,pottarg(1,istartt),gradtarg(1,1,istartt), &
     &    hesstarg(1,1,istartt),thresh)   
              enddo
            enddo
!$OMP END PARALLEL DO
          endif

          if(ifcharge.eq.0.and.ifdipole.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,istartt,iendt,npts0,i,jbox,jstart,jend,npts) &
!$OMP SCHEDULE(DYNAMIC)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istartt = itargse(1,ibox)
              iendt = itargse(2,ibox)
              npts0 = iendt-istartt+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectdh(nd,sourcesort(1,jstart), &
     &    dipvecsort(1,1,jstart),npts,targsort(1,istartt), &
     &    npts0,pottarg(1,istartt),gradtarg(1,1,istartt), &
     &    hesstarg(1,1,istartt),thresh)          
              enddo
            enddo
!$OMP END PARALLEL DO
          endif

          if(ifcharge.eq.1.and.ifdipole.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ibox,istartt,iendt,npts0,i,jbox,jstart,jend,npts) &
!$OMP SCHEDULE(DYNAMIC)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istartt = itargse(1,ibox)
              iendt = itargse(2,ibox)
              npts0 = iendt-istartt+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectcdh(nd,sourcesort(1,jstart), &
     &    chargesort(1,jstart), &
     &    dipvecsort(1,1,jstart),npts,targsort(1,istartt), &
     &    npts0,pottarg(1,istartt),gradtarg(1,1,istartt), &
     &    hesstarg(1,1,istartt),thresh)          
              enddo
            enddo
!$OMP END PARALLEL DO
          endif
        endif
      enddo
 1000 continue      
 
      call cpu_time(time2)
!$        time2=omp_get_wtime()
      timeinfo(6) = time2-time1
      if(ifprint.ge.1) call prin2('timeinfo=*',timeinfo,6)
      d = 0
      do i = 1,6
         d = d + timeinfo(i)
      enddo

      if(ifprint.ge.1) call prin2('sum(timeinfo)=*',d,1)

      return
      end
!------------------------------------------------
      subroutine lfmm3dexpc_direct_tree(nd,istart,iend,jstart,jend, &
     &     source,ifcharge,charge,ifdipole, &
     &     dipvec,expc,texps,scj,ntj,wlege,nlege)
!--------------------------------------------------------------------
!     This subroutine adds the local expansions due to sources
!     istart to iend in the source array at the expansion centers
!     jstart to jend in the expansion center array to the existing 
!     local expansions at the corresponding expansion centers.
!
!     INPUT arguments
!-------------------------------------------------------------------
!     nd           in: integer
!                  number of charge densities
!
!     istart       in:Integer
!                  Starting index in source array whose expansions
!                  we wish to add
!
!     iend         in:Integer
!                  Last index in source array whose expansions
!                  we wish to add
!
!     jstart       in: Integer
!                  First index in the expansion center array at 
!                  which we  wish to compute the expansions
! 
!     jend         in:Integer
!                  Last index in expansion center array at 
!                  which we wish to compute the expansions
! 
!     scjsort      in: double precision(*)
!                  Scale of expansions formed at the expansion centers
!
!     source       in: double precision(3,ns)
!                  Source locations
!
!     ifcharge     in: Integer
!                  flag for including expansions due to charges
!                  The expansion due to charges will be included
!                  if ifcharge == 1
!
!     charge       in: double precision
!                  Charge at the source locations
!
!     ifdipole     in: Integer
!                 flag for including expansions due to dipoles
!                 The expansion due to dipoles will be included
!                 if ifdipole == 1
!
!     dipvec      in: double precision(3,ns)
!                 Dipole orientation vector at the source locations
!
!     expc        in: double precision(3,nexpc)
!                 Expansion center locations
!
!     ntj         in: Integer
!                 Number of terms in expansion
!
!     wlege       in: double precision(0:nlege,0:nlege)
!                 precomputed array of recurrence relation
!                 coeffs for Ynm calculation.
!
!    nlege        in: integer
!                 dimension parameter for wlege
!------------------------------------------------------------
!     OUTPUT
!
!   Updated expansions at the targs
!   texps       out: double complex(0:ntj,-ntj:ntj,expc) 
!                 coeffs for local expansions
!-------------------------------------------------------               
        implicit none
!
        integer istart,iend,jstart,jend,ns,j, nlege
        integer ifcharge,ifdipole,ier,nd
        double precision source(3,*)
        double precision scj(*)
        double precision wlege(*)
        double precision charge(nd,*)
        double precision dipvec(nd,3,*)
        double precision expc(3,*)

        integer nlevels,ntj
!
        double complex texps(nd,0:ntj,-ntj:ntj,*)
        
!
        ns = iend - istart + 1
        if(ifcharge.eq.1.and.ifdipole.eq.0) then
          do j=jstart,jend
            call l3dformtac(nd,scj(j), &
     &    source(1,istart),charge(1,istart),ns, &
     &    expc(1,j),ntj,texps(1,0,-ntj,j),wlege,nlege)
           enddo
         endif

         if(ifcharge.eq.0.and.ifdipole.eq.1) then
          do j=jstart,jend
            call l3dformtad(nd,scj(j), &
     &    source(1,istart), &
     &    dipvec(1,1,istart),ns,expc(1,j),ntj,texps(1,0,-ntj,j), &
     &    wlege,nlege)
           enddo
         endif

         if(ifcharge.eq.1.and.ifdipole.eq.1) then
          do j=jstart,jend
            call l3dformtacd(nd,scj(j), &
     &    source(1,istart),charge(1,istart), &
     &    dipvec(1,1,istart),ns,expc(1,j),ntj,texps(1,0,-ntj,j), &
     &    wlege,nlege)
           enddo
         endif

!
        return
        end