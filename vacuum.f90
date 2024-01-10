! 8/24/23
! This file contains taylor_vaccum_solver(), a subroutine that solves a system 
! of integral equations to compute the magnetic field corresponding with Taylor
! states in a toroidal geometry. Other subroutines here are helpers. 

subroutine lpcomp_virtualcasing_addsub_complex(npatches, norders, ixyzs, &
    iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
    eps, nnz, row_ptr, col_ind, iquad, nquad, wnear, rjvec, rho, &
    novers, nptso, ixyzso, srcover, whtsover, curlj, gradrho)
    ! 
    ! Wrapper for lpcomp_virtualcasing_addsub that accommodates complex 
    ! input/output.
    !
    implicit none

    integer npatches, npts, ndtarg, ntarg, nnz, nquad, nptso
    integer norders(npatches), ixyzs(npatches+1), iptype(npatches)
    integer row_ptr(npts+1), col_ind(nnz), iquad(nnz+1)
    integer novers(npatches), ixyzso(npatches+1)
    
    integer i, j

    real *8 eps
    real *8 srccoefs(9,npts),srcvals(12,npts),targs(ndtarg,ntarg)
    real *8 wnear(nquad,3),srcover(12,nptso),whtsover(nptso)
    real *8 curljre(3,ntarg),curljim(3,ntarg)
    real *8 gradrhore(3,ntarg),gradrhoim(3,ntarg)

    complex *16 rjvec(3,npts), rho(npts), curlj(3,ntarg), gradrho(3,ntarg)
    complex *16 ima
    data ima/(0.0d0,1.0d0)/

    ! real part
    call lpcomp_virtualcasing_addsub(npatches,norders,ixyzs,&
        iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
        eps,nnz,row_ptr,col_ind,iquad,nquad,wnear,real(rjvec),real(rho),&
        novers,nptso,ixyzso,srcover,whtsover,curljre,gradrhore)
    ! imaginary part
    call lpcomp_virtualcasing_addsub(npatches,norders,ixyzs,&
        iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
        eps,nnz,row_ptr,col_ind,iquad,nquad,wnear,aimag(rjvec),aimag(rho),&
        novers,nptso,ixyzso,srcover,whtsover,curljim,gradrhoim)
    ! $OMP PARALLEL DO DEFAULT(SHARED)
    do i = 1,3
        do j = 1,ntarg
            curlj(i,j) = curljre(i,j) + ima*curljim(i,j)
        enddo
    enddo
    ! $OMP END PARALLEL DO 
    ! $OMP PARALLEL DO DEFAULT(SHARED)
    do i = 1,3
        do j = 1,ntarg
            gradrho(i,j) = gradrhore(i,j) + ima*gradrhoim(i,j)
        enddo
    enddo
    ! $OMP END PARALLEL DO 
    return 
end subroutine lpcomp_virtualcasing_addsub_complex

subroutine lpcomp_lap_comb_dir_addsub_complex_vec(nd, npatches, norders, ixyzs,&
    iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
    eps, dpars, nnz, row_ptr, col_ind, iquad, nquad, wnear, sigma, &
    novers, nptso, ixyzso, srcover, whtsover, pot)

    implicit none 
    integer, intent(in) :: nd,npatches,npts,ndtarg,ntarg
    integer, intent(in) :: norders(npatches),ixyzs(npatches+1)
    integer, intent(in) :: ixyzso(npatches+1),iptype(npatches)
    real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
    real *8, intent(in) :: targs(ndtarg,ntarg)
    real *8, intent(in) :: dpars(2)
    integer, intent(in) :: nnz,row_ptr(ntarg+1),col_ind(nnz),nquad
    integer, intent(in) :: iquad(nnz+1)
    real *8, intent(in) :: wnear(nquad)
    complex *16, intent(in) :: sigma(nd,npts)
    integer, intent(in) :: novers(npatches+1)
    integer, intent(in) :: nptso
    real *8, intent(in) :: srcover(12,nptso),whtsover(nptso)
    complex *16, intent(out) :: pot(nd,ntarg)
    
    integer i,j
    real *8 pot_real_temp(ntarg),pot_imag_temp(ntarg)

    complex *16 ima
    data ima/(0.0d0,1.0d0)/

    do i = 1,nd
        ! pot_real_temp = 0
        ! pot_imag_temp = 0
        call lpcomp_lap_comb_dir_addsub(npatches,norders,ixyzs,&
            iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,eps,&
            dpars,nnz,row_ptr,col_ind,iquad,nquad,wnear,real(sigma(i,1:npts)),&
            novers,nptso,ixyzso,srcover,whtsover,pot_real_temp)
        call lpcomp_lap_comb_dir_addsub(npatches,norders,ixyzs,&
            iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,eps,&
            dpars,nnz,row_ptr,col_ind,iquad,nquad,wnear,aimag(sigma(i,1:npts)),&
            novers,nptso,ixyzso,srcover,whtsover,pot_imag_temp)
        do j = 1,ntarg
            pot(i,j) = pot_real_temp(j) + ima*pot_imag_temp(j)
        enddo
    enddo

    return
end subroutine lpcomp_lap_comb_dir_addsub_complex_vec

subroutine fun_surf_interp_complex(nd,npatches,norders,ixyzs,iptype,npts,&
    f,na,apatches,auv,finterp)
    !
    ! Wrapper for virtual-casing/src/surf_routs.f90:fun_surf_interp() that 
    ! handles complex I/O
    !

    implicit none
    integer, intent(in) :: nd,npatches,npts
    integer, intent(in) :: norders(npatches),ixyzs(npatches+1),iptype(npatches)
    complex *16, intent(in) :: f(nd,npts)
    integer, intent(in) :: na,apatches(na)
    real *8, intent(in) :: auv(2,na)
    complex *16, intent(out) :: finterp(nd,na)

    real *8 frealintp(nd,na),fimagintp(nd,na)
    integer i,j

    complex *16 ima
    data ima/(0.0d0,1.0d0)/

    call fun_surf_interp(nd,npatches,norders,ixyzs,iptype,npts,real(f),na,&
        apatches,auv,frealintp)
    call fun_surf_interp(nd,npatches,norders,ixyzs,iptype,npts,aimag(f),na,&
        apatches,auv,fimagintp)
    do i=1,na
        do j=1,nd
            finterp(j,i) = frealintp(j,i) + ima*fimagintp(j,i)
        enddo
    enddo
    return
end subroutine fun_surf_interp_complex

subroutine taylor_vaccuum_solver(npatches,norders,ixyzs,iptype,npts,&
    srccoefs,srcvals,ipars,eps,numit,flux,eps_gmres,solvetest,B0,&
    niter1,niter2,errs1,errs2,rres1,rres2,sigma,alpha,B)
    ! 
    ! Solves the Taylor state BIE 
    !     
    !          -sigma     dS_0                 ___
    !     0 = -------- - ------[sigma] + i(n • \ / x S_0[m_H])
    !             2        dn                   v
    ! 
    ! subject to the flux constraint 
    ! 
    !      /\
    !      |
    !      O S_0[n x B] • dl = flux,
    !      |
    !     \/
    ! 
    ! with m_H represented as 
    ! 
    !           N_s
    !           ---
    !           \
    !     m_H = /   alpha_k m_H^k,    
    !           ---
    !           k=1
    ! 
    ! where {m_H^1, ..., m_H^(N_s)} is a basis of harmonic surface vector 
    ! fields, and B on the surface defined as 
    ! 
    !          -sigma      m_H    ___                ___
    !     B = --------n + ----- - \ / S_0[sigma] + i \ / x S_0[m_H].
    !             2         2      v                  v  
    !
    ! The solution is provided as a vector of sigma evaluations and the 
    ! coefficients alpha_k. The linear system is solved iteratively using GMRES.
    ! For now, we deal with the genus-one case (N_s = 1), e.g. a torus. 
    ! 
    ! References:
    !   [1] D. Malhotra et. al., "Taylor states in stellarators: A fast 
    !       high-order boundary integral solver"
    ! 
    ! ==========================================================================
    ! 
    ! Input: 
    !   npatches - integer
    !     number of patches 
    ! 
    !   norders- integer(npatches)
    !     order of discretization on each patch 
    !
    !   ixyzs - integer(npatches+1)
    !     ixyzs(i) denotes the starting location in srccoefs,
    !     and srcvals array corresponding to patch i
    !
    !   iptype - integer(npatches)
    !     type of patch
    !     iptype = 1, triangular patch discretized using RV nodes
    !
    !   npts - integer
    !     total number of discretization points on the boundary
    !
    !   srccoefs - real *8 (9,npts)
    !     koornwinder expansion coefficients of xyz, dxyz/du, and dxyz/dv on 
    !     each patch. 
    !     For each point srccoefs(1:3,i) is xyz info
    !                    srccoefs(4:6,i) is dxyz/du info
    !                    srccoefs(7:9,i) is dxyz/dv info
    !
    !   srcvals - real *8 (12,npts)
    !     xyz(u,v) and derivative info sampled at the 
    !     discretization nodes on the surface
    !     srcvals(1:3,i) - xyz info
    !     srcvals(4:6,i) - dxyz/du info
    !     srcvals(7:9,i) - dxyz/dv info
    !     srcvals(10:12,i) - nxyz info
    ! 
    !   ipars - integer(2)
    !     parameters for surface discretization
    !
    !   eps - real *8
    !     precision requested for computing quadrature and fmm
    !     tolerance
    !
    !   numit - integer
    !     max number of gmres iterations
    ! 
    !   flux - real *8
    !     right-hand side of the flux condition
    !
    !   eps_gmres - real *8
    !     gmres tolerance requested
    ! 
    !   solvetest - integer
    !     if 0, use the solver as intended
    !     if 1, test the solver using a reference field B0 (see 4.3 in [1])
    ! 
    !   B0 - complex *16(3,npts)
    !     reference field for testing (see 4.3 in [1])
    ! 
    ! ==========================================================================
    !
    ! Output:
    !   niter1 - integer
    !     number of GMRES iterations required for relative residual when 
    !       computing A^-1 B0 . n
    !     ignore this if solvetest = 0
    ! 
    !   niter2 - integer
    !     number of GMRES iterations required for relative residual when
    !       computing A^-1 u_1 (see 4.3.1 in [1])
    !      
    !   errs1 - real *8(numit+1)
    !     relative residual as a function of iteration number for first GMRES
    ! 
    !   errs2 - real *8(numit+1)
    !     relative residual as a function of iteration number for second GMRES
    !
    !   rres1 - real *8
    !     relative residual for computed solution for first GMRES
    ! 
    !   rres2 - real *8
    !     relative residual for computed solution for second GMRES
    !          
    !   sigma - complex *16(npts)
    !     solution density 
    !   
    !   alpha - complex *16
    !     solution coefficient on harmonic surface vector field
    ! 
    !   B - complex *16
    !     magnetic field
    ! 
    !   fluxcheck - complex *16
    !     computed flux, for testing purposes
    !
    implicit none

    ! constants
    real *8 done,pi

    ! surface discretization
    integer ndtarg,ntarg
    real *8, allocatable :: targs(:,:),uvs_targ(:,:)
    integer, allocatable :: ipatch_id(:)
    integer iptype(npatches)
    integer iptype_avg,norder_avg
    real *8 rfac,rfac0
    real *8, allocatable :: cms(:,:),rads(:),rad_near(:) 

    ! near quadrature correction
    integer nnz, ikerorder
    integer, allocatable :: row_ptr(:),col_ind(:),iquad(:) 
    real *8 dpars(2)
    integer, allocatable :: ixyzso(:),novers(:)
    complex *16 ztmp
    integer nptso
    real *8, allocatable :: srcover(:,:),wover(:)
    integer nquad
    real *8, allocatable :: wnear(:,:)
    integer iquadtype
    real *8 t1,t2

    ! circulation (A-cycle) integral 
    integer m, na, nb 
    real *8, allocatable :: avals(:,:),bvals(:,:)
    real *8, allocatable :: awts(:),bwts(:)
    real *8, allocatable :: auv(:,:),buv(:,:)
    integer, allocatable :: apatches(:),bpatches(:)

    ! test solver with a reference field B0
    integer solvetest
    complex *16 B0(3,npts)
    real *8 nxB0(3,npts),nxB0temp(npts),SnxB0(3,npts),SnxB0temp(npts)
    real *8, allocatable :: SnxB0intp(:,:)
    real *8 AcycSnxB0, flux
    complex *16 rhs(npts) 

    complex *16, allocatable :: SnxgradSDintp(:,:)
    complex *16, allocatable :: SnxgradSwintp(:,:)
    complex *16, allocatable :: SmHnxcurlSmHintp(:,:)

    integer npatches,npts
    integer norders(npatches),ixyzs(npatches+1)
    real *8 srccoefs(9,npts),srcvals(12,npts),eps,eps_gmres
    complex *16 soln(npts)

    real *8 errs1(numit+1),errs2(numit+1)
    real *8 rres1,rres2
    integer niter1,niter2

    integer i,j

    integer ipars(2)

    ! gmres variables
    real *8 did,rb,wnrm2
    integer numit,it,iind,it1,k,l
    real *8 rmyerr
    complex *16 dtmp,temp,w(npts)
    complex *16 vmat(npts,numit+1),hmat(numit,numit)
    complex *16 cs(numit),sn(numit),svec(numit+1),yvec(numit+1)
    complex *16 wtmp(npts),gradSsigma(3,npts),gradSw(3,npts),curlSmH(3,npts)
    complex *16 mH(3,npts),mHcomps(1,2,npts)
    complex *16 jtmp(3,npts),curlSjtemp(3,npts),rhs_gmres(npts)

    ! post-gmres variables
    complex *16 nxgradSD(3,npts),SnxgradSD(3,npts),AcycSnxgradSD
    complex *16 nxgradSw(3,npts),SnxgradSw(3,npts),AcycSnxgradSw
    complex *16 mHnxcurlSmH(3,npts),SmHnxcurlSmH(3,npts),AcycSmHnxcurlSmH

    complex *16 sigma(npts),alpha,B(3,npts),nxB(3,npts),SnxB(3,npts),fluxcheck
    complex *16, allocatable :: AcycSnxBintp(:,:)

    complex *16 ima
    data ima/(0.0d0,1.0d0)/

    done = 1
    pi = atan(done)*4

    ! set up targets as on-surface discretization points
    ndtarg = 3
    ntarg = npts
    allocate(targs(ndtarg,npts),uvs_targ(2,ntarg),ipatch_id(ntarg))

    ! $OMP PARALLEL DO DEFAULT(SHARED)
    do i=1,ntarg
        targs(1,i) = srcvals(1,i)
        targs(2,i) = srcvals(2,i)
        targs(3,i) = srcvals(3,i)
        ipatch_id(i) = -1
        uvs_targ(1,i) = 0
        uvs_targ(2,i) = 0
    enddo
    ! $OMP END PARALLEL DO 

    ! initialize patch_id and uv_tar for on-surface targets
    call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts,ipatch_id,uvs_targ)
    iptype_avg = floor(sum(iptype)/(npatches+0.0d0))
    norder_avg = floor(sum(norders)/(npatches+0.0d0))
    call get_rfacs(norder_avg,iptype_avg,rfac,rfac0)
    allocate(cms(3,npatches),rads(npatches),rad_near(npatches))
    call get_centroid_rads(npatches,norders,ixyzs,iptype,npts,srccoefs,cms,rads)
    
    ! $OMP PARALLEL DO DEFAULT(SHARED) 
    do i=1,npatches
        rad_near(i) = rads(i)*rfac
    enddo
    ! $OMP END PARALLEL DO

    ! find near quadrature correction interactions
    print *, "entering find near mem"
    call findnearmem(cms,npatches,rad_near,ndtarg,targs,npts,nnz)
    print *, "nnz=",nnz
    allocate(row_ptr(npts+1),col_ind(nnz))
    call findnear(cms,npatches,rad_near,ndtarg,targs,npts,row_ptr,col_ind)
    allocate(iquad(nnz+1)) 
    call get_iquad_rsc(npatches,ixyzs,npts,nnz,row_ptr,col_ind,iquad)
    ikerorder = -1
    ! set kernel parameters
    dpars(1) = 1.0d0
    dpars(2) = 0
    if(abs(dpars(2)).gt.1.0d-16) ikerorder = 0

    ! estimate oversampling for far field and oversample geometry
    allocate(novers(npatches),ixyzso(npatches+1))
    print *, "beginning far order estimation"
    ztmp = 0
    call get_far_order(eps,npatches,norders,ixyzs,iptype,cms,rads,npts,&
        srccoefs,ndtarg,npts,targs,ikerorder,ztmp,nnz,row_ptr,col_ind,rfac,&
        novers,ixyzso)
    nptso = ixyzso(npatches+1)-1
    print *, "nptso=",nptso
    allocate(srcover(12,nptso),wover(nptso))
    call oversample_geom(npatches,norders,ixyzs,iptype,npts,srccoefs,srcvals,&
        novers,ixyzso,nptso,srcover)
    call get_qwts(npatches,novers,ixyzso,iptype,nptso,srcover,wover)

    ! compute near quadrature correction
    nquad = iquad(nnz+1)-1
    print *, "nquad=",nquad
    allocate(wnear(nquad,3))
    
    ! $OMP PARALLEL DO DEFAULT(SHARED)      
    do i=1,nquad
        do j = 1,3
            wnear(i,j) = 0
        enddo
    enddo
    ! $OMP END PARALLEL DO

    iquadtype = 1
    print *, "starting to generate near quadrature"
    call cpu_time(t1)
    call getnearquad_magnetostatics(npatches,norders,ixyzs, &
    iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,ipatch_id, &
    uvs_targ,eps,iquadtype,nnz,row_ptr,col_ind,iquad,rfac0,nquad, &
    wnear)

    ! reading near quadrature from file 
    ! open(50,form='unformatted',file='wnear_fat.dat',action='read')
    ! read(50) wnear
    ! close(50)

    call cpu_time(t2)
    print *, 'quadrature generation time=',t2-t1
    print *, "done generating near quadrature"

    ! writing near quadrature to file 
    ! open(50,file='wnear_fat.dat',form='unformatted')
    ! write(50) wnear
    ! close(50)
    ! print *, 'done writing near quadrature to file'
    
    ! get A-cycle parametrization
    m = 64
    na = ipars(2)*m
    nb = ipars(1)*m 
    allocate(avals(9,na),awts(na),auv(2,na),apatches(na))
    allocate(bvals(9,nb),bwts(nb),buv(2,nb),bpatches(nb))
    call get_ab_cycles_torusparam(npatches,norders,ixyzs,iptype,npts,&
        srccoefs,srcvals,ipars,m,na,avals,awts,apatches,auv,nb,bvals,&
        bwts,bpatches,buv)
    
    ! rhs computation for testing
    if(solvetest.eq.1) then
        do i=1,npts
            nxB0(1,i) = srcvals(11,i)*real(B0(3,i)) - srcvals(12,i)*real(B0(2,i))
            nxB0(2,i) = srcvals(12,i)*real(B0(1,i)) - srcvals(10,i)*real(B0(3,i))
            nxB0(3,i) = srcvals(10,i)*real(B0(2,i)) - srcvals(11,i)*real(B0(1,i))
        enddo
        ! compute S[n x B0]
        do i = 1,3
            do j = 1,npts
                nxB0temp(j) = nxB0(i,j)
            enddo
            call lpcomp_lap_comb_dir_addsub(npatches,norders,&
                ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,eps,&
                dpars,nnz,row_ptr,col_ind,iquad,nquad,wnear,nxB0temp,&
                novers,nptso,ixyzso,srcover,wover,SnxB0temp)
            do k=1,npts
                SnxB0(i,k) = SnxB0temp(k)
            enddo
        enddo
        ! evaluate on A-cycle
        allocate(SnxB0intp(3,na))
        call fun_surf_interp(3,npatches,norders,ixyzs,iptype,npts,&
            SnxB0,na,apatches,auv,SnxB0intp)
        ! integrate along A-cycle
        AcycSnxB0 = 0
        do i = 1,na
            AcycSnxB0 = AcycSnxB0 + ((SnxB0intp(1,i))*avals(4,i) + &
                SnxB0intp(2,i)*avals(5,i) + & 
                SnxB0intp(3,i)*avals(6,i)) * awts(i)
        enddo
        flux = AcycSnxB0
        do i=1,npts
            rhs(i) = -(B0(1,i)*srcvals(10,i) + B0(2,i)*srcvals(11,i) + &
                B0(3,i)*srcvals(12,i))  
        enddo
    endif
    print *, "now starting gmres 1"

    ! start GMRES for w = A11 \ rhs

    ! $OMP PARALLEL DO DEFAULT(SHARED)
    do i=1,npts
        vmat(i,1) = 0
    enddo
    ! $OMP END PARALLEL DO 

    did = -1d0/2
    niter1 = 0
    rb = 0
    do i = 1,numit
        cs(i) = 0
        sn(1) = 0
    enddo
    ! $OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:rb)
    do i=1,npts
        rb = rb + abs(rhs(i))**2
    enddo
    ! $OMP END PARALLEL DO
    rb = sqrt(rb)
    ! $OMP PARALLEL DO DEFAULT(SHARED)
    do i=1,npts
        w(i) = 0
    enddo
    ! $OMP END PARALLEL DO
    ! only do first GMRES if rhs is nonzero
    if (abs(rb).gt.1.0d-16) then 
        ! $OMP PARALLEL DO DEFAULT(SHARED)
        do i=1,npts
            vmat(i,1) = rhs(i)/rb
        enddo
        ! $OMP END PARALLEL DO 
        svec(1) = rb

        ! dummy argument (current) for layer potential evaluation
        ! $OMP PARALLEL DO DEFAULT(SHARED)
        do i=1,3
            do j=1,npts
                jtmp(i,j) = 0
                curlSjtemp(i,j) = 0
            enddo
        enddo
        ! $OMP PARALLEL DO DEFAULT(SHARED)

        ! GMRES iteration
        do it = 1,numit
            it1 = it+1
            ! layer potential evaluation
            call lpcomp_virtualcasing_addsub_complex(npatches,norders,ixyzs,&
            iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,eps,nnz,&
            row_ptr,col_ind,iquad,nquad,wnear,jtmp,vmat(1:npts,it),&
            novers,nptso,ixyzso,srcover,wover,curlSjtemp,gradSsigma)
            ! $OMP PARALLEL DO DEFAULT(SHARED)
            do i=1,npts
                wtmp(i) = -gradSsigma(1,i)*srcvals(10,i)
                wtmp(i) = wtmp(i) - gradSsigma(2,i)*srcvals(11,i)
                wtmp(i) = wtmp(i) - gradSsigma(3,i)*srcvals(12,i)
            enddo
            ! $OMP END PARALLEL DO

            do k = 1,it
                dtmp = 0
                ! $OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:dtmp)
                do j = 1,npts
                    dtmp = dtmp + wtmp(j)*vmat(j,k)
                enddo
                ! $OMP END PARALLEL DO
                hmat(k,it) = dtmp
                ! $OMP PARALLEL DO DEFAULT(SHARED) 
                do j = 1,npts
                    wtmp(j) = wtmp(j) - hmat(k,it)*vmat(j,k)
                enddo
                ! $OMP END PARALLEL DO 
            enddo

            hmat(it,it) = hmat(it,it) + did
            wnrm2 = 0
            ! $OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:dtmp)
            do j = 1,npts
                wnrm2 = wnrm2 + abs(wtmp(j))**2
            enddo
            ! $OMP END PARALLEL DO
            wnrm2 = sqrt(wnrm2)

            ! $OMP PARALLEL DO DEFAULT(SHARED) 
            do j = 1,npts
                vmat(j,it1) = wtmp(j)/wnrm2
            enddo
            ! $OMP END PARALLEL DO 

            ! adapted for complex system
            do k=1,it-1
                temp = cs(k)*hmat(k,it)+sn(k)*hmat(k+1,it)
                hmat(k+1,it) = -sn(k)*hmat(k,it)+cs(k)*hmat(k+1,it)
                hmat(k,it) = temp
            enddo
    
            dtmp = wnrm2
    
            call rotmat_complex_gmres(hmat(it,it),dtmp,cs(it),sn(it))
                
            hmat(it,it) = cs(it)*hmat(it,it)+sn(it)*wnrm2
            svec(it1) = -sn(it)*svec(it)
            svec(it) = cs(it)*svec(it)
            rmyerr = abs(svec(it1))/rb
            errs1(it) = rmyerr
            print *, "iter=",it,errs1(it)
    
            if(rmyerr.le.eps_gmres.or.it.eq.numit) then
                ! solve the linear system corresponding to upper triangular part 
                ! of hmat to obtain yvec
                ! y = triu(H(1:it,1:it))\s(1:it)
                do j = 1,it
                    iind = it-j+1
                    yvec(iind) = svec(iind)
                    do l=iind+1,it
                        yvec(iind) = yvec(iind) - hmat(iind,l)*yvec(l)
                    enddo
                    yvec(iind) = yvec(iind)/hmat(iind,iind)
                enddo
                ! estimate w
                ! $OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
                do j = 1,npts
                    w(j) = 0
                    do i=1,it
                        w(j) = w(j) + yvec(i)*vmat(j,i)
                    enddo
                enddo
                ! $OMP END PARALLEL DO          
                rres1 = 0
                ! $OMP PARALLEL DO DEFAULT(SHARED) 
                do j = 1,npts
                    wtmp(j) = 0
                enddo
                ! $OMP END PARALLEL DO 
                ! compute error
                call lpcomp_virtualcasing_addsub_complex(npatches,norders,ixyzs,&
                iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,eps,nnz,row_ptr,&
                col_ind,iquad,nquad,wnear,jtmp,w,novers,nptso,ixyzso,srcover,&
                wover,curlSjtemp,gradSw)
                ! curlSmH = curl S[mH], gradSw = grad S[w]
                ! compute -n . grad S[w]
                ! $OMP PARALLEL DO DEFAULT(SHARED)
                do i=1,npts
                    wtmp(i) = -gradSw(1,i)*srcvals(10,i)
                    wtmp(i) = wtmp(i) - gradSw(2,i)*srcvals(11,i)
                    wtmp(i) = wtmp(i) - gradSw(3,i)*srcvals(12,i)
                enddo
                ! $OMP END PARALLEL DO
                ! $OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:rres)
                do i = 1,npts
                    rres1 = rres1 + abs(did*w(i) + wtmp(i) - rhs(i))**2
                enddo
                ! $OMP END PARALLEL DO       
                rres1 = sqrt(rres1)/rb   
                niter1 = it
                ! Now w = rhs \ A11, ending the first GMRES 
                exit
            endif
        enddo
    endif
    print *,'now starting gmres 2'

    ! start GMRES for D = A11 \ -A12

    ! Find A_12 = i n . curl S_0[m_H] (see [1], section 3.1.4)

    ! compute harmonic surface vector field 
    call axi_surf_harm_field(npts,srcvals,mH,mHcomps)
    ! compute curl S_0[m_H]
    ! This subroutine computes grad S and curl S at the same time, 
    ! so a dummy argument is used for the grad S computation.

    call lpcomp_virtualcasing_addsub_complex(npatches,norders,ixyzs,iptype,&
        npts,srccoefs,srcvals,ndtarg,ntarg,targs,eps,nnz,row_ptr,col_ind,&
        iquad,nquad,wnear,mH,vmat(1:npts,1),novers,nptso,ixyzso,srcover,&
        wover,curlSmH,gradSsigma)
    ! compute -A_12 = rhs_gmres
    ! $OMP PARALLEL DO DEFAULT(SHARED)
    do i = 1,npts
        rhs_gmres(i) = -curlSmH(1,i)*srcvals(10,i) 
        rhs_gmres(i) = rhs_gmres(i) - curlSmH(2,i)*srcvals(11,i) 
        rhs_gmres(i) = rhs_gmres(i) - curlSmH(3,i)*srcvals(12,i) 
        rhs_gmres(i) = ima*rhs_gmres(i) 
    enddo
    ! $OMP END PARALLEL DO 

    niter2 = 0
    ! compute norm of RHS and initialize v
    rb = 0
    do i = 1,numit
        cs(i) = 0
        sn(1) = 0
    enddo
    ! $OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:rb)
    do i=1,npts
        rb = rb + abs(rhs_gmres(i))**2
    enddo
    ! $OMP END PARALLEL DO
    rb = sqrt(rb)
    ! $OMP PARALLEL DO DEFAULT(SHARED)
    do i=1,npts
        vmat(i,1) = rhs_gmres(i)/rb
    enddo
    ! $OMP END PARALLEL DO 
    svec(1) = rb

    ! dummy argument (current) for layer potential evaluation
    ! $OMP PARALLEL DO DEFAULT(SHARED)
    do i=1,3
        do j=1,npts
            jtmp(i,j) = 0
            curlSjtemp(i,j) = 0
        enddo
    enddo
    ! $OMP PARALLEL DO DEFAULT(SHARED)

    ! GMRES iteration
    do it = 1,numit
        it1 = it+1
        ! layer potential evaluation
        call lpcomp_virtualcasing_addsub_complex(npatches,norders,ixyzs,&
        iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,eps,nnz,&
        row_ptr,col_ind,iquad,nquad,wnear,jtmp,vmat(1:npts,it),&
        novers,nptso,ixyzso,srcover,wover,curlSjtemp,gradSsigma)
        ! $OMP PARALLEL DO DEFAULT(SHARED)
        do i=1,npts
            wtmp(i) = -gradSsigma(1,i)*srcvals(10,i)
            wtmp(i) = wtmp(i) - gradSsigma(2,i)*srcvals(11,i)
            wtmp(i) = wtmp(i) - gradSsigma(3,i)*srcvals(12,i)
        enddo
        ! $OMP END PARALLEL DO

        do k = 1,it
            dtmp = 0
            ! $OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:dtmp)
            do j = 1,npts
                dtmp = dtmp + conjg(wtmp(j))*vmat(j,k)
            enddo
            ! $OMP END PARALLEL DO
            hmat(k,it) = dtmp
            ! $OMP PARALLEL DO DEFAULT(SHARED) 
            do j = 1,npts
                wtmp(j) = wtmp(j) - hmat(k,it)*vmat(j,k)
            enddo
            ! $OMP END PARALLEL DO 
        enddo

        hmat(it,it) = hmat(it,it) + did
        wnrm2 = 0
        ! $OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:dtmp)
        do j = 1,npts
            wnrm2 = wnrm2 + abs(wtmp(j))**2
        enddo
        ! $OMP END PARALLEL DO
        wnrm2 = sqrt(wnrm2)

        ! $OMP PARALLEL DO DEFAULT(SHARED) 
        do j = 1,npts
            vmat(j,it1) = wtmp(j)/wnrm2
        enddo
        ! $OMP END PARALLEL DO 

        ! adapted for complex system
        do k=1,it-1
            temp = conjg(cs(k))*hmat(k,it)+conjg(sn(k))*hmat(k+1,it)
            hmat(k+1,it) = -sn(k)*hmat(k,it)+cs(k)*hmat(k+1,it)
            hmat(k,it) = temp
        enddo
  
        dtmp = wnrm2
  
        call rotmat_complex_gmres(hmat(it,it),dtmp,cs(it),sn(it))

        hmat(it,it) = conjg(cs(it))*hmat(it,it)+sn(it)*wnrm2
        svec(it1) = -sn(it)*svec(it)
        svec(it) = cs(it)*svec(it)
        rmyerr = abs(svec(it1))/rb
        errs2(it) = rmyerr
        print *, "iter=",it,errs2(it)
  
        if(rmyerr.le.eps_gmres.or.it.eq.numit) then
            ! solve the linear system corresponding to upper triangular part of 
            ! hmat to obtain yvec
            ! y = triu(H(1:it,1:it))\s(1:it)
            do j = 1,it
                iind = it-j+1
                yvec(iind) = svec(iind)
                do l=iind+1,it
                  yvec(iind) = yvec(iind) - hmat(iind,l)*yvec(l)
                enddo
                yvec(iind) = yvec(iind)/hmat(iind,iind)
            enddo
            ! estimate x
            ! $OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
            do j = 1,npts
                soln(j) = 0
                do i=1,it
                  soln(j) = soln(j) + yvec(i)*vmat(j,i)
                enddo
            enddo
            ! $OMP END PARALLEL DO          
            rres2 = 0
            ! $OMP PARALLEL DO DEFAULT(SHARED) 
            do j = 1,npts
                wtmp(j) = 0
            enddo
            ! $OMP END PARALLEL DO 
            ! compute error
            call lpcomp_virtualcasing_addsub_complex(npatches,norders,ixyzs,&
            iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,eps,nnz,row_ptr,&
            col_ind,iquad,nquad,wnear,mH,soln,novers,nptso,ixyzso,srcover,&
            wover,curlSmH,gradSsigma)
            ! curlSmH = curl S[mH], gradSsigma = grad S[soln]
            ! $OMP PARALLEL DO DEFAULT(SHARED)
            do i=1,npts
                wtmp(i) = -gradSsigma(1,i)*srcvals(10,i)
                wtmp(i) = wtmp(i) - gradSsigma(2,i)*srcvals(11,i)
                wtmp(i) = wtmp(i) - gradSsigma(3,i)*srcvals(12,i)
            enddo
            ! $OMP END PARALLEL DO
            ! $OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:rres)
            do i = 1,npts
                rres2 = rres2 + abs(did*soln(i) + wtmp(i) - rhs_gmres(i))**2
            enddo
            ! $OMP END PARALLEL DO          
            rres2 = sqrt(rres2)/rb
            niter2 = it

            ! Now soln = D (see [1]). We solve for alpha and sigma using (31) 
            ! and (32). 
            ! If there is a rhs, we do the same computations on w as we do on D
            
            ! compute -n x grad S[D]
            ! $OMP PARALLEL DO DEFAULT(SHARED)
            do i = 1,npts
                nxgradSD(1,i) = -srcvals(11,i)*gradSsigma(3,i)&
                    + srcvals(12,i)*gradSsigma(2,i)
                nxgradSD(2,i) = -srcvals(12,i)*gradSsigma(1,i)&
                    + srcvals(10,i)*gradSsigma(3,i)
                nxgradSD(3,i) = -srcvals(10,i)*gradSsigma(2,i)&
                    + srcvals(11,i)*gradSsigma(1,i)
            enddo
            ! $OMP END PARALLEL DO
            ! compute -n x grad S[w]
            ! $OMP PARALLEL DO DEFAULT(SHARED)
            do i = 1,npts
                nxgradSw(1,i) = -srcvals(11,i)*gradSw(3,i)&
                    + srcvals(12,i)*gradSw(2,i)
                nxgradSw(2,i) = -srcvals(12,i)*gradSw(1,i)&
                    + srcvals(10,i)*gradSw(3,i)
                nxgradSw(3,i) = -srcvals(10,i)*gradSw(2,i)&
                    + srcvals(11,i)*gradSw(1,i)
            enddo
            ! $OMP END PARALLEL DO

            ! compute -S[n x grad S[D]]
            call lpcomp_lap_comb_dir_addsub_complex_vec(3,npatches,norders,&
                ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,eps,&
                dpars,nnz,row_ptr,col_ind,iquad,nquad,wnear,nxgradSD,novers,&
                nptso,ixyzso,srcover,wover,SnxgradSD)
            ! compute -S[n x grad S[w]]
            call lpcomp_lap_comb_dir_addsub_complex_vec(3,npatches,norders,&
                ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,eps,&
                dpars,nnz,row_ptr,col_ind,iquad,nquad,wnear,nxgradSw,novers,&
                nptso,ixyzso,srcover,wover,SnxgradSw)

            allocate(SnxgradSDintp(3,na),SnxgradSwintp(3,na))
            ! evaluate -S[n x grad S[D]] on A-cycle
            call fun_surf_interp_complex(3,npatches,norders,ixyzs,iptype,npts,&
                SnxgradSD,na,apatches,auv,SnxgradSDintp)
            ! evaluate -S[n x grad S[w]] on A-cycle
            call fun_surf_interp_complex(3,npatches,norders,ixyzs,iptype,npts,&
                SnxgradSw,na,apatches,auv,SnxgradSwintp)
            AcycSnxgradSD = 0
            ! integrate along A-cycle
            do i = 1,na
                AcycSnxgradSD = AcycSnxgradSD + (SnxgradSDintp(1,i)*avals(4,i) &
                    + SnxgradSDintp(2,i)*avals(5,i) &
                    + SnxgradSDintp(3,i)*avals(6,i))*awts(i)
            enddo
            AcycSnxgradSw = 0
            do i = 1,na
                AcycSnxgradSw = AcycSnxgradSw + (SnxgradSwintp(1,i)*avals(4,i) &
                    + SnxgradSwintp(2,i)*avals(5,i) &
                    + SnxgradSwintp(3,i)*avals(6,i))*awts(i)
            enddo

            ! compute -m_H/2 + n x curl S[mH]
            ! $OMP PARALLEL DO DEFAULT(SHARED)
            do i = 1,npts
                mHnxcurlSmH(1,i) = -mH(1,i)/2 + &
                    srcvals(11,i)*curlSmH(3,i) - srcvals(12,i)*curlSmH(2,i)
                mHnxcurlSmH(2,i) = -mH(2,i)/2 + &
                    srcvals(12,i)*curlSmH(1,i) - srcvals(10,i)*curlSmH(3,i)
                mHnxcurlSmH(3,i) = -mH(3,i)/2 + &
                    srcvals(10,i)*curlSmH(2,i) - srcvals(11,i)*curlSmH(1,i)
            enddo
            ! $OMP END PARALLEL DO
            ! compute S[-m_H/2 + n x curl S[mH]]
            call lpcomp_lap_comb_dir_addsub_complex_vec(3,npatches,norders,&
                ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,eps,&
                dpars,nnz,row_ptr,col_ind,iquad,nquad,wnear,mHnxcurlSmH,&
                novers,nptso,ixyzso,srcover,wover,SmHnxcurlSmH)

            ! evaluate S[-m_H/2 + n x curl S[mH]] on A-cycle
            allocate(SmHnxcurlSmHintp(3,na))
            call fun_surf_interp_complex(3,npatches,norders,ixyzs,iptype,npts,&
                SmHnxcurlSmH,na,apatches,auv,SmHnxcurlSmHintp)

            ! integrate along A-cycle
            AcycSmHnxcurlSmH = 0
            do i = 1,na
                AcycSmHnxcurlSmH = AcycSmHnxcurlSmH + awts(i)*(&
                    SmHnxcurlSmHintp(1,i)*avals(4,i) + &
                    SmHnxcurlSmHintp(2,i)*avals(5,i) + & 
                    SmHnxcurlSmHintp(3,i)*avals(6,i))
            enddo
            ! multiply by i
            AcycSmHnxcurlSmH = ima*AcycSmHnxcurlSmH

            alpha = (flux + AcycSnxgradSw)/(AcycSnxgradSD + AcycSmHnxcurlSmH)
            ! $OMP PARALLEL DO DEFAULT(SHARED)
            do i = 1,npts
                sigma(i) = alpha*soln(i) - w(i)
            enddo
            ! $OMP END PARALLEL DO

            ! compute the magnetic field from sigma and alpha
            ! added 11/2/23
            call lpcomp_virtualcasing_addsub_complex(npatches,norders,ixyzs,&
            iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,eps,nnz,row_ptr,&
            col_ind,iquad,nquad,wnear,mH,sigma,novers,nptso,ixyzso,srcover,&
            wover,curlSmH,gradSsigma)
            ! $OMP PARALLEL DO DEFAULT(SHARED)
            do i = 1,npts
                do j = 1,3
                    B(j,i) = (-sigma(i)*srcvals(j+9,i) + alpha*mH(j,i))/2 &
                        - gradSsigma(j,i) + alpha*ima*curlSmH(j,i)
                enddo
            enddo 
            ! $OMP END PARALLEL DO

            ! compute the flux using the A-cycle integral representation
            ! $OMP PARALLEL DO DEFAULT(SHARED)
            do i = 1,npts
                nxB(1,i) = srcvals(11,i)*B(3,i) - srcvals(12,i)*B(2,i)
                nxB(2,i) = srcvals(12,i)*B(1,i) - srcvals(10,i)*B(3,i)
                nxB(3,i) = srcvals(10,i)*B(2,i) - srcvals(11,i)*B(1,i)
            enddo
            ! $OMP END PARALLEL DO
            call lpcomp_lap_comb_dir_addsub_complex_vec(3,npatches,norders,&
                ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,eps,&
                dpars,nnz,row_ptr,col_ind,iquad,nquad,wnear,nxB,novers,nptso,&
                ixyzso,srcover,wover,SnxB)
            allocate(AcycSnxBintp(3,na))    
            call fun_surf_interp_complex(3,npatches,norders,ixyzs,iptype,npts,&
                SnxB,na,apatches,auv,AcycSnxBintp)
            fluxcheck = 0
            do i = 1,na
                fluxcheck = fluxcheck + awts(i)*(&
                    AcycSnxBintp(1,i)*avals(4,i) + &
                    AcycSnxBintp(2,i)*avals(5,i) + &
                    AcycSnxBintp(3,i)*avals(6,i))
            enddo
            print *,'fluxcheck',fluxcheck
            return 
        endif
    enddo
    return
end subroutine taylor_vaccuum_solver

subroutine axi_surf_harm_field(npts, srcvals, mH, mHcomps)
    !
    ! Computes a particular harmonic surface vector field on an axisymmetric, 
    ! genus-one surface. Harmonicity requires 
    ! 
    !     ___                   ___
    !     \ / _ . m_H = 0  and  \ / _ . n x m_H = 0, 
    !      v |                   v |
    ! 
    ! and we enforce the additional condition
    ! 
    !     n x m_H = im_H. 
    ! 
    ! Such fields are scalar multiples of 
    !              ^          ^
    !            theta     i tau
    !     m_H = ------- + -------. 
    !              r         r
    ! This subroutine computes the above.
    ! 
    ! ==========================================================================
    ! 
    ! Input:
    !   npts - integer
    !     total number of discretization points on the boundary
    ! 
    !   srccoefs - real *8 (9,npts)
    !     koornwinder expansion coefficients of xyz, dxyz/du, and dxyz/dv on 
    !     each patch. 
    !     For each point srccoefs(1:3,i) is xyz info
    !                    srccoefs(4:6,i) is dxyz/du info
    !                    srccoefs(7:9,i) is dxyz/dv info
    !
    !   srcvals - real *8 (12,npts)
    !     xyz(u,v) and derivative info sampled at the 
    !     discretization nodes on the surface
    !     srcvals(1:3,i) - xyz info
    !     srcvals(4:6,i) - dxyz/du info
    !     srcvals(7:9,i) - dxyz/dv info 
    !     srcvals(10:12,i) - nxyz info 
    ! 
    ! ==========================================================================
    ! 
    ! Output:
    !   mH - real *8 (3,npts)
    !     the harmonic surface vector field
    ! 
    !   mHcomps - real *8 (1,2,npts)
    !     theta- and tau-components of mH
    ! 
    implicit none 
    integer npts
    real *8 srcvals(12,npts)

    complex *16 mH(3,npts), mHcomps(1,2,npts)

    integer i   
    real *8 overr, cosphi, sinphi, overu2, overv2
    real *8 tauhat(3)
    complex *16 ima
    data ima/(0.0d0,1.0d0)/

    do i = 1,npts
        overr = srcvals(1,i)**2 + srcvals(2,i)**2
        overr = 1/sqrt(overr)
        cosphi = srcvals(1,i)*overr 
        sinphi = srcvals(2,i)*overr 
        tauhat(1) = -srcvals(12,i)*cosphi
        tauhat(2) = -srcvals(12,i)*sinphi
        tauhat(3) = srcvals(10,i)*cosphi + srcvals(11,i)*sinphi
        ! mHcomps(1,1,i) = overr
        ! mHcomps(1,2,i) = ima*overr
        mH(1,i) = (-sinphi + ima*tauhat(1))*overr
        mH(2,i) = (cosphi + ima*tauhat(2))*overr
        mH(3,i) = ima*tauhat(3)*overr
        mHcomps(1,1,i) = mH(1,i)*srcvals(4,i) + mH(2,i)*srcvals(5,i) 
        mHcomps(1,1,i) = mHcomps(1,1,i) + mH(3,i)*srcvals(6,i)
        overu2 = srcvals(4,i)**2 + srcvals(5,i)**2 + srcvals(6,i)**2
        overu2 = 1/overu2
        mHcomps(1,1,i) = mHcomps(1,1,i)*overu2
        mHcomps(1,2,i) = mH(1,i)*srcvals(7,i) + mH(2,i)*srcvals(8,i) 
        mHcomps(1,2,i) = mHcomps(1,2,i) + mH(3,i)*srcvals(9,i)
        overv2 = srcvals(7,i)**2 + srcvals(8,i)**2 + srcvals(9,i)**2
        overv2 = 1/overv2
        mHcomps(1,2,i) = mHcomps(1,2,i)*overv2
    enddo

    return
end subroutine axi_surf_harm_field

subroutine rotmat_complex_gmres(a,b,c,s)
    implicit none
    complex *16 a,b,c,s,temp,aphase,bphase

    if (a.eq.0) then
        c = 0
        s = 1
    else
        temp = abs(b/a)**2
        aphase = a/abs(a)
        bphase = b/abs(b)
        c = aphase/sqrt(1.0d0+temp)
        s = bphase/sqrt(1.0d0+1/temp)
    endif
end subroutine rotmat_complex_gmres
