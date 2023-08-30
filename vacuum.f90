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
    real *8 curljre(3,npts),curljim(3,npts),gradrhore(3,npts),gradrhoim(3,npts)

    complex *16 rjvec(3,npts), rho(npts), curlj(3,npts), gradrho(3,npts)
    complex *16 ima
    data ima/(0.0d0,1.0d0)/

    ! real part
    call lpcomp_virtualcasing_addsub(npatches, norders, ixyzs, &
    iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
    eps, nnz, row_ptr, col_ind, iquad, nquad, wnear, real(rjvec), real(rho), &
    novers, nptso, ixyzso, srcover, whtsover, curljre, gradrhore)
    ! imaginary part
    call lpcomp_virtualcasing_addsub(npatches, norders, ixyzs, &
    iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
    eps, nnz, row_ptr, col_ind, iquad, nquad, wnear, aimag(rjvec), aimag(rho), &
    novers, nptso, ixyzso, srcover, whtsover, curljim, gradrhoim)
    ! $OMP PARALLEL DO DEFAULT(SHARED)
    do i = 1,3
        do j = 1,npts
            curlj(i,j) = curljre(i,j) + ima*curljim(i,j)
        enddo
    enddo
    ! $OMP END PARALLEL DO 
    ! $OMP PARALLEL DO DEFAULT(SHARED)
    do i = 1,3
        do j = 1,npts
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
        pot_real_temp = 0
        pot_imag_temp = 0
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

subroutine taylor_vaccuum_solver(npatches, norders, ixyzs, iptype, npts, &
    srccoefs, srcvals, ipars, eps, numit, flux, eps_gmres, niter, errs, rres, &
    sigma, alpha, B, fluxcheck)
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
    ! ==========================================================================
    !
    ! Output:
    !   niter - integer
    !     number of gmres iterations required for relative residual
    !      
    !   errs(1:iter) - relative residual as a function of iteration
    !     number
    !
    !   rres - real *8
    !     relative residual for computed solution
    !          
    !   sigma - copmlex *16(npts)
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
    integer npatches,norder,npols,npts
    integer ifinout
    integer norders(npatches),ixyzs(npatches+1)
    integer iptype(npatches)
    real *8 srccoefs(9,npts),srcvals(12,npts),eps,eps_gmres
    real *8 dpars(2)
    real *8 flux
    complex *16 soln(npts)

    real *8, allocatable :: targs(:,:)
    integer, allocatable :: ipatch_id(:)
    real *8, allocatable :: uvs_targ(:,:)
    integer ndtarg,ntarg

    real *8 errs(numit+1)
    real *8 rres,eps2
    integer niter

    integer nover,npolso,nptso
    integer nnz,nquad
    integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)
    real *8, allocatable :: wnear(:,:)

    real *8, allocatable :: srcover(:,:),wover(:)
    integer, allocatable :: ixyzso(:),novers(:)
    real *8, allocatable :: cms(:,:),rads(:),rad_near(:) 

    integer i,j,jpatch,jquadstart,jstart

    integer ipars(2)
    complex *16 zpars
    real *8 timeinfo(10),t1,t2,omp_get_wtime

    real *8 ttot,done,pi
    real *8 rfac,rfac0
    integer iptype_avg,norder_avg
    integer ikerorder,iquadtype

    ! gmres variables
    real *8 did,rb,wnrm2
    integer numit,it,iind,it1,k,l
    real *8 rmyerr
    complex *16 dtmp,temp
    complex *16 vmat(npts,numit+1),hmat(numit,numit)
    complex *16 cs(numit),sn(numit),svec(numit+1),yvec(numit+1)
    complex *16 wtmp(npts),gradSsigma(3,npts),curlSmH(3,npts)
    complex *16 mH(3,npts),mHcomps(1,2,npts)
    complex *16 jtmp(3,npts),curlSjtemp(3,npts),rhs_gmres(npts)

    ! post-gmres variables
    complex *16 nxgradSD(3,npts), SnxgradSD(3,npts), AcycSnxgradSD
    complex *16 mHnxcurlSmH(3,npts), SmHnxcurlSmH(3,npts), AcycSmHnxcurlSmH
    ! circulation (b-cycle) integral 
    integer m, na, nb 
    real *8, allocatable :: avals(:,:),bvals(:,:)
    real *8, allocatable :: awts(:),bwts(:)
    real *8, allocatable :: auv(:,:),buv(:,:)
    integer, allocatable :: apatches(:),bpatches(:)
    real *8, allocatable :: SnxgradSDrealintp(:,:),SnxgradSDimagintp(:,:)
    real *8, allocatable :: SmHnxcurlSmHrintp(:,:),SmHnxcurlSmHiintp(:,:)

    complex *16 sigma(npts),alpha,B(3,npts),nxB(3,npts),SnxB(3,npts),fluxcheck
    real *8, allocatable :: AcycSnxBrealintp(:,:), AcycSnxBimagintp(:,:)

    complex *16 ima
    data ima/(0.0d0,1.0d0)/

    complex *16 ztmp

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
    call cpu_time(t2)
    call prin2('quadrature generation time=*',t2-t1,1)
    print *, "done generating near quadrature, now starting gmres"

    ! start GMRES 

    ! Find A_12 = i n . curl S_0[m_H] (see [1], section 3.1.4)

    ! compute harmonic surface vector field 
    call axi_surf_harm_field(npts,srcvals,mH,mHcomps)
    ! compute curl S_0[m_H]
    ! This subroutine computes grad S and curl S at the same time, 
    ! so a dummy argument is used for the grad S computation.
    ! $OMP PARALLEL DO DEFAULT(SHARED)
    do i=1,npts
        vmat(i,1) = 0
    enddo
    ! $OMP END PARALLEL DO 
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

    did = -1d0/2
    niter = 0
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
            
        hmat(it,it) = cs(it)*hmat(it,it)+sn(it)*wnrm2
        svec(it1) = -sn(it)*svec(it)
        svec(it) = cs(it)*svec(it)
        rmyerr = abs(svec(it1))/rb
        errs(it) = rmyerr
        print *, "iter=",it,errs(it)
  
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
            rres = 0
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
                rres = rres + abs(did*soln(i) + wtmp(i) - rhs_gmres(i))**2
            enddo
            ! $OMP END PARALLEL DO          
            rres = sqrt(rres)/rb
            niter = it

            ! Now soln = D (see [1]). We solve for alpha and sigma using (31) 
            ! and (32). 
            
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
            ! compute -S[n x grad S[D]]
            call lpcomp_lap_comb_dir_addsub_complex_vec(3,npatches,norders,&
                ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,eps,&
                dpars,nnz,row_ptr,col_ind,iquad,nquad,wnear,nxgradSD,novers,&
                npts,ixyzso,srcover,wover,SnxgradSD)

            ! get A-cycle parametrization
            m = 16 
            na = ipars(2)*m
            nb = ipars(1)*m 
            allocate(avals(9,na),awts(na),auv(2,na),apatches(na))
            allocate(bvals(9,nb),bwts(nb),buv(2,nb),bpatches(nb))
            call get_ab_cycles_torusparam(npatches,norders,ixyzs,iptype,npts,&
                srccoefs,srcvals,ipars,m,na,avals,awts,apatches,auv,nb,bvals,&
                bwts,bpatches,buv)

            ! evaluate -S[n x grad S[D]] on A-cycle
            allocate(SnxgradSDrealintp(3,na),SnxgradSDimagintp(3,na))
            call fun_surf_interp(3,npatches,norders,ixyzs,iptype,npts,&
                real(SnxgradSD),na,apatches,auv,SnxgradSDrealintp)
            call fun_surf_interp(3,npatches,norders,ixyzs,iptype,npts,&
                aimag(SnxgradSD),na,apatches,auv,SnxgradSDimagintp)
            
            ! integrate along A-cycle
            AcycSnxgradSD = 0
            do i = 1,na
                AcycSnxgradSD = AcycSnxgradSD + (&
                    (SnxgradSDrealintp(1,i) + ima*SnxgradSDimagintp(1,i)) * &
                        avals(4,i) + &
                    (SnxgradSDrealintp(2,i) + ima*SnxgradSDimagintp(2,i)) * & 
                        avals(5,i) + & 
                    (SnxgradSDrealintp(3,i) + ima*SnxgradSDimagintp(3,i)) * & 
                        avals(6,i)) * awts(i)
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
            ! compute S[-m_H/2 + curl S[mH]]
            call lpcomp_lap_comb_dir_addsub_complex_vec(3,npatches,norders,&
                ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,eps,&
                dpars,nnz,row_ptr,col_ind,iquad,nquad,wnear,mHnxcurlSmH,&
                novers,npts,ixyzso,srcover,wover,SmHnxcurlSmH)

            ! evaluate S[-m_H/2 + curl S[mH]] on A-cycle
            allocate(SmHnxcurlSmHrintp(3,na),SmHnxcurlSmHiintp(3,na))
            call fun_surf_interp(3,npatches,norders,ixyzs,iptype,npts,&
                real(SmHnxcurlSmH),na,apatches,auv,SmHnxcurlSmHrintp)
            call fun_surf_interp(3,npatches,norders,ixyzs,iptype,npts,&
                aimag(SmHnxcurlSmH),na,apatches,auv,SmHnxcurlSmHiintp)

            ! integrate along A-cycle
            AcycSmHnxcurlSmH = 0
            do i = 1,na
                AcycSmHnxcurlSmH = AcycSmHnxcurlSmH + (&
                    (SmHnxcurlSmHrintp(1,i) + ima*SmHnxcurlSmHiintp(1,i)) * &
                        avals(4,i) + &
                    (SmHnxcurlSmHrintp(2,i) + ima*SmHnxcurlSmHiintp(2,i)) * & 
                        avals(5,i) + & 
                    (SmHnxcurlSmHrintp(3,i) + ima*SmHnxcurlSmHiintp(3,i)) * & 
                        avals(6,i)) * awts(i)
            enddo
            AcycSmHnxcurlSmH = ima*AcycSmHnxcurlSmH

            print *,"A21D + A22 = ",AcycSnxgradSD + AcycSmHnxcurlSmH
            alpha = flux/(AcycSnxgradSD + AcycSmHnxcurlSmH)
            ! $OMP PARALLEL DO DEFAULT(SHARED)
            do i = 1,npts
                sigma(i) = alpha*soln(i)
            enddo
            ! $OMP END PARALLEL DO

            ! compute the magnetic field from sigma and alpha
            ! $OMP PARALLEL DO DEFAULT(SHARED)
            do i = 1,npts
                do j = 1,3
                    B(j,i) = alpha*(-soln(i)*srcvals(j+9,i)/2 + mH(j,i)/2 &
                        - gradSsigma(j,i) + ima*curlSmH(j,i))
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
                dpars,nnz,row_ptr,col_ind,iquad,nquad,wnear,nxB,novers,npts,&
                ixyzso,srcover,wover,SnxB)
            allocate(AcycSnxBrealintp(3,na),AcycSnxBimagintp(3,na))
            call fun_surf_interp(3,npatches,norders,ixyzs,iptype,npts,&
                real(SnxB),na,apatches,auv,AcycSnxBrealintp)
            call fun_surf_interp(3,npatches,norders,ixyzs,iptype,npts,&
                aimag(SnxB),na,apatches,auv,AcycSnxBimagintp)
            fluxcheck = 0
            do i = 1,na
                fluxcheck = fluxcheck + (&
                    (AcycSnxBrealintp(1,i) + ima*AcycSnxBimagintp(1,i))&
                        *avals(4,i) + &
                    (AcycSnxBrealintp(2,i) + ima*AcycSnxBimagintp(2,i))&
                        *avals(5,i) + &
                    (AcycSnxBrealintp(3,i) + ima*AcycSnxBimagintp(3,i))&
                        *avals(6,i)) * awts(i)
            enddo
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

! subroutine test_axi_surf_harm_field(npatches, norders, ixyzs, iptype, npts, &
!     srccoefs, srcvals, f, surfdivf, sum)
!     ! 
!     ! Tests that the vector field generated by axi_surf_harm_field() satisfies 
!     ! the conditions for harmonicity. 
!     ! 
!     ! ==========================================================================
!     ! 
!     ! Input: (TODO)
!     !   npts - integer
!     !      total number of discretiztion points on the boundary
!     !   
!     !   srcvals - real *8 (12,npts)
!     !     xyz(u,v) and derivative info sampled at the 
!     !     discretization nodes on the surface
!     !     srcvals(1:3,i) - xyz info
!     !     srcvals(4:6,i) - dxyz/du info
!     !     srcvals(7:9,i) - dxyz/dv info 
!     !     srcvals(10:12,i) - nxyz info 
!     ! 
!     !   field - real *8 (3,npts)
!     !     field generated by axi_surf_harm_field()
!     ! 
!     implicit none
!     integer, intent(in) :: npatches, norders(npatches), ixyzs(npatches+1)
!     integer, intent(in) :: iptype(npatches), npts
!     real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts), f(3,npts)
!     real *8, intent(out) :: surfdivf(npts), sum

!     integer :: i, istart, npols
!     real *8 :: dfuv(3,2,npts), ginv(2,2,npts), temp(6,npts)
!     ! integer i, j
!     ! real *8 gdet
!     ! real *8 umet(3,npts), vmet(3,npts)

!     ! compute metric quantities, see Dan's paper
!     ! do i = 1,npts
!     !     ginv(1,1,i) = 0
!     !     do j = 1,npts
!     !         ginv(1,1,i) = ginv(1,1,i) + srcvals(6+j,i)**2
!     !     enddo
!     !     ginv(1,2,i) = 0
!     !     do j = 1,npts
!     !         ginv(1,2,i) = ginv(1,2,i) - srcvals(3+j,i)*srcvals(6+j,i)
!     !     enddo
!     !     ginv(2,1,i) = 0
!     !     do j = 1,npts
!     !         ginv(2,1,i) = ginv(2,1,i) - srcvals(6+j,i)*srcvals(3+j,i)
!     !     enddo
!     !     ginv(2,2,i) = 0
!     !     do j = 1,npts
!     !         ginv(2,2,i) = ginv(2,2,i) + srcvals(3+j,i)**2
!     !     enddo
!     ! enddo
!     ! do i = 1,npts
!     !     gdet = abs(ginv(1,1,i)*ginv(2,2,i) - ginv(1,2,i)*ginv(2,1,i))
!     !     ginv(1,1,i) = ginv(1,1,i)/gdet
!     !     ginv(1,2,i) = ginv(1,2,i)/gdet
!     !     ginv(2,1,i) = ginv(2,1,i)/gdet
!     !     ginv(2,2,i) = ginv(2,2,i)/gdet
!     !     do j = 1,npts
!     !         umet(j,i) = ginv(1,1,i)*srcvals(3+j,i) + ginv(2,1,i)*srcvals(6+j,i)
!     !         vmet(j,i) = ginv(1,2,i)*srcvals(3+j,i) + ginv(2,2,i)*srcvals(6+j,i)
!     !     enddo
!     ! enddo
!     ! NEED DERIVATIVES OF m_H

!     call get_inv_first_fundamental_form(npatches, norders, ixyzs, iptype, &
!     npts, srccoefs, srcvals, ginv)
!     print *,"dfuv before"
!     print *,dfuv(:,:,npts-1) !write(13,*) dfuv(:,:,1)
!     !write(13,*) "before"
!     do i = 1,npatches
!         istart = ixyzs(i)
!         npols = ixyzs(i+1)-ixyzs(i)
!         call get_surf_uv_grad_guru(3, norders(i), npols, &
!         iptype(i), f(1,istart), dfuv(1,1,istart))
!     enddo
!     ! call get_surf_uv_grad(3, npatches, norders, ixyzs, iptype, npts, srccoefs, &
!     ! srcvals, f, dfuv)
!     print *,"dfuv after"
!     print *,dfuv(:,:,npts-1) !write(13,*) dfuv(:,:,1)
!     !write(13,*) "after"
!     do i = 1,npts
!         temp(1,i) = ginv(1,1,i)*srcvals(4,i) + ginv(1,2,i)*srcvals(7,i)
!         temp(2,i) = ginv(1,1,i)*srcvals(5,i) + ginv(1,2,i)*srcvals(8,i)
!         temp(3,i) = ginv(1,1,i)*srcvals(6,i) + ginv(1,2,i)*srcvals(9,i)
!         temp(4,i) = ginv(2,1,i)*srcvals(4,i) + ginv(2,2,i)*srcvals(7,i)
!         temp(5,i) = ginv(2,1,i)*srcvals(5,i) + ginv(2,2,i)*srcvals(8,i)
!         temp(6,i) = ginv(2,1,i)*srcvals(6,i) + ginv(2,2,i)*srcvals(9,i)
!     enddo
!     do i = 1,npts
!         surfdivf(i) = temp(1,i)*dfuv(1,1,i)
!         surfdivf(i) = surfdivf(i) + temp(2,i)*dfuv(2,1,i)
!         surfdivf(i) = surfdivf(i) + temp(3,i)*dfuv(3,1,i)
!         surfdivf(i) = surfdivf(i) + temp(4,i)*dfuv(1,2,i)
!         surfdivf(i) = surfdivf(i) + temp(5,i)*dfuv(2,2,i)
!         surfdivf(i) = surfdivf(i) + temp(6,i)*dfuv(3,2,i)
!         ! if (i.eq.26) then
!         !     print *,dfuv(3,2,i)
!         !     ! print *,surfdivf(i)
!         ! endif
!     enddo
!     ! testing
!     ! do i = 1,npts
!     !     surfdivf(i) = ginv(2,2,i)
!     ! enddo
!     ! end testing
!     sum = 0
!     do i = 1,npts
!         sum = sum + abs(surfdivf(i))
!     enddo

!     return
! end subroutine test_axi_surf_harm_field

subroutine setup_geom_vacuum(igeomtype, norder, npatches, ipars, srcvals, &
    srccoefs, ifplot, fname)
    ! 
    ! Geometry setup function lifted from superconductor-type1/test/
    ! test_ab_cycle.f

    implicit none
    integer npols,itype,npmax,ntri,nover
    real *8 done,pi,umin,umax,vmin,vmax
    integer igeomtype,norder,npatches,ipars(*),ifplot
    character (len=*) fname
    real *8 srcvals(12,*), srccoefs(9,*)
    real *8, allocatable :: uvs(:,:),umatr(:,:),vmatr(:,:),wts(:)

    real *8, pointer :: ptr1,ptr2,ptr3,ptr4
    integer, pointer :: iptr1,iptr2,iptr3,iptr4
    real *8, target :: p1(10),p2(10),p3(10),p4(10)
    real *8, allocatable, target :: triaskel(:,:,:)
    real *8, allocatable, target :: deltas(:,:)
    integer, allocatable :: isides(:)
    integer, target :: nmax,mmax

    procedure (), pointer :: xtri_geometry

    external xtri_stell_eval,xtri_sphere_eval,xtri_wtorus_eval
      
    npols = (norder+1)*(norder+2)/2
    allocate(uvs(2,npols),umatr(npols,npols),vmatr(npols,npols))
    allocate(wts(npols))

    call vioreanu_simplex_quad(norder,npols,uvs,umatr,vmatr,wts)

    if(igeomtype.eq.1) then
        itype = 2
        allocate(triaskel(3,3,npatches))
        allocate(isides(npatches))
        npmax = npatches
        ntri = 0
        call xtri_platonic(itype, ipars(1), npmax, ntri, triaskel, isides)

        xtri_geometry => xtri_sphere_eval
        ptr1 => triaskel(1,1,1)
        ptr2 => p2(1)
        ptr3 => p3(1)
        ptr4 => p4(1)


        if(ifplot.eq.1) then
            call xtri_vtk_surf(fname,npatches,xtri_geometry,ptr1,ptr2,ptr3,&
                ptr4, norder,'Triangulated surface of the sphere')
        endif


        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,ptr3,ptr4,npols,uvs,&
            umatr,srcvals,srccoefs)
      endif

    if(igeomtype.eq.2) then
        done = 1
        pi = atan(done)*4
        umin = 0
        umax = 2*pi
        vmin = 2*pi
        vmax = 0
        allocate(triaskel(3,3,npatches))
        nover = 0
        call xtri_rectmesh_ani(umin,umax,vmin,vmax,ipars(1),ipars(2),nover,&
            npatches,npatches,triaskel)

        mmax = 2
        nmax = 1
        xtri_geometry => xtri_stell_eval

        allocate(deltas(-1:mmax,-1:nmax))
        deltas(-1,-1) = 0.17d0
        deltas(0,-1) = 0
        deltas(1,-1) = 0
        deltas(2,-1) = 0

        deltas(-1,0) = 0.11d0
        deltas(0,0) = 1
        deltas(1,0) = 4.5d0
        deltas(2,0) = -0.25d0

        deltas(-1,1) = 0
        deltas(0,1) = 0.07d0
        deltas(1,1) = 0
        deltas(2,1) = -0.45d0

        ptr1 => triaskel(1,1,1)
        ptr2 => deltas(-1,-1)
        iptr3 => mmax
        iptr4 => nmax

        if(ifplot.eq.1) then
            call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2,iptr3,&
            iptr4,norder,'Triangulated surface of the stellarator')
        endif

        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,iptr3,iptr4,npols,&
            uvs,umatr,srcvals,srccoefs)
    endif

    if(igeomtype.eq.3) then
        done = 1
        pi = atan(done)*4
        umin = 0
        umax = 2*pi
        vmin = 0
        vmax = 2*pi
        allocate(triaskel(3,3,npatches))
        nover = 0
        call xtri_rectmesh_ani(umin,umax,vmin,vmax,ipars(1),ipars(2),nover,&
            npatches,npatches,triaskel)
        call prinf('npatches=*',npatches,1)
         
        p1(1) = 1
        p1(2) = 2
        p1(3) = 0.25d0

        p2(1) = 1.0d0
        p2(2) = 1.0d0
        p2(3) = 1.0d0

        ! number of oscillations
        p4(1) = 3.0d0


        ptr1 => triaskel(1,1,1)
        ptr2 => p1(1)
        ptr3 => p2(1)
        ptr4 => p4(1)
        xtri_geometry => xtri_wtorus_eval
        if(ifplot.eq.1) then
            call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2,ptr3,&
                ptr4, norder,'Triangulated surface of the wtorus')
        endif


        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,ptr3,ptr4,npols,uvs,&
            umatr,srcvals,srccoefs)
    endif
      
    if(igeomtype.eq.4) then
        done = 1
        pi = atan(done)*4
        umin = 0
        umax = 2*pi
        vmin = 0
        vmax = 2*pi
        allocate(triaskel(3,3,npatches))
        nover = 0
        call xtri_rectmesh_ani(umin,umax,vmin,vmax,ipars(1),ipars(2),nover,&
            npatches,npatches,triaskel)
        call prinf('npatches=*',npatches,1)
         
        p1(1) = 1.0d0
        p1(2) = 1.75d0
        p1(3) = 0.25d0

        p2(1) = 1.0d0
        p2(2) = 1.0d0
        p2(3) = 1.0d0
        ! number of oscillations
        p4(1) = 0.0d0


        ptr1 => triaskel(1,1,1)
        ptr2 => p1(1)
        ptr3 => p2(1)
        ptr4 => p4(1)
        xtri_geometry => xtri_wtorus_eval
        if(ifplot.eq.1) then
            call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2,ptr3,&
                ptr4, norder,'Triangulated surface of the torus')
        endif


        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,ptr3,ptr4,npols,uvs,&
            umatr,srcvals,srccoefs)
    endif

    if(igeomtype.eq.5) then
        done = 1
        pi = atan(done)*4
        umin = 0
        umax = 2*pi
        vmin = 2*pi
        vmax = 0
        allocate(triaskel(3,3,npatches))
        nover = 0
        call xtri_rectmesh_ani(umin,umax,vmin,vmax,ipars(1),ipars(2),nover,&
            npatches,npatches,triaskel)

        mmax = 2
        nmax = 1
        xtri_geometry => xtri_stell_eval

        allocate(deltas(-1:mmax,-1:nmax))
        deltas(-1,-1) = 0.17d0*0
        deltas(0,-1) = 0
        deltas(1,-1) = 0
        deltas(2,-1) = 0

        deltas(-1,0) = 0.11d0*0
        deltas(0,0) = 1*0
        deltas(1,0) = 2.0d0
        deltas(2,0) = -0.25d0*0

        deltas(-1,1) = 0
        deltas(0,1) = 1.0d0
        deltas(1,1) = 0
        deltas(2,1) = -0.45d0*0


        ptr1 => triaskel(1,1,1)
        ptr2 => deltas(-1,-1)
        iptr3 => mmax
        iptr4 => nmax

        if(ifplot.eq.1) then
            call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2,iptr3,&
                iptr4, norder,'Triangulated surface of the stellarator')
        endif

        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,iptr3,iptr4,npols,&
            uvs,umatr,srcvals,srccoefs)
    endif
      
    return  
end subroutine setup_geom_vacuum

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

program vacuum
    implicit none

    integer :: norder,npols,npatches,npts,m,na,nb,ipars(2)
    integer, allocatable :: norders(:),ixyzs(:),iptype(:)
    real *8, allocatable :: srccoefs(:,:),srcvals(:,:),surfdivf(:)
    real *8, allocatable :: mHrealinterp(:,:),mHimaginterp(:,:)
    real *8, allocatable :: avals(:,:),bvals(:,:),awts(:),bwts(:)
    real *8, allocatable :: auv(:,:),buv(:,:)
    integer, allocatable :: apatches(:),bpatches(:)
    complex *16 :: sum
    complex *16, allocatable :: mH(:,:),mHcomps(:,:,:),nxmH(:,:),dmH(:,:)
    real *8, allocatable :: ffform(:,:,:),dmHreal(:,:),dmHimag(:,:)
    integer i,igeomtype,ifplot,j
    character *300 fname

    ! add'l variables for taylor_vacuum_solver()
    real *8 eps,flux,eps_gmres,rres
    integer numit,niter
    real *8, allocatable :: errs(:)
    complex *16, allocatable :: sigma(:),B(:,:),Bdotn(:)
    real *8, allocatable :: AcycBrealintp(:,:), AcycBimagintp(:,:)
    complex *16 alpha, fluxcheck

    complex *16 ima
    data ima/(0.0d0,1.0d0)/

    ! if the number of sample points is too low, get a seg fault
    igeomtype = 5 ! set to 4 for axisymmetric geometry
    norder = 8 
    npols = (norder+1)*(norder+2)/2
    ifplot = 0
    fname = 'torus.vtk'

    ipars(1) = 6
    if(igeomtype.eq.2) ipars(1) = 10
    ipars(2) = ipars(1)*3
    npatches = 2*ipars(1)*ipars(2)
    
    npts = npols*npatches
    print *,"npts = ",npts

    allocate(norders(npatches),ixyzs(npatches+1),iptype(npatches))
    do i = 1,npatches
        norders(i) = norder
        ixyzs(i) = 1 + (i-1)*npols 
        iptype(i) = 1
    enddo
    ixyzs(npatches+1) = npts+1

    allocate(srccoefs(9,npts),srcvals(12,npts),surfdivf(npts))

    call setup_geom_vacuum(igeomtype,norder,npatches,ipars,srcvals,srccoefs,&
        ifplot,fname)
    
    allocate(mH(3,npts),mHcomps(1,2,npts),nxmH(3,npts))
    call axi_surf_harm_field(npts,srcvals,mH,mHcomps)
    ! open(13, file="torusmH.txt")
    ! do i = 1,npts
    !     do j = 1,3
    !         write(13,*) mH(j,i)
    !     enddo
    ! enddo
    ! close(13)

    allocate(dmH(1,npts))
    allocate(ffform(2,2,npts),dmHreal(1,npts),dmHimag(1,npts))
    call get_first_fundamental_form(npatches,norders,ixyzs,iptype,npts,&
        srccoefs,srcvals,ffform)
    print *,ffform(1,1,2),ffform(1,2,2)
    print *,ffform(2,1,2),ffform(2,2,2)
    call get_surf_div(1,npatches,norders,ixyzs,iptype,npts,srccoefs,srcvals,&
        real(mHcomps),dmHreal)
    call get_surf_div(1,npatches,norders,ixyzs,iptype,npts,srccoefs,srcvals,&
        aimag(mHcomps),dmHimag)
    sum = 0
    do i = 1,npts
        dmH(1,i) = dmHreal(1,i) + ima*dmHimag(1,i)
        sum = sum + abs(dmH(1,i))
    enddo
    print *,"dmH",sum

    ! do i = 1,npts
    !     nxmH(1,i) = srcvals(11,i)*mH(3,i) - srcvals(12,i)*mH(2,i) 
    !     nxmH(2,i) = srcvals(12,i)*mH(1,i) - srcvals(10,i)*mH(3,i)
    !     nxmH(3,i) = srcvals(10,i)*mH(2,i) - srcvals(11,i)*mH(1,i)
    ! enddo
    ! open(40,file="nxmHplusimH.txt")
    ! do i = 1,npts
    !     do j = 1,3
    !         write(40,*) nxmH(j,i) + ima*mH(j,i)
    !     enddo
    ! enddo
    ! close(40)

    ! This is mH = phihat, used for testing. 
    ! do i = 1,npts
    !     mH(1,i) = -srcvals(2,i)
    !     mH(2,i) = srcvals(1,i)
    !     mH(3,i) = 0
    ! enddo

    ! ===================================== 
    
    ! m = 16
    ! na = ipars(2)*m
    ! nb = ipars(1)*m
    ! allocate(avals(9,na),awts(na),auv(2,na),apatches(na))
    ! allocate(auv(2,na)) ! why does gmres blow up if this is commented out?
    ! allocate(bvals(9,nb),bwts(nb),buv(2,nb),bpatches(nb))
    ! call get_ab_cycles_torusparam(npatches,norders,ixyzs,iptype,npts,srccoefs,&
    !     srcvals,ipars,m,na,avals,awts,apatches,auv,nb,bvals,bwts,bpatches,buv)
    ! open(41,file="torusavals.txt")
    ! write(41,*) avals
    ! close(41)
    
    ! allocate(mHrealinterp(3,nb),mHimaginterp(3,nb))    
    ! call fun_surf_interp(3,npatches,norders,ixyzs,iptype,npts,real(mH),nb,&
    !     bpatches,buv,mHrealinterp)
    ! call fun_surf_interp(3,npatches,norders,ixyzs,iptype,npts,aimag(mH),nb,&
    !     bpatches,buv,mHimaginterp)
    ! open(14,file="torusmHrealinterpbcycle.txt")
    ! write(14,*) mHrealinterp
    ! close(14)
    ! open(15,file="torusmHimaginterpbcycle.txt")
    ! write(15,*) mHimaginterp
    ! close(15)

    ! =====================================

    ! allocate(mHrealinterp(3,nb),mHimaginterp(3,nb))    
    ! call fun_surf_interp(3,npatches,norders,ixyzs,iptype,npts,real(mH),nb,&
    !     bpatches,buv,mHrealinterp)
    ! call fun_surf_interp(3,npatches,norders,ixyzs,iptype,npts,aimag(mH),nb,&
    !     bpatches,buv,mHimaginterp)
    ! open(14,file="mHrealinterpbcycle.txt")
    ! write(14,*) mHrealinterp
    ! close(14)
    ! open(15,file="mHimaginterpbcycle.txt")
    ! write(15,*) mHimaginterp
    ! close(15)

    ! call test_axi_surf_harm_field(npatches,norders,ixyzs,iptype,npts,srccoefs,&
    ! srcvals,real(mH),surfdivf,sum)
    ! open(16, file='surfdivf.txt')
    ! do i = 1,npts
    !     write(16,*) surfdivf(i)
    ! enddo
    ! close(16)

    ! ===========================================

    eps = 0.51d-4
    eps_gmres = 1d-8
    flux = 0.3d0
    numit = 200
    niter = 0

    allocate(errs(numit+1))
    allocate(sigma(npts),B(3,npts))

    call taylor_vaccuum_solver(npatches,norders,ixyzs,iptype,npts,&
        srccoefs,srcvals,ipars,eps,numit,flux,eps_gmres,niter,&
        errs,rres,sigma,alpha,B,fluxcheck)

    open(16,file="sigma.txt")
    do i = 1,npts
        write(16,*) sigma(i)
    enddo
    close(16)

    open(17,file="alpha.txt")
    write(17,*) alpha
    close(17)

    open(18,file="B.txt")
    do i = 1,npts
        do j = 1,3
            write(18,*) B(j,i)
        enddo
    enddo
    close(18)

    ! check that B . n = 0
    allocate(Bdotn(npts))
    do i = 1,npts
        Bdotn(i) = 0
        do j = 1,3
            Bdotn(i) = Bdotn(i) + B(j,i)*srcvals(j+9,i)
        enddo 
    enddo
    open(19,file="Bdotn.txt")
    do i = 1,npts
        write(19,*) Bdotn(i)
    enddo
    close(19)

    ! check flux
    print *,"A-cyc integral:",fluxcheck
end program vacuum 