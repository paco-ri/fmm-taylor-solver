subroutine Sprime_matrix(npatches, norders, ixyzs, iptype, npts, &
    srccoefs, srcvals, ipars, Sprime)
    ! 
    ! Computes the matrix associated with -n . S_0 (negative sign included!)
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
    ! ==========================================================================
    !
    ! Output:
    !
    !
    implicit none
    integer npatches,npts
    integer norders(npatches),ixyzs(npatches+1)
    integer iptype(npatches)
    real *8 srccoefs(9,npts),srcvals(12,npts),eps
    real *8 dpars(2)

    real *8, allocatable :: targs(:,:)
    integer, allocatable :: ipatch_id(:)
    real *8, allocatable :: uvs_targ(:,:)
    integer ndtarg,ntarg

    integer nptso
    integer nnz,nquad
    integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)
    real *8, allocatable :: wnear(:,:)

    real *8, allocatable :: srcover(:,:),wover(:)
    integer, allocatable :: ixyzso(:),novers(:)
    real *8, allocatable :: cms(:,:),rads(:),rad_near(:) 

    integer i,j

    integer ipars(2)
    real *8 t1,t2

    real *8 done,pi
    real *8 rfac,rfac0
    integer iptype_avg,norder_avg
    integer ikerorder,iquadtype

    ! operator evaluation variables
    real *8 did
    complex *16 w(npts),unitvec(npts)
    complex *16 wtmp(npts),gradSsigma(3,npts)
    complex *16 jtmp(3,npts),curlSjtemp(3,npts)

    ! output
    complex *16 Sprime(npts,npts)

    ! circulation (A-cycle) integral 
    integer m, na, nb 
    real *8, allocatable :: avals(:,:),bvals(:,:)
    real *8, allocatable :: awts(:),bwts(:)
    real *8, allocatable :: auv(:,:),buv(:,:)
    integer, allocatable :: apatches(:),bpatches(:)

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
    print *, "starting to generate/load near quadrature"
    call cpu_time(t1)
    ! trying to read wnear from file 11/6/23
    ! call getnearquad_magnetostatics(npatches,norders,ixyzs, &
    ! iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,ipatch_id, &
    ! uvs_targ,eps,iquadtype,nnz,row_ptr,col_ind,iquad,rfac0,nquad, &
    ! wnear)
    open(50,form='unformatted',file='wnear_fat.dat',action='read')
    read(50) wnear
    close(50)
    call cpu_time(t2)
    print *, 'quadrature generation time=',t2-t1
    print *, "done generating near quadrature"

    ! open(50,file='wnear.dat',form='unformatted')
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

    did = -1d0/2

    ! $OMP PARALLEL DO DEFAULT(SHARED)
    do i=1,npts
        w(i) = 0
    enddo
    ! $OMP END PARALLEL DO

    ! dummy argument (current) for layer potential evaluation
    ! $OMP PARALLEL DO DEFAULT(SHARED)
    do i=1,3
        do j=1,npts
            jtmp(i,j) = 0
            curlSjtemp(i,j) = 0
        enddo
    enddo
    ! $OMP PARALLEL DO DEFAULT(SHARED)

    ! layer potential evaluation
    do j = 1,npts
        print *,'iter',j
        do i = 1,npts
            unitvec(i) = 0
        enddo
        unitvec(j) = 1.0d0
        call lpcomp_virtualcasing_addsub_complex(npatches,norders,ixyzs,&
        iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,eps,nnz,&
        row_ptr,col_ind,iquad,nquad,wnear,jtmp,unitvec,&
        novers,nptso,ixyzso,srcover,wover,curlSjtemp,gradSsigma)

        ! $OMP PARALLEL DO DEFAULT(SHARED)
        do i=1,npts
            wtmp(i) = -gradSsigma(1,i)*srcvals(10,i)
            wtmp(i) = wtmp(i) - gradSsigma(2,i)*srcvals(11,i)
            wtmp(i) = wtmp(i) - gradSsigma(3,i)*srcvals(12,i)
        enddo
        ! $OMP END PARALLEL DO

        do i = 1,npts
            Sprime(i,j) = wtmp(i)
        enddo
    enddo
end subroutine Sprime_matrix

program Sprimematrix
    implicit none

    integer norder,npols,npatches,npts,ipars(2)
    integer, allocatable :: norders(:),ixyzs(:),iptype(:)
    real *8, allocatable :: srccoefs(:,:),srcvals(:,:),surfdivf(:)
    complex *16, allocatable :: mH(:,:),mHcomps(:,:,:),nxmH(:,:)
    integer i,igeomtype,ifplot,j
    character *300 fname

    ! testing 
    real *8 pi,rmin,rmaj
    integer Nphi,Ntheta,solvetest
    real *8, allocatable :: targ(:,:)
    complex *16, allocatable :: B0(:,:),B0dotn(:)

    ! Sprime
    complex *16, allocatable :: Sprime(:,:)

    complex *16 ima
    data ima/(0.0d0,1.0d0)/

    ! if the number of sample points is too low, get a seg fault
    igeomtype = 4 ! set to 4 for axisymmetric geometry
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

    ! Generate field for testing, B0
    pi = 4*atan(1.0d0)
    rmaj = 2.0d0 ! 10.0d0 ! 1.75d0
    rmin = 1.25d0 ! 2.0d0 ! 1.5d0
    Nphi = 128
    Ntheta = 128 ! for B0 from torus
    ! Ntheta = 2**14
    allocate(targ(3,npts))
    allocate(B0(3,npts))
    do i=1,npts
        do j=1,3
            targ(j,i) = srcvals(j,i)
        enddo
    enddo
    
    call curl_vector_potential_torus(Nphi,Ntheta,rmin,rmaj,1.0d0,npts,targ,B0)

    ! get B0 . n
    allocate(B0dotn(npts))
    do i = 1,npts
        B0dotn(i) = B0(1,i)*srcvals(10,i) + B0(2,i)*srcvals(11,i) + &
            B0(3,i)*srcvals(12,i)
    enddo

    allocate(Sprime(npts,npts))
    call Sprime_matrix(npatches,norders,ixyzs,iptype,npts,&
        srccoefs,srcvals,ipars,Sprime)

    open(21,form='unformatted',file='Sprime.dat')
    write(21) Sprime
    close(21)
end program Sprimematrix