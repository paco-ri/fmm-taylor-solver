program vacuum_example
    implicit none

    integer norder,npols,npatches,npts,ipars(2)
    integer, allocatable :: norders(:),ixyzs(:),iptype(:)
    real *8, allocatable :: srccoefs(:,:),srcvals(:,:),surfdivf(:)
    complex *16, allocatable :: mH(:,:),mHcomps(:,:,:),nxmH(:,:)
    integer i,igeomtype,ifplot,j
    character *300 fname

    ! add'l variables for taylor_vacuum_solver()
    real *8 eps,flux,eps_gmres,rres1,rres2
    integer numit,niter1,niter2
    real *8, allocatable :: errs1(:),errs2(:)
    complex *16, allocatable :: sigma(:),B(:,:),Bdotn(:)
    complex *16 alpha, fluxcheck
    complex *16, allocatable :: Bint(:,:)

    ! testing 
    real *8 pi,rmin,rmaj
    integer Nphi,Ntheta,solvetest
    real *8, allocatable :: targ(:,:)
    complex *16, allocatable :: B0(:,:),B0dotn(:)
    integer nintpts
    real *8, allocatable :: intpts(:,:)
    complex *16 B0surfint

    ! outfile naming
    character(len=2) nor
    character(len=8) npa
    character(len=16) tag

    complex *16 ima
    data ima/(0.0d0,1.0d0)/

    igeomtype = 4 ! set to 4 for toroidal geometry
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

    write(nor,'(I2.2)') norder
    write(npa,'(I8.8)') npatches
    tag = 'nor'//nor//'npa'//npa

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

    ! setting up GMRES
    eps = 1d-8 ! 0.51d-4
    eps_gmres = 1d-8
    flux = 0.3d0
    numit = 200
    niter1 = 0
    niter2 = 0

    allocate(errs1(numit+1),errs2(numit+1))
    allocate(sigma(npts),B(3,npts))

    ! Generate field for testing, B0
    pi = 4*atan(1.0d0)
    rmaj = 2.0d0 ! 10.0d0 ! 1.75d0
    rmin = 1.25d0 ! 2.0d0 ! 1.5d0
    Nphi = 128  

    ! Ntheta = 128 ! for B0 from torus
    Ntheta = 2**14 ! for B0 from loop
    allocate(targ(3,npts))
    allocate(B0(3,npts))
    do i=1,npts
        do j=1,3
            targ(j,i) = srcvals(j,i)
        enddo
    enddo
    
    ! call curl_vector_potential_torus(Nphi,Ntheta,rmin,rmaj,1.0d0,npts,targ,B0)
    call curl_vector_potential_loop(Ntheta,rmin,rmaj,1.0d0,npts,targ,B0)
    ! do i = 1,npts
    !     B0(3,i) = 0
    ! enddo
    open(15,file='B0'//tag//'.txt')
    do i = 1,npts
        do j = 1,3
            write(15,*) B0(j,i)
        enddo
    enddo
    close(15)

    ! get B0 . n
    allocate(B0dotn(npts))
    do i = 1,npts
        B0dotn(i) = B0(1,i)*srcvals(10,i) + B0(2,i)*srcvals(11,i) + &
            B0(3,i)*srcvals(12,i)
    enddo

    ! generate points interior to Omega
    nintpts = 100
    allocate(intpts(3,nintpts))
    do i = 1,nintpts
        intpts(1,i) = rmaj*cos(2*pi*(i-1)/nintpts)
        intpts(2,i) = rmaj*sin(2*pi*(i-1)/nintpts)
        intpts(3,i) = 0
    enddo
    allocate(Bint(3,nintpts))

    solvetest = 1
    call taylor_vaccuum_solver(npatches,norders,ixyzs,iptype,npts,&
        srccoefs,srcvals,ipars,eps,numit,flux,eps_gmres,solvetest,B0,&
        niter1,niter2,errs1,errs2,rres1,rres2,sigma,alpha)
    print *,"flux=",flux

    ! open(16,file="sigma.txt")
    ! do i = 1,npts
    !     write(16,*) sigma(i)
    ! enddo
    ! close(16)

    ! open(17,file="alpha.txt")
    ! write(17,*) alpha
    ! close(17)

    open(18,file='B'//tag//'.txt')
    do i = 1,npts
        do j = 1,3
            write(18,*) B(j,i)
        enddo
    enddo
    close(18)

    ! check that B . n = B0 . n
    allocate(Bdotn(npts))
    do i = 1,npts
        Bdotn(i) = 0
        do j = 1,3
            Bdotn(i) = Bdotn(i) + B(j,i)*srcvals(j+9,i)
        enddo 
    enddo
    open(19,file='Bdotn'//tag//'.txt')
    do i = 1,npts
        write(19,*) Bdotn(i)
    enddo
    close(19)
    open(23,file='B0dotn'//tag//'.txt')
    do i = 1,npts
        write(23,*) B0dotn(i)
    enddo
    close(23)

    ! integrate B0dotn over boundary of Omega
    B0surfint = 0
    do i = 1,npts
        B0surfint = B0surfint + B0dotn(i)*(srcvals(1,i)**2 + srcvals(2,i)**2)
    enddo
    B0surfint= B0surfint*rmin/(4*pi**2*npts)
    print *,'B0 . n surface integral', B0surfint

    open(20,file='Bdotnerr'//tag//'.txt')
    do i = 1,npts
        write(20,*) Bdotn(i)-B0dotn(i)
    enddo
    close(20)

    open(21,file='Berrrel'//tag//'.txt')
    do i = 1,npts
        do j = 1,3
            write(21,*) abs(B(j,i)-B0(j,i))/abs(B0(j,i))
        enddo
    enddo
    close(21)

    ! 20.28425343*
    ! open(22,file="B0ampere.txt")
    ! do i = 1,npts
    !     write(22,*) 127.45268097113586*srcvals(2,i)/(2*pi*(srcvals(1,i)**2+srcvals(2,i)**2))
    !     write(22,*) -127.45268097113586*srcvals(1,i)/(2*pi*(srcvals(1,i)**2+srcvals(2,i)**2))
    !     write(22,*) 0
    ! enddo
    ! close(22)

    open(23,file='Berrabs'//tag//'.txt')
    do i = 1,npts
        do j = 1,3
            write(23,*) abs(B(j,i)-B0(j,i))
        enddo
    enddo
    close(23)

    ! open(24,file='B0intanalytic.txt')
    ! do i = 1,nintpts
    !     write(24,*) 20.28425343*intpts(2,i)/(2*pi*rmaj**2)
    !     write(24,*) 20.28425343*intpts(1,i)/(2*pi*rmaj**2)
    !     write(24,*) 0 
    ! enddo
    ! close(24)

    ! open(25,file='Bint.txt')
    ! do i = 1,nintpts
    !     do j = 1,3
    !         write(25,*) Bint(j,i)
    !     enddo
    ! enddo
    ! close(25)

    ! check flux
    print *,"A-cyc integral:",fluxcheck
end program vacuum_example
