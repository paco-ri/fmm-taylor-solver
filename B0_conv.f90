program B0_conv
    implicit none
    integer igeomtype,norder,npols,npatches,ipars(2),ifplot,i,k
    real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
    character *300 :: fname

    integer Nphi,Ntheta,npts
    real *8 rmin,rmaj,jmag
    real *8, allocatable :: targ(:,:)
    complex *16, allocatable :: B01(:,:),B02(:,:)

    real *8 pi

    pi = 4*atan(1.0d0)

    ! set up geometry

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

    allocate(srcvals(12,npts),srccoefs(9,npts),targ(3,npts))
    call setup_geom_vacuum(igeomtype,norder,npatches,ipars,srcvals,srccoefs,&
        ifplot,fname)
    do i = 1,npts
        do k = 1,3
            targ(k,i) = srcvals(k,i)
        enddo
    enddo

    ! set up B0

    rmaj = 10.0d0 ! 1.75d0
    rmin = 2.0d0 ! 1.5d0
    Nphi = 128
    Ntheta = 128
    jmag = 1.0d0

    allocate(B01(3,npts),B02(3,npts))
    call curl_vector_potential_torus(Nphi,Ntheta,rmin,rmaj,jmag,npts,targ,B01)
    print *,real(B01(1:3,1:4))
    ! Nphi = 256
    ! Ntheta = 256
    ! call curl_vector_potential_torus(Nphi,Ntheta,rmin,rmaj,jmag,npts,targ,B02)
    ! print *,real(B02(1:3,1:4))
    ! print *,maxval(real(B01-B02))

    open(22,file="B0ampere.txt")
    do i = 1,npts
        write(22,*) 20.28425343*srcvals(2,i)/(srcvals(1,i)**2+srcvals(2,i)**2)
        write(22,*) -20.28425343*srcvals(1,i)/(srcvals(1,i)**2+srcvals(2,i)**2)
        write(22,*) 0
    enddo
    close(22)
end program B0_conv