subroutine curl_vector_potential(nsource, source, j, ntarg, targ, curlA)
    !
    ! Computes the curl of the vector potential A, where
    ! 
    !         /\
    !         |       
    !     A = |   g_0(x,x') J(x') ds(x'), 
    !         |    
    !        \/ S
    ! 
    ! J is a current on an axisymmetric surface S that surrounds the solution 
    ! domain, and g_0 is the Laplace kernel |x-x'|^-1. 
    ! 
    ! ==========================================================================
    ! 
    ! Input:
    !   nsource - integer
    !     total number of discretization points on S
    ! 
    !   source - real *8 (3,nsource)
    !     discretization nodes on S
    ! 
    !   j - complex *16 (3,nsource)
    !     current on S
    ! 
    !   ntarg - integer
    !     total number of discretization points on solution domain boundary
    ! 
    !   targ - real *8 (3,ntarg)
    !     discretization nodes on solution domain boundary
    ! 
    ! Output:
    !   curlA - complex *16 (3,ntarg)
    ! 
    implicit none 
    integer nsource,ntarg
    real *8 source(3,nsource),targ(3,ntarg)
    complex *16 j(3,nsource),curlA(3,ntarg)

    integer i,k
    real *8 gradlapker(3,nsource),targpt(3),temp,pi

    pi = 4*atan(1.0d0)

    do i = 1,ntarg
        do k = 1,3
            curlA = (0d0,0d0)
        enddo
    enddo
    do i = 1,ntarg
        targpt(1) = targ(1,i)
        targpt(2) = targ(2,i)
        targpt(3) = targ(3,i)

        do k = 1,nsource
            gradlapker(1,k) = targpt(1) - source(1,k)
            gradlapker(2,k) = targpt(2) - source(2,k)
            gradlapker(3,k) = targpt(3) - source(3,k)
            temp = sqrt(gradlapker(1,k)**2 + gradlapker(2,k)**2 + &
                gradlapker(3,k)**2)
            temp = temp**3
            gradlapker(1,k) = gradlapker(1,k)/temp
            gradlapker(2,k) = gradlapker(2,k)/temp
            gradlapker(3,k) = gradlapker(3,k)/temp
        enddo

        do k = 1,nsource
            curlA(1,i) = curlA(1,i) + &
                j(2,k)*dcmplx(gradlapker(3,k)) - j(3,k)*dcmplx(gradlapker(2,k)) 
            curlA(2,i) = curlA(2,i) + &
                j(3,k)*dcmplx(gradlapker(1,k)) - j(1,k)*dcmplx(gradlapker(3,k)) 
            curlA(3,i) = curlA(3,i) + & 
                j(1,k)*dcmplx(gradlapker(2,k)) - j(2,k)*dcmplx(gradlapker(1,k)) 
        enddo
        
        curlA(1,i) = pi*curlA(1,i)/nsource
        curlA(2,i) = pi*curlA(2,i)/nsource
        curlA(3,i) = pi*curlA(3,i)/nsource
    enddo

    return
end subroutine curl_vector_potential

subroutine curl_vector_potential_torus(Nphi, Ntheta, rmin, rmaj, jmag, ntarg, &
    targ, curlA)
    !
    ! Computes the curl of the vector potential A, where
    ! 
    !         /\
    !         |       
    !     A = |   g_0(x,x') J(x') ds(x'), 
    !         |    
    !        \/ S
    ! 
    ! J is a current on a toroidal surface S with minor and major radii <rmin> 
    ! and <rmaj>, respectively, discretized in the toroidal and poloidal angles 
    ! with <Nphi> and <Ntheta> points, respectively, that surrounds the solution 
    ! domain, and g_0 is the Laplace kernel |x-x'|^-1. J is a surface current on
    ! S with magnitude <jmag>.
    ! 
    ! ==========================================================================
    ! 
    ! Input:
    !   Nphi - integer
    !     number of discretization points in the toroidal direction
    ! 
    !   Ntheta - integer
    !     number of discretization points in the poloidal direction
    ! 
    !   rmin - real *8
    !     minor radius of S
    ! 
    !   rmaj - real *8
    !     major radius of S
    ! 
    !   jmag - complex *16
    !     magnitude of current J
    ! 
    !   ntarg - integer
    !     total number of discretization points on solution domain boundary
    ! 
    !   targ - real *8 (3,ntarg)
    !     discretization nodes on solution domain boundary
    ! 
    ! Output:
    !   curlA - complex *16 (3,ntarg)
    ! 
    implicit none 
    integer Nphi,Ntheta,nsource,ntarg
    real *8 rmin,rmaj,phi(Nphi),theta(Ntheta),targ(3,ntarg),jmag
    complex *16 curlA(3,ntarg)
    real *8, allocatable :: source(:,:),j(:,:),gradlapker(:,:)
    complex *16, allocatable :: integrand(:,:)

    integer i,k,l
    real *8 targpt(3),temp,pi

    pi = 4*atan(1.0d0)

    nsource = Nphi*Ntheta

    ! set up phi and theta grids
    do i=1,Nphi
        phi(i) = 2*pi*(i-1)/Nphi
    enddo
    do i=1,Ntheta
        theta(i) = 2*pi*(i-1)/Ntheta
    enddo

    ! set up current geometry
    allocate(source(3,nsource))
    do i = 1,Ntheta
        do k = 1,Nphi
            source(1,Nphi*(i-1)+k) = (rmaj+rmin*cos(theta(i)))*cos(phi(k))
            source(2,Nphi*(i-1)+k) = (rmaj+rmin*cos(theta(i)))*sin(phi(k))
            source(3,Nphi*(i-1)+k) = rmin*sin(theta(i))
        enddo
    enddo
    allocate(j(3,nsource))
    do i=1,Ntheta
        do k=1,Nphi
            j(1,Nphi*(i-1)+k) = -jmag*rmin*sin(theta(i))*cos(phi(k))
            j(2,Nphi*(i-1)+k) = -jmag*rmin*sin(theta(i))*sin(phi(k))
            j(3,Nphi*(i-1)+k) = jmag*rmin*cos(theta(i))
        enddo
    enddo
    do i = 1,ntarg
        do k = 1,3
            curlA(k,i) = 0
        enddo
    enddo

    allocate(gradlapker(3,nsource))
    allocate(integrand(3,nsource))
    ! perform integral for each target point
    do i = 1,ntarg
        targpt(1) = targ(1,i)
        targpt(2) = targ(2,i)
        targpt(3) = targ(3,i)

        ! construct the gradient of the Laplace kernel
        do k = 1,nsource
            gradlapker(1,k) = targpt(1) - source(1,k)
            gradlapker(2,k) = targpt(2) - source(2,k)
            gradlapker(3,k) = targpt(3) - source(3,k)
            temp = sqrt(gradlapker(1,k)**2 + gradlapker(2,k)**2 + &
                gradlapker(3,k)**2)
            temp = temp**3
            gradlapker(1,k) = gradlapker(1,k)/temp
            gradlapker(2,k) = gradlapker(2,k)/temp
            gradlapker(3,k) = gradlapker(3,k)/temp
        enddo

        ! form integrand including part of the change-of-variables weight 
        ! do k = 1,Ntheta
        !     do l = 1,Nphi
        !         integrand(1,Nphi*(k-1)+l) = j(2,k)*gradlapker(3,k) &
        !             - j(3,k)*gradlapker(2,k)
        !         integrand(1,Nphi*(k-1)+l) = integrand(1,Nphi*(k-1)+l)*&
        !             (rmaj+rmin*cos(theta(k)))
        !         integrand(2,Nphi*(k-1)+l) = j(3,k)*gradlapker(1,k) &
        !             - j(1,k)*gradlapker(3,k)
        !         integrand(2,Nphi*(k-1)+l) = integrand(2,Nphi*(k-1)+l)*&
        !             (rmaj+rmin*cos(theta(k)))
        !         integrand(3,Nphi*(k-1)+l) = j(1,k)*gradlapker(2,k) &
        !             - j(2,k)*gradlapker(1,k)
        !         integrand(3,Nphi*(k-1)+l) = integrand(3,Nphi*(k-1)+l)*&
        !             (rmaj+rmin*cos(theta(k)))
        !     enddo
        ! enddo

        do k = 1,Ntheta
            do l = 1,Nphi
                integrand(1,Nphi*(k-1)+l) = & 
                    j(2,Nphi*(k-1)+l)*gradlapker(3,Nphi*(k-1)+l) &
                    - j(3,Nphi*(k-1)+l)*gradlapker(2,Nphi*(k-1)+l)
                integrand(1,Nphi*(k-1)+l) = integrand(1,Nphi*(k-1)+l)*&
                    (rmaj+rmin*cos(theta(k)))
                integrand(2,Nphi*(k-1)+l) = &
                    j(3,Nphi*(k-1)+l)*gradlapker(1,Nphi*(k-1)+l) &
                    - j(1,Nphi*(k-1)+l)*gradlapker(3,Nphi*(k-1)+l)
                integrand(2,Nphi*(k-1)+l) = integrand(2,Nphi*(k-1)+l)*&
                    (rmaj+rmin*cos(theta(k)))
                integrand(3,Nphi*(k-1)+l) = &
                    j(1,Nphi*(k-1)+l)*gradlapker(2,Nphi*(k-1)+l) &
                    - j(2,Nphi*(k-1)+l)*gradlapker(1,Nphi*(k-1)+l)
                integrand(3,Nphi*(k-1)+l) = integrand(3,Nphi*(k-1)+l)*&
                    (rmaj+rmin*cos(theta(k)))
            enddo
        enddo

        do k = 1,nsource
            curlA(1,i) = curlA(1,i) + integrand(1,k)
            curlA(2,i) = curlA(2,i) + integrand(2,k)
            curlA(3,i) = curlA(3,i) + integrand(3,k)
        enddo
        
        curlA(1,i) = rmin*pi*curlA(1,i)/nsource
        curlA(2,i) = rmin*pi*curlA(2,i)/nsource
        curlA(3,i) = rmin*pi*curlA(3,i)/nsource
    enddo

    return
end subroutine curl_vector_potential_torus

subroutine curl_vector_potential_loop(Ntheta, rmin, rmaj, jmag, ntarg, targ, &
    curlA)
    ! 
    ! Same as above, but J is a current on a loop in the same plane as a 
    ! cross-section of Omega.
    ! 
    ! =================================================================
    ! 
    ! Input: 
    !   Ntheta - integer
    !     number of discretization points along loop 
    ! 
    !   rmin - real *8
    !     radius of loop
    ! 
    !   rmaj - real *8
    !     distance between origin and loop center
    ! 
    !   jmag - complex *16
    !     magnitude of current J
    ! 
    !   ntarg - integer
    !     total number of discretization points on solution domain boundary
    ! 
    !   targ - real *8 (3,ntarg)
    !     discretization nodes on solution domain boundary
    ! 
    ! Output:
    !   curlA - complex *16 (3,ntarg)
    ! 

    implicit none
    integer Ntheta,ntarg
    real *8 rmin,rmaj,theta(Ntheta),targ(3,ntarg),jmag,source(3,Ntheta)
    complex *16 curlA(3,ntarg),integrand(3,Ntheta)
    real *8 j(3,Ntheta),gradlapker(3,Ntheta)

    integer i,k,l
    real *8 targpt(3),temp,pi

    pi = 4*atan(1.0d0)

    ! set up theta grid
    do i=1,Ntheta
        theta(i) = 2*pi*(i-1)/Ntheta
    enddo

    ! set up current geometry
    do i = 1,Ntheta
        source(1,i) = rmaj+rmin*cos(theta(i))
        source(2,i) = 0
        source(3,i) = rmin*sin(theta(i))
        j(1,i) = -jmag*rmin*sin(theta(i))
        j(2,i) = 0
        j(3,i) = jmag*rmin*cos(theta(i))
    enddo
    do i = 1,ntarg
        do k = 1,3
            curlA(k,i) = 0
        enddo
    enddo

    ! perform integral for each target point
    do i = 1,ntarg
        do k = 1,3
            targpt(k) = targ(k,i)
        enddo

        ! construct the gradient of the Laplace kernel
        do k = 1,Ntheta
            gradlapker(1,k) = targpt(1) - source(1,k)
            gradlapker(2,k) = targpt(2) - source(2,k)
            gradlapker(3,k) = targpt(3) - source(3,k)
            temp = sqrt(gradlapker(1,k)**2 + gradlapker(2,k)**2 + &
                gradlapker(3,k)**2)
            temp = temp**3
            gradlapker(1,k) = gradlapker(1,k)/temp
            gradlapker(2,k) = gradlapker(2,k)/temp
            gradlapker(3,k) = gradlapker(3,k)/temp
        enddo

        ! form integrand without change-of-variables weight
        do k = 1,Ntheta
            integrand(1,k) = j(2,k)*gradlapker(3,k) - j(3,k)*gradlapker(2,k)
            integrand(2,k) = j(3,k)*gradlapker(1,k) - j(1,k)*gradlapker(3,k)
            integrand(3,k) = j(1,k)*gradlapker(2,k) - j(2,k)*gradlapker(1,k)
        enddo

        do k = 1,Ntheta
            do l = 1,3
                curlA(l,i) = curlA(l,i) + integrand(l,k)
            enddo
        enddo

        do k = 1,3
            curlA(k,i) = rmin*curlA(k,i)/(2*Ntheta)
        enddo
    enddo

end subroutine curl_vector_potential_loop

! program test_vacuum
!     implicit none
!     real *8 pi,rmin,rmaj
!     integer Nphi,Ntheta,i,k
!     real *8, allocatable :: phi(:),theta(:)

!     ! vars for setup_geom_vacuum()
!     integer igeomtype,norder,npols,npatches,ipars(2),ifplot
!     real *8, allocatable :: targvals(:,:),targcoefs(:,:)
!     character (len=300) fname

!     ! vars for curl_vector_potential()
!     integer nsource,ntarg
!     real *8, allocatable :: source(:,:),targ(:,:)
!     complex *16, allocatable :: j(:,:),curlA(:,:)

!     pi = 4*atan(1.0d0)
!     rmaj = 1.5d0
!     rmin = 1.0d0
!     Nphi = 10
!     Ntheta = 16

!     ! set up phi and theta grids
!     allocate(phi(Nphi),theta(Ntheta))
!     do i=1,Nphi
!         phi(i) = 2*pi*(i-1)/(Nphi+1)
!     enddo
!     do i=1,Ntheta
!         theta(i) = 2*pi*(i-1)/(Ntheta+1)
!     enddo

!     ! get solution domain geometry
!     igeomtype = 4 ! set to 4 for axisymmetric geometry
!     norder = 8 
!     npols = (norder+1)*(norder+2)/2
!     ifplot = 0
!     fname = 'torus.vtk'
!     ipars(1) = 6
!     if(igeomtype.eq.2) ipars(1) = 10
!     ipars(2) = ipars(1)*3
!     npatches = 2*ipars(1)*ipars(2)
!     nsource = Nphi*Ntheta
!     ntarg = npols*npatches
!     allocate(targvals(12,ntarg),targcoefs(9,ntarg))
!     call setup_geom_vacuum(igeomtype,norder,npatches,ipars,targvals,&
!         targcoefs,ifplot,fname)

!     ! set up current geometry
!     allocate(source(3,nsource),targ(3,ntarg))
!     allocate(j(3,nsource),curlA(3,ntarg))
!     do i=1,Nphi
!         do k=1,Ntheta
!             source(1,Nphi*(i-1)+k) = (rmaj+rmin*cos(theta(k)))*cos(phi(i))
!             source(2,Nphi*(i-1)+k) = (rmaj+rmin*cos(theta(k)))*sin(phi(i))
!             source(3,Nphi*(i-1)+k) = rmin*sin(theta(k))
!         enddo
!     enddo
!     do i=1,ntarg
!         do k=1,3
!             targ(k,i) = targvals(k,i)
!         enddo
!     enddo
!     do i=1,nsource
!         j(1,i) = 1.0
!         j(2,i) = 1.0
!         j(3,i) = 1.0
!     enddo

!     call curl_vector_potential(nsource,source,j,ntarg,targ,curlA)
!     ! open(11,file='curlA.txt')
!     ! open(12,file='targ.txt')
!     ! do i=1,ntarg
!     !     do k=1,3
!     !         write(11,*) curlA(k,i)
!     !         write(12,*) targ(k,i)
!     !     enddo
!     ! enddo
!     ! close(12)
!     ! close(11)

!     print *,"ntarg = ",ntarg
! end program test_vacuum