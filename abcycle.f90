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

program abcycle
    implicit none

    integer :: norder, npols, npatches, npts, m, na, nb, ipars(2)
    integer, allocatable :: norders(:), ixyzs(:), iptype(:)
    real *8, allocatable :: srccoefs(:,:), srcvals(:,:), surfdivf(:)
    real *8, allocatable :: avals(:,:),bvals(:,:)
    real *8, allocatable :: awts(:),bwts(:)
    real *8, allocatable :: auv(:,:),buv(:,:)
    integer, allocatable :: apatches(:),bpatches(:)
    
    integer i, igeomtype, ifplot
    character *300 fname

    igeomtype = 3
    norder = 8
    npols = (norder+1)*(norder+2)/2
    ifplot = 0
    fname = 'torus.vtk'

    ipars(1) = 6
    if(igeomtype.eq.2) ipars(1) = 10
    ipars(2) = ipars(1)*3
    npatches = 2*ipars(1)*ipars(2)

    npts = npatches*npols
    allocate(norders(npatches),ixyzs(npatches+1),iptype(npatches))
    allocate(srccoefs(9,npts),srcvals(12,npts),surfdivf(npts))

    do i = 1,npatches
        norders(i) = norder
        ixyzs(i) = 1 + (i-1)*npols 
        iptype(i) = 1
    enddo
    ixyzs(npatches+1) = npts+1

    call setup_geom_vacuum(igeomtype,norder,npatches,ipars,srcvals,srccoefs,&
        ifplot,fname)
    open(11, file="absrcvals.txt")
    write(11,*) srcvals
    close(11)

    m = 16
    na = ipars(2)*m
    nb = ipars(1)*m
    allocate(avals(9,na),awts(na),auv(2,na),apatches(na))
    allocate(bvals(9,nb),bwts(nb),buv(2,nb),bpatches(nb))
    call get_ab_cycles_torusparam(npatches,norders,ixyzs,iptype,npts,srccoefs,&
        srcvals,ipars,m,na,avals,awts,apatches,auv,nb,bvals,bwts,bpatches,buv)

    open(14, file="avals.txt")
    write(14,*) avals
    close(14)
    open(15, file="bvals.txt")
    write(15,*) bvals
    close(15)
end program abcycle