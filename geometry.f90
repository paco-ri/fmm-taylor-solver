subroutine setup_geom_vacuum(igeomtype, norder, npatches, ipars, srcvals, &
    srccoefs, ifplot, fname)
    ! 
    ! Geometry setup function lifted from superconductor-type1/test/
    ! test_ab_cycle.f
    !
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
    external xtri_wtorus_star_eval
      
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
      
    ! fat torus
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
         
        p1(1) = 0.75d0 ! 1.0d0
        p1(2) = 2.0d0 ! 10.0d0 ! 1.75d0
        p1(3) = 0.0d0 ! 0.25d0

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

    if(igeomtype.eq.6) then
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
         
        p1(1) = 0.5d0 ! 1.0d0
        p1(2) = 1.0d0 ! 1.75d0
        p1(3) = 0.0d0 ! 0.3d0 ! 0.25d0
        p1(4) = 0.2d0

        p2(1) = 1.0d0
        p2(2) = 1.0d0
        p2(3) = 1.0d0
        ! number of oscillations
        p4(1) = 0.0d0
        p4(2) = 2.0d0

        ptr1 => triaskel(1,1,1)
        ptr2 => p1(1)
        ptr3 => p2(1)
        ptr4 => p4(1)
        xtri_geometry => xtri_wtorus_star_eval
        if(ifplot.eq.1) then
            call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2,ptr3,&
                ptr4, norder,'Triangulated surface of the torus')
        endif
        
        
        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,ptr3,ptr4,npols,uvs,&
            umatr,srcvals,srccoefs)
    endif

    return  
end subroutine setup_geom_vacuum

subroutine xtri_wtorus_star_eval(itri, u, v, xyz, dxyzduv, triainfo, &
    radii, scales, p4)
  implicit real *8 (a-h,o-z)
  real *8 :: xyz(3), dxyzduv(3,2), triainfo(3,3,*), scales(3)
  real *8 :: radii(4), p4(2)

  !
  ! project the triangle itri in triainfo onto a torus
  !
  !    Input:
  ! itri - triangle number to map
  ! u,v - local uv coordinates on triangle itri
  ! triainfo - flat skeleton triangle info
  ! radii - minor radius, major radius, radius of toroidal oscillation,
  !   radius of poloidal oscillation 
  ! scales - scaling for x,y,z components from the standard torus
  ! p4 - number of oscillations (must be an integer currently recast
  !   as a double precision number)
  !
  !    Output:
  ! xyz - point on the sphere
  ! dxyzduv - first and second derivative information
  !
  !

  x0=triainfo(1,1,itri)
  y0=triainfo(2,1,itri)
  z0=triainfo(3,1,itri)

  x1=triainfo(1,2,itri)
  y1=triainfo(2,2,itri)
  z1=triainfo(3,2,itri)

  x2=triainfo(1,3,itri)
  y2=triainfo(2,3,itri)
  z2=triainfo(3,3,itri)

  rminor = radii(1)
  rmajor = radii(2)
  rwave = radii(3)
  rwavep = radii(4)

  a = scales(1)
  b = scales(2)
  c = scales(3)

  nosc = p4(1)
  noscp = p4(2)


  !
  ! ... process the geometry, return the point location on the almond
  ! and the derivatives with respect to u and v
  !
  s = x0+u*(x1-x0)+v*(x2-x0)
  t = y0+u*(y1-y0)+v*(y2-y0)

  rr = rmajor+rminor*cos(t)*(1+rwavep*cos(noscp*t))+rwave*cos(nosc*s)


  xyz(1) = a*rr*cos(s)
  xyz(2) = b*rr*sin(s)
  xyz(3) = c*rminor*sin(t)

  dsdu = (x1-x0)
  dsdv = (x2-x0)
  dtdu = (y1-y0)
  dtdv = (y2-y0)

  drrds = -nosc*rwave*sin(nosc*s)
  drrdt = rminor*(-sin(t)*(1+rwavep*cos(noscp*t))-cos(t)*rwavep*noscp*sin(noscp*t))

  dxds = a*drrds*cos(s) - a*rr*sin(s)
  dyds = b*drrds*sin(s) + b*rr*cos(s)
  dzds = 0

  dxdt = a*drrdt*cos(s)
  dydt = b*drrdt*sin(s)
  dzdt = c*rminor*cos(t)
  
  dxdu = dxds*dsdu + dxdt*dtdu
  dydu = dyds*dsdu + dydt*dtdu
  dzdu = dzds*dsdu + dzdt*dtdu

  dxdv = dxds*dsdv + dxdt*dtdv
  dydv = dyds*dsdv + dydt*dtdv
  dzdv = dzds*dsdv + dzdt*dtdv

  dxyzduv(1,1) = dxdu
  dxyzduv(2,1) = dydu
  dxyzduv(3,1) = dzdu

  dxyzduv(1,2) = dxdv
  dxyzduv(2,2) = dydv
  dxyzduv(3,2) = dzdv

  return
end subroutine xtri_wtorus_star_eval

