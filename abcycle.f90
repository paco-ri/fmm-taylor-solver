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

    igeomtype = 6
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