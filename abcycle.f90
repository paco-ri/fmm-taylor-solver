program abcycle
    implicit none

    integer :: norder, npols, npatches, npts, m, na, nb, ipars(2)
    integer, allocatable :: norders(:), ixyzs(:), iptype(:)
    real *8, allocatable :: srccoefs(:,:), srcvals(:,:), surfdivf(:)
    real *8, allocatable :: avals(:,:),bvals(:,:)
    real *8, allocatable :: awts(:),bwts(:)
    real *8, allocatable :: auv(:,:),buv(:,:)
    integer, allocatable :: apatches(:),bpatches(:)
    
    integer i, igeomtype

    igeomtype = 4
    norder = 8
    npols = (norder+1)*(norder+2)/2

    if(igeomtype.eq.2.or.igeomtype.eq.4) then
        ipars(1) = 6
        if(igeomtype.eq.2) ipars(1) = 10
        ipars(2) = ipars(1)*3
        npatches = 2*ipars(1)*ipars(2)
    endif
    npts = npatches*npols
    allocate(norders(npatches),ixyzs(npatches+1),iptype(npatches))
    allocate(srccoefs(9,npts),srcvals(12,npts),surfdivf(npts))

    if(igeomtype.eq.2) then 
        open(12, file="axisrccoefs.txt")
        read(12,*) srccoefs
        close(12)
        open(13, file="axisrcvals.txt")
        read(13,*) srcvals
        close(13)
    endif
    if(igeomtype.eq.4) then
        open(12, file="torussrccoefs.txt")
        read(12,*) srccoefs
        close(12)
        open(13, file="torussrcvals.txt")
        read(13,*) srcvals
        close(13)
    endif

    do i = 1,npatches
        norders(i) = norder
        ixyzs(i) = 1 + (i-1)*npols 
        iptype(i) = 1
    enddo
    ixyzs(npatches+1) = 1 + npts

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