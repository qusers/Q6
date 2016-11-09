! (C) 2014 Uppsala Molekylmekaniska HB, Uppsala, Sweden
! simprep.f90
! based on
! md.f90
! by Johan Åqvist, John Marelius, Anders Kaplan, Isabella Feierberg, Martin Nervall & Martin Almlöf
! general simulation setup 

module SIMPREP

! used modules
!use PROFILING
use SIZES
use TRJ
use MPIGLOB
use QATOM
use VERSIONS
use GLOBALS
use EXC
use QMATH
use QALLOC
use QCP
use NONBONDENE
use MPIGLOB
!$ use omp_lib
implicit none



!----START OF PUBLIC SUBROUTINES
contains


!----------------------------------------------------------------------


subroutine simprep_startup
! initialise used modules
call qatom_startup
call trj_startup


! initialise constants
pi = 4.0_prec*atan(one)
deg2rad = pi/180.0_prec

end subroutine simprep_startup


!----------------------------------------------------------------------

subroutine simprep_shutdown
! call used modules' shutdown subroutines
if (nodeid .eq. 0) then
if (use_excluded_groups) then
call excluded_shutdown(ngroups_gc)
end if
call qatom_shutdown
call index_shutdown
call trj_shutdown
endif
end subroutine simprep_shutdown

!----------------------------------------------------------------------

#ifdef USE_GRID
subroutine create_grid_pp
! *** local variables
integer						:: i,j,k,ndim,num
integer						:: i3
real(kind=prec)					:: xtmp,ytmp,ztmp
real(kind=prec)					:: xmax,ymax,zmax
real(kind=prec)					:: xmin,ymin,zmin
!have 1.75 times the cutoff length to account for large charge groups in solute residues
real(kind=prec)					:: gridRc
integer						:: li,lj,lk,ui,uj,uk

gridRc = 1.75_prec * Rcpp

if(use_PBC) then
! when using period boundary conditions, the grids are generated accoring to the box centre and the desired  system size
! this depends on the volume, so we need to be careful with the coordinates to not miss any atom
! so the coordinate information needs to be checked/reset in populate grids
! Paul Bauer 2015

! one problem could be the large movement of molecules if they are not put back in the box
! this needs to be checked for, too
! also, the coordinates of each grid need to be corrected after a box movement

! number of one dimensional spacings cubed
pp_gridnum=int(ceiling(((boxlength(1)+boxlength(2)+boxlength(3))/3)/gridRc)**3)
ndim=ceiling(((boxlength(1)+boxlength(2)+boxlength(3))/3)/gridRc)
pp_ndim=ndim
if (ncgp_solute .gt. 500) then
gridstor_pp = (ncgp_solute*8/pp_gridnum)+100
else
gridstor_pp = ncgp_solute
endif
call allocate_grid_pp

! set starting coordinates for first grid
xtmp=boxcentre(1)-boxlength(1)/2
ytmp=boxcentre(2)-boxlength(2)/2
ztmp=boxcentre(3)-boxlength(3)/2

num=0
do i=1,ndim
 ytmp=boxcentre(2)-boxlength(2)/2
 do j=1,ndim
  ztmp=boxcentre(3)-boxlength(3)/2
  do k=1,ndim
   num=num+1

   grid_pp_int(num,:,:,:)=.false.
   grid_pp_int(num,i,j,k)=.true.

! stuff to account for PBC wraparound for the interaction matrix
   if(i+1.le.ndim) then
    ui=i+1
   else
    ui=1
   end if
   if(i-1.ge.1) then
    li=i-1
   else
    li=ndim
   end if
   if(j+1.le.ndim) then
    uj=j+1
   else
    uj=1
   end if
   if(j-1.ge.1) then
    lj=j-1
   else
    lj=ndim
   end if
   if(k+1.le.ndim) then
    uk=k+1
   else
    uk=1
   end if
   if(k-1.ge.1) then
    lk=k-1
   else
    lk=ndim
   end if
   
! generating interaction matrix
   grid_pp_int(num,ui,uj,uk)=.true.
   grid_pp_int(num,ui,uj,lk)=.true.
   grid_pp_int(num,ui,lj,uk)=.true.
   grid_pp_int(num,ui,lj,lk)=.true.
   grid_pp_int(num,li,uj,uk)=.true.
   grid_pp_int(num,li,uj,lk)=.true.
   grid_pp_int(num,li,lj,uk)=.true.
   grid_pp_int(num,li,lj,lk)=.true.
   grid_pp_int(num,ui,uj,k)=.true.
   grid_pp_int(num,ui,lj,k)=.true.
   grid_pp_int(num,li,uj,k)=.true.
   grid_pp_int(num,li,lj,k)=.true.
   grid_pp_int(num,ui,j,uk)=.true.
   grid_pp_int(num,ui,j,lk)=.true.
   grid_pp_int(num,li,j,uk)=.true.
   grid_pp_int(num,li,j,lk)=.true.
   grid_pp_int(num,i,uj,uk)=.true.
   grid_pp_int(num,i,uj,lk)=.true.
   grid_pp_int(num,i,lj,uk)=.true.
   grid_pp_int(num,i,lj,lk)=.true.
   grid_pp_int(num,ui,j,k)=.true.
   grid_pp_int(num,li,j,k)=.true.
   grid_pp_int(num,i,j,uk)=.true.
   grid_pp_int(num,i,j,lk)=.true.
   grid_pp_int(num,i,uj,k)=.true.
   grid_pp_int(num,i,lj,k)=.true.

   grid_pp(num)%x=xtmp
   grid_pp(num)%y=ytmp
   grid_pp(num)%z=ztmp

   grid_pp(num)%xend=xtmp+gridRc
   grid_pp(num)%yend=ytmp+gridRc
   grid_pp(num)%zend=ztmp+gridRc

   if (k .eq. ndim) then
    if (grid_pp(num)%zend .lt. (boxcentre(3)+(boxlength(3)/2))) then
     grid_pp(num)%zend=boxcentre(3)+(boxlength(3)/2)
    end if
   end if
   ztmp=ztmp+gridRc
  end do
  if (j .eq. ndim) then
   if(grid_pp(num)%yend .lt. (boxcentre(2)+(boxlength(2)/2))) then
    grid_pp(num)%yend=boxcentre(2)+(boxlength(2)/2)
   end if
  end if
  ytmp=ytmp+gridRc
 end do
 if (i .eq. ndim) then
  if(grid_pp(num)%xend .lt. (boxcentre(1)+(boxlength(1)/2))) then
   grid_pp(num)%xend=boxcentre(1)+(boxlength(1)/2)
  end if
 end if
 xtmp=xtmp+gridRc
end do

else
! we are using spherical boundaries, meaning that we now have some problems
! we need to know the largest and smallest values for all coordiantes to generate the grid
i = 1
i3 = i*3-3
xmax=x(i3+1)
ymax=x(i3+2)
zmax=x(i3+3)
xmin=xmax
ymin=ymax
zmin=zmax
! now iterate over all cooridnates to find min and max
! atom one is already read in, don't need to read it again
do i=2,natom
 i3=i*3-3
 if (xmin .gt. x(i3+1)) xmin=x(i3+1)
 if (ymin .gt. x(i3+2)) ymin=x(i3+2)
 if (zmin .gt. x(i3+3)) zmin=x(i3+3)
 if (xmax .lt. x(i3+1)) xmax=x(i3+1)
 if (ymax .lt. x(i3+2)) ymax=x(i3+2)
 if (zmax .lt. x(i3+3)) zmax=x(i3+3)
end do
! now add some buffer to each of them -> half of Rcpp
xmin=xmin-(gridRc/3)
ymin=ymin-(gridRc/3)
zmin=zmin-(gridRc/3)
xmax=xmax+(gridRc/3)
ymax=ymax+(gridRc/3)
zmax=zmax+(gridRc/3)
! the system of boxes will be cubic in any case, so now find the largest number of 
! grid spaces one any side
ndim=ceiling((xmax-xmin)/gridRc)
if (ceiling((ymax-ymin)/gridRc) .gt. ndim ) then
 ndim=ceiling((ymax-ymin)/gridRc)
else if (ceiling((zmax-zmin)/gridRc) .gt. ndim ) then
 ndim=(ceiling((zmax-zmin)/gridRc))
end if
pp_gridnum=ndim**3
if (ncgp_solute .gt. 500) then
gridstor_pp = (ncgp_solute*8/pp_gridnum)+100
else
gridstor_pp = ncgp_solute
endif
pp_ndim=ndim
call allocate_grid_pp

xtmp=xmin
ytmp=ymin
ztmp=zmin

num=1
do i=1,ndim
 ytmp=ymin
 do j=1,ndim
  ztmp=zmin
  do k=1,ndim
! now fill the interaction matrix
! no wraparound in this case, so be careful
! when there would be wraparound, just point to my own cell again

grid_pp_int(num,:,:,:)=.false.
grid_pp_int(num,i,j,k)=.true.


   if((i+1).le.ndim) then
    ui=i+1
   else
    ui=i
   end if
   if((i-1).ge.1) then
    li=i-1
   else
    li=i
   end if
   if((j+1).le.ndim) then
    uj=j+1
   else
    uj=j
   end if
   if((j-1).ge.1) then
    lj=j-1
   else
    lj=j
   end if
   if((k+1).le.ndim) then
    uk=k+1
   else
    uk=k
   end if
   if((k-1).ge.1) then
    lk=k-1
   else
    lk=k
   end if

   grid_pp_int(num,ui,uj,uk)=.true.
   grid_pp_int(num,ui,uj,lk)=.true.
   grid_pp_int(num,ui,lj,uk)=.true.
   grid_pp_int(num,ui,lj,lk)=.true.
   grid_pp_int(num,li,uj,uk)=.true.
   grid_pp_int(num,li,uj,lk)=.true.
   grid_pp_int(num,li,lj,uk)=.true.
   grid_pp_int(num,li,lj,lk)=.true.
   grid_pp_int(num,ui,uj,k)=.true.
   grid_pp_int(num,ui,lj,k)=.true.
   grid_pp_int(num,li,uj,k)=.true.
   grid_pp_int(num,li,lj,k)=.true.
   grid_pp_int(num,ui,j,uk)=.true.
   grid_pp_int(num,ui,j,lk)=.true.
   grid_pp_int(num,li,j,uk)=.true.
   grid_pp_int(num,li,j,lk)=.true.
   grid_pp_int(num,i,uj,uk)=.true.
   grid_pp_int(num,i,uj,lk)=.true.
   grid_pp_int(num,i,lj,uk)=.true.
   grid_pp_int(num,i,lj,lk)=.true.
   grid_pp_int(num,ui,j,k)=.true.
   grid_pp_int(num,li,j,k)=.true.
   grid_pp_int(num,i,j,uk)=.true.
   grid_pp_int(num,i,j,lk)=.true.
   grid_pp_int(num,i,uj,k)=.true.
   grid_pp_int(num,i,lj,k)=.true.

   grid_pp(num)%x=xtmp
   grid_pp(num)%y=ytmp
   grid_pp(num)%z=ztmp

   grid_pp(num)%xend=xtmp+gridRc
   grid_pp(num)%yend=ytmp+gridRc
   grid_pp(num)%zend=ztmp+gridRc

   num=num+1
   ztmp=ztmp+gridRc
  end do
  ytmp=ytmp+gridRc
 end do
 xtmp=xtmp+gridRc
end do

end if

end subroutine create_grid_pp 

subroutine create_grid_pw
! *** local variables
integer                                         :: i,j,k,ndim,num
integer                                         :: i3
real(kind=prec)                                 :: xtmp,ytmp,ztmp
real(kind=prec)                                 :: xmax,ymax,zmax
real(kind=prec)                                 :: xmin,ymin,zmin
! again make grid larger than cutoff to account for charge groups
! can be smaller this time because we only have water as one of them
real(kind=prec)					:: gridRc
integer                                         :: li,lj,lk,ui,uj,uk

gridRc = 1.5_prec * Rcpw

if(use_PBC) then
! when using period boundary conditions, the grids are generated accoring to the box centre and the desired  system size
! this depends on the volume, so we need to be careful with the coordinates to not miss any atom
! so the coordinate information needs to be checked/reset in populate grids
! Paul Bauer 2015

! one problem could be the large movement of molecules if they are not put back in the box
! this needs to be checked for, too
! also, the coordinates of each grid need to be corrected after a box movement

! number of one dimensional spacings cubed
pw_gridnum=int(ceiling(((boxlength(1)+boxlength(2)+boxlength(3))/3)/gridRc)**3)
ndim=ceiling(((boxlength(1)+boxlength(2)+boxlength(3))/3)/gridRc)
pw_ndim=ndim
if (nwat .gt. 500) then
gridstor_pw = (nwat*8/pw_gridnum)+100
else
gridstor_pw = nwat
endif
call allocate_grid_pw

! set starting coordinates for first grid
xtmp=boxcentre(1)-(boxlength(1)/2)
ytmp=boxcentre(2)-(boxlength(2)/2)
ztmp=boxcentre(3)-(boxlength(3)/2)


num=0
do i=1,ndim
 ytmp=boxcentre(2)-(boxlength(2)/2)
 do j=1,ndim
  ztmp=boxcentre(3)-(boxlength(3)/2)
  do k=1,ndim
   num=num+1

   grid_pw_int(num,:,:,:)=.false.
   grid_pw_int(num,i,j,k)=.true.

! stuff to account for PBC wraparound for the interaction matrix
   if(i+1.le.ndim) then
    ui=i+1
   else
    ui=1
   end if
   if(i-1.ge.1) then
    li=i-1
   else
    li=ndim
   end if
   if(j+1.le.ndim) then
    uj=j+1
   else
    uj=1
   end if
   if(j-1.ge.1) then
    lj=j-1
   else
    lj=ndim
   end if
   if(k+1.le.ndim) then
    uk=k+1
   else
    uk=1
   end if
   if(k-1.ge.1) then
    lk=k-1
   else
    lk=ndim
   end if

! generating interaction matrix
   grid_pw_int(num,ui,uj,uk)=.true.
   grid_pw_int(num,ui,uj,lk)=.true.
   grid_pw_int(num,ui,lj,uk)=.true.
   grid_pw_int(num,ui,lj,lk)=.true.
   grid_pw_int(num,li,uj,uk)=.true.
   grid_pw_int(num,li,uj,lk)=.true.
   grid_pw_int(num,li,lj,uk)=.true.
   grid_pw_int(num,li,lj,lk)=.true.
   grid_pw_int(num,ui,uj,k)=.true.
   grid_pw_int(num,ui,lj,k)=.true.
   grid_pw_int(num,li,uj,k)=.true.
   grid_pw_int(num,li,lj,k)=.true.
   grid_pw_int(num,ui,j,uk)=.true.
   grid_pw_int(num,ui,j,lk)=.true.
   grid_pw_int(num,li,j,uk)=.true.
   grid_pw_int(num,li,j,lk)=.true.
   grid_pw_int(num,i,uj,uk)=.true.
   grid_pw_int(num,i,uj,lk)=.true.
   grid_pw_int(num,i,lj,uk)=.true.
   grid_pw_int(num,i,lj,lk)=.true.
   grid_pw_int(num,ui,j,k)=.true.
   grid_pw_int(num,li,j,k)=.true.
   grid_pw_int(num,i,j,uk)=.true.
   grid_pw_int(num,i,j,lk)=.true.
   grid_pw_int(num,i,uj,k)=.true.
   grid_pw_int(num,i,lj,k)=.true.


   grid_pw(num)%x=xtmp
   grid_pw(num)%y=ytmp
   grid_pw(num)%z=ztmp

   grid_pw(num)%xend=xtmp+gridRc
   grid_pw(num)%yend=ytmp+gridRc
   grid_pw(num)%zend=ztmp+gridRc

   if (k .eq. ndim) then
    if (grid_pw(num)%zend .lt. (boxcentre(3)+(boxlength(3)/2))) then
     grid_pw(num)%zend=boxcentre(3)+(boxlength(3)/2)
    end if
   end if
   ztmp=ztmp+gridRc
  end do
  if (j .eq. ndim) then
   if(grid_pw(num)%yend .lt. (boxcentre(2)+(boxlength(2)/2))) then
    grid_pw(num)%yend=boxcentre(2)+(boxlength(2)/2)
   end if
  end if
  ytmp=ytmp+gridRc
 end do
 if (i .eq. ndim) then
  if(grid_pw(num)%xend .lt. (boxcentre(1)+(boxlength(1)/2))) then
   grid_pw(num)%xend=boxcentre(1)+(boxlength(1)/2)
  end if
 end if
 xtmp=xtmp+gridRc
end do

else
! we are using spherical boundaries, meaning that we now have some problems
! we need to know the largest and smallest values for all coordiantes to generate the grid
i = 1
i3 = i*3-3
xmax=x(i3+1)
ymax=x(i3+2)
zmax=x(i3+3)
xmin=xmax
ymin=ymax
zmin=zmax
! now iterate over all cooridnates to find min and max
! atom one is already read in, don't need to read it again
do i=2,natom
 i3=i*3-3
 if (xmin .gt. x(i3+1)) xmin=x(i3+1)
 if (ymin .gt. x(i3+2)) ymin=x(i3+2)
 if (zmin .gt. x(i3+3)) zmin=x(i3+3)
 if (xmax .lt. x(i3+1)) xmax=x(i3+1)
 if (ymax .lt. x(i3+2)) ymax=x(i3+2)
 if (zmax .lt. x(i3+3)) zmax=x(i3+3)
end do
! now add some buffer to each of them -> half of Rcpp
xmin=xmin-(gridRc/3)
ymin=ymin-(gridRc/3)
zmin=zmin-(gridRc/3)
xmax=xmax+(gridRc/3)
ymax=ymax+(gridRc/3)
zmax=zmax+(gridRc/3)
! the system of boxes will be cubic in any case, so now find the largest number of 
! grid spaces one any side
ndim=ceiling((xmax-xmin)/gridRc)
if (ceiling((ymax-ymin)/gridRc) .gt. ndim ) then
 ndim=ceiling((ymax-ymin)/gridRc)
else if (ceiling((zmax-zmin)/gridRc) .gt. ndim ) then
 ndim=(ceiling((zmax-zmin)/gridRc))
end if
pw_gridnum=ndim**3
pw_ndim=ndim
if (nwat .gt. 500) then
gridstor_pw = (nwat*8/pw_gridnum)+100
else
gridstor_pw = nwat
endif
call allocate_grid_pw

xtmp=xmin
ytmp=ymin
ztmp=zmin

num=1
do i=1,ndim
 ytmp = ymin
 do j=1,ndim
  ztmp = zmin
  do k=1,ndim
! now fill the interaction matrix
! no wraparound in this case, so be careful
! when there would be wraparound, just point to my own cell again

grid_pw_int(num,:,:,:)=.false.
grid_pw_int(num,i,j,k)=.true.


   if(i+1.le.ndim) then
    ui=i+1
   else
    ui=i
   end if
   if(i-1.ge.1) then
    li=i-1
   else
    li=i
   end if
   if(j+1.le.ndim) then
    uj=j+1
   else
    uj=j
   end if
   if(j-1.ge.1) then
    lj=j-1
   else
    lj=j
   end if
   if(k+1.le.ndim) then
    uk=k+1
   else
    uk=k
   end if
   if(k-1.ge.1) then
    lk=k-1
   else
    lk=k
   end if

   grid_pw_int(num,ui,uj,uk)=.true.
   grid_pw_int(num,ui,uj,lk)=.true.
   grid_pw_int(num,ui,lj,uk)=.true.
   grid_pw_int(num,ui,lj,lk)=.true.
   grid_pw_int(num,li,uj,uk)=.true.
   grid_pw_int(num,li,uj,lk)=.true.
   grid_pw_int(num,li,lj,uk)=.true.
   grid_pw_int(num,li,lj,lk)=.true.
   grid_pw_int(num,ui,uj,k)=.true.
   grid_pw_int(num,ui,lj,k)=.true.
   grid_pw_int(num,li,uj,k)=.true.
   grid_pw_int(num,li,lj,k)=.true.
   grid_pw_int(num,ui,j,uk)=.true.
   grid_pw_int(num,ui,j,lk)=.true.
   grid_pw_int(num,li,j,uk)=.true.
   grid_pw_int(num,li,j,lk)=.true.
   grid_pw_int(num,i,uj,uk)=.true.
   grid_pw_int(num,i,uj,lk)=.true.
   grid_pw_int(num,i,lj,uk)=.true.
   grid_pw_int(num,i,lj,lk)=.true.
   grid_pw_int(num,ui,j,k)=.true.
   grid_pw_int(num,li,j,k)=.true.
   grid_pw_int(num,i,j,uk)=.true.
   grid_pw_int(num,i,j,lk)=.true.
   grid_pw_int(num,i,uj,k)=.true.
   grid_pw_int(num,i,lj,k)=.true.

   grid_pw(num)%x=xtmp
   grid_pw(num)%y=ytmp
   grid_pw(num)%z=ztmp

   grid_pw(num)%xend=xtmp+gridRc
   grid_pw(num)%yend=ytmp+gridRc
   grid_pw(num)%zend=ztmp+gridRc

   num=num+1
   ztmp=ztmp+gridRc
  end do
  ytmp=ytmp+gridRc
 end do
 xtmp=xtmp+gridRc
end do

end if
end subroutine create_grid_pw

subroutine create_grid_ww
! *** local variables
integer                                         :: i,j,k,ndim,num
integer                                         :: i3
real(kind=prec)                                 :: xtmp,ytmp,ztmp
real(kind=prec)                                 :: xmax,ymax,zmax
real(kind=prec)                                 :: xmin,ymin,zmin
! again grids are slightly larger than cutoff, can be onle 1.25 because we only have water molecules
real(kind=prec)					:: gridRc 
integer                                         :: li,lj,lk,ui,uj,uk

gridRc = 1.25_prec * Rcww

if(use_PBC) then
! when using period boundary conditions, the grids are generated accoring to the box centre and the desired  system size
! this depends on the volume, so we need to be careful with the coordinates to not miss any atom
! so the coordinate information needs to be checked/reset in populate grids
! Paul Bauer 2015

! one problem could be the large movement of molecules if they are not put back in the box
! this needs to be checked for, too
! also, the coordinates of each grid need to be corrected after a box movement

! number of one dimensional spacings cubed
ww_gridnum=int(ceiling(((boxlength(1)+boxlength(2)+boxlength(3))/3)/gridRc)**3)
ndim=ceiling(((boxlength(1)+boxlength(2)+boxlength(3))/3)/gridRc)
ww_ndim=ndim
if (nwat .gt. 500) then
gridstor_ww = (nwat*8/ww_gridnum)+100
else
gridstor_ww = nwat
endif
call allocate_grid_ww
! set starting coordinates for first grid
xtmp=boxcentre(1)-(boxlength(1)/2)
ytmp=boxcentre(2)-(boxlength(2)/2)
ztmp=boxcentre(3)-(boxlength(3)/2)


num=0
do i=1,ndim
 ytmp=boxcentre(2)-(boxlength(2)/2)
 do j=1,ndim
  ztmp=boxcentre(3)-(boxlength(3)/2)
  do k=1,ndim
   num=num+1

   grid_ww_int(num,:,:,:)=.false.
   grid_ww_int(num,i,j,k)=.true.

! stuff to account for PBC wraparound for the interaction matrix
   if(i+1.le.ndim) then
    ui=i+1
   else
    ui=1
   end if
   if(i-1.ge.1) then
    li=i-1
   else
    li=ndim
   end if
   if(j+1.le.ndim) then
    uj=j+1
   else
    uj=1
   end if
   if(j-1.ge.1) then
    lj=j-1
   else
    lj=ndim
   end if
   if(k+1.le.ndim) then
    uk=k+1
   else
    uk=1
   end if
   if(k-1.ge.1) then
    lk=k-1
   else
    lk=ndim
   end if

! generating interaction matrix
   grid_ww_int(num,ui,uj,uk)=.true.
   grid_ww_int(num,ui,uj,lk)=.true.
   grid_ww_int(num,ui,lj,uk)=.true.
   grid_ww_int(num,ui,lj,lk)=.true.
   grid_ww_int(num,li,uj,uk)=.true.
   grid_ww_int(num,li,uj,lk)=.true.
   grid_ww_int(num,li,lj,uk)=.true.
   grid_ww_int(num,li,lj,lk)=.true.
   grid_ww_int(num,ui,uj,k)=.true.
   grid_ww_int(num,ui,lj,k)=.true.
   grid_ww_int(num,li,uj,k)=.true.
   grid_ww_int(num,li,lj,k)=.true.
   grid_ww_int(num,ui,j,uk)=.true.
   grid_ww_int(num,ui,j,lk)=.true.
   grid_ww_int(num,li,j,uk)=.true.
   grid_ww_int(num,li,j,lk)=.true.
   grid_ww_int(num,i,uj,uk)=.true.
   grid_ww_int(num,i,uj,lk)=.true.
   grid_ww_int(num,i,lj,uk)=.true.
   grid_ww_int(num,i,lj,lk)=.true.
   grid_ww_int(num,ui,j,k)=.true.
   grid_ww_int(num,li,j,k)=.true.
   grid_ww_int(num,i,j,uk)=.true.
   grid_ww_int(num,i,j,lk)=.true.
   grid_ww_int(num,i,uj,k)=.true.
   grid_ww_int(num,i,lj,k)=.true.

   grid_ww(num)%x=xtmp
   grid_ww(num)%y=ytmp
   grid_ww(num)%z=ztmp
   grid_ww(num)%xend=xtmp+gridRc
   grid_ww(num)%yend=ytmp+gridRc
   grid_ww(num)%zend=ztmp+gridRc

   if (k .eq. ndim) then
    if (grid_ww(num)%zend .lt. (boxcentre(3)+(boxlength(3)/2))) then
     grid_ww(num)%zend=boxcentre(3)+(boxlength(3)/2)
    end if
   end if
   ztmp=ztmp+gridRc
  end do
  if (j .eq. ndim) then
   if(grid_ww(num)%yend .lt. (boxcentre(2)+(boxlength(2)/2))) then
    grid_ww(num)%yend=boxcentre(2)+(boxlength(2)/2)
   end if
  end if
  ytmp=ytmp+gridRc
 end do
 if (i .eq. ndim) then
  if(grid_ww(num)%xend .lt. (boxcentre(1)+(boxlength(1)/2))) then
   grid_ww(num)%xend=boxcentre(1)+(boxlength(1)/2)
  end if
 end if
 xtmp=xtmp+gridRc
end do

else
! we are using spherical boundaries, meaning that we now have some problems
! we need to know the largest and smallest values for all coordiantes to generate the grid
i = 1
i3 = i*3-3
xmax=x(i3+1)
ymax=x(i3+2)
zmax=x(i3+3)
xmin=xmax
ymin=ymax
zmin=zmax
! now iterate over all cooridnates to find min and max
! atom one is already read in, don't need to read it again
do i=2,natom
 i3=i*3-3
 if (xmin .gt. x(i3+1)) xmin=x(i3+1)
 if (ymin .gt. x(i3+2)) ymin=x(i3+2)
 if (zmin .gt. x(i3+3)) zmin=x(i3+3)
 if (xmax .lt. x(i3+1)) xmax=x(i3+1)
 if (ymax .lt. x(i3+2)) ymax=x(i3+2)
 if (zmax .lt. x(i3+3)) zmax=x(i3+3)
end do
! now add some buffer to each of them -> half of Rcpp
xmin=xmin-(gridRc/3)
ymin=ymin-(gridRc/3)
zmin=zmin-(gridRc/3)
xmax=xmax+(gridRc/3)
ymax=ymax+(gridRc/3)
zmax=zmax+(gridRc/3)
! the system of boxes will be cubic in any case, so now find the largest number of 
! grid spaces one any side
ndim=ceiling((xmax-xmin)/gridRc)
if (ceiling((ymax-ymin)/gridRc) .gt. ndim ) then
 ndim=ceiling((ymax-ymin)/gridRc)
else if (ceiling((zmax-zmin)/gridRc) .gt. ndim ) then
 ndim=(ceiling((zmax-zmin)/gridRc))
end if
ww_gridnum=ndim**3
ww_ndim=ndim
if (nwat .gt. 500) then
gridstor_ww = (nwat*8/ww_gridnum)+100
else
gridstor_ww = nwat
endif
call allocate_grid_ww
xtmp=xmin
ytmp=ymin
ztmp=zmin

num=1
do i=1,ndim
 ytmp = ymin
 do j=1,ndim
  ztmp = zmin
  do k=1,ndim
! now fill the interaction matrix
! no wraparound in this case, so be careful
! when there would be wraparound, just point to my own cell again

grid_ww_int(num,:,:,:)=.false.
grid_ww_int(num,i,j,k)=.true.

   if(i+1.le.ndim) then
    ui=i+1
   else
    ui=i
   end if
   if(i-1.ge.1) then
    li=i-1
   else
    li=i
   end if
   if(j+1.le.ndim) then
    uj=j+1
   else
    uj=j
   end if
   if(j-1.ge.1) then
    lj=j-1
   else
    lj=j
   end if
   if(k+1.le.ndim) then
    uk=k+1
   else
    uk=k
   end if
   if(k-1.ge.1) then
    lk=k-1
   else
    lk=k
   end if

   grid_ww_int(num,ui,uj,uk)=.true.
   grid_ww_int(num,ui,uj,lk)=.true.
   grid_ww_int(num,ui,lj,uk)=.true.
   grid_ww_int(num,ui,lj,lk)=.true.
   grid_ww_int(num,li,uj,uk)=.true.
   grid_ww_int(num,li,uj,lk)=.true.
   grid_ww_int(num,li,lj,uk)=.true.
   grid_ww_int(num,li,lj,lk)=.true.
   grid_ww_int(num,ui,uj,k)=.true.
   grid_ww_int(num,ui,lj,k)=.true.
   grid_ww_int(num,li,uj,k)=.true.
   grid_ww_int(num,li,lj,k)=.true.
   grid_ww_int(num,ui,j,uk)=.true.
   grid_ww_int(num,ui,j,lk)=.true.
   grid_ww_int(num,li,j,uk)=.true.
   grid_ww_int(num,li,j,lk)=.true.
   grid_ww_int(num,i,uj,uk)=.true.
   grid_ww_int(num,i,uj,lk)=.true.
   grid_ww_int(num,i,lj,uk)=.true.
   grid_ww_int(num,i,lj,lk)=.true.
   grid_ww_int(num,ui,j,k)=.true.
   grid_ww_int(num,li,j,k)=.true.
   grid_ww_int(num,i,j,uk)=.true.
   grid_ww_int(num,i,j,lk)=.true.
   grid_ww_int(num,i,uj,k)=.true.
   grid_ww_int(num,i,lj,k)=.true.

   grid_ww(num)%x=xtmp
   grid_ww(num)%y=ytmp
   grid_ww(num)%z=ztmp

   grid_ww(num)%xend=xtmp+gridRc
   grid_ww(num)%yend=ytmp+gridRc
   grid_ww(num)%zend=ztmp+gridRc

   num=num+1
   ztmp=ztmp+gridRc
  end do
  ytmp=ytmp+gridRc
 end do
 xtmp=xtmp+gridRc
end do

end if
end subroutine create_grid_ww
#endif

!-----------------------------------------------------------------------


subroutine open_files
! --> restart file (2)
if(restart) then
open (unit=2, file=restart_file, status='old', form='unformatted', action='read', err=2)
end if

! --> final coords (3)
open (unit=3, file=xfin_file, status='unknown', form='unformatted', action='write', err=3)

! --> energy output file (11)
if ( iene_cycle .gt. 0 ) then
open (unit=11, file=ene_file, status='unknown', form='unformatted', action='write', err=11)
end if

! --> external file for implicit position restraints (12)
if ( implicit_rstr_from_file .eq. 1 ) then
open (unit=12, file=exrstr_file, status='old', form='unformatted', action='read', err=12)
end if


return

! crude error handling
2 call die('error opening restart file.')
3 call die('error opening final coordinates file.')
11 call die('error opening energy output file.')
12 call die('error opening position restraints file.')

end subroutine open_files

!-----------------------------------------------------------------------

subroutine gauss (am,sd,v,ig)
! arguments
real(kind=prec)					::	am,sd,v
integer					::	ig

! local variables
integer					::	i
real(kind=prec)					::	a,y

a=zero
do i=1,12
y=randm(ig)
a=a+y
end do
v=(a-6.0_prec)*sd+am
end subroutine gauss

!-----------------------------------------------------------------------

subroutine q_gauss (am,sd,v,ig)
! same as gauss, but returns a qr_vec object instead of single real number
! actually just calls gauss three times :D
! arguments
real(kind=prec)         :: am,sd
integer                 :: ig
TYPE(qr_vec)            :: v
! locals

call gauss(am,sd,v%x,ig)
call gauss(am,sd,v%y,ig)
call gauss(am,sd,v%z,ig)
end subroutine q_gauss


!-----------------------------------------------------------------------
subroutine get_fep
! local variables
character(len=200)			::	libtext
character(len=2)			::	qaname
integer					::	i,j,k,iat
!temp. array for reallocating long-range exclusion list
integer(AI), pointer	::	tempexlong(:,:)

! --- # states, # q-atoms
if(.not. qatom_load_atoms(fep_file)) then
        call die('failure to load Q-atoms from FEP file.')
end if

! set flags
do i=1,nqat
        if(iqseq(i) > 0 .and. iqseq(i) <= nat_solute)  then
                iqatom(iqseq(i)) = i
        else if(iqseq(i) == 0) then
                write(*,10) i
        else
                write(*,20) i, iqseq(i)
                call die('invalid q-atom data')
        end if
end do
10	format('>>> WARNING: Q-atom no. ',i2,' is not associated with a topology atom.')
20	format('>>>>> ERROR: Q-atom no. ',i2,' has invalid topology number ',i5)
!allocate memory for qatom charges
allocate(qcrg(nqat,nstates), stat=alloc_status)
call check_alloc('Qatom charges')

! --- copy topology charges

do i=1,nqat
        do j=1,nstates
                qcrg(i,j)=crg(iqseq(i))
        end do
end do

!initialize softcore lookup array
allocate (sc_lookup(nqat,natyps+nqat,nstates))
sc_lookup(:,:,:)=zero

!load rest of fep file
if(.not. qatom_load_fep(fep_file)) then
        call die('failure to load FEP file.')
end if

!Adapt LJ parameters to topology
!If arithmetic combination rule take sqrt(epsilon) now
if (qvdw_flag .and. ivdw_rule .eq. 2 ) then
        qbvdw(:,1) = sqrt( qbvdw(:,1) )
        qbvdw(:,3) = sqrt( qbvdw(:,3) )
end if

!remove redefined bonded interactions from topology
if(nqbond > 0 .or. nqangle > 0 .or. nqtor > 0 .or. nqimp > 0 ) then
        write(*,*)
        call centered_heading('Removing redefined interactions from topology','-')
230		format('type',t10,' atom1 atom2 atom3 atom4')
        write(*,230)
231		format(a,t10,4i6)
        !remove bonds that were redefined
        do i=1,nbonds
                do j=1,nqbond
                        if ( (bnd(i)%i==qbnd(j)%i .and. bnd(i)%j==qbnd(j)%j) .or. &
                                (bnd(i)%i==qbnd(j)%j .and. bnd(i)%j==qbnd(j)%i) ) then
                                bnd(i)%cod = 0
                                write (*,231) 'bond',bnd(i)%i,bnd(i)%j
                        end if
                end do
        end do

        !remove angles that were redefined
        do i=1,nangles
                do j=1,nqangle
                        if((ang(i)%i.eq.qang(j)%i .and. ang(i)%j.eq.qang(j)%j .and. &
                                ang(i)%k.eq.qang(j)%k)                          .or. &
                                (ang(i)%i.eq.qang(j)%k .and. ang(i)%j.eq.qang(j)%j .and. &
                                ang(i)%k.eq.qang(j)%i) )                         then

                                ang(i)%cod = 0
                                write (*,231) 'angle',ang(i)%i,ang(i)%j,ang(i)%k
                        end if
                end do
        end do

        !remove torsions that were redefined
        do i=1,ntors
                do j=1,nqtor
			if(( (tor(i)%i.eq.qtor(j)%i .and. tor(i)%j.eq.qtor(j)%j .and. &
				tor(i)%k.eq.qtor(j)%k .and. tor(i)%l.eq.qtor(j)%l) .or. &
				(tor(i)%i.eq.qtor(j)%l .and. tor(i)%j.eq.qtor(j)%k .and. &
				tor(i)%k.eq.qtor(j)%j .and. tor(i)%l.eq.qtor(j)%i )) .and. &
				tor(i)%cod .ne. 0) then
                                tor(i)%cod = 0
                                write (*,231) 'torsion', tor(i)%i,tor(i)%j,tor(i)%k,tor(i)%l
						end if
                end do
        end do


        !remove impropers that were redefined
        select case(ff_type)
        case(FF_CHARMM) !special code for CHARMM
                do i=1,nimps
                        do j=1,nqimp
				if ( ((	(imp(i)%i.eq.qimp(j)%i) .or. &
					(imp(i)%i.eq.qimp(j)%l) .or. &
					(imp(i)%l.eq.qimp(j)%i)) .and. &
				     (	(imp(i)%j.eq.qimp(j)%l) .or. &
					(imp(i)%j.eq.qimp(j)%i) .or. &
					(imp(i)%j.eq.qimp(j)%j) .or. &
					(imp(i)%j.eq.qimp(j)%k)) .and. &
				     (	(imp(i)%k.eq.qimp(j)%i) .or. &
					(imp(i)%k.eq.qimp(j)%j) .or. &
					(imp(i)%k.eq.qimp(j)%k) .or. &
					(imp(i)%k.eq.qimp(j)%l))) .and. &
				     (	imp(i)%cod .ne. 0)) then
                                        imp(i)%cod = 0
                                        write (*,231) &
                                        'improper',imp(i)%i,imp(i)%j,imp(i)%k,imp(i)%l
                                end if
                        end do
                end do

        case default
        do i=1,nimps
                do j=1,nqimp
                        if(( (imp(i)%i.eq.qimp(j)%i) .and. (imp(i)%j.eq.qimp(j)%j) .and. &
                                (imp(i)%k.eq.qimp(j)%k) .and. (imp(i)%l.eq.qimp(j)%l)) .or. &
                                ((imp(i)%i.eq.qimp(j)%k) .and. (imp(i)%j.eq.qimp(j)%j) .and. &
                                (imp(i)%k.eq.qimp(j)%i) .and. (imp(i)%l.eq.qimp(j)%l)) .or. &
                                ((imp(i)%i.eq.qimp(j)%k) .and. (imp(i)%j.eq.qimp(j)%j) .and. &
                                (imp(i)%k.eq.qimp(j)%l) .and. (imp(i)%l.eq.qimp(j)%i)) .or. &
                                ((imp(i)%i.eq.qimp(j)%l) .and. (imp(i)%j.eq.qimp(j)%j) .and. &
                                (imp(i)%k.eq.qimp(j)%i) .and. (imp(i)%l.eq.qimp(j)%k)) .or. &
                                ((imp(i)%i.eq.qimp(j)%l) .and. (imp(i)%j.eq.qimp(j)%j) .and. &
                                (imp(i)%k.eq.qimp(j)%k) .and. (imp(i)%l.eq.qimp(j)%i))) then
                                imp(i)%cod = 0
                                write(*,231)'improper',imp(i)%i,imp(i)%j,imp(i)%k,imp(i)%l
                        end if
                end do
        end do
        end select
end if

!check special exclusions
!modify exclusion lists to inclue special exclusions between Q and non-Q
if(nexspec > 0) then
        allocate(tempexlong(2,nexlong+nexspec))
        tempexlong(:, 1:nexlong) = listexlong(:, 1:nexlong)
        deallocate(listexlong)
        listexlong => tempexlong
end if

do k = 1, nexspec
        i = exspec(k)%i
        j = exspec(k)%j
        if(i < 1 .or. i > nat_pro .or. j < 1 .or. j > nat_pro) then
                write(*, 592) k, i, j
                call die('invalid special exclusion data')
        end if
        !if one or more non-Q-atoms modify exclusion lists
        if(iqatom(i)==0 .or. iqatom(j)==0) then
                !With non-Q-atoms involved only accept all or no states
                if(any(exspec(k)%flag(1:nstates))) then
                        if(.not. all(exspec(k)%flag(1:nstates))) then
                                write(*,594) k
                                call die('invalid special exclusion data')
                        else !exlcude in all states
                                if(abs(j-i) <= max_nbr_range) then
                                        if(i < j) then
                                                listex(j-i,i) = .true.
                                        else
                                                listex(i-j,j) = .true.
                                        end if
                                else
                                        nexlong = nexlong + 1
                                        listexlong(1, nexlong) = i
                                        listexlong(2, nexlong) = j
                                end if
                        end if
                end if
        end if
end do
592	format('>>>>> ERROR: Special exclusion pair ',i2,' (',i5,1x,i5,') is invalid')
594	format('>>>>> ERROR: Non-Q-atom special excl. pair ',i2,' must be on in all or no states')
end subroutine get_fep

!-----------------------------------------------------------------------

subroutine get_fname (text,length,filnam)
! arguments
character(len=200)			::	text,filnam
integer					::	length
! local variables


integer					::	i

length=200
do i=1,200
if ( text(i:i) .eq. ' ' ) then
length=i-1
goto 10
end if
end do
10 filnam(1:length)=text(1:length)

end subroutine get_fname

!-----------------------------------------------------------------------

#if defined (USE_MPI)
!Defines and allocates variables needed in the md-calculations
!The node initiation is written for AI = 4. If changes are made to any size in
! sizes.f90 the MPI-code must be changed accordingly. It is not dynamically
! implemented yet.
subroutine init_nodes
!
! initialise slave nodes, sending to slaves:
!
! variables:
!  natom,nwat,nsteps,use_LRF,NBcycle,crg_ow,crg_hw,Rcpp,Rcww,Rcpw,Rcq,xpcent
!  nat_solute,ncgp,ncgp_solute,ivdw_rule,iuse_switch_atom,el14_scale,n14long
!  nexlong,natyps,nljtyp,rexcl_o,nstates,nqat,qvdw_flag,nqlib,RcLRF,
!  use_PBC, qswitch, nmol, nat_pro
!
! arrays:
!  x,v,iqatom,ljcod,qconn,iwhich_cgp,lrf,excl,iac,crg,cgpatom,cgp,iaclib
!  list14,listex,list14long,listexlong,iqseq,qiac,qcrg,qavdw,qbvdw,EQ(:)%lambda,
!  boxlength, inv_boxl, boxcentre, sc_lookup 
!

integer, parameter			:: vars = 40    !increment this var when adding data to broadcast in batch 1
integer				   	:: blockcnt(vars), ftype(vars)
integer(kind=MPI_ADDRESS_KIND)			   	:: fdisp(vars)
integer					:: mpitype_batch,mpitype_batch2
integer					:: nat3,j,jj,ii
!temp for shake
integer,allocatable                     :: shakeconsttmp(:)
integer,allocatable                     :: shakebondi(:)
integer,allocatable                     :: shakebondj(:)
integer                                 :: shakenumconst
real(kind=prec),allocatable             :: shakedist(:)
real(kind=prec), allocatable		:: temp_lambda(:)
integer, parameter                      :: maxint=2147483647
real(kind=prec), parameter              :: maxreal=1E35_prec
integer  :: MPI_AI_INTEGER, MPI_TINY_INTEGER, i_loop

!external MPI_Address
!external MPI_Bcast

!**********
!2002-11-28 
!MN-> This will work with new implementations of MPI standard >= 2
!The MPI library at PDC does not support these definitions when I tried to use them.
!Using these routines will allow a change made to the sizes in sizes.f90 to
! affect the mpi. Without them the variables below marked (AI) and (TINY) will have to 
! be changed manually.
!When using this part make sure the vars marked with comments (AI) and (TINY) are 
! changed to MPI_AI_INTEGER and MPI_TINY_INTEGER.

!external MPI_Type_Create_F90_Integer
!external MPI_SizeOf

!Define data types
! This is wrong, the 1:st param is "Precision, in decimal digits", not bits
!call MPI_Type_Create_F90_Integer((8*AI-1),MPI_AI_INTEGER,ierr)
!call MPI_Type_Create_F90_Integer((8*TINY-1),MPI_TINY_INTEGER,ierr)
!To check the size in bytes of the new types use
!call MPI_SizeOf(MPI_AI_INTEGER,size,ierr)
!call MPI_SizeOf(MPI_TINY_INTEGER,size,ierr)



if (nodeid .eq. 0) call centered_heading('Distributing data to slave nodes', '-')
!***************************

! --- mandatory data, first batch ---

if (nodeid .eq. 0) write (*,'(80a)') 'MD data, first batch'

! initialise custom MPI data type; default is an integer scalar
! The variables below are sorted in modular order.
! Make sure to add next variable to the right module and 
! add +1 to 'vars'.
blockcnt(:) = 1
ftype(:) = MPI_INTEGER

! run control constants: natom, nwat, nsteps, NBmethod, NBcycle
call MPI_Bcast(natom, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast natom')
call MPI_Bcast(nwat, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast nwat')
call MPI_Bcast(nsteps, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast nsteps')
call MPI_Bcast(use_LRF, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast use_LRF')
call MPI_Bcast(NBcycle, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast NBcycle')
call MPI_Bcast(iene_cycle,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast iene_cycle')
!some more new stuff for the thermostats/integrators
!now done in parallel
call MPI_Bcast(friction,1,MPI_REAL8,0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast friction')
call MPI_Bcast(thermostat,1,MPI_INTEGER,0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast thermostat')
call MPI_Bcast(dt,1,MPI_REAL8,0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast dt')
call MPI_Bcast(dt2,1,MPI_REAL8,0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast dt2')
call MPI_Bcast(Temp0,1,MPI_REAL8,0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast Temp0')
call MPI_Bcast(iseed,1,MPI_INTEGER,0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast dt2')
call MPI_Bcast(Temp0,1,MPI_REAL8,0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast iseed')

! water parameters: chg_solv array and solv_atom (used by nonbond_?w)
! also data type and internal interaction array
ftype(1) = MPI_INTEGER4
ftype(2) = MPI_INTEGER4
ftype(3) = MPI_REAL8
ftype(4) = MPI_REAL8
ftype(5) = MPI_REAL8
blockcnt(1) = 1
blockcnt(2) = 1
blockcnt(3) = 1
blockcnt(4) = 1
blockcnt(5) = 1
fdisp(1) = 0
fdisp(2) = 4
fdisp(3) = 4+4
fdisp(4) = 4+4+8
fdisp(5) = 4+4+8+8
call MPI_Type_create_struct(5, blockcnt, fdisp, ftype, mpitype_batch_solv_int, ierr)
if (ierr .ne. 0) call die('init_nodes MPI_Type_create_struct solv_int')
call MPI_Type_commit(mpitype_batch_solv_int, ierr)
if (ierr .ne. 0) call die('init_nodes MPI_Type_commit solv_int')

if (nwat .gt. 0 ) then
call MPI_Bcast(solvent_type,1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
if (ierr .ne. 0 ) call die('init_nodes/MPI_Bcast solvent_type')
call MPI_Bcast(solv_atom,1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
if (ierr .ne. 0 ) call die('init_nodes/MPI_Bcast solv_atom')
call MPI_Bcast(num_solv_int,1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
if (ierr .ne. 0 ) call die('init_nodes/MPI_Bcast num_solv_int')

if (nodeid .ne. 0 ) then
allocate(chg_solv(solv_atom),stat=alloc_status)
call check_alloc('chg_solv')
allocate(aLJ_solv(solv_atom,3),stat=alloc_status)
call check_alloc('aLJ_solv')
allocate(bLJ_solv(solv_atom,3),stat=alloc_status)
call check_alloc('bLJ_solv')
allocate(nonbnd_solv_int(num_solv_int),stat=alloc_status)
call check_alloc('nonbnd_solv_int')
end if
call MPI_Bcast(chg_solv,solv_atom,MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast chg_solv array')
call MPI_Bcast(aLJ_solv,3*solv_atom,MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast chg_solv array')
call MPI_Bcast(bLJ_solv,3*solv_atom,MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast chg_solv array')
call MPI_Bcast(nonbnd_solv_int,num_solv_int,mpitype_batch_solv_int, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast nonbnd_solv_int struct')
end if

! cutoffs: Rcpp, Rcww, Rcpw, Rcq, RcLRF (used by pair list generating functions)
call MPI_Bcast(Rcpp, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast Rcpp')
call MPI_Bcast(Rcww, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast Rcww')
call MPI_Bcast(Rcpw, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast Rcpw')
call MPI_Bcast(Rcq, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast Rcq')
call MPI_Bcast(RcLRF, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast RcLRF')




!Periodic Boudary Condition
call MPI_Bcast(use_PBC, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast use_PBC')
call MPI_Bcast(boxcentre, 3, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast boxcentre')
call MPI_Bcast(boxlength, 3, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast boxlength')
call MPI_Bcast(inv_boxl, 3, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast inv_boxl')
call MPI_Bcast(qswitch, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast qswitch')
call MPI_Bcast(constant_pressure, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast constant_pressure')
call MPI_Bcast(ivolume_cycle, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast ivolume_cycle')
call MPI_Bcast(rigid_box_centre, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast rigid_box_centre')
call MPI_Bcast(put_solvent_back_in_box, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast put_solvent_back_in_box')
call MPI_Bcast(put_solute_back_in_box, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast put_solute_back_in_box')

! xpcent            from TOPO, needed for listgeneration
call MPI_Bcast(xpcent, 3, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast xpcent')


! a bunch of vars from the TOPO module
call MPI_Bcast(nat_solute, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast nat_solute')
call MPI_Bcast(ncgp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast ncgp')
call MPI_Bcast(ncgp_solute, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast ncgp_solute')
call MPI_Bcast(ivdw_rule, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast ivdw_rule')
call MPI_Bcast(iuse_switch_atom, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast iuse_switch_atom')
call MPI_Bcast(el14_scale, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast el14_scale')
call MPI_Bcast(n14long, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast n14long')
call MPI_Bcast(nexlong, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast nexlong')
call MPI_Bcast(natyps, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast natyps')
call MPI_Bcast(rexcl_o, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast rexcl')
call MPI_Bcast(nmol, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast nmol')
call MPI_Bcast(nres_solute, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast nres_solute')
call MPI_Bcast(nat_pro, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast nat_pro')

! water pol stuff
call MPI_Bcast(wpol_restr, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast wpol_restr')
call MPI_Bcast(n_max_insh, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast n_max_insh')
call MPI_Bcast(nwpolr_shell, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast nwpolr_shell')
call MPI_Bcast(fk_pshell, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast fk_pshell')
call MPI_Bcast(fk_wsphere, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast fk_wsphere')
call MPI_Bcast(Dwmz, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast Dwmz')
call MPI_Bcast(awmz, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast awmz')
call MPI_Bcast(fkwpol, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast fkwpol')
call MPI_Bcast(rwat, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast rwat')



#ifdef USE_GRID
!first few grid variables
call MPI_Bcast(gridstor_pp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast gridstor_pp')
call MPI_Bcast(gridstor_pw, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast gridstor_pw')
call MPI_Bcast(gridstor_ww, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast gridstor_ww')

! water grid lookup
if (nwat .gt. 0 ) then
if (nodeid .ne. 0 ) then
allocate(ww_igrid(nwat),stat=alloc_status)
call check_alloc('ww_igrid slave nodes')
allocate(pw_igrid(ncgp_solute),stat=alloc_status)
call check_alloc('pw_igrid slave nodes')
end if
call MPI_Bcast(ww_igrid,nwat,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast ww_igrid')
call MPI_Bcast(pw_igrid,ncgp_solute,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast pw_igrid')
end if
! protein grid lookup
if (ncgp_solute .gt. 0) then
if (nodeid .ne. 0 ) then
allocate(pp_igrid(ncgp_solute),stat=alloc_status)
call check_alloc('pp_igrid slave nodes')
end if
call MPI_Bcast(pp_igrid,ncgp_solute,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast pp_igrid')
end if ! ncgp_solute .gt. 0
#endif

!vars from QATOM
call MPI_Bcast(nstates, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast nstates')
call MPI_Bcast(use_excluded_groups,1,MPI_LOGICAL, 0,MPI_COMM_WORLD,ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast use_excluded_groups')

if (use_excluded_groups) then
	call MPI_Bcast(ngroups_gc, 1, MPI_INTEGER, 0,MPI_COMM_WORLD,ierr)
	if (ierr .ne. 0) call die('init_nodes/MPI_Bcast ngroups_gc')
end if
call MPI_Bcast(nqat, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast nqat')
call MPI_Bcast(qvdw_flag, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast qvdw_flag')
call MPI_Bcast(nqlib, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast nqlib')

if (use_excluded_groups) then
        allocate(tempmask(nat_pro,ngroups_gc),stat=alloc_status)
        call check_alloc('temp mask arrays')
        if (nodeid .ne. 0) then
                allocate(ST_gc(ngroups_gc),stat=alloc_status)
                call check_alloc('node gc array')
        end if
        if (nodeid.eq.0) then
                do jj=1,ngroups_gc
                tempmask(1:nat_pro,jj)=ST_gc(jj)%gcmask%mask(1:nat_pro)
                end do
        end if
! Wait first for Master here
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_Bcast(tempmask,nat_pro*ngroups_gc,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
        if (ierr .ne. 0) call die('init_nodes/MPI_Bcast gc atom masks')
        if (nodeid.ne.0) then
                do jj=1,ngroups_gc
                call mask_initialize(ST_gc(jj)%gcmask)
                ST_gc(jj)%gcmask%mask(1:nat_pro)=tempmask(1:nat_pro,jj)
                end do
        end if
        deallocate(tempmask,stat=alloc_status)
	call MPI_Bcast(ST_gc(1:ngroups_gc)%caltype,ngroups_gc,MPI_INTEGER, 0,MPI_COMM_WORLD,ierr)
	if (ierr .ne. 0) call die('init_nodes/MPI_Bcast group contrib calc types')
end if

call MPI_Bcast(use_qcp,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
if (ierr .ne.0) call die('init_nodes/MPI_BCast use_qcp')
! send all qcp stuff we might be needing
if (use_qcp) then
call MPI_Bcast(qcp_size,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast qcp_size')
call MPI_Bcast(qcp_level,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast qcp_level')
call MPI_Bcast(qcp_steps,2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast qcp_steps')
call MPI_Bcast(qcp_atnum,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast qcp_atnum')
if (nodeid .ne. 0) then
allocate(qcp_atom(qcp_atnum))
end if
call MPI_Bcast(qcp_atom,qcp_atnum,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast qcp_atom')
end if


call MPI_Bcast(xwcent,3,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast wxcent')

! --- MD data, second batch ---

if (nodeid .eq. 0) write (*,'(80a)') 'MD data, second batch'

! allocate arrays
if (nodeid .ne. 0) then
 call allocate_natom_arrays
 if (thermostat == NOSEHOOVER) then
  call allocate_nhchain_arrays
 end if
end if

! broadcast x, v and winv
nat3 = 3*natom

call MPI_Bcast(x, nat3, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast x')
call MPI_Bcast(v, nat3, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast v')


call MPI_Bcast(winv,natom,MPI_REAL8,0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast winv')

!Broadcast iqatom
call MPI_Bcast(iqatom, natom, MPI_INTEGER2, 0, MPI_COMM_WORLD, ierr) !(TINY)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast iqatom')

call MPI_Bcast(num_atyp,1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast num_atyp')

!Broadcast ljcod
if (nodeid .ne. 0) then
allocate(ljcod(num_atyp,num_atyp))
end if
call MPI_Bcast(ljcod, size(ljcod), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast ljcod')

!Broadcast qconn(nstates,nat_solute, nqat)
if (nodeid .ne. 0) then
allocate(qconn(nstates,nat_solute, nqat),stat=alloc_status)
call check_alloc('qconn')
end if
call MPI_Bcast(qconn, size(qconn), MPI_INTEGER2, 0, MPI_COMM_WORLD, ierr) ! (TINY)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast qconn')


! --- Periodic boundary condition data ---


if (nodeid.ne.0) then
allocate(mvd_mol(nmol))
end if


! --- lrf data ---

if (use_LRF) then
! lrf stuff

if (nodeid .eq. 0) write (*,'(80a)') 'lrf data'

! allocate arrays
if (nodeid .ne. 0) call allocate_lrf_arrays

!MPI_INTEGER4 is used instead of MPI_AI_INTEGER
!Change to mpi_type_create, see note above or note2 in sizes.f90
! iwhich_cgp
call MPI_Bcast(iwhich_cgp, natom, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr) !(AI)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast lrf parameters')


! lrf
! we now keep the data type of the lrf for later usage
ftype(:) = MPI_REAL8
blockcnt(1) = 3					! real(kind=prec) cgp_cent(3)
fdisp(1) = 0
blockcnt(2) = 1					! real(kind=prec) phi0
fdisp(2) = 3*8
blockcnt(3) = 3					! real(kind=prec) phi1(3)
fdisp(3) = 3*8 + 8
blockcnt(4) = 9					! real(kind=prec) phi2(9)
fdisp(4) = 3*8 + 8 + 3*8
blockcnt(5) = 27				! real(kind=prec) phi3(27)
fdisp(5) = 3*8 + 8 + 3*8 + 9*8
call MPI_Type_create_struct(5, blockcnt, fdisp, ftype, mpitype_batch_lrf, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Type_create_struct')
call MPI_Type_commit(mpitype_batch_lrf, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Type_commit')
call MPI_Bcast(lrf, ncgp, mpitype_batch_lrf, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast lrf parameters')
call MPI_Op_create(lrf_add,.false.,mpi_lrf_add,ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Op_create mpi_lrf_add')
call MPI_Op_create(lrf_cgp_rep,.false.,mpi_lrf_cgp_rep,ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Op_create mpi_lrf_cgp_rep')

end if !(use_LRF)

! --- data from the TOPO module ---

if (nodeid .eq. 0) write (*,'(80a)') 'TOPO data'

! allocate topology arrays
if (nodeid .ne. 0) then
! don't allocate memory for stuff we don't need
! these array size variables are actually used
max_cgp=ncgp
max_atyps = natyps
max_14long = n14long
max_exlong = nexlong
max_atom = natom

call topo_allocate_atom(alloc_status)
call check_alloc('topology arrays')
call topo_allocate_potential(alloc_status)
call check_alloc('topology arrays')
allocate(istart_mol(nmol+1), &
stat=alloc_status)
call check_alloc('topology arrays')
end if

! broadcast excl
call MPI_Bcast(excl, natom, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast excl')
! broadcast istart_mol
call MPI_Bcast(istart_mol, nmol+1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast istart_mol')

! Bcast iac, crg and cgpatom 
call MPI_Bcast(iac, natom, MPI_INTEGER2, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast iac')
call MPI_Bcast(crg, natom, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast crg')
call MPI_Bcast(cgpatom, natom, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr) !(AI)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast cgpatom')

call MPI_Bcast(xtop,nat3,MPI_REAL8,0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast xtop')
call MPI_Bcast(shell,natom, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast shell')


! cgp
!Use MPI_Type_create_struct here too
ftype(:) = MPI_INTEGER4 !(AI)
blockcnt(:) = 1
fdisp(1) = 0				! integer(AI) iswitch
fdisp(2) = AI				! integer(AI) first
fdisp(3) = AI + AI			! integer(AI) last
call MPI_Type_create_struct(3, blockcnt, fdisp, ftype, mpitype_batch, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Type_create_struct')
call MPI_Type_commit(mpitype_batch, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Type_commit')
call MPI_Bcast(cgp, ncgp, mpitype_batch, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast cgp')
call MPI_Type_free(mpitype_batch, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Type_free')

! iaclib
ftype(:) = MPI_REAL8
blockcnt(1) = 1					! real(kind=prec) mass
fdisp(1) = 0
blockcnt(2) = nljtyp				! real(kind=prec) avdw(nljtyp)
fdisp(2) = 8
blockcnt(3) = nljtyp				! real(kind=prec) bvdw(nljtyp)
fdisp(3) = 8 + 8*nljtyp
call MPI_Type_create_struct(3, blockcnt, fdisp, ftype, mpitype_batch, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Type_create_struct')
call MPI_Type_commit(mpitype_batch, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Type_commit')
call MPI_Bcast(iaclib, max_atyps, mpitype_batch, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast iaclib')
call MPI_Type_free(mpitype_batch, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Type_free')

! list14 and listex share the same format: logical listxx(max_nbr_range,max_atom)
call MPI_Bcast(list14, size(list14), MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast list14')
call MPI_Bcast(listex, size(listex), MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast listex')

! list14long and listexlong share the same format: integer(AI) listxxlong(2,max_nxxlong)
call MPI_Bcast(list14long, 2*n14long, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr) !(AI)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast list14long')
call MPI_Bcast(listexlong, 2*nexlong, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast listexlong')

#ifdef USE_GRID

! data for the different grids
call MPI_Bcast(pp_gridnum,1,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast pp_gridnum')
call MPI_Bcast(pp_ndim,1,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast pp_ndim')
call MPI_Bcast(pw_gridnum,1,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast pw_gridnum')
call MPI_Bcast(pw_ndim,1,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast pw_ndim')
call MPI_Bcast(ww_gridnum,1,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast ww_gridnum')
call MPI_Bcast(ww_ndim,1,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast ww_ndim')

if (nodeid.ne.0) then
!prepare grids on the different nodes using the data from master
if (ncgp_solute .gt. 0) then
allocate(grid_pp(pp_gridnum),stat=alloc_status)
call check_alloc('MPI pp_grid')
allocate(grid_pp_grp(pp_gridnum,gridstor_pp),stat=alloc_status)
call check_alloc('MPI grid_pp_grp')
allocate(grid_pp_ngrp(pp_gridnum),stat=alloc_status)
call check_alloc('MPI grid_pp_ngrp')
allocate(grid_pp_int(pp_gridnum,pp_ndim,pp_ndim,pp_ndim),stat=alloc_status)
call check_alloc('MPI grid_pp_int')
end if ! ncgp_solute .gt. 0
if (nwat .gt. 0) then
allocate(grid_pw(pw_gridnum),stat=alloc_status)
call check_alloc('MPI pw_grid')
allocate(grid_ww(ww_gridnum),stat=alloc_status)
call check_alloc('MPI ww_grid')
allocate(grid_pw_grp(pw_gridnum,gridstor_pw),stat=alloc_status)
call check_alloc('MPI grid_pw_grp')
allocate(grid_pw_ngrp(pw_gridnum),stat=alloc_status)
call check_alloc('MPI grid_pw_ngrp')
allocate(grid_pw_int(pw_gridnum,pw_ndim,pw_ndim,pw_ndim),stat=alloc_status)
call check_alloc('MPI grid_pw_int')
allocate(grid_ww_grp(ww_gridnum,gridstor_ww),stat=alloc_status)
call check_alloc('MPI grid_ww_grp')
allocate(grid_ww_ngrp(ww_gridnum),stat=alloc_status)
call check_alloc('MPI grid_ww_ngrp')
allocate(grid_ww_int(ww_gridnum,ww_ndim,ww_ndim,ww_ndim),stat=alloc_status)
call check_alloc('MPI grid_ww_int')
end if !nwat > 0
end if ! nodeid .ne. 0

!make first process wait until nodes have allocated grid arrays
call MPI_BARRIER(MPI_COMM_WORLD, ierr)

! also use make struct for grids
! loop based assignment not working
ftype(1) = MPI_REAL8            ! real(kind=prec)       x
ftype(2) = MPI_REAL8            ! real(kind=prec)       y
ftype(3) = MPI_REAL8            ! real(kind=prec)       z
ftype(4) = MPI_REAL8            ! real(kind=prec)       xend
ftype(5) = MPI_REAL8            ! real(kind=prec)       yend
ftype(6) = MPI_REAL8            ! real(kind=prec)       zend
blockcnt(1) = 1
blockcnt(2) = 1
blockcnt(3) = 1
blockcnt(4) = 1
blockcnt(5) = 1
blockcnt(6) = 1
fdisp(1) = 0     
fdisp(2) = 1*8
fdisp(3) = 2*8
fdisp(4) = 3*8
fdisp(5) = 4*8
fdisp(6) = 5*8
call MPI_Type_create_struct(6, blockcnt, fdisp, ftype, mpitype_batch, ierr)
if (ierr .ne. 0) call die('init_nodes MPI_Type_create_struct grid')
call MPI_Type_commit(mpitype_batch, ierr)
if (ierr .ne. 0) call die('init_nodes MPI_Type_commit grid')
if (ncgp_solute .gt. 0) then
call MPI_Bcast(grid_pp, pp_gridnum, mpitype_batch, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes MPI_Bcast grid_pp')
end if
if (nwat .gt. 0 ) then
call MPI_Bcast(grid_pw, pw_gridnum, mpitype_batch, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes MPI_Bcast grid_pw')
call MPI_Bcast(grid_pw, pw_gridnum, mpitype_batch, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes MPI_Bcast grid_ww')
end if ! nwat > 0
call MPI_Type_free(mpitype_batch, ierr)
if (ierr .ne. 0) call die('init_nodes MPI_Type_free grid')

! now the rest of the individual arrays
! only need to send interaction matrix, rest is not constant and updated in
! populate grids
if (ncgp_solute .gt. 0 ) then
call MPI_Bcast(grid_pp_int,pp_gridnum*pp_ndim**3,MPI_LOGICAL,0, MPI_COMM_WORLD,ierr)
if (ierr .ne. 0) call die('init_nodes MPI_Bcast grid_pp_int')
! make grid_pp_grp mpi data structure here
call MPI_Type_contiguous(gridstor_pp, MPI_INTEGER, mpitype_batch_ppgrid, ierr)
if (ierr .ne. 0) call die('init_nodes MPI_Type_contiguous gridstor_pp')

endif
if (nwat .gt. 0) then
call MPI_Bcast(grid_pw_int,pw_gridnum*pw_ndim**3,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
if (ierr .ne. 0) call die('init_nodes MPI_Bcast grid_pw_int')
call MPI_Bcast(grid_ww_int,ww_gridnum*ww_ndim**3,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
if (ierr .ne. 0) call die('init_nodes MPI_Bcast grid_ww_int')

call MPI_Type_contiguous(gridstor_pw, MPI_INTEGER, mpitype_batch_pwgrid, ierr)
if (ierr .ne. 0) call die('init_nodes MPI_Type_contiguous gridstor_pw')
call MPI_Type_contiguous(gridstor_ww, MPI_INTEGER, mpitype_batch_wwgrid, ierr)
if (ierr .ne. 0) call die('init_nodes MPI_Type_contiguous gridstor_ww')

endif

call MPI_Op_create(grid_add,.false.,mpi_grid_add,ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Op_create mpi_grid_add')


#endif

! --- data from the QATOM module ---

if (nodeid .eq. 0) write (*,'(80a)') 'QATOM data'

! allocate memory
if (nodeid .ne. 0) then
allocate(iqseq(nqat), &
qiac(nqat,nstates), &
qcrg(nqat,nstates), &
qavdw(nqlib,nljtyp), &
qbvdw(nqlib,nljtyp), &
EQ(nstates), &
sc_lookup(nqat,natyps+nqat,nstates), &
stat=alloc_status)
!if we use PBC, we allocate the storage for the old energies here
if( use_PBC ) then
allocate(old_EQ(nstates))
end if
call check_alloc('Q-atom arrays')
end if
!Broadcast sc_lookup(nqat,natyps+nqat,nstates)
if (nstates.ne.0) then
call MPI_Bcast(sc_lookup, size(sc_lookup), MPI_REAL8, 0, MPI_COMM_WORLD,ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast sc_lookup')
else
call MPI_Bcast(sc_lookup,  nstates, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast sc_lookup')
end if

! integer(AI) ::  iqseq(nqat)
!Change to mpi_type_create  (AI)
call MPI_Bcast(iqseq, nqat, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast iqseq')

!  integer ::  qiac(nqat,nstates)
if (nstates.ne.0) then
call MPI_Bcast(qiac, size(qiac), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast qiac')
else
call MPI_Bcast(qiac, nstates, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast qiac')
end if
! real(4) ::  qcrg(nqat,nstates)
if (nstates.ne.0) then
call MPI_Bcast(qcrg, size(qcrg), MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast qcrg')
else
call MPI_Bcast(qcrg, nstates, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast qcrg')
end if

if(qvdw_flag) then
!MN20030409-> Havn't tried with qvdw_flag == .true.
! qavdw and qbvdw share the same format: real(kind=prec) qxvdw(nqlib,nljtyp)
if (nstates.ne.0) then
call MPI_Bcast(qavdw, size(qavdw), MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast qavdw')
call MPI_Bcast(qbvdw, size(qbvdw), MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast qbvdw')
else
call MPI_Bcast(qavdw, nljtyp, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast qavdw')
call MPI_Bcast(qbvdw, nljtyp, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast qbvdw')
end if
end if
if (nstates .gt. 0) then
! Broadcast EQ(:)%lambda
allocate(temp_lambda(1:nstates), stat=alloc_status)
call check_alloc('Q-atom energy array')
if (nodeid .eq. 0) temp_lambda(1:nstates) = EQ(1:nstates)%lambda
call MPI_Bcast(temp_lambda, nstates, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast EQ%lambda')
if (nodeid .ne. 0) EQ(1:nstates)%lambda = temp_lambda(1:nstates)
deallocate(temp_lambda)
end if
!Final Bcast for the group contribution header data type
call MPI_Bcast(ene_header%arrays,1,MPI_INTEGER,0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast group contrib arrays')
call MPI_Bcast(ene_header%totresid,1,MPI_INTEGER,0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast group contrib totresid')

if (nodeid .ne. 0 ) then
allocate(ene_header%types(ene_header%arrays),ene_header%numres(ene_header%arrays),&
	ene_header%resid(ene_header%totresid),ene_header%gcnum(ene_header%arrays),&
	stat=alloc_status)
call check_alloc('Group contrib type arrays')
end if
call MPI_Bcast(ene_header%types,ene_header%arrays,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast group contrib arrays')
call MPI_Bcast(ene_header%numres,ene_header%arrays,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast group contrib arrays')
call MPI_Bcast(ene_header%gcnum,ene_header%arrays,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast group contrib arrays')
call MPI_Bcast(ene_header%resid,ene_header%totresid,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast group contrib arrays')

if (use_qcp .and. nodeid .ne. 0) then
allocate(qcp_EQ(nstates),stat=alloc_status)
call check_alloc('QCP node EQ arrays')
end if
!send the shake data here now, deleted old reference to this above
!shake now done on each node for own set of atoms
!PB 2015
if (nodeid .eq. 0) write (*,'(80a)') 'Shake control data'
call MPI_Bcast(shake_hydrogens,1,MPI_LOGICAL,0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast shake hydrogens')
call MPI_Bcast(shake_heavy,1,MPI_LOGICAL,0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast shake heavy')
call MPI_Bcast(shake_solvent,1,MPI_LOGICAL,0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast shake solvent')
call MPI_Bcast(shake_solute,1,MPI_LOGICAL,0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast shake solute')

call MPI_Bcast(Ndegf, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)    !bara i div init_
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast Ndegf')
call MPI_Bcast(Ndegfree, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)    !bara i div init_
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast Ndegfree')

call MPI_Bcast(shake_molecules,1,MPI_INTEGER,0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast shake molecules')
call MPI_Bcast(shake_constraints,1,MPI_INTEGER,0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast shake constraints')

!now we have to do some nasty stuff to get the shake data from master to the
!slaves, we can not just send the shake arrays, as they contain pointers and are
!thus not MPI transferable
!so we first make temp arrays for shake on master and send those to the nodes
if (shake_molecules .gt. 0) then
allocate(shakeconsttmp(shake_molecules), stat=alloc_status)
call check_alloc('Temporary shake arrays')
if (nodeid .eq. 0) then
shakenumconst = 0
do jj=1,shake_molecules
shakeconsttmp(jj)=shake_mol(jj)%nconstraints
shakenumconst = shakenumconst + shake_mol(jj)%nconstraints
end do
end if
if (nodeid .ne.0 ) then
allocate(shake_mol(nmol),stat = alloc_status)
call check_alloc('Slave node shake main array')
end if
call MPI_Bcast(shakenumconst,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast shakenumconst')
call MPI_Bcast(shakeconsttmp,shake_molecules,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast shakeconsttmp')

allocate(shakebondi(shakenumconst),shakebondj(shakenumconst),shakedist(shakenumconst),&
                stat=alloc_status)
call check_alloc('Temporary shake constraint arrays')
if (nodeid .eq. 0 ) then
shakenumconst = 0
do jj = 1,shake_molecules
        ii = 1
        do while (ii .le. shake_mol(jj)%nconstraints)
        shakenumconst = shakenumconst + 1
        shakebondi(shakenumconst) = shake_mol(jj)%bond(ii)%i
        shakebondj(shakenumconst) = shake_mol(jj)%bond(ii)%j
        shakedist(shakenumconst) = shake_mol(jj)%bond(ii)%dist2
        ii = ii +1
        end do
end do
end if
call MPI_Bcast(shakebondi,shakenumconst,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast shakebondi')
call MPI_Bcast(shakebondj,shakenumconst,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast shakebondj')
call MPI_Bcast(shakedist,shakenumconst,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast shakedist')

!yeah, now the slaves have all the shake information, need to put it now back
!into the data structure
if (nodeid .ne.0 ) then
do jj = 1, shake_molecules
shakenumconst = 0
shake_mol(jj)%nconstraints = shakeconsttmp(jj)
allocate(shake_mol(jj)%bond(shakeconsttmp(jj)),stat=alloc_status)
call check_alloc('slave node shkae bond array')
        ii = 1
        do while (ii .le. shake_mol(jj)%nconstraints)
        shakenumconst = shakenumconst + 1
        shake_mol(jj)%bond(ii)%i = shakebondi(shakenumconst)
        shake_mol(jj)%bond(ii)%j = shakebondj(shakenumconst)
        shake_mol(jj)%bond(ii)%dist2 = shakedist(shakenumconst)
        ii = ii + 1
        end do
end do
end if

deallocate(shakeconsttmp,shakebondi,shakebondj,shakedist)

end if ! shake_mol .gt. 0


!And the last -> length of bytes to recieve from the EQ arrays
call MPI_Bcast(reclength,1,MPI_INTEGER,0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast EQ array byte length')

if (nodeid .eq. 0) then 
call centered_heading('End of initiation', '-')
print *
end if

!and we sync before the end :)
call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!Finally allocate for  slaves:E_send, EQ_send
!For master :E_recv,d_recv
call allocate_mpi  

end subroutine init_nodes
#endif


!-----------------------------------------------------------------------

subroutine init_shake
!
! initialize shake constraints
!
!locals
integer						::	mol, b, ia, ja, constr, angle, i
real(kind=prec)						:: exclshk
integer						::	src, trg
integer						::	solute_shake_constraints
logical                                         ::      shaken

!allocate molecule list
allocate(shake_mol(nmol), stat=alloc_status)
call check_alloc('shake molecule array')

shake_mol(:)%nconstraints = 0
mol = 0
exclshk = zero

!count bonds to be constrained in each molecule
!also count shake constraints involving excluded atoms
do b=1,nbonds
        ia = bnd(b)%i
        ja = bnd(b)%j
        do while(ia >= istart_mol(mol+1))
                !new molecule
                mol = mol +1
        end do
        !skip redefined bonds
        if(bnd(b)%cod == 0) cycle
!new shake flag to spec if people want all atoms, all solvent hydrogens, solute
!hydrogens, solute atoms ...
        shaken = .false.
        if(shake_solute .and. (ia.le.nat_solute)) then
                if(.not. heavy(ia) .or. .not. heavy(ja)) then
                        if(shake_hydrogens) then
                                shake_mol(mol)%nconstraints = shake_mol(mol)%nconstraints + 1
                                shaken = .true.
                        end if
                else
                        if(shake_heavy) then
                                shake_mol(mol)%nconstraints = shake_mol(mol)%nconstraints + 1
                                shaken = .true.
                        end if
                end if
        end if
        if(shake_solvent .and. (ia.gt.nat_solute)) then
                if(.not. heavy(ia) .or. .not. heavy(ja)) then
                        if(shake_hydrogens) then
                                shake_mol(mol)%nconstraints = shake_mol(mol)%nconstraints + 1
                                shaken = .true.
                        end if
                else
                        if(shake_heavy) then
                                shake_mol(mol)%nconstraints = shake_mol(mol)%nconstraints + 1
                                shaken = .true.
                        end if
                end if
        end if

           if(shaken .and. ( .not. use_PBC) ) then
                if(excl(ia)) exclshk = exclshk + 0.5_prec
                if(excl(ja)) exclshk = exclshk + 0.5_prec
           end if

end do
!count extra shake constraints from fep file in appropriate molecule
do b = 1, nqshake
    ia=iqshake(b)
        mol = 1
        do while(ia >= istart_mol(mol+1))
                mol = mol + 1
        end do
        shake_mol(mol)%nconstraints = shake_mol(mol)%nconstraints + 1
end do

!allocate bond lists for each molecule
do mol = 1, nmol
        !allocate(sbtemp(nconstr(mol), status = alloc_status)
        allocate(shake_mol(mol)%bond(shake_mol(mol)%nconstraints), stat = alloc_status)
        call check_alloc('shake bond array')
        !shake_mol(mol)%bonds => sbtemp
end do

mol = 0
!add the constraint
do b=1,nbonds
        shaken = .false.
        ia = bnd(b)%i
        ja = bnd(b)%j
        do while(ia >= istart_mol(mol+1)) 
                !new molecule
                mol = mol +1
                shake_mol(mol)%nconstraints = 0
        end do
        !skip redefined bonds
        if(bnd(b)%cod == 0) cycle
        if(shake_solute .and. (ia.le.nat_solute)) then
                if(.not. heavy(ia) .or. .not. heavy(ja)) then
                        if(shake_hydrogens) then
                                shaken = .true.
                        end if
                else
                        if(shake_heavy) then
                                shaken = .true.
                        end if
                end if
        end if
        if(shake_solvent .and. (ia.gt.nat_solute)) then
                if(.not. heavy(ia) .or. .not. heavy(ja)) then
                        if(shake_hydrogens) then
                                shaken = .true.
                        end if
                else
                        if(shake_heavy) then
                                shaken = .true.
                        end if
                end if
        end if


        if(shaken) then
                shake_mol(mol)%nconstraints = shake_mol(mol)%nconstraints + 1
                shake_mol(mol)%bond(shake_mol(mol)%nconstraints)%i = ia
                shake_mol(mol)%bond(shake_mol(mol)%nconstraints)%j = ja
                shake_mol(mol)%bond(shake_mol(mol)%nconstraints)%dist2 = &
                        bondlib(bnd(b)%cod)%bnd0**2
                !set the bond code to -1 for shaken bonds
                !bnd(b) will be deleted by shrink_topology
                bnd(b)%cod = -1
        end if
end do

!add extra shake constraints from fep file to appropriate molecule
do b = 1, nqshake
    ia=iqshake(b)
        ja=jqshake(b)
        mol = 1
        do while(ia >= istart_mol(mol+1))
                mol = mol + 1
        end do
        !see if already shaken
        do constr = 1, shake_mol(mol)%nconstraints
                if((ia == shake_mol(mol)%bond(constr)%i .and. &
                        ja == shake_mol(mol)%bond(constr)%j) .or. &
                   (ja == shake_mol(mol)%bond(constr)%i .and. &
                        ia == shake_mol(mol)%bond(constr)%j)) then
                        !found it: will overwrite 
                        !also decrement number of constraints
                        shake_mol(mol)%nconstraints = shake_mol(mol)%nconstraints - 1
                        exit
                end if
        end do
        !constr now contains the right index
        shake_mol(mol)%bond(constr)%i = ia
        shake_mol(mol)%bond(constr)%j = ja
        shake_mol(mol)%bond(constr)%dist2 = &
                dot_product(EQ(1:nstates)%lambda,qshake_dist(b,1:nstates))**2
        shake_mol(mol)%nconstraints = shake_mol(mol)%nconstraints + 1
end do

!get total number of shake constraints in solute (used for separate scaling of temperatures)
solute_shake_constraints = sum(shake_mol(1:nmol-nwat)%nconstraints)


!remove molecules with zero constraints from list
trg = 1
src = 2
do while(src <= nmol)
        if(shake_mol(trg)%nconstraints == 0) then
                shake_mol(trg) = shake_mol(src)
                !clear source
                shake_mol(src)%nconstraints = 0
                src = src + 1
        else
                trg = trg + 1
                if(trg == src) src = src + 1
        end if
end do
shake_molecules = trg

!total number of constraints
shake_constraints = sum(shake_mol(1:shake_molecules)%nconstraints)
write(*,100) shake_constraints
write(*,101) shake_molecules
100	format(/,'Number of shake constraints             = ',i10)
101	format('No. molecules with shake constraints    = ',i10)
! calculate #degrees of freedom
Ndegf=3*natom-shake_constraints    !changed from Ndegf=3*natom-3-shake_constraints, center of mass position is NOT CONstrained in the simulation, but IS constrained for initial temperatures....
Ndegfree=Ndegf-3*nexats+exclshk

Ndegf_solvent = Ndegf - 3*nat_solute + solute_shake_constraints
Ndegf_solute = Ndegf - Ndegf_solvent

Ndegfree_solvent = 3*(natom - nat_solute) - (shake_constraints - solute_shake_constraints)
Ndegfree_solute = Ndegfree - Ndegfree_solvent

if (Ndegfree_solvent*Ndegfree_solute .eq. 0) then    ! if either solvent or solute have 0 degrees of freedom, turn off separate scaling (in case it's on) and do not print detailed temperatures
	detail_temps = .false.
	separate_scaling = .false.
else
	detail_temps = .true.
end if


!clear angles which are shaken (i and k atoms shaken)
do mol=1, shake_molecules
        do constr = 1, shake_mol(mol)%nconstraints
        ia = shake_mol(mol)%bond(constr)%i
        ja = shake_mol(mol)%bond(constr)%j
                do angle = 1, nangles
                        if((ang(angle)%i == ia .and. ang(angle)%k == ja) .or. &
                           (ang(angle)%i == ja .and. ang(angle)%k == ia)) then
                                ang(angle)%cod = 0
                                exit
                        end if
                end do
        end do
end do
! we can not initilaize the nose-hoover variables before
! as they need the shake information
! so we set it after shake has been started
if ( thermostat == NOSEHOOVER ) then
        do i=1,numchain
                qnh(i)=nhq
                vnh(i)=0
                xnh(i)=0
        end do
        qnh(1)=Ndegf*nhq
end if

end subroutine init_shake

!-----------------------------------------------------------------------

subroutine maxwell
! *** local variables
integer						:: i,j,k
real(kind=prec)					:: sd,kT
TYPE(qr_vec)                                    :: vg
!	Generate Maxwellian velocities 
kT = Boltz*Tmaxw

do i=1,natom
        sd = sqrt (kT/iaclib(iac(i))%mass)
        call q_gauss(zero,sd,vg,iseed)
        v(i) = vg
end do

end subroutine maxwell

!-----------------------------------------------------------------------
! Subroutine for the Nosé-Hoover chain propagation
subroutine nh_prop

real(kind=prec)				:: dt4, dt8, expf, s
integer				:: i, i3

!initializing all the constants to be used
dt4=0.5_prec*dt2
dt8=0.5_prec*dt4

Gnh(1) = Ndegf*Boltz*(Temp-Temp0)/qnh(1)

do i=2,numchain
	Gnh(i)=(qnh(i-1)*vnh(i-1)**2 - kbT)/qnh(i)
end do

vnh(numchain)=vnh(numchain)+Gnh(numchain)*dt4
xnh(numchain)=xnh(numchain)+vnh(numchain)*dt2

do i=numchain-1,1,-1
	expf=q_exp(-dt8*vnh(i+1))
	vnh(i)=(vnh(i)*expf+Gnh(i)*dt4)*expf
	xnh(i)=xnh(i)+vnh(i)*dt2
end do

s=q_exp(-vnh(1)*dt2)
do i=1,natom
	v(i)=v(i)*s
end do

Gnh(1) = Ndegf*Boltz*(s*s*Temp - Temp0)/qnh(1)
vnh(1) = (vnh(1)*expf+Gnh(1)*dt4)*expf

do i=2,numchain-1
        expf=q_exp(-dt8*vnh(i+1))
        Gnh(i)=(qnh(i-1)*vnh(i-1)**2 - kbT)/qnh(i)
        vnh(i)=(vnh(i)*expf+Gnh(i)*dt4)*expf
end do

Gnh(numchain)=(qnh(numchain-1)*vnh(numchain-1)**2 - kbT)/qnh(numchain)
vnh(numchain)=vnh(numchain)+Gnh(numchain)*dt4

end subroutine nh_prop

!-----------------------------------------------------------------------
subroutine prep_coord


! local variables
integer(4)          :: i,nat3
logical             :: old_restart = .false.
TYPE(qr_vec),allocatable :: x_old(:),v_old(:)
TYPE(qr_vec)        :: old_boxlength, old_boxcentre
integer             :: headercheck,myprec
!new variables for differen tprecision restart files
TYPE(qr_vecs),allocatable :: x_single(:),v_single(:)
TYPE(qr_vecs)             :: boxl_single,boxc_single
TYPE(qr_vecd),allocatable :: x_double(:),v_double(:)
TYPE(qr_vecd)             :: boxl_double,boxc_double
#ifndef PGI
TYPE(qr_vecq),allocatable   :: x_quad(:),v_quad(:)
TYPE(qr_vecq)               :: boxl_quad,boxc_quad
#endif
if (prec .eq. singleprecision) then
myprec = -137
elseif (prec .eq. doubleprecision) then
myprec = -1337
#ifndef PGI
elseif (prec .eq. quadprecision) then
myprec = -13337
#endif
else
call die('No such precision')
end if



! --- Refresh topology coords. if needed (external restraints file)
if ( implicit_rstr_from_file .eq. 1 ) then
write (*,'(/,a,/)') 'Refreshing topology coords for restraining...'
read(12) headercheck
  if ((headercheck .ne. -137).and.(headercheck.ne.-1337).and.(headercheck.ne.-13337)) then
!old restart file without header canary
      rewind(12)
      old_restart = .true.
  end if
      read (12) nat3
      rewind(12)
        if(nat3 /= 3*natom) then
          write(*,100) nat3/3, natom
          call die('wrong number of atoms in restart file')
       end if
  if (.not.old_restart) then
     read (12) headercheck
     if (myprec .ne. headercheck) then
       write(*,*) '>>> WARNING: Using mismatched precision in restart file'
       if (headercheck .eq. -137) then
         allocate(x_single(natom))
         read (12,err=112,end=112) nat3, (x_single(i),i=1,natom)
         xtop(1:natom)%x = x_single(1:natom)%x
         xtop(1:natom)%y = x_single(1:natom)%y
         xtop(1:natom)%z = x_single(1:natom)%z
        deallocate(x_single)
       else if (headercheck .eq. -1337) then
         allocate(x_double(natom))
         read (12,err=112,end=112) nat3, (x_double(i),i=1,natom)
         xtop(1:natom)%x = x_double(1:natom)%x
         xtop(1:natom)%y = x_double(1:natom)%y
         xtop(1:natom)%z = x_double(1:natom)%z
         deallocate(x_double)
       else if (headercheck .eq. -13337) then
#ifndef PGI
         allocate(x_quad(natom))
         read (12,err=112,end=112) nat3, (x_quad(i),i=1,natom)
         xtop(1:natom)%x = x_quad(1:natom)%x
         xtop(1:natom)%y = x_quad(1:natom)%y
         xtop(1:natom)%z = x_quad(1:natom)%z
         deallocate(x_quad)
#else
         call die('Quadruple precision not supported in PGI')
#endif
       end if
     else
       read (12,err=112,end=112) nat3, (xtop(i),i=1,natom)
     end if
  else
     allocate(x_old(natom))
     read (12,err=112,end=112) nat3, (x_old(i),i=1,natom)
     xtop(1:natom)%x = x_old(1:natom)%x
     xtop(1:natom)%y = x_old(1:natom)%y
     xtop(1:natom)%z = x_old(1:natom)%z
     deallocate(x_old)
  end if

end if
old_restart =.false.
headercheck=0
!Assign restraints of kind res:atom their numerical atom numbers
do i=1,nrstr_dist
  if(rstdis(i)%itext .ne. 'nil') then
    if (scan(rstdis(i)%itext,':') .ne. 0) then
      rstdis(i)%i=get_atom_from_resnum_atnum(rstdis(i)%itext)
	else
      read(rstdis(i)%itext,*) rstdis(i)%i
    end if
    if (scan(rstdis(i)%jtext,':') .ne. 0) then
      rstdis(i)%j=get_atom_from_resnum_atnum(rstdis(i)%jtext)
    else
      read(rstdis(i)%jtext,*) rstdis(i)%j
    end if
  end if
end do

! --- Make spherical restraining shell lists based on
!     the xtop coords.
if (.not. use_PBC) then

	if(rexcl_i > rexcl_o) then
	  call die('inner radius of restrained shell must be < exclusion radius')
	end if
	!first find atoms in shell 
	if (rexcl_i >= 0) then      !if rexcl_i is defined...
	  if (rexcl_i <= 1.00) then   !if rexcl_i is defined as fraction of rexcl_o
	    rexcl_i = rexcl_i * rexcl_o  !translate to Angstrom
	  end if 
	  if(iuse_switch_atom == 1) then
	     call make_shell
	  else
	     call make_shell2
	  end if
	else
	  write (*,'(/,a,/)') 'Restrained shell not defined!'
	end if
else
	shell(:) = .false.
end if ! .not. use_PBC
! --- read restart file

call allocate_natom_arrays
 if (thermostat == NOSEHOOVER) then
  call allocate_nhchain_arrays
 end if
if(restart) then
        ! topology routine has determined nwat, natom and allocated storage
        call centered_heading('Reading restart file','-')
        read(2) headercheck
        if ((headercheck .ne. -137).and.(headercheck.ne.-1337).and.(headercheck.ne.-13337)) then
!old restart file without header canary
          rewind(2)
          old_restart = .true.
        end if
          read (2) nat3
          rewind(2)
        if(nat3 /= 3*natom) then
                write(*,100) nat3/3, natom
100			format('>>>>> ERROR:',i5,' atoms in restart file not equal to',i5,&
                        ' in topology.')
                call die('wrong number of atoms in restart file')
        end if
        if (.not.old_restart) then
          read(2) headercheck
          if (myprec .ne. headercheck) then
            write(*,*) '>>> WARNING: Using mismatched precision in restart file'
            if (headercheck .eq. -137) then
              allocate(x_single(natom),v_single(natom))
              read (2,err=112,end=112) nat3, (x_single(i),i=1,natom)
              read (2,err=112,end=112) nat3, (v_single(i),i=1,natom)
              x(1:natom)%x = x_single(1:natom)%x
              x(1:natom)%y = x_single(1:natom)%y
              x(1:natom)%z = x_single(1:natom)%z
              v(1:natom)%x = v_single(1:natom)%x
              v(1:natom)%y = v_single(1:natom)%y
              v(1:natom)%z = v_single(1:natom)%z
              deallocate(x_single,v_single)
              if( use_PBC) then
                 read(2,err=112,end=112) boxl_single
                 read(2,err=112,end=112) boxc_single
                 boxlength%x = boxl_single%x
                 boxlength%y = boxl_single%y
                 boxlength%z = boxl_single%z
                 boxcentre%x = boxc_single%x
                 boxcentre%y = boxc_single%y
                 boxcentre%z = boxc_single%z
              end if
            else if (headercheck .eq. -1337) then
              allocate(x_double(natom),v_double(natom))
              read (2,err=112,end=112) nat3, (x_double(i),i=1,natom)
              read (2,err=112,end=112) nat3, (v_double(i),i=1,natom)
              x(1:natom)%x = x_double(1:natom)%x
              x(1:natom)%y = x_double(1:natom)%y
              x(1:natom)%z = x_double(1:natom)%z
              v(1:natom)%x = v_double(1:natom)%x
              v(1:natom)%y = v_double(1:natom)%y
              v(1:natom)%z = v_double(1:natom)%z
              deallocate(x_double,v_double)
              if( use_PBC) then
                 read(2,err=112,end=112) boxl_double
                 read(2,err=112,end=112) boxc_double
                 boxlength%x = boxl_double%x
                 boxlength%y = boxl_double%y
                 boxlength%z = boxl_double%z
                 boxcentre%x = boxc_double%x
                 boxcentre%y = boxc_double%y
                 boxcentre%z = boxc_double%z
              end if
            else if (headercheck .eq. -13337) then
#ifndef PGI
              allocate(x_quad(natom),v_quad(natom))
              read (2,err=112,end=112) nat3, (x_quad(i),i=1,natom)
              read (2,err=112,end=112) nat3, (v_quad(i),i=1,natom)
              x(1:natom)%x = x_quad(1:natom)%x
              x(1:natom)%y = x_quad(1:natom)%y
              x(1:natom)%z = x_quad(1:natom)%z
              v(1:natom)%x = v_quad(1:natom)%x
              v(1:natom)%y = v_quad(1:natom)%y
              v(1:natom)%z = v_quad(1:natom)%z
              deallocate(x_quad,v_quad)
              if( use_PBC) then
                 read(2,err=112,end=112) boxl_quad
                 read(2,err=112,end=112) boxc_quad
                 boxlength%x = boxl_quad%x
                 boxlength%y = boxl_quad%y
                 boxlength%z = boxl_quad%z
                 boxcentre%x = boxc_quad%x
                 boxcentre%y = boxc_quad%y
                 boxcentre%z = boxc_quad%z
              end if
#else
              call die('Quadruple precision not supported in PGI')
#endif
            end if
          else
              read (2,err=112,end=112) nat3, (x(i),i=1,natom)
              read (2,err=112,end=112) nat3, (v(i),i=1,natom)
              if( use_PBC) then
                read(2,err=112,end=112) boxlength
                read(2,err=112,end=112) boxcentre
              end if
          end if
        else
          allocate(x_old(natom),v_old(natom))
          read (2,err=112,end=112) nat3, (x_old(i),i=1,natom)
          read (2,err=112,end=112) nat3, (v_old(i),i=1,natom)
          write(*,*) 'Read coordinates and velocities from previous version of qdyn'
          x(1:natom)%x = x_old(1:natom)%x
          x(1:natom)%y = x_old(1:natom)%y
          x(1:natom)%z = x_old(1:natom)%z
          v(1:natom)%x = v_old(1:natom)%x
          v(1:natom)%y = v_old(1:natom)%y
          v(1:natom)%z = v_old(1:natom)%z
          deallocate(x_old,v_old)
          if( use_PBC) then
             read(2,err=112,end=112) old_boxlength
             read(2,err=112,end=112) old_boxcentre
             write(*,*) 'Read boxlength and center from previous version of qdyn'
             boxlength%x = old_boxlength%x
             boxlength%y = old_boxlength%y
             boxlength%z = old_boxlength%z
             boxcentre%x = old_boxcentre%x
             boxcentre%y = old_boxcentre%y
             boxcentre%z = old_boxcentre%z
          end if
        end if
        write (*,'(a30,i8)')   'Total number of atoms        =',natom
        write (*,'(a30,i8,/)') 'Number of waters encountered =',nwat

        if( use_PBC) then
                write(*,*)
                write(*,'(a16,3f8.3)') 'Boxlength     =', boxlength
                write(*,'(a16,3f8.3)') 'Centre of box =', boxcentre
        end if
        !water polarisation data will be read from restart file in wat_shells
else
        x(1:natom) = xtop(1:natom)
end if

! clear iqatom atom array
iqatom(:) = 0

return
#if defined(USE_MPI)
112 call MPI_Abort(MPI_COMM_WORLD,1,ierr)
#else
112 stop 'Aborting due to errors reading restart file.'
#endif

end subroutine prep_coord

!-----------------------------------------------------------------------

!Sort out heavy atoms in restrained shell. Use protein center to calculate distance.
!Uses coordinates from topology unless 'implicit_rstr_from_file' is specified.
subroutine make_shell
! *** Local variables
	integer						::	i,ig,i3
	real(kind=prec)						::	rin2,r2

	nshellats = 0
	rin2  = rexcl_i**2

	shell(:) = .false.

	do ig=1,ncgp_solute
       if (.not. excl(cgp(ig)%iswitch) .and. heavy(cgp(ig)%iswitch)) then 
                r2 = q_dist4(xtop(cgp(ig)%iswitch),xpcent)
		if(r2 > rin2) then
			do i=cgp(ig)%first, cgp(ig)%last
        nshellats = nshellats + 1 
				shell(cgpatom(i)) = .true.
			end do
		end if
   end if
	end do
	write(*,105) nshellats, rexcl_i, rexcl_o
105	format('Found   ',i6,' solute atoms in the restrained shell region (',f6.2,' to ',f6.2,')')
end subroutine make_shell

!------------------------------------------------------------------------

!Sort out heavy atoms in restrained shell. Use protein center to calculate distance.
!Use coordinates from topology unless 'implicit_rstr_from_file' is specified
subroutine make_shell2
! *** Local variables
	integer						::	i,ig,k
	real(kind=prec)						::	rout2,rin2,r2
        TYPE(qr_vec), allocatable		::	cgp_cent(:)
	nshellats = 0
	rin2  = rexcl_i**2

	shell(:) = .false.

	allocate(cgp_cent(ncgp+nwat))

	cgp_cent = cgp_cent * zero

	do ig=1,ncgp_solute
    if (.not. excl(cgp(ig)%iswitch) .and. heavy(cgp(ig)%iswitch)) then
                do i = cgp(ig)%first,cgp(ig)%last
			cgp_cent(ig) = cgp_cent(ig) + xtop(cgpatom(i))
		end do
        cgp_cent(ig) = cgp_cent(ig)/real(cgp(ig)%last - cgp(ig)%first +1, kind=prec)
		r2 = q_dist4(cgp_cent(ig),xpcent)

		if ( r2 .gt. rin2 ) then
			do i=cgp(ig)%first, cgp(ig)%last
				nshellats = nshellats + 1
				shell(cgpatom(i)) = .true.
			end do
		end if
   end if
	end do

	deallocate(cgp_cent) 
	write(*,105) nshellats, rexcl_i, rexcl_o
105	format('Found   ',i6,' solute atoms in the restrained shell region (',f6.2,' to ',f6.2,'Å)')
end subroutine make_shell2

!-----------------------------------------------------------------------

subroutine init_trj
!locals
integer						::	trj_atoms

!initialise trajectory atom mask
if(itrj_cycle > 0) then
        call trj_initialize(frames=nsteps/itrj_cycle, steps_before=itrj_cycle,&
                interval=itrj_cycle, steps=nsteps,	degf=ndegfree, &
                topfile=top_file)

        trj_atoms = trj_commit_mask()
        write(*,100) trj_atoms
        if(.not. trj_create(trj_file)) then
                call die('failure to open trajectory file')
        end if
end if

100	format('Coordinates for',i6,' atoms will be written to the trajectory.')
end subroutine init_trj

!-----------------------------------------------------------------------
subroutine precompute_interactions
! nice and small routine to call all the other stuff below
if (nat_solute .ne. 0) call pp_int_comp
if ((nwat .gt. 0) .and. (nat_solute .ne.0)) call pw_int_comp
if ((nat_solute .ne. 0) .and. (nqat .ne.0)) call qp_int_comp
if (nqat .ne. 0) call qq_int_comp
if ((nwat .gt. 0) .and. (nqat.ne.0)) call qw_int_comp
if (nwat .gt. 0) call ww_int_comp
if (use_excluded_groups) call exc_assign

end subroutine precompute_interactions

!-----------------------------------------------------------------------

subroutine exc_assign
! assigns the interactions between excluded groups and q atoms
! the information is gained from the previous generated lists
! but needs another pass over all atoms because we have to allocate the arrays first
!locals
integer                         :: ig,iq,jq,istate
integer                         :: igc
integer                         :: exc_qqlist,exc_qqplist,exc_qplist
! allocate all arrays
allocate(exc_nbqq_list(ngroups_gc,nstates,maxval(exc_nbqq)),exc_nbqqp_list(ngroups_gc,nstates,maxval(exc_nbqqp)),&
        exc_nbqp_list(ngroups_gc,nstates,maxval(exc_nbqp)))


exc_qqlist  = 1
exc_qplist  = 1
exc_qqplist = 1

do ig = 1, nat_solute
        do igc = 1, ngroups_gc
                if (ST_gc(igc)%gcmask%mask(ig)) then
                        do iq = 1, nqat
                        do istate = 1, nstates
                        if (iqatom(ig) .ne. 0) then
                        ! atom ig is a qatom, add to exc_qq list with state info
                                if (.not.qq_precomp(iq,iqatom(ig),istate)%set) cycle
                                ! if those two atoms are allowed to interact
                                exc_nbqq_list(igc,istate,exc_qqlist)%iq = iq
                                exc_nbqq_list(igc,istate,exc_qqlist)%jq = jq
                                exc_nbqq_list(igc,istate,exc_qqlist)%vdWA = qq_precomp(iq,jq,istate)%vdWA
                                exc_nbqq_list(igc,istate,exc_qqlist)%vdWB = qq_precomp(iq,jq,istate)%vdWB
                                exc_nbqq_list(igc,istate,exc_qqlist)%elec = qq_precomp(iq,jq,istate)%elec
                                exc_nbqq_list(igc,istate,exc_qqlist)%soft = qq_precomp(iq,jq,istate)%soft
                                exc_nbqq_list(igc,istate,exc_qqlist)%score = qq_precomp(iq,jq,istate)%score
                                exc_qqlist = exc_qqlist + 1
                        else if(any(qconn(:,ig,:) <= 3)) then
                        ! it is bonded somehow to a q atom        
                                if (.not.qp_precomp(ig,iq,istate)%set) cycle
                                ! ! if atoms can interact, should not be needed but who knows
                                exc_nbqqp_list(igc,istate,exc_qqplist)%i = iq
                                exc_nbqqp_list(igc,istate,exc_qqplist)%j = ig
                                exc_nbqqp_list(igc,istate,exc_qqplist)%vdWA = qp_precomp(ig,iq,istate)%vdWA
                                exc_nbqqp_list(igc,istate,exc_qqplist)%vdWB = qp_precomp(ig,iq,istate)%vdWB
                                exc_nbqqp_list(igc,istate,exc_qqplist)%elec = qp_precomp(ig,iq,istate)%elec
                                exc_nbqqp_list(igc,istate,exc_qqplist)%score = qp_precomp(ig,iq,istate)%score
                                exc_qqplist = exc_qqplist + 1
                        else
                        ! plain old nonbonded qp
                                if (.not.qp_precomp(iq,iq,istate)%set) cycle
                                exc_nbqp_list(igc,istate,exc_qplist)%i = iq
                                exc_nbqp_list(igc,istate,exc_qplist)%j = ig
                                exc_nbqp_list(igc,istate,exc_qplist)%vdWA = qp_precomp(ig,iq,istate)%vdWA
                                exc_nbqp_list(igc,istate,exc_qplist)%vdWB = qp_precomp(ig,iq,istate)%vdWB
                                exc_nbqp_list(igc,istate,exc_qplist)%elec = qp_precomp(ig,iq,istate)%elec
                                exc_nbqp_list(igc,istate,exc_qplist)%score = qp_precomp(ig,iq,istate)%score
                                exc_qplist = exc_qplist + 1
                        end if
                        end do ! nstates
                        end do ! nqat
                end if
        end do
end do

end subroutine exc_assign

subroutine pp_int_comp
! here comes the real stuff
! this one takes info from the list update routines to 
! make the full list of all possible interactions and have them ready for
! immidiate lookup

! locals
integer                         :: ig,jg,nl
integer                         :: ia,ja,i,j
allocate(pp_precomp(nat_solute,nat_solute),stat=alloc_status)
call check_alloc('Protein-Protein precomputation array')

! 
pp_precomp(:,:)%set   = .false.
pp_precomp(:,:)%score = zero
pp_precomp(:,:)%elec  = zero
pp_precomp(:,:)%vdWA  = zero
pp_precomp(:,:)%vdWB  = zero

igloop: do ig = calculation_assignment%pp%start, calculation_assignment%pp%end
        ia = cgp(ig)%iswitch
        if ( excl(ia) ) cycle igloop

jgloop: do jg = 1, ncgp_solute
                ja = cgp(jg)%iswitch
                if ( excl(ja) ) cycle jgloop
! count each charge group pair once only
                if ( ((ig .gt. jg) .and. (mod(ig+jg,2) .eq. 0)) .or. &
                        ((ig .lt. jg) .and. (mod(ig+jg,2) .eq. 1)) ) &
                        cycle jgloop
ialoop:         do ia = cgp(ig)%first, cgp(ig)%last
                        i = cgpatom(ia)
!             --- q-atom ? ---
                        if ( iqatom(i) .ne. 0 ) cycle ialoop

jaloop:                 do ja = cgp(jg)%first, cgp(jg)%last
                                j = cgpatom(ja)
!             --- q-atom ? ---
                                if ( iqatom(j).ne.0 ) cycle jaloop
! count once
                                if ( ig .eq. jg .and. i .ge. j ) cycle jaloop

                                if ( abs(j-i) .le. max_nbr_range ) then
                                        if ( i .lt. j ) then
                                                if ( listex(j-i,i) ) then
                                                        cycle jaloop
                                                else if (list14(j-i,i)) then
                                                        call precompute_set_values_pp(i,j,3)
                                                        cycle jaloop
                                                end if
                                        else
                                                if ( listex(i-j,j) ) then
                                                        cycle jaloop
                                                else if ( list14(i-j,j)) then
                                                        call precompute_set_values_pp(i,j,3)
                                                        cycle jaloop
                                                end if
                                        end if
                                else
                                        do nl = 1, nexlong
                                                if ( (listexlong(1,nl) .eq. i .and. &
                                                        listexlong(2,nl) .eq. j      ) .or. &
                                                        (listexlong(1,nl) .eq. j .and. &
                                                        listexlong(2,nl) .eq. i      ) ) then
                                                        cycle jaloop
                                                end if
                                        end do
                                        do nl = 1, n14long
                                                if ( (list14long(1,nl) .eq. i .and. &
                                                        list14long(2,nl) .eq. j      ) .or. &
                                                        (list14long(1,nl) .eq. j .and. &
                                                        list14long(2,nl) .eq. i      ) ) then
                                                        call precompute_set_values_pp(i,j,3)
                                                        cycle jaloop
                                                end if
                                        end do
                                end if
                                call precompute_set_values_pp(i,j,ljcod(iac(i),iac(j)))
                        end do jaloop
                end do ialoop
        end do jgloop
end do igloop
end subroutine pp_int_comp

!-----------------------------------------------------------------------

subroutine pw_int_comp
! locals
integer                         :: ig,jg,ia,a_ind,b_ind

allocate(pw_precomp(nat_solute,solv_atom),stat=alloc_status)
call check_alloc('Protein-Solvent precomputation array')

pw_precomp(:,:)%set = .false.
pw_precomp(:,:)%score = zero
pw_precomp(:,:)%elec  = zero
pw_precomp(:,:)%vdWA  = zero
pw_precomp(:,:)%vdWB  = zero

igloop:do ig = calculation_assignment%pw%start, calculation_assignment%pw%end
jgloop: do jg = 1, solv_atom
ialoop:         do ia = cgp(ig)%first, cgp(ig)%last
                        a_ind = cgpatom(ia)
                        b_ind = jg
                        if ( iqatom(a_ind) .ne. 0 ) cycle ialoop
                        call precompute_set_values_pw(a_ind,b_ind,ljcod(iac(a_ind),iac(nat_solute+b_ind)))
                end do ialoop
        end do jgloop
end do igloop
end subroutine pw_int_comp

!-----------------------------------------------------------------------

subroutine qp_int_comp
! locals
integer                         :: ig,jq,ia,a_ind,is,vdw,j,k

allocate(qp_precomp(nat_solute,nqat,nstates),stat=alloc_status)
call check_alloc('Protein-QAtom precomputation array')
if(use_excluded_groups) then
        allocate(exc_nbqp(ngroups_gc),exc_nbqqp(ngroups_gc))
end if
qp_precomp(:,:,:)%set   = .false.
qp_precomp(:,:,:)%score = zero
qp_precomp(:,:,:)%elec  = zero
qp_precomp(:,:,:)%vdWA  = zero
qp_precomp(:,:,:)%vdWB  = zero

igloop:do ig = 1, nat_solute
        if(any(qconn(:,ig,:) .le. 3)) cycle igloop
        do jq = 1, nqat
                vdw = ljcod(iac(ig),iac(iqseq(jq)))
                do is = 1 ,nstates
                        if(qconn(is, ig, jq) .eq. 4) vdw = 3
                        call precompute_set_values_qp(jq,ig,is,vdw)
                        if (use_excluded_groups) then
                                do k = 1, ngroups_gc
                                if(ST_gc(k)%gcmask%mask(ig)) then
                                        exc_nbqp(k) = exc_nbqp(k) + 1
                                end if
                                end do
                        end if
                end do
        end do
end do igloop
end subroutine qp_int_comp

!-----------------------------------------------------------------------

subroutine qq_int_comp
! locals
integer                         :: iq,jq,ia,ja,k,l,vdw,i,is
real(kind=prec)                 :: tmp_elscale
logical                         :: found
allocate(qq_precomp(nqat,nqat,nstates),stat=alloc_status)
call check_alloc('QAtom-QAtom precomputation array')
if (use_excluded_groups) then
        allocate(exc_nbqq(ngroups_gc))
end if

qq_precomp(:,:,:)%set   = .false.
qq_precomp(:,:,:)%soft  = .false.
qq_precomp(:,:,:)%score = zero
qq_precomp(:,:,:)%elec  = zero
qq_precomp(:,:,:)%vdWA  = zero
qq_precomp(:,:,:)%vdWB  = zero

do iq = 1, nqat - 1
        ia = iqseq(iq)
        do jq = iq + 1, nqat
                ja = iqseq(jq)
                do is = 1, nstates
                        if(qconn(is, ja, iq) .ge. 4) then
                                tmp_elscale = one
                                if (nel_scale .ne. 0) then
                                        i = 1
                                        found = .false.
                                        do while(( i .le. nel_scale) .and. &
                                                (.not.found))
                                                k=qq_el_scale(i)%iqat
                                                l=qq_el_scale(i)%jqat
                                                if ((iq .eq. k .and. jq .eq. l) .or. &
                                                        (iq .eq. l .and. jq .eq. k)) then
                                                        tmp_elscale = qq_el_scale(i)%el_scale(is)
                                                        found = .true.
                                                        ! from Masoud
                                                end if
                                                i = i + 1
                                        end do
                                end if
                                if(qconn(is, ja, iq) .eq. 4) then
                                        vdw = 3
                                elseif(.not. qvdw_flag) then
                                        vdw = ljcod(iac(ia),iac(ja))
                                else
                                        vdw = 1
                                        i = 1
                                        found = .false.
                                        do while(( i .le. nqexpnb) .and. &
                                                (.not.found))
                                                if ((iq .eq. iqexpnb(i) .and.  &
                                                        jq .eq. jqexpnb(i))  .or. &
                                                        ( jq .eq. iqexpnb(i) .and. &
                                                        iq .eq. jqexpnb(i)))  then
                                                        vdw  = 2
                                                        found = .true.
                                                end if
                                                i = i + 1
                                        end do
                                end if ! (qconn = 4)
                                call precompute_set_values_qq(iq,jq,is,vdw,tmp_elscale)
                                if(use_excluded_groups) then
                                        do k = 1, ngroups_gc
                                        if(ST_gc(k)%gcmask%mask(ia)) then
                                                exc_nbqq(k) = exc_nbqq(k) + 1
                                                !call exc_add_int_qq(k,iq,jq)
                                        else if(ST_gc(k)%gcmask%mask(ja)) then
                                                exc_nbqq(k) = exc_nbqq(k) + 1
                                        end if
                                        end do
                                end if
                        end if
                end do
        end do
end do
! now we put those that were before also on q-q, but are actually q-p, on the
! q-p list where they belong!
do ja = 1, nat_solute
        if(iqatom(ja) .ne. 0) cycle
        if(any(qconn(:,ja,:) .le. 3)) then
                !bonded or angled to at least one Q-atom
                do iq = 1, nqat
                        do is = 1, nstates
                                if(qconn(is, ja, iq) .ge. 4) then
                                        if(qconn(is, ja, iq) .eq. 4) then
                                                vdw = 3
                                        elseif(qvdw_flag) then
                                                vdw = 1
                                        else
                                                vdw = ljcod(iac(ia),iac(ja))
                                        end if
                                        call precompute_set_values_qp(iq,ja,is,vdw)
                                        if (use_excluded_groups) then
                                                do k = 1, ngroups_gc
                                                if (ST_gc(k)%gcmask%mask(ja)) then
                                                        exc_nbqqp(k) = exc_nbqqp(k) + 1
                                                end if
                                                end do
                                        end if
                                end if
                        end do
                end do
        end if
end do
! prepare q-atom nonbond lists that do not need updating
#ifdef USE_MPI
if(nodeid .eq.0) then
#endif
call nbqqlist
#ifdef USE_MPI
end if
#endif
end subroutine qq_int_comp

!-----------------------------------------------------------------------

subroutine qw_int_comp
! locals
integer                         :: iq,jg,ia,a_ind,b_ind

allocate(qw_precomp(nqat,solv_atom,nstates),stat=alloc_status)
call check_alloc('QAtom-Solvent precomputation array')

qw_precomp(:,:,:)%set = .false.
qw_precomp(:,:,:)%score = zero
qw_precomp(:,:,:)%elec  = zero
qw_precomp(:,:,:)%vdWA  = zero
qw_precomp(:,:,:)%vdWB  = zero

igloop:do iq = 1, nqat
jgloop: do jg = 1, solv_atom
                a_ind = iq
                b_ind = jg
                call precompute_set_values_qw(a_ind,b_ind,ljcod(iac(nat_solute+b_ind),iac(iqseq(iq))))
        end do jgloop
end do igloop
end subroutine qw_int_comp

!-----------------------------------------------------------------------

subroutine ww_int_comp
! locals
integer                         :: ig,jg,ia,a_ind,b_ind

allocate(ww_precomp(solv_atom,solv_atom),stat=alloc_status)
call check_alloc('Solvent-Solvent precomputation array')

ww_precomp(:,:)%set = .false.
ww_precomp(:,:)%score = zero
ww_precomp(:,:)%elec  = zero
ww_precomp(:,:)%vdWA  = zero
ww_precomp(:,:)%vdWB  = zero

igloop:do ig = 1, solv_atom
jgloop: do jg = 1, solv_atom
                a_ind = ig
                b_ind = jg
                call precompute_set_values_ww(a_ind,b_ind,ljcod(iac(nat_solute+a_ind),iac(nat_solute+b_ind)))
        end do jgloop
end do igloop

end subroutine ww_int_comp

!-----------------------------------------------------------------------

subroutine precompute_set_values_pp(i,j,vdw)
! arguments
integer                         :: i,j,vdw
! locals
real(kind=prec)                 :: tempA,tempB
select case(ivdw_rule)
case(VDW_GEOMETRIC)
tempA = iaclib(iac(i))%avdw(vdw) * iaclib(iac(j))%avdw(vdw)
tempB = iaclib(iac(i))%bvdw(vdw) * iaclib(iac(j))%bvdw(vdw)
pp_precomp(i,j)%vdWA = tempA
pp_precomp(i,j)%vdWB = tempB

case(VDW_ARITHMETIC)
tempA = iaclib(iac(i))%avdw(vdw) + iaclib(iac(j))%avdw(vdw)
tempA = tempA**2
tempA = tempA * tempA * tempA
tempB = iaclib(iac(i))%bvdw(vdw) * iaclib(iac(j))%bvdw(vdw)

pp_precomp(i,j)%vdWA = (tempA**2) * tempB
pp_precomp(i,j)%vdWB = 2.0_prec * tempA * tempB

end select
pp_precomp(i,j)%elec = crg(i) * crg(j)

if ( (pp_precomp(i,j)%vdWA.ne.zero) .or. (pp_precomp(i,j)%vdWB.ne.zero) &
        .or. (pp_precomp(i,j)%elec.ne.zero)) pp_precomp(i,j)%set  = .true.

if (vdw .eq. 3 ) pp_precomp(i,j)%elec = pp_precomp(i,j)%elec * el14_scale

! make sure both possible ways to get there are prepared
pp_precomp(j,i) = pp_precomp(i,j)

end subroutine precompute_set_values_pp

!-----------------------------------------------------------------------

subroutine precompute_set_values_pw(i,j,vdw)
! arguments
integer                         :: i,j,vdw
! locals
real(kind=prec)                 :: tempA,tempB
select case(ivdw_rule)
case(VDW_GEOMETRIC)
tempA = iaclib(iac(i))%avdw(vdw) * aLJ_solv(j,vdw)
tempB = iaclib(iac(i))%bvdw(vdw) * bLJ_solv(j,vdw)
pw_precomp(i,j)%vdWA = tempA
pw_precomp(i,j)%vdWB = tempB

case(VDW_ARITHMETIC)
tempA = iaclib(iac(i))%avdw(vdw) + aLJ_solv(j,vdw)
tempA = tempA**2
tempA = tempA * tempA * tempA
tempB = iaclib(iac(i))%bvdw(vdw) * bLJ_solv(j,vdw)

pw_precomp(i,j)%vdWA = (tempA**2) * tempB
pw_precomp(i,j)%vdWB = 2.0_prec * tempA * tempB

end select
pw_precomp(i,j)%elec = crg(i) * chg_solv(j)

if ( (pw_precomp(i,j)%vdWA.ne.zero) .or. (pw_precomp(i,j)%vdWB.ne.zero) &
        .or. (pw_precomp(i,j)%elec.ne.zero)) pw_precomp(i,j)%set  = .true.

end subroutine precompute_set_values_pw

!-----------------------------------------------------------------------

subroutine precompute_set_values_qp(iq,j,istate,vdw)
! arguments
integer                         :: iq,i,j,vdw,qvdw,istate
! locals
real(kind=prec)                 :: tempA,tempB

i = iqseq(iq)
if (vdw .eq. 2 ) then
qvdw = 1
else
qvdw = vdw
end if
! for reference, we put all (!!!!!!!!!einself) qp interactions here
! even those listed for some reason under nbqqlist
! now nbqq only (!!!einself) has the q-q interactions

select case(ivdw_rule)
case(VDW_GEOMETRIC)
if (qvdw_flag) then
tempA = qavdw(qiac(iq,istate),qvdw) * iaclib(iac(j))%avdw(vdw)
tempB = qbvdw(qiac(iq,istate),qvdw) * iaclib(iac(j))%bvdw(vdw)
else
tempA = iaclib(iac(i))%avdw(vdw) * iaclib(iac(j))%avdw(vdw)
tempB = iaclib(iac(i))%bvdw(vdw) * iaclib(iac(j))%bvdw(vdw)
end if
qp_precomp(j,iq,istate)%vdWA = tempA
qp_precomp(j,iq,istate)%vdWB = tempB
case(VDW_ARITHMETIC)
if (qvdw_flag) then
tempA = qavdw(qiac(iq,istate),qvdw) + iaclib(iac(j))%avdw(vdw)
tempB = qbvdw(qiac(iq,istate),qvdw) * iaclib(iac(j))%bvdw(vdw)
else
tempA = iaclib(iac(i))%avdw(vdw) + iaclib(iac(j))%avdw(vdw)
tempB = iaclib(iac(i))%bvdw(vdw) * iaclib(iac(j))%bvdw(vdw)
end if
tempA = tempA**2
tempA = tempA * tempA * tempA
qp_precomp(j,iq,istate)%vdWA = (tempA**2) * tempB
qp_precomp(j,iq,istate)%vdWB = 2.0_prec * tempA * tempB
end select
if(.not. qq_use_library_charges) then
qp_precomp(j,iq,istate)%elec = qcrg(iq,istate) * crg(j)
else
qp_precomp(j,iq,istate)%elec = crg(i) * crg(j)
end if

if ( (qp_precomp(j,iq,istate)%vdWA.ne.zero) .or. (qp_precomp(j,iq,istate)%vdWB.ne.zero) &
        .or. (qp_precomp(j,iq,istate)%elec.ne.zero)) qp_precomp(j,iq,istate)%set  = .true.

qp_precomp(j,iq,istate)%score = sc_lookup(iq,iac(j),istate)
if (vdw .eq. 3 ) qp_precomp(j,iq,istate)%elec = qp_precomp(j,iq,istate)%elec * el14_scale

end subroutine precompute_set_values_qp

!-----------------------------------------------------------------------

subroutine precompute_set_values_qq(iq,jq,istate,vdw,q_elscale)
! arguments
integer                         :: iq,jq,i,j,vdw,istate
real(kind=prec)                 :: q_elscale ! from Masoud
! locals
real(kind=prec)                 :: tempA,tempB

i = iqseq(iq)
j = iqseq(jq)

! this one needs special treatment, as the nbqqlist routine makes a list of both
! q-q and q-p interactions for some reason (because of bonded/angled
! interactions to at least one q atom)
! so we just split it because else this will lead to branch points later on

select case(ivdw_rule)
case(VDW_GEOMETRIC)
if (qvdw_flag) then
tempA = qavdw(qiac(iq,istate),vdw) * qavdw(qiac(jq,istate),vdw)
tempB = qbvdw(qiac(iq,istate),vdw) * qbvdw(qiac(jq,istate),vdw)
else
tempA = iaclib(iac(i))%avdw(vdw) * iaclib(iac(j))%avdw(vdw)
tempB = iaclib(iac(i))%bvdw(vdw) * iaclib(iac(j))%bvdw(vdw)
end if
qq_precomp(iq,jq,istate)%vdWA = tempA
qq_precomp(iq,jq,istate)%vdWB = tempB
case(VDW_ARITHMETIC)
if (qvdw_flag) then
! soft pair case needs special treatment !
if (vdw .eq. 2 ) then
tempA = qavdw(qiac(iq,istate),vdw) * qavdw(qiac(jq,istate),vdw)
else
tempA = qavdw(qiac(iq,istate),vdw) + qavdw(qiac(jq,istate),vdw)
end if
tempB = qbvdw(qiac(iq,istate),vdw) * qbvdw(qiac(jq,istate),vdw)
else
tempA = iaclib(iac(i))%avdw(vdw) + iaclib(iac(j))%avdw(vdw)
tempB = iaclib(iac(i))%bvdw(vdw) * iaclib(iac(j))%bvdw(vdw)
end if
! soft pair case needs special treatment !
if ((vdw .eq. 2 ) .and. (qvdw_flag)) then
qq_precomp(iq,jq,istate)%vdWA = tempA
qq_precomp(iq,jq,istate)%vdWB = tempB
else
tempA = tempA**2
tempA = tempA * tempA * tempA
qq_precomp(iq,jq,istate)%vdWA = (tempA**2) * tempB
qq_precomp(iq,jq,istate)%vdWB = 2.0_prec * tempA * tempB
end if
end select
if(.not. qq_use_library_charges) then
qq_precomp(iq,jq,istate)%elec = qcrg(iq,istate) * qcrg(jq,istate)
else
qq_precomp(iq,jq,istate)%elec = crg(i) * crg(j)
end if
qq_precomp(iq,jq,istate)%elec  = qq_precomp(iq,jq,istate)%elec * q_elscale
qq_precomp(iq,jq,istate)%score = sc_lookup(iq,natyps+jq,istate)

if ( (qq_precomp(iq,jq,istate)%vdWA.ne.zero) .or. (qq_precomp(iq,jq,istate)%vdWB.ne.zero) &
        .or. (qq_precomp(iq,jq,istate)%elec.ne.zero)) qq_precomp(iq,jq,istate)%set = .true.

if (vdw .eq. 3 ) qq_precomp(iq,jq,istate)%elec = qq_precomp(iq,jq,istate)%elec * el14_scale
if ((vdw .eq. 2 ).and.(qvdw_flag)) qq_precomp(iq,jq,istate)%soft = .true.
end subroutine precompute_set_values_qq

!-----------------------------------------------------------------------

subroutine precompute_set_values_qw(iq,j,vdw)
! arguments
integer                         :: iq,i,j,vdw,istate,iacj,qvdw
! locals
real(kind=prec)                 :: tempA,tempB
! need to check first if vdw is right, qvdw is 1,1,3 instead of normal 1,2,3
if (vdw .eq. 2 ) then
qvdw = 1 
else
qvdw = vdw
end if

i = iqseq(iq)
iacj = iac(nat_solute+j)
do istate = 1, nstates
select case(ivdw_rule)
case(VDW_GEOMETRIC)
if (qvdw_flag) then
tempA = qavdw(qiac(iq,istate),qvdw) * aLJ_solv(j,vdw)
tempB = qbvdw(qiac(iq,istate),qvdw) * bLJ_solv(j,vdw)
else
tempA = iaclib(iac(i))%avdw(vdw) * aLJ_solv(j,vdw)
tempB = iaclib(iac(i))%bvdw(vdw) * bLJ_solv(j,vdw)
end if
qw_precomp(iq,j,istate)%vdWA = tempA
qw_precomp(iq,j,istate)%vdWB = tempB
case(VDW_ARITHMETIC)
if (qvdw_flag) then
tempA = qavdw(qiac(iq,istate),qvdw) + aLJ_solv(j,vdw)
tempB = qbvdw(qiac(iq,istate),qvdw) * bLJ_solv(j,vdw)
else
tempA = iaclib(iac(i))%avdw(vdw) + aLJ_solv(j,vdw)
tempB = iaclib(iac(i))%bvdw(vdw) * bLJ_solv(j,vdw)
end if
tempA = tempA**2
tempA = tempA * tempA * tempA
qw_precomp(iq,j,istate)%vdWA = (tempA**2) * tempB
qw_precomp(iq,j,istate)%vdWB = 2.0_prec * tempA * tempB
end select
if(.not. qq_use_library_charges) then
qw_precomp(iq,j,istate)%elec = qcrg(iq,istate) * chg_solv(j)
else
qw_precomp(iq,j,istate)%elec = crg(i) * chg_solv(j)
end if

qw_precomp(iq,j,istate)%score = sc_lookup(iq,iacj,istate)
end do

end subroutine precompute_set_values_qw

!-----------------------------------------------------------------------

subroutine precompute_set_values_ww(i,j,vdw)
! arguments
integer                         :: i,j,vdw
! locals
real(kind=prec)                 :: tempA,tempB
select case(ivdw_rule)
case(VDW_GEOMETRIC)
tempA = aLJ_solv(i,vdw) * aLJ_solv(j,vdw)
tempB = bLJ_solv(i,vdw) * bLJ_solv(j,vdw)
ww_precomp(i,j)%vdWA = tempA
ww_precomp(i,j)%vdWB = tempB

case(VDW_ARITHMETIC)
tempA = aLJ_solv(i,vdw) + aLJ_solv(j,vdw)
tempA = tempA**2
tempA = tempA * tempA * tempA
tempB = bLJ_solv(i,vdw) * bLJ_solv(j,vdw)

ww_precomp(i,j)%vdWA = (tempA**2) * tempB
ww_precomp(i,j)%vdWB = 2.0_prec * tempA * tempB

end select

ww_precomp(i,j)%set  = .true.
ww_precomp(i,j)%elec = chg_solv(i) * chg_solv(j)

end subroutine precompute_set_values_ww

!-----------------------------------------------------------------------
subroutine prep_sim
! local variables
integer						:: i, j, ig, istate,runvar,iw,irc_solvent
type(Q_ENE_HEAD)				:: tempheader

if (nodeid .eq. 0) then	
        write(*,*)
        call centered_heading('Initialising dynamics', '-')
end if

! Set parameters (bonds, angles, charges,...) & restraints for water   
! also get the number of atoms in each solvent molecule
! just take all non solute atoms and divide by nwat
! allocate chg_solv array here after getting number of solvent atoms/molecule
if(nwat > 0) then
! solv atom is read in during topo read now :)
	allocate(chg_solv(solv_atom))
	allocate(aLJ_solv(solv_atom,3))
	allocate(bLJ_solv(solv_atom,3))
        select case (solvent_type)
	 case (SOLVENT_SPC,SOLVENT_ALLATOM)
		do iw=1,solv_atom
			chg_solv(iw)=crg(nat_solute+iw)
			aLJ_solv(iw,1:3)=iaclib(iac(nat_solute+iw))%avdw(1:3)
			bLJ_solv(iw,1:3)=iaclib(iac(nat_solute+iw))%bvdw(1:3)
		end do
        case(SOLVENT_GENERAL)
                !add appropriate code here
                call die('Topology contains mixed solvent. This feature is not implemented yet.')
        end select

        if( .not. use_PBC ) then
! sovent density also set in wat_sphere to value from topology
! default will be  0.0335, water density
                call wat_sphere
                if (wpol_restr) call wat_shells

        else !compute charges of the system for box case 
        !(done in subroutine wat_sphere for sphere case)
                !calc. total charge of non-Q-atoms
                crgtot = zero
                do i = 1, nat_solute
                        if ( iqatom(i)==0 ) crgtot = crgtot + crg(i)
                end do
                write (*,60) crgtot
60 format ('Total charge of non-Q atoms             = ',f10.2)


        !calc effective charge of whole system at this lambda
                crgQtot = zero
                do i = 1, nqat
                        do istate = 1, nstates
                        crgtot = crgtot + qcrg(i,istate)*EQ(istate)%lambda
                        crgQtot = crgQtot + qcrg(i,istate)*EQ(istate)%lambda
                        end do
                end do

                write (*,70) crgtot

70 format ('Total charge of system                  = ',f10.2)

        end if
end if
! set the charge group membership for every topology atom only if using LRF or PBC
if(use_LRF .or. use_PBC) then
        call allocate_lrf_arrays

        do ig = 1, ncgp
                do i = cgp(ig)%first, cgp(ig)%last
                        iwhich_cgp(cgpatom(i)) = ig
                end do
        end do
!set nsolvent and nsolute to values from topology
!for box control later
	nsolvent = nres - nres_solute
	nsolute = nmol - nsolvent
end if

!	Prepare an array of inverse masses
winv(:) = one/iaclib(iac(:))%mass

if(use_PBC) allocate(mvd_mol(nmol))

if( use_PBC ) then 
        !compute masses of all molecules
        allocate(mol_mass(nmol))
        mol_mass(:) = zero

        do i = 1,nmol-1 !all molecules but the last
                do j = istart_mol(i), istart_mol(i+1)-1 !all atoms of molecule
                        mol_mass(i) = mol_mass(i) + iaclib(iac(j))%mass
                end do
        end do

        do j = istart_mol(nmol), natom !last molecule
                mol_mass(nmol) = mol_mass(nmol) + iaclib(iac(j))%mass
        end do

        mol_mass(:) = one/mol_mass(:)

        !prepare array of masses
        allocate(mass(natom))
!make this a duplicate of the iaclib array to avoid truncation errors
!        mass(:) = one/winv(:)
        mass(:)=iaclib(iac(:))%mass
! moved here or it will segfault :(
        if (control_box) then
        boxlength = new_boxl
        if ( put_solute_back_in_box .or. put_solvent_back_in_box ) then !only put back in box if either solute or solvent should be put back (qdyn input option)
                call put_back_in_box
        end if
        write(*,'(a)') 'Boxsize changed. Equilibration may be needed'
        end if
end if

!initialization for the Nosé-Hoover chain thermostat
!is done after shake has started, here shake_constraints are still zero

!scale charges by sqrt(coulomb_constant) 
crg(:) = crg(:) * sqrt(coulomb_constant)
if(nwat > 0) then
chg_solv(:) = chg_solv(:) * sqrt(coulomb_constant)
end if
if(nqat > 0) then
        qcrg(:,:) = qcrg(:,:) * sqrt(coulomb_constant)
end if
if (use_excluded_groups) then
!start writing the energy file header with the needed information
!one array for each group + 1 for the default environment
! needs to be higher if QCP also on
	write(*,*)
	write(*,*) 'Preparing residue groups for group contribution calculation'
	ene_header%arrays = ngroups_gc + 1 + QCP_N
	allocate(ene_header%types(ngroups_gc+1),ene_header%numres(ngroups_gc+1),&
			ene_header%gcnum(ngroups_gc+1))
	runvar = 2
        do i=1,ngroups_gc
                call mask_initialize(ST_gc(i)%gcmask)
                call gc_make_array(ST_gc(i))
		ene_header%numres(runvar)=ST_gc(i)%count
		select case (trim(ST_gc(i)%caltypen))
			case ('full')
			ene_header%types(runvar)=FULL
			ene_header%gcnum(runvar)=i
			case ('electro')
			ene_header%types(runvar)=ELECTRO
			ene_header%gcnum(runvar)=i
			case ('vdw')
			ene_header%types(runvar)=VDW
			ene_header%gcnum(runvar)=i
			case ('all')
!special case, do all the calculations together
!needs reallocation of header arrays
			allocate(tempheader%types(runvar),tempheader%numres(runvar),tempheader%gcnum(runvar))
			tempheader%types(1:runvar)=ene_header%types(1:runvar)
			tempheader%numres(1:runvar)=ene_header%numres(1:runvar)
			tempheader%gcnum(1:runvar)=ene_header%gcnum(1:runvar)
			deallocate(ene_header%types,ene_header%numres,ene_header%gcnum)
!add two more calculations
!			runvar = ene_header%arrays + 2
			allocate(ene_header%types(ene_header%arrays+2),&
			ene_header%numres(ene_header%arrays+2),&
			ene_header%gcnum(ene_header%arrays+2))
			ene_header%types(1:runvar)=tempheader%types(1:runvar)
			ene_header%types(runvar)=FULL
			ene_header%types(runvar+1)=ELECTRO
			ene_header%types(runvar+2)=VDW
			ene_header%numres(1:runvar)=tempheader%numres(1:runvar)
			ene_header%numres(runvar+1)=ene_header%numres(runvar)
			ene_header%numres(runvar+2)=ene_header%numres(runvar)
			ene_header%gcnum(1:runvar)=tempheader%gcnum(1:runvar)
			ene_header%gcnum(runvar)=i
			ene_header%gcnum(runvar+1)=i
			ene_header%gcnum(runvar+2)=i
			ene_header%arrays=ene_header%arrays+2
			runvar = runvar + 2
			deallocate(tempheader%types,tempheader%numres,tempheader%gcnum)


			case default
!should never ever reach this, end execution now!
			write(*,*) 'Unrecognized statement in group contribution case list'
			call die('group contribution case list')
		end select
                do j=1,ST_gc(i)%count
                	stat = mask_add(ST_gc(i)%gcmask,ST_gc(i)%sendtomask(j))
                end do
		runvar = runvar + 1
	write (*,*) 'Finished group ',i
        end do
	do i=2,ene_header%arrays
		ene_header%totresid = ene_header%totresid + ene_header%numres(i)
	end do
else
!no extra groups, write minimal info to header
!only need to know if QCP is being performed
	ene_header%arrays = 1 + QCP_N
	ene_header%totresid = 0 ! to prevent allocation errors
	allocate(ene_header%types(1+QCP_N),ene_header%numres(1+QCP_N),ene_header%gcnum(1+QCP_N))
end if
	ene_header%totresid = ene_header%totresid + 1
	allocate(ene_header%resid(ene_header%totresid))
	ene_header%types(1)=NOGC
	ene_header%numres(1)=1
	ene_header%resid(1)=-1
	ene_header%gcnum(1)=-1
	runvar=2
	do i=2,ene_header%arrays - QCP_N
!give each residue its topology number
		do j=1,ene_header%numres(i)
			ene_header%resid(runvar)=get_from_mask(ST_gc(ene_header%gcnum(i)),j)
			runvar = runvar + 1
		end do
	end do
! set qcp_pos to last array entry
if (use_qcp) qcp_pos = ene_header%arrays
!now that we know how many arrays we need -> allocate them in EQ_save for file writing
        allocate(EQ_save(nstates))
	do i=1,nstates
	allocate(EQ_save(i)%qx(ene_header%arrays),EQ_save(i)%qq(ene_header%arrays),&
		EQ_save(i)%qp(ene_header%arrays),EQ_save(i)%qw(ene_header%arrays),&
		EQ_save(i)%total(ene_header%arrays))
	end do

!if using PBC, alredy allocate the old_EQ arrays so we don't have to do it
!every time we update the volume
if( use_PBC ) then
        allocate(old_EQ(nstates))
end if

if (use_qcp) call qcp_init

!TODO write QCP info to header

! the last thing to prepare are the MD grids
! so we call them here
! they will be filled later when the pair lists are created the first time
#ifdef USE_GRID
call create_grid_pp
if (nwat .gt. 0) then
call create_grid_pw
call create_grid_ww
end if

! now allocate the array needed to store the grid lookup information for
! the charge groups
allocate(pp_igrid(ncgp_solute))
allocate(pw_igrid(ncgp_solute))
allocate(ww_igrid(nwat))
#endif

! prepare internal nonbonded interactions for solvent
call prep_sim_precompute_solvent

#if defined (USE_MPI)
!	reclength=nstates*((2*ene_header%arrays)+(2*ene_header%arrays))
        reclength= nstates*(2+2)
        !qp interactions and qw interactions give 4*nstates real
#endif
end subroutine prep_sim

!-----------------------------------------------------------------------

subroutine prep_sim_precompute_solvent
! locals
integer,allocatable		:: interaction(:,:)
integer				:: i,j,na,nb,counter,last
real(kind=prec)			:: tempA,tempB
type(NB_TYPE),allocatable	:: tmp_solv_int(:)

if (ntors .gt. ntors_solute) then
! we need a way to troll the solvent interactions without having to build a list for the 1-4
! interactions and save it in the topology, this could kill topology reading otherwise
! so we start trolling them here before shake and other stuff is assigned
! This means we know all solvent internal stuff
! the maximum possible number is n^2, but we exclude self interaction and nearest neighbors

! fill array with all possible interaction codes
allocate(interaction(solv_atom,solv_atom),stat=alloc_status)
call check_alloc('solvent self interaction array')
do i = 1 , solv_atom
do j = 1 , solv_atom
interaction(i,j) = ljcod(iac(nat_solute+i),iac(nat_solute+j))
end do
end do

! first, exclude all self interactions
do i = 1 , solv_atom
interaction(i,i) = 0
end do

! need to know kast atom of first solvent molecule, because we can stop searching there :)
last = nat_solute + solv_atom

! now, troll the bond list for interactions to exclude
do i = nbonds_solute+1 , nbonds
na = bnd(i)%i
nb = bnd(i)%j
if (na.gt.last .or. nb.gt.last) cycle
! get the actual index of the solvent atom
! and set this index to zero
na = na - nat_solute
nb = nb - nat_solute
interaction(na,nb) = 0
interaction(nb,na) = 0
end do
! same for angle list
do i = nangles_solute+1, nangles
na = ang(i)%i
nb = ang(i)%k
if (na.gt.last .or. nb.gt.last) cycle
na = na - nat_solute
nb = nb - nat_solute
interaction(na,nb) = 0
interaction(nb,na) = 0
end do
! different now for torsion
! set to type 3 (1-4 interaction) if not previously set to zero
! and set inverted interaction to zero to prevent double counting
do i = ntors_solute+1, ntors
na = tor(i)%i
nb = tor(i)%l
if (na.gt.last .or. nb.gt.last) cycle
na = na - nat_solute
nb = nb - nat_solute
if ((interaction(na,nb) .ne. 0) .or. (interaction(nb,na) .ne. 0) ) then
interaction(na,nb) = 3
interaction(nb,na) = 0
end if
end do
counter = 0
! now we know the remaining interactions that we need to precalculate
! allocate full size array and later shrink it to the number needed
allocate(tmp_solv_int(solv_atom**2),stat=alloc_status)
call check_alloc('tmp solvent internal nonbonds')
do i = 1, solv_atom
do j = 1, solv_atom
if (interaction(i,j) .ne. 0 ) then
counter = counter + 1
tmp_solv_int(counter)%i=i
tmp_solv_int(counter)%j=j
select case(ivdw_rule)
case(VDW_GEOMETRIC)
tmp_solv_int(counter)%vdWA = aLJ_solv(i,interaction(i,j)) * aLJ_solv(j,interaction(i,j))
tmp_solv_int(counter)%vdWB = bLJ_solv(i,interaction(i,j)) * bLJ_solv(j,interaction(i,j))
case(VDW_ARITHMETIC)
tempA = aLJ_solv(i,interaction(i,j)) + aLJ_solv(j,interaction(i,j))
tempA = tempA**2
tempA = tempA * tempA * tempA
tempB = bLJ_solv(i,interaction(i,j)) * bLJ_solv(j,interaction(i,j))

tmp_solv_int(counter)%vdWA = (tempA**2) * tempB
tmp_solv_int(counter)%vdWB = 2.0_prec * tempA * tempB

end select

tmp_solv_int(counter)%elec = chg_solv(i) * chg_solv(j)
if (interaction(i,j) .eq. 3) &
tmp_solv_int(counter)%elec = tmp_solv_int(counter)%elec * el14_scale
end if
end do
end do
if (counter .le. 0) then
num_solv_int = 1
tmp_solv_int(:)%i=-1
tmp_solv_int(:)%j=-1
tmp_solv_int(:)%vdWA=-1E35_prec
tmp_solv_int(:)%vdWB=-1E35_prec
tmp_solv_int(:)%elec=-1E35_prec
else
num_solv_int = counter
end if

! now make real array to shrink done to the actual number of interactions
allocate(nonbnd_solv_int(num_solv_int),stat=alloc_status)
call check_alloc('solvent internal nonbonds')
nonbnd_solv_int(1:num_solv_int) = tmp_solv_int(1:num_solv_int)

! done, clean up
deallocate(tmp_solv_int)
deallocate(interaction)

else
! needed to prevent MPI from crashing during broadcast
num_solv_int = 1
allocate(nonbnd_solv_int(num_solv_int),stat=alloc_status)
call check_alloc('solvent internal nonbonds')
nonbnd_solv_int(:)%i=-1
nonbnd_solv_int(:)%j=-1
nonbnd_solv_int(:)%vdWA=-1E35_prec
nonbnd_solv_int(:)%vdWB=-1E35_prec
nonbnd_solv_int(:)%elec=-1E35_prec
end if


end subroutine prep_sim_precompute_solvent

!-----------------------------------------------------------------------

subroutine prep_sim_version(version)
! arguments
character(*)	:: version

! local values
integer :: canary = 1337

if (prec .eq. singleprecision) then
canary = 137
else if (prec .eq. doubleprecision) then
canary = 1337
#ifndef PGI
else if (prec .eq. quadprecision) then
canary = 13337
#endif
else
call die('No such precision')
end if
! hard coded value
!now write version number so Miha stops complaining
ene_header%version = trim(version)

if ( iene_cycle .gt. 0 ) then
write (11) canary,ene_header%arrays,ene_header%totresid,ene_header%types(1:ene_header%arrays),&
ene_header%numres(1:ene_header%arrays),ene_header%resid(1:ene_header%totresid),&
ene_header%gcnum(1:ene_header%arrays),ene_header%version
end if

end subroutine prep_sim_version

!-----------------------------------------------------------------------
!Put molecules back in box for nice visualisation.
!Change boxcentre if rigid_box_centre is off.
!Update cgp_centers for LRF.
!-----------------------------------------------------------------------
subroutine put_back_in_box

TYPE(qr_vec)				::	boxc,old_boxc
integer				::	i, j, starten, slutet
!the borders of the periodic box
TYPE(qr_vec)                            :: maxn,minn
TYPE(qr_vec)                            :: cm
integer				::  k, ig
integer				::  pbib_start, pbib_stop

!the function can be run exclusivly on node 0
!no gain from doing it on each node
if (nodeid.eq.0) then
old_boxc = boxcentre
if( .not. rigid_box_centre ) then !if the box is allowed to float around, center on solute if present, otherwise center on solvent
        if( nat_solute > 0) then
                slutet = nat_solute
                !starten = ncgp_solute + 1
                
        else !if no solute present, centre box around solvent
                slutet = natom
                !starten = 1
        end if
        !find center
        boxc = boxc * zero
        do i = 1,slutet
                boxc = boxc + x(i)
        end do
        boxc = boxc / real(slutet, kind=prec)
        boxcentre = boxc !store new boxcentre
        ! starten = ncgp_solute + 1 !solute can not move around the corner
else
        !use boxcenter given in topology, ie. 'moving' solute
        boxc = boxcentre
        !starten = 1
end if

!calculate the borders of the periodic box
maxn = boxc + boxlength/2.0_prec
minn = boxc - boxlength/2.0_prec

mvd_mol(:) = 0 

!pbib_start and pbib_stop are the starting and stopping molecule indexes of which molecules to Put Back In Box
pbib_start = 1
pbib_stop  = nmol
if ( .not. put_solute_back_in_box ) then !we're not putting solute back in box
!we now have nsolute to keep track of this stuff, don't use this old way
	pbib_start = nsolute + 1
end if
if ( .not. put_solvent_back_in_box ) then !we're not putting solvent back in box
!same here, we know the molecule numbers
	pbib_stop = nsolute
end if          


do i=pbib_start,pbib_stop
        cm = cm * zero
        do j = istart_mol(i),istart_mol(i+1)-1 !loop over all atoms in molecule i
                cm = cm + x(j)*mass(j)
        end do
        cm = cm * mol_mass(i) !centre of mass of molecule i
!mol_mass is actually 1/mol_mass
        !x-direction
        if( cm%x .gt. maxn%x) then !position of centre of mass

                        do j= istart_mol(i),istart_mol(i+1)-1 !move the molecule
						x(j)%x = x(j)%x - boxlength%x
				        end do
						mvd_mol(i) = 1
		else if ( cm%x .lt. minn%x ) then
                        do j= istart_mol(i),istart_mol(i+1)-1 !move the molecule
						x(j)%x = x(j)%x + boxlength%x
				        end do
						mvd_mol(i) = 1
        end if

       ! y-direction
        if( cm%y .gt. maxn%y) then !position of centre of mass
                        do j= istart_mol(i),istart_mol(i+1)-1 !move the molecule
						x(j)%y = x(j)%y - boxlength%y
				        end do
						mvd_mol(i) = 1
        else if ( cm%y .lt. minn%y ) then
                        do j= istart_mol(i),istart_mol(i+1)-1 !move the molecule
						x(j)%y = x(j)%y + boxlength%y
				        end do
						mvd_mol(i) = 1
        end if

        !z-direction
        if( cm%z .gt. maxn%z) then !position of centre of mass
                        do j= istart_mol(i),istart_mol(i+1)-1 !move the molecule
						x(j)%z = x(j)%z - boxlength%z
				        end do
						mvd_mol(i) = 1
        else if ( cm%z .lt. minn%z ) then
                        do j= istart_mol(i),istart_mol(i+1)-1 !move the molecule
						x(j)%z = x(j)%z + boxlength%z
				        end do
						mvd_mol(i) = 1
  		end if
end do !over molecules
end if !if(nodeid .eq. 0)
!now that node 0 has done the job, broadcast the changes to the slaves
#if defined(USE_MPI)
!only during MD, skip for first setup, mpitypes are 
!set after this, so save to use as proxy
if(mpitype_batch_lrf.ne.0) then
	call MPI_Bcast(x, nat3, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
	if (ierr .ne. 0) call die('init_nodes/MPI_Bcast x')
	call MPI_Bcast(boxcentre, 3, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
	if (ierr .ne. 0) call die('init_nodes/MPI_Bcast x')
	call MPI_Bcast(old_boxc, 3, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
	if (ierr .ne. 0) call die('init_nodes/MPI_Bcast x')
end if
#endif

#ifdef USE_GRID
! we now need to update the coordinates of all grid boxes
! to shift by the amount the new box has shifted
grid_pp(:)%x = grid_pp(:)%x - (old_boxc-boxcentre)
grid_pp(:)%y = grid_pp(:)%y - (old_boxc-boxcentre)
grid_pp(:)%z = grid_pp(:)%z - (old_boxc-boxcentre)

grid_pp(:)%xend = grid_pp(:)%xend - (old_boxc-boxcentre)
grid_pp(:)%yend = grid_pp(:)%yend - (old_boxc-boxcentre)
grid_pp(:)%zend = grid_pp(:)%zend - (old_boxc-boxcentre)

grid_pw(:)%x = grid_pw(:)%x - (old_boxc-boxcentre)
grid_pw(:)%y = grid_pw(:)%y - (old_boxc-boxcentre)
grid_pw(:)%z = grid_pw(:)%z - (old_boxc-boxcentre)

grid_pw(:)%xend = grid_pw(:)%xend - (old_boxc-boxcentre)
grid_pw(:)%yend = grid_pw(:)%yend - (old_boxc-boxcentre)
grid_pw(:)%zend = grid_pw(:)%zend - (old_boxc-boxcentre)

grid_ww(:)%x = grid_ww(:)%x - (old_boxc-boxcentre)
grid_ww(:)%y = grid_ww(:)%y - (old_boxc-boxcentre)
grid_ww(:)%z = grid_ww(:)%z - (old_boxc-boxcentre)

grid_ww(:)%xend = grid_ww(:)%xend - (old_boxc-boxcentre)
grid_ww(:)%yend = grid_ww(:)%yend - (old_boxc-boxcentre)
grid_ww(:)%zend = grid_ww(:)%zend - (old_boxc-boxcentre)

#endif
!LRF: if molecule moved update all cgp_centers from first charge group to the last one
if (use_LRF) then

if (nodeid .eq.0) then
do k=pbib_start,pbib_stop
		if (mvd_mol(k) == 1) then  
						do ig=iwhich_cgp(istart_mol(k)),iwhich_cgp(istart_mol(k+1)-1)
								lrf(ig)%cgp_cent = lrf(ig)%cgp_cent * zero
								do i  = cgp(ig)%first, cgp(ig)%last
										lrf(ig)%cgp_cent = lrf(ig)%cgp_cent + x(cgpatom(i))
								end do
						lrf(ig)%cgp_cent = lrf(ig)%cgp_cent/real(cgp(ig)%last - cgp(ig)%first +1, kind=prec)
						end do
		end if 
end do
end if ! nodeid .eq. 0
#if defined(USE_MPI)
!for now just bcast the whole thing
!only during MD, skip for first setup, mpitypes are 
!set after this, so save to use as proxy
if(mpitype_batch_lrf.ne.0) then
	call MPI_Bcast(lrf,ncgp,mpitype_batch_lrf,0,MPI_COMM_WORLD, ierr)
	if (ierr .ne. 0) call die('put_back_in_box MPI_Bcast lrf')
end if
! for later use of PBC update on each node
!call lrf_gather
#endif
end if

end subroutine put_back_in_box

!-----------------------------------------------------------------------


real(kind=prec) function randm (ig)
! arguments
integer					::	ig

! local variables
integer, parameter		::	m = 100000000
integer, parameter		::  m1 = 10000
integer, parameter		::	mult=31415821
integer					::	irandh,irandl,multh,multl
real(kind=prec)					::	r
integer, save				::	irand = 0
integer, save				::  new = 0

if (new .eq. 0) then
new = 1
irand = mod (iabs(ig),m)
end if

! --- multiply irand by mult, but take into account that overflow must
! --- be discarded, and do not generate an error.
irandh = irand / m1
irandl = mod(irand, m1)
multh = mult / m1
multl = mod(mult, m1)

irand = mod(irandh*multl + irandl*multh, m1) * m1 + irandl*multl
irand = mod(irand + 1, m)

! --- convert irand to a real random number between 0 and 1.
r = real(irand / 10, kind=prec) * 10 / real(m, kind=prec)
if ((r .le. 0.e0_prec) .or. (r .gt. 1.e0_prec)) r = 0.e0_prec
randm = r
ig = irand

end function randm

!-----------------------------------------------------------------------

subroutine shrink_topology
!get rid of bonds and angles where all atoms are excluded
!or where the code has been set to 0 due to q-[bonds|angles|...]

!locals
integer						::	i, removed

if(exclude_bonded) then
        call centered_heading &
                ('Eliminating torsions & impropers for excluded atoms', '-')
end if

10	format('Reduced number of ',a,t31,'from ',i8,' to ')
12	format(i8)

i = 1
removed = 0
do while(i <= nbonds)
        !if all atoms excluded
        if(bnd(i)%cod <= 0) then
                !bond code either 0 (bond redefined in FEP file) 
                !or -1 (bond removed by shake)
                if(i <= nbonds_solute) then
                        bnd(i) = bnd(nbonds_solute)
                        bnd(nbonds_solute) = bnd(nbonds)
                        nbonds_solute = nbonds_solute - 1
                else
                        bnd(i) = bnd(nbonds)
                endif
                nbonds = nbonds - 1
                cycle !don't change i now, 
        end if
        i = i + 1 
end do

i = 1
do while(i <= nangles)
        !if all atoms excluded
        if(ang(i)%cod == 0) then
                !move last angle to current position
                if(i <= nangles_solute) then
                        ang(i) = ang(nangles_solute)
                        ang(nangles_solute) = ang(nangles)
                        nangles_solute = nangles_solute - 1
                else
                        ang(i) = ang(nangles)
                endif
                nangles = nangles - 1
                cycle !don't change i now, 
        end if
        i = i + 1 
end do

if(exclude_bonded) write(*,10, advance='no') 'torsions', ntors
i = 1
do while(i <= ntors)
        !if all atoms excluded
        if((exclude_bonded .and. excl(tor(i)%i) .and. excl(tor(i)%j) &
                .and. excl(tor(i)%k) .and. excl(tor(i)%l)) &
                .or. tor(i)%cod == 0) then
                !move last bond to current position
                if(i <= ntors_solute) then
                        tor(i) = tor(ntors_solute)
                        tor(ntors_solute) = tor(ntors)
                        ntors_solute = ntors_solute - 1
                else
                        tor(i) = tor(ntors)
                endif
                ntors = ntors - 1
                cycle !don't change i now, 
        end if
        i = i + 1 
end do
if(exclude_bonded) write(*, 12) ntors

if(exclude_bonded) write(*,10, advance='no') 'impropers', nimps
i = 1
do while(i <= nimps)
        !if all atoms excluded
        if(exclude_bonded .and. excl(imp(i)%i) .and. excl(imp(i)%j) &
                .and. excl(imp(i)%k) .and. excl(imp(i)%l) &
                .or. imp(i)%cod == 0) then
                if(i <= nimps_solute) then
                        imp(i) = imp(nimps_solute)
                        imp(nimps_solute) = imp(nimps)
                        nimps_solute = nimps_solute - 1
                else
                        imp(i) = imp(nimps)
                endif
                nimps = nimps - 1
                cycle !don't change i now, 
        end if
        i = i + 1 
end do
if(exclude_bonded) write(*, 12) nimps

end subroutine shrink_topology

!--------------------------------------------------------------------


subroutine stop_cm_translation
! local variables
integer						::	i
real(kind=prec)						::	rmass,totmass
TYPE(qr_vec)						::	vcm

! calculate totmass and vcm
totmass = zero
vcm = vcm * zero
do i=1,natom
rmass = iaclib(iac(i))%mass
totmass=totmass+rmass
vcm = vcm +  (v(i) * rmass)
end do

! scale vcm
vcm = vcm / totmass

! update v
do i=1,natom
v(i) = v(i) - vcm
end do
end subroutine stop_cm_translation

!-----------------------------------------------------------------------
subroutine topology
! local variables
integer					::	nat3
integer					::	i
real(kind=prec)					::	box_min, vtemp, vtemp1

!
! read topology
!
! will init:
!  natom
!  lots of stuff from topo_load
!  nwat
!  anglib, torlib, implib (conversion)
!  ljcod
!  [iaclib%bvdw] (conversion)

if(.not. topo_load(top_file, require_version=4.15_prec)) then
        call die('Failed to load topology.')
end if
natom = nat_pro

nwat = (natom - nat_solute) / solv_atom
!add extra element to molecule start atom array to keep track of last atom
istart_mol(nmol+1) = nat_pro + 1

! abort if no atoms
if (natom .eq. 0) call die('zero particles to simulate')

! convert libraries from degrees to radians
anglib(1:nangcod)%ang0 = deg2rad*anglib(1:nangcod)%ang0
torlib(1:ntorcod)%paths = one/torlib(1:ntorcod)%paths
torlib(1:ntorcod)%deltor = deg2rad*torlib(1:ntorcod)%deltor
implib(1:nimpcod)%imp0 = deg2rad*implib(1:nimpcod)%imp0

num_atyp = maxval(iac)

allocate(ljcod(num_atyp,num_atyp))

ljcod(:,:) = 1
do i=1,nlj2
ljcod(lj2(i)%i, lj2(i)%j) = 2
ljcod(lj2(i)%j, lj2(i)%i) = 2
end do


!	If arithmetic combination rule (ivdw_rule=2) take sqrt(epsilon) now
if ( ivdw_rule .eq. 2 ) then
do i=1,natyps
  iaclib(i)%bvdw(:) = sqrt(abs(iaclib(i)%bvdw(:)))
end do
end if

!check if same boundary in topology and Qdyn input file
if(  ((.not. use_PBC) .and.  box) .or. ( use_PBC .and. (.not. box)) ) then
call die('Must have same boundary (sphere or box) in topology and input file')
end if

if( use_PBC ) then
if( (boxlength%x == zero .or. boxlength%y == zero .or. boxlength%z == zero ) ) then
        inv_boxl = inv_boxl * zero
else
        inv_boxl = one/boxlength
end if

!check cut-offs if periodic box used
box_min = min( boxlength%x, boxlength%y, boxlength%z )

!Solute-solute cut-off radii
if( .not. (box_min .gt. Rcpp*2) ) then
        call die('Solute-solute cut-off radii too large')

!Solvent-solvent
else if( .not. (box_min .gt. Rcww*2) ) then
        call die('Solvent-solvent cut-off radii too large')

!Solute-solvent
else if( .not. (box_min .gt. Rcpw*2) ) then
        call die('Solute-solvent cut-off radii too large')

!Q-atom
else if( .not. (box_min .gt. rcq*2) ) then
        call die('Q-atom cut-off radii too large')
!LRF
else if( .not. (box_min .gt. RcLRF*2) ) then
        call die('LRF cut-off radii too large')
end if
end if
end subroutine topology

!-----------------------------------------------------------------------


subroutine wat_sphere
! local variables
integer					:: i,i3,kr,isort,int_wat,istate
real(kind=prec)					:: rc,rnwat
real(kind=prec)					:: crgexcl

!possibly override target sphere radius from topology
if(rwat_in > zero) rwat = rwat_in


!calc. total charge of non-excluded non-Q-atoms and excluded atoms
crgtot = zero
crgexcl = zero
rho_wat = topo_rho_wat
do i = 1, nat_solute
	if ( .not. excl(i) ) then
		if ( iqatom(i)==0 ) then
			crgtot = crgtot + crg(i)
		end if
	else
		crgexcl = crgexcl + crg(i)
	end if
end do
write (*,60) 'non-Q atoms', crgtot
write (*,60) 'excluded atoms', crgexcl
60 format ('Total charge of ',a,t41,'= ',f10.2)

!calc effective charge of simulation sphere at this lambda
crgQtot = zero
do i = 1, nqat
	do istate = 1, nstates
		crgtot = crgtot + qcrg(i,istate)*EQ(istate)%lambda
		crgQtot = crgQtot + qcrg(i,istate)*EQ(istate)%lambda
	end do
end do
write (*,70) crgtot
70 format ('Total charge of system                  = ',f10.2)

if (.not. wpol_born) crgQtot = zero !ignore total Q if Born corr. is off
if ( nwat .eq. 0 ) return


!	Set default values for unspecified optional parameters
if(fk_wsphere == -1) then
        !
        ! To be replaced by function of rc giving appropriate default for any sphere
        !
        fk_wsphere = fk_wsphere_default
end if
if(fkwpol == -1) then
        !
        ! To be replaced by function of rc giving appropriate default for any sphere
        !
        fkwpol = fkwpol_default
end if
if(Dwmz == -1) then !Use magic function to get suitable Dwmz
        Dwmz = 0.26_prec*exp(-0.19_prec*(rwat-15.0_prec))+0.74_prec
end if
if(awmz == -1) then !use magic for the reach of the Morse potential
        awmz = 0.2_prec/(one+exp(0.4_prec*(rwat-25.0_prec)))+0.3_prec
end if

write (*,90) rwat, fk_wsphere, Dwmz, awmz
90	format ('Target water sphere radius              = ',f10.2,/,&
                'Surface inward harmonic force constant  = ',f10.2,/,&
                'Surface attraction well depth           = ',f10.2,/,&
                'Surface attraction well width           = ',f10.2)
92	format ('Water polarisation restraints           : ',a)
if(.not. wpol_restr) then
        write(*,92) 'OFF'
else if(wpol_born) then
        write(*,92) 'ON, Born correction enabled'
        write(*, 100) fkwpol
else 
        write(*,92) 'ON, Born correction disabled'
        write(*, 100) fkwpol
end if
100	format('Radial polarisation force constant      = ',f10.2)

end subroutine wat_sphere

!-----------------------------------------------------------------------

subroutine wat_shells
! set up the shells for polarisation restraining

! local variables
real(kind=prec)						::	rout, dr, ri, Vshell, rshell, drs, eps_diel
integer						::	is, n_insh


integer						::	nwpolr_shell_restart, filestat
integer						::	bndcodw, angcodw

logical :: old_restart = .false.
real(8),allocatable :: x_old(:), v_old(:)
real(kind=singleprecision),allocatable :: x_1(:), v_1(:)
real(kind=doubleprecision),allocatable :: x_2(:), v_2(:)
#ifndef PGI
real(kind=quadprecision),allocatable :: x_3(:), v_3(:)
#endif
real(8) :: old_boxlength(3),old_boxcentre(3)
real(kind=singleprecision) :: boxlength_1(3),boxcentre_1(3)
real(kind=doubleprecision) :: boxlength_2(3),boxcentre_2(3)
#ifndef PGI
real(kind=quadprecision) :: boxlength_3(3),boxcentre_3(3)
#endif
integer :: headercheck,i,myprec,dummy

if (prec .eq. singleprecision) then
myprec = -137
elseif (prec .eq. doubleprecision) then
myprec = -1337
#ifndef PGI
elseif (prec .eq. quadprecision) then
myprec = -13337
#endif
else
call die('No such precision')
end if

!calc mu_w
!look up bond code for water
bndcodw = bnd(nbonds)%cod
angcodw = ang(nangles)%cod
!find charge of water O = charge of 1st solvent atom

mu_w = -chg_solv(1)*bondlib(bndcodw)%bnd0*cos(anglib(angcodw)%ang0/2)

! shell widths are drout, 2drout, 3drout
drs = wpolr_layer / drout !number of drouts 

! calc number of shells based on arithmetic series sum formula
nwpolr_shell = int(-0.5_prec + sqrt(2*drs + 0.25_prec)) 
allocate(wshell(nwpolr_shell), stat=alloc_status)
call check_alloc('water polarisation shell array')

write(*, 100) nwpolr_shell
100	format(/,'Setting up ', i1, ' water shells for polarisation restraints.')

if(restart) then !try to load theta_corr from restart file
!first, rewind file to find out if we have an old or new restart
   rewind(2)
   read(2) headercheck
   if ((headercheck .ne. -137).and.(headercheck .ne. -1337).and.(headercheck .ne. -13337)) then
     old_restart = .true.
     rewind(2)
     allocate(x_old(3*natom),v_old(3*natom))
     read(2) dummy,(x_old(i),i=1,nat3)
     read(2) dummy,(v_old(i),i=1,nat3)
     if( use_PBC) then
       read(2) old_boxlength(:)
       read(2) old_boxcentre(:)
     end if
     deallocate(x_old,v_old)
   else if (headercheck .eq. -137) then
     allocate(x_1(3*natom),v_1(3*natom))
     read(2) dummy,(x_1(i),i=1,nat3)
     read(2) dummy,(v_1(i),i=1,nat3)
     if( use_PBC) then
       read(2) boxlength_1(:)
       read(2) boxcentre_1(:)
     end if
     deallocate(x_1,v_1)
   else if (headercheck .eq. -1337) then
     allocate(x_2(3*natom),v_2(3*natom))
     read(2) dummy,(x_2(i),i=1,nat3)
     read(2) dummy,(v_2(i),i=1,nat3)
     if( use_PBC) then
       read(2) boxlength_2(:)
       read(2) boxcentre_2(:)
     end if
     deallocate(x_2,v_2)
   else if (headercheck .eq. -13337) then
#ifndef PGI
     allocate(x_3(3*natom),v_3(3*natom))
     read(2) dummy,(x_3(i),i=1,nat3)
     read(2) dummy,(v_3(i),i=1,nat3)
     if( use_PBC) then
       read(2) boxlength_3(:)
       read(2) boxcentre_3(:)
     end if
     deallocate(x_3,v_3)
#else
     call die('Quadruple precision not supported in PGI')
#endif
   end if
        read(2, iostat=filestat) nwpolr_shell_restart
        if(filestat .ne. 0 .or. nwpolr_shell_restart /= nwpolr_shell) then
                write(*,102) 
                wshell(:)%theta_corr = zero
        else
                backspace(2)
                if (old_restart) then
                   allocate(old_wshell(nwpolr_shell), stat=alloc_status)
                   read(2) nwpolr_shell_restart,old_wshell(:)%theta_corr
                   wshell(:)%theta_corr = old_wshell(:)%theta_corr
                   write(*,*) 'Loaded water polarization data from old qdyn restart file'
                   deallocate(old_wshell)
                else
                   if (headercheck .ne. myprec) then
                   write(*,*) '>>> WARNING: Using water polarisation from restart with mismatched precision'
                     if (headercheck .eq. -137) then
                       allocate(wshell_single(nwpolr_shell), stat=alloc_status)
                       read(2) nwpolr_shell_restart,wshell_single(:)%theta_corr
                       wshell(:)%theta_corr = wshell_single(:)%theta_corr
                       deallocate(wshell_single)
                     else if (headercheck .eq. -1337) then
                       allocate(wshell_double(nwpolr_shell), stat=alloc_status)
                       read(2) nwpolr_shell_restart,wshell_double(:)%theta_corr
                       wshell(:)%theta_corr = wshell_double(:)%theta_corr
                       deallocate(wshell_double)
                     else if (headercheck .eq. -13337) then
#ifndef PGI
                       allocate(wshell_quad(nwpolr_shell), stat=alloc_status)
                       read(2) nwpolr_shell_restart,wshell_quad(:)%theta_corr
                       wshell(:)%theta_corr = wshell_quad(:)%theta_corr
                       deallocate(wshell_quad)
#else
                       call die('Quadruple precision not supported in PGI')
#endif
                     end if
                   else
                     read(2) nwpolr_shell_restart,wshell(:)%theta_corr
                   end if
                end if
                write(*,103)
        end if
else
        wshell(:)%theta_corr = zero
end if

102	format('>>> WARNING: Failed to read polarisation restraint data from restart file.')
103	format('Loaded polarisation restraint data from restart file.')

write(*,'(a)') 'Shell #    outer radius    inner radius'
110	format(i7, 2f16.2)

rout = rwat
n_max_insh = 0
do is = 1, nwpolr_shell 
        wshell(is)%avtheta = 0
        wshell(is)%avn_insh = 0
        wshell(is)%rout = rout
    dr = drout*is
        ri = rout - dr
        wshell(is)%dr = dr
        Vshell = rout**3 - ri**3
        n_insh = int(4 * pi/3 * Vshell * rho_wat)
        if (n_insh > n_max_insh) n_max_insh = n_insh
rshell = (0.5_prec*(rout**3+ri**3))**(one/3.0_prec)


        ! --- Note below: 0.98750 = (1-1/epsilon) for water
        ! now deprecated because we actually read the dielectric from the
        ! solvent files as dielectric = epsilon * 1000
        eps_diel = real(dielectric/1000,kind=prec)
        eps_diel = one-(one/eps_diel)
!wshell(is)%cstb = crgQtot*0.98750_prec/(rho_wat*mu_w*4.0_prec*pi*rshell**2)
        wshell(is)%cstb = crgQtot*eps_diel/(rho_wat*mu_w*4.0_prec*pi*rshell**2)
        write(*, 110) is, rout, ri
        rout = rout - dr
end do

n_max_insh = n_max_insh * 1.5_prec !take largest and add some extra
call allocate_watpol_arrays

end subroutine wat_shells

!----------------------------------------------------------------------------------------
!*******************************************************
!Will find and return the xtop atom number from 
!  residue number and atom number in residue from
!  library sequence.
! Uses global variables: xtop,nres,res
!*******************************************************

integer function get_atom_from_resnum_atnum(aid)
!arguments
character(*), intent(in)	::	aid	!string=residue:atom
	
!locals
integer						::	separator_pos
character(len=20)			::	res_str
character(len=5)			::	atom_str
integer						::	filestat
integer						::	resnum, atnum

get_atom_from_resnum_atnum = 0

separator_pos = scan(aid, ':')
if(separator_pos < 2 .or. separator_pos == len_trim(aid)) return !no valid colon found
res_str = aid(1:separator_pos-1)
atom_str = aid(separator_pos+1:len_trim(aid))
read(res_str, *, iostat=filestat) resnum
read(atom_str, *, iostat=filestat) atnum
if(filestat > 0) return

!Residue must be in topology
if(resnum < 1 .or. resnum > nres) then                     
  return                                                 
end if

if(atnum .le. (res(resnum+1)%start - res(resnum)%start)) then
  get_atom_from_resnum_atnum = res(resnum)%start + atnum - 1
return
end if

!we have an error: 
write(*, 120) atnum, resnum
call die('error in finding atom number from resnum:atnum.')

120	format('>>>>> ERROR: There is no atom number ',i4,' in residue ',i4,'.')
end function get_atom_from_resnum_atnum

!----------------------------------------------------------------------------------------

end module SIMPREP




