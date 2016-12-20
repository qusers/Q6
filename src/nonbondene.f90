! (C) 2014 Uppsala Molekylmekaniska HB, Uppsala, Sweden
! nonbondene.f90 
! based on md.f90
! by Johan Åqvist, John Marelius, Anders Kaplan, Isabella Feierberg, Martin Nervall & Martin Almlöf
! calculation of nonbonded energies and list generation
! by Paul Bauer

module NONBONDENE

! used modules
use SIZES
use QATOM
use QALLOC
use GLOBALS
use QMATH
use NONBONDED
use EXC
use MPIGLOB
!$ use omp_lib
implicit none


contains

subroutine cgp_centers
! *** local variables
integer						::	ig,i
#ifdef _OPENMP
integer :: quotient, remainder


threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (ncgp)/threads_num
remainder = MOD(ncgp, threads_num)
mp_start = thread_id * quotient + MIN(thread_id, remainder) + 1
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

do ig = mp_start,mp_end
#else
do ig = 1, ncgp
#endif
        lrf(ig)%cgp_cent  = zero
        lrf(ig)%phi0      = zero
        lrf(ig)%phi1      = zero
        lrf(ig)%phi2(:)   = zero
        lrf(ig)%phi3(:)   = zero

        do i  = cgp(ig)%first, cgp(ig)%last
                lrf(ig)%cgp_cent = lrf(ig)%cgp_cent + x(cgpatom(i))
        end do

lrf(ig)%cgp_cent = lrf(ig)%cgp_cent/real(cgp(ig)%last - cgp(ig)%first +1, kind=prec)

end do

end subroutine cgp_centers

subroutine distribute_nonbonds
!locals
integer					:: npp, npw, nqp, nww, nqw
type(NODE_ASSIGNMENT_TYPE),allocatable  :: node_assignment(:)
real					:: avgload, old_avgload
integer					:: i, last_cgp, last_pair
integer					:: mpitype_pair_assignment, mpitype_node_assignment
integer                                 :: average_pairs,inode,icgp,sum,less_than_sum
integer                                 :: n_bonded, n_nonbonded, master_assign
real                                    :: percent
integer                                 :: master_sum
!!!!Tmp vars för allokering
integer,parameter			:: vars = 5
integer     :: mpi_batch, i_loop
integer     :: blockcnt(vars),type(vars)
integer(8)  :: disp(vars)
!!!


! count the number of nonbonded interactions and distribute them among the nodes

if (nodeid .eq. 0) then
! nice header
call centered_heading('Distribution of charge groups','-')

!Allocate node_assignment
allocate(node_assignment(0:numnodes-1),stat=alloc_status)
call check_alloc('node_assignment')

!Allocate arrays that hold no. pairs per chargegroup.
call allocate_nbxx_per_cgp

!Count stuff for balancing nodes and allocating nonbond arrays nbxx()
nbpp_per_cgp = 0
nbpw_per_cgp = 0
nbqp_per_cgp = 0
nbqw_per_cgp = 0
nbww_per_cgp = 0

call nbpp_count(npp, nbpp_per_cgp,Rcpp**2)    !Only for switching atoms!!!?
call nbpw_count(npw, nbpw_per_cgp,Rcpw**2)
call nbqp_count(nqp, nbqp_per_cgp,Rcq**2)
call nbqw_count(nqw, nbqw_per_cgp,Rcq**2) 
call nbww_count(nww, nbww_per_cgp,Rcww**2) 

!For keeping track of actual # of nonbonded pairs
totnbpp = npp
totnbpw = npw
totnbww = nww*solv_atom**2
totnbqp = nqp
totnbqw = nqw*solv_atom*nqat

if (numnodes .eq. 1) then
! only one node: no load balancing

! make the master node handle everything
calculation_assignment%pp%start = 1
calculation_assignment%pp%end = ncgp_solute
calculation_assignment%pw%start = 1
calculation_assignment%pw%end = ncgp_solute
calculation_assignment%qp%start = 1
calculation_assignment%qp%end = ncgp_solute
calculation_assignment%qw%start = 1
calculation_assignment%qw%end = nwat
calculation_assignment%ww%start = 1
calculation_assignment%ww%end = nwat

calculation_assignment%shake%start        = 1
calculation_assignment%shake%end          = shake_molecules
calculation_assignment%natom%start        = 1
calculation_assignment%natom%end          = natom
calculation_assignment%nat_solute%start   = 1
calculation_assignment%nat_solute%end     = nat_solute
calculation_assignment%nat_solvent%start  = nat_solute + 1
calculation_assignment%nat_solvent%end    = natom
calculation_assignment%nmol%start         = 1
calculation_assignment%nmol%end           = nmol
calculation_assignment%ncgp%start         = 1
calculation_assignment%ncgp%end           = ncgp

#if defined (USE_MPI)
else ! i.e. slave nodes exists


! A simple solution to avoid parallelising the bonded
! Calculate n_bonded and n_nonbonded 
! Approximate time of computing one bonded with one nonbonded
! The number of qq-interactions are neglected
n_bonded = nbonds + nangles + ntors + nimps
n_nonbonded = totnbpp + totnbpw + totnbww + totnbqw + totnbqp

! Compare to determine how many nonbonded master should get
! A bonded is faster, so this favours an early completion for master
master_assign =  n_nonbonded/numnodes - n_bonded * numnodes

! calculate the assignments ********************

!Calculate balanced assignment for p-p pairs
icgp=0
sum=0
 !First assign master a small part
node_assignment(0)%pp%start=icgp+1
percent=REAL(totnbpp, kind=prec)/n_nonbonded
less_than_sum = master_assign*percent  ! No. of pp-type to assign master
do while((icgp .lt. ncgp_solute) .and. (sum .lt. less_than_sum))
   icgp=icgp+1
   sum=sum + nbpp_per_cgp(icgp)
end do
node_assignment(0)%pp%end=icgp
master_sum=sum
 !Now assign slaves
average_pairs=(totnbpp-sum)/(numnodes-1)
do inode=1,numnodes-2
  node_assignment(inode)%pp%start=icgp+1
  less_than_sum=average_pairs*inode+master_sum
  do while (sum .lt. less_than_sum) 
     icgp=icgp+1
     sum=sum + nbpp_per_cgp(icgp)
  end do
  node_assignment(inode)%pp%end=icgp
end do
node_assignment(numnodes-1)%pp%start=icgp+1
node_assignment(numnodes-1)%pp%end=ncgp_solute

!Calculate balanced assignment for p-w pairs
icgp=0
sum=0
node_assignment(0)%pw%start=icgp+1
percent=REAL(totnbpw, kind=prec)/n_nonbonded
less_than_sum = master_assign*percent
do while((icgp .lt. ncgp_solute) .and. (sum .lt. less_than_sum))
   icgp=icgp+1
   sum=sum + nbpw_per_cgp(icgp)
end do
node_assignment(0)%pw%end=icgp
master_sum=sum
average_pairs=(totnbpw-sum)/(numnodes-1)
do inode=1,numnodes-2
  node_assignment(inode)%pw%start=icgp+1
  less_than_sum=average_pairs*inode+master_sum
  do while (sum .lt. less_than_sum)
     icgp=icgp+1
     sum=sum + nbpw_per_cgp(icgp)
  end do
  node_assignment(inode)%pw%end=icgp
end do
node_assignment(numnodes-1)%pw%start=icgp+1
node_assignment(numnodes-1)%pw%end=ncgp_solute

!Calculate balanced assignment for q-p pairs
icgp=0
sum=0
node_assignment(0)%qp%start=icgp+1
percent=REAL(totnbqp, kind=prec)/n_nonbonded
less_than_sum = master_assign*percent
do while((icgp .lt. ncgp_solute) .and. (sum .lt. less_than_sum))
   icgp=icgp+1
   sum=sum + nbqp_per_cgp(icgp)
end do
node_assignment(0)%qp%end=icgp
master_sum=sum
average_pairs=(totnbqp-sum)/(numnodes-1)
do inode=1,numnodes-2
  node_assignment(inode)%qp%start=icgp+1
  less_than_sum=average_pairs*inode+master_sum
  do while(sum .lt. less_than_sum)
     icgp=icgp+1
     sum=sum + nbqp_per_cgp(icgp)
  end do
  node_assignment(inode)%qp%end=icgp
end do
node_assignment(numnodes-1)%qp%start=icgp+1
node_assignment(numnodes-1)%qp%end=ncgp_solute

!Calculate balanced assignment for w-w pairs
icgp=0
sum=0
node_assignment(0)%ww%start=icgp+1
percent=REAL(totnbww, kind=prec)/n_nonbonded
less_than_sum = master_assign*percent
do while((icgp .lt. nwat) .and. (sum .lt. less_than_sum))
   icgp=icgp+1
   sum=sum + nbww_per_cgp(icgp)
end do
node_assignment(0)%ww%end=icgp
master_sum=sum
average_pairs=(totnbww-sum)/(numnodes-1)
do inode=1,numnodes-2
  node_assignment(inode)%ww%start=icgp+1
  less_than_sum=average_pairs*inode+master_sum
  do while(sum .lt. less_than_sum)
     icgp=icgp+1
     sum=sum + nbww_per_cgp(icgp)
  end do
  node_assignment(inode)%ww%end=icgp
end do
node_assignment(numnodes-1)%ww%start=icgp+1
node_assignment(numnodes-1)%ww%end=nwat

!Calculate balanced assignment for q-w pairs
icgp=0
sum=0
node_assignment(0)%qw%start=icgp+1
percent=REAL(totnbqw, kind=prec)/n_nonbonded
less_than_sum = master_assign*percent
do while((icgp .lt. nwat) .and. (sum .lt. less_than_sum))
   icgp=icgp+1
   sum=sum + nbqw_per_cgp(icgp)
end do
node_assignment(0)%qw%end=icgp
master_sum=sum
average_pairs=(totnbqw-sum)/(numnodes-1)
do inode=1,numnodes-2
  node_assignment(inode)%qw%start=icgp+1
  less_than_sum=average_pairs*inode+master_sum
  do while(sum .lt. less_than_sum)
     icgp=icgp+1
     sum=sum + nbqw_per_cgp(icgp)
  end do
  node_assignment(inode)%qw%end=icgp
end do
node_assignment(numnodes-1)%qw%start=icgp+1
node_assignment(numnodes-1)%qw%end=nwat

!now search the shake constrained molecules and give a balanced assignment for
!each node, based on the number of shake constraints in each searched molecule
icgp=0
sum=0
average_pairs=REAL((shake_constraints/numnodes),kind=prec)
do inode=0,numnodes-2
        node_assignment(inode)%shake%start=icgp+1
        less_than_sum = (inode+1)*average_pairs
        do while ((sum .lt. less_than_sum).and.(icgp.lt.shake_molecules))
             icgp=icgp+1
             sum=sum + shake_mol(icgp)%nconstraints
        end do
        node_assignment(inode)%shake%end=icgp
end do
node_assignment(numnodes-1)%shake%start=icgp+1
node_assignment(numnodes-1)%shake%end=shake_molecules

!assignment of atom numbers
sum=0
average_pairs=REAL((natom/numnodes),kind=prec)
do inode = 0,numnodes-2
        node_assignment(inode)%natom%start = sum + 1
        less_than_sum = (inode+1)*average_pairs
        do while (sum.lt.less_than_sum)
                sum = sum + average_pairs
        end do
        node_assignment(inode)%natom%end   = sum
end do
node_assignment(numnodes-1)%natom%start = node_assignment(numnodes-2)%natom%end + 1
node_assignment(numnodes-1)%natom%end   = natom

!assignment of solute atoms
sum=0
average_pairs=REAL((nat_solute/numnodes),kind=prec)
do inode = 0,numnodes-2
        node_assignment(inode)%nat_solute%start = sum + 1
        less_than_sum = (inode+1)*average_pairs
        do while (sum.lt.less_than_sum)
                sum = sum + average_pairs
        end do
        node_assignment(inode)%nat_solute%end   = sum
end do
node_assignment(numnodes-1)%nat_solute%start = & 
node_assignment(numnodes-2)%nat_solute%end + 1
node_assignment(numnodes-1)%nat_solute%end   = nat_solute

!assignment of solvent atoms
sum=nat_solute
average_pairs=REAL(((natom-nat_solute+1)/numnodes),kind=prec)
do inode = 0,numnodes-2
        node_assignment(inode)%nat_solvent%start = sum + 1
        less_than_sum = (inode+1)*average_pairs + nat_solute
        do while (sum.lt.less_than_sum)
                sum = sum + average_pairs
        end do
        node_assignment(inode)%nat_solvent%end   = sum
end do
node_assignment(numnodes-1)%nat_solvent%start = & 
node_assignment(numnodes-2)%nat_solvent%end + 1
node_assignment(numnodes-1)%nat_solvent%end   = natom

sum=0
average_pairs=REAL((nmol/numnodes),kind=prec)
do inode = 0, numnodes-2
        node_assignment(inode)%nmol%start = sum + 1
        less_than_sum = (inode+1)*average_pairs
        do while (sum.lt.less_than_sum)
                sum = sum + average_pairs
        end do
        node_assignment(inode)%nmol%end   = sum
end do
node_assignment(numnodes-1)%nmol%start = &   
node_assignment(numnodes-2)%nmol%end + 1
node_assignment(numnodes-1)%nmol%end   = nmol

sum=0
average_pairs=REAL((ncgp/numnodes),kind=prec)
do inode = 0, numnodes-2
        node_assignment(inode)%ncgp%start = sum + 1
        less_than_sum = (inode+1)*average_pairs
        do while (sum.lt.less_than_sum)
                sum = sum + average_pairs
        end do
        node_assignment(inode)%ncgp%end   = sum
end do
node_assignment(numnodes-1)%ncgp%start = &    
node_assignment(numnodes-2)%ncgp%end + 1
node_assignment(numnodes-1)%ncgp%end   = ncgp


#endif
end if    !if (numnodes .eq. 1)

! deallocate bookkeeping arrays
!deallocate(nppcgp, npwcgp, nqpcgp, nwwmol)

end if   !if (nodeid .eq. 0)

! distribute assignments to the nodes
#if defined (USE_MPI)
if (numnodes .gt. 1) then
    if (nodeid .ne. 0) then
	! Dummy allocation to avoid runtime errors when using pointer checking
	allocate(node_assignment(1),stat=alloc_status)
    endif 
! register data types
call MPI_Type_contiguous(3, MPI_INTEGER, mpitype_pair_assignment, ierr)
if (ierr .ne. 0) call die('failure while creating custom MPI data type')
call MPI_Type_commit(mpitype_pair_assignment, ierr)
if (ierr .ne. 0) call die('failure while creating custom MPI data type')

call MPI_Type_contiguous(11, mpitype_pair_assignment, mpitype_node_assignment, ierr)
if (ierr .ne. 0) call die('failure while creating custom MPI data type')
call MPI_Type_commit(mpitype_node_assignment, ierr)
if (ierr .ne. 0) call die('failure while creating custom MPI data type')

! distribute
call MPI_Scatter(node_assignment, 1, mpitype_node_assignment, &
    calculation_assignment, 1, mpitype_node_assignment, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('failure while sending node assignments')

! free data type
call MPI_Type_free(mpitype_node_assignment, ierr)
call MPI_Type_free(mpitype_pair_assignment, ierr)
    if (nodeid .ne. 0) then
	deallocate(node_assignment)
    endif 
end if
#endif

if (nodeid .eq. 0) then
 ! print a status report
 write(*,98) 'solute-solute', 'solute-water', 'water-water', 'Q-solute', 'Q-water', '#Shake', '#Atoms', '#Molecules'
 write(*,99) 'total',ncgp_solute,ncgp_solute,nwat,ncgp_solute,nwat,shake_molecules,natom,nmol
if (numnodes .gt. 1) then
 do i=0,numnodes-1
    write(*,100) i, 'assigned cgps', &
         node_assignment(i)%pp%end-node_assignment(i)%pp%start+1, &
         node_assignment(i)%pw%end-node_assignment(i)%pw%start+1, &
         node_assignment(i)%ww%end-node_assignment(i)%ww%start+1, &
         node_assignment(i)%qp%end-node_assignment(i)%qp%start+1, &
         node_assignment(i)%qw%end-node_assignment(i)%qw%start+1, &
         node_assignment(i)%shake%end-node_assignment(i)%shake%start+1, &
         node_assignment(i)%natom%end-node_assignment(i)%natom%start+1, &
         node_assignment(i)%nmol%end-node_assignment(i)%nmol%start+1
 end do
end if
end if

#if defined (USE_MPI)
blockcnt(:) = 1
type(:) = MPI_INTEGER
call MPI_Bcast(totnbpp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
call MPI_Bcast(totnbpw, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
call MPI_Bcast(totnbqp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
call MPI_Bcast(totnbqw, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
call MPI_Bcast(totnbww, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif

! allocate
calculation_assignment%pp%max = totnbpp/numnodes + 0.20*totnbpp
allocate(nbpp(calculation_assignment%pp%max), stat=alloc_status)
call check_alloc('solute-solute non-bond list')

calculation_assignment%pw%max = (totnbpw+3)/numnodes + 0.20*totnbpw
allocate(nbpw(calculation_assignment%pw%max), stat=alloc_status)
call check_alloc('solute-solvent non-bond list')

calculation_assignment%ww%max = (totnbww+nwat)/numnodes + 0.20*totnbww
allocate(nbww(calculation_assignment%ww%max), stat=alloc_status)
call check_alloc('solvent-solvent non-bond list')

calculation_assignment%qp%max = totnbqp/numnodes + 0.20*totnbqp
allocate(nbqp(calculation_assignment%qp%max,nstates), stat=alloc_status)
call check_alloc('Qatom - solute non-bond list')

calculation_assignment%qw%max = totnbqw/numnodes + 0.20*totnbqw
allocate(nbqw(calculation_assignment%qw%max,nstates), stat=alloc_status)
call check_alloc('Qatom - water non-bond list')



  if (use_PBC) then
	!allocate array to keep track of chargegroups
	!approximate with one half of the number of atompairs
	allocate(nbpp_cgp(calculation_assignment%pp%max / 2), stat=alloc_status)
	call check_alloc('solute-solute non-bond charge group pair list')
	allocate(nbpw_cgp(calculation_assignment%pw%max / 2), stat=alloc_status)
	call check_alloc('solute-solvent non-bond charge group pair list')
	allocate(nbqp_cgp(calculation_assignment%qp%max / 2), stat=alloc_status)
	call check_alloc('qatom-solute non-bond charge group pair list')
!we now allocate the storage arrays for MC_volume update here, so we don't have to realloc every time 
!the function is called
        allocate(old_nbpp_cgp(calculation_assignment%pp%max / 2), stat=alloc_status)
        call check_alloc('solute-solute non-bond charge group pair list')
        allocate(old_nbpw_cgp(calculation_assignment%pw%max / 2), stat=alloc_status)
        call check_alloc('solute-solvent non-bond charge group pair list')
        allocate(old_nbqp_cgp(calculation_assignment%qp%max / 2), stat=alloc_status)
        call check_alloc('qatom-solute non-bond charge group pair list')

allocate(old_nbpp(calculation_assignment%pp%max), stat=alloc_status)
call check_alloc('solute-solute non-bond list')

allocate(old_nbpw(calculation_assignment%pw%max), stat=alloc_status)
call check_alloc('solute-solvent non-bond list')

allocate(old_nbww(calculation_assignment%ww%max), stat=alloc_status)
call check_alloc('solvent-solvent non-bond list')

allocate(old_nbqp(calculation_assignment%qp%max,nstates), stat=alloc_status)
call check_alloc('Qatom - solute non-bond list')
  end if	

!Kanske deallokera nbxx_per_cgp TODO

98 format('node value ',8a13)
99 format(a10,1x,8(1x,i12))
!99 format(a4,2x,a,t18,i13,3x,i13,3x,i13,3x,i13)
100 format(i4,1x,a5,1x,8(1x,i12))
!100 format(i4,2x,a,t18,i13,3x,i13,3x,i13,3x,i13)

if (nodeid .eq. 0)  call centered_heading('End of distribution', '-')

end subroutine distribute_nonbonds
!-----------------------------------------------------------------------
subroutine lrf_taylor(LRF_loc)
! arguments
real(kind=prec)                                 :: LRF_loc
! *** local variables
integer						::	i,ic
real(kind=prec)						::	Vij, q
TYPE(qr_vec)                                            :: dr,df,tmp

! global variables used:
!  E%LRF, natom, excl, iqatom, iwhich_cgp, lrf, x, crg, d
! change to use the calculation_assignment%natom variable
! to only calculate the stuff each node should do
do i = calculation_assignment%natom%start, calculation_assignment%natom%end
! for every atom on this node:

if ( ( use_PBC .and. (iqatom(i)==0) ) .or. ( (.not. excl(i) ) .and. (iqatom(i)==0) ) ) then
! unless excluded atom or q-atom:

! find the displacement dr from the center of the charge group
ic = iwhich_cgp(i)

dr = lrf(ic)%cgp_cent - x(i)

! --- Electric potential
tmp%x = q_dotprod(dr,lrf(ic)%phi2(1))
tmp%y = q_dotprod(dr,lrf(ic)%phi2(2))
tmp%z = q_dotprod(dr,lrf(ic)%phi2(3))

Vij=lrf(ic)%phi0 + &
        q_dotprod(dr,lrf(ic)%phi1) + &
        0.5_prec*q_dotprod(dr,tmp)

LRF_loc = LRF_loc + 0.5_prec * crg(i) * Vij

! --- Electric field
tmp%x = q_dotprod(dr,lrf(ic)%phi3(1))
tmp%y = q_dotprod(dr,lrf(ic)%phi3(2))
tmp%z = q_dotprod(dr,lrf(ic)%phi3(3))

df%x = lrf(ic)%phi1%x + &
        q_dotprod(dr,lrf(ic)%phi2(1)) + &
        0.5_prec*q_dotprod(dr,tmp)

tmp%x = q_dotprod(dr,lrf(ic)%phi3(4))
tmp%y = q_dotprod(dr,lrf(ic)%phi3(5))
tmp%z = q_dotprod(dr,lrf(ic)%phi3(6))

df%y = lrf(ic)%phi1%y + &
        q_dotprod(dr,lrf(ic)%phi2(2)) + &
        0.5_prec*q_dotprod(dr,tmp)

tmp%x = q_dotprod(dr,lrf(ic)%phi3(7))
tmp%y = q_dotprod(dr,lrf(ic)%phi3(8))
tmp%z = q_dotprod(dr,lrf(ic)%phi3(9))

df%z = lrf(ic)%phi1%z + &
        q_dotprod(dr,lrf(ic)%phi2(3)) + &
        0.5_prec*q_dotprod(dr,tmp)

! update d
d(i) = d(i) - df * crg(i)
end if
end do

end subroutine lrf_taylor

!-----------------------------------------------------------------------
! function to gather lrf information from nodes and set it the values
! gathered there
#ifdef USE_MPI
integer function lrf_add(lrf1,lrf2,length,MPI_Datatype)
type(LRF_TYPE),INTENT (IN)      :: lrf1(length)
type(LRF_TYPE)                  :: lrf2(length)
integer                         :: length,MPI_Datatype

lrf2(1:length)%cgp_cent           = lrf1(1:length)%cgp_cent
lrf2(1:length)%phi0               = lrf1(1:length)%phi0+lrf2(1:length)%phi0
lrf2(1:length)%phi1               = lrf1(1:length)%phi1   +lrf2(1:length)%phi1    
lrf2(1:length)%phi2(1)            = lrf1(1:length)%phi2(1)+lrf2(1:length)%phi2(1)
lrf2(1:length)%phi2(2)            = lrf1(1:length)%phi2(2)+lrf2(1:length)%phi2(2)
lrf2(1:length)%phi2(3)            = lrf1(1:length)%phi2(3)+lrf2(1:length)%phi2(3)
lrf2(1:length)%phi3(1)            = lrf1(1:length)%phi3(1)+lrf2(1:length)%phi3(1)
lrf2(1:length)%phi3(2)            = lrf1(1:length)%phi3(2)+lrf2(1:length)%phi3(2)
lrf2(1:length)%phi3(3)            = lrf1(1:length)%phi3(3)+lrf2(1:length)%phi3(3)
lrf2(1:length)%phi3(4)            = lrf1(1:length)%phi3(4)+lrf2(1:length)%phi3(4)
lrf2(1:length)%phi3(5)            = lrf1(1:length)%phi3(5)+lrf2(1:length)%phi3(5)
lrf2(1:length)%phi3(6)            = lrf1(1:length)%phi3(6)+lrf2(1:length)%phi3(6)
lrf2(1:length)%phi3(7)            = lrf1(1:length)%phi3(7)+lrf2(1:length)%phi3(7)
lrf2(1:length)%phi3(8)            = lrf1(1:length)%phi3(8)+lrf2(1:length)%phi3(8)
lrf2(1:length)%phi3(9)            = lrf1(1:length)%phi3(9)+lrf2(1:length)%phi3(9)

lrf_add = -1

end function lrf_add
!-----------------------------------------------------------------------
integer function lrf_cgp_rep(lrf1,lrf2,length,MPI_Datatype)
type(LRF_TYPE),INTENT (IN)      :: lrf1(length)
type(LRF_TYPE)                  :: lrf2(length)
integer                         :: length,MPI_Datatype
lrf2(1:length)%cgp_cent           = lrf1(1:length)%cgp_cent
lrf_cgp_rep = -1

end function lrf_cgp_rep
#endif
!-----------------------------------------------------------------------
! this one gets the lrf info from all nodes and puts it on master after a lrf
! update has been completed
! only used with MPI
#ifdef USE_MPI
subroutine lrf_gather
call MPI_Allreduce(MPI_IN_PLACE,LRF,ncgp,mpitype_batch_lrf,mpi_lrf_add,&
MPI_COMM_WORLD,ierr)
if (ierr .ne. 0) call die('lrf_gather reduce')
!if (nodeid .eq.0 ) LRF = lrf_reduce
end subroutine lrf_gather
#endif

!-----------------------------------------------------------------------
!DEC$ ATTRIBUTES FORCEINLINE :: lrf_update
subroutine lrf_update(group1,group2)
! --- input variables
integer				:: group1,group2

! --- local variables
integer				:: i,ia,i_sw,j,k,l,m,n
real(kind=prec)			:: r2,field0,field1,field2
TYPE(qr_vec)                    :: shift,dr,tmp(9),uvec
! get switch atoms for PBC
if (use_PBC) then
i_sw = cgp(group1)%iswitch
shift = x(i_sw) - lrf(group2)%cgp_cent
shift = boxlength*q_nint(shift*inv_boxl)
else
shift = zero
end if
! unit vector for later assignments
  uvec = one


iloop:        do ia = cgp(group1)%first, cgp(group1)%last
  ! skip if q-atom
  i = cgpatom(ia)
  if ( iqatom(i)/=0 ) cycle iloop

  dr = x(i) - lrf(group2)%cgp_cent - shift
  r2 = qvec_square(dr)

  field0=crg(i)/(r2*q_sqrt(r2))
  field1=3.0_prec*field0/r2
  field2=-field1/r2
!$omp critical
  lrf(group2)%phi0    = lrf(group2)%phi0+field0*r2
  lrf(group2)%phi1    = lrf(group2)%phi1 - dr * field0
  tmp(1) = dr * field1
  lrf(group2)%phi2(1) = lrf(group2)%phi2(1) + dr * tmp(1)%x 
  lrf(group2)%phi2(2) = lrf(group2)%phi2(2) + dr * tmp(1)%y 
  lrf(group2)%phi2(3) = lrf(group2)%phi2(3) + dr * tmp(1)%z 

  lrf(group2)%phi2(1)%x = lrf(group2)%phi2(1)%x - field0
  lrf(group2)%phi2(2)%y = lrf(group2)%phi2(2)%y - field0
  lrf(group2)%phi2(3)%z = lrf(group2)%phi2(3)%z - field0
!$omp end critical
!$omp critical

  tmp(1) = uvec * dr%x
  tmp(2) = uvec * dr%y
  tmp(3) = uvec * dr%z

  do n = 1,3

  do j = 1,3
  k = j+3
  l = k+3
  m = j + (n-1)*3
  tmp(k) = tmp(n) * tmp(j)
  tmp(l) = tmp(k) * dr
  tmp(l) = tmp(l) * 5.0_prec
  tmp(l) = tmp(l) * field2
  lrf(group2)%phi3(m) = lrf(group2)%phi3(m) + tmp(l)
  end do

  end do

  tmp(1) = dr * r2
  lrf(group2)%phi3(1)%x = lrf(group2)%phi3(1)%x - field2*(3.0_prec*tmp(1)%x)
  lrf(group2)%phi3(1)%y = lrf(group2)%phi3(1)%y - field2*(tmp(1)%y)
  lrf(group2)%phi3(1)%z = lrf(group2)%phi3(1)%z - field2*(tmp(1)%z)
  lrf(group2)%phi3(2)%x = lrf(group2)%phi3(2)%x - field2*(tmp(1)%y)
  lrf(group2)%phi3(2)%y = lrf(group2)%phi3(2)%y - field2*(tmp(1)%x)

  lrf(group2)%phi3(3)%x = lrf(group2)%phi3(3)%x - field2*(tmp(1)%z)

  lrf(group2)%phi3(3)%z = lrf(group2)%phi3(3)%z - field2*(tmp(1)%x)
  lrf(group2)%phi3(4)%x = lrf(group2)%phi3(4)%x - field2*(tmp(1)%y)
  lrf(group2)%phi3(4)%y = lrf(group2)%phi3(4)%y - field2*(tmp(1)%x)

  lrf(group2)%phi3(5)%x = lrf(group2)%phi3(5)%x - field2*(tmp(1)%x)
  lrf(group2)%phi3(5)%y = lrf(group2)%phi3(5)%y - field2*(3.0_prec*tmp(1)%y)
  lrf(group2)%phi3(5)%z = lrf(group2)%phi3(5)%z - field2*(tmp(1)%z)

  lrf(group2)%phi3(6)%y = lrf(group2)%phi3(6)%y - field2*(tmp(1)%z)
  lrf(group2)%phi3(6)%z = lrf(group2)%phi3(6)%z - field2*(tmp(1)%y)
  lrf(group2)%phi3(7)%x = lrf(group2)%phi3(7)%x - field2*(tmp(1)%z)

  lrf(group2)%phi3(7)%z = lrf(group2)%phi3(7)%z - field2*(tmp(1)%x)

  lrf(group2)%phi3(8)%y = lrf(group2)%phi3(8)%y - field2*(tmp(1)%z)
  lrf(group2)%phi3(8)%z = lrf(group2)%phi3(8)%z - field2*(tmp(1)%y)
  lrf(group2)%phi3(9)%x = lrf(group2)%phi3(9)%x - field2*(tmp(1)%x)
  lrf(group2)%phi3(9)%y = lrf(group2)%phi3(9)%y - field2*(tmp(1)%y)
  lrf(group2)%phi3(9)%z = lrf(group2)%phi3(9)%z - field2*(3.0_prec*tmp(1)%z)


!$omp end critical
end do iloop

end subroutine lrf_update

!-----------------------------------------------------------------------

subroutine make_nbqqlist
!locals
integer						::	is

call make_qconn
nbqq_max = nbqq_count()
allocate(nbqq(nbqq_max, nstates),nbqqp(nbqq_max, nstates),stat=alloc_status)
call check_alloc('Qatom-Qatom non-bond list')
do is =1, nstates
        write (*,200) nbqq_pair(is),is
end do
write (*,*)

200 format ('No. of Rcq indep. nb pairs involving q-atoms = ',i5, &
' in state :',i3)
!qqlist now done after precomputation
end subroutine make_nbqqlist

!-----------------------------------------------------------------------

subroutine make_pair_lists(Rq,Rcq2,RcLRF2,Rcpp2,Rcpw2,Rcww2)
! arguments
real(kind=prec)                 :: Rq,Rcq2,Rcpp2,Rcpw2,Rcww2,RcLRF2
#if defined (PROFILING)
real(kind=prec)                                         :: start_loop_time
start_loop_time = rtime()
#endif
! start with populating the grid
! each node populates its assigned ncgp_solute
! and nwat and then intercommunication does the rest

! we make this whole thing now omp parallel
! start opening region here so we don't have to pay the penalty
! to do this in every function

!$omp parallel default(none) shared(use_PBC,use_LRF,iuse_switch_atom) 

#ifdef USE_GRID
call populate_grids
#endif
if( use_PBC ) then
if (.not. use_LRF)then
        !cutoff
   if(iuse_switch_atom == 1) then
        call nbpplist_box(Rcpp2)
        call nbpwlist_box(Rcpw2)
        call nbqplist_box(Rcq2,Rq)
    else
        call nbpplis2_box(Rcpp2)
        call nbpwlis2_box(Rcpw2)
        call nbqplis2_box(Rcq2,Rq)
    end if
    call nbwwlist_box(Rcww2)
else
    call cgp_centers
                      if ( iuse_switch_atom == 1 ) then 
                        call nbpplist_box_lrf(Rcpp2,RcLRF2)
                              call nbpwlist_box_lrf(Rcpw2,RcLRF2)
                              call nbqplist_box(Rcq2,Rq)
                      else
                        call nbpplis2_box_lrf(Rcpp2,RcLRF2)
                              call nbpwlis2_box_lrf(Rcpw2,RcLRF2)
                              call nbqplis2_box(Rcq2,Rq)
    endif
    call nbwwlist_box_lrf(Rcww2,RcLRF2)
endif
        call nbqwlist_box(Rcq2,Rq)
else !spherical case
if(.not. use_LRF) then
        ! cutoff
        if( iuse_switch_atom .eq. 1 ) then
                call nbpplist(Rcpp2)
                call nbpwlist(Rcpw2)
                call nbqplist(Rcq2)
        else
                call nbpplis2(Rcpp2)
                call nbpwlis2(Rcpw2)
                call nbqplis2(Rcq2)
        end if
        call nbwwlist(Rcww2)
else 
        ! cutoff with lrf
        call cgp_centers ! *** måste anropas av alla noder (nollställer lrf)
        if( iuse_switch_atom .eq. 1 ) then
                call nbpplist_lrf(Rcpp2,RcLRF2)
                call nbpwlist_lrf(Rcpw2,RcLRF2)
                call nbqplist(Rcq2)
        else
                call nbpplis2_lrf(Rcpp2,RcLRF2)
                call nbpwlis2_lrf(Rcpw2,RcLRF2)
                call nbqplis2(Rcq2)
        end if
        call nbwwlist_lrf(Rcww2,RcLRF2)
end if

call nbqwlist(Rcq2)
end if

!$omp end parallel

! we are gathering the LRF info here
#ifdef USE_MPI
if (use_LRF) call lrf_gather
#endif
#if defined (PROFILING)
profile(1)%time = profile(1)%time + rtime() - start_loop_time
#endif

end subroutine make_pair_lists
!-----------------------------------------------------------------------


subroutine nbpp_count(npp, nppcgp,rcut)
! arguments
integer						:: npp
integer						:: nppcgp(:)
real(kind=prec)                                 :: rcut
! local variables
integer						:: i,j,ig,jg,ia,ja,nl
TYPE(qr_vec)					:: shift
real(kind=prec)						:: rcut2,r2
integer						:: LJ_code
! This routine counts non-bonded solute-solute atom pairs 
! excluding any Q-atoms.

! uses the global variables:
!  Rcpp, ncgp, cgp, excl, x, cgpatom, iqatom, ljcod, crg, 
!  iaclib, max_nbr_range, listex, nexlong, listexlong


npp = 0
rcut2 = rcut**2

igloop: do ig = 1, ncgp_solute
nppcgp(ig) = 0

ia = cgp(ig)%iswitch
 
! skip if excluded group
if ( .not. use_PBC .and. excl(ia) ) cycle igloop


jgloop:	do jg = 1, ncgp_solute
  ja = cgp(jg)%iswitch
 
   ! skip if excluded group
  if ( .not. use_PBC .and. excl(ja) ) cycle jgloop

  ! count each charge group pair once only
  if ( ((ig .gt. jg) .and. (mod(ig+jg,2) .eq. 0)) .or. &
           ((ig .lt. jg) .and. (mod(ig+jg,2) .eq. 1)) ) &
           cycle jgloop


  !******PWadded if-statement 2001-10-01

  if( .not. use_PBC ) then
  r2 = q_dist4(x(ia),x(ja))
  else
	shift = x(ia) - x(ja)
        r2 = q_dist4(shift,boxlength*q_nint(shift*inv_boxl))
  end if


  ! skip if outside cutoff
  if ( r2 .gt. rcut2 ) cycle jgloop

ialoop: do ia = cgp(ig)%first, cgp(ig)%last
        i = cgpatom(ia)

        ! skip if q-atom
        if ( iqatom(i)/=0 ) cycle ialoop

jaloop:	do ja = cgp(jg)%first, cgp(jg)%last
          j = cgpatom(ja)

          ! skip if q-atom
          if ( iqatom(j)/=0 ) cycle jaloop

          if ( ig .eq. jg .and. i .ge. j ) cycle jaloop

          ! skip if all interactions zero
          LJ_code = ljcod(iac(i),iac(j))
          if((crg(i) * crg(j) == zero) &
                .and. &
                (iaclib(iac(i))%avdw(LJ_code)*iaclib(iac(j))%avdw(LJ_code) == zero) &
                .and. &
                (iaclib(iac(i))%bvdw(LJ_code)*iaclib(iac(j))%bvdw(LJ_code) == zero)) &
                cycle jaloop

          ! check bonded exclusions and 1-4 nbors
          if ( abs(j-i) .le. max_nbr_range ) then
                if ( i .lt. j ) then
                  if ( listex(j-i,i) ) cycle jaloop
                else
                  if ( listex(i-j,j) ) cycle jaloop
                end if
          else
                do nl = 1, nexlong
                  if ( (listexlong(1,nl) .eq. i .and. &
                        listexlong(2,nl) .eq. j      ) .or. &
                        (listexlong(1,nl) .eq. j .and. &
                        listexlong(2,nl) .eq. i      ) ) cycle jaloop
                end do
          end if

          ! passed all tests -- count the pair
          npp = npp + 1
          nppcgp(ig) = nppcgp(ig) + 1

        end do jaloop
  end do ialoop
end do jgloop
end do igloop

end subroutine nbpp_count

!-----------------------------------------------------------------------


subroutine nbpplis2(Rcut)
! args
real(kind=prec)					:: Rcut
! local variables
integer						:: i,j,ig,jg,ia,ja,nl,inside,jgr
real(kind=prec)						:: rcut2,r2

integer						::iagrid, igrid, jgrid, kgrid, gridnum

! for spherical boundary 
!	This routine makes a list of non-bonded solute-solute atom pairs 
!	excluding any Q-atoms.
! now rewritten to take into account that we already know all the stuff from the
! precomputation

#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)                                         :: start_loop_time
start_loop_time = rtime()
#endif
!$omp single 
nbpp_pair = 0
!$omp end single 
!$omp barrier
rcut2 = Rcut

! new stuff because we are now having a grid 
! check first if the atoms are in one the adjecent grids
! if not, cycle them
! grid information taken from new costum structure, updated every nb steps
! after getting the first group, we cycle through the grids to find those that interact
! and only search further in them

! interaction matrix is saved in grid_pp_int array 

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%pp%end - calculation_assignment%pp%start + 1)/threads_num
remainder = MOD(calculation_assignment%pp%end - calculation_assignment%pp%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%pp%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

igloop: do ig = mp_start,mp_end
#else
igloop:  do ig = calculation_assignment%pp%start, calculation_assignment%pp%end
#endif

ia = cgp(ig)%iswitch
if ( excl(ia) ) cycle igloop
#ifdef USE_GRID
iagrid = pp_igrid(ig)
!start at the first grid
gridnum = 0
igridloop:  do igrid = 1, pp_ndim
jgridloop:   do jgrid = 1, pp_ndim
kgridloop:    do kgrid = 1, pp_ndim
		gridnum = gridnum + 1
		if (.not. (grid_pp_int(iagrid,igrid,jgrid,kgrid))) cycle kgridloop
! else, we have interacting grids
! now search for the charge groups in there
! and cycle over all of them
! save the real number in jgr
jgloop: do jg = 1, grid_pp_ngrp(gridnum)
                jgr = grid_pp_grp(gridnum,jg)
#else
jgloop: do jgr = 1, ncgp_solute
#endif
  ! count each charge group pair once only
                if ( ((ig .gt. jgr) .and. (mod(ig+jgr,2) .eq. 0)) .or. &
                        ((ig .lt. jgr) .and. (mod(ig+jgr,2) .eq. 1)) ) &
                        cycle jgloop
                ja = cgp(jgr)%iswitch
                if ( excl(ja) ) cycle jgloop

!	      --- outside cutoff ? ---
                inside = 0
                ia = cgp(ig)%first
                do while ((ia .le. cgp(ig)%last) .and. (inside .eq. 0))
                        i = cgpatom(ia)
                        ja = cgp(jgr)%first
                        do while ((ja .le. cgp(jgr)%last) .and. (inside .eq. 0))
                                j = cgpatom(ja)
                                r2 = q_dist4(x(i),x(j))
                                if ( r2 .le. rcut2 ) then
! one atom pair is within cutoff: set inside
                                        inside = 1
                                end if
                                ja = ja + 1
                        end do
                        ia = ia + 1
                end do
                if (inside .eq. 0) cycle jgloop
ialoop:         do ia = cgp(ig)%first, cgp(ig)%last
                i = cgpatom(ia)
!	         --- q-atom ? ---
                if ( iqatom(i).ne.0 ) cycle ialoop

jaloop:                 do ja = cgp(jgr)%first, cgp(jgr)%last
                                j = cgpatom(ja)
!                --- q-atom ? ---
                                if ( iqatom(j).ne.0 ) cycle jaloop
! make sure each pair is only counted once
                                if ( ig .eq. jgr .and. i .ge. j ) cycle jaloop
! do they interact?
                                if(pp_map(i-pp_low,j).eq.0) cycle jaloop
                                if(.not.pp_precomp(pp_map(i-pp_low,j))%set) cycle jaloop
! if out of space then make more space
!$omp critical
                                if (nbpp_pair .eq. calculation_assignment%pp%max) call reallocate_nonbondlist_pp
                                nbpp_pair = nbpp_pair + 1
                                nbpp(nbpp_pair)%i = i
                                nbpp(nbpp_pair)%j = j
                                nbpp(nbpp_pair)%vdWA = pp_precomp(pp_map(i-pp_low,j))%vdWA
                                nbpp(nbpp_pair)%vdWB = pp_precomp(pp_map(i-pp_low,j))%vdWB
                                nbpp(nbpp_pair)%elec = pp_precomp(pp_map(i-pp_low,j))%elec
!$omp end critical

                        end do jaloop
                end do ialoop
        end do jgloop
#ifdef USE_GRID
	end do kgridloop
	end do jgridloop
	end do igridloop
#endif
end do igloop
#if defined (PROFILING) 
profile(3)%time = profile(3)%time + rtime() - start_loop_time
#endif 

end subroutine nbpplis2


!--------------------------------------------------------------------

subroutine nbpplis2_box(Rcut)
! args
real(kind=prec)					:: Rcut
! local variables
integer						:: i,j,ig,jg,ia,ja,nl,inside,jgr
real(kind=prec)					:: rcut2,r2
TYPE(qr_vec)					:: shift

integer                                         ::iagrid, igrid, jgrid, kgrid, gridnum
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)                                         :: start_loop_time
start_loop_time = rtime()
#endif
 
  ! for periodic boundary conditions
  !	This routine makes a list of non-bonded solute-solute atom pairs 
  !	excluding any Q-atoms.
!$omp single
  nbpp_pair = 0
  nbpp_cgp_pair = 0
!$omp end single
!$omp barrier
  rcut2 = Rcut
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%pp%end - calculation_assignment%pp%start + 1)/threads_num
remainder = MOD(calculation_assignment%pp%end - calculation_assignment%pp%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%pp%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

igloop: do ig = mp_start,mp_end
#else
igloop: do ig = calculation_assignment%pp%start, calculation_assignment%pp%end
#endif
#ifdef USE_GRID
        iagrid = pp_igrid(ig)
        gridnum = 0
igridloop: do igrid = 1, pp_ndim
jgridloop:  do jgrid = 1, pp_ndim
kgridloop:   do kgrid = 1, pp_ndim
		gridnum = gridnum + 1
		if (.not. (grid_pp_int(iagrid,igrid,jgrid,kgrid))) cycle kgridloop
! else, we have interacting grids
! now search for the charge groups in there
! and cycle over all of them
! save the real number in jgr
jgloop:         do jg = 1, grid_pp_ngrp(gridnum)
                        jgr = grid_pp_grp(gridnum,jg)
#else
jgloop:         do jgr = 1, ncgp_solute
#endif
! count each charge group pair once only
                        if ( ((ig .gt. jgr) .and. (mod(ig+jgr,2) .eq. 0)) .or. &
                                ((ig .lt. jgr) .and. (mod(ig+jgr,2) .eq. 1)) ) &
                                cycle jgloop

!      --- outside cutoff ? ---
                        inside = 0
                        ia = cgp(ig)%first
                        do while ((ia .le. cgp(ig)%last) .and. (inside .eq. 0))
                                i = cgpatom(ia)
                                ja = cgp(jgr)%first
                                do while ((ja .le. cgp(jgr)%last) .and. (inside .eq. 0))
                                        j = cgpatom(ja)
					shift = x(i) -x(j)
                                        r2 = q_dist4(shift,boxlength*q_nint(shift*inv_boxl))
                                        if ( r2 .le. rcut2 ) then
! one atom pair is within cutoff: set inside
                                                inside = 1
!$omp critical
                                                if (nbpp_cgp_pair .eq. size(nbpp_cgp, 1) )  call reallocate_nbpp_cgp
                                                nbpp_cgp_pair = nbpp_cgp_pair + 1
                                                nbpp_cgp(nbpp_cgp_pair)%i = i
                                                nbpp_cgp(nbpp_cgp_pair)%j = j
!$omp end critical
                                        end if
                                        ja = ja + 1
                                end do
                                ia = ia + 1
                        end do
                        if (inside .eq. 0) cycle jgloop

ialoop:                 do ia = cgp(ig)%first, cgp(ig)%last
                                i = cgpatom(ia)
!	         --- q-atom ? ---
                                if ( iqatom(i).ne.0 ) cycle ialoop
jaloop:                         do ja = cgp(jgr)%first, cgp(jgr)%last
                                        j = cgpatom(ja)
!                --- q-atom ? ---
                                        if ( iqatom(j).ne.0 ) cycle jaloop
! make sure each pair is only counted once
                                        if ( ig .eq. jgr .and. i .ge. j ) cycle jaloop
! do they interact?
                                        if(pp_map(i-pp_low,j).eq.0) cycle jaloop
                                        if(.not.pp_precomp(pp_map(i-pp_low,j))%set) cycle jaloop
! if out of space then make more space
!$omp critical
                                        if (nbpp_pair .eq. calculation_assignment%pp%max) call reallocate_nonbondlist_pp
                                        nbpp_pair = nbpp_pair + 1
                                        nbpp(nbpp_pair)%i = i
                                        nbpp(nbpp_pair)%j = j
                                        nbpp(nbpp_pair)%vdWA = pp_precomp(pp_map(i-pp_low,j))%vdWA
                                        nbpp(nbpp_pair)%vdWB = pp_precomp(pp_map(i-pp_low,j))%vdWB
                                        nbpp(nbpp_pair)%elec = pp_precomp(pp_map(i-pp_low,j))%elec
                                        nbpp(nbpp_pair)%cgp_pair = nbpp_cgp_pair
!$omp end critical
                                end do jaloop
                        end do ialoop
                end do jgloop
#ifdef USE_GRID
	end do kgridloop
	end do jgridloop
	end do igridloop
#endif
end do igloop
#if defined (PROFILING) 
profile(3)%time = profile(3)%time + rtime() - start_loop_time
#endif 

end subroutine nbpplis2_box
!-----------------------------------------------------------------------------
subroutine nbpplis2_box_lrf(Rcut,RLRF)
! args
real(kind=prec)						:: Rcut,RLRF
! local variables
integer						:: i,j,ig,jg,ia,ja,nl,inside,jgr
TYPE(qr_vec)					:: shift
real(kind=prec)						:: rcut2,r2

real(kind=prec)						::RcLRF2
integer						::inside_LRF, is3
integer						::iagrid, igrid, jgrid, kgrid, gridnum
! for periodic boundary conditions
!	This routine makes a list of non-bonded solute-solute atom pairs 
!	excluding any Q-atoms.
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)                                         :: start_loop_time
start_loop_time = rtime()
#endif
!$omp single
  nbpp_pair = 0
  nbpp_cgp_pair = 0
!$omp end single
!$omp barrier
  rcut2 = Rcut
  RcLRF2 = RLRF
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%pp%end - calculation_assignment%pp%start + 1)/threads_num
remainder = MOD(calculation_assignment%pp%end - calculation_assignment%pp%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%pp%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

igloop: do ig = mp_start,mp_end
#else
igloop: do ig = calculation_assignment%pp%start, calculation_assignment%pp%end
#endif
#ifdef USE_GRID
         iagrid = pp_igrid(ig)
	  gridnum = 0
igridloop:  do igrid = 1, pp_ndim
jgridloop:   do jgrid = 1, pp_ndim
kgridloop:    do kgrid = 1, pp_ndim
		gridnum = gridnum + 1
		if (.not. (grid_pp_int(iagrid,igrid,jgrid,kgrid))) then
ggloop:		 do jg = 1, grid_pp_ngrp(gridnum)
		  jgr = grid_pp_grp(gridnum,jg)
		   if ( ((ig .gt. jgr) .and. (mod(ig+jgr,2) .eq. 0)) .or. &
                   ((ig .lt. jgr) .and. (mod(ig+jgr,2) .eq. 1)) .or.&
                   (ig.eq.jgr)) &
                   cycle ggloop
		    call lrf_update(ig,jgr)
		    call lrf_update(jgr,ig)
		   end do ggloop
		 else 
! we have real interacting grids , do the rest of the calculations
jgloop: do jg = 1, grid_pp_ngrp(gridnum)
                jgr = grid_pp_grp(gridnum,jg)
#else
jgloop: do jgr = 1, ncgp_solute
#endif
! count each charge group pair once only
                if ( ((ig .gt. jgr) .and. (mod(ig+jgr,2) .eq. 0)) .or. &
                        ((ig .lt. jgr) .and. (mod(ig+jgr,2) .eq. 1)) ) &
                        cycle jgloop
!      --- outside cutoff ? ---
                inside = 0
                inside_LRF = 0
                ia = cgp(ig)%first
                do while ((ia .le. cgp(ig)%last) .and. (inside .eq. 0))
                        i = cgpatom(ia)
                        ja = cgp(jgr)%first
                        do while ((ja .le. cgp(jgr)%last) .and. (inside .eq. 0))
                                j = cgpatom(ja)
                                shift = x(i) - x(j)
                                r2 = q_dist4(shift,boxlength*q_nint(shift*inv_boxl))
                                if ( r2 .le. rcut2 ) then
! one atom pair is within cutoff: set inside
                                        inside = 1
!$omp critical
                                        if (nbpp_cgp_pair .eq. size(nbpp_cgp, 1) )  call reallocate_nbpp_cgp
                                        nbpp_cgp_pair = nbpp_cgp_pair + 1
                                        nbpp_cgp(nbpp_cgp_pair)%i = i
                                        nbpp_cgp(nbpp_cgp_pair)%j = j
!$omp end critical
                                elseif ((r2 .le. RcLRF2) .or. (RLRF .eq. -one)) then
                                        inside_LRF = 1
                                end if
                                ja = ja + 1
                        end do ! ia .le. cgp(ig)%last
                        ia = ia + 1
                end do ! ja .le. cgp(jgr)%last
		
                if (inside .eq. 1) then
ialoop:                 do ia = cgp(ig)%first, cgp(ig)%last
                                i = cgpatom(ia)
!	         --- q-atom ? ---
                                if ( iqatom(i).ne.0 ) cycle ialoop
jaloop:                         do ja = cgp(jgr)%first, cgp(jgr)%last
                                        j = cgpatom(ja)
!                --- q-atom ? ---
                                        if ( iqatom(j).ne.0 ) cycle jaloop
! make sure each pair is only counted once
                                        if ( ig .eq. jgr .and. i .ge. j ) cycle jaloop
! do they interact?
                                        if(pp_map(i-pp_low,j).eq.0) cycle jaloop
                                        if(.not.pp_precomp(pp_map(i-pp_low,j))%set) cycle jaloop
! if out of space then make more space
!$omp critical
                                        if (nbpp_pair .eq. calculation_assignment%pp%max) call reallocate_nonbondlist_pp
                                        nbpp_pair = nbpp_pair + 1
                                        nbpp(nbpp_pair)%i = i
                                        nbpp(nbpp_pair)%j = j
                                        nbpp(nbpp_pair)%vdWA = pp_precomp(pp_map(i-pp_low,j))%vdWA
                                        nbpp(nbpp_pair)%vdWB = pp_precomp(pp_map(i-pp_low,j))%vdWB
                                        nbpp(nbpp_pair)%elec = pp_precomp(pp_map(i-pp_low,j))%elec
                                        nbpp(nbpp_pair)%cgp_pair = nbpp_cgp_pair
!$omp end critical
                                end do jaloop
                        end do ialoop
                elseif((inside_LRF .eq. 1) .and. (inside .eq. 0)) then
! outside pp-cutoff but inside LRF cut-off use LRF
!ig : jgr calculation
                        call lrf_update(ig,jgr)
!jgr : ig calculations
                        call lrf_update(jgr,ig)
                end if ! outside cutoff
        end do jgloop
#ifdef USE_GRID
	end if ! interacting boxes
	end do kgridloop
	end do jgridloop
	end do igridloop
#endif
end do igloop
#if defined (PROFILING) 
profile(3)%time = profile(3)%time + rtime() - start_loop_time
#endif 

end subroutine nbpplis2_box_lrf
!-----------------------------------------------------------------------
subroutine nbpplis2_lrf(Rcut,RLRF)
! args
real(kind=prec)					:: Rcut,RLRF
! local variables
integer						:: i,j,ig,jg,ia,ja,nl,is,jgr
logical						::	inside
real(kind=prec)						:: rcut2,r2
real(kind=prec)						::	RcLRF2

integer							::iagrid, igrid, jgrid, kgrid, gridnum
!	This routine makes a list of non-bonded solute-solute atom pairs 
!	excluding any Q-atoms.
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)                                         :: start_loop_time
start_loop_time = rtime()
#endif
!$omp single
nbpp_pair = 0
!$omp end single
!$omp barrier
rcut2 = Rcut
RcLRF2 = RLRF

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%pp%end - calculation_assignment%pp%start + 1)/threads_num
remainder = MOD(calculation_assignment%pp%end - calculation_assignment%pp%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%pp%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

igloop: do ig = mp_start,mp_end
#else
igloop: do ig = calculation_assignment%pp%start, calculation_assignment%pp%end
#endif

        ! skip if excluded group
        is = cgp(ig)%iswitch
        if ( excl(is) ) cycle igloop
#ifdef USE_GRID
        gridnum = 0
	iagrid = pp_igrid(ig)
igridloop:  do igrid = 1, pp_ndim
jgridloop:   do jgrid = 1, pp_ndim
kgridloop:    do kgrid = 1, pp_ndim
		gridnum = gridnum + 1
                if (.not. (grid_pp_int(iagrid,igrid,jgrid,kgrid))) then
ggloop:          do jg = 1, grid_pp_ngrp(gridnum)
                  jgr = grid_pp_grp(gridnum,jg)
		  ja = cgp(jgr)%iswitch
                  if ( excl(ja) ) cycle ggloop
                   if ( ((ig .gt. jgr) .and. (mod(ig+jgr,2) .eq. 0)) .or. &
                   ((ig .lt. jgr) .and. (mod(ig+jgr,2) .eq. 1)) .or. &
                   (ig.eq.jgr)) &
                   cycle ggloop
                    call lrf_update(ig,jgr)
                    call lrf_update(jgr,ig)
                   end do ggloop
                 else
! we have real interacting grids , do the rest of the calculations
jgloop: do jg = 1, grid_pp_ngrp(gridnum)
                jgr = grid_pp_grp(gridnum,jg)
#else
jgloop: do jgr = 1, ncgp_solute
#endif
! count each charge group pair once only
                if( ((ig .gt. jgr) .and. (mod(ig+jgr,2) .eq. 0)) .or. &
                        ((ig .lt. jgr) .and. (mod(ig+jgr,2) .eq. 1)) ) &
                        cycle jgloop

!	      --- excluded group ? ---
                ja = cgp(jgr)%iswitch
                if ( excl(ja) ) cycle jgloop

!	      --- outside cutoff ? ---
                inside = .false.
                ia = cgp(ig)%first
pairloop:	do while( (ia .le. cgp(ig)%last) .and. (inside.eqv..false.) )
                        i = cgpatom(ia)
			ja = cgp(jgr)%first
                        do while (( ja.le.cgp(jgr)%last) .and. (inside.eqv..false.))
                                j = cgpatom(ja)
				r2 = q_dist4(x(i),x(j))

                                if ( r2 <= rcut2 ) then
                                        ! one atom pair is within cutoff: set inside
                                        inside = .true.
                                end if
                        ja = ja + 1
                        end do
                        ia = ia + 1
                end do pairloop
!	      --- inside cutoff ? ---
                if (inside) then
ialoop:                 do ia = cgp(ig)%first, cgp(ig)%last
                                i = cgpatom(ia)

                                !	             --- q-atom ? ---
                                if ( iqatom(i).ne.0 ) cycle ialoop
jaloop:				do ja = cgp(jgr)%first, cgp(jgr)%last
                                        j = cgpatom(ja)
!                --- q-atom ? ---
                                        if ( iqatom(j).ne.0 ) cycle jaloop
! make sure each pair is only counted once
                                        if ( ig .eq. jgr .and. i .ge. j ) cycle jaloop
! do they interact?
                                        if(pp_map(i-pp_low,j).eq.0) cycle jaloop
                                        if(.not.pp_precomp(pp_map(i-pp_low,j))%set) cycle jaloop

                                        ! if out of space then make more space
!$omp critical
                                        if (nbpp_pair == calculation_assignment%pp%max) call reallocate_nonbondlist_pp

                                        nbpp_pair = nbpp_pair + 1
                                        nbpp(nbpp_pair)%i = i
                                        nbpp(nbpp_pair)%j = j
                                        nbpp(nbpp_pair)%vdWA = pp_precomp(pp_map(i-pp_low,j))%vdWA
                                        nbpp(nbpp_pair)%vdWB = pp_precomp(pp_map(i-pp_low,j))%vdWB
                                        nbpp(nbpp_pair)%elec = pp_precomp(pp_map(i-pp_low,j))%elec
!$omp end critical
                                end do jaloop
                        end do ialoop
                elseif(r2 .le. RcLRF2) then
                        ! outside pp-cutoff but inside LRF cut-off: use LRF
			call lrf_update(ig,jgr)
			call lrf_update(jgr,ig)
                end if
        end do jgloop
#ifdef USE_GRID
	end if ! interacting boxes
	end do kgridloop
	end do jgridloop
	end do igridloop
#endif
end do igloop
#if defined (PROFILING) 
profile(3)%time = profile(3)%time + rtime() - start_loop_time
#endif 

end subroutine nbpplis2_lrf

!-----------------------------------------------------------------------

subroutine nbpplist(Rcut)
! args
real(kind=prec)					:: Rcut
! local variables
integer						:: i,j,ig,jg,ia,ja,nl,jgr
real(kind=prec)						:: rcut2,r2
integer						:: LJ_code
integer                                                 ::iagrid, igrid, jgrid, kgrid, gridnum

! For use with spherical boundary   
!	This routine makes a list of non-bonded solute-solute atom pairs 
!	excluding any Q-atoms.

! uses the global variables:
!  Rcpp, ncgp, cgp, excl, x, cgpatom, iqatom, ljcod, crg, iaclib, max_nbr_range, listex
!  nexlong, listexlong, calculation_assignment%pp%max, alloc_status, list14, n14long, list14long

! reset #pairs
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)                                         :: start_loop_time
start_loop_time = rtime()
#endif
!$omp single
nbpp_pair = 0
!$omp end single
!$omp barrier
rcut2 = Rcut
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%pp%end - calculation_assignment%pp%start + 1)/threads_num
remainder = MOD(calculation_assignment%pp%end - calculation_assignment%pp%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%pp%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

igloop: do ig = mp_start,mp_end
#else
igloop: do ig = calculation_assignment%pp%start, calculation_assignment%pp%end
#endif
! for every assigned charge group:

! skip if excluded group
        ia = cgp(ig)%iswitch
        if ( excl(ia) ) cycle igloop
#ifdef USE_GRID
iagrid = pp_igrid(ig)
gridnum = 0

igridloop:	do igrid = 1, pp_ndim
jgridloop:	 do jgrid = 1, pp_ndim
kgridloop:	  do kgrid = 1, pp_ndim
			gridnum = gridnum + 1
			if (.not. (grid_pp_int(iagrid,igrid,jgrid,kgrid))) cycle kgridloop


jgloop: do jg = 1, grid_pp_ngrp(gridnum) 
                jgr = grid_pp_grp(gridnum,jg)
#else
jgloop: do jgr = 1, ncgp_solute
#endif
! for every charge group:
! count each charge group pair once only
                if ( ((ig .gt. jgr) .and. (mod(ig+jgr,2) .eq. 0)) .or. &
                        ((ig .lt. jgr) .and. (mod(ig+jgr,2) .eq. 1)) ) &
                        cycle jgloop
! skip if excluded group
                ja = cgp(jgr)%iswitch
                if ( excl(ja) ) cycle jgloop
		r2 = q_dist4(x(ia),x(ja))

! skip if outside cutoff
                if ( r2 .gt. rcut2 ) cycle jgloop
ialoop:         do ia = cgp(ig)%first, cgp(ig)%last
! for every atom in the charge group (of the outermost loop):
                        i = cgpatom(ia)
! skip if q-atom
                        if ( iqatom(i).ne.0 ) cycle ialoop
jaloop:                 do ja = cgp(jgr)%first, cgp(jgr)%last
! for every atom in the charge group (innermost loop)
                                j = cgpatom(ja)
!                --- q-atom ? ---
                                if ( iqatom(j).ne.0 ) cycle jaloop
! make sure each pair is only counted once
                                if ( ig .eq. jgr .and. i .ge. j ) cycle jaloop
! do they interact?
                                if(pp_map(i-pp_low,j).eq.0) cycle jaloop
                                if(.not.pp_precomp(pp_map(i-pp_low,j))%set) cycle jaloop
! if out of space then make more space
!$omp critical
                                if (nbpp_pair .eq. calculation_assignment%pp%max) call reallocate_nonbondlist_pp
! all tests passed, add the pair
                                nbpp_pair = nbpp_pair + 1
                                nbpp(nbpp_pair)%i = i
                                nbpp(nbpp_pair)%j = j 
                                nbpp(nbpp_pair)%vdWA = pp_precomp(pp_map(i-pp_low,j))%vdWA
                                nbpp(nbpp_pair)%vdWB = pp_precomp(pp_map(i-pp_low,j))%vdWB
                                nbpp(nbpp_pair)%elec = pp_precomp(pp_map(i-pp_low,j))%elec
!$omp end critical

                        end do jaloop
                end do ialoop
        end do jgloop
#ifdef USE_GRID
end do kgridloop
end do jgridloop
end do igridloop
#endif
end do igloop
#if defined (PROFILING) 
profile(3)%time = profile(3)%time + rtime() - start_loop_time
#endif 

end subroutine nbpplist

!-----------------------------------------------------------------------

subroutine nbpplist_box(Rcut)
! args
real(kind=prec)						:: Rcut
! local variables
integer						:: i,j,ig,jg,ia,ja,nl,ig_sw, jg_sw,jgr
TYPE(qr_vec)					:: shift
real(kind=prec)						:: rcut2,r2
integer						:: LJ_code
integer                                                 ::iagrid, igrid, jgrid, kgrid, gridnum

! For use with periodic boundary conditions
!	This routine makes a list of non-bonded solute-solute atom pairs 
!	excluding any Q-atoms.

! uses the global variables:
!  Rcpp, ncgp, cgp, excl, x, cgpatom, iqatom, ljcod, crg, iaclib, max_nbr_range, listex
!  nexlong, listexlong, calculation_assignment%pp%max, alloc_status, list14, n14long, list14long

! reset #pairs
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)                                         :: start_loop_time
start_loop_time = rtime()
#endif
!$omp single
  nbpp_pair = 0 !atom pairs
  nbpp_cgp_pair = 0 !chargegroup pairs
!$omp end single
!$omp barrier
  rcut2 = Rcut

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%pp%end - calculation_assignment%pp%start + 1)/threads_num
remainder = MOD(calculation_assignment%pp%end - calculation_assignment%pp%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%pp%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

igloop: do ig = mp_start,mp_end
#else
igloop: do ig = calculation_assignment%pp%start, calculation_assignment%pp%end
#endif
	! for every assigned charge group:
	ig_sw = cgp(ig)%iswitch !switching atom in charge group ig
#ifdef USE_GRID
        iagrid = pp_igrid(ig)
	gridnum = 0
igridloop: do igrid = 1, pp_ndim
jgridloop:  do jgrid = 1, pp_ndim
kgridloop:   do kgrid = 1, pp_ndim
		gridnum = gridnum + 1
		if (.not. (grid_pp_int(iagrid,igrid,jgrid,kgrid))) cycle kgridloop

jgloop: do jg = 1, grid_pp_ngrp(gridnum)
		jgr = grid_pp_grp(gridnum,jg)
#else
jgloop: do jgr = 1, ncgp_solute
#endif
	  ! for every charge group:

	  ! count each charge group pair once only
	  if ( ((ig .gt. jgr) .and. (mod(ig+jgr,2) .eq. 0)) .or. &
		   ((ig .lt. jgr) .and. (mod(ig+jgr,2) .eq. 1)) ) &
		   cycle jgloop

	  jg_sw = cgp(jgr)%iswitch !switching atom in charge group jg
	  shift = x(ig_sw) - x(jg_sw)
          r2 = q_dist4(shift,boxlength*q_nint(shift*inv_boxl))

	  ! skip if outside cutoff
	  if ( r2 .gt. rcut2 ) cycle jgloop

	  !inside cutoff

	  !check if more memory is needed
!$omp critical
	  if(nbpp_cgp_pair .eq. size(nbpp_cgp, 1)) call reallocate_nbpp_cgp
	  !add the charge group pair
	  nbpp_cgp_pair = nbpp_cgp_pair + 1
	  nbpp_cgp(nbpp_cgp_pair)%i = ig_sw !the switching atoms of the charge groups in the pair
	  nbpp_cgp(nbpp_cgp_pair)%j = jg_sw
!$omp end critical

ialoop:         do ia = cgp(ig)%first, cgp(ig)%last
! for every atom in the charge group ig (of the outermost loop):
                        i = cgpatom(ia)
! skip if q-atom
                        if ( iqatom(i).ne.0 ) cycle ialoop
jaloop:                 do ja = cgp(jgr)%first, cgp(jgr)%last
! for every atom in the charge group jg (innermost loop)
                                j = cgpatom(ja)
!                --- q-atom ? ---
                                if ( iqatom(j).ne.0 ) cycle jaloop
! make sure each pair is only counted once
                                if ( ig .eq. jgr .and. i .ge. j ) cycle jaloop
! do they interact?
                                if(pp_map(i-pp_low,j).eq.0) cycle jaloop
                                if(.not.pp_precomp(pp_map(i-pp_low,j))%set) cycle jaloop
! if out of space then make more space
!$omp critical
                                if (nbpp_pair .eq. calculation_assignment%pp%max) call reallocate_nonbondlist_pp
! all tests passed, add the pair
                                nbpp_pair = nbpp_pair + 1
                                nbpp(nbpp_pair)%i = i
                                nbpp(nbpp_pair)%j = j 
                                nbpp(nbpp_pair)%vdWA = pp_precomp(pp_map(i-pp_low,j))%vdWA
                                nbpp(nbpp_pair)%vdWB = pp_precomp(pp_map(i-pp_low,j))%vdWB
                                nbpp(nbpp_pair)%elec = pp_precomp(pp_map(i-pp_low,j))%elec
                                nbpp(nbpp_pair)%cgp_pair = nbpp_cgp_pair !which pair of charge groups the atom pair belongs to
!$omp end critical		
                        end do jaloop
                end do ialoop
	end do jgloop
#ifdef USE_GRID
	end do kgridloop
	end do jgridloop
	end do igridloop
#endif
end do igloop
#if defined (PROFILING) 
profile(3)%time = profile(3)%time + rtime() - start_loop_time
#endif 

end subroutine nbpplist_box

!--------------------------------------------------------------------------------------

subroutine nbpplist_lrf(Rcut,RLRF)
! args
real(kind=prec)					:: Rcut,RLRF
! local variables
integer						:: i,j,ig,jg,ia,ja,nl,is,is3,jgr
real(kind=prec)						:: rcut2,r2
integer						:: LJ_code
real(kind=prec)						::	RcLRF2
integer                                                 ::iagrid, igrid, jgrid, kgrid, gridnum
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)						:: start_loop_time
start_loop_time = rtime()
#endif


!	This routine makes a list of non-bonded solute-solute atom pairs 
!	excluding any Q-atoms.

! uses the global variables:
!  nbpp_pair, Rcpp, RcLRF, cgp, excl, ncgp, x, cgpatom, iqatom, ljcod, crg, iaclib, max_nbr_range,
!  listex, nexlong, listexlong, nbpp, alloc_status, list14, n14long, list14long, lrf

!$omp single
nbpp_pair = 0
!$omp end single
!$omp barrier
rcut2 = Rcut
RcLRF2 = RLRF
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%pp%end - calculation_assignment%pp%start + 1)/threads_num
remainder = MOD(calculation_assignment%pp%end - calculation_assignment%pp%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%pp%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

igloop: do ig = mp_start,mp_end
#else
igloop: do ig = calculation_assignment%pp%start, calculation_assignment%pp%end
#endif
! for every assigned charge group:

! skip if excluded group
        is = cgp(ig)%iswitch
        if ( excl(is) ) cycle igloop
#ifdef USE_GRID
iagrid = pp_igrid(ig)
gridnum = 0
igridloop:	do igrid = 1, pp_ndim
jgridloop:	 do jgrid = 1, pp_ndim
kgridloop:	  do kgrid = 1, pp_ndim
			gridnum = gridnum + 1
			if (.not. (grid_pp_int(iagrid,igrid,jgrid,kgrid))) then
ggloop:			 do jg = 1, grid_pp_ngrp(gridnum)
			  jgr = grid_pp_grp(gridnum,jg)
                          ja = cgp(jgr)%iswitch
                          if ( excl(ja) ) cycle ggloop
			  if ( ((ig .gt. jgr) .and. (mod(ig+jgr,2) .eq. 0)) .or. &
			  ((ig .lt. jgr) .and. (mod(ig+jgr,2) .eq. 1)) .or. &
                          (ig.eq.jgr)) &
			  cycle ggloop
			  call lrf_update(ig,jgr)
			  call lrf_update(jgr,ig)
			 end do ggloop
			else
jgloop: do jg = 1, grid_pp_ngrp(gridnum)
                jgr = grid_pp_grp(gridnum,jg)
#else
jgloop: do jgr = 1, ncgp_solute
#endif
! count each charge group pair once only
                if ( ((ig .gt. jgr) .and. (mod(ig+jgr,2) .eq. 0)) .or. &
                        ((ig .lt. jgr) .and. (mod(ig+jgr,2) .eq. 1)) ) &
                        cycle jgloop
! skip if excluded group
                ja = cgp(jgr)%iswitch
                if ( excl(ja) ) cycle jgloop
		r2 = q_dist4(x(is),x(ja))

! inside cutoff?
                if ( r2 .le. rcut2 ) then
ialoop:                 do ia = cgp(ig)%first, cgp(ig)%last
! skip if q-atom
                                i = cgpatom(ia)
                                if ( iqatom(i).ne.0 ) cycle ialoop
jaloop:                         do ja = cgp(jgr)%first, cgp(jgr)%last
                                        j = cgpatom(ja)
! skip if q-atom
                                        if ( iqatom(j).ne.0 ) cycle jaloop
! make sure each pair is only counted once
                                        if ( ig .eq. jgr .and. i .ge. j ) cycle jaloop
! do they interact?
                                        if(pp_map(i-pp_low,j).eq.0) cycle jaloop
                                        if(.not.pp_precomp(pp_map(i-pp_low,j))%set) cycle jaloop
! if out of space then make more space
!$omp critical
                                        if (nbpp_pair .eq. calculation_assignment%pp%max) call reallocate_nonbondlist_pp
! all tests passed, add the pair
                                        nbpp_pair = nbpp_pair + 1
                                        nbpp(nbpp_pair)%i = i
                                        nbpp(nbpp_pair)%j = j
                                        nbpp(nbpp_pair)%vdWA = pp_precomp(pp_map(i-pp_low,j))%vdWA
                                        nbpp(nbpp_pair)%vdWB = pp_precomp(pp_map(i-pp_low,j))%vdWB
                                        nbpp(nbpp_pair)%elec = pp_precomp(pp_map(i-pp_low,j))%elec
!$omp end critical
                                end do jaloop
                        end do ialoop
                elseif(r2 .le. RcLRF2) then
! outside pp-cutoff but inside LRF cut-off use LRF
                        call lrf_update(ig,jgr)
                        call lrf_update(jgr,ig)
                end if ! LRF cutoff
        end do jgloop
#ifdef USE_GRID
end if ! interaction
end do kgridloop
end do jgridloop
end do igridloop
#endif
end do igloop
#if defined (PROFILING)
profile(3)%time = profile(3)%time + rtime() - start_loop_time
#endif

end subroutine nbpplist_lrf
!----------------LRF version of PW PBC-----------------------
subroutine nbpplist_box_lrf(Rcut,RLRF)
! args
real(kind=prec)						:: Rcut,RLRF
! local variables
integer						:: i,j,ig,jg,ia,ja,nl,ig_sw, jg_sw, is3, jgr
real(kind=prec)						:: rcut2,r2
TYPE(qr_vec)					:: shift
integer						:: LJ_code

real(kind=prec)						::RcLRF2

integer                                                 ::iagrid, igrid, jgrid, kgrid, gridnum
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)                                         :: start_loop_time
start_loop_time = rtime()
#endif

  ! For use with periodic boundary conditions
  !	This routine makes a list of non-bonded solute-solute atom pairs 
  !	excluding any Q-atoms.

  ! uses the global variables:
  !  Rcpp, ncgp, cgp, excl, x, cgpatom, iqatom, ljcod, crg, iaclib, max_nbr_range, listex
  !  nexlong, listexlong, calculation_assignment%pp%max, alloc_status, list14, n14long, list14long

  ! reset #pairs
!$omp single
  nbpp_pair = 0 !atom pairs
  nbpp_cgp_pair = 0 !chargegroup pairs
!$omp end single
!$omp barrier
  rcut2 = Rcut
  RcLRF2 = RLRF
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%pp%end - calculation_assignment%pp%start + 1)/threads_num
remainder = MOD(calculation_assignment%pp%end - calculation_assignment%pp%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%pp%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

igloop: do ig = mp_start,mp_end
#else
igloop: do ig = calculation_assignment%pp%start, calculation_assignment%pp%end
#endif
	! for every assigned charge group:
	ig_sw = cgp(ig)%iswitch !switching atom in charge group ig
#ifdef USE_GRID
	iagrid = pp_igrid(ig)
	gridnum = 0
igridloop:	do igrid = 1, pp_ndim
jgridloop:	 do jgrid = 1, pp_ndim
kgridloop:	  do kgrid = 1, pp_ndim
			gridnum = gridnum + 1
			if (.not. (grid_pp_int(iagrid,igrid,jgrid,kgrid))) then
ggloop:			 do jg = 1, grid_pp_ngrp(gridnum)
			  jgr = grid_pp_grp(gridnum,jg)
			  if ( ((ig .gt. jgr) .and. (mod(ig+jgr,2) .eq. 0)) .or. &
	                   ((ig .lt. jgr) .and. (mod(ig+jgr,2) .eq. 1)) .or. &
                           (ig.eq.jgr)) &
        	           cycle ggloop
			  call lrf_update(ig,jgr)
			  call lrf_update(jgr,ig)
			 end do ggloop
			 else
jgloop:	do jg = 1, grid_pp_ngrp(gridnum)
	        jgr = grid_pp_grp(gridnum,jg)
#else
jgloop: do jgr = 1, ncgp_solute
#endif
! for every charge group:
! count each charge group pair once only
                if ( ((ig .gt. jgr) .and. (mod(ig+jgr,2) .eq. 0)) .or. &
                        ((ig .lt. jgr) .and. (mod(ig+jgr,2) .eq. 1)) ) &
                        cycle jgloop
                jg_sw = cgp(jgr)%iswitch !switching atom in charge group jg
	 	shift = x(ig_sw) -  x(jg_sw)
                r2 = q_dist4(shift,boxlength*q_nint(shift*inv_boxl))
	
! skip if outside cutoff
                if ( r2 .le. rcut2 ) then 
!inside cutoff
!check if more memory is needed
!$omp critical
                        if(nbpp_cgp_pair .eq. size(nbpp_cgp, 1)) call reallocate_nbpp_cgp
!add the charge group pair
                        nbpp_cgp_pair = nbpp_cgp_pair + 1
                        nbpp_cgp(nbpp_cgp_pair)%i = ig_sw !the switching atoms of the charge groups in the pair
                        nbpp_cgp(nbpp_cgp_pair)%j = jg_sw
!$omp end critical
ialoop:                 do ia = cgp(ig)%first, cgp(ig)%last
! for every atom in the charge group ig (of the outermost loop):
                                i = cgpatom(ia)
!skip if q-atom
                                if ( iqatom(i).ne.0 ) cycle ialoop

jaloop:                         do ja = cgp(jgr)%first, cgp(jgr)%last
! for every atom in the charge group jg (innermost loop)
                                        j = cgpatom(ja)
!                --- q-atom ? ---
                                        if ( iqatom(j).ne.0 ) cycle jaloop
! make sure each pair is only counted once
                                        if ( ig .eq. jgr .and. i .ge. j ) cycle jaloop
! do they interact?
                                        if(pp_map(i-pp_low,j).eq.0) cycle jaloop
                                        if(.not.pp_precomp(pp_map(i-pp_low,j))%set) cycle jaloop
! if out of space then make more space
!$omp critical
                                        if (nbpp_pair .eq. calculation_assignment%pp%max) call reallocate_nonbondlist_pp
! all tests passed, add the pair
                                        nbpp_pair = nbpp_pair + 1
                                        nbpp(nbpp_pair)%i = i
                                        nbpp(nbpp_pair)%j = j 
                                        nbpp(nbpp_pair)%vdWA = pp_precomp(pp_map(i-pp_low,j))%vdWA
                                        nbpp(nbpp_pair)%vdWB = pp_precomp(pp_map(i-pp_low,j))%vdWB
                                        nbpp(nbpp_pair)%elec = pp_precomp(pp_map(i-pp_low,j))%elec
                                        nbpp(nbpp_pair)%cgp_pair = nbpp_cgp_pair !which pair of charge groups the atom pair belongs to
!$omp end critical
                                end do jaloop
                        end do ialoop
                elseif ((r2 .le. RcLRF2) .or. (RLRF .eq. -one)) then
! outside pp-cutoff but inside LRF cut-off use LRF
!ig : jg calculation
                        call lrf_update(ig,jgr)
!jg : ig calculations
                        call lrf_update(jgr,ig)
                end if ! outside cutoff

        end do jgloop
#ifdef USE_GRID
	end if ! interaction
	end do kgridloop
	end do jgridloop
	end do igridloop
#endif
end do igloop
#if defined (PROFILING)
profile(3)%time = profile(3)%time + rtime() - start_loop_time
#endif

end subroutine nbpplist_box_lrf
!-----------------------------------------------------------------------
!******PWchanged 2001-10-01
subroutine nbpw_count(npw, npwcgp,Rcut)
! arguments
integer						:: npw
integer						:: npwcgp(:)
real(kind=prec)                                 :: rcut
! local variables
integer						:: i,ig,jg,ia,ja
real(kind=prec)						:: rcut2,r2
integer						:: LJ_code
TYPE(qr_vec)					:: shift

!	This routine makes a list of non-bonded solute-solvent atom pairs
!	excluding Q-atoms.

! uses the global variables:
!  Rcpw, ncgp, cgp, excl, nwat, nat_solute, x, cgpatom, iqatom, ljcod, crg, iaclib

npw = 0
rcut2 = Rcut**2

igloop: do ig = 1, ncgp_solute
! for each charge group of the protein:
npwcgp(ig) = 0

! skip if excluded charge group
ia = cgp(ig)%iswitch
if ( .not.use_PBC .and. excl(ia) ) cycle igloop


jgloop: do jg = 1, nwat
  ! for each water molecule:
ja = nat_solute + solv_atom*jg-(solv_atom-1)
if(.not. use_PBC .and. excl(ja) ) cycle jgloop ! skip excluded waters


  if( .not. use_PBC ) then
	r2 = q_dist4(x(ia),x(ja))
  else
	shift = x(ia) -  x(ja)
        r2 = q_dist4(shift,boxlength*q_nint(shift*inv_boxl))
  end if

! skip if outside cutoff
if ( r2 .gt. rcut2 ) cycle jgloop

ialoop:  do ia = cgp(ig)%first, cgp(ig)%last
! for each atom in the protein charge group:
i = cgpatom(ia)

! skip if q-atom
if ( iqatom(i)/=0 ) cycle ialoop

jaloop:	do ja = nat_solute + solv_atom*jg-(solv_atom-1), nat_solute + solv_atom*jg
          ! for every atom of the water molecule:

          ! calculate LJ_code for the pair
          LJ_code = ljcod(iac(i),iac(ja))

          ! skip pairs with zero interaction
          if((crg(i) * crg(ja) == zero) &
                  .and. &
                  (iaclib(iac(i))%avdw(LJ_code)*iaclib(iac(ja))%avdw(LJ_code) == zero) &
                  .and. &
                  (iaclib(iac(i))%bvdw(LJ_code)*iaclib(iac(ja))%bvdw(LJ_code) == zero)) &
                  cycle jaloop

          ! count the pair
          npw = npw + 1
          npwcgp(ig) = npwcgp(ig) + 1

        end do jaloop
end do ialoop
end do jgloop
end do igloop

end subroutine nbpw_count

!-----------------------------------------------------------------------

subroutine nbpwlis2(Rcut)
! args
real(kind=prec)					:: Rcut
! local variables
integer						:: i,ig,jg,ia,ja, inside,jgr,j
real(kind=prec)						:: rcut2,r2
integer                                                 ::iagrid, igrid, jgrid, kgrid, gridnum, LJ_code
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)						:: start_loop_time
start_loop_time = rtime()
#endif

!	This routine makes a list of non-bonded solute-solvent atom pairs
!	excluding Q-atoms.

!$omp single
nbpw_pair = 0
!$omp end single
!$omp barrier
rcut2 = Rcut

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%pw%end - calculation_assignment%pw%start + 1)/threads_num
remainder = MOD(calculation_assignment%pw%end - calculation_assignment%pw%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%pw%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

igloop: do ig = mp_start,mp_end
#else
igloop: do ig = calculation_assignment%pw%start, calculation_assignment%pw%end
#endif
!	   --- excluded group ? ---
        ia = cgp(ig)%iswitch
        if ( excl(ia) ) cycle igloop
#ifdef USE_GRID
        iagrid = pp_igrid(ig)
igridloop:	do igrid = 1, pw_ndim
jgridloop:	 do jgrid = 1, pw_ndim
kgridloop:	  do kgrid = 1, pw_ndim
			gridnum = gridnum + 1
			if (.not. (grid_pw_int(iagrid,igrid,jgrid,kgrid))) cycle kgridloop
jgloop: do jg = 1, grid_pw_ngrp(gridnum)
                jgr = grid_pw_grp(gridnum,jg)
#else
jgloop: do jgr = 1, nwat
#endif
                ja = nat_solute + solv_atom*jgr-(solv_atom-1)
                if ( excl(ja) ) cycle jgloop ! skip excluded waters

!	      --- outside cutoff ? ---
                inside = 0
                ia = cgp(ig)%first
                do while ((ia .le. cgp(ig)%last) .and. (inside .eq. 0))
                        i = cgpatom(ia)
			r2 = q_dist4(x(ia),x(ja))

                        if ( r2 .le. rcut2 ) then
! inside cutoff, raise the flag
                                inside = 1
                        end if
                        ia = ia + 1
                end do
                if (inside .eq. 0) cycle jgloop

ialoop:         do ia = cgp(ig)%first, cgp(ig)%last
!	         --- q-atom ? ---
                        i = cgpatom(ia)
                        if ( iqatom(i)/=0 ) cycle ialoop
! if out of space then make more space
!$omp critical
                        if (nbpw_pair .gt. calculation_assignment%pw%max-solv_atom) call reallocate_nonbondlist_pw
jaloop:                 do j = 1, solv_atom
                                ja = nat_solute + (solv_atom*jgr) - solv_atom + j
! add the pair
                                nbpw_pair = nbpw_pair + 1
                                nbpw(nbpw_pair)%i = i
                                nbpw(nbpw_pair)%j = ja
                                nbpw(nbpw_pair)%vdWA = pw_precomp(i-pw_low,j)%vdWA
                                nbpw(nbpw_pair)%vdWB = pw_precomp(i-pw_low,j)%vdWB
                                nbpw(nbpw_pair)%elec = pw_precomp(i-pw_low,j)%elec
                                nbpw(nbpw_pair)%cgp_pair = nbpw_cgp_pair
                        end do jaloop
!$omp end critical
                end do ialoop
        end do jgloop
#ifdef USE_GRID
end do kgridloop
end do jgridloop
end do igridloop
#endif
end do igloop
#if defined (PROFILING)
profile(4)%time = profile(4)%time + rtime() - start_loop_time
#endif

end subroutine nbpwlis2

!------------------------------------------------------------------------------

!******PWadded 2001-10-18
subroutine nbpwlis2_box(Rcut)
! args
real(kind=prec)						:: Rcut
! local variables
integer						:: i,ig,jg,ia,ja,ig_atom, inside, jgr, LJ_code,j
TYPE(qr_vec)					:: shift
real(kind=prec)						:: rcut2,r2
integer                                                 ::iagrid, igrid, jgrid, kgrid, gridnum
! for periodic boundary conditions
!	This routine makes a list of non-bonded solute-solvent atom pairs
!	excluding Q-atoms.
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)						:: start_loop_time
start_loop_time = rtime()
#endif
!$omp single
nbpw_pair = 0
nbpw_cgp_pair = 0
!$omp end single
!$omp barrier
rcut2 = Rcut

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%pw%end - calculation_assignment%pw%start + 1)/threads_num
remainder = MOD(calculation_assignment%pw%end - calculation_assignment%pw%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%pw%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

igloop: do ig = mp_start,mp_end
#else
igloop: do ig = calculation_assignment%pw%start, calculation_assignment%pw%end
#endif
#ifdef USE_GRID
	iagrid = pw_igrid(ig)
	gridnum = 0
igridloop: do igrid = 1, pw_ndim
jgridloop:  do jgrid = 1, pw_ndim
kgridloop:   do kgrid = 1, pw_ndim
		gridnum = gridnum + 1
		if (.not. (grid_pw_int(iagrid,igrid,jgrid,kgrid))) cycle kgridloop

jgloop:	do jg = 1, grid_pw_ngrp(gridnum)
                jgr = grid_pw_grp(gridnum,jg)
#else
jgloop: do jgr = 1, nwat
#endif
                ja = nat_solute + solv_atom*jgr-(solv_atom-1)
!	      --- outside cutoff ? ---
		inside = 0
		ig_atom = cgp(ig)%first
		do while ((ig_atom .le. cgp(ig)%last) .and. (inside .eq. 0))
                   	i = cgpatom(ig_atom)
			shift = x(i) - x(ja)
                        r2 = q_dist4(shift,boxlength*q_nint(shift*inv_boxl))

                        if ( r2 .le. rcut2 ) then
! inside cutoff, raise the flag
                                inside = 1
!$omp critical
                                if (nbpw_cgp_pair .eq. size(nbpw_cgp, 1)) call reallocate_nbpw_cgp
                                nbpw_cgp_pair = nbpw_cgp_pair + 1
                                nbpw_cgp(nbpw_cgp_pair)%i = i
                                nbpw_cgp(nbpw_cgp_pair)%j = ja
!$omp end critical
                        end if
                        ig_atom = ig_atom + 1 !ia = ia + 1
                end do
		if (inside .eq. 0) cycle jgloop
ialoop:         do ia = cgp(ig)%first, cgp(ig)%last
!	         --- q-atom ? ---
                        i = cgpatom(ia)
                        if ( iqatom(i)/=0 ) cycle ialoop
! if out of space then make more space
!$omp critical
                        if (nbpw_pair .gt. calculation_assignment%pw%max-solv_atom) call reallocate_nonbondlist_pw
jaloop:                 do j = 1, solv_atom
                                ja = nat_solute + (solv_atom*jgr) - solv_atom + j
! for every atom of the water molecule:
! add the pair
                                nbpw_pair = nbpw_pair + 1
                                nbpw(nbpw_pair)%i = i
                                nbpw(nbpw_pair)%j = ja
                                nbpw(nbpw_pair)%vdWA = pw_precomp(i-pw_low,j)%vdWA
                                nbpw(nbpw_pair)%vdWB = pw_precomp(i-pw_low,j)%vdWB
                                nbpw(nbpw_pair)%elec = pw_precomp(i-pw_low,j)%elec
                                nbpw(nbpw_pair)%cgp_pair = nbpw_cgp_pair
                        end do jaloop
!$omp end critical
                end do ialoop
        end do jgloop
#ifdef USE_GRID
	end do kgridloop
	end do jgridloop
	end do igridloop
#endif
end do igloop
#if defined (PROFILING)
profile(4)%time = profile(4)%time + rtime() - start_loop_time
#endif

end subroutine nbpwlis2_box
!-----------------------------------------------------------------------
subroutine nbpwlis2_box_lrf(Rcut,RLRF)
! args
real(kind=prec)					:: Rcut,RLRF
! local variables
integer						:: i,ig,jg,ia,ja,ig_atom, inside, jgr, LJ_code, j
real(kind=prec)					:: rcut2,r2
TYPE(qr_vec)					:: shift
!LRF
real(kind=prec)					:: RcLRF2
integer						:: jg_cgp, inside_LRF, is3

integer                                         ::iagrid, igrid, jgrid, kgrid, gridnum

! for periodic boundary conditions
!	This routine makes a list of non-bonded solute-solvent atom pairs
!	excluding Q-atoms.
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)						:: start_loop_time
start_loop_time = rtime()
#endif
!$omp single
  nbpw_pair = 0
  nbpw_cgp_pair = 0
!$omp end single
!$omp barrier
  rcut2 = Rcut
  RcLRF2 = RLRF

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%pw%end - calculation_assignment%pw%start + 1)/threads_num
remainder = MOD(calculation_assignment%pw%end - calculation_assignment%pw%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%pw%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

igloop: do ig = mp_start,mp_end
#else
igloop: do ig = calculation_assignment%pw%start, calculation_assignment%pw%end
#endif
#ifdef USE_GRID
	iagrid = pw_igrid(ig)
igridloop:	do igrid = 1, pw_ndim
jgridloop:	 do jgrid = 1, pw_ndim
kgridloop:	  do kgrid = 1, pw_ndim
			gridnum = gridnum + 1
			if (.not. (grid_pw_int(iagrid,igrid,jgrid,kgrid))) then
ggloop:			 do jg = 1, grid_pw_ngrp(gridnum)
			  jgr = grid_pw_grp(gridnum,jg)
			  ja = nat_solute + solv_atom*jgr-(solv_atom-1)
			  jg_cgp = iwhich_cgp(ja)
			  call lrf_update(ig,jg_cgp)
			  call lrf_update(jg_cgp,ig)
			 end do ggloop
			else
jgloop:	do jg = 1, grid_pw_ngrp(gridnum)
                jgr = grid_pw_grp(gridnum,jg)
#else
jgloop: do jgr = 1, nwat
#endif
                ja = nat_solute + solv_atom*jgr-(solv_atom-1)
                jg_cgp = iwhich_cgp(ja)
!	      --- outside cutoff ? ---
                inside = 0
                inside_LRF = 0
                ig_atom = cgp(ig)%first
                do while ((ig_atom .le. cgp(ig)%last) .and. (inside .eq. 0))
                        i = cgpatom(ig_atom)
                        shift = x(i) - x(ja)
                        r2 = q_dist4(shift,boxlength*q_nint(shift*inv_boxl))
                        if ( r2 .le. rcut2 ) then
! inside cutoff, raise the flag
                                inside = 1
!$omp critical
                                if (nbpw_cgp_pair .eq. size(nbpw_cgp, 1)) call reallocate_nbpw_cgp
                                nbpw_cgp_pair = nbpw_cgp_pair + 1
                                nbpw_cgp(nbpw_cgp_pair)%i = i
                                nbpw_cgp(nbpw_cgp_pair)%j = ja
!$omp end critical
                        elseif ((r2 .le. RcLRF2) .or. (RLRF .eq. -one)) then
                                inside_LRF = 1
                        end if
                        ig_atom = ig_atom + 1 !ia = ia + 1
                end do
                if (inside .eq. 1) then
ialoop:                 do ia = cgp(ig)%first, cgp(ig)%last
!	         --- q-atom ? ---
                        i = cgpatom(ia)
                        if ( iqatom(i)/=0 ) cycle ialoop
! if out of space then make more space
!$omp critical
                                if (nbpw_pair .gt. calculation_assignment%pw%max - 3) call reallocate_nonbondlist_pw
jaloop:                         do j = 1, solv_atom
                                        ja = nat_solute + (solv_atom*jgr) - solv_atom + j 
! for every atom of the water molecule:
! add the pair
                                        nbpw_pair = nbpw_pair + 1
                                        nbpw(nbpw_pair)%i = i
                                        nbpw(nbpw_pair)%j = ja 
                                        nbpw(nbpw_pair)%vdWA = pw_precomp(i-pw_low,j)%vdWA
                                        nbpw(nbpw_pair)%vdWB = pw_precomp(i-pw_low,j)%vdWB
                                        nbpw(nbpw_pair)%elec = pw_precomp(i-pw_low,j)%elec
                                        nbpw(nbpw_pair)%cgp_pair = nbpw_cgp_pair
                                end do jaloop
!$omp end critical

                        end do ialoop
                elseif((inside_LRF .eq. 1) .and. (inside .eq. 0)) then   
! outside pw-cutoff but inside LRF cut-off: use LRF
!solut : solvent
                        call lrf_update(ig,jg_cgp)
!solvent : solut	
                        call lrf_update(jg_cgp,ig)
                end if
        end do jgloop
#ifdef USE_GRID
	end if ! interaction
end do kgridloop
end do jgridloop
end do igridloop
#endif
end do igloop
#if defined (PROFILING)
profile(4)%time = profile(4)%time + rtime() - start_loop_time
#endif

end subroutine nbpwlis2_box_lrf
!-----------------------------------------------------------------------
subroutine nbpwlis2_lrf(Rcut,RLRF)
! args
real(kind=prec)					:: Rcut, RLRF
! local variables
integer						:: i,j,ig,jg,jg_cgp,ia,ja,inside,is,jgr
real(kind=prec)						:: rcut2,r2
real(kind=prec)						::	RcLRF2
integer                                                 ::iagrid, igrid, jgrid, kgrid, gridnum, LJ_code
#ifdef _OPENMP
integer :: quotient, remainder
#endif

!	This routine makes a list of non-bonded solute-solvent atom pairs
!	excluding Q-atoms.
#if defined (PROFILING)
real(kind=prec)						:: start_loop_time
start_loop_time = rtime()
#endif
!$omp single
nbpw_pair = 0
!$omp end single
!$omp barrier
rcut2 = Rcut
RcLRF2 = RLRF
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%pw%end - calculation_assignment%pw%start + 1)/threads_num
remainder = MOD(calculation_assignment%pw%end - calculation_assignment%pw%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%pw%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

igloop: do ig = mp_start,mp_end
#else
igloop: do ig = calculation_assignment%pw%start, calculation_assignment%pw%end
#endif

!	   --- excluded group ? ---
        is = cgp(ig)%iswitch
        if ( excl(is) ) cycle igloop
#ifdef USE_GRID
        iagrid = pw_igrid(ig)
        gridnum = 0
igridloop:	do igrid = 1, pw_ndim
jgridloop:	 do jgrid = 1, pw_ndim
kgridloop:	  do kgrid = 1, pw_ndim
			gridnum = gridnum + 1
			if (.not. (grid_pw_int(iagrid,igrid,jgrid,kgrid))) then
ggloop:			do jg = 1, grid_pw_ngrp(gridnum)
			 jgr = grid_pw_grp(gridnum,jg)
			 ja = nat_solute + solv_atom*jgr-(solv_atom-1)
			 if(excl(ja)) cycle ggloop ! skip excluded waters
			 jg_cgp = iwhich_cgp(ja)
			 call lrf_update(ig,jg_cgp)
			 call lrf_update(jg_cgp,ig)
			end do ggloop
			else
jgloop:	do jg = 1, grid_pw_ngrp(gridnum)
	        jgr = grid_pw_grp(gridnum,jg)
#else
jgloop: do jgr = 1, nwat
#endif
                ja = nat_solute + solv_atom*jgr-(solv_atom-1)
                if(excl(ja)) cycle jgloop
                jg_cgp = iwhich_cgp(ja)

!	      --- outside cutoff ? ---
                inside = 0
                ia = cgp(ig)%first
                do while ((ia .le. cgp(ig)%last) .and. (inside .eq. 0 ))
                        i = cgpatom(ia)
			r2 = q_dist4(x(is),x(ja))
                        if ( r2 .le. rcut2 ) then
! inside cutoff, raise the flag
                                inside = 1
                        end if
                        ia = ia + 1
                end do

                if ( inside .eq. 1 ) then
ialoop:                 do ia = cgp(ig)%first, cgp(ig)%last
!	            --- q-atom ? ---
                                i = cgpatom(ia)
                                if ( iqatom(i)/=0 ) cycle ialoop
! if out of space then make more space
!$omp critical
                                if (nbpw_pair .gt. calculation_assignment%pw%max-solv_atom) call reallocate_nonbondlist_pw
jaloop:                         do j = 1, solv_atom 
                                        ja = nat_solute + (solv_atom*jgr) - solv_atom + j
! for every atom of the water molecule:
! add the pair
                                        nbpw_pair = nbpw_pair + 1
                                        nbpw(nbpw_pair)%i = i
                                        nbpw(nbpw_pair)%j = ja
                                        nbpw(nbpw_pair)%vdWA = pw_precomp(i-pw_low,j)%vdWA
                                        nbpw(nbpw_pair)%vdWB = pw_precomp(i-pw_low,j)%vdWB
                                        nbpw(nbpw_pair)%elec = pw_precomp(i-pw_low,j)%elec
                                        nbpw(nbpw_pair)%cgp_pair = nbpw_cgp_pair
                                end do jaloop
!$omp end critical
                        end do ialoop

                elseif(r2 .le. RcLRF2) then   
! outside pw-cutoff but inside LRF cut-off: use LRF
                	call lrf_update(ig,jg_cgp)
	                call lrf_update(jg_cgp,ig)
                end if
        end do jgloop
#ifdef USE_GRID
end if ! interaction
end do kgridloop
end do jgridloop
end do igridloop
#endif
end do igloop
#if defined (PROFILING)
profile(4)%time = profile(4)%time + rtime() - start_loop_time
#endif

end subroutine nbpwlis2_lrf

!-----------------------------------------------------------------------

subroutine nbpwlist(Rcut)
! args
real(kind=prec)					:: Rcut
! local variables
integer						:: i,ig,jg,ia,ja,jgr,j
real(kind=prec)						:: rcut2,r2
integer						:: LJ_code
integer                                                 ::iagrid, igrid, jgrid, kgrid, gridnum

! for use with spherical boundary

!	This routine makes a list of non-bonded solute-solvent atom pairs
!	excluding Q-atoms.
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)						:: start_loop_time
start_loop_time = rtime()
#endif

! reset nbpw_pair
!$omp single
nbpw_pair = 0
!$omp end single
!$omp barrier
rcut2 = Rcut
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%pw%end - calculation_assignment%pw%start + 1)/threads_num
remainder = MOD(calculation_assignment%pw%end - calculation_assignment%pw%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%pw%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

igloop: do ig = mp_start,mp_end
#else
igloop: do ig = calculation_assignment%pw%start, calculation_assignment%pw%end
#endif
! for every charge group:

! skip excluded groups
        ia = cgp(ig)%iswitch
        if ( excl(ia) ) cycle igloop
#ifdef USE_GRID
        iagrid = pw_igrid(ig)
        gridnum = 0
igridloop:	do igrid = 1, pw_ndim
jgridloop:	 do jgrid = 1, pw_ndim
kgridloop:	  do kgrid = 1, pw_ndim
			gridnum = gridnum + 1
			if (.not. (grid_pw_int(iagrid,igrid,jgrid,kgrid))) cycle kgridloop
jgloop:	do jg = 1, grid_pw_ngrp(gridnum)
	        jgr = grid_pw_grp(gridnum,jg)
#else
jgloop: do jgr = 1, nwat
#endif
! for every water molecule in grid:
                ja = nat_solute + solv_atom*jgr-(solv_atom-1)
                if(excl(ja)) cycle jgloop ! skip excluded waters
		r2 = q_dist4(x(ia),x(ja))

! skip water outside cutoff
                if ( r2 .gt. rcut2 ) cycle jgloop

ialoop:         do ia = cgp(ig)%first, cgp(ig)%last
! for every atom in the charge group:

! find the atom index of the atom in the charge group
                        i = cgpatom(ia)
!	skip q-atoms
                        if ( iqatom(i)/=0 ) cycle ialoop
!$omp critical
! if out of space then make more space
                        if (nbpw_pair .gt. calculation_assignment%pw%max - solv_atom) call reallocate_nonbondlist_pw
jaloop:                 do j = 1, solv_atom
                                ja = nat_solute + (solv_atom*jgr) - solv_atom + j
! for every atom of the water molecule:
! add the pair
                                nbpw_pair = nbpw_pair + 1
                                nbpw(nbpw_pair)%i = i
                                nbpw(nbpw_pair)%j = ja 
                                nbpw(nbpw_pair)%vdWA = pw_precomp(i-pw_low,j)%vdWA
                                nbpw(nbpw_pair)%vdWB = pw_precomp(i-pw_low,j)%vdWB
                                nbpw(nbpw_pair)%elec = pw_precomp(i-pw_low,j)%elec
                        end do jaloop
!$omp end critical
                end do ialoop
        end do jgloop
#ifdef USE_GRID
end do kgridloop
end do jgridloop
end do igridloop
#endif
end do igloop
#if defined (PROFILING)
profile(4)%time = profile(4)%time + rtime() - start_loop_time
#endif

end subroutine nbpwlist

!-----------------------------------------------------------------------
!******PWadded 2001-10-18

subroutine nbpwlist_box(Rcut)
! args
real(kind=prec)					:: Rcut
! local variables
integer						:: i,ig,jg,ia,ja,ig_sw,jg_sw,jgr,j
real(kind=prec)						:: rcut2,r2
TYPE(qr_vec)					:: shift
integer						:: LJ_code
integer                                                 ::iagrid, igrid, jgrid, kgrid, gridnum
! For use with periodic boundary conditions
!	This routine makes a list of non-bonded solute-solvent atom pairs
!	excluding Q-atoms.
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)						:: start_loop_time
start_loop_time = rtime()
#endif

! reset nbpw_pair
!$omp single
nbpw_pair = 0
nbpw_cgp_pair = 0
!$omp end single
!$omp barrier
rcut2 = Rcut
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%pw%end - calculation_assignment%pw%start + 1)/threads_num
remainder = MOD(calculation_assignment%pw%end - calculation_assignment%pw%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%pw%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

igloop: do ig = mp_start,mp_end
#else
igloop: do ig = calculation_assignment%pw%start, calculation_assignment%pw%end
#endif
! for every charge group:
        ig_sw = cgp(ig)%iswitch
#ifdef USE_GRID
        iagrid = pw_igrid(ig)
        gridnum = 0
igridloop:	do igrid = 1, pw_ndim
jgridloop:	 do jgrid = 1, pw_ndim
kgridloop:	  do kgrid = 1, pw_ndim
			gridnum = gridnum + 1
			if (.not. (grid_pw_int(iagrid,igrid,jgrid,kgrid))) cycle kgridloop
jgloop: do jg = 1, grid_pw_ngrp(gridnum)
	        jgr = grid_pw_grp(gridnum,jg)
#else
jgloop: do jgr = 1, nwat
#endif
! for every water molecule:
                jg_sw = nat_solute + solv_atom*jgr-(solv_atom-1)
		shift = x(ig_sw) - x(jg_sw)
                r2 = q_dist4(shift,boxlength*q_nint(shift*inv_boxl))

! skip water outside cutoff
                if ( r2 .gt. rcut2 ) cycle jgloop

	  !inside cut-off
!check if charge group pair list is big enough
!$omp critical
                if(nbpw_cgp_pair .eq. size(nbpw_cgp, 1) ) call reallocate_nbpw_cgp
                nbpw_cgp_pair = nbpw_cgp_pair + 1
                nbpw_cgp(nbpw_cgp_pair)%i = ig_sw  !solute
                nbpw_cgp(nbpw_cgp_pair)%j = jg_sw  !water
!$omp end critical
ialoop:         do ia = cgp(ig)%first, cgp(ig)%last
! for every atom in the charge group:

! find the atom index of the atom in the charge group
                        i = cgpatom(ia)
!	skip q-atoms
                        if ( iqatom(i)/=0 ) cycle ialoop
! if out of space then make more space
!$omp critical
                        if (nbpw_pair .gt. calculation_assignment%pw%max - solv_atom) call reallocate_nonbondlist_pw
jaloop:                 do j = 1, solv_atom
                                ja = nat_solute + (solv_atom*jgr) - solv_atom + j
! for every atom of the water molecule:
! add the pair
                                nbpw_pair = nbpw_pair + 1
                                nbpw(nbpw_pair)%i = i
                                nbpw(nbpw_pair)%j = ja 
                                nbpw(nbpw_pair)%vdWA = pw_precomp(i-pw_low,j)%vdWA
                                nbpw(nbpw_pair)%vdWB = pw_precomp(i-pw_low,j)%vdWB
                                nbpw(nbpw_pair)%elec = pw_precomp(i-pw_low,j)%elec
                                nbpw(nbpw_pair)%cgp_pair = nbpw_cgp_pair
                        end do jaloop
!$omp end critical
                end do ialoop
        end do jgloop
#ifdef USE_GRID
end do kgridloop
end do jgridloop
end do igridloop
#endif
  end do igloop
#if defined (PROFILING)
profile(4)%time = profile(4)%time + rtime() - start_loop_time
#endif

end subroutine nbpwlist_box
!-------------------------------------------------------------------------------------

subroutine nbpwlist_lrf(Rcut,RLRF)
! args
real(kind=prec)					:: Rcut,RLRF
! local variables
integer						:: i,j,ig,jg,jg_cgp,ia,ja,is,jgr
real(kind=prec)						:: rcut2,r2
integer						:: LJ_code
real(kind=prec)						::	RcLRF2
integer                                                 ::iagrid, igrid, jgrid, kgrid, gridnum
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)						:: start_loop_time
start_loop_time = rtime()
#endif

!	This routine makes a list of non-bonded solute-solvent atom pairs
!	excluding Q-atoms.

! reset nbpw_pair
!$omp single
nbpw_pair = 0
!$omp end single
!$omp barrier
rcut2 = Rcut
RcLRF2 = RLRF
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%pw%end - calculation_assignment%pw%start + 1)/threads_num
remainder = MOD(calculation_assignment%pw%end - calculation_assignment%pw%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%pw%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

igloop: do ig = mp_start,mp_end
#else
igloop: do ig = calculation_assignment%pw%start, calculation_assignment%pw%end
#endif
! for every charge group:

! skip excluded groups
        is = cgp(ig)%iswitch
        if ( excl(is) ) cycle igloop
#ifdef USE_GRID
        iagrid = pw_igrid(ig)
        gridnum = 0
igridloop:	do igrid = 1, pw_ndim
jgridloop:	 do jgrid = 1, pw_ndim
kgridloop:	  do kgrid = 1, pw_ndim
			gridnum = gridnum + 1
			if (.not. (grid_pw_int(iagrid,igrid,jgrid,kgrid))) then
ggloop:			do jg = 1, grid_pw_ngrp(gridnum)
			 jgr = grid_pw_grp(gridnum,jg)
			 ja = nat_solute + solv_atom*jgr-(solv_atom-1)
			 if(excl(ja)) cycle ggloop ! skip excluded waters
			 jg_cgp = iwhich_cgp(ja)
			 call lrf_update(ig,jg_cgp)
			 call lrf_update(jg_cgp,ig)
			end do ggloop
			else
jgloop:	do jg = 1, grid_pw_ngrp(gridnum)
	        jgr = grid_pw_grp(gridnum,jg)
#else
jgloop: do jgr = 1,nwat
#endif
                ja = nat_solute + solv_atom*jgr-(solv_atom-1)
                jg_cgp = iwhich_cgp(ja)
                if(excl(ja)) cycle jgloop
		r2 = q_dist4(x(is),x(ja))
                if ( r2 .le. rcut2 ) then
! within the cutoff radix:

ialoop:                 do ia = cgp(ig)%first, cgp(ig)%last
                                i = cgpatom(ia)
! skip q-atoms
                                if ( iqatom(i)/=0 ) cycle ialoop
!$omp critical
! if out of space then make more space
                                if (nbpw_pair .gt. calculation_assignment%pw%max - solv_atom) call reallocate_nonbondlist_pw
jaloop:                         do j = 1, solv_atom
                                        ja = nat_solute + (solv_atom*jgr) - solv_atom + j
! add the pair
                                        nbpw_pair = nbpw_pair + 1
                                        nbpw(nbpw_pair)%i = i
                                        nbpw(nbpw_pair)%j = ja 
                                        nbpw(nbpw_pair)%vdWA = pw_precomp(i-pw_low,j)%vdWA
                                        nbpw(nbpw_pair)%vdWB = pw_precomp(i-pw_low,j)%vdWB
                                        nbpw(nbpw_pair)%elec = pw_precomp(i-pw_low,j)%elec
                                end do jaloop
!$omp end critical
                        end do ialoop
                elseif(r2 .le. RcLRF2) then   
! outside pw-cutoff but inside LRF cut-off: use LRF
                        call lrf_update(ig,jg_cgp)
                        call lrf_update(jg_cgp,ig)
                end if
        end do jgloop
#ifdef USE_GRID
end if ! interaction
end do kgridloop
end do jgridloop
end do igridloop
#endif
end do igloop
#if defined (PROFILING)
profile(4)%time = profile(4)%time + rtime() - start_loop_time
#endif

end subroutine nbpwlist_lrf
!---------------LRF version of PW PBC-----------------------
subroutine nbpwlist_box_lrf(Rcut,RLRF)
! args
real(kind=prec)					:: Rcut,RLRF
! local variables
integer						:: i,ig,jg,ia,ja,ig_sw,jg_sw,jgr
real(kind=prec)						:: rcut2,r2
integer						:: LJ_code
TYPE(qr_vec)					:: shift
! LRF
real(kind=prec)						:: RcLRF2
integer						:: jg_cgp, j, is3
integer                                                 ::iagrid, igrid, jgrid, kgrid, gridnum
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)                                         :: start_loop_time
start_loop_time = rtime()
#endif

  ! For use with periodic boundary conditions
  !	This routine makes a list of non-bonded solute-solvent atom pairs
  !	excluding Q-atoms.

  ! reset nbpw_pair
!$omp single
  nbpw_pair = 0
  nbpw_cgp_pair = 0
!$omp end single
!$omp barrier
  rcut2 = Rcut
  RcLRF2 = RLRF
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%pw%end - calculation_assignment%pw%start + 1)/threads_num
remainder = MOD(calculation_assignment%pw%end - calculation_assignment%pw%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%pw%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

igloop: do ig = mp_start,mp_end
#else
igloop: do ig = calculation_assignment%pw%start, calculation_assignment%pw%end
#endif
	! for every charge group:
        ig_sw = cgp(ig)%iswitch
#ifdef USE_GRID
        iagrid = pw_igrid(ig)
        gridnum = 0
igridloop:	do igrid = 1, pw_ndim
jgridloop:	 do jgrid = 1, pw_ndim
kgridloop:	  do kgrid = 1, pw_ndim
			gridnum = gridnum + 1
			if (.not. (grid_pw_int(iagrid,igrid,jgrid,kgrid))) then
ggloop:			do jg = 1, grid_pw_ngrp(gridnum)
			 jgr = grid_pw_grp(gridnum,jg)
			 jg_sw = nat_solute + (solv_atom*jgr-(solv_atom-1))
			 jg_cgp = iwhich_cgp(jg_sw)
			 call lrf_update(ig,jg_cgp)
			 call lrf_update(jg_cgp,ig)
			end do ggloop
			else
jgloop: do jg = 1, grid_pw_ngrp(gridnum)
	        jgr = grid_pw_grp(gridnum,jg)
#else
jgloop: do jgr = 1, nwat
#endif
                jg_sw = nat_solute + (solv_atom*jgr-(solv_atom-1))
                jg_cgp = iwhich_cgp(jg_sw)
		shift = x(ig_sw) - x(jg_sw)
                r2 = q_dist4(shift,boxlength*q_nint(shift*inv_boxl))

! skip water outside cutoff
                if ( r2 .le. rcut2 ) then
!inside cut-off
!check if charge group pair list is big enough
!$omp critical
                        if(nbpw_cgp_pair .eq. size(nbpw_cgp, 1) ) call reallocate_nbpw_cgp
                        nbpw_cgp_pair = nbpw_cgp_pair + 1
                        nbpw_cgp(nbpw_cgp_pair)%i = ig_sw  !solute
                        nbpw_cgp(nbpw_cgp_pair)%j = jg_sw  !water
!$omp end critical
ialoop:                 do ia = cgp(ig)%first, cgp(ig)%last
! for every atom in the charge group:
! find the atom index of the atom in the charge group
                                i = cgpatom(ia)
!	skip q-atoms
                                if ( iqatom(i)/=0 ) cycle ialoop
! if out of space then make more space
!$omp critical
                                if (nbpw_pair .gt. calculation_assignment%pw%max - solv_atom) call reallocate_nonbondlist_pw
jaloop:                         do j = 1, solv_atom
                                        ja = nat_solute + (solv_atom*jgr) - solv_atom + j
! for every atom of the solvent molecule:
! add the pair
                                        nbpw_pair = nbpw_pair + 1
                                        nbpw(nbpw_pair)%i = i
                                        nbpw(nbpw_pair)%j = ja 
                                        nbpw(nbpw_pair)%vdWA = pw_precomp(i-pw_low,j)%vdWA
                                        nbpw(nbpw_pair)%vdWB = pw_precomp(i-pw_low,j)%vdWB
                                        nbpw(nbpw_pair)%elec = pw_precomp(i-pw_low,j)%elec
                                        nbpw(nbpw_pair)%cgp_pair = nbpw_cgp_pair
                                end do jaloop
!$omp end critical
                        end do ialoop
                elseif((r2 .le. RcLRF2) .or. (RcLRF .eq. -one)) then   
! outside pw-cutoff but inside LRF cut-off: use LRF
!solut : solvent
                        call lrf_update(ig,jg_cgp)
!solvent : solut	
                        call lrf_update(jg_cgp,ig)
                end if
        end do jgloop
#ifdef USE_GRID
end if !interaction
end do kgridloop
end do jgridloop
end do igridloop
#endif
  end do igloop
#if defined (PROFILING)
profile(4)%time = profile(4)%time + rtime() - start_loop_time
#endif    

end subroutine nbpwlist_box_lrf
!-------------------------------------------------------------------------------------
subroutine make_qconn

integer						::	 i, iq, is, ia

allocate(qconn(nstates,nat_solute, nqat))
qconn(:,:,:) = 9

do iq = 1, nqat
        qconn(:, iqseq(iq), iq) = 1
end do

do iq = 1, nqat
        do is = 1, nstates
                i = iqseq(iq)
                call find_bonded(origin=i, current=i, level=1, state=is)
        end do
end do

!modify matrix to take special exclusions into account
do i = 1, nexspec
        iq = iqatom(exspec(i)%i)
        if(iq > 0) then
                do is = 1, nstates
                        if(exspec(i)%flag(is)) then
                                qconn(is, exspec(i)%j, iq) = 0 !exclude by setting to 0
                        end if
                end do
        end if
        iq = iqatom(exspec(i)%j)
        if(iq > 0) then
                do is = 1, nstates
                        if(exspec(i)%flag(is)) then
                                qconn(is, exspec(i)%i, iq) = 0 !exclude by setting to 0
                        end if
                end do
        end if
end do

end subroutine make_qconn


!------------------------------------------------------------------------------


recursive subroutine find_bonded(origin, current, level, state)
!args
integer, intent(in)			::	origin, current, level, state
!locals
integer						::	b, newcurrent, newlevel

!find q-atom connectivity using the bond list and the q-bond list
!shaken bonds (code -1) must be taken into account, but not 
!redefined bonds in the topology
do b = 1, nbonds_solute
        if(bnd(b)%cod == 0) cycle !skip redefined (but not shaken)
        if(bnd(b)%i  == current) then
                newlevel = level + 1 
                newcurrent = bnd(b)%j
                if(qconn(state, newcurrent, iqatom(origin)) > newlevel)  then
                        qconn(state, newcurrent, iqatom(origin)) = newlevel
                        if(newlevel < 4) then
                                call find_bonded(origin, newcurrent, newlevel, state)
                        end if
                end if
        elseif(bnd(b)%j  == current) then
                newlevel = level + 1 
                newcurrent = bnd(b)%i
                if(qconn(state, newcurrent, iqatom(origin)) > newlevel)  then
                        qconn(state, newcurrent, iqatom(origin)) = newlevel
                        if(newlevel < 4) then
                                call find_bonded(origin, newcurrent, newlevel, state)
                        end if
                end if
        end if
end do
do b = 1, nqbond
        if(qbnd(b)%cod(state) > 0) then
                if(qbnd(b)%i  == current) then
                        newlevel = level + 1 
                        newcurrent = qbnd(b)%j
                        if(qconn(state, newcurrent, iqatom(origin)) > newlevel)  then
                                qconn(state, newcurrent, iqatom(origin)) = newlevel
                                if(newlevel < 4) then
                                        call find_bonded(origin, newcurrent, newlevel, state)
                                end if
                        end if
                elseif(qbnd(b)%j  == current) then
                        newlevel = level + 1 
                        newcurrent = qbnd(b)%i
                        if(qconn(state, newcurrent, iqatom(origin)) > newlevel)  then
                                qconn(state, newcurrent, iqatom(origin)) = newlevel
                                if(newlevel < 4) then
                                        call find_bonded(origin, newcurrent, newlevel, state)
                                end if
                        end if
                end if
        end if
end do

end subroutine find_bonded


!------------------------------------------------------------------------------


integer function nbqq_count()

integer						::	iq, j, jq, is

allocate(nbqq_pair(nstates),nbqqp_pair(nstates))

nbqq_pair(:) = 0

!count Q-Q
do iq = 1, nqat - 1
        do jq = iq + 1, nqat
                do is = 1, nstates
                        if(qconn(is, iqseq(jq), iq) > 3) then
                                nbqq_pair(is) = nbqq_pair(is)+1
                        end if
                end do
        end do
end do
!count Q-non-Q
do j = 1, nat_solute
        if(iqatom(j) > 0) cycle
        if(any(qconn(:,j,:) <= 3)) then 
                !bonded or angled to at least one Q-atom in any state
                do iq = 1, nqat
                        do is = 1, nstates
                                if(qconn(is, j, iq) >= 4) then
                                        nbqq_pair(is) = nbqq_pair(is)+1
                                end if
                        end do
                end do
        end if
end do

nbqq_count = maxval(nbqq_pair(:))
end function nbqq_count

!---------------------------------------------------------------------------------

subroutine nbqqlist
integer						::	iq, j, jq, is, i, k,l, ia , ja
real(kind=prec)                     :: el_scale
logical                     :: set

nbqq_pair(:) = 0

nbqqp_pair(:) = 0
!list Q-Q
do iq = 1, nqat - 1
        ia = iqseq(iq)
        do jq = iq + 1, nqat
                ja = iqseq(jq)
                do is = 1, nstates
                        if(.not.qq_precomp(iq,jq,is)%set) cycle
                        nbqq_pair(is) = nbqq_pair(is)+1
                        nbqq(nbqq_pair(is),is)%iq = iq
                        nbqq(nbqq_pair(is),is)%jq = jq
                        nbqq(nbqq_pair(is),is)%vdWA = qq_precomp(iq,jq,is)%vdWA
                        nbqq(nbqq_pair(is),is)%vdWB = qq_precomp(iq,jq,is)%vdWB
                        nbqq(nbqq_pair(is),is)%elec = qq_precomp(iq,jq,is)%elec
                        nbqq(nbqq_pair(is),is)%score = qq_precomp(iq,jq,is)%score
                        nbqq(nbqq_pair(is),is)%soft = qq_precomp(iq,jq,is)%soft
                end do
        end do
end do

!list variable Q-P atoms
do ja = 1, nat_solute
        if(iqatom(ja) .ne. 0) cycle
        if(any(qconn(:,ja,:) <= 3)) then
                !bonded or angled to at least one Q-atom
                do iq = 1, nqat
                        do is = 1, nstates
                                if(.not.qp_precomp(ja,iq,is)%set) cycle
! the list style for the qp interactions has now been changed to be state
! dependent, so we can use a single list to calculate all qp interactions
! as this one here will be static, the minimum number of qp pairs will be kept
! the same always and is saved under the new variable nbqqp_pair(:)
! the lists will now be allocated based on the maximum size needed that is
! determined after running the precomputation
                                nbqqp_pair(is) = nbqqp_pair(is)+1
                                nbqqp(nbqqp_pair(is),is)%i = iq
                                nbqqp(nbqqp_pair(is),is)%j = ja
                                nbqqp(nbqqp_pair(is),is)%vdWA = qp_precomp(ja,iq,is)%vdWA
                                nbqqp(nbqqp_pair(is),is)%vdWB = qp_precomp(ja,iq,is)%vdWB
                                nbqqp(nbqqp_pair(is),is)%elec = qp_precomp(ja,iq,is)%elec
                                nbqqp(nbqqp_pair(is),is)%score = qp_precomp(ja,iq,is)%score
                        end do
                end do
        end if
end do

end subroutine nbqqlist

!-----------------------------------------------------------------------
!******PWchanged 2001-10-01
subroutine nbqp_count(nqp, nqpcgp,rcut)
! arguments
integer						:: nqp
integer						:: nqpcgp(:)
real(kind=prec)                                 :: rcut

! local variables
integer						:: ig,ia,i,j,iq
real(kind=prec)					:: rcut2,r2
TYPE(qr_vec)					:: shift

!	This routine counts non-bonded atom pairs involving
!	*one* Q-atom and *one* non-Q-atom, where the latter is *not connected*
!	(meaning not bonded or angled) to any Q-atom.
!
!	( i , j ) pairs correspond to
!	( iq, j ) with first index being a Q-atom with the *Q-atom numbering*, &
!	          and the second index is the non-Q-atom.

!	Note that for PBC the default is to count all atoms as interacting with each other
!	if the user did not specify a q - q cutoff of > 0

nqp = 0
rcut2 = Rcut**2



if(nqat .eq. 0) return

! --- solute - Q-atoms

igloop: do ig = 1, ncgp_solute

nqpcgp(ig) = 0

! skip if excluded group
ia = cgp(ig)%iswitch
if ( .not. use_PBC .and. excl(ia) ) cycle igloop

!******PWadded if 2001-10-01

if( .not. use_PBC ) then
	r2 = q_dist4(x(ia),xpcent)
else
	if (Rcq .gt. zero ) then
		shift = x(ia) - x(qswitch)
                r2 = q_dist4(shift,boxlength*q_nint(shift*inv_boxl))
	end if
end if

! skip if outside cutoff
if ( ( ( use_PBC ) .and. (r2 .gt. rcut2) .and. (Rcq .gt. zero ) ) .or. & 
	( ( r2 .gt. rcut2 ) ) ) cycle igloop

ialoop: do ia = cgp(ig)%first, cgp(ig)%last
  i = cgpatom(ia)

  if(iqatom(i) .ne. 0) cycle
  ! check if already on qq list
  if(any(qconn(:,i,:) <= 3)) cycle ialoop

  ! count the pairs
  nqp = nqp + nqat
  nqpcgp(ig) = nqpcgp(ig) + nqat

end do ialoop
end do igloop

end subroutine nbqp_count

!-----------------------------------------------------------------------
!******PWchanged 2001-10-01
subroutine nbqw_count(nqw, nqwmol,rcut)
! arguments
integer						:: nqw
integer						:: nqwmol(:)
real(kind=prec)                                 :: rcut

! local variables
integer						:: ig,ia,i,j,iq
real(kind=prec)					:: rcut2,r2
TYPE(qr_vec)					:: shift

!	This routine counts water molecules that interact with q-atoms
!	Note that in PBC the default is no cutoff (rcq < 0), so we just count all of them
nqw = 0
rcut2 = Rcut**2


if(nqat .eq. 0) return

! --- solvent - Q-atoms

iwloop: do ig = 1, nwat
nqwmol(ig) = 0
ia = nat_solute + solv_atom*ig-(solv_atom-1)
if(.not. use_PBC .and. excl(ia)) cycle iwloop ! skip excluded waters

!******PWadded if-statement 2001-10-01
if( .not. use_PBC ) then
	r2 = q_dist4(x(ia),xpcent)
else
	if (Rcq .gt. zero) then
		shift = x(ia) - x(qswitch)
                r2 = q_dist4(shift,boxlength*q_nint(shift*inv_boxl))
	end if
end if

! skip if outside cutoff
if ( (r2 .lt. rcut2 ) .or. (use_PBC.and.(rcq.lt.zero))) then
        nqw = nqw + 1
        nqwmol(ig) = solv_atom*nqat
end if
end do iwloop

end subroutine nbqw_count

!-----------------------------------------------------------------------

subroutine nbqplis2(Rcut)
! args
real(kind=prec)					:: Rcut
! local variables
integer						:: ig,ia,i,j,iq,nl,inside,is
real(kind=prec)						:: rcut2,r2
integer						:: xspec
logical, save					:: list_done
#ifdef _OPENMP
integer :: quotient, remainder
#endif

!	This routine makes a list of non-bonded atom pairs involving
!	*one* Q-atom and *one* non-Q-atom, where the latter is *not connected*
!	(meaning not bonded or angled) to any Q-atom.
!
!	( i , j ) pairs correspond to
!	( iq, j ) with first index being a Q-atom with the *Q-atom numbering*, &
!	          and the second index is the non-Q-atom.

! uses the global variables:
!  nbqs_pair, Rcq, cgp, excl, cgpatom, x, xpcent, nqat, iqseq, 
!  qconn, calculation_assignment%qs%max, nbqs, ljcod, nwat, nat_solute

!don't remake pair list if q-atoms interact with all solute atoms (rcq>rexcl_o)
!and list already made
#if defined (PROFILING)
real(kind=prec)                                         :: start_loop_time
start_loop_time = rtime()
#endif



! we don't remake the q - p list in PBC either now, by just including all freaking atoms into it
! if the q-q cutoff is set to < 0
! this is made default for the standart inputs, you can still set it to something else if you want to
if(list_done .and. (Rcut > rexcl_o**2)) return
!$omp single
nbqp_pair = 0
!$omp end single
!$omp barrier
rcut2 = Rcut


if(nqat .eq. 0) return


! --- solute - Q-atoms

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%qp%end - calculation_assignment%qp%start + 1)/threads_num
remainder = MOD(calculation_assignment%qp%end - calculation_assignment%qp%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%qp%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif
igloop: do ig = mp_start,mp_end
#else
igloop: do ig = calculation_assignment%qp%start, calculation_assignment%qp%end
#endif
! for every assigned charge group:

! skip if excluded group
        ia = cgp(ig)%iswitch
        if ( excl(ia) ) cycle igloop

! check cutoff
        inside = 0
        ia = cgp(ig)%first
        do while ((ia .le. cgp(ig)%last) .and. (inside .eq. 0))
                i = cgpatom(ia)

		r2 = q_dist4(x(i),xpcent)
                if ( r2 .le. rcut2 ) then
                        inside = 1
                end if
                ia = ia + 1
        end do
        if (inside .eq. 0) cycle igloop

ialoop: do ia = cgp(ig)%first, cgp(ig)%last
                i = cgpatom(ia)
! see if it is already on the list after nbqq
                if(any(qconn(:,i,:) <= 3)) cycle ialoop
!$omp critical
                ! if out of space then make more space
                if (nbqp_pair .ge. calculation_assignment%qp%max-nqat) call reallocate_nonbondlist_qp
qaloop:         do iq = 1, nqat
                        nbqp_pair = nbqp_pair + 1
                        do is = 1, nstates
! store the pair
                                nbqp(nbqp_pair,is)%i = iq
                                nbqp(nbqp_pair,is)%j = i
                                nbqp(nbqp_pair,is)%vdWA = qp_precomp(i,iq,is)%vdWA
                                nbqp(nbqp_pair,is)%vdWB = qp_precomp(i,iq,is)%vdWB
                                nbqp(nbqp_pair,is)%elec = qp_precomp(i,iq,is)%elec
                                nbqp(nbqp_pair,is)%score = qp_precomp(i,iq,is)%score
                        end do
                end do qaloop
!$omp end critical
        end do ialoop
end do igloop
!$omp single
list_done = .true. !save this value
!$omp end single
#if defined (PROFILING)
profile(5)%time = profile(5)%time + rtime() - start_loop_time
#endif

end subroutine nbqplis2



!----------------------------------------------------------------------------

subroutine nbqplis2_box(Rcut,Rq)
! args
real(kind=prec)					:: Rcut,Rq
!  ! local variables
integer						:: ig,ia,i,j,iq,nl,inside,ig_atom,is
real(kind=prec)					:: rcut2,r2
integer						:: xspec
logical,save					:: list_done = .false.
TYPE(qr_vec)					:: shift

! for periodic boundary conditions
!	This routine makes a list of non-bonded atom pairs involving
!	*one* Q-atom and *one* non-Q-atom, where the latter is *not connected*
!	(meaning not bonded or angled) to any Q-atom.
!
!	( i , j ) pairs correspond to
!	( iq, j ) with first index being a Q-atom with the *Q-atom numbering*, &
!	          and the second index is the non-Q-atom.

! uses the global variables:
!  nbqs_pair, Rcq, cgp, excl, cgpatom, x, xpcent, nqat, iqseq, 
!  qconn, calculation_assignment%qs%max, nbqs, ljcod, nwat, nat_solute

! now with more parallel
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)                                         :: start_loop_time
start_loop_time = rtime()
#endif

! we don't remake the q - p list in PBC either now, by just including all freaking atoms into it
! if the q-q cutoff is set to < 0
! this is made default for the standart inputs, you can still set it to something else if you want to
if(list_done .and. (((use_PBC).and.(Rq.lt.zero)))) return
if(nqat .eq. 0) return

!$omp single
nbqp_pair = 0
nbqp_cgp_pair = 0
!$omp end single
!$omp barrier
rcut2 = Rcut

 ! --- solute - Q-atoms
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%qp%end - calculation_assignment%qp%start + 1)/threads_num
remainder = MOD(calculation_assignment%qp%end - calculation_assignment%qp%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%qp%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

igloop: do ig = mp_start,mp_end
#else
igloop: do ig = calculation_assignment%qp%start, calculation_assignment%qp%end
#endif
    ! for every assigned charge group:

	! check cutoff
	! if Rcq < 0 -> no cutoff, skip distance calculation here
	if (Rq.lt.zero) then
	        inside = 1
	else
	        inside = 0
	end if
	ig_atom = cgp(ig)%first
	do while ((ig_atom .le. cgp(ig)%last) .and. (inside .eq. 0))
                i = cgpatom(ig_atom)
                shift = x(i) - x(qswitch)
                r2 = q_dist4(shift,boxlength*q_nint(shift*inv_boxl))
                if ( r2 .le. rcut2 ) then
		        inside = 1
!$omp critical
                        if(nbqp_cgp_pair .eq. size(nbqp_cgp, 1) ) call reallocate_nbqp_cgp
                        nbqp_cgp_pair = nbqp_cgp_pair + 1
                        nbqp_cgp(nbqp_cgp_pair)%i = i !leave %j empty, equals qswitch
!$omp end critical
                end if

                ig_atom = ig_atom + 1 !ia = ia + 1
        end do
	if (inside .eq. 0) cycle igloop

ialoop: do ia = cgp(ig)%first, cgp(ig)%last
                i = cgpatom(ia)

! check if already on qq list
                if(any(qconn(:,i,:) <= 3)) cycle ialoop

!$omp critical
                ! if out of space then make more space
                if (nbqp_pair .ge. calculation_assignment%qp%max-nqat) call reallocate_nonbondlist_qp
qaloop:         do iq = 1, nqat
                        nbqp_pair = nbqp_pair + 1
                        do is = 1, nstates
! store the pair
                                nbqp(nbqp_pair,is)%i = iq
                                nbqp(nbqp_pair,is)%j = i
                                nbqp(nbqp_pair,is)%vdWA = qp_precomp(i,iq,is)%vdWA
                                nbqp(nbqp_pair,is)%vdWB = qp_precomp(i,iq,is)%vdWB
                                nbqp(nbqp_pair,is)%elec = qp_precomp(i,iq,is)%elec
                                nbqp(nbqp_pair,is)%score = qp_precomp(i,iq,is)%score
                                nbqp(nbqp_pair,is)%cgp_pair = nbqp_cgp_pair
                        end do
                end do qaloop
!$omp end critical
        end do ialoop
end do igloop
#if defined (PROFILING)
profile(5)%time = profile(5)%time + rtime() - start_loop_time
#endif
!$omp single
list_done = .true. !save this value
!$omp end single

end subroutine nbqplis2_box
!-----------------------------------------------------------------------

subroutine nbqplist(Rcut)
! args
real(kind=prec)					:: Rcut
! local variables
integer						:: ig,ia,i,j,iq,nl,is
real(kind=prec)						:: rcut2,r2
integer						::	xspec
logical, save				::	list_done = .false.
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)						:: start_loop_time
start_loop_time = rtime()
#endif


!For spherical boundary	
!	This routine makes a list of non-bonded atom pairs involving
!	*one* Q-atom and *one* non-Q-atom, where the latter is *not connected*
!	(meaning not bonded or angled) to any Q-atom.
!
!	( i , j ) pairs correspond to
!	( iq, j ) with first index being a Q-atom with the *Q-atom numbering*, &
!	          and the second index is the non-Q-atom.

! global variables used:
!  nbqs_pair, Rcq, cgp, excl, x, xpcent, cgpatom, nqat, iqseq,
!  calculation_assignment%qp%max, nbqs, ljcod, nwat, nat_solute

!don't remake pair list if q-atoms interact with all solute atoms (rcq>rexcl_o)
!and list already made
if(list_done .and. (Rcut .gt. rexcl_o**2)) return
if(nqat .eq. 0) return
!$omp single
nbqp_pair = 0
!$omp end single
!$omp barrier
rcut2 = Rcut


! --- solute - Q-atoms

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%qp%end - calculation_assignment%qp%start + 1)/threads_num
remainder = MOD(calculation_assignment%qp%end - calculation_assignment%qp%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%qp%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif
igloop: do ig = mp_start,mp_end
#else
igloop: do ig = calculation_assignment%qp%start, calculation_assignment%qp%end
#endif

! skip if excluded group
        ia = cgp(ig)%iswitch
        if ( excl(ia) ) cycle igloop
	r2 = q_dist4(x(ia),xpcent)

! skip if outside cutoff
        if ( r2 .gt. rcut2 ) cycle igloop

ialoop: do ia = cgp(ig)%first, cgp(ig)%last

                i = cgpatom(ia)

  ! check if already on qq list
                if(any(qconn(:,i,:) <= 3)) cycle ialoop
! if out of space then make more space
!$omp critical
                if (nbqp_pair .ge. calculation_assignment%qp%max-nqat) call reallocate_nonbondlist_qp
qaloop:         do iq = 1, nqat
                        nbqp_pair = nbqp_pair + 1
                        do is = 1, nstates
! store the pair
                                nbqp(nbqp_pair,is)%i = iq
                                nbqp(nbqp_pair,is)%j = i
                                nbqp(nbqp_pair,is)%vdWA = qp_precomp(i,iq,is)%vdWA
                                nbqp(nbqp_pair,is)%vdWB = qp_precomp(i,iq,is)%vdWB
                                nbqp(nbqp_pair,is)%elec = qp_precomp(i,iq,is)%elec
                                nbqp(nbqp_pair,is)%score = qp_precomp(i,iq,is)%score
                        end do
                end do qaloop
!$omp end critical
        end do ialoop
end do igloop
#if defined (PROFILING)
profile(5)%time = profile(5)%time + rtime() - start_loop_time
#endif
!$omp single
list_done = .true. !save this value
!$omp end single
end subroutine nbqplist

!-----------------------------------------------------------------------

!******PWadded 2001-10-18
subroutine nbqplist_box(Rcut,Rq)
! args
real(kind=prec)					:: Rq,Rcut
! local variables
integer						:: ig,ia,i,j,iq,nl,ig_atom,is
real(kind=prec)					:: rcut2,r2
integer						::	xspec
integer						:: ga, gb, inside
logical,save					:: list_done = .false.
TYPE(qr_vec)					:: shift
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)                                         :: start_loop_time
start_loop_time = rtime()
#endif

!For periodic boundary conditions	
!	This routine makes a list of non-bonded atom pairs involving
!	*one* Q-atom and *one* non-Q-atom, where the latter is *not connected*
!	(meaning not bonded or angled) to any Q-atom.
!
!	( i , j ) pairs correspond to
!	( iq, j ) with first index being a Q-atom with the *Q-atom numbering*, &
!	          and the second index is the non-Q-atom.

! global variables used:
!  nbqs_pair, Rcq, cgp, x, xpcent, cgpatom, nqat, iqseq,
!  calculation_assignment%qs%max, nbqs, ljcod, nwat, nat_solute

! we don't remake the q - p list in PBC either now, by just including all freaking atoms into it
! if the q-q cutoff is set to < 0
! this is made default for the standart inputs, you can still set it to something else if you want to
if(list_done .and. (((use_PBC).and.(Rq.lt.zero)))) return

if(nqat .eq. 0) return
!$omp single
nbqp_pair = 0
!$omp end single
!$omp barrier
rcut2 = Rcut
! --- solute - Q-atoms
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%qp%end - calculation_assignment%qp%start + 1)/threads_num
remainder = MOD(calculation_assignment%qp%end - calculation_assignment%qp%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%qp%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

igloop: do ig = mp_start,mp_end
#else
igloop: do ig = calculation_assignment%qp%start, calculation_assignment%qp%end
#endif
        if (Rcq.lt.zero) then
                inside = 1
        else
                inside = 0
        end if
        ig_atom = cgp(ig)%first
        do while ((ig_atom .le. cgp(ig)%last) .and. (inside .eq. 0))
                i = cgpatom(ig_atom)
                shift = x(i) - x(qswitch)
                r2 = q_dist4(shift,boxlength*q_nint(shift*inv_boxl))
                if ( r2 .le. rcut2 ) then
                        inside = 1

!$omp critical
                        if( nbqp_cgp_pair .eq. size(nbqp_cgp, 1) ) call reallocate_nbqp_cgp
                        nbqp_cgp_pair = nbqp_cgp_pair + 1
                        nbqp_cgp(nbqp_cgp_pair)%i = i !leave %j empty, equals qswitch
!$omp end critical
                end if
                ig_atom = ig_atom + 1 !ia = ia + 1
        end do
        if (inside .eq. 0) cycle igloop

ialoop: do ia = cgp(ig)%first, cgp(ig)%last

                i = cgpatom(ia)
! check if already on qq list
                if(any(qconn(:,i,:) <= 3)) cycle ialoop
! if out of space then make more space
!$omp critical
                ! if out of space then make more space
                if (nbqp_pair .ge. calculation_assignment%qp%max-nqat) call reallocate_nonbondlist_qp
qaloop:         do iq = 1, nqat
                        nbqp_pair = nbqp_pair + 1
                        do is = 1, nstates
! store the pair
                                nbqp(nbqp_pair,is)%i = iq
                                nbqp(nbqp_pair,is)%j = i
                                nbqp(nbqp_pair,is)%vdWA = qp_precomp(i,iq,is)%vdWA
                                nbqp(nbqp_pair,is)%vdWB = qp_precomp(i,iq,is)%vdWB
                                nbqp(nbqp_pair,is)%elec = qp_precomp(i,iq,is)%elec
                                nbqp(nbqp_pair,is)%score = qp_precomp(i,iq,is)%score
                        end do
                end do qaloop
!$omp end critical
        end do ialoop
end do igloop

#if defined (PROFILING)
profile(5)%time = profile(5)%time + rtime() - start_loop_time
#endif
!$omp single
list_done = .true.
!$omp end single

end subroutine nbqplist_box
!------------------------------------------------------------------------------------

subroutine nbqwlist(Rcut)
! args
real(kind=prec)					:: Rcut
! local variables
integer						:: ig,ia,i,is,iq,istate
real(kind=prec)						:: rcut2,r2
logical,save					:: list_done = .false.

! we do not know the number of molecules in the solvent, so we need to get them first before allocating
! this is done in prep_sim and stored in solv_atom

!	This routine makes a list of water molecules within rcq from xpcent,
! i.e. the q-atom - water non-bond lists which implicitly includes all
! q-atoms with all atoms of the listed water
! waters may not have bonded interactions with q-atoms !
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)                                         :: start_loop_time
start_loop_time = rtime()
#endif



!We don't have to remake the list if q-atoms interact with all waters
!(rcq > rexcl_o) and list already made (nbwq_pair >0)
if( list_done .and.  (Rcut .gt. rexcl_o**2)) return

if(nqat .eq. 0) return
!$omp single
nbqw_pair = 0
!$omp end single
!$omp barrier
rcut2 = Rcut



#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%qw%end - calculation_assignment%qw%start + 1)/threads_num
remainder = MOD(calculation_assignment%qw%end - calculation_assignment%qw%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%qw%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

iwloop: do ig = mp_start,mp_end
#else
iwloop: do ig = calculation_assignment%qw%start, calculation_assignment%qw%end
#endif
        ia = nat_solute + solv_atom*ig-(solv_atom-1)
        if( excl(ia) ) cycle iwloop ! skip excluded waters
	r2 = q_dist4(x(ia),xpcent)

! store if inside cutoff
        if ( r2 <= rcut2 ) then
!$omp critical
qaloop:         do iq = 1, nqat
solvloop:               do is = 1, solv_atom
                                nbqw_pair = nbqw_pair + 1
                                nbqw(nbqw_pair,1:nstates)%i=iq
                                nbqw(nbqw_pair,1:nstates)%j=ia +(is-1)
                                nbqw(nbqw_pair,1:nstates)%elec  = qw_precomp(iq,is,1:nstates)%elec 
                                nbqw(nbqw_pair,1:nstates)%vdWA  = qw_precomp(iq,is,1:nstates)%vdWA
                                nbqw(nbqw_pair,1:nstates)%vdWB  = qw_precomp(iq,is,1:nstates)%vdWB
                                nbqw(nbqw_pair,1:nstates)%score = qw_precomp(iq,is,1:nstates)%score
                        end do solvloop
                end do qaloop
!$omp end critical
        end if
end do iwloop
#if defined (PROFILING)
profile(6)%time = profile(6)%time + rtime() - start_loop_time
#endif
! save this to prevent remaking of the list
!$omp single
list_done = .true.
!$omp end single

end subroutine nbqwlist

!-----------------------------------------------------------------------

!******PWadded 2001-10-18
subroutine nbqwlist_box(Rcut,Rq)
! args
real(kind=prec)					:: Rq,Rcut
! local variables
integer						:: ig,ia,i,is,iq
real(kind=prec)					:: rcut2,r2
logical,save					:: list_done = .false.
TYPE(qr_vec)					:: shift
! we do not know the number of molecules in the solvent, so we need to get them first before allocating
! this is done in prep_sim and stored in solv_atom
!	This routine makes a list of water molecules within rcq from xswitch,
! i.e. the q-atom - water non-bond lists which implicitly includes all
! q-atoms with all atoms of the listed water
! waters may not have bonded interactions with q-atoms !
!now with more parallel
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)                                         :: start_loop_time
start_loop_time = rtime()
#endif

if(list_done .and. (Rq.lt.zero)) return
if(nqat .eq. 0) return
!$omp single
nbqw_pair = 0
!$omp end single
!$omp barrier
rcut2 = Rcut

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%qw%end - calculation_assignment%qw%start + 1)/threads_num
remainder = MOD(calculation_assignment%qw%end - calculation_assignment%qw%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%qw%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

iwloop: do ig = mp_start,mp_end
#else
iwloop: do ig = calculation_assignment%qw%start, calculation_assignment%qw%end
#endif
	ia = nat_solute + solv_atom*ig-(solv_atom-1)
	if (Rq.gt.zero) then
	shift = x(ia) - x(qswitch)
        r2 = q_dist4(shift,boxlength*q_nint(shift*inv_boxl))
	end if
	! store if inside cutoff
	if ( ( r2 .le. rcut2 ).or. (Rq.lt.zero)) then
!$omp critical
qaloop:	        do iq = 1, nqat
solvloop:	        do is = 1, solv_atom
                                nbqw_pair = nbqw_pair + 1
                                nbqw(nbqw_pair,:)%i = iq
                                nbqw(nbqw_pair,:)%j = ia + is - 1
                                nbqw(nbqw_pair,:)%elec  = qw_precomp(iq,is,:)%elec 
                                nbqw(nbqw_pair,:)%vdWA  = qw_precomp(iq,is,:)%vdWA
                                nbqw(nbqw_pair,:)%vdWB  = qw_precomp(iq,is,:)%vdWB
                                nbqw(nbqw_pair,:)%score = qw_precomp(iq,is,:)%score
		        end do solvloop
		end do qaloop
!$omp end critical
        end if
end do iwloop
#if defined (PROFILING)
profile(6)%time = profile(6)%time + rtime() - start_loop_time
#endif
!$omp single
list_done = .true.
!$omp end single
end subroutine nbqwlist_box
!-------------------------------------------------------------------------------------
!******PWchanged 2001-10-01
subroutine nbww_count(nww, nwwmol,rcut)
! arguments
integer						:: nww
integer						:: nwwmol(:)
real(kind=prec)                                 :: rcut

! local variables
integer						:: iw,jw,ia,ja
real(kind=prec)					:: rcut2,r2
TYPE(qr_vec)					:: shift

! This routine counts non-bonded solvent-solvent atom pairs.

nww = 0
rcut2 = Rcut**2

iwloop: do iw = 1, nwat
nwwmol(iw) = 0

ia = nat_solute + solv_atom*iw-(solv_atom-1)
if(.not. use_PBC .and. excl(ia)) cycle iwloop ! skip excluded waters


jwloop: do jw = 1, nwat
  ja = nat_solute + solv_atom*jw-(solv_atom-1)
  if(.not. use_PBC .and. excl(ja)) cycle jwloop ! skip excluded waters


  ! count each w-w pair once only
  if ( ((iw .gt. jw) .and. (mod(iw+jw,2) .eq. 0)) .or. &
           ((iw .lt. jw) .and. (mod(iw+jw,2) .eq. 1)) .or. &
           (iw .eq. jw)) &
           cycle jwloop

  if( use_PBC ) then
  	shift = x(ia) - x(ja)
        r2 = q_dist4(shift,boxlength*q_nint(shift*inv_boxl))
  else
        r2 = q_dist4(x(ia),x(ja))
  end if
  ! count the pair if inside cutoff
  if ( r2 .le. rcut2 ) then
        nww = nww + 1
        nwwmol(iw) = nwwmol(iw) + solv_atom**2
  end if

end do jwloop
end do iwloop
end subroutine nbww_count

!-----------------------------------------------------------------------
subroutine nbwwlist(Rcut)
! args
real(kind=prec)					:: Rcut
! local variables
integer						:: iw,jw,ia,ja,la,ka,jwr
real(kind=prec)						:: rcut2,r2
integer                                                 ::iagrid, igrid, jgrid, kgrid, gridnum
! This routine makes a list of non-bonded solvent-solvent atom pairs
! The number of atoms in each solvent molecule is stored in the global solv_atom
! uses the global variables:
!  nbww_pair, Rcww, nat_solute, excl, nwat, x, nbww, calculation_assignment%ww%max
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)                                         :: start_loop_time
start_loop_time = rtime()
#endif

!$omp single 
nbww_pair = 0
!$omp end single
!$omp barrier
rcut2 = Rcut
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%ww%end - calculation_assignment%ww%start + 1)/threads_num
remainder = MOD(calculation_assignment%ww%end - calculation_assignment%ww%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%ww%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

iwloop: do iw = mp_start,mp_end
#else
iwloop: do iw = calculation_assignment%ww%start, calculation_assignment%ww%end
#endif

	ia = nat_solute + solv_atom*iw-(solv_atom-1)
        if (excl(ia)) cycle iwloop
#ifdef USE_GRID
        iagrid = ww_igrid(iw)
	gridnum = 0
igridloop:	do igrid = 1, ww_ndim
jgridloop:	 do jgrid = 1, ww_ndim
kgridloop:	  do kgrid = 1, ww_ndim
			gridnum = gridnum + 1
			if (.not. (grid_ww_int(iagrid,igrid,jgrid,kgrid))) cycle kgridloop
jwloop: do jw=1, grid_ww_ngrp(gridnum)
	        jwr = grid_ww_grp(gridnum,jw)
#else
jwloop:	do jwr =1, nwat
#endif
                ja = nat_solute + solv_atom*jwr-(solv_atom-1)
                if (excl(ja)) cycle jwloop ! skip excluded waters
! count each w-w pair once only
                if ( ((iw .gt. jwr) .and. (mod(iw+jwr,2) .eq. 0)) .or. &
                     ((iw .lt. jwr) .and. (mod(iw+jwr,2) .eq. 1)) .or. (iw .eq. jwr))  cycle jwloop

		r2 = q_dist4(x(ia),x(ja))

                if ( r2 .le. rcut2 ) then
! inside cutoff: add the pair

! check if there is enough space
!$omp critical
			if (nbww_pair .ge. (calculation_assignment%ww%max - solv_atom**2)) call reallocate_nonbondlist_ww
laloop:			do la = 1, solv_atom
kaloop:			        do ka = 1, solv_atom
                                        nbww_pair = nbww_pair + 1
                                        nbww(nbww_pair)%i = nat_solute +solv_atom*iw-(solv_atom-la)
        				nbww(nbww_pair)%j = nat_solute +solv_atom*jwr-(solv_atom-ka)
                                        nbww(nbww_pair)%elec = ww_precomp(la,ka)%elec
                                        nbww(nbww_pair)%vdWA = ww_precomp(la,ka)%vdWA
                                        nbww(nbww_pair)%vdWB = ww_precomp(la,ka)%vdWB
			        end do kaloop
			end do laloop
!$omp end critical
                end if
        end do jwloop
#ifdef USE_GRID
end do kgridloop
end do jgridloop
end do igridloop
#endif
end do iwloop
#if defined (PROFILING)
profile(2)%time = profile(2)%time + rtime() - start_loop_time
#endif
end subroutine nbwwlist
!------------PWadded 2001-10-18------------------------------------------
subroutine nbwwlist_box(Rcut)
! args
real(kind=prec)					:: Rcut
! local variables
integer						:: iw,jw,ia,ja,la,ka,jwr
real(kind=prec)					:: rcut2,r2
integer                                         :: iagrid, igrid, jgrid, kgrid, gridnum
TYPE(qr_vec)					:: shift
! for periodic boundary conditions
! This routine makes a list of non-bonded solvent-solvent atom pairs
! uses the global variables:
!  nbww_pair, Rcww, nat_solute, excl, nwat, x, nbww, calculation_assignment%ww%max

! now with more parallel
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)                                         :: start_loop_time
start_loop_time = rtime()
#endif


!$omp single
nbww_pair = 0
!$omp end single
!$omp barrier
rcut2 = Rcut

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%ww%end - calculation_assignment%ww%start + 1)/threads_num
remainder = MOD(calculation_assignment%ww%end - calculation_assignment%ww%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%ww%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif
iwloop: do iw = mp_start,mp_end
#else
iwloop: do iw = calculation_assignment%ww%start, calculation_assignment%ww%end
#endif
        ia = nat_solute + solv_atom*iw-(solv_atom-1)
#ifdef USE_GRID
        iagrid = ww_igrid(iw)
	gridnum = 0
igridloop:	do igrid = 1, ww_ndim
jgridloop:	 do jgrid = 1, ww_ndim
kgridloop:	  do kgrid = 1, ww_ndim
			gridnum = gridnum + 1
			if (.not. (grid_ww_int(iagrid,igrid,jgrid,kgrid))) cycle kgridloop
jwloop: do jw=1, grid_ww_ngrp(gridnum)
                jwr = grid_ww_grp(gridnum,jw)
#else
jwloop:	do jwr= 1, nwat
#endif
                ja = nat_solute + solv_atom*jwr-(solv_atom-1)
! count each w-w pair once only
                if ( ((iw .gt. jwr) .and. (mod(iw+jwr,2) .eq. 0)) .or. &
                        ((iw .lt. jwr) .and. (mod(iw+jwr,2) .eq. 1)) .or. &
                        (iw .eq. jwr)) &
                        cycle jwloop

                shift = x(ia) - x(ja)
                r2 = q_dist4(shift,boxlength*q_nint(shift*inv_boxl))
                if ( r2 .le. rcut2 ) then
! inside cutoff: add the pair
! check if there is enough space
!$omp critical
                        if (nbww_pair .ge. (calculation_assignment%ww%max - solv_atom**2)) call reallocate_nonbondlist_ww
laloop:                 do la = 1, solv_atom
kaloop:                         do ka = 1, solv_atom
                                        nbww_pair = nbww_pair + 1
                                        nbww(nbww_pair)%i = nat_solute+solv_atom*iw-(solv_atom-la)
                                        nbww(nbww_pair)%j = nat_solute+solv_atom*jwr-(solv_atom-ka)
                                        nbww(nbww_pair)%elec = ww_precomp(la,ka)%elec
                                        nbww(nbww_pair)%vdWA = ww_precomp(la,ka)%vdWA
                                        nbww(nbww_pair)%vdWB = ww_precomp(la,ka)%vdWB
                                end do kaloop
                        end do laloop
!$omp end critical
                end if
        end do jwloop
#ifdef USE_GRID
end do kgridloop
end do jgridloop
end do igridloop
#endif
end do iwloop
#if defined (PROFILING)
profile(2)%time = profile(2)%time + rtime() - start_loop_time
#endif

end subroutine nbwwlist_box

!---------------------------------------------------------------------------

subroutine nbwwlist_lrf(Rcut,RLRF)
! args
real(kind=prec)					:: Rcut,RLRF
! local variables
integer						:: i,j,ig,jg,iw,jw,ia,ja,is,is3,la,ka,jwr
real(kind=prec)						:: rcut2,r2
real(kind=prec)						::	RcLRF2
integer                                                 ::iagrid, igrid, jgrid, kgrid, gridnum
! now with more parallel -> OMP added
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#if defined (PROFILING)
real(kind=prec)						:: start_loop_time
start_loop_time = rtime()
#endif
!	This routine makes a list of non-bonded solvent-solvent atom pairs.

! uses the global variables:
!  nbww_pair, Rcww, nwat, nat_solute, excl, x, nbww, ncgp, lrf, crg, calculation_assignment%ww%max



!$omp single
nbww_pair = 0
!$omp end single
!$omp barrier
rcut2 = Rcut
RcLRF2 = RLRF

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%ww%end - calculation_assignment%ww%start + 1)/threads_num
remainder = MOD(calculation_assignment%ww%end - calculation_assignment%ww%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%ww%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

iwloop: do iw = mp_start,mp_end
#else
iwloop: do iw = calculation_assignment%ww%start, calculation_assignment%ww%end
#endif
        is  = nat_solute + solv_atom*iw-(solv_atom-1)
        if( excl(is)) cycle iwloop
        ig = iwhich_cgp(is)
#ifdef USE_GRID
        iagrid = ww_igrid(iw)
        gridnum = 0
igridloop:	do igrid = 1, ww_ndim
jgridloop:	 do jgrid = 1, ww_ndim
kgridloop:	  do kgrid = 1, ww_ndim
			gridnum = gridnum + 1
			if (.not. (grid_ww_int(iagrid,igrid,jgrid,kgrid))) then
ggloop:			 do jw = 1, grid_ww_ngrp(gridnum)
			  jwr = grid_ww_grp(gridnum,jw)
			  ja = nat_solute + solv_atom*jwr-(solv_atom-1)
			  if(excl(ja)) cycle ggloop
			  jg = iwhich_cgp(ja)
                          if ( ((iw .gt. jwr) .and. (mod(iw+jwr,2) .eq. 0)) .or. &
                                  ((iw .lt. jwr) .and. (mod(iw+jwr,2) .eq. 1)).or. &
                                  (iw .eq. jwr) ) &
                                  cycle ggloop
			  call lrf_update(ig,jg)
			  call lrf_update(jg,ig)
			 end do ggloop
			else
jwloop:	do jw = 1, grid_ww_ngrp(gridnum)
	        jwr = grid_ww_grp(gridnum,jw)
#else
jwloop: do jwr = 1, nwat
#endif
                ja = nat_solute + solv_atom*jwr-(solv_atom-1)
                if(excl(ja)) cycle jwloop
                jg = iwhich_cgp(ja)
! count each w-w pair once only
                if ( ((iw .gt. jwr) .and. (mod(iw+jwr,2) .eq. 0)) .or. &
                        ((iw .lt. jwr) .and. (mod(iw+jwr,2) .eq. 1)) .or. &
                        (iw .eq. jwr) ) &
                        cycle jwloop

		r2 = q_dist4(x(is),x(ja))
                if ( r2 .le. rcut2 ) then
! inside cutoff: add the pair
! this now has to be OMP critical so the memory is not getting corrupted
!$omp critical 
! check if there is enough space
                        if (nbww_pair .ge. (calculation_assignment%ww%max - solv_atom**2)) call reallocate_nonbondlist_ww
laloop:                 do la = 1, solv_atom
kaloop:                         do ka = 1, solv_atom
                                        nbww_pair = nbww_pair + 1
                                        nbww(nbww_pair)%i = nat_solute+solv_atom*iw-(solv_atom-la)
                                        nbww(nbww_pair)%j = nat_solute+solv_atom*jwr-(solv_atom-ka)
                                        nbww(nbww_pair)%elec = ww_precomp(la,ka)%elec
                                        nbww(nbww_pair)%vdWA = ww_precomp(la,ka)%vdWA
                                        nbww(nbww_pair)%vdWB = ww_precomp(la,ka)%vdWB
                                end do kaloop
                        end do laloop
!$omp end critical

                elseif (r2 .le. RcLRF2) then
! outside ww-cutoff but inside LRF cut-off: use LRF
                        call lrf_update(ig,jg)
                        call lrf_update(jg,ig)
                end if
        end do jwloop
#ifdef USE_GRID
end if ! interaction
end do kgridloop
end do jgridloop
end do igridloop
#endif
end do iwloop
#if defined (PROFILING)
profile(2)%time = profile(2)%time + rtime() - start_loop_time
#endif

end subroutine nbwwlist_lrf
!--------------LRF version of PW PBC-----------------------------------
subroutine nbwwlist_box_lrf(Rcut,RLRF)
! args
real(kind=prec)					:: Rcut,RLRF
! local variables
integer						:: i,j,ig,jg,iw,jw,ia,ja,is,is3,la,ka,jwr
real(kind=prec)					:: rcut2,r2
real(kind=prec)					::	RcLRF2
integer                                         :: iagrid, igrid, jgrid, kgrid, gridnum
TYPE(qr_vec)					:: shift
#ifdef _OPENMP
integer :: quotient, remainder
#endif

#if defined (PROFILING)
real(kind=prec)                                         :: start_loop_time
start_loop_time = rtime()
#endif


!	This routine makes a list of non-bonded solvent-solvent atom pairs.

! uses the global variables:
!  nbww_pair, Rcww, nwat, nat_solute, excl, x, nbww, ncgp, lrf, crg, calculation_assignment%ww%max


!$omp single
nbww_pair = 0
!$omp end single
!$omp barrier
rcut2 = Rcut
RcLRF2 = RLRF

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%ww%end - calculation_assignment%ww%start + 1)/threads_num
remainder = MOD(calculation_assignment%ww%end - calculation_assignment%ww%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%ww%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

iwloop: do iw = mp_start,mp_end
#else
iwloop: do iw = calculation_assignment%ww%start, calculation_assignment%ww%end
#endif
        is  = nat_solute + solv_atom*iw-(solv_atom-1)
        ig = iwhich_cgp(is)
#ifdef USE_GRID
        iagrid = ww_igrid(iw)
        gridnum = 0
igridloop:      do igrid = 1, ww_ndim
jgridloop:       do jgrid = 1, ww_ndim
kgridloop:        do kgrid = 1, ww_ndim
                        gridnum = gridnum + 1
                        if (.not. (grid_ww_int(iagrid,igrid,jgrid,kgrid))) then
ggloop:                  do jw = 1, grid_ww_ngrp(gridnum)
                          jwr = grid_ww_grp(gridnum,jw)
                          ja = nat_solute + solv_atom*jwr-(solv_atom-1)
                          jg = iwhich_cgp(ja)
                          call lrf_update(ig,jg)
                          call lrf_update(jg,ig)
                         end do ggloop
                        else
jwloop: do jw = 1, grid_ww_ngrp(gridnum)
                jwr = grid_ww_grp(gridnum,jw)
#else
jwloop: do jwr = 1, nwat
#endif
                ja = nat_solute + solv_atom*jwr-(solv_atom-1)

! count each w-w pair once only
                if ( ((iw .gt. jwr) .and. (mod(iw+jwr,2) .eq. 0)) .or. &
                        ((iw .lt. jwr) .and. (mod(iw+jwr,2) .eq. 1)) .or. &
                        (iw .eq. jwr) ) &
                        cycle jwloop

                shift = x(is) - x(ja)
                r2 = q_dist4(shift,boxlength*q_nint(shift*inv_boxl))
                jg = iwhich_cgp(ja)
                if ( r2 .le. rcut2 ) then
! inside cutoff: add the pair
! check if there is enough space
!$omp critical
                        if (nbww_pair .ge. (calculation_assignment%ww%max - solv_atom**2)) call reallocate_nonbondlist_ww
laloop:                 do la = 1, solv_atom
kaloop:                         do ka = 1, solv_atom
                                        nbww_pair = nbww_pair + 1
                                        nbww(nbww_pair)%i = nat_solute+solv_atom*iw-(solv_atom-la)
                                        nbww(nbww_pair)%j = nat_solute+solv_atom*jwr-(solv_atom-ka)
                                        nbww(nbww_pair)%elec = ww_precomp(la,ka)%elec
                                        nbww(nbww_pair)%vdWA = ww_precomp(la,ka)%vdWA
                                        nbww(nbww_pair)%vdWB = ww_precomp(la,ka)%vdWB
                                end do kaloop
                        end do laloop
!$omp end critical

                elseif((r2 .le. RcLRF2) .or. (RcLRF .eq. -one)) then
! outside ww-cutoff but inside LRF cut-off: use LRF
!iw interaction     
                        call lrf_update(ig,jg)	
!jw interaction
                        call lrf_update(jg,ig)
                end if
        end do jwloop
#ifdef USE_GRID
end if ! interaction
end do kgridloop
end do jgridloop
end do igridloop
#endif
end do iwloop
#if defined (PROFILING)
profile(2)%time = profile(2)%time + rtime() - start_loop_time
#endif
end subroutine nbwwlist_box_lrf
!---------------------------------------------------------------------------
subroutine nbmonitorlist
!precalculate interactions for later use here, too
!old stuff was not usefull any more
! local variables
integer         :: i,j,ig,jg,ia,ja,nl,istate,LJ_code,maxingroups,par, atomnri
integer         :: max_int
integer,allocatable :: num_int(:)
integer         :: grpi,grpj,atomi,atomj,qq_pair,aLJ,bLJ


if (monitor_group_pairs == 0) return

!check the size of the largest group
maxingroups=maxval(monitor_atom_group(:)%n)
max_int = sum(monitor_atom_group(:)%n)

allocate(monitor_group_int(max_int,nstates))
monitor_group_int(:,:)%score = zero
monitor_group_int(:,:)%soft  = .false.
monitor_group_int(:,:)%iq    = -1
monitor_group_int(:,:)%jq    = -1
monitor_group_int(:,:)%elec  = -1E35_prec
monitor_group_int(:,:)%vdWA  = -1E35_prec
monitor_group_int(:,:)%vdWB  = -1E35_prec
allocate(num_int(nstates))
num_int(:) = 0

do par=1,monitor_group_pairs
        monitor_group_pair(par)%lstart(:) = num_int(:) + 1
        grpi=monitor_group_pair(par)%i      
        grpj=monitor_group_pair(par)%j 
        do i=1,monitor_atom_group(grpi)%n
                atomi=monitor_atom_group(grpi)%atom(i)
                do j=1,monitor_atom_group(grpj)%n
                        atomj=monitor_atom_group(grpj)%atom(j)    
!assert that atoms are different
                        if(atomi == atomj) then 
                                call die('two paired monitor atom groups contain the same atom')
                        end if
!check which list we need to check
                        if((iqatom(atomi).ne.0) .and. (iqatom(atomj).ne.0)) then
!both q_atoms, are on qq list
!check if those bastards are allowed to interact
!we now do the real checking to see if the atoms actually interact
!to make sure we get the interaction energy right 
                                do istate = 1 , nstates
                                if (.not.(qq_precomp(iqatom(atomi),iqatom(atomj),istate)%set)) cycle
                                        num_int(istate) = num_int(istate) + 1
                                        monitor_group_int(num_int(istate),istate)%iq= atomi
                                        monitor_group_int(num_int(istate),istate)%jq= atomj
                                        monitor_group_int(num_int(istate),istate)%elec = &
                                                qq_precomp(iqatom(atomi),iqatom(atomj),istate)%elec
                                        monitor_group_int(num_int(istate),istate)%vdWA = &
                                                qq_precomp(iqatom(atomi),iqatom(atomj),istate)%vdWA
                                        monitor_group_int(num_int(istate),istate)%vdWB = &
                                                qq_precomp(iqatom(atomi),iqatom(atomj),istate)%vdWB
! for q-q atoms we also need the info on soft pairs and soft core stuff
! so the group_int stuff needs to be able to store this, too
                                        monitor_group_int(num_int(istate),istate)%soft = &
                                                qq_precomp(iqatom(atomi),iqatom(atomj),istate)%soft
                                        monitor_group_int(num_int(istate),istate)%score = &
                                                qq_precomp(iqatom(atomi),iqatom(atomj),istate)%score
                                end do

                        else if(iqatom(atomi).ne.0) then
!only atom i is Q atom, atom j will be on q-p list then
!so we check the q-p precompute table for them, to also account for qqp atoms
                                do istate = 1 , nstates
                                if(.not.(qp_precomp(iqatom(atomi),atomj,istate)%set)) cycle
                                        num_int(istate) = num_int(istate) + 1
                                        monitor_group_int(num_int(istate),istate)%iq= atomi
                                        monitor_group_int(num_int(istate),istate)%jq= atomj
                                        monitor_group_int(num_int(istate),istate)%elec = &
                                                qp_precomp(iqatom(atomi),atomj,istate)%elec
                                        monitor_group_int(num_int(istate),istate)%vdWA = &
                                                qp_precomp(iqatom(atomi),atomj,istate)%vdWA
                                        monitor_group_int(num_int(istate),istate)%vdWB = &
                                                qp_precomp(iqatom(atomi),atomj,istate)%vdWB
! for q-q atoms we also need the info on soft pairs and soft core stuff
! so the group_int stuff needs to be able to store this, too
                                        monitor_group_int(num_int(istate),istate)%score = &
                                                qp_precomp(iqatom(atomi),atomj,istate)%score
                                end do
!other way around, j is Q atom, same code but inverted
                        else if(iqatom(atomj).ne.0) then
                                do istate = 1 , nstates
                                if(.not.(qp_precomp(iqatom(atomj),atomi,istate)%set)) cycle
                                        num_int(istate) = num_int(istate) + 1
                                        monitor_group_int(num_int(istate),istate)%iq= atomi
                                        monitor_group_int(num_int(istate),istate)%jq= atomj
                                        monitor_group_int(num_int(istate),istate)%elec = &
                                                qp_precomp(iqatom(atomj),atomi,istate)%elec
                                        monitor_group_int(num_int(istate),istate)%vdWA = &
                                                qp_precomp(iqatom(atomj),atomi,istate)%vdWA
                                        monitor_group_int(num_int(istate),istate)%vdWB = &
                                                qp_precomp(iqatom(atomj),atomi,istate)%vdWB
! for q-q atoms we also need the info on soft pairs and soft core stuff
! so the group_int stuff needs to be able to store this, too
                                        monitor_group_int(num_int(istate),istate)%score = &
                                                qp_precomp(iqatom(atomj),atomi,istate)%score
                                end do
! and if neither is a qatom, we use the pp-list
                        else
                                if(pp_map(i-pp_low,j).eq.0) cycle
                                if(.not.(pp_precomp(pp_map(i-pp_low,j))%set)) cycle
                                num_int(:) = num_int(:) + 1
                                monitor_group_int(num_int(:),:)%iq= atomi
                                monitor_group_int(num_int(:),:)%jq= atomj
                                monitor_group_int(num_int(:),:)%elec = &
                                        pp_precomp(pp_map(i-pp_low,j))%elec
                                monitor_group_int(num_int(:),:)%vdWA = &
                                        pp_precomp(pp_map(i-pp_low,j))%vdWA
                                monitor_group_int(num_int(:),:)%vdWB = &
                                        pp_precomp(pp_map(i-pp_low,j))%vdWB
                        end if
                end do
        end do
        monitor_group_pair(par)%lend(:) = num_int(:)
end do
! cleaning up again
deallocate(num_int)                                
end subroutine nbmonitorlist

!---------------------------------------------------------------------------------------------------------		   

subroutine nonbond_monitor
!monitor nonbonded energies between selected groups of atoms
!rewritten tu use precomputed interaction energies

real(kind=prec)  :: r,r2,r6,r12,r6_hc
integer                 :: istate,par
integer                 :: cstart,cend,calcgroup
real(kind=prec)         :: V_a,V_b,Vel,Vvdw,Vwel,Vwvdw,Vwsum
TYPE(ENERET_TYPE)       :: nb_ene

do par=1,monitor_group_pairs
        monitor_group_pair(par)%Vel(:)=zero
        monitor_group_pair(par)%Vlj(:)=zero
        monitor_group_pair(par)%Vwel = zero
        monitor_group_pair(par)%Vwlj = zero
        monitor_group_pair(par)%Vwsum= zero  

        do istate = 1, nstates   
! loop over all our nice pairs for each state
! defined under lstart and lend
! that then point to the right number in monitor_int
        cstart = monitor_group_pair(par)%lstart(istate)
        cend   = monitor_group_pair(par)%lend(istate)

                do calcgroup = cstart, cend
                        if(use_PBC) then
                                nb_ene = nbe_qqb(monitor_group_int(calcgroup,istate),EQ(istate)%lambda)
                        else
                                nb_ene = nbe_qq(monitor_group_int(calcgroup,istate),EQ(istate)%lambda)
                        end if
                        Vel  = nb_ene%Vel
                        Vvdw = nb_ene%V_a - nb_ene%V_b
                        monitor_group_pair(par)%Vel(istate) = &
                                monitor_group_pair(par)%Vel(istate)+Vel
                        monitor_group_pair(par)%Vlj(istate) = &
                                monitor_group_pair(par)%Vlj(istate)+Vvdw
                end do ! groups
                monitor_group_pair(par)%Vwel = monitor_group_pair(par)%Vwel + &
                        monitor_group_pair(par)%Vel(istate) * EQ(istate)%lambda
                monitor_group_pair(par)%Vwlj = monitor_group_pair(par)%Vwlj + &
                        monitor_group_pair(par)%Vlj(istate) * EQ(istate)%lambda
        end do ! nstates
        monitor_group_pair(par)%Vwsum= monitor_group_pair(par)%Vwlj+monitor_group_pair(par)%Vwel
end do !par

end subroutine nonbond_monitor

!-------------------------------------------------------------------------
! nonbond routines are now consolidated to just sphere/box versions
! because combination rules are applied earlier
! because it is its own module, they get passed the list they need


subroutine nonbond_pp(E_loc)
! arguments
TYPE(NB_ENERGIES)                  :: E_loc
! local variables
integer                         :: ip,i,j
real(kind=prec)                 :: Vel,V_a,V_b,dv
TYPE(ENERET_TYPE)               :: nb_ene
#ifdef _OPENMP
integer :: quotient, remainder
real(kind=prec) :: Vel_omp, Vvdw_omp
#endif

! global variables used:
! x,d,E


#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (nbpp_pair)/threads_num
remainder = MOD(nbpp_pair, threads_num)
mp_start = thread_id * quotient + 1 + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif
Vel_omp = zero
Vvdw_omp = zero
do ip = mp_start,mp_end
#else
do ip = 1, nbpp_pair 
#endif

! do this for every pair and use the compiler to unroll the loop
i      = nbpp(ip)%i
j      = nbpp(ip)%j
nb_ene = nbe(nbpp(ip))

! calculate Vel and dv
Vel  = nb_ene%Vel  
V_a  = nb_ene%V_a
V_b  = nb_ene%V_b
dv   = nb_ene%dv

! update d
d(i) = d(i) - nb_ene%vec*dv
d(j) = d(j) + nb_ene%vec*dv

! update energies
#ifdef _OPENMP
Vel_omp  = Vel_omp + Vel
Vvdw_omp = Vvdw_omp +  V_a - V_b
#else
E_loc%el  = E_loc%el + Vel 
E_loc%vdw = E_loc%vdw + V_a - V_b
#endif
end do

#ifdef _OPENMP
!$omp atomic update
E%el  = E%el  + Vel_omp
!$omp atomic update
E%vdw = E%vdw + Vvdw_omp
#endif


end subroutine nonbond_pp

!------------------------------------------------------------------------
subroutine nonbond_pp_box(E_loc)
! arguments
TYPE(NB_ENERGIES)                       :: E_loc
! local variables
integer                                 :: ip, ga, gb, group,i,j
real(kind=prec)                         :: Vel,V_a,V_b,dv
TYPE(qr_vec)                            :: shift
TYPE(ENERET_TYPE)                       :: nb_ene
#ifdef _OPENMP
integer :: quotient, remainder
real(kind=prec) :: Vel_omp, Vvdw_omp
#endif
! global variables used:
! x,d,E

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (nbpp_cgp_pair)/threads_num
remainder = MOD(nbpp_cgp_pair, threads_num)
mp_start = thread_id * quotient + 1 + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

do group = mp_start,mp_end
#else
do group = 1, nbpp_cgp_pair
#endif
ga = nbpp_cgp(group)%i !atom index for the two switching atoms
gb = nbpp_cgp(group)%j

!the distance between the two switching atoms
shift = x(ga) - x(gb)

nbpp_cgp(group)%shift = boxlength*q_nint(shift*inv_boxl)

end do
!$omp barrier
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (nbpp_pair)/threads_num
remainder = MOD(nbpp_pair, threads_num)
mp_start = thread_id * quotient + 1 + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif
Vel_omp = zero
Vvdw_omp = zero
do ip = mp_start,mp_end
#else
do ip = 1, nbpp_pair !- 1, 2
#endif
! for every pair, not every second pair

i      = nbpp(ip)%i
j      = nbpp(ip)%j
group  = nbpp(ip)%cgp_pair
nb_ene = nbe_b(nbpp(ip),nbpp_cgp(group)%shift)

Vel  = nb_ene%Vel  
V_a  = nb_ene%V_a  
V_b  = nb_ene%V_b  
dv   = nb_ene%dv  

! update d
d(i) = d(i) - nb_ene%vec*dv
d(j) = d(j) + nb_ene%vec*dv

! update energies
#ifdef _OPENMP
Vel_omp  = Vel_omp  + Vel
Vvdw_omp = Vvdw_omp + V_a - V_b
#else
E_loc%el  = E_loc%el  + Vel
E_loc%vdw = E_loc%vdw + V_a - V_b
#endif
end do

#ifdef _OPENMP
!$omp atomic update
E%el  = E%el  + Vel_omp
!$omp atomic update
E%vdw = E%vdw + Vvdw_omp
#endif


end subroutine nonbond_pp_box
!----------------------------------------------------------------------

subroutine nonbond_pw(E_loc)
! arguments
TYPE(NB_ENERGIES)                               :: E_loc
! local variables
integer                                         :: ip,i,j
real(kind=prec)                                 :: Vel,V_a,V_b,dv
TYPE(ENERET_TYPE)                               :: nb_ene
#ifdef _OPENMP
integer :: quotient, remainder
real(kind=prec) :: Vel_omp, Vvdw_omp
#endif
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (nbpw_pair)/threads_num
remainder = MOD(nbpw_pair, threads_num)
mp_start = thread_id * quotient + 1 + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif
Vel_omp = zero
Vvdw_omp = zero
do ip = mp_start,mp_end
#else
! global variables used:
! x, d, E

do ip = 1, nbpw_pair
#endif
! for every assigned pair:
i      = nbpw(ip)%i
j      = nbpw(ip)%j
nb_ene = nbe(nbpw(ip))

! calculate Vel and dv
Vel  = nb_ene%Vel 
V_a  = nb_ene%V_a 
V_b  = nb_ene%V_b 
dv   = nb_ene%dv 

! update forces
d(i) = d(i) - nb_ene%vec*dv
d(j) = d(j) + nb_ene%vec*dv

! update energies
#ifdef _OPENMP
Vel_omp  = Vel_omp  + Vel
Vvdw_omp = Vvdw_omp + V_a - V_b
#else
E_loc%el  = E_loc%el  + Vel       
E_loc%vdw = E_loc%vdw + V_a - V_b
#endif
end do
#ifdef _OPENMP
!$omp atomic update
E_loc%el  = E_loc%el  + Vel_omp
!$omp atomic update
E_loc%vdw = E_loc%vdw + Vvdw_omp
#endif
end subroutine nonbond_pw

!-----------------------------------------------------------------------

subroutine nonbond_pw_box(E_loc)
! arguments
TYPE(NB_ENERGIES)                               :: E_loc
! local variables
integer                                         :: ip,i,j
real(kind=prec)                                 :: Vel,V_a,V_b,dv
integer                                         :: group, ga, gb
TYPE(qr_vec)                                    :: shift
TYPE(ENERET_TYPE)                               :: nb_ene
#ifdef _OPENMP
integer :: quotient, remainder
real(kind=prec) :: Vel_omp, Vvdw_omp
#endif
  ! global variables used:
  !  iac, crg, iaclib, x, d, E

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (nbpw_cgp_pair)/threads_num
remainder = MOD(nbpw_cgp_pair, threads_num)
mp_start = thread_id * quotient + 1 + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
mp_end = mp_end + 1
endif

do group = mp_start,mp_end
#else
!compute the peridocal shift for every charge group pair
do group = 1, nbpw_cgp_pair
#endif
ga = nbpw_cgp(group)%i !atom index for solute switching atom
gb = nbpw_cgp(group)%j  !atom index for the solvent switching atom

!the distance between the two switching atoms
shift = x(ga) - x(gb)
nbpw_cgp(group)%shift = boxlength*q_nint(shift*inv_boxl)

end do
!$omp barrier
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (nbpw_pair)/threads_num
remainder = MOD(nbpw_pair, threads_num)
mp_start = thread_id * quotient + 1 + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
mp_end = mp_end + 1
endif
Vel_omp = zero
Vvdw_omp = zero
do ip = mp_start,mp_end
#else
do ip = 1, nbpw_pair
#endif
! for every assigned pair:
i      = nbpw(ip)%i
j      = nbpw(ip)%j
group  = nbpw(ip)%cgp_pair
nb_ene = nbe_b(nbpw(ip),nbpw_cgp(group)%shift)


! calculate Vel and dv
Vel  = nb_ene%Vel 
V_a  = nb_ene%V_a 
V_b  = nb_ene%V_b 
dv   = nb_ene%dv 
! update forces
d(i) = d(i) - nb_ene%vec*dv
d(j) = d(j) + nb_ene%vec*dv

! update energies
#ifdef _OPENMP
Vel_omp  = Vel_omp  + Vel
Vvdw_omp = Vvdw_omp + V_a - V_b
#else
E_loc%el  = E_loc%el  + Vel       
E_loc%vdw = E_loc%vdw + V_a - V_b
#endif
end do
#ifdef _OPENMP
!$omp atomic update
E_loc%el  = E_loc%el  + Vel_omp
!$omp atomic update
E_loc%vdw = E_loc%vdw + Vvdw_omp
#endif

end subroutine nonbond_pw_box

!-----------------------------------------------------------------------

subroutine nonbond_qq(EQ_loc,lambda)
! arguments
TYPE(NB_ENERGIES)                               :: EQ_loc(:)
real(kind=prec)                                 :: lambda(:)
! local variables
integer                                         :: ip,i,j,iq,jq,istate
real(kind=prec)                                 :: Vel,V_a,V_b,dv
TYPE(ENERET_TYPE)                               :: nb_ene
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#ifdef _OPENMP
qomp_elec  = zero
qomp_vdw   = zero
#endif
do istate = 1 , nstates
! for every state:
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (nbqq_pair(istate))/threads_num
remainder = MOD(nbqq_pair(istate), threads_num)
mp_start = thread_id * quotient + 1 + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
mp_end = mp_end + 1
endif

do ip = mp_start,mp_end
#else
do ip = 1, nbqq_pair(istate)
#endif
  ! for every pair:

iq = nbqq(ip,istate)%iq
jq = nbqq(ip,istate)%jq
i  = iqseq(iq)
j  = iqseq(jq)


nb_ene = nbe_qq(nbqq(ip,istate),lambda(istate))


! calculate Vel, V_a, V_b and dv
Vel = nb_ene%Vel 
V_b = nb_ene%V_b 
V_a = nb_ene%V_a 
dv  = nb_ene%dv 

! update forces
d(i) = d(i) - nb_ene%vec*dv
d(j) = d(j) + nb_ene%vec*dv
! update energies
#ifdef _OPENMP
qomp_elec(istate)  = qomp_elec(istate) + Vel
qomp_vdw(istate)   = qomp_vdw(istate)  + V_a - V_b
#else
EQ_loc(istate)%el  = EQ_loc(istate)%el  + Vel
EQ_loc(istate)%vdw = EQ_loc(istate)%vdw + V_a - V_b
#endif
end do ! ip

end do ! istate
#ifdef _OPENMP
!$omp atomic update
EQ_loc%el  = EQ_loc%el  + qomp_elec
!$omp atomic update
EQ_loc%vdw = EQ_loc%vdw + qomp_vdw
#endif
end subroutine nonbond_qq

!-----------------------------------------------------------------------
subroutine nonbond_qqp(EQ_loc,lambda)
! arguments
TYPE(NB_ENERGIES)                               :: EQ_loc(:)
real(kind=prec)                                 :: lambda(:)
! local variables
integer                                         :: ip,iq,i,j,istate
real(kind=prec)                                 :: Vel,V_a,V_b,dv
TYPE(ENERET_TYPE)                               :: nb_ene
TYPE(qr_dist3)                                  :: dist
#ifdef _OPENMP
integer :: quotient, remainder
#endif
#ifdef _OPENMP
qomp_elec  = zero
qomp_vdw   = zero
#endif
do istate = 1, nstates
! for every state:
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (nbqqp_pair(istate))/threads_num
remainder = MOD(nbqqp_pair(istate), threads_num)
mp_start = thread_id * quotient + 1 + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
mp_end = mp_end + 1
endif

do ip = mp_start,mp_end
#else
do ip = 1, nbqqp_pair(istate)
#endif
  ! for every pair:
iq = nbqqp(ip,istate)%i
i  = iqseq(iq)
j  = nbqqp(ip,istate)%j


dist   = q_dist3(x(i),x(j))
nb_ene = nbe_qx(nbqqp(ip,istate),lambda(istate),dist)

! calculate Vel, V_a, V_b and dv
Vel = nb_ene%Vel 
V_a = nb_ene%V_a 
V_b = nb_ene%V_b 
dv  = nb_ene%dv 

! update forces
d(i) = d(i) - dist%vec*dv
d(j) = d(j) + dist%vec*dv
! update energies
#ifdef _OPENMP
qomp_elec(istate) = qomp_elec(istate) + Vel
qomp_vdw(istate)  = qomp_vdw(istate)  + V_a - V_b
#else
EQ_loc(istate)%el  = EQ_loc(istate)%el  + Vel
EQ_loc(istate)%vdw = EQ_loc(istate)%vdw + V_a - V_b
#endif
end do ! ip

end do ! istate
#ifdef _OPENMP
!$omp atomic update
EQ_loc%el  = EQ_loc%el  + qomp_elec
!$omp atomic update
EQ_loc%vdw = EQ_loc%vdw + qomp_vdw
#endif
end subroutine nonbond_qqp

!-----------------------------------------------------------------------

subroutine nonbond_qp(EQ_loc,lambda)
! arguments
TYPE(NB_ENERGIES)                               :: EQ_loc(:)
real(kind=prec)                                 :: lambda(:)
! local variables
integer                                         :: istate
integer                                         :: ip,iq,i,j
real(kind=prec)                                 :: Vel,V_a,V_b,dv
TYPE(ENERET_TYPE)                               :: nb_ene
TYPE(qr_dist3)                                  :: dist
#ifdef _OPENMP
integer :: quotient, remainder
#endif
! global variables used:
!  iqseq, iac, crg, x, nstates, qvdw_flag, iaclib, qiac, qavdw, qbvdw, qcrg, el14_scale, EQ, d, nat_solute
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (nbqp_pair)/threads_num
remainder = MOD(nbqp_pair, threads_num)
mp_start = thread_id * quotient + 1 + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

qomp_elec  = zero
qomp_vdw   = zero
do ip = mp_start,mp_end
#else
do ip = 1, nbqp_pair
#endif

iq = nbqp(ip,1)%i
i  = iqseq(iq)
j  = nbqp(ip,1)%j

dist = q_dist3(x(i),x(j))

do istate = 1, nstates
! for every assigned q-s pair:
nb_ene = nbe_qx(nbqp(ip,istate),lambda(istate),dist)

! calculate qi, Vel, V_a, V_b and dv
Vel  = nb_ene%Vel 
V_a  = nb_ene%V_a 
V_b  = nb_ene%V_b 
dv   = nb_ene%dv 

! update forces
d(i) = d(i) - dist%vec*dv
d(j) = d(j) + dist%vec*dv

  ! update q-protein or q-water energies
#ifdef _OPENMP
qomp_elec(istate) = qomp_elec(istate) + Vel
qomp_vdw(istate)  = qomp_vdw(istate)  + V_a - V_b
#else
EQ_loc(istate)%el  = EQ_loc(istate)%el  + Vel
EQ_loc(istate)%vdw = EQ_loc(istate)%vdw + V_a - V_b
#endif
end do ! istate

end do

#ifdef _OPENMP
do istate=1, nstates
!$omp atomic update
EQ_loc(istate)%el  = EQ_loc(istate)%el  + qomp_vdw(istate)
!$omp atomic update
EQ_loc(istate)%vdw = EQ_loc(istate)%vdw + qomp_elec(istate)
end do
#endif
end subroutine nonbond_qp
!-----------------------------------------------------------------------

!******PWadded 2001-10-23
subroutine nonbond_qp_box(EQ_loc,lambda)
! arguments
TYPE(NB_ENERGIES)                               :: EQ_loc(:)
real(kind=prec)                                 :: lambda(:)
! local variables
integer                                         :: ip,iq,i,j
integer                                         :: istate
real(kind=prec)                                 :: Vel,V_a,V_b,dv
integer                                         :: group, gr, ia
TYPE(qr_vec)                                    :: shift
TYPE(ENERET_TYPE)                               :: nb_ene
TYPE(qr_dist3)                                  :: dist
#ifdef _OPENMP
integer :: quotient, remainder
#endif
  ! global variables used:
  !  iqseq, iac, crg, x, nstates, qvdw_flag, iaclib, qiac, qavdw, qbvdw, qcrg, el14_scale, EQ, d, nat_solute

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (nbqp_cgp_pair)/threads_num
remainder = MOD(nbqp_cgp_pair, threads_num)
mp_start = thread_id * quotient + 1 + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
mp_end = mp_end + 1
endif

do gr = mp_start,mp_end
#else
!compute the peridocal shift for every charge group pair
do gr = 1, nbqp_cgp_pair
#endif
ia = nbqp_cgp(gr)%i !atom index for the atom

!the distance between the two switching atoms
shift = x(ia) - x(qswitch)
nbqp_cgp(gr)%shift = boxlength*q_nint(shift*inv_boxl)

end do
!$omp barrier
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (nbqp_pair)/threads_num
remainder = MOD(nbqp_pair, threads_num)
mp_start = thread_id * quotient + 1 + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
mp_end = mp_end + 1
endif
qomp_elec = zero
qomp_vdw  = zero
do ip = mp_start,mp_end
#else
do ip = 1, nbqp_pair
#endif


! init state-invariant variables:
group = nbqp(ip,1)%cgp_pair
iq    = nbqp(ip,1)%i
i     = iqseq(iq)
j     = nbqp(ip,1)%j

shift = x(i) - x(j)
dist  = q_dist3(shift,nbqp_cgp(group)%shift)


do istate = 1, nstates

nb_ene = nbe_qx(nbqp(ip,istate),lambda(istate),dist)

! calculate qi, Vel, V_a, V_b and dv
Vel  = nb_ene%Vel 
V_a  = nb_ene%V_a 
V_b  = nb_ene%V_b 
dv   = nb_ene%dv 

! update forces
d(i) = d(i) - dist%vec*dv
d(j) = d(j) + dist%vec*dv

! update q-protein or q-water energies
#ifdef _OPENMP
qomp_elec(istate) = qomp_elec(istate) + Vel
qomp_vdw(istate)  = qomp_vdw(istate)  + V_a - V_b
#else
EQ_loc(istate)%el  = EQ_loc(istate)%el  + Vel
EQ_loc(istate)%vdw = EQ_loc(istate)%vdw + V_a - V_b
#endif
end do ! istate

end do
#ifdef _OPENMP
!$omp atomic update
EQ_loc(istate)%el  = EQ_loc(istate)%el  + qomp_elec(istate)
!$omp atomic update
EQ_loc(istate)%vdw = EQ_loc(istate)%vdw + qomp_vdw(istate)
#endif
end subroutine nonbond_qp_box

!-----------------------------------------------------------------------

subroutine nonbond_qw(EQ_loc,lambda)
! arguments
TYPE(NB_ENERGIES)                               :: EQ_loc(:)
real(kind=prec)                                 :: lambda(:)
! local variables
integer                                         :: jw,iq,iw,i,j,jp
integer                                         :: istate
real(kind=prec)                                 :: dv,Vel,V_a,V_b
TYPE(qr_dist3)                                  :: dist
TYPE(ENERET_TYPE)                               :: nb_ene
#ifdef _OPENMP
integer :: quotient, remainder
#endif
! global variables used:
!  iqseq, iac, crg, x, nstates, qvdw_flag, iaclib, qiac, qavdw, qbvdw, qcrg, el14_scale, EQ, d, nat_solute


#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (nbww_pair/solv_atom)/threads_num
remainder = MOD(nbww_pair/solv_atom, threads_num)
mp_start = thread_id * quotient + 1 + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
mp_start = (mp_start*solv_atom) - solv_atom + 1

if (remainder .gt. thread_id) then
mp_end = mp_end + 1
endif
mp_end = (mp_end*solv_atom) - solv_atom + 1
qomp_elec = zero
qomp_vdw  = zero

do iw = mp_start,mp_end,solv_atom
#else
do iw = 1, nbqw_pair, solv_atom
#endif
! init state-invariant variables:
iq   = nbqw(iw,1)%i
i    = iqseq(iq)


do jw = 0, solv_atom-1
! for every assigned q-w pair:
! init state-invariant variables:
jp   = iw+jw
j    = nbqw(jp,1)%j

dist = q_dist3(x(i),x(j))

do istate = 1, nstates

nb_ene = nbe_qx(nbqw(jp,istate),lambda(istate),dist)
Vel = nb_ene%Vel
V_a = nb_ene%V_a
V_b = nb_ene%V_b
dv  = nb_ene%dv
! update forces
d(i) = d(i) - dist%vec*dv
d(j) = d(j) + dist%vec*dv

#ifdef _OPENMP
qomp_elec(istate) = qomp_elec(istate) + Vel
qomp_vdw(istate)  = qomp_vdw(istate)  + V_a - V_b
#else
EQ_loc(istate)%el  = EQ_loc(istate)%el  + Vel
EQ_loc(istate)%vdw = EQ_loc(istate)%vdw + V_a - V_b 
#endif

end do ! nstates
end do !jw

end do ! nbqw
#ifdef _OPENMP
!$omp atomic update
EQ_loc(istate)%el  = EQ_loc(istate)%el  + qomp_elec(istate)
!$omp atomic update
EQ_loc(istate)%vdw = EQ_loc(istate)%vdw + qomp_vdw(istate)
#endif
end subroutine nonbond_qw

!-----------------------------------------------------------------------
!******PWadded 2001-10-23
subroutine nonbond_qw_box(EQ_loc,lambda)
! arguments
TYPE(NB_ENERGIES)                               :: EQ_loc(:)
real(kind=prec)                                 :: lambda(:)
! local variables
integer                                         :: jw,iw,jp,i,j,iq
integer						:: istate
real(kind=prec)                                 :: dv,Vel,V_a,V_b
TYPE(qr_dist3)					:: dist
TYPE(qr_vec)					:: shift
TYPE(ENERET_TYPE)                               :: nb_ene
#ifdef _OPENMP
integer :: quotient, remainder
#endif
	! global variables used:
	!  iqseq, iac, crg, x, nstates, qvdw_flag, iaclib, qiac, qavdw, qbvdw, qcrg, el14_scale, EQ, d, nat_solute
  
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (nbww_pair/solv_atom)/threads_num
remainder = MOD(nbww_pair/solv_atom, threads_num)
mp_start = thread_id * quotient + 1 + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
mp_start = (mp_start*solv_atom) - solv_atom + 1

if (remainder .gt. thread_id) then
mp_end = mp_end + 1
endif
mp_end = (mp_end*solv_atom) - solv_atom + 1

qomp_elec = zero
qomp_vdw  = zero

do iw = mp_start,mp_end,solv_atom
#else
do iw = 1, nbqw_pair, solv_atom
#endif
! for every assigned q-s pair:
! init state-invariant variables:
iq    = nbqw(iw,1)%i
i     = iqseq(iq)
do jw = 0, solv_atom-1
! init state-invariant variables:
jp    = iw+jw
j     = nbqw(jp,1)%j

!compute the periodical shift
shift = x(qswitch) - x(j)
shift = boxlength*q_nint(shift*inv_boxl)
dist  = q_dist3(x(i)-x(j),shift)

do istate = 1, nstates

nb_ene = nbe_qx(nbqw(jp,istate),EQ(istate)%lambda,dist)
! calculate qi, Vel, V_a, V_b and dv
Vel  = nb_ene%Vel 
V_a  = nb_ene%V_a 
V_b  = nb_ene%V_b 
dv   = nb_ene%dv 

! update forces
d(i) = d(i) - dist%vec*dv
d(j) = d(j) + dist%vec*dv

#ifdef _OPENMP
qomp_elec(istate) = qomp_elec(istate) + Vel
qomp_vdw(istate)  = qomp_vdw(istate)  + V_a - V_b
#else
EQ_loc(istate)%el  = EQ_loc(istate)%el  + Vel
EQ_loc(istate)%vdw = EQ_loc(istate)%vdw + V_a - V_b
#endif

end do ! nstates
end do ! jw

end do ! nbqw
#ifdef _OPENMP
!$omp atomic update
EQ_loc(istate)%el  = EQ_loc(istate)%el  + qomp_elec(istate)
!$omp atomic update
EQ_loc(istate)%vdw = EQ_loc(istate)%vdw + qomp_vdw(istate)
#endif
end subroutine nonbond_qw_box

!----------------------------------------------------------------------------------------------------

subroutine nonbond_ww(E_loc)
! arguments
TYPE(NB_ENERGIES)                               :: E_loc
! local variables
integer						:: iw
integer						:: i,j
real(kind=prec)					:: Vel,V_a,V_b,dv
TYPE(ENERET_TYPE)                               :: nb_ene
#ifdef _OPENMP
integer :: quotient, remainder
real(kind=prec) :: Vel_omp,Vvdw_omp
#endif
! global variables used:
!  nat_solute, iac, crg, ljcod, iaclib, x, d, E

! totally rewritten now
! we are now using the nbww_pair index and nbww true list
! this function before was garbage
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()

quotient = (nbww_pair/(solv_atom**2))/threads_num
remainder = MOD(nbww_pair/(solv_atom**2), threads_num)

mp_start = thread_id * quotient + 1 + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
mp_start = (mp_start*(solv_atom**2)) - (solv_atom**2) + 1
if (remainder .gt. thread_id) then
mp_end = mp_end + 1
endif
mp_end = (mp_end*(solv_atom**2))

Vel_omp=zero
Vvdw_omp=zero
do iw = mp_start,mp_end,solv_atom
#else
do iw = 1, nbww_pair, solv_atom**2
#endif
i     = nbww(iw)%i
j     = nbww(iw)%j
nb_ene= nbe(nbww(iw))

Vel   = nb_ene%Vel 
V_a   = nb_ene%V_a 
V_b   = nb_ene%V_b 
dv    = nb_ene%dv 

d(i) = d(i) - nb_ene%vec*dv
d(j) = d(j) + nb_ene%vec*dv
#ifdef _OPENMP
Vel_omp  = Vel_omp  + Vel
Vvdw_omp = Vvdw_omp + V_a - V_b
#else
E_loc%el  = E_loc%el  + Vel       
E_loc%vdw = E_loc%vdw + V_a - V_b
#endif

end do ! iw
#ifdef _OPENMP
!$omp atomic update
E_loc%el  = E_loc%el + Vel_omp
!$omp atomic update
E_loc%vdw = E_loc%vdw + Vvdw_omp
#endif

end subroutine nonbond_ww

!----------------------------------------------------------------------------------------------------
subroutine nonbond_ww_box(E_loc)
! arguments
TYPE(NB_ENERGIES)                               :: E_loc
! local variables
integer						:: iw,ip,i,j,ia
real(kind=prec)					:: Vel,V_a,V_b,dv
TYPE(qr_vec)					:: shift
TYPE(ENERET_TYPE)                               :: nb_ene
#ifdef _OPENMP
integer :: quotient, remainder
real(kind=prec) :: Vel_omp,Vvdw_omp
#endif
! global variables used:
!  nat_solute, iac, crg, ljcod, iaclib, x, d, E

! totally rewritten now
! we are now using the nbww_pair index and nbww true list
! this function before was garbage
#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (nbww_pair/(solv_atom**2))/threads_num
remainder = MOD(nbww_pair/(solv_atom**2), threads_num)
mp_start = thread_id * quotient + 1 + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
mp_start = (mp_start*(solv_atom**2)) - (solv_atom**2) + 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif
mp_end = (mp_end*(solv_atom**2))

Vel_omp=zero
Vvdw_omp=zero
do iw = mp_start,mp_end,solv_atom**2
#else
do iw = 1, nbww_pair, solv_atom**2
#endif
i    = nbww(iw)%i
j    = nbww(iw)%j

!distance between this oxygen atom and the oxygen atom of the above watermolecule, iw
!this is always the first interaction :P
! periodic shift
shift = x(i) - x(j)
shift = boxlength*q_nint(shift*inv_boxl)

do ia = 1, (solv_atom**2) - 1
ip    = iw+ia

i     = nbww(ip)%i
j     = nbww(ip)%j

nb_ene= nbe_b(nbww(ip),shift)

Vel   = nb_ene%Vel 
V_a   = nb_ene%V_a 
V_b   = nb_ene%V_b 
dv    = nb_ene%dv 

d(i) = d(i) - nb_ene%vec*dv
d(j) = d(j) + nb_ene%vec*dv

#ifdef _OPENMP
Vel_omp  = Vel_omp  + Vel
Vvdw_omp = Vvdw_omp + V_a - V_b
#else
E_loc%el  = E_loc%el  + Vel       
E_loc%vdw = E_loc%vdw + V_a - V_b
#endif
end do ! ia
end do ! iw
#ifdef _OPENMP
!$omp atomic update
E_loc%el  = E_loc%el  + Vel_omp
!$omp atomic update
E_loc%vdw = E_loc%vdw + Vvdw_omp
#endif

end subroutine nonbond_ww_box

!-----------------------------------------------------------------------

subroutine nonbond_qw_spc(EQ_loc,lambda)
!calculate non-bonded interactions between Q-atoms and SPC water molecules
!(optimisations rely on LJ params = 0 for water H) using geometric comb. rule
! arguments
TYPE(NB_ENERGIES)                               :: EQ_loc(:)
real(kind=prec)                                 :: lambda(:)
! local variables
integer						:: iw,jw,iq,i,j,ip,ja,jp
integer						:: istate
real(kind=prec)					:: Vel, dv
real(kind=prec)					:: V_a, V_b
TYPE(qr_dist3)					:: dist
TYPE(qr_dist)					:: eldist
TYPE(ENERET_TYPE)                               :: nb_ene
#ifdef _OPENMP
integer :: quotient, remainder
#endif
! has been extensivly rewritten to use n atom spc solvent
! and the new nbqw list
! old function made Paul very sad

! global variables used:
!  iqseq, iac, crg, x, nstates, qvdw_flag, iaclib, qiac, qavdw, qbvdw, qcrg, el14_scale, EQ, d, nat_solute

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (nbqw_pair/solv_atom)/threads_num
remainder = MOD(nbqw_pair/solv_atom, threads_num)
mp_start = thread_id * quotient + 1 + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
mp_start = (mp_start*solv_atom) - solv_atom + 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif
mp_end = (mp_end*solv_atom) 

qomp_elec = zero
qomp_vdw  = zero

do iw = mp_start,mp_end,solv_atom
#else
do iw = 1, nbqw_pair, solv_atom
#endif
! init state-invariant variables:
iq   = nbqw(iw,1)%i
i    = iqseq(iq)
j    = nbqw(iw,1)%j

! only q - O distance first, only one with vdW
dist = q_dist3(x(i),x(j))

dv = zero
do istate = 1, nstates
nb_ene = nbe_qx(nbqw(iw,istate),lambda(istate),dist)
! calculate qi, Vel, V_a, V_b and dv
V_a    = nb_ene%V_a 
V_b    = nb_ene%V_b 
Vel    = nb_ene%Vel 
dv     = dv + nb_ene%dv 
! update q-water energies
#ifdef _OPENMP
qomp_elec(istate) = qomp_elec(istate) + Vel
qomp_vdw(istate)  = qomp_vdw(istate)  + V_a - V_b
#else
EQ_loc(istate)%el  = EQ_loc(istate)%el  + Vel
EQ_loc(istate)%vdw = EQ_loc(istate)%vdw + V_a - V_b
#endif
end do !istate
! update forces on q atom
d(i) = d(i) - dist%vec*dv
! update forces on water
d(j) = d(j) + dist%vec*dv

!now calculate only charge-charge interaction for hydrogens
do ip = 1, solv_atom-1

jp   = iw + ip
ja   = j +  ip

eldist = q_dist(x(i),x(ja))

dv = zero
do istate = 1, nstates
nb_ene = nbe_qspc(nbqw(jp,istate),lambda(istate),eldist)

Vel = nb_ene%Vel
#ifdef _OPENMP
qomp_elec(istate) = qomp_elec(istate) + Vel
#else
EQ_loc(istate)%el = EQ_loc(istate)%el + Vel
#endif
dv = dv + nb_ene%dv
end do ! istate
! update forces on q atom
d(i)  = d(i)  - eldist%vec*dv
! update forces on water
d(ja) = d(ja) + eldist%vec*dv
end do ! ip

end do ! nbqw
#ifdef _OPENMP
do istate = 1, nstates
!$omp atomic update
EQ_loc(istate)%vdw = EQ_loc(istate)%vdw + qomp_vdw(istate)
!$omp atomic update
EQ_loc(istate)%el  = EQ_loc(istate)%el  + qomp_elec(istate)
end do ! nstates
#endif
end subroutine nonbond_qw_spc

!-----------------------------------------------------------------------
!******PWadded 2001-10-23
subroutine nonbond_qw_spc_box(EQ_loc,lambda)
!calculate non-bonded interactions between Q-atoms and SPC water molecules
!(optimisations rely on LJ params = 0 for water H) using geometric comb. rule
! aeguments
TYPE(NB_ENERGIES)                               :: EQ_loc(:)
real(kind=prec)                                 :: lambda(:)
! local variables
integer                                         :: iw,iq,i,j,ip,ja,jp
integer						:: istate
real(kind=prec)                                 :: Vel, dv
real(kind=prec)                                 :: V_a, V_b
TYPE(qr_dist3)					:: dist
TYPE(qr_dist)					:: eldist
TYPE(qr_vec)					:: shift
TYPE(ENERET_TYPE)                               :: nb_ene
#ifdef _OPENMP
integer :: quotient, remainder
#endif

! has been extensivly rewritten to use n atom spc solvent
! and the new nbqw list

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (nbqw_pair/solv_atom)/threads_num
remainder = MOD(nbqw_pair/solv_atom, threads_num)
mp_start = thread_id * quotient + 1 + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
mp_start = (mp_start*solv_atom) - solv_atom + 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif
mp_end = (mp_end*solv_atom) 
qomp_elec = zero
qomp_vdw  = zero

do iw = mp_start,mp_end,solv_atom
#else

do iw = 1, nbqw_pair, solv_atom
#endif
! init state-invariant variables:
iq   = nbqw(iw,1)%i
i    = iqseq(iq)
j    = nbqw(iw,1)%j

shift = x(qswitch)-x(j)
shift = boxlength*q_nint(shift*inv_boxl)
dist  = q_dist3(x(i)-x(j),shift)


dv = zero
do istate = 1, nstates
nb_ene = nbe_qx(nbqw(iw,istate),lambda(istate),dist)

! calculate qi, Vel, V_a, V_b and dv
V_a = nb_ene%V_a 
V_b = nb_ene%V_b 
Vel = nb_ene%Vel 
dv  = dv + nb_ene%dv 
!softcore r6*r6_hc is (r^6/(r^6+alpha))
! update q-water energies
#ifdef _OPENMP
qomp_elec(istate) = qomp_elec(istate) + Vel
qomp_vdw(istate)  = qomp_vdw(istate)  + V_a - V_b
#else
EQ_loc(istate)%el  = EQ_loc(istate)%el  + Vel
EQ_loc(istate)%vdw = EQ_loc(istate)%vdw + V_a - V_b
#endif
end do !istate			
! update forces on q atom
d(i) = d(i) - dist%vec*dv
! update forces on water
d(j) = d(j) + dist%vec*dv

!now calculate only charge-charge interaction for hydrogens
do ip = 1, solv_atom-1

ja   = j  + ip
jp   = iw + ip
eldist = q_dist(x(i)-x(ja),shift)

dv = zero
do istate = 1, nstates

nb_ene = nbe_qspc(nbqw(jp,istate),lambda(istate),eldist)
Vel    = nb_ene%Vel
#ifdef _OPENMP
qomp_elec(istate) = qomp_elec(istate) + Vel
#else
EQ_loc(istate)%el = EQ_loc(istate)%el + Vel
#endif
dv = dv + nb_ene%dv
end do ! istate
! update forces on q atom
d(i)  = d(i)  - eldist%vec*dv
! update forces on water
d(ja) = d(ja) + eldist%vec*dv
end do ! ip

end do ! nbqw
#ifdef _OPENMP
do istate = 1, nstates
!$omp atomic update
EQ_loc(istate)%vdw = EQ_loc(istate)%vdw + qomp_vdw(istate)
!$omp atomic update
EQ_loc(istate)%el  = EQ_loc(istate)%el  + qomp_elec(istate)
end do ! nstates
#endif
end subroutine nonbond_qw_spc_box


!-----------------------------------------------------------------------

subroutine nonbond_ww_spc(E_loc)
! function has been totally rewritten to be used for n - atom spc solvent
! only vdW interactions between first two atoms are calculated -> those are the heavy atoms
! rest of the atoms only interact via coloumb interactions
! old function gave me nightmares
! arguments
TYPE(NB_ENERGIES)                               :: E_loc
! local variables
integer                                         :: iw,jw,i,j,ip
real(kind=prec)                                 :: Vel,V_a,V_b,dv
TYPE(ENERET_TYPE)                               :: nb_ene
#ifdef _OPENMP
integer :: quotient, remainder
real(kind=prec) :: Vel_omp, Vvdw_omp
#endif

! the loop is now organised in the following way (as before for the arithmetic case)
! we always jump solv_atom positions, to calculate atomAn-atomB(1-solv_atom)
! and iterate a separate counter to get the current position of the atomA(1-solv_atom) position
! vdW is only calculated if atomA1 - atomB1

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (nbww_pair/(solv_atom**2))/threads_num
remainder = MOD(nbww_pair/(solv_atom**2), threads_num)
mp_start = thread_id * quotient + 1 + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
mp_start = (mp_start*(solv_atom**2)) - (solv_atom**2) + 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif
mp_end = (mp_end*(solv_atom**2)) 

Vel_omp = zero
Vvdw_omp = zero
do iw = mp_start,mp_end,solv_atom**2
#else
do iw = 1, nbww_pair, solv_atom**2
#endif

i    = nbww(iw)%i
j    = nbww(iw)%j

nb_ene = nbe(nbww(iw))
V_a = nb_ene%V_a 
V_b = nb_ene%V_b 
Vel = nb_ene%Vel 

#ifdef _OPENMP
Vvdw_omp = Vvdw_omp + V_a - V_b
Vel_omp = Vel_omp + Vel
#else
E_loc%vdw = E_loc%vdw + V_a - V_b
E_loc%el  = E_loc%el  + Vel
#endif
dv = nb_ene%dv
d(i) = d(i) - nb_ene%vec*dv
d(j) = d(j) + nb_ene%vec*dv


! now calc the rest of the interactions
! remember that we need to switch the atoms every solv_atom steps

do jw = 1, (solv_atom**2) - 1

ip    = iw+jw
i     = nbww(ip)%i
j     = nbww(ip)%j

nb_ene = nbe_spc(nbww(ip))

Vel  = nb_ene%Vel
#ifdef _OPENMP
Vel_omp  = Vel_omp  + Vel
#else
E_loc%el = E_loc%el + Vel
#endif

dv = nb_ene%dv
d(i) = d(i) - nb_ene%vec*dv
d(j) = d(j) + nb_ene%vec*dv

end do
end do ! iw
#ifdef _OPENMP
!$omp atomic update
E_loc%el  = E_loc%el  + Vel_omp
!$omp atomic update
E_loc%vdw = E_loc%vdw + Vvdw_omp
#endif
end subroutine nonbond_ww_spc

!-----------------------------------------------------------------------
!******PWadded 2001-10-23
subroutine nonbond_ww_spc_box(E_loc)

! function has been totally rewritten to be used for n - atom spc solvent
! only vdW interactions between first two atoms are calculated -> those are the heavy atoms
! rest of the atoms only interact via coloumb interactions
! old function gave me nightmares
! arguments
TYPE(NB_ENERGIES)                               :: E_loc
  ! local variables
integer                                         :: iw,jw,ip,i,j
real(kind=prec)					:: Vel,V_a,V_b,dv
TYPE(qr_vec)					:: shift
TYPE(ENERET_TYPE)                               :: nb_ene
#ifdef _OPENMP
integer :: quotient, remainder
real(kind=prec) :: Vel_omp, Vvdw_omp
#endif

! the loop is now organised in the following way (as before for the arithmetic case)
! we always jump solv_atom positions, to calculate atomAn-atomB(1-solv_atom)
! and iterate a separate counter to get the current position of the atomA(1-solv_atom) position
! vdW is only calculated if atomA1 - atomB1

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (nbww_pair/(solv_atom**2))/threads_num
remainder = MOD(nbww_pair/(solv_atom**2), threads_num)
mp_start = thread_id * quotient + 1 + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
mp_start = (mp_start*(solv_atom**2)) - (solv_atom**2) + 1
if (remainder .gt. thread_id) then
mp_end = mp_end + 1
endif
mp_end = (mp_end*(solv_atom**2)) 

Vel_omp = zero
Vvdw_omp = zero

do iw = mp_start,mp_end,solv_atom**2
#else
do iw = 1, nbww_pair, solv_atom**2
#endif
i    = nbww(iw)%i
j    = nbww(iw)%j


! very first one is the only calc with vdW
! and the one defining the box shift
shift = x(i)-x(j)
shift = boxlength*q_nint(shift*inv_boxl)

nb_ene = nbe_b(nbww(iw),shift)

V_a = nb_ene%V_a
V_b = nb_ene%V_b 
Vel = nb_ene%Vel 

#ifdef _OPENMP
Vvdw_omp = Vvdw_omp + V_a - V_b
Vel_omp  = Vel_omp  + Vel
#else
E_loc%vdw = E_loc%vdw + V_a - V_b
E_loc%el  = E_loc%el  + Vel
#endif
dv = nb_ene%dv
d(i) = d(i) - nb_ene%vec*dv
d(j) = d(j) + nb_ene%vec*dv


! now calc the rest of the interactions

do jw = 1, (solv_atom**2) - 1

ip    = iw + jw
i     = nbww(ip)%i
j     = nbww(ip)%j

nb_ene = nbe_spcb(nbww(ip),shift)

Vel  = nb_ene%Vel
#ifdef _OPENMP
Vel_omp  = Vel_omp  + Vel
#else
E_loc%el = E_loc%el + Vel
#endif
dv = nb_ene%dv
d(i) = d(i) - nb_ene%vec*dv
d(j) = d(j) + nb_ene%vec*dv

end do
end do ! iw
#ifdef _OPENMP
!$omp atomic update
E_loc%el  = E_loc%el  + Vel_omp
!$omp atomic update
E_loc%vdw = E_loc%vdw + Vvdw_omp
#endif
end subroutine nonbond_ww_spc_box


!-----------------------------------------------------------------------

subroutine nonbond_solvent_internal(E_loc)
! new function that handles non bonded interactions within a solvent molecule
! only called for solvents with several heavy atoms and if the solvent has torsions
! all interactions are precomputed, only need distances
! each node only works on its own solvent atoms
! arguments
TYPE(NB_ENERGIES)                               :: E_loc
! local variables
integer                                         :: iw,is,ip,ia,ib
real(kind=prec)					:: r2,r6,r,r12
real(kind=prec)					:: Vel,V_a,V_b,dv
TYPE(qr_dist2)					:: distance
#ifdef _OPENMP
integer :: quotient, remainder
real(kind=prec) :: Vel_omp, Vvdw_omp
#endif

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%ww%end - calculation_assignment%ww%start + 1)/threads_num
remainder = MOD(calculation_assignment%ww%end - calculation_assignment%ww%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%ww%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

Vel_omp = zero
Vvdw_omp = zero

do iw = mp_start,mp_end
#else
do iw = calculation_assignment%ww%start, calculation_assignment%ww%end
#endif
is  = nat_solute + solv_atom*iw-(solv_atom)

! loop over all interactions within this solvent molecule
do ip = 1 , num_solv_int 
ia   = is + nonbnd_solv_int(ip)%i
ib   = is + nonbnd_solv_int(ip)%j

distance = q_dist2(x(ia),x(ib))

r2  = distance%r2 
r   = distance%r 
r6  = distance%r6 
r12 = distance%r12 

Vel  = nonbnd_solv_int(ip)%elec * r
V_a  = nonbnd_solv_int(ip)%vdWA * r12
V_b  = nonbnd_solv_int(ip)%vdWB * r6

dv   = r2*( -Vel) + r2*(-12.0_prec*V_a +6.0_prec*V_b )
#ifdef _OPENMP
Vel_omp  = Vel_omp  + Vel
Vvdw_omp = Vvdw_omp + V_a - V_b
#else
E_loc%el  = E_loc%el  + Vel
E_loc%vdw = E_loc%vdw + V_a - V_b
#endif
d(ia) = d(ia) - distance%vec*dv
d(ib) = d(ib) + distance%vec*dv
end do
end do
#ifdef _OPENMP
!$omp atomic update
E_loc%el  = E_loc%el  + Vel_omp
!$omp atomic update
E_loc%vdw = E_loc%vdw + Vvdw_omp
#endif

end subroutine nonbond_solvent_internal

!-----------------------------------------------------------------------

subroutine nonbond_solvent_internal_box(E_loc)
! new function that handles non bonded interactions within a solvent molecule
! only called for solvents with several heavy atoms and if the solvent has torsions
! all interactions are precomputed, only need distances
! each node only works on its own solvent atoms
! arguments
TYPE(NB_ENERGIES)                               :: E_loc
! local variables
integer                                         :: iw,is,ip,ia,ib
real(kind=prec)                                 :: r2,r6,r,r12
real(kind=prec)                                 :: Vel,V_a,V_b,dv
TYPE(qr_dist2)					:: distance
TYPE(qr_vec)					:: shift
#ifdef _OPENMP
integer :: quotient, remainder
real(kind=prec) :: Vel_omp, Vvdw_omp
#endif

#ifdef _OPENMP
threads_num = omp_get_num_threads()
thread_id = omp_get_thread_num()
quotient = (calculation_assignment%ww%end - calculation_assignment%ww%start + 1)/threads_num
remainder = MOD(calculation_assignment%ww%end - calculation_assignment%ww%start + 1, threads_num)
mp_start = thread_id * quotient + calculation_assignment%ww%start + MIN(thread_id, remainder)
mp_end = mp_start + quotient - 1
if (remainder .gt. thread_id) then
    mp_end = mp_end + 1
endif

Vel_omp = zero
Vvdw_omp = zero

do iw = mp_start,mp_end
#else
do iw = calculation_assignment%ww%start, calculation_assignment%ww%end
#endif
is  = nat_solute + solv_atom*iw-(solv_atom)

! loop over all interactions within this solvent molecule
do ip = 1 , num_solv_int
ia   = is + nonbnd_solv_int(ip)%i
ib   = is + nonbnd_solv_int(ip)%j

shift = x(ia)-x(ib)
shift = boxlength*q_nint(shift*inv_boxl)
distance = q_dist2(x(ia)-x(ib),shift)

r2  = distance%r2
r   = distance%r
r6  = distance%r6
r12 = distance%r12

Vel  = nonbnd_solv_int(ip)%elec * r
V_a  = nonbnd_solv_int(ip)%vdWA * r12
V_b  = nonbnd_solv_int(ip)%vdWB * r6

dv   = r2*( -Vel) + r2*(-12.0_prec*V_a +6.0_prec*V_b )
#ifdef _OPENMP
Vel_omp  = Vel_omp  + Vel
Vvdw_omp = Vvdw_omp + V_a - V_b
#else
E_loc%el  = E_loc%el  + Vel
E_loc%vdw = E_loc%vdw + V_a - V_b
#endif
d(ia) = d(ia) - distance%vec*dv
d(ib) = d(ib) + distance%vec*dv
end do
end do
#ifdef _OPENMP
!$omp atomic update
E_loc%el  = E_loc%el  + Vel_omp
!$omp atomic update
E_loc%vdw = E_loc%vdw + Vvdw_omp
#endif


end subroutine nonbond_solvent_internal_box

!-----------------------------------------------------------------------

subroutine offdiag
! local variables
integer						:: io,i,j,k,l
real(kind=prec)					:: r

! global variables used:
!  offd, noffd, iqseq, x, Hij, offd2

do io = 1, noffd
! for every offd:

i  = offd(io)%i
j  = offd(io)%j
k  = iqseq(offd2(io)%k)
l  = iqseq(offd2(io)%l)

r = q_sqrt(q_dist4(x(k),x(l)))

Hij(i,j) = offd2(io)%A * exp(-offd2(io)%mu*r)
offd(io)%Hij = Hij(i,j)	! store for save
offd(io)%rkl = r
end do
end subroutine offdiag

!-----------------------------------------------------------------------
#ifdef USE_GRID
subroutine populate_grids
! this function searches the number of charge groups and number of water molecules
! and places them into the new md grids
! this is done every NB steps in make_pair_lists as the first thing
! in the end master broadcasts the new grids to the slave nodes
! needs some trickery so that we don't overflow any arrays

integer					:: ichg,igrid
integer					:: i
TYPE(qr_vec)                            :: boxmin,boxmax,xc
logical                                 :: set = .false.
#ifdef USE_MPI
integer, parameter                      :: vars = 40    !increment this var when adding data to broadcast in batch 1
integer                                 :: blockcnt(vars), ftype(vars)
integer(kind=MPI_ADDRESS_KIND)                          :: fdisp(vars)
#endif

!each node now does its share

if (use_PBC) then
boxmin = boxcentre-boxlength/2.0_prec
boxmax = boxcentre+boxlength/2.0_prec
end if

! reset the group information for all grids to zero
! and null all arrays
grid_pp_ngrp(:)=0
grid_pw_ngrp(:)=0
grid_ww_ngrp(:)=0
ww_igrid(:)=0
pw_igrid(:)=0
pp_igrid(:)=0
grid_pp_grp(:,:)=0
grid_pw_grp(:,:)=0
grid_ww_grp(:,:)=0

do ichg=calculation_assignment%pp%start,calculation_assignment%pp%end
! we decide the grid by the position of the switch atom of each charge group
i   = cgp(ichg)%iswitch
xc = x(i)
! this is to account for the possibility of atoms to move out of the PBC box if they are not put back in
! in this case we just shift the atom implicitly by the distance
! this is just a bookkeeping issue
do while
(use_PBC.and.((xc%x.gt.(boxmax%x)).or.(xc%x.lt.(boxmin%x)).or.(xc%y.gt.(boxmax%y)).or.(xc%y.lt.(boxmin%y)).or.(xc%z.gt.(boxmax%z)).or.(xc%z.lt.(boxmin%z))))
if (xc%x.gt.boxmax%x) then
xc%x = boxmin%x + (xc%x-boxmax%x)
elseif (xc%x.lt.boxmin%x) then
xc%x = boxmax%x - (xc%x-boxmin%x)
end if
if (xc%y.gt.boxmax%y) then
xc%y = boxmin%y + (xc%y-boxmax%y)
elseif (xc%y.lt.boxmin%y) then
xc%y = boxmax%y - (xc%y-boxmin%y)
end if
if (xc%z.gt.boxmax%z) then
xc%z = boxmin%z + (xc%z-boxmax%z)
elseif (xc%z.lt.boxmin%z) then
xc%z = boxmax%z - (xc%z-boxmin%z)
end if
end do

set   = .false.
igrid = 1
! now loop first over the pp_grid and then over the pw_grid for ncgp_solute
pploop: do while (( igrid .le. pp_gridnum).and.(.not.set))
		if ((xc%x.ge.grid_pp(igrid)%gstart%x).and.(xc%x.lt.grid_pp(igrid)%gend%x)) then
			if ((xc%y.ge.grid_pp(igrid)%gstart%y).and.(xc%y.lt.grid_pp(igrid)%gend%y)) then
				if ((xc%z.ge.grid_pp(igrid)%gstart%z).and.(xc%z.lt.grid_pp(igrid)%gend%z)) then  
! yeah, we found our grid, now place the atom in here and set all the bookkeeping stuff
! then cycle the grid and continue with the next one
! needs temp storage because of the OMP part!!!!!
				pp_igrid(ichg) = igrid
				grid_pp_ngrp(igrid) = grid_pp_ngrp(igrid) + 1
				grid_pp_grp(igrid,grid_pp_ngrp(igrid)) = ichg
                                set = .true.
				cycle pploop
! this should (!!!!) ensure we place all groups, but needs to be tested
				end if
			end if
		end if
                igrid = igrid + 1
	end do pploop
! now the pw_grid, same procedure
set   = .false.
igrid = 1
pwloop1: do while (( igrid .le. pw_gridnum).and.(.not.set))
		if ((xc%x.ge.grid_pw(igrid)%gstart%x).and.(xc%x.lt.grid_pw(igrid)%gend%x)) then
			if ((xc%y.ge.grid_pw(igrid)%gstart%y).and.(xc%y.lt.grid_pw(igrid)%gend%y)) then
				if ((xc%z.ge.grid_pw(igrid)%gstart5z).and.(xc%z.lt.grid_pw(igrid)%gend%z)) then
				pw_igrid(ichg) = igrid
                                set = .true.
				cycle pwloop1
				end if
			end if
		end if
                igrid = igrid + 1
	end do pwloop1
end do ! ncgp_solute

! now do the same for nwat
do ichg=calculation_assignment%ww%start,calculation_assignment%ww%end
! we decide the grid by the position of the first atom of each solvent molecule
i   = nat_solute+(ichg*solv_atom)-(solv_atom-1)
xc = x(i)
! this is to account for the possibility of atoms to move out of the PBC box if they are not put back in
! in this case we just shift the atom implicitly by the distance
! this is just a bookkeeping issue
do while
(use_PBC.and.((xc%x.gt.(boxmax%x)).or.(xc%x.lt.(boxmin%x)).or.(xc%y.gt.(boxmax%y)).or.(xc%y.lt.(boxmin%y)).or.(xc%z.gt.(boxmax%z)).or.(xc%z.lt.(boxmin%z))))
if (xc%x.gt.boxmax%x) then
xc%x = boxmin%x + (xc%x-boxmax%x)
elseif (xc%x.lt.boxmin%x) then
xc%x = boxmax%x - (xc%x-boxmin%x)
end if
if (xc%y.gt.boxmax%y) then
xc%y = boxmin%y + (xc%y-boxmax%y)
elseif (xc%y.lt.boxmin%y) then
xc%y = boxmax%y - (xc%y-boxmin%y)
end if
if (xc%z.gt.boxmax%z) then
xc%z = boxmin%z + (xc%z-boxmax%z)
elseif (xc%z.lt.boxmin%z) then
xc%z = boxmax%z - (xc%z-boxmin%z)
end if
end do

set   = .false.
igrid = 1
pwloop2: do while (( igrid .le. pw_gridnum).and.(.not.set)) 
                if ((xc%x.ge.grid_pw(igrid)%gstart%x).and.(xc%x.lt.grid_pw(igrid)%gend%x)) then
                        if ((xc%y.ge.grid_pw(igrid)%gstart%y).and.(xc%y.lt.grid_pw(igrid)%gend%y)) then
                                if ((xc%z.ge.grid_pw(igrid)%gstart%z).and.(xc%z.lt.grid_pw(igrid)%gend%z)) then
				grid_pw_ngrp(igrid) = grid_pw_ngrp(igrid) + 1
				grid_pw_grp(igrid,grid_pw_ngrp(igrid)) = ichg
                                set = .true.
                                cycle pwloop2
                                end if
                        end if
                end if
                igrid = igrid + 1
        end do pwloop2
set   = .false.
igrid = 1
wwloop: do while (( igrid .le. ww_gridnum).and.(.not.set))
		if ((xc%x.ge.grid_ww(igrid)%gstart%x).and.(xc%x.lt.grid_ww(igrid)%gend%x)) then
			if ((xc%y.ge.grid_ww(igrid)%gstart%y).and.(xc%y.lt.grid_ww(igrid)%gend%y)) then
				if ((xc%z.ge.grid_ww(igrid)%gstart%z).and.(xc%z.lt.grid_ww(igrid)%gend%z)) then
				grid_ww_ngrp(igrid) = grid_ww_ngrp(igrid) + 1
				grid_ww_grp(igrid,grid_ww_ngrp(igrid)) = ichg
				ww_igrid(ichg) = igrid
                                set = .true.
				cycle wwloop
				end if
			end if
		end if
                igrid = igrid + 1
	end do wwloop
end do ! nwat


! now collect and distribute all the info to all nodes
#ifdef USE_MPI

if (ncgp_solute .ne. 0 ) then

call MPI_Allreduce(MPI_IN_PLACE,grid_pp_ngrp,pp_gridnum,MPI_INTEGER,MPI_SUM, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('populate grids MPI_Allreduce grid_pp_ngrp')
call MPI_Allreduce(MPI_IN_PLACE,pp_igrid,ncgp_solute,MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('populate grids MPI_Allreduce pp_igrid')

do igrid = 1, pp_gridnum
call MPI_Allreduce(MPI_IN_PLACE,grid_pp_grp(igrid,:),gridstor_pp,MPI_INTEGER,&
                        mpi_grid_add, MPI_COMM_WORLD, ierr)
end do
if (ierr .ne. 0) call die('populate grids MPI_Allreduce grid_pp_grp')

end if
if (nwat .gt. 0 ) then
call MPI_Allreduce(MPI_IN_PLACE,grid_pw_ngrp,pw_gridnum,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('populate grids MPI_Allreduce grid_pw_ngrp')
call MPI_Allreduce(MPI_IN_PLACE,pw_igrid,ncgp_solute,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('populate grids MPI_Allreduce pw_igrid')
call MPI_Allreduce(MPI_IN_PLACE,grid_ww_ngrp,ww_gridnum,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('populate grids MPI_Allreduce grid_ww_ngrp')
call MPI_Allreduce(MPI_IN_PLACE,ww_igrid,nwat,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('populate grids MPI_Allreduce ww_igrid')

do igrid = 1, pw_gridnum
call MPI_Allreduce(MPI_IN_PLACE,grid_pw_grp(igrid,:),gridstor_pw,MPI_INTEGER,&
                        mpi_grid_add, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('populate grids MPI_Allreduce grid_pw_grp')
end do
do igrid = 1, ww_gridnum
call MPI_Allreduce(MPI_IN_PLACE,grid_ww_grp(igrid,:),gridstor_ww,MPI_INTEGER,&
                        mpi_grid_add, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('populate grids MPI_Allreduce grid_ww_grp')
end do
end if
#endif
 

end subroutine populate_grids
#endif

!-----------------------------------------------------------------------


subroutine restrain_solvent(E_loc)
! arguments
TYPE(RESTRAINT_ENERGIES)                        :: E_loc
! local variables
integer						:: iw,i,isolv,jsolv
real(kind=prec)					:: db,erst,dv,fexp,b
TYPE(qr_vec)					:: dcent
real(kind=prec)					:: shift
TYPE(qr_dist5)                                  :: distance

! global variables used:
!  E, Boltz, Tfree, fk_wsphere, nwat, nat_pro, x, xwcent, rwat, Dwmz, awmz, d

if(fk_wsphere .ne. zero) then
        shift = sqrt (Boltz*Tfree/fk_wsphere)
else
        shift = zero
end if

! make decision here to use old or new code, aviod branch points during loop
if(solv_atom.eq.3) then
do iw = ncgp_solute + 1, ncgp
        i  = cgp(iw)%iswitch
        if (excl(i)) cycle ! skip excluded topology waters
        distance = q_dist5(xwcent,x(i))
        b  = q_sqrt(distance%r2)
        db = b - (rwat - shift)
        ! calculate erst and dv
        if ( db .gt. zero ) then
                erst = 0.5_prec * fk_wsphere * db**2 - Dwmz
                dv = fk_wsphere*db/b
        else
                if (b .gt. zero) then
                  fexp = q_exp ( awmz*db )
                  erst = Dwmz*(fexp*fexp-2.0_prec*fexp)
                  dv = -2.0_prec*Dwmz*awmz*(fexp-fexp*fexp)/b
                else
                  dv = zero
                  erst = zero
                end if
        end if
        d(i) = d(i) + distance%vec*dv
        E_loc%solvent_radial = E_loc%solvent_radial + erst
        end do
else
do iw = ncgp_solute + 1, ncgp
        i  = cgp(iw)%iswitch
        if (excl(i)) cycle ! skip excluded topology waters
        dcent = zero
        jsolv = iw - ncgp_solute
        i = nat_solute + solv_atom*jsolv - (solv_atom-1)
        do isolv = 0 , solv_atom - 1
                dcent = dcent + x(i+isolv)
        end do
        dcent = dcent/real(solv_atom,kind=prec)
        distance = q_dist5(xwcent,dcent)
        b = q_sqrt(distance%r2)
        db = b - (rwat - shift)

        ! calculate erst and dv
        if ( db .gt. zero ) then
                erst = 0.5_prec * fk_wsphere * db**2 - Dwmz
                dv = fk_wsphere*db/b
        else
                if (b .gt. zero) then
                  fexp = q_exp ( awmz*db )
                  erst = Dwmz*(fexp*fexp-2.0_prec*fexp)
                  dv = -2.0_prec*Dwmz*awmz*(fexp-fexp*fexp)/b
                else
                  dv = zero
                  erst = zero
                end if
        end if
        dv = dv / solv_atom

        ! update energy and forces
        E_loc%solvent_radial = E_loc%solvent_radial + erst
        jsolv = iw - ncgp_solute
        i = nat_solute + solv_atom*jsolv - (solv_atom-1)
        do isolv = 0 , solv_atom - 1
               d(i+isolv) = d(i+isolv) + distance%vec*dv 
        end do
end do
end if
end subroutine restrain_solvent

!-----------------------------------------------------------------------

subroutine watpol(E_loc,md)
! arguments
TYPE(RESTRAINT_ENERGIES)                :: E_loc
logical                                 :: md
! local variables
integer					:: iw,is,i,il,jl,jw,imin,jmin,isolv,jsolv,j3
real(kind=prec)				:: dr,rw,rshell,rm,rc,scp
real(kind=prec)				:: tmin,arg,avtdum,dv,f0
TYPE(qr_vec)                            :: f1,f2,f3,rmu,rcu,rmc

! global variables used:
!  E, wshell, bndw0, deg2rad, angw0, nwat, theta, theta0, nat_pro, x, xwcent,
!  tdum, nwpolr_shell, list_sh, pi, nsort, istep, itdis_update, fkwpol, d

! reset wshell%n_insh
if (md) then
 wshell(:)%n_insh = 0

! calculate theta(:), tdum(:), wshell%n_insh
do iw = 1, nwat

theta(iw)  = zero
theta0(iw) = zero

i  = nat_solute + solv_atom*iw - (solv_atom-1) 
if(excl(i)) cycle ! skip excluded topology solvent

! function needs to be rewritten as we no longer just have 3-point water
! molecules, meaning that the simple geometry won't work any longer
! we will only keep the old code for SPC/TIP3P like solvent

rmu = zero
rmc = zero

if (solv_atom .eq. 3) then
!SPC/TIP3P solvent, keep old code
do isolv = 1, solv_atom - 1
jsolv = i + isolv
rmu = rmu + x(jsolv)
end do
rmu = rmu - x(i)*2.0_prec
else
! no longer 3-point solvent, create solvent vector using individual atom charges
! of the solvent molecule
do isolv = 0, solv_atom - 1
jsolv = i + isolv
rmc = rmc + x(jsolv)
rmu = rmu + x(jsolv) * chg_solv(isolv+1)
end do
rmc = rmc/real(solv_atom,kind=prec)
end if

rm = q_sqrt(qvec_square(rmu))
rmu = rmu/rm

if (solv_atom .eq. 3) then
!SPC/TIP3P solvent, keep old code
rcu = x(i) - xwcent !Radial vector to center solvent atom
else
! different solvent, use new code for geometric center
rcu = rmc - xwcent
end if
rc = q_sqrt(qvec_square(rcu))
rcu = rcu/rc
scp = q_dotprod(rmu,rcu) !Calculate angle between solvent vector and radial vector
if ( scp .gt.  one ) scp =  one
if ( scp .lt. -one ) scp = -one
theta(iw) = q_acos( scp )
tdum(iw) = theta(iw)

if ( rc > wshell(nwpolr_shell)%rout-wshell(nwpolr_shell)%dr ) then
  do is = nwpolr_shell, 2, -1
        if(rc <= wshell(is)%rout) exit
  end do
  wshell(is)%n_insh = wshell(is)%n_insh + 1
  if (wshell(is)%n_insh .ge. n_max_insh) then
	call reallocate_watpol_arrays
  end if
  list_sh(wshell(is)%n_insh,is) = iw
end if
end do

! sort the solvent molecules according to theta
do is = 1, nwpolr_shell
imin = 0
do il = 1, wshell(is)%n_insh
tmin = 2.0_prec*pi
do jl = 1, wshell(is)%n_insh
jw = list_sh(jl,is)
if ( tdum(jw) .lt. tmin ) then
  jmin = jw
  tmin = theta(jw)
end if
end do
imin = imin+1
nsort(imin,is) = jmin
  tdum(jmin) = 99999.0_prec

end do

end do

end if ! if md

! calculate energy and force
if ( istep .ne. 0 .and. mod(istep,itdis_update) .eq. 0 .and. md) then
call centered_heading('Solvent polarisation restraint data', '-')
write(*,'(a)') 'shell    <n>    <theta>    theta_0 theta_corr'
do is = 1, nwpolr_shell
wshell(is)%avtheta = wshell(is)%avtheta / real (itdis_update, kind=prec)
wshell(is)%avn_insh = wshell(is)%avn_insh / real (itdis_update, kind=prec)
wshell(is)%theta_corr = wshell(is)%theta_corr + wshell(is)%avtheta-acos(wshell(is)%cstb)
write (*,10) is,wshell(is)%avn_insh,wshell(is)%avtheta/deg2rad, &
     acos(wshell(is)%cstb)/deg2rad,wshell(is)%theta_corr/deg2rad
10	  format(i5,1x,f6.1,3x,f8.3,3x,f8.3,3x,f8.3)
wshell(is)%avtheta = zero
wshell(is)%avn_insh = zero
end do
end if

do is = 1, nwpolr_shell
if(wshell(is)%n_insh == 0) cycle !skip empty shell
avtdum = zero
do il = 1, wshell(is)%n_insh
iw = nsort(il,is)
arg = one + (one - 2.0_prec*real(il, kind=prec))/real(wshell(is)%n_insh, kind=prec)
theta0(il) = acos ( arg )
theta0(il) = theta0(il)-3.0_prec*sin(theta0(il))*wshell(is)%cstb/2.0_prec
if ( theta0(il) .lt. zero ) theta0(il) = zero
if ( theta0(il) .gt. pi)   theta0(il) = pi

avtdum = avtdum + theta(iw)

E_loc%water_pol = E_loc%water_pol + 0.5_prec*fkwpol* &
     (theta(iw)-theta0(il)+wshell(is)%theta_corr)**2

dv = fkwpol*(theta(iw)-theta0(il)+wshell(is)%theta_corr)

i  = nat_solute + solv_atom*iw - (solv_atom-1)
rmu= zero ! solvent vector
rmc= zero
if (solv_atom .eq. 3) then
!SPC/TIP3P solvent, keep old code
do isolv = 1, solv_atom - 1
jsolv = i + isolv
rmu = rmu + x(jsolv)
end do
rmu = rmu - x(i)* 2.0_prec
else
! no longer 3-point solvent, create solvent vector based on individual charges of
! the solvent molecule
do isolv = 0, solv_atom - 1
jsolv = i + isolv
rmc = rmc + x(jsolv)
rmu = rmu + x(jsolv) * chg_solv(isolv+1)
end do
rmc = rmc/real(solv_atom,kind=prec)
end if

rm = q_sqrt(qvec_square(rmu))
rmu = rmu / rm

if (solv_atom .eq. 3) then
!SPC/TIP3P solvent, keep old code
rcu = x(i) - xwcent

else
! different solvent, use new code for geometric center
rcu = rmc - xwcent
end if
rc = q_sqrt(qvec_square(rcu))
rcu = rcu / rc

scp = q_dotprod(rmu,rcu)
if ( scp .gt.  one ) scp =  one
if ( scp .lt. -one ) scp = -one
f0 = q_sin ( q_acos(scp) )
if ( abs(f0) .lt. 1.e-12_prec ) f0 = 1.e-12_prec
f0 = -one / f0
f0 = dv*f0

f1 = ((rcu - (rmu*scp)) * (-2.0_prec))/rm
f3 = (rcu - (rmu*scp))/rm
f2 = (rmu - (rcu*scp))/rc

d(i) = d(i) + (f1+f2)*f0

do isolv = 1, solv_atom - 1
jsolv = i + isolv
d(jsolv) = d(jsolv) + (f3*f0)

end do
end do
if (md) then
wshell(is)%avtheta = wshell(is)%avtheta + avtdum/real(wshell(is)%n_insh, kind=prec)
wshell(is)%avn_insh = wshell(is)%avn_insh + wshell(is)%n_insh
end if
end do
end subroutine watpol

!----------------------------------------------------------------------------


#if defined (USE_MPI)
!***********************
!subroutine handling summation of nonbonded energies from slave nodes.
!***********************
! Use the global vars
!  request_recv, E_send,EQ_send,E_recv,EQ_Recv,d_recv
! Allocate  - status


subroutine gather_nonbond(E_loc,EQ_loc)
! arguments
TYPE(ENERGIES)                          :: E_loc
TYPE(OQ_ENERGIES)                       :: EQ_loc(:)
integer,parameter                       :: vars=3
integer,dimension(3,numnodes-1)         :: tag
integer,dimension(vars)	                :: blockcnt,ftype 
integer(kind=MPI_ADDRESS_KIND), dimension(vars)	:: fdisp, base
integer                                 :: mpitype_package,mpitype_send
integer                                 :: i,j,istate

do i=1,numnodes-1
tag(1,i)=numnodes*100+i
tag(2,i)=numnodes*200+i
tag(3,i)=numnodes*300+i
end do

if (nodeid .eq. 0) then        !master

!reclength is defined in mpiglob and set in prep_sim
!used for the now not previously known length of the EQ array
!for sending same length distributed to nodes

! Post receives for each of the d/E/EQ_recv structures
! E/EQ_Recv should really be handled with MPI_Type_create_struct
! and d_recv's type should be handled correctly (it's KIND=wp8)
! should preferably use size(d_recv, 1) for count
do i = 1,numnodes-1
  call MPI_IRecv(d_recv(1,i), natom*3, MPI_REAL8, i, tag(1,i), MPI_COMM_WORLD, &
       request_recv(i,1),ierr)
  if (ierr .ne. 0) call die('gather_nonbond/MPI_IRecv d_recv')
  call MPI_IRecv(E_recv(i), 3*2+1, MPI_REAL8, i, tag(2,i), MPI_COMM_WORLD, &
       request_recv(i,2),ierr)
  if (ierr .ne. 0) call die('gather_nonbond/MPI_IRecv E_recv')
  call MPI_IRecv(EQ_recv(1,i), reclength, MPI_REAL8, i, tag(3,i), MPI_COMM_WORLD, &
	request_recv(i,3),ierr)
  if (ierr .ne. 0) call die('gather_nonbond/MPI_IRecv EQ_recv')
end do

else                  !slave nodes
E_send%pp%el  = E_loc%pp%el
E_send%pp%vdw = E_loc%pp%vdw
E_send%pw%el  = E_loc%pw%el
E_send%pw%vdw = E_loc%pw%vdw
E_send%ww%el  = E_loc%ww%el
E_send%ww%vdw = E_loc%ww%vdw
E_send%lrf    = E_loc%lrf
do i=1,nstates
EQ_send(i)%qp%el  = EQ_loc(i)%qp%el
EQ_send(i)%qp%vdw = EQ_loc(i)%qp%vdw
EQ_send(i)%qw%el  = EQ_loc(i)%qw%el
EQ_send(i)%qw%vdw = EQ_loc(i)%qw%vdw
end do

! See comments above on the IRecv part
call MPI_Send(d, natom*3, MPI_REAL8, 0, tag(1,nodeid), MPI_COMM_WORLD,ierr) 
if (ierr .ne. 0) call die('gather_nonbond/Send d')
call MPI_Send(E_send, 3*2+1, MPI_REAL8, 0, tag(2,nodeid), MPI_COMM_WORLD,ierr) 
if (ierr .ne. 0) call die('gather_nonbond/Send E_send')
call MPI_Send(EQ_send, reclength, MPI_REAL8, 0, tag(3,nodeid), MPI_COMM_WORLD,ierr) 
if (ierr .ne. 0) call die('gather_nonbond/Send EQ_send')

end if
end subroutine gather_nonbond

#endif
!----------------------------------------------------------------------------------------

subroutine calculate_exc
! nasty piece of work to only calculate the interactions between the excluded atoms
! and subsequently substract them from the previously calculate energies of the whole system
! better than having n2 if cases during the actual nb calculations
! locals
integer                         :: gcnum,istate,aindex,num
integer                         :: iq,jq,i,j,ia,ic
TYPE(qr_dist3)                  :: dist
TYPE(qr_vec)                    :: shift
TYPE(ENERET_TYPE)               :: nb_ene
TYPE(OQ_ENERGIES)               :: EQ_loc(size(EQ))


do aindex = 2, ene_header%arrays - QCP_N
do istate = 1, nstates

! first set total energies
EQ_loc(istate) = EQ(istate)
! get actual index of calculation type
gcnum = ene_header%gcnum(aindex)

! start looping over all the excluded interactions
! first qq stuff, substract from qq index
do num = 1, exc_nbqq(gcnum,istate)%inter
iq = exc_nbqq(gcnum,istate)%list(num)%iq
jq = exc_nbqq(gcnum,istate)%list(num)%jq

nb_ene = nbe_qq(exc_nbqq(gcnum,istate)%list(num),EQ(istate)%lambda)
select case (ene_header%types(aindex))
case (FULL)
        EQ_loc(istate)%qq%el  = EQ_loc(istate)%qq%el  - nb_ene%Vel
        EQ_loc(istate)%qq%vdw = EQ_loc(istate)%qq%vdw - (nb_ene%V_a - nb_ene%V_b)
case (ELECTRO)
        EQ_loc(istate)%qq%el  = EQ_loc(istate)%qq%el  - nb_ene%Vel
case (VDW)
        EQ_loc(istate)%qq%vdw = EQ_loc(istate)%qq%vdw - (nb_ene%V_a - nb_ene%V_b)
end select
end do ! num
do num = 1, exc_nbqqp(gcnum,istate)%inter
iq = exc_nbqqp(gcnum,istate)%list(num)%i
i  = iqseq(iq)
j  = exc_nbqqp(gcnum,istate)%list(num)%j

dist = q_dist3(x(i),x(j))
nb_ene = nbe_qx(exc_nbqqp(gcnum,istate)%list(num),EQ(istate)%lambda,dist)
select case (ene_header%types(aindex))
case (FULL)
        EQ_loc(istate)%qp%el  = EQ_loc(istate)%qp%el  - nb_ene%Vel
        EQ_loc(istate)%qp%vdw = EQ_loc(istate)%qp%vdw - (nb_ene%V_a - nb_ene%V_b)
case (ELECTRO)
        EQ_loc(istate)%qp%el  = EQ_loc(istate)%qp%el  - nb_ene%Vel
case (VDW)
        EQ_loc(istate)%qp%vdw = EQ_loc(istate)%qp%vdw - (nb_ene%V_a - nb_ene%V_b)
end select
end do ! num
do num = 1 , exc_nbqp(gcnum,istate)%inter
iq = exc_nbqp(gcnum,istate)%list(num)%i
i  = iqseq(iq)
j  = exc_nbqp(gcnum,istate)%list(num)%j

if (use_PBC) then
        ia = iwhich_cgp(j)
        ic = cgpatom(cgp(ia)%first)
        shift = x(ic) - x(qswitch)
        dist  = q_dist3(x(i)-x(j),boxlength*q_nint(shift*inv_boxl))
else
        dist = q_dist3(x(i),x(j))
end if
nb_ene = nbe_qx(exc_nbqp(gcnum,istate)%list(num),EQ(istate)%lambda,dist)
select case (ene_header%types(aindex))
case (FULL)
        EQ_loc(istate)%qp%el  = EQ_loc(istate)%qp%el  - nb_ene%Vel
        EQ_loc(istate)%qp%vdw = EQ_loc(istate)%qp%vdw - (nb_ene%V_a - nb_ene%V_b)
case (ELECTRO)
        EQ_loc(istate)%qp%el  = EQ_loc(istate)%qp%el  - nb_ene%Vel
case (VDW)
        EQ_loc(istate)%qp%vdw = EQ_loc(istate)%qp%vdw - (nb_ene%V_a - nb_ene%V_b)
end select
end do ! num

! update total EQ for this state and index
EQ_loc(istate)%qx%el  = EQ_loc(istate)%qq%el  + EQ_loc(istate)%qp%el  + EQ_loc(istate)%qw%el
EQ_loc(istate)%qx%vdw = EQ_loc(istate)%qq%vdw + EQ_loc(istate)%qp%vdw + EQ_loc(istate)%qw%vdw
EQ_loc(istate)%total  = EQ_loc(istate)%q%bond + EQ_loc(istate)%q%angle + EQ_loc(istate)%q%torsion + &
        EQ_loc(istate)%q%improper + EQ_loc(istate)%qx%el + EQ_loc(istate)%qx%vdw + EQ_loc(istate)%restraint

end do ! nstates

call qatom_savetowrite(EQ_loc,aindex)

end do ! aindex

end subroutine calculate_exc

end module NONBONDENE




