!	(C) 2000 Uppsala Molekylmekaniska HB, Uppsala, Sweden
!	exc.f90
!	by Paul Bauer, based on trj.f90 file by John Marelius
!	Q subroutine to allow energy calculations excluding 
!	defined groups from the full energy
!	to get the real contribution to the energy profile

module EXC

use ATOM_MASK
use MISC
use QATOM
implicit none

	character(*), private, parameter	:: MODULE_VERSION = '0.01'
	character(*), private, parameter	:: MODULE_DATE    = '2014-04-21'

	integer, private						::	mask_rows=0
	integer,parameter						:: 	max_rows=20
	character*80							::	mask_row(max_rows)

		
	!mask 	
	integer, private							::	ncoords
	type ::EXC_PARAMS
		character*80		:: seltype
		integer			:: count
		integer			:: curcount=0
		character,allocatable	:: maskarray(:)
		type(MASK_TYPE)		:: gcmask
		integer			:: fileunit
		character		:: filename
	end type EXC_PARAMS

	type(EXC_PARAMS),allocatable		:: ST_gc(:)

	integer,allocatable			:: gc_mask_array(:)
!	integer,allocatable			:: gc_mask_map(:)



contains


subroutine gc_alloc_array(atoms)
	!arguments
	integer					:: atoms,alloc_status

	allocate(gc_mask_array(atoms),stat=alloc_status)
	call check_alloc_general(alloc_status,'gc mmain array')
	gc_mask_array(:)=0
	

end subroutine gc_alloc_array

logical function gc_store_mask(ST_gc_local,line)
	!arguments
	type(EXC_PARAMS)			::	ST_gc_local
	character(*)				::	line
	
	if(ST_gc_local%curcount < ST_gc_local%count) then
		ST_gc_local%curcount = ST_gc_local%curcount + 1
		ST_gc_local%maskarray(ST_gc_local%curcount) = line
		gc_store_mask = .true.
	else
		write(*,900) ST_gc_local%count
900		format('>>>>> ERROR: Too many atom mask rows (max',i3,').')
		gc_store_mask = .false.
	end if
end function gc_store_mask


subroutine gc_make_array(ST_gc_local)
	!arguments
	type(EXC_PARAMS)				:: ST_gc_local
	!locals
	integer						:: i,j,length,stat,alloc_status
	integer,parameter				:: max_length=80
	character*80,allocatable			:: tmp(:)
	character*80					:: str1,str2
!	type(MASK_TYPE)				        ::      masks
!Here we make the real magic, telling the program what atoms we want
!to have in the different groups to be excluded
!        call mask_initialize(masks)

	allocate(tmp(ST_gc_local%count),stat=alloc_status)
	call check_alloc_general(alloc_status,'gc tmp mask arrays')

	do i=1,ST_gc_local%count
!to prevent overflow
	length=(len_trim(ST_gc_local%seltype)+len_trim(ST_gc_local%maskarray(i))+1)
	if ( length.lt.max_length) then
		str1 = trim(ST_gc_local%seltype)//' '
		write(str2,'(i10)') ST_gc_local%maskarray(i)
		tmp(i) = str1//trim(str2)
	end if
	end do
	do i=1,ST_gc_local%count
		stat=mask_add(ST_gc_local%gcmask,tmp(i))
	end do

!	call mask_initialize(ST_gc_local%gcmask)

	do i=1,nat_pro
		if (ST_gc_local%gcmask%MASK(i)) then
			gc_mask_array(i)=1
		end if
	end do	
	deallocate(tmp, stat=alloc_status)
end subroutine gc_make_array		
		
subroutine set_gc_energies(atom,ngroups,Vel,VvdW,EQtot,EQ_exc,ST_local)
	!parameters
	integer				:: atom, ngroups
	real*8				:: Vel, VvdW
	type(EXC_PARAMS)		:: ST_local(:)
	type(Q_ENERGIES)		:: EQtot
	type(Q_ENERGIES)		:: EQ_exc(:)
	!locals
	integer				:: iii,jjj,maskatoms
        do iii=1,ngroups
                EQ_exc(iii)%qp%el         =       EQtot%qp%el
                EQ_exc(iii)%qp%vdw        =       EQtot%qp%vdw
                if (gc_mask_array(atom).eq.1) then
                        maskatoms = ST_local(iii)%gcmask%included
                        if(ST_local(iii)%gcmask%MASK(atom)) then
                                EQ_exc(iii)%qp%el         = EQ_exc(iii)%qp%el  - Vel
                                EQ_exc(iii)%qp%vdw        = EQ_exc(iii)%qp%vdw - VvdW
                        end if
                end if
        end do
end subroutine set_gc_energies

end module EXC

