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
use GLOBALS

implicit none

	character(*), private, parameter	:: MODULE_VERSION = '0.01'
	character(*), private, parameter	:: MODULE_DATE    = '2014-04-21'

	integer, private						::	mask_rows=0
	integer,parameter						:: 	max_rows=20
	character*80							::	mask_row(max_rows)

		
	!mask 	
	integer, private							::	ncoords
        logical,allocatable             :: tempmask(:,:)
contains


logical function gc_store_mask(ST_gc_local,line)!,p_mask)
	!arguments
	type(EXC_PARAMS)        	::	ST_gc_local
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
	type(EXC_PARAMS)			:: ST_gc_local
	!locals
	integer						:: i,j,length,stat,alloc_status
	integer,parameter				:: max_length=80
	character*80					:: str1,str2

	allocate(ST_gc_local%sendtomask(ST_gc_local%count))
	do i=1,ST_gc_local%count
!to prevent overflow
	length=(len_trim(ST_gc_local%seltype)+len_trim(ST_gc_local%maskarray(i))+1)
	if ( length.lt.max_length) then
		str1 = trim(ST_gc_local%seltype)//' '
		str2 = trim(ST_gc_local%maskarray(i))
                if ( trim(str1) == 'atom' ) then
                        ST_gc_local%sendtomask(i) = trim(str2)
                else
		        ST_gc_local%sendtomask(i) = trim(str1)//' '//trim(str2)
                end if
	end if
	end do
end subroutine gc_make_array

integer function get_from_mask(ST_gc_local,nres)
	!arguments
	type(EXC_PARAMS)                        :: ST_gc_local
	integer					:: nres
	!locals
!	        integer                                         :: i,j,length,
!        integer,parameter                               :: max_length=80
        character*10                                    :: str1,str2
	str1 = trim(ST_gc_local%maskarray(nres))
	read(str1,*) get_from_mask
end function get_from_mask

subroutine set_gc_energies(atomi,atomj,Vel,VvdW,totVel,totVvdW,mask,gctype)
	!parameters
	integer				:: atomi,atomj,gctype
	real(kind=prec)				:: Vel, VvdW
        logical(1),pointer                 :: mask(:)
	real(kind=prec)			        :: totVel,totVvdw
	
        if(((mask(atomi)).or.(mask(atomj)))) then
		select case (gctype)
			case (FULL)
			!Do nothing
			case (ELECTRO)
			totVvdW = totVvdW + VvdW
			case (VDW)
			totVel = totVel  + Vel
		end select
	else
                totVel = totVel  + Vel
                totVvdW = totVvdW + VvdW
        end if
end subroutine set_gc_energies
subroutine	excluded_shutdown(nexc)
!arguments
integer			:: nexc
!locals
integer			:: i
do i=1,nexc
        call mask_deallocate(ST_gc(i)%gcmask)
end do
if(allocated(exc_nbqq_list)) then
        deallocate(exc_nbqq_list)
end if
if(allocated(exc_nbqp_list)) then
        deallocate(exc_nbqp_list)
end if

end subroutine excluded_shutdown

end module EXC

