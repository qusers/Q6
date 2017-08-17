! Q6: A comprehensive simulation package for molecular dynamics simulations and 
! free energy calculations, including empirical valence bond simulations, 
! linear interaction energy calculations, and free energy perturbation.
! 
! Copyright © 2017 Johan Åqvist, John Marelius, Shina Caroline Lynn Kamerlin and Paul Bauer
! 
! This program is free software; you can redistribute it and/or modify it under the 
! terms of the GNU General Public License as published by the Free 
! Software Foundation; either version 2 of the License, or any later version.
! 
! This program is distributed in the hope that it will be useful, 
! but WITHOUT ANY WARRANTY; without even the implied warranty of 
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
! See the GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License along with 
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin 
! Street, Fifth Floor, Boston, MA  02110-1301, USA. Also add information on 
! how to contact you by electronic and paper mail.
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

	integer,parameter						:: 	max_rows=20
	character*80							::	mask_row(max_rows)

		
	!mask 	
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
	integer						:: i,length
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
        character*10                                    :: str1
	str1 = trim(ST_gc_local%maskarray(nres))
	read(str1,*) get_from_mask
end function get_from_mask

subroutine excluded_shutdown(nexc)
!arguments
integer                 :: nexc
!locals
integer                 :: i
do i=1,nexc
        call mask_deallocate(ST_gc(i)%gcmask)
end do
if(allocated(exc_nbqq))  deallocate(exc_nbqq)
if(allocated(exc_nbqp))  deallocate(exc_nbqp)
if(allocated(exc_nbqqp)) deallocate(exc_nbqqp)

end subroutine excluded_shutdown

end module EXC

