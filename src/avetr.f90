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

module AVETR
  use PREP
  implicit none

	integer, parameter             :: AVE_PDB = 11
	integer(4), private            :: ncoords, N_sets = 0
! changed to use q vector structure to store coordinates
! one element has all three coordinates
	TYPE(qr_vec), allocatable, private  :: x_in(:), x_sum(:), x2_sum(:)
	real(kind=prec), private               :: rmsd
contains
!TODO: *choose which frames, add more trajectories, divide x_sum every 100 steps

!******************************************************
!Main subroutine
!******************************************************
subroutine avetr_calc
  integer :: allocation_status
  character(len=1) :: ans
  logical :: fin
  N_sets = 0
  call trajectory
! number of coordinates is still natom*3
! but store is now size of natom so we can
! access everything by the atom index instead of 
! having to play around with numbers
  ncoords = trj_get_ncoords() 
  allocate(x_in(ncoords/3), x_sum(ncoords/3), x2_sum(ncoords/3), &
				stat=allocation_status)
  if (allocation_status .ne. 0) then
    write(*,*) 'Out of memory!'
    return
  end if
  do while(trj_read_masked(x_in))  !add from first file
    call add_coordinates
  end do

  !add from multiple files
  fin = .false.
  do while(.not. fin)
    CALL get_string_arg(ans, '-----> Add more frames? (y or n): ')
    if (ans .eq. 'y') then
	  call trajectory
	  do while(trj_read_masked(x_in))  !add from additional files
		call add_coordinates
	  end do
	else
	  fin = .true.
	end if
  end do


  call average
  call write_average
  deallocate(x_in, x_sum, x2_sum, stat=allocation_status)
end subroutine avetr_calc

!******************************************************
!Sum the coordinates and the sqared coordinates
!******************************************************
subroutine add_coordinates
	x_sum = x_sum + x_in
	x2_sum = x2_sum + (x_in*x_in)
	N_sets = N_sets +1
end subroutine add_coordinates

!******************************************************
!Make average and rmsd
!******************************************************
subroutine average
TYPE(qr_vec)    :: temp
integer         :: i
	x_sum = x_sum / real(N_sets,kind=prec)
	x2_sum = x2_sum / real(N_sets,kind=prec)
        temp = temp * zero
        do i=1, ncoords/3
        temp = x2_sum(i) - (x_sum(i)*x_sum(i))
        end do
        temp = temp / real(ncoords,kind=prec)
        rmsd = q_sqrt(temp%x+temp%y+temp%z)
end subroutine average

!******************************************************
!Write average coords to pdb file.
!Variables used from prep: mask
!Variables used from topo: xtop
!******************************************************
subroutine write_average
	!assign masked coordinates to right atom in topology
	call mask_put(mask, xtop, x_sum)
    call writepdb
	write(*,'(a,f6.3,a)') 'Root mean square co-ordinate deviation ', rmsd, ' A'
	x_sum  = x_sum * zero
	x2_sum = x2_sum * zero
end subroutine write_average

end module AVETR
