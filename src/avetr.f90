! (C) 2014 Uppsala Molekylmekaniska HB, Uppsala, Sweden
! average co-ordinates from Qdyn trajectory files and write pdb-structure
! Added to Qprep March 2004 by Martin Nervall
! Tested to reproduce average structures from vmd

module AVETR
  use PREP
  implicit none

	integer, parameter             :: AVE_PDB = 11
	integer(4), private            :: ncoords, N_sets = 0
	real(4), allocatable, private  :: x_in(:), x_sum(:), x2_sum(:)
	real(8), private               :: rmsd
contains
!TODO: *choose which frames, add more trajectories, divide x_sum every 100 steps

!******************************************************
!Main subroutine
!******************************************************
subroutine avetr_calc
  integer :: i, allocation_status
  character(len=1) :: ans
  logical :: fin
  N_sets = 0
  call trajectory
  ncoords = trj_get_ncoords()
  allocate(x_in(ncoords), x_sum(ncoords), x2_sum(ncoords), &
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
	x2_sum = x2_sum + x_in**2
	N_sets = N_sets +1
end subroutine add_coordinates

!******************************************************
!Make average and rmsd
!******************************************************
subroutine average
	x_sum = x_sum / N_sets
	x2_sum = x2_sum / N_sets
	rmsd = sqrt(sum(x2_sum - x_sum**2)/ncoords)
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
	x_sum = 0
	x2_sum = 0
end subroutine write_average

end module AVETR
