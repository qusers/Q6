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
!  calc_rms.f90
!  by John Marelius
!  root mean square coordinate deviation
module calc_rms
  use calc_base
  use maskmanip
  implicit none

!constants
  integer, parameter              :: MAX_MASKS = 10

!module variables
  type(MASK_TYPE), private, target :: masks(MAX_MASKS)
  integer, private                 :: Nmasks = 0
  type RMS_COORD_TYPE
    real(kind=prec), pointer       :: x(:), x0(:)
  end type RMS_COORD_TYPE
  type(RMS_COORD_TYPE), private    :: coords(MAX_MASKS)
  character(80)                    :: rms_basename(MAX_MASKS)
contains

subroutine rms_initialize
end subroutine rms_initialize

subroutine rms_finalize(i)
  integer                         :: i
  call mask_finalize(masks(i))
end subroutine rms_finalize

integer function rms_add(desc)
  !arguments
  character(*)                    :: desc
  integer                         :: ats,fn
  if(Nmasks == MAX_MASKS) then
    write(*,10) MAX_MASKS
    return
  end if
10      format('Sorry, the maximum number of RMSD calculations is ',i2)
  !add a new RMSD mask
  Nmasks = Nmasks + 1
  call mask_initialize(masks(Nmasks))
  ats =  maskmanip_make(masks(Nmasks))
  !discard if no atoms in mask
  if(ats == 0) then
    call mask_finalize(masks(Nmasks))
    Nmasks = Nmasks - 1
    RMS_add = 0
    return
  end if

  allocate(coords(Nmasks)%x(3*ats), coords(Nmasks)%x0(3*ats))

  call rms_make_ref(Nmasks)
  RMS_add = Nmasks
  ! get name for file print out
  call getlin(rms_basename(Nmasks),'--> RMSD print out filename: ')
  ! open and close file to overwrite current saved data
  ! this needs to be changed
  fn = freefile()
  open(fn,file=rms_basename(Nmasks),status='replace',action='write')
21      format('RMSD calculation ',i6)
  write(fn,20) masks(Nmasks)%included
  write(fn,21) Nmasks
  close(fn)
  write(desc, 20) masks(Nmasks)%included
20      format('Root Mean Square Deviation for ',i6,' atoms')
 end function rms_add


subroutine rms_calc(i)
  !arguments
  integer, intent(in)             :: i

  !locals
  real(kind=prec)                 :: r
  integer                         :: fn,ats

  fn = freefile()
  if(i < 1 .or. i > Nmasks) return
  open(fn,file=rms_basename(i),status='old',action='write',position='append')
  ats = masks(i)%included
  allocate(tmp_coords(ats))
  tmp_coords(:) = translate_coords(coords(i)%x)
  call mask_get(masks(i), translate_coords(xin), tmp_coords)
  coords(i)%x = btranslate_coords(tmp_coords)
  deallocate(tmp_coords)
  r = sqrt(  sum((coords(i)%x-coords(i)%x0)**2)/(masks(i)%included) )  ! removed 3 in front of mask(i)
  write(*,100, advance='no') r
  write(fn, 100, advance='yes') r
  close(fn)
100 format(f8.6)
end subroutine rms_calc


subroutine rms_make_ref(i)
  integer                         :: i,at,ats
  if(i < 1 .or. i > Nmasks) return
  ats = masks(i)%included
  allocate(tmp_coords(ats))
  tmp_coords(:) = translate_coords(coords(i)%x0)
  call mask_get(masks(i), xtop, tmp_coords)
  coords(i)%x0 = btranslate_coords(tmp_coords)
  deallocate(tmp_coords)
end subroutine rms_make_ref


subroutine rms_heading(i)
  integer                         :: i
  write(*,'(a)', advance='no') ' RMSD(A)'
end subroutine rms_heading


end module calc_rms
