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
! version.f90
! Initial date: 2015
! by Ireneusz Szeler
! Q version and help print info

module VERSIONS

  implicit none

contains

subroutine version_check(Q_PROGRAM, Q_VERSION, Q_DATE, Q_SUFFIX)

! arguments
  character(*)  :: Q_PROGRAM
  character(*)  :: Q_VERSION
  character(*)  :: Q_DATE
  character(*)  :: Q_SUFFIX

! local
  logical  :: fin
  integer  :: num_arg
  character(200) :: flag

!  if (nodeid .eq. 0) then
    
    fin = .false.
    num_arg = command_argument_count()
    if (num_arg .gt. 0) then
      call getarg(1,flag)
      call lowcase(flag)
      if ( flag .eq. '-v' .or. flag .eq. '-version' .or. flag .eq. '--version') fin = .true. 
    end if

    call version_print(Q_PROGRAM, Q_VERSION, Q_SUFFIX)
    
    call lowcase(Q_PROGRAM)
    
    if ( flag .eq. '-h' .or. flag .eq. '-help' .or. flag .eq. '--help') then
       fin = .true.
       write(*,*) 
       write(*,'(a,a)') trim(Q_PROGRAM), ' help information'
       write(*,*) 'Q6: A comprehensive simulation package for molecular dynamics simulations and'
       write(*,*) 'free energy calculations, including empirical valence bond simulations,'
       write(*,*) 'linear interaction energy calculations, and free energy perturbation.'
       write(*,*) 
       write(*,*) 'Copyright © 2017 Johan Åqvist, John Marelius, Shina Caroline Lynn Kamerlin and Paul Bauer'
       write(*,*) 
       write(*,*) 'This program is free software; you can redistribute it and/or modify it under the'
       write(*,*) 'terms of the GNU General Public License as published by the Free'
       write(*,*) 'Software Foundation; either version 2 of the License, or any later version.'
       write(*,*) 
       write(*,*) 'This program is distributed in the hope that it will be useful,'
       write(*,*) 'but WITHOUT ANY WARRANTY; without even the implied warranty of'
       write(*,*) 'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.'
       write(*,*) 'See the GNU General Public License for more details.'
       write(*,*) 
       write(*,*) 'You should have received a copy of the GNU General Public License along with'
       write(*,*) 'this program; if not, write to the Free Software Foundation, Inc., 51 Franklin'
       write(*,*) 'Street, Fifth Floor, Boston, MA  02110-1301, USA. Also add information on'
       write(*,*) 'how to contact you by electronic and paper mail.'
       write(*,*) 
       write(*,*) 'Information about copying and warrantry can be found in the file LICENCE'
       write(*,*) 
       select case (Q_PROGRAM)
         case ('qdyn')
             write(*,*) 
             write(*,'(a)') 'To run calculations use: '
             write(*,'(a,a,a)') '    ',trim(Q_PROGRAM), '6 inputfile.inp > output.file'
             write(*,'(a)') ' or for parallel version'
             write(*,'(a,a,a)') '    mpienvironment ', trim(Q_PROGRAM), '6p inputfile.inp > output.file'
             write(*,*) 
             write(*,'(a)') 'where:'
             write(*,'(a)') 'mpienvironment - e.q. mpirun -n 4, for more info check cluster informations '
             write(*,*) 
         case ('qfep')
             write(*,*) 
             write(*,'(a,a)') 'In this moment no info available for ', trim(Q_PROGRAM) 
             write(*,*) 
         case ('qprep')
             write(*,*) 
             write(*,'(a,a)') 'In this moment no info available for ', trim(Q_PROGRAM) 
             write(*,*) 
         case ('qcalc')
             write(*,*) 
             write(*,'(a,a)') 'In this moment no info available for ', trim(Q_PROGRAM) 
             write(*,*) 
         case ('qpi')
              write(*,*)
              write(*,'(a,a)') 'At this moment no info available for ', trim(Q_PROGRAM)
         case default
             write(*,*) 
             write(*,'(a)') 'It is some new program added to Q?'
             write(*,*) 
         end select
    end if
  
    if(fin) STOP

  select case (Q_PROGRAM)
    case ('qdyn')
      write(*,*) 
      write(*,'(a,a,a,a)') 'Welcome in ', trim(Q_PROGRAM), ' modification date ', trim(Q_DATE)
      write(*,*)
    case ('qfep')
      write(*,*) 
      write(*,'(a,a,a,a)') 'Welcome in ', trim(Q_PROGRAM), ' modification date ', trim(Q_DATE) 
      write(*,*) 
    case ('qprep')
      write(*,*) 
      write(*,'(a,a,a,a)') 'Welcome in ', trim(Q_PROGRAM), ' modification date ', trim(Q_DATE) 
      write(*,*) 
    case ('qcalc')
      write(*,*) 
      write(*,'(a,a,a,a)') 'Welcome in ', trim(Q_PROGRAM), ' modification date ', trim(Q_DATE) 
      write(*,*)
    case ('qpi')
      write(*,*)
      write(*,'(a,a,a,a)') 'Welcome in ', trim(Q_PROGRAM), ' modification date ', trim(Q_DATE)
    case default
      write(*,*) 
      write(*,'(a)') 'Welcome in ...' 
      write(*,*) 
  end select
  write(*,*) 'Q6, Copyright © 2017 Johan Åqvist, John Marelius, Shina Caroline Lynn Kamerlin and Paul Bauer'
  write(*,*) 'Q6 comes with ABSOLUTELY NO WARRANTY.  This is free software, and you are welcome'
  write(*,*) 'to redistribute it under certain conditions. For details, add the --help flag.'


!  end if ! node .eq. 0

end subroutine version_check

elemental subroutine lowcase(word)

 character(*), intent(in out) :: word
 integer   :: i, ic, nlen

 nlen = len(word) 
 do i=1, nlen
   ic = ichar(word(i:i))
   if (ic >= 65 .and. ic < 90) word(i:i) = char(ic+32)
 end do

end subroutine lowcase

subroutine version_print(Q_PROGRAM, Q_VERSION, Q_SUFFIX)

! arguments
  character(*)  :: Q_PROGRAM
  character(*)  :: Q_VERSION
  character(*)  :: Q_SUFFIX

! local
  integer  :: datum(8)

! start-of-header
  write(*,*)
  write(*,'(a)') 'Build and version information'
  write(*,*)
!  write (*,'(79a)') ('#',i=1,79)
  if ( Q_PROGRAM .eq. 'Qdyn') then
#if defined (DUM)
    write(*,'(a,a,a)') 'QDum input checker version ', trim(Q_VERSION), ' initialising'
#elif defined(EVAL)
    write(*,'(a,a,a)') 'QDyn evaluation version ', trim(Q_VERSION), ' initialising'
    write(*,'(a)') 'This version is for evaluation purposes only.'
    write(*,'(a)') 'Optimisations are disabled - runs at <20% of maximum speed.'
#endif
  end if
#if defined (BUILD_USERNAME) && defined (BUILD_HOSTNAME) && defined (BUILD_DATE) && defined (BUILD_SOURCE) && defined (BUILD_NUMBER) && defined(BUILD_COMPILER)
  write(*,'(a,a)') 'Build number ', BUILD_NUMBER
  write(*,'(a,a)') 'Build date   ', BUILD_DATE
  write(*,'(a)')   'Built:       '
  write(*,'(a,a)') '      by     ', BUILD_USERNAME
  write(*,'(a,a)') '      on     ', BUILD_HOSTNAME
  write(*,'(a,a)') '      git id ', BUILD_SOURCE
  write(*,'(a,a)') '      with   ', BUILD_COMPILER
#endif
  write(*,'(a,a,a,a,a)')  trim(Q_PROGRAM), ' version ', trim(Q_VERSION), trim(Q_SUFFIX),' initialising'
  call date_and_time(values=datum)
  write(*,130) datum(1),datum(2),datum(3),datum(5),datum(6),datum(7)
130	format('Current date ',i4,'-',i2,'-',i2,' and time ',i2,':',i2,':',i2)

end subroutine version_print

!new function to pass essential version information to other routines
character(80) function version_pass()
! local	
 

#if defined (BUILD_SOURCE) && defined (BUILD_NUMBER)
 version_pass = trim(BUILD_NUMBER)//' , git id='//trim(BUILD_SOURCE)
#else
 version_pass = trim(BUILD_NUMBER)
#endif

end function version_pass

character(80) function date_pass()
! local 

 date_pass = trim(BUILD_DATE)

end function date_pass

end module VERSIONS
