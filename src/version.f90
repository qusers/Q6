! (C) 2014 Uppsala Molekylmekaniska HB, Uppsala, Sweden
! version.f90
! Initial date: 2015
! by Ireneusz Szeler
! Q version and halp print info

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

    call version_print(Q_PROGRAM, Q_VERSION, Q_DATE, Q_SUFFIX)
    
    call lowcase(Q_PROGRAM)
    
    if ( flag .eq. '-h' .or. flag .eq. '-help' .or. flag .eq. '--help') then
       fin = .true.
       write(*,*) 
       write(*,'(a,a)') trim(Q_PROGRAM), ' help information'
       select case (Q_PROGRAM)
         case ('qdyn')
             write(*,*) 
             write(*,'(a)') 'To run calculations use: '
             write(*,'(a,a,a)') '    ',trim(Q_PROGRAM), '5 inputfile.inp > output.file'
             write(*,'(a)') ' or for parallel version'
             write(*,'(a,a,a)') '    mpienvirment ', trim(Q_PROGRAM), '5p inputfile.inp > output.file'
             write(*,*) 
             write(*,'(a)') 'where:'
             write(*,'(a)') 'mpienvirment - e.q. mpirun -n 4, for more info check cluster informations '
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
    case default
      write(*,*) 
      write(*,'(a)') 'Welcome in ...' 
      write(*,*) 
  end select


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

subroutine version_print(Q_PROGRAM, Q_VERSION, Q_DATE, Q_SUFFIX)

! arguments
  character(*)  :: Q_PROGRAM
  character(*)  :: Q_VERSION
  character(*)  :: Q_DATE
  character(*)  :: Q_SUFFIX

! local
  integer  :: i
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
#else
  write(*,'(a,a,a,a,a)')  trim(Q_PROGRAM), ' version ', trim(Q_VERSION), trim(Q_SUFFIX),' initialising'
#endif
  call date_and_time(values=datum)
  write(*,130) datum(1),datum(2),datum(3),datum(5),datum(6),datum(7)
130	format('Current date ',i4,'-',i2,'-',i2,' and time ',i2,':',i2,':',i2)

end subroutine version_print

end module VERSIONS
