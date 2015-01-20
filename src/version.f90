! (C) 2014 Uppsala Molekylmekaniska HB, Uppsala, Sweden
! version.f90
! Initial date: 2015
! by Ireneusz Szeler
! Q version and halp print info

module VERSION

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
  integer  :: num_arg, i
  integer  :: datum(8)
  character(200) :: flag

!  if (nodeid .eq. 0) then
    
    fin = .false.
    num_arg = command_argument_count()
    if (num_arg .gt. 0) then
      call getarg(1,flag)
      call lowcase(flag)
      if ( flag .eq. '-v' .or. flag .eq. '-version' .or. flag .eq. '--version') fin = .true. 
    end if

! start-of-header
    write (*,'(79a)') ('#',i=1,79)
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
    
    if ( flag .eq. '-h' .or. flag .eq. '-help' .or. flag .eq. '--help') then
       fin = .true.
    end if
  
    if(fin) STOP


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

end module VERSION
