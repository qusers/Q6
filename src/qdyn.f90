! (C) 2014 Uppsala Molekylmekaniska HB, Uppsala, Sweden
! qdyn.f90
! Initial date: 2000
! by Johan Åqvist, John Marelius, Anders Kaplan & Martin Nervall
! Qdyn molecular dynamics main program

program Qdyn5
  use MD
  use VERSIONS
  use MPIGLOB ! use MPI global data
#if defined (_DF_VERSION_)
  use dfport  ! portability lib for signals, used by Compaq Visual Fortran compiler
#endif
!$ use omp_lib

  implicit none

  ! version data
  character(10)					:: QDYN_NAME = 'Qdyn'
  character(10)					:: QDYN_VERSION = '5.06'
  character(12)					:: QDYN_DATE    = '2014-03-14'
#if defined (USE_MPI)
  character(10)					:: QDYN_SUFFIX = '_parallel'
#else
  character(10)					:: QDYN_SUFFIX = ''
#endif


#if defined (USE_MPI)
  ! MPI error code
  integer						:: qdyn_ierr
#endif

  ! signal handler data and declarations
  integer(4)					:: sigret
#if defined (_DF_VERSION_)
  ! nothing
#else
  integer(4), parameter			:: SIGINT  = 2 ! CTRL-C signal
  integer(4), parameter			:: SIGKILL = 9 ! kill/CTRL-BREAK signal
  integer(4), parameter			:: SIGABRT = 6 ! kill/CTRL-BREAK signal
#endif
  external sigint_handler
  external sigkill_handler
  external sigabrt_handler

#if defined (USE_MPI)
  ! initialize MPI
  call MPI_Init(qdyn_ierr)
  if (qdyn_ierr .ne. MPI_SUCCESS) call die('failure at MPI init')
  call MPI_Comm_rank(MPI_COMM_WORLD, nodeid, qdyn_ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, numnodes, qdyn_ierr)
#else
	nodeid = 0
	numnodes = 1
#endif

  ! initialize signal handlers
  sigret = qsignal(SIGINT, sigint_handler, -1_4)
  sigret = qsignal(SIGKILL, sigkill_handler, -1_4)
  sigret = qsignal(SIGABRT, sigabrt_handler, -1_4)

  ! initialise static data, display banner etc
  call startup

  if (nodeid .eq. 0) then
#ifdef _OPENMP
!$omp parallel
  threads_num = omp_get_num_threads()
!$omp end parallel
#endif

	! master node: read input and initialise

	if(.not. initialize()) call die('Invalid data in input file')						! read input data
	call open_files						! open necessary files
	call topology				! read topology
	call prep_coord						! read coords, solvates etc
	if ( nstates > 0 ) call get_fep	! read fep/evb strategy
	!remove things with code 0 and maybe excluded bonded interactions

	call prep_sim						! prepare for simulation (calc. inv. mass, total charge,...)
	call close_input_files				! close input files

	call init_shake
	call make_nbqqlist
	call shrink_topology
	call nbmonitorlist
	call init_trj

    ! generate Maxwellian velocities and shake initial x and v if necessary
    if ( iseed > 0 ) then
      call maxwell
	  call initial_shaking
 	  call stop_cm_translation
    end if
  end if

#if defined (USE_MPI)
  ! initialise slave nodes
  if (numnodes .gt. 1) call init_nodes
#endif
	! count non-bonded pairs to get the maximum number, then distribute them among the nodes

call distribute_nonbonds

 ! do the work!
  call md_run

  if (nodeid .eq. 0) then
	! master node: close output files
	call close_output_files
  end if

  ! deallocate memory etc.
  call shutdown

#if defined (USE_MPI)
  ! shut down MPI
  call MPI_Finalize(qdyn_ierr)
#endif

contains

!-----------------------------------------------------------------------

! startup/shutdown code

subroutine startup

  if (nodeid .eq. 0) then
    call version_check(QDYN_NAME, QDYN_VERSION, QDYN_DATE, QDYN_SUFFIX) ! print version and chack for flags
  end if

  ! initialise used modules
  call md_startup

end subroutine startup

!-----------------------------------------------------------------------

subroutine shutdown
  integer						:: i

  if (nodeid .eq. 0) then
#if defined (DUM)
	write(*,*) 'QDum input checker version ', QDYN_VERSION, ' terminated normally.'
#else
	write(*,*) 'QDyn version ', trim(QDYN_VERSION), trim(QDYN_SUFFIX), ' terminated normally.'
#endif
	write (*,'(79a)') ('#',i=1,79)
  end if

  ! call md's shutdown
  call md_shutdown
end subroutine shutdown

!-----------------------------------------------------------------------
INTEGER(4) FUNCTION qsignal( signum, proc, sigflag )
!       use MD
       implicit none
       INTEGER(4)                               :: signum, sigflag
       external proc
       qsignal = 1
END FUNCTION qsignal

end program Qdyn5

! signal handlers

INTEGER(4) FUNCTION sigint_handler(sig_num)
  use MD
  implicit none
  INTEGER(4)					:: sig_num

  call die('user request (control-C)')
  sigint_handler = 1
END FUNCTION sigint_handler

INTEGER(4) FUNCTION sigkill_handler(sig_num)
  use MD
  implicit none
  INTEGER(4)					:: sig_num

  call die('kill signal')
  sigkill_handler = 1
END FUNCTION sigkill_handler

INTEGER(4) FUNCTION sigabrt_handler(sig_num)
  use MD
  implicit none
  INTEGER(4)					:: sig_num

  call die('kill signal')
  sigabrt_handler = 1
END FUNCTION sigabrt_handler
