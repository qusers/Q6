!	(C) 2000 Uppsala Molekylmekaniska HB, Uppsala, Sweden

!	q_prep.f90
!	by Johan Åqvist & John Marelius

!	Qprep topology preparation main program

PROGRAM qprep5
	use PREP
	use AVETR
	IMPLICIT none

	character(*), parameter	::	PROGRAM_VERSION = '5.04'
	character(*), parameter	::	PROGRAM_DATE = '2005-04-14'
	logical                 ::  use_inputfile
	character(200)          ::  fileName=''
!.......................................................................

	call startup

!	allocate memory for libraries
	call topo_set_max(0, max_lib, 0)
	
!	reset residue count
	CALL clearlib

	use_inputfile = check_inputfile(fileName)
	if (use_inputfile) then
	  call qprep_from_inputfile(fileName)
	else
	  call qprep_from_commandline
	endif
		
	call shutdown

contains


!*************************************************************************
!*  Read input from file and execute commands
!*************************************************************************
subroutine qprep_from_inputfile(filename)

	character(200)	::	filename
	character(200)	::	command
	logical         ::  readable
	integer			:: INPF_U = 0
	integer			:: stat


!Fixa så att filen läses per line och parsas på samma sätt som med command line

	if (.not. parse_open_file(filename)) then
		write( * , '(a,a)') 'Could not open ',filename
		return
	endif

	do
		call get_string_arg(command)
!		if(command == '') then
!			return
!		endif

		select case (command)
		case('quit', 'q') 
			exit
		case default
			call parse_command(command)
		end select
	enddo

end subroutine qprep_from_inputfile

!*************************************************************************
!*  Read input from command line and execute commands
!*************************************************************************
subroutine qprep_from_commandline

	character(200)					::	command

	do
		call parse_reset
		call get_string_arg(command, 'Qprep> ')
		call locase(command) !make all lower case
		select case (command)
		case('quit', 'q') 
			exit
		case default
			call parse_command(command)
		end select
	enddo

end subroutine qprep_from_commandline

   
!*************************************************************************
!*  Parse a command and call corresponding subroutine
!*************************************************************************
subroutine parse_command(command)

character(*), intent(IN)  ::  command
! --- Command loop
!
		select case (command)
		case ('average', 'av')
			CALL avetr_calc
		case ('readlib', 'rl')
			CALL readlib
		case('clearlib', 'cl')
			CALL clearlib
		case('readpdb', 'rp')
			CALL readpdb
		case('readprm', 'readff', 'rff', 'rprm')
			CALL readff
		case('addbond', 'ab')
			CALL addbond
		case('clearbond', 'clearbonds')
			CALL clearbond
		case('maketop', 'mt')
			CALL maketop
		case('cleartop')
			CALL cleartop
		case('listseq', 'ls')
			CALL listseq
		case('listres', 'lr')
			CALL listres
		case('writetop', 'wt')
			CALL writetop
		case('checkbonds', 'cb')
			CALL checkbonds
		case('checkangs', 'ca')
			CALL checkangs
		case('checktors', 'ct')
			CALL checktors
		case('checkimps', 'ci')
			CALL checkimps
		case('changeimp')
			CALL changeimp
		case('readtop', 'rt')
			CALL readtop
		case('readx', 'rx')
			CALL readx
		case('makeshell', 'ms')
			CALL make_shell2
		case('mask', 'ma')
			CALL modify_mask
		case('trajectory', 'trj', 'tr')
			CALL trajectory
		case('readtrajectory', 'readframe','rf')
			CALL readframe
		case('readnext', 'rn')
			CALL readnext
		case('solvate', 'so')
			CALL solvate
		case('writepdb', 'wp')
			CALL writepdb 
		case('writemol2', 'wm')
			CALL writemol2
		case('xlink', 'crosslink', 'xl')
			CALL xlink
		case('prefs', 'preferences')
			CALL listprefs
		case('set', 's')
			CALL set
		case('help', '?', 'h')
			call help
		case('boundary', 'bc')
			call define_boundary_condition
		case default
			WRITE( * , '(/,a,a,a,/)') 'unrecognized command "', trim(command), '" (type ? for help)'
		END select

end subroutine parse_command


!-----------------------------------------------------------------------

subroutine help
		WRITE( *, * )
		WRITE( * , '(a)') &
		'command      argument            decription', &
		'----------- ------------------- -------------------------------------------------', &
		'addbond                          adds extra bonds(e.g. S-S)',&
		'average                          Calulates an average structure from a trajectory file.',&
		'boundary    [boundry condition]  set boundary condition',&
		'            [centre] ',&
		'            (box) [boxlengths]',&
		'            (sphere) [radius]',&
		'            (sphere) [inner radius]',&
		'changeimp                        redefine(specified) improper torsions',&
		'checkangs   [energy_threshold]   check angle energies',&
		'checkbonds  [energy_threshold]   check bonds energies',&
		'checkimps   [energy_threshold]   check improper torsion energies',&
		'checktors   [energy_threshold]   check torsion energies',&
		'clearbond                        clears extra bonds',&
		'clearlib                         unloads all libaries',&
		'cleartop                         clears topology & parameters',&
		'help                             shows command list',&
		'listseq                          lists the residue sequence',&
		'listres     [residue_number]     lists atoms in residue',&
		'makeshell                        fix the mask of the atoms in the restrained shell',&
		'maketop     [name]               generates the topology',&
		'mask        [mask_def|none]      add to or clear atom mask',&
		'preferences                      list preferences', &
		'quit                             quits the program',&
		'readframe   [frame]              reads coordinates from trajectory',&
		'readlib     [library_file]       reads library file',&
		'readnext                         reads next frame from trajectory', &
		'readpdb     [pdb_file]           reads pdb file',&
		'readprm     [param.file]         reads FF parameter file',&
		'readtop     [topology_file]      reads topology file',&
		'readx       [restart_file]       reads coord. file',&
		'set                              set preferences',&
		'solvate(boundry=sphere)          solvate sphere with specified options',&	
		'            [centre]',&
		'            [radius]',&
		'            [grid|file|restart] ',&
		'            [solvent name] ',&
		'            [file name] ',&
		'',&
		'solvate(boundry=box)             solvate box with specified options',&
		'            [grid|file|restart]  ',&	
		'            [solvent name] ',&
		'            [file name] ',&
		'',&
		'trajectory  [trajectory_file]    open trajectory file',&
		'writetop    [topology_file]      writes topology file',&
		'writepdb    [pdb_file]           writes pdb file', &
		'            [hydrogen_flag [gap_flag [water_flag]]]',&
		'',&                                
		'writemol2   [mol_no]             writes molecule mol_no (0 for all) to',&
		'            [mol2_file]          SYBYL mol2 file', &
		'            [hydrogen_flag ',& 
		'			  [water_flag]]',&
		'',&                                
		'xlink                            add crosslink bonds interactively', &
		'', &
		'short form meaning' ,&
		'---------- -------', &
		'ab          addbond', &
		'av          average', &
		'bc          boundary', &
		'cb,ca,ct,ci checkbonds, checkangs, checktors, checkimps', &
		'cl          clearlib', &
		'h, ?        help', &
		'lr          listres', &
		'ls          listseq', &
		'ma          mask', &
		'ma 0        clear mask', &
		'ms          makeshell', &
		'mt          maketop', &
		'prefs       preferences', &
		'q           quit', &
		'rf          readframe', &
		'rff         readprm', &
		'rl          readlib', &
		'rn          readnext', &
		'rp          readpdb', &
		'rt          readtop', &
		'rx          readx', &
		's           set',&
		'so          solvate', &
		'trj, tr     trajectory', &
		'wp          writepdb', &
		'wt          writetop', &
		'xl          xlink'

		WRITE( *, * )
end subroutine help

!-----------------------------------------------------------------------

subroutine startup
	integer						::	i

	write(*,'(79a)')('#',i=1,79)
	write(*,'(a,a)') 'Welcome to Qprep version ',PROGRAM_VERSION
	call prep_startup
	write(*,*)

end subroutine startup

!-----------------------------------------------------------------------

subroutine shutdown
	call prep_shutdown
	stop 'Qprep terminated'
end subroutine shutdown

!*************************************************************************
!*  Determine if Qprep is to be run from command line or from input file
!*************************************************************************

logical function check_inputfile(infilename)
!local variables
integer :: num_args
character(200), intent(OUT)  :: infilename
character(300)  :: text


! read name of input file from the command line
num_args = command_argument_count()
if (num_args .lt. 1) then
  check_inputfile = .false.
  return
endif

call getarg(num_args, infilename)

text = 'Reading input from '//infilename
call centered_heading(trim(text), '-')


check_inputfile = .true.

end function check_inputfile




END PROGRAM qprep5
