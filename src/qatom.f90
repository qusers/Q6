! (C) 2014 Uppsala Molekylmekaniska HB, Uppsala, Sweden
! qatom.f90
! by John Marelius, Johan Åqvist & Martin Almlöf
! Q-atom force field data and FEP file reading

module QATOM
use SIZES
use NRGY
use MISC
use PRMFILE
use INDEXER
!use MPIGLOB
use TOPO

implicit none

!constants
	character(*), private, parameter	:: MODULE_NAME    = 'Q-atom'
	character(*), private, parameter	:: MODULE_VERSION = '5.06'
	character(*), private, parameter	:: MODULE_DATE    = '2014-01-01'

!	Constants
	real(8), private					:: pi, deg2rad	!set in sub startup

!-----------------------------------------------------------------------
!	fep/evb information
!-----------------------------------------------------------------------

	integer, parameter			::	max_states		= 10
	integer, parameter			::	max_qat			= 100
	integer, parameter			::	max_link		= 10

	integer						::	nstates, nqat

	!******PWadded this variable
	!Topology atom number of the switching atom of the Q-atoms
	!Needed if periodic boundaries are used
	integer						::  qswitch

	integer 	::	offset !offset number for topology atom numbers
	logical						::	qvdw_flag
	logical						::	qq_use_library_charges
	integer(AI), allocatable	::	iqseq(:)
	integer, allocatable			::	qiac(:,:)
	character(len=KEYLENGTH), allocatable	::qtac(:)
	integer																::	nqexpnb
	integer(AI), allocatable							::	iqexpnb(:), jqexpnb(:)

	integer									::	nqlib
	real(4), allocatable		::	qcrg(:,:)
	real(8), allocatable		::	qmass(:)
	real(8), allocatable		::	qavdw(:,:), qbvdw(:,:)

	integer				::	nqbond

	type QBOND_TYPE
		integer(AI)		::	i,j
		integer(TINY)		::	cod(max_states)
	end type QBOND_TYPE

	type QBONDLIB_TYPE
		!Morse diss. en, Morse Alpha, Morse/Harm. r0
		real(8)					::	Dmz, amz, r0
		!Harm. force const
		real(8)					::	fk
	end type QBONDLIB_TYPE

	type(QBOND_TYPE), allocatable	 :: qbnd(:)
	type(QBONDLIB_TYPE), allocatable :: qbondlib(:)


	type QANGLE_TYPE
		integer(AI)				::	i,j,k
		integer(TINY)			::	cod(max_states)
	end type QANGLE_TYPE

	integer						::	nqangle

	type(QANGLE_TYPE), allocatable::	qang(:)
	type(ANGLIB_TYPE), allocatable::	qanglib(:)

	integer												::	nqtor
	integer(AI), allocatable			::	iqtor(:),jqtor(:),kqtor(:),lqtor(:)
	integer(TINY), allocatable		::	qtorcod(:,:)
	real(8), allocatable					::	qfktor(:),qrmult(:),qdeltor(:)

	integer						::	nqimp
	integer(AI), allocatable	::	iqimp(:),jqimp(:),kqimp(:),lqimp(:)
	integer(TINY), allocatable	::	qimpcod(:,:)
	real(8), allocatable		::	qfkimp(:),qimp0(:)

	integer						::	nang_coupl,ntor_coupl,nimp_coupl
	integer(AI)					::	iang_coupl(3,max_qat)
	integer(AI)					::	itor_coupl(3,max_qat)
	integer(AI)					::	iimp_coupl(3,max_qat)

	integer						::	nqshake
	integer(AI)				::	iqshake(max_qat),jqshake(max_qat)
	real(8)						::	qshake_dist(max_qat,max_states)

	integer						::	noffd
	type(OFFDIAG_SAVE), allocatable::	offd(:)
	type(OFFDIAG_AUX), allocatable::	offd2(:)

    integer						::	nexspec
	type SPECEX_TYPE
		integer(AI)				::	i,j
		logical					::	flag(max_states)
	end type SPECEX_TYPE
	type(SPECEX_TYPE), allocatable :: exspec(:)

	! Monitoring of nonbonded interactions between selected groups of atoms
	type monitor_group_pair_TYPE
		integer					::	i,j
		real(8)					::	Vel(max_states), Vlj(max_states)
		real(8)					::	Vwel, Vwlj, Vwsum
	end type  monitor_group_pair_TYPE

	type(monitor_group_pair_TYPE), allocatable::monitor_group_pair(:)

	type monitor_atom_group_TYPE
	  integer, pointer			::	atom(:)	! the atoms
	  integer								::	n		! #atoms
	end type monitor_atom_group_TYPE

	type (monitor_atom_group_TYPE), allocatable::monitor_atom_group(:)

	!the maximal number of atoms in a group to be monitored
	!this is limited by the line length used in the prmfile module!

	integer, parameter			::	MAX_ATOMS_IN_SPECIAL_GROUP=30

	integer						::	monitor_group_pairs, monitor_groups

	!type holding all scaling factors for electostatic interactions in qq-pairs
	type qq_el_scale_TYPE
		integer(AI)      ::iqat,jqat
		real(8)          ::el_scale(max_states)! to hold the el_scale for different states "masoud Oct_2013"
	end type qq_el_scale_TYPE

	type(qq_el_scale_TYPE),allocatable             :: qq_el_scale(:)

	integer               ::nel_scale !number of defined scale factors

    	integer						::	tmpindex,numsoftlines,i2
    	real(8), allocatable				::	sc_lookup(:,:,:), alpha_max(:,:)
		real(8)						::	sc_aq,sc_bq,sc_aj,sc_bj,alpha_max_tmp
		logical						::  softcore_use_max_potential

!-----------------------------------------------------------------------
!	fep/evb energies
!-----------------------------------------------------------------------
	type(Q_ENERGIES), allocatable::	EQ(:)
	real(8)						::	Hij(max_states,max_states)
	real(8)						::	EMorseD(max_qat)
	real(8)						::	dMorse_i(3,max_qat)
	real(8)						::	dMorse_j(3,max_qat)

!miscellany
	logical						::	use_new_fep_format

contains

!-------------------------------------------------------------------------

subroutine qatom_startup
	! initialise used modules
	call prmfile_startup
	call nrgy_startup

	! initialise constants
	pi = 4.0*atan(1.0)
	deg2rad = pi/180.0

100	format(a,' module',t30,'version ',a,t50,'(modified on ',a,')')
end subroutine qatom_startup

!-------------------------------------------------------------------------

subroutine qatom_shutdown
	integer						::	alloc_status
	deallocate(EQ, stat=alloc_status)
	deallocate(iqseq, qiac, iqexpnb, jqexpnb, qcrg, stat=alloc_status)
	deallocate(qmass, stat=alloc_status)
	deallocate(qavdw, qbvdw, stat=alloc_status)
	deallocate(qbnd, qbondlib, stat=alloc_status)
	deallocate(qang, stat=alloc_status)
	deallocate(qanglib, stat=alloc_status)
	deallocate(iqtor, jqtor, kqtor, lqtor, qtorcod, stat=alloc_status)
	deallocate(qfktor, qrmult, qdeltor, stat=alloc_status)
	deallocate(iqimp, jqimp, kqimp, lqimp, qimpcod, stat=alloc_status)
	deallocate(qfkimp,qimp0, stat=alloc_status)
	deallocate(exspec, stat=alloc_status)
	deallocate(offd, offd2, stat=alloc_status)
	if (allocated(qq_el_scale)) deallocate(qq_el_scale)
end subroutine qatom_shutdown

!-------------------------------------------------------------------------

logical function qatom_old_load_atoms(fep_file)
!arguments
	character(*), intent(in)	::	fep_file
  ! *** local variables
  integer					::	i,j,k,iat
  !.......................................................................
	qatom_old_load_atoms = .false.

  call centered_heading('Reading Q atom list','-')
  ! --> fep file (4)

  open (unit=4,file=fep_file,status='old',form='formatted', &
       action='read')

  ! --- # states, # q-atoms

  read (4,*) nstates,nqat
  write (*,20) nstates,nqat
20 format ('No. of fep/evb states    = ',i5,5x, &
       'No. of fep/evb atoms     = ',i5)
  !allocate memory for qatom list
  allocate(iqseq(nqat))


  read (4,*) (iqseq(i),i=1,nqat)
  write (*,40) (iqseq(i),i=1,nqat)
40 format ('Atom nos.:',10i6)
	qatom_old_load_atoms = .true.
end function qatom_old_load_atoms

!-------------------------------------------------------------------------

logical function qatom_load_atoms(fep_file)
	!arguments
	character(*), intent(in)	::	fep_file
	! *** local variables
	integer					::	i,j,k,iat,ires,iatq
	integer					::	s, topno, icase, stat, last, add_res
	integer					::	qflag(max_states)
	integer					::	type_count !counts number of parameters
	logical					::	yes
	character(len=4)		::	offset_name
	character(len=4)		::	res_name
	character(len=7)		::	res_str
	character(len=50)		::	line
	character(len=50)		::	word
	character(len=4), allocatable :: names(:)
	integer					::	offset_residue, max_res, resno
	integer, allocatable	::  residues(:)
  !.......................................................................

	use_new_fep_format = .true.
	qatom_load_atoms = .true.				!assume this for a start
	softcore_use_max_potential = .false.	!default

	call centered_heading('Reading Q atom list','-')

	if(.not. prm_open_section('atoms', fep_file)) then !it's not a new file
		write(*,'(a)') &
			'>>> WARNING: No [atoms] section in fep file. Trying old format.'
		call prm_close
		use_new_fep_format = .false.
		qatom_load_atoms = qatom_old_load_atoms(fep_file)
		return
	end if

	if(.not. prm_open_section('FEP')) then
		nstates = 1
		write(*,21, advance='no') nstates
	else
		if(.not.prm_get_integer_by_key('states', nstates)) then
			nstates = 1
			write(*,21) nstates
		else
			write(*,20) nstates
		end if
		yes = prm_get_logical_by_key('qq_use_library_charges', qq_use_library_charges, .false.)
		yes = prm_get_logical_by_key('softcore_use_max_potential', softcore_use_max_potential, .false.)

		offset = -1
		!should an offset be applied to topology atom numbers?
		if(prm_get_integer_by_key('offset', offset)) then
			!got an atom number offset
			if(offset < 1 .or. offset > nat_solute) then
				!it's invalid
				write(*, '(a,i5)') '>>>>> ERROR: Invalid topology atom number offset value:', offset
				qatom_load_atoms = .false.
				return
			end if
			ligand_offset = offset
		elseif(prm_get_string_by_key('offset_name', offset_name)) then
			do i = 1, nres_solute
				if(offset_name == res(i)%name) then
					offset = res(i)%start - 1
					exit
				end if
			end do
			if(offset == -1) then !not found
				write(*, '(a,a4,a)') '>>>>> ERROR: Residue name ',offset_name, &
					'not found.'
				qatom_load_atoms = .false.
				return
			end if
			ligand_offset = offset
		elseif(prm_get_integer_by_key('offset_residue', offset_residue)) then
			if(offset_residue < 1 .or. offset_residue > nres_solute) then
				write(*, '(a,i5)') '>>>>> ERROR: Invalid residue number for offset:',offset_residue
				qatom_load_atoms = .false.
				return
			end if
			offset = res(offset_residue)%start - 1
		else
			offset = 0
		end if
	end if

	if(nstates == 0) then
		write(*,'(/,a)') &
			'>>>>> ERROR: Number of states must be at least 1. Aborting.'
		qatom_load_atoms = .false.
		return
	end if

	yes = prm_open_section('atoms') !by now we know it's there
	type_count = prm_count('atoms')	!count number of q-atom lines
	if(type_count == 0) then
		write(*,'(a)') &
			'>>> WARNING: Number of Q-atoms is zero. Fep file will not be loaded.'
		call prm_close
		!qatom_load_atoms = .false. !this condition IS OK.
		return
	end if

	if (.not. prm_get_line(line)) then
		write(*, '(a)') ">>>>> ERROR: reading line in 'atoms' section"
		qatom_load_atoms = .false.
		return
	end if
	icase=0
	read(line, fmt=*, iostat=stat) word, i !read integer from value
	if (stat == 0) then
		read(word, fmt=*, iostat=stat) i !read integer from value
		if (stat == 0) then
			icase = 1   !int_int
		else
			icase = 2 !string_int
		end if
	else
		icase = 3  !string_string
	end if

	select case(icase)
	case (1)
		nqat = prm_max_enum('atoms', type_count) !count number of q-atoms & get highest q-atom number
		!allocate memory for qatom list
		allocate(iqseq(nqat))
		yes = prm_open_section('atoms') !rewind section
		do i = 1, type_count
			if(prm_get_int_int(s, topno)) then
				if(topno + offset < 1 .or. topno + offset > nat_solute) then
					write(*, '(a,i5,a, i2)') '>>>>> ERROR: invalid topology atom number', &
						topno + offset, ' for Q-atom',i
					qatom_load_atoms = .false.
				end if
				iqseq(s) = topno + offset
			else
				write(*,'(a,i2)') &
					'>>> WARNING: Failed to read Q-atom ',i
				qatom_load_atoms = .false.
			end if
		end do
	!Input is of type 'res xx'
	case(2)
		yes = prm_open_section('atoms') !rewind section
		if (offset > 0) then
			write(*,'(a)') &
				'>>> WARNING: Offset can only be used when defining q-atoms with atom nubers. \n Offset will be set to zero.'
				offset=0
		end if

		max_res = prm_max_enum2('atoms', type_count) !count number of residues & get highest residue number
		if (max_res > nres_solute) then
			write(*, '(a,i5)') '>>>>> ERROR: invalid topology solute residue number', max_res
			qatom_load_atoms = .false.
			return
		end if
		allocate(residues(type_count))
		nqat=0
		!Count number of q-atoms
		do i = 1, type_count
			if(prm_get_string_int(res_str, resno)) then
				call upcase(res_str)
				if(res_str == 'RES' .or. res_str == 'RESIDUE') then
					residues(i)=resno
					if (resno < nres_solute) then
						nqat = nqat + (res(resno+1)%start - res(resno)%start)
					else
						nqat = nqat + ((nat_solute +1) - res(resno)%start)
					end if
				else
					write(*, '(a)') ">>>>> ERROR: invalid selection syntax in fep file, section 'atoms'"
					qatom_load_atoms = .false.
					return
				end if
			else
				write(*,'(a,i2)') &
					'>>> WARNING: Failed to read Q-residue ',i
				qatom_load_atoms = .false.
			end if
		end do
		allocate(iqseq(nqat))
		iatq=1
		!Now assign iqseq with proper topology atom numbers
		do i = 1, type_count
			ires=residues(i)
			if (ires < nres_solute) then
				last=res(ires+1)%start-1
				do iat = res(ires)%start,last
					iqseq(iatq)= iat
					iatq=iatq+1
				end do
			else
				do iat = res(ires)%start, nat_solute
					iqseq(iatq)= iat
					iatq=iatq+1
				end do
			end if
		end do
		offset = 0; !offset, offset not compatible with this definition
		write (*,31) (residues(i),i=1,type_count)

	!it is string_string, 'all TYPE'
    case(3)
		yes = prm_open_section('atoms') !rewind
		if (offset > 0) then
			write(*,'(a)') &
				'>>> WARNING: Offset can only be used when defining q-atoms with atom nubers. \n Offset will be set to zero.'
				offset=0
		end if

		allocate(names(type_count))
		do i=1, type_count
			if(.not. prm_get_string_string(res_str, res_name)) then
				write(*, '(a)') ">>>>> ERROR: invalid syntax in fep file, section 'atoms'"
				qatom_load_atoms = .false.
				return
			end if
			call upcase(res_str)
			if(.not. (res_str == 'ALL')) then
				write(*, '(a)', advance='no') ">>>>> ERROR: invalid syntax in fep file, &
					&section 'atoms': ",res_str
				qatom_load_atoms = .false.
				return
			end if
			call upcase(res_name)
			names(i)=res_name
		end do
		add_res=0
		allocate(residues(nres_solute))  !allocate for worst case
		do ires=1,nres_solute
			do j=1,type_count
				if(names(j) == res(ires)%name) then
					add_res=add_res+1
					residues(add_res)=ires
					exit
				end if
			end do
		end do
		if (add_res == 0) then
			write(*, '(a)') ">>>>> ERROR: No matching residues found in topology. Could not assign qatoms."
			qatom_load_atoms = .false.
			return
		end if
		!Count number of q-atoms
		do i = 1, add_res
			ires=residues(i)
			if (ires < nres_solute) then
				nqat = nqat + (res(ires+1)%start - res(ires)%start)
			else
				nqat = nqat + (nat_solute +1 - res(ires)%start)
			end if
		end do
		allocate(iqseq(nqat))
		iatq=1
		!Now assign iqseq with proper topology atom numbers
		do i = 1, add_res
			ires=residues(i)
			if (ires < nres_solute) then
				last=res(ires+1)%start-1
				do iat = res(ires)%start,last
					iqseq(iatq)= iat
					iatq=iatq+1
				end do
			else
				do iat = res(ires)%start, nat_solute
					iqseq(iatq)= iat
					iatq=iatq+1
				end do
			end if
		end do
		offset = 0; !offset, offset not compatible with this definition
		write (*,31) (residues(i),i=1,add_res)

	case default
		write(*,'(a,i2)') &
			">>> ERROR: Failed to read fep file. Syntax error in 'atoms' section."
		qatom_load_atoms = .false.
		return
	end select

  if (allocated(residues)) deallocate(residues)
  if (allocated(names)) deallocate(names)

  write (*,30) nqat
  write (*,40) (iqseq(i),i=1,nqat)

20 format ('No. of fep/evb states    = ',i5)
21 format ('Default fep/evb states   = ',i5)
25	format('Offset for topology atom numbers = ',i5)
30 format ('No. of fep/evb atoms     = ',i5)
31 format ('Assigning q-atoms from residues: ',5i6)
40 format ('Atom nos.:',10i6)

end function qatom_load_atoms


!-------------------------------------------------------------------------

logical function qatom_old_load_fep()
  ! *** local variables
  character					::	libtext*80,qaname*2
  integer					::	i,j,k,iat
  integer					::	nqcrg,nqcod
	integer					::	qflag(max_states)
	!temp. array to read integer flags before switching to logicals
	integer					::	exspectemp(max_states)

	qatom_old_load_fep  = .false.

  call centered_heading('Reading fep/evb strategy','-')

  !allocate memory for qatom arrays
	allocate(qiac(nqat,nstates))
	allocate(iqexpnb(nqat))
	allocate(jqexpnb(nqat))
	!qcrg may be allocated in MD to copy topology charges
	if(.not. allocated(qcrg)) allocate(qcrg(nqat,nstates))

  ! --- Set new charges
  read (4,*) nqcrg
	if(nqcrg > 0) then
		write (*,60) nqcrg

60		format (/,'No. of changing charges  = ',i5)
80		format ('q-atom state_1 state_2 ...')
		do i=1,nqcrg
			read (4,*) iat,(qcrg(iat,j),j=1,nstates)
	        write (*,'(i6,8f8.4)') iat,(qcrg(iat,j),j=1,nstates)
		end do
	end if


  ! --- Set new vdw params / only if qvdw_flag is true

  read (4,*) iat
  if(iat > 0) qvdw_flag = .true.
  if(qvdw_flag) then
	write (*,100)
100 format (/,'Q-atom vdW parameters are to be redefined:')


     write (*,80)
     do i=1,nqat
        read (4,*) iat,(qiac(iat,j),j=1,nstates)
        write (*,'(i6,8i8)') iat,(qiac(iat,j),j=1,nstates)
     end do

  ! Read Qatom type library
     read (4,*) nqlib
     write (*,120) nqlib
120  format (/,'No. fep/evb lib entries  = ',i5)
     read (4,'(a80)') libtext
     write (*,'(a80)') libtext

	 allocate(qmass(nqlib), qavdw(nqlib,nljtyp), qbvdw(nqlib,nljtyp))

     do i=1,nqlib
        read (4,*) j,qaname,(qavdw(j,k),qbvdw(j,k),k=1,3),qmass(j)
        write (*,140)j,qaname,(qavdw(j,k),qbvdw(j,k),k=1,3),qmass(j)
     end do
140  format (i5,3x,a2,1x,f9.2,2f8.2,f6.2,f9.2,2f8.2)
     !

     read (4,*) nqexpnb
     write (*,144) nqexpnb
144  format (/,'No. C*exp(-ar) nb pairs  = ',i5)
     do i=1,nqexpnb
        read (4,*) iqexpnb(i),jqexpnb(i)
        write (*,146) iqexpnb(i),jqexpnb(i)
     end do
146  format ('atom_i, atom_j             : ',2i5)

  end if

  ! --- Set new bonds

  read (4,*) nqbond
  if(nqbond > 0) write (*,160) nqbond
160 format (/,'No. of changing bonds    = ',i5)
  allocate(qbnd(nqbond))
  do i=1,nqbond
     read (4,*) qbnd(i)%i,qbnd(i)%j
     read (4,*) (qflag(j),j=1,nstates)
     read (4,*) (qbnd(i)%cod(j),j=1,nstates)
     write (*,180) qbnd(i)%i,qbnd(i)%j
     write (*,200) (qflag(j),j=1,nstates)
     write (*,220) (qbnd(i)%cod(j),j=1,nstates)
		!set bond type = 0 where bond presence flag = 0
		do j=1,nstates
			if(qflag(j) == 0) qbnd(i)%cod(j) = 0
		end do
  end do
180 format (/,'atom_i -- atom_j              : ',2i5)
200 format ('exists in state_1 state_2 ... : ',10i5)
220 format ('codes  in state_1 state_2 ... : ',10i5)


  read (4,*) nqcod
  allocate(qbondlib(nqcod))
	qbondlib(:)%fk = 0. !clear Harmonic f.c, not available with old FEP file
  if ( nqcod .gt. 0 ) then
     write (*,*)
     write (*,'(a)') 'Morse  E_diss   alpha      b0'
     do i=1,nqcod
        read (4,*) qbondlib(i)%Dmz,qbondlib(i)%amz,qbondlib(i)%r0
        write (*,225) i,qbondlib(i)%Dmz,qbondlib(i)%amz,qbondlib(i)%r0
     end do
225  format (i5,3f8.2)
     write (*,*)
  end if

  ! --- Set new angles

  read (4,*) nqangle
  if(nqangle > 0) write (*,260) nqangle
260 format (/,'No. of changing angles   = ',i5)
  allocate(qang(nqangle))

  do i=1,nqangle
     read (4,*) qang(i)%i,qang(i)%j,qang(i)%k
     read (4,*) qflag(1:nstates)
     read (4,*) qang(i)%cod(1:nstates)
     write (*,280) qang(i)%i,qang(i)%j,qang(i)%k
     write (*,200) qflag(1:nstates)
     write (*,220) qang(i)%cod(1:nstates)
		!set code to 0 where presence flag = 0
		do j=1,nstates
			if(qflag(j) == 0) qang(i)%cod(j) = 0
		end do
  end do
280 format (/,1x,'atom_i -- atom_j -- atom_k    : ',3i5)

  read (4,*) nqcod
  allocate(qanglib(nqcod))
	!set new features to 0.
	qanglib(:)%ureyfk = 0.
	qanglib(:)%ureyr0 = 0.
  if ( nqcod .gt. 0 ) then
     write (*,*)
     write (*,'(a)') 'Angle force-k  theta0'
     do i=1,nqcod
        read (4,*) qanglib(i)%fk,qanglib(i)%ang0
        write (*,225) i,qanglib(i)%fk,qanglib(i)%ang0
        qanglib(i)%ang0 = deg2rad*qanglib(i)%ang0
     end do
     write (*,*)
  end if

  ! --- Set new torsions

  read (4,*) nqtor
  if(nqtor > 0) write (*,360) nqtor
360 format (/,'No. of changing torsions = ',i5)
  allocate(iqtor(nqtor), &
	jqtor(nqtor), &
	kqtor(nqtor), &
	lqtor(nqtor), &
	qtorcod(nqtor,nstates))

  do i=1,nqtor
     read (4,*) iqtor(i),jqtor(i),kqtor(i),lqtor(i)
     read (4,*) (qflag(j),j=1,nstates)
     read (4,*) (qtorcod(i,j),j=1,nstates)
     write (*,380) iqtor(i),jqtor(i),kqtor(i),lqtor(i)
     write (*,200) (qflag(j),j=1,nstates)
     write (*,220) (qtorcod(i,j),j=1,nstates)
		do j=1, nstates
			if(qflag(j) == 0) qtorcod(i,j) =0
		end do
  end do
380 format (/,'at_i -- at_j -- at_k -- at_l  : ',4i5)

  read (4,*) nqcod
  allocate(qfktor(nqcod), &
	qrmult(nqcod), &
	qdeltor(nqcod))
  if ( nqcod .gt. 0 ) then
     write (*,*)
     write (*,'(a)') ' Tors force-k    mult   delta'
     do i=1,nqcod
        read (4,*) qfktor(i),qrmult(i),qdeltor(i)
        write (*,225) i,qfktor(i),qrmult(i),qdeltor(i)
        qdeltor(i) = deg2rad*qdeltor(i)
     end do
     write (*,*)
  end if

  ! --- Set new impropers

  read (4,*) nqimp
  if(nqimp > 0) write (*,460) nqimp
460 format (/,'No. of changing impropers= ',i5)

  allocate(iqimp(nqimp), &
	jqimp(nqimp), &
	kqimp(nqimp), &
	lqimp(nqimp), &
	qimpcod(nqimp,nstates))

  do i=1,nqimp
     read (4,*) iqimp(i),jqimp(i),kqimp(i),lqimp(i)
     read (4,*) (qflag(j),j=1,nstates)
     read (4,*) (qimpcod(i,j),j=1,nstates)
     write (*,380) iqimp(i),jqimp(i),kqimp(i),lqimp(i)
     write (*,200) (qflag(j),j=1,nstates)
     write (*,220) (qimpcod(i,j),j=1,nstates)
		!set type = 0 if presence flag = 0
		do j=1,nstates
			if(qflag(j) == 0) qimpcod(i,j) = 0
		end do
  end do

  read (4,*) nqcod
  allocate(qfkimp(nqcod), qimp0(nqcod))

  if ( nqcod .gt. 0 ) then
     write (*,*)
     write (*,'(a)') ' Impr force-k    imp0'
     do i=1,nqcod
        read (4,*) qfkimp(i),qimp0(i)
        write (*,225) i,qfkimp(i),qimp0(i)
        qimp0(i) = deg2rad*qimp0(i)
     end do
     write (*,*)
  end if


  ! --- Read angle,torsion and improper couplings to Morse bonds

  read (4,*) nang_coupl
  if(nang_coupl > 0) write (*,544) nang_coupl
544 format (/,'No. of angle-Morse couplings = ',i5)
  do i=1,nang_coupl
     read (4,*) (iang_coupl(j,i),j=1,2)
     write (*,545) (iang_coupl(j,i),j=1,2)
  end do
545 format ('angle_i, bond_j              : ',2i5)

  read (4,*) ntor_coupl
  if( ntor_coupl > 0) write (*,548) ntor_coupl
548 format ('No. of tors.-Morse couplings = ',i5)
  do i=1,ntor_coupl
     read (4,*) (itor_coupl(j,i),j=1,2)
     write (*,549) (itor_coupl(j,i),j=1,2)
  end do
549 format ('torsion_i, bond_j            : ',2i5)

  read (4,*) nimp_coupl
  if(nimp_coupl > 0 ) write (*,552) nimp_coupl
552 format ('No. of impr.-Morse couplings = ',i5)
  do i=1,nimp_coupl
     read (4,*) (iimp_coupl(j,i),j=1,2)
     write (*,553) (iimp_coupl(j,i),j=1,2)
  end do
553 format ('improper_i, bond_j           : ',2i5)

  ! --- Read Q-atom shake constraints

  read (4,*) nqshake
  if(nqshake > 0) write (*,560) nqshake
560 format ('No. fep/evb shake contraints = ',i5)
  do i=1,nqshake
     read (4,*) iqshake(i),jqshake(i),(qshake_dist(i,j),j=1,nstates)
     write (*,580) iqshake(i),jqshake(i),(qshake_dist(i,j),j=1,nstates)
  end do
580 format ('i -- j, dist in state_1, state_2, ... : ',2i5,10f6.2)

  ! --- Read off-diagonal hamiltonial matrix functions

	read (4,*) noffd
	if(noffd > 0) write (*,600) noffd
600 format ('No. offdiagonal (Hij) funcs. = ',i5)
	if(noffd > max_states) then
		write(*,'(a,i2,a)') &
			'>>> Error: Too many off-diagonal functions, (max_states is ', &
				max_states, ' Aborting.'
		return
	end if
	allocate(offd(noffd), offd2(noffd))
	if ( noffd .gt. 0 ) then
		write (*,'(a)') 'state_i state_j atom_k atom_l     Aij mu_ij'
	end if
	do i=1,noffd
		read (4,*) offd(i)%i,offd(i)%j,offd2(i)%k,offd2(i)%l,offd2(i)%A, &
          offd2(i)%mu
		write (*,620) offd(i)%i,offd(i)%j,offd2(i)%k,offd2(i)%l,offd2(i)%A, &
          offd2(i)%mu
	end do
620 format (i7,1x,3i7,f8.2,f6.2)

! --- Read special exclusions among quantum atoms

      read (unit=4,fmt=*,end=999,err=999) nexspec
      if(nexspec > 0) write (*,585) nexspec
	  allocate(exspec(nexspec))

 585  format (/,'No. special exclusions       = ',i5)
      do i=1,nexspec
         read (4,*) exspec(i)%i,exspec(i)%j
         read (4,*) exspectemp(1:nstates)
         write (*,590) exspec(i)%i,exspec(i)%j
         write (*,200) exspectemp(1:nstates)
			do j = 1, nstates
				if(exspectemp(j) == 0) then
					exspec(i)%flag(j) = .false.
				else if(exspectemp(j) == 1) then
					exspec(i)%flag(j) = .true.
				else
					write(*,592)
					return
				end if
			end do
       end do
590  format ('i -- j  special exclusion pair  , ... : ',2i5)
592	format('>>>>> ERROR: Special exclusion state flags are invalid.')


999  close(4)
	qatom_old_load_fep = .true.
  !.......................................................................
end function qatom_old_load_fep

!-----------------------------------------------------------------------
logical function qatom_load_fep(fep_file)
!arguments
	character(*), intent(in)	::	fep_file
  ! *** local variables
  character(len=200)		::	line
  integer					::	i,j,k,l,iat,st
  integer					::	nqcrg,nqcod
  character(len=40)			::	section
  logical, allocatable, dimension(:)		::	type_read
  integer					::	type_count, filestat
  real(8)                   ::  el_scale(max_states) !local variable for scaling of different states "masoud Oct_2013"
  integer                   ::  stat

	!temp. array to read integer flags before switching to logicals
	integer					::	exspectemp(max_states)
	character(len=keylength)	:: qtac_tmp(max_states)

	!temp array for reading special atom group members
	integer                   ::  temp_atom(MAX_ATOMS_IN_SPECIAL_GROUP)


	if(.not. use_new_fep_format) then
		qatom_load_fep = qatom_old_load_fep()
		return
	end if

	qatom_load_fep = .true.

	if(.not. prm_open(fep_file)) then
		write(*,'(a,a)') '>>>>> ERROR: Could not open fep file ',trim(fep_file)
		qatom_load_fep = .false.
		return
	end if

	call centered_heading('Reading fep/evb strategy','-')

	!When PBC is used, a switching atom for the Q-atoms has to be defined.
	if( use_PBC ) then
		if( .not. prm_open_section('PBC') ) then
			write(*,51)
			qatom_load_fep = .false.
			return
		else
			if( .not. prm_get_integer_by_key('switching_atom', qswitch) ) then
				write(*,50)
				qatom_load_fep = .false.
			else
				qswitch = qswitch + offset
				write (*,'(a,i6,a)') 'Using topology atom number ',qswitch,' when creating q-atom based nonbond lists.'
			end if
		end if
	end if
51		format('>>>>> ERROR: Section PBC is required when using periodic boundary.')
50		format('>>>>> ERROR: Switching atom could not be read.')



  !allocate memory for qatom arrays
	allocate(qiac(nqat,nstates))
	allocate(iqexpnb(nqat))
	allocate(jqexpnb(nqat))
	!qcrg may be allocated in MD to copy topology charges
	if(.not. allocated(qcrg)) allocate(qcrg(nqat,nstates))

	!read flag for use of library charges in qq-nonbond
	!assume section FEP is open
	! --- Set new charges
	section = 'change_charges'
	nqcrg = prm_count(section)
	if(nqcrg > 0) then
		write (*,60) nqcrg
		if(qq_use_library_charges) then
			write(*,'(a)') 'Intra-Q-atom electrostatic interactions will not be changed.'
		end if

60		format (/,'No. of changing charges  = ',i5)
		do i=1,nqcrg
			if(.not. prm_get_line(line)) goto 1000
			!read atom to check
			read(line,*,err=1000) iat
			if(iat < 1 .or. iat > nqat) then
				write(*,82) iat
				qatom_load_fep = .false.
				cycle !dont even try to read charges
			else if(iqseq(iat) == 0) then
				write(*,82) iat
				qatom_load_fep = .false.
				cycle !dont even try to read charges
			end if
			!read atom and all charges

			read (line,*, err=1000) iat,qcrg(iat,:)
		end do

		!list Q-atom charges
		write(*,29)
		write(*,30) ('state',i,i=1,nstates)
		do i=1,nqat
			write(*,31) i, qcrg(i,:)
		end do
		write(*,32) sum(qcrg(:,:), dim=1)

	end if
29	format('Effective Q-atom charges for all Q-atoms')
30	format('Q atom    charge in',7(1x,a5,i2))
31	format(     i6,t20,10f8.3)
32	format(/,'   SUM',t20,10f8.3)

82	format('>>>>> ERROR: ',i2,' is not a valid q-atom number')

	! Read Qatom type library only if iqvdw_flag is true
	section = 'atom_types'
	nqlib = prm_count(section)
	if(nqlib > 0) then
		write (*,120) nqlib
		if(ivdw_rule==VDW_GEOMETRIC) then
			write(*,130)
		else
			write(*,131)
		end if
		allocate(qtac(nqlib))
		call index_create(nqlib) !set up atom type name lookup table
		!clear read flag for all parameters
120		format (/,'No. of Q-atom types  = ',i5)
130		format('Name            Ai      Bi      Ci     ai Ai(1-4) Bi(1-4)    Mass')
131		format('Name            R*i     ei      Ci     ai R*i(1-4)ei(1-4)    Mass')

		allocate(qmass(nqlib), qavdw(nqlib,nljtyp), qbvdw(nqlib,nljtyp))

		do i=1,nqlib
			if(.not. prm_get_line(line)) goto 1000
			read (line,*, err=1000) qtac(i),(qavdw(i,k),qbvdw(i,k),k=1,nljtyp),qmass(i)
			write (*,140) qtac(i),(qavdw(i,k),qbvdw(i,k),k=1,nljtyp),qmass(i)
			if(.not. index_add(qtac(i), i)) then
				write(*,83) qtac(i)
				qatom_load_fep = .false.
			end if
		end do
140		format (a8,1x,f9.2,2f8.2,f6.2,f9.2,2f8.2)
83		format ('>>>>> ERROR: Could not enumerate q-atom type ',a,' Duplicate name?')
	end if

	! --- Set new vdw params
	section = 'change_atoms'
	iat = prm_max_enum(section, type_count)
	if(iat == 0) then
		qvdw_flag = .false.
	else if(iat /= nqat .or. type_count /= nqat) then
		write(*,'(a)') '>>>>> ERROR: Atom types of Q-atoms must be given for every Q-atom!'
		qatom_load_fep = .false.
		return
	else
		qvdw_flag = .true.
		write (*,100)
100		format (/,'Assigning Q-atom types to all Q atoms:')
80		format('Q atom    atom type in',6(2x,a5,i2))
		write (*,80) ('state',i,i=1,nstates)
		do i=1,nqat !we know nqat = number of type_count
			!reading should not fail at this stage - get_line was called by prm_count!
			if(.not. prm_get_line(line)) goto 1000
        	read (line,*, err=1000) iat,qtac_tmp(1:nstates)
			write (*,'(i6,t25,6(a8,1x))') iat,qtac_tmp(1:nstates)
			do j = 1, nstates
				!check that atom type is defined and assign qiac(iat,j)
				if(.not. index_get(qtac_tmp(j), qiac(iat,j))) then
					write(*,113) qtac_tmp(j)
					qatom_load_fep = .false.
				end if
			end do !checking
		end do
	end if !if redefine atom types
112	format('>>>>> ERROR: ',a,' type ',i2,' has not been defined.')
113	format('>>>>> ERROR: Q-atom type ',a,' has not been defined.')

! --- Soft repulsion pairs
	section = 'soft_pairs'
	nqexpnb = prm_count(section)
	if(nqexpnb > 0) then
		write (*,144) nqexpnb
		write(*,145)
		do i=1,nqexpnb
			if(.not. prm_get_line(line)) goto 1000
			read(line,*, err=1000) j, k
			if(j < 1 .or. j > nqat .or. k < 1 .or. k > nqat) then
				write(*,148) j,k
				qatom_load_fep = .false.
				cycle !dont even try to read charges
			else if(iqseq(j) == 0 .or. iqseq(k) == 0) then
				write(*,148) j,k
				qatom_load_fep = .false.
				cycle !dont even try to read
			end if
			iqexpnb(i) = j
			jqexpnb(i) = k
			write (*,146) iqexpnb(i),jqexpnb(i)
		end do
	end if
144	format (/,'No. of soft repulsion non-bonded pairs = ',i5)
145	format('q-atom_i q-atom_j')
146	format (i8,1x,i8)
148	format('>>>>> ERROR: Invalid q-atom number in this group: ',4i5)

! --- Read scaling factors for el. interactions between qq-atoms
	section = 'el_scale'
	nel_scale = prm_count(section)
	if(nel_scale > 0) then
	allocate(qq_el_scale(nel_scale))
		write (*,154) nel_scale
		write(*,155) ('state',i, i=1,nstates) !heading for qq_el_scale at different states "masoud Oct_2013"
		do i=1,nel_scale
			if(.not. prm_get_line(line)) goto 1000
			read(line,*, err=1000) j, k, el_scale(1:nstates)
			if(j < 1 .or. j > nqat .or. k < 1 .or. k > nqat) then !Check if qatoms...
				write(*,148) j,k
				qatom_load_fep = .false.
				cycle
			else if(iqseq(j) == 0 .or. iqseq(k) == 0) then
				write(*,148) j,k
				qatom_load_fep = .false.
				cycle
			end if
            qq_el_scale(i)%iqat=j
            qq_el_scale(i)%jqat=k
            qq_el_scale(i)%el_scale(1:nstates)=el_scale(1:nstates) !assigning scale factor to qq_el_scale variable "masoud Oct_2013"
!			write (*,156) qq_el_scale(i)%iqat,qq_el_scale(i)%jqat,qq_el_scale(i)%el_scale
!			write (*,156) qq_el_scale(i)%iqat,qq_el_scale(i)%jqat,(qq_el_scale(i)%el_scale(h), h=1,nstates)
   			write (*,156) qq_el_scale(i)%iqat,qq_el_scale(i)%jqat,(qq_el_scale(i)%el_scale(1:nstates))   
		end do
	end if
154	format (/,'No. of electrostatic scaling factors between q-q-atoms = ',i5)
155	format('q-atom_i q-atom_j el_scale' , 7(1x,a5,i2))
156	format (i8,1x,i8,7x,7f8.2) !print out up to 7 states "masoud Oct_2013"

! --- Read special exclusions among quantum atoms
	section='excluded_pairs'
	nexspec = prm_count(section)
    if(nexspec > 0) then
		write (*,585) nexspec
		write(*,586) ('state',i, i=1,nstates)
		allocate(exspec(nexspec))

585		format (/,'No. of excluded non-bonded pairs = ',i5)
		do i=1,nexspec
			if(.not. prm_get_line(line)) goto 1000
			read(line,*, err=1000) exspec(i)%i, exspec(i)%j, exspectemp(1:nstates)
			exspec(i)%i = exspec(i)%i + offset
			exspec(i)%j = exspec(i)%j + offset
			write (*,587) exspec(i)%i,exspec(i)%j,exspectemp(1:nstates)
			do j = 1, nstates
				if(exspectemp(j) == 0) then
					exspec(i)%flag(j) = .false.
				elseif(exspectemp(j) == 1) then
					exspec(i)%flag(j) = .true.
				else
					write(*,592)
					qatom_load_fep = .false.
				end if
			end do
		end do
586	format('atom_i atom_j    excluded in', 7(1x,a5,i2))
587	format(i6,1x,i6,t29,7i8)
592		format('>>>>> ERROR: Special exclusion state flags are invalid.')
	end if


	! --- Set new bonds
	section='bond_types'
	nqcod = prm_max_enum(section, type_count)
	if(nqcod > 0 ) then
		allocate(type_read(nqcod))
	type_read(:) = .false.
		allocate(qbondlib(nqcod))
		write (*,150)
		write (*,'(a)') 'type #  Morse E_diss    alpha       b0  Harmonic force_k '

		do i=1,type_count
			if(.not. prm_get_line(line)) goto 1000
			read(line, *, iostat=filestat) j, qbondlib(j)%Dmz,qbondlib(j)%amz, &
				qbondlib(j)%r0
			if(filestat /= 0) then !could not read Morse - try harmonic
				read(line, *, iostat=filestat) j,qbondlib(j)%fk,qbondlib(j)%r0
				if(filestat /= 0) then
					goto 1000
				else
					!Harmonic OK - clear Morse
					qbondlib(j)%Dmz = 0
					qbondlib(j)%amz = 0
					type_read(j) =.true.
					write (*,226) j, qbondlib(j)%r0, qbondlib(j)%fk
				end if
			else
				!Morse OK - clear harmonic
				qbondlib(j)%fk = 0
				type_read(j) =.true.
				write (*,225) j,qbondlib(j)%Dmz,qbondlib(j)%amz, qbondlib(j)%r0
			end if
		end do
150		format (/,'Q-bond types:')
225		format(i6,6x,f8.2,1x,f8.2,1x,f8.2)
226		format(i6,24x,f8.2,10x,f8.2)
		write (*,*)
	end if

	section = 'change_bonds'
	nqbond = prm_count(section)
	if(nqbond > 0)  then
		write (*,160) nqbond
		write(*,161) ('state',i,i=1,nstates)
		allocate(qbnd(nqbond))
		do i=1,nqbond
			if(.not. prm_get_line(line)) goto 1000
			read(line,*, err=1000) qbnd(i)%i, qbnd(i)%j, qbnd(i)%cod(1:nstates)
			qbnd(i)%i = qbnd(i)%i + offset
			qbnd(i)%j = qbnd(i)%j + offset
			write (*,162) qbnd(i)%i,qbnd(i)%j,qbnd(i)%cod(1:nstates)
			!check types
			do j=1,nstates
				if(qbnd(i)%cod(j) >0) then
					if (allocated(type_read)) then
					if(.not. type_read(qbnd(i)%cod(j))) then
						write(*,112) 'Q-bond', qbnd(i)%cod(j)
						qatom_load_fep = .false.
					end if
					else
					    write(*,112) 'Q-bond', qbnd(i)%cod(j)
					    qatom_load_fep = .false.
					end if
				end if
			end do
		end do

160		format (/,'No. of changing bonds    = ',i5)
161		format('atom_i atom_j   bond type in',6(1x,a5,i2))
162		format(i6,1x,i6,t29,6i8)
	end if
	if (allocated(type_read)) then
	    deallocate(type_read)
	end if


	! --- Set new angles
	section = 'angle_types'
	nqcod = prm_max_enum(section, type_count)
	if ( nqcod .gt. 0 ) then
		allocate(type_read(nqcod))
	type_read(:) = .false.
		allocate(qanglib(nqcod))
		!set optional params to 0
		qanglib(:)%ureyfk = 0.
		qanglib(:)%ureyr0 = 0.
		write (*,'(/,a)') 'Q-angle types:'
		write (*,'(a)', advance='no') 'type #       force-k   theta0'
		do i=1,type_count
			if(.not. prm_get_line(line)) goto 1000
			read(line,*, iostat=filestat) j, qanglib(j)
			if(filestat == 0) then
				if(i == 1) write(*,'(a)') '  Urey-Bradley force-k      r0'
				type_read(j) = .true.
				write (*,230) j,qanglib(j)
			elseif(filestat < 0) then
				if(i == 1) write(*,*)
				type_read(j) = .true.
				write (*,225) j,qanglib(j)%fk, qanglib(j)%ang0
			else
				goto 1000
			end if
			qanglib(j)%ang0 = deg2rad*qanglib(j)%ang0
		end do
		write (*,*)
	end if

230	format(i6,2f8.2,14x,2f8.2)

	section = 'change_angles'
	nqangle = prm_count(section)
	if(nqangle > 0) then
		write (*,260) nqangle
260		format (/,'No. of changing angles   = ',i5)
		write(*,261) ('state',i,i=1,nstates)
		allocate(qang(nqangle))

		do i=1,nqangle
			if(.not. prm_get_line(line)) goto 1000
			read(line,*, err=1000) qang(i)%i,qang(i)%j,qang(i)%k, qang(i)%cod(1:nstates)
			qang(i)%i = qang(i)%i + offset
			qang(i)%j = qang(i)%j + offset
			qang(i)%k = qang(i)%k + offset
			write (*,262) qang(i)%i,qang(i)%j,qang(i)%k, qang(i)%cod(1:nstates)
			!check types
			do j=1,nstates
				if(qang(i)%cod(j) > 0) then
					if (allocated(type_read)) then
					if(.not. type_read(qang(i)%cod(j))) then
						write(*,112) 'Q-angle', qang(i)%cod(j)
						qatom_load_fep = .false.
					end if
					else
					    write(*,112) 'Q-angle', qang(i)%cod(j)
					    qatom_load_fep = .false.
					end if
				end if
			end do !j
		end do !i
261		format('atom_i atom_j atom_k    angle type in',5(1x,a5,i2))
262		format(i6,1x,i6,1x,i6,t38,5i8)
	end if
	if (allocated(type_read)) then
	    deallocate(type_read)
	end if

  ! --- Set new torsions
	section = 'torsion_types'
	nqcod = prm_max_enum(section, type_count)
	allocate(qfktor(nqcod), qrmult(nqcod), qdeltor(nqcod))
	if ( nqcod .gt. 0 ) then
		allocate(type_read(nqcod))
		type_read(:) = .false.
		write (*,350)
		do i=1,type_count
			if(.not. prm_get_line(line)) goto 1000
			read(line,*, err=1000) j, qfktor(j),qrmult(j),qdeltor(j)
			type_read(j) = .true.
			write (*,225) j,qfktor(j),qrmult(j),qdeltor(j)
			qdeltor(j) = deg2rad*qdeltor(j)
		end do
		write (*,*)
	end if
350	format(/,'Q-torsion types:',/,'type #       force-k     mult    delta')
	section='change_torsions'
	nqtor = prm_count(section)
	if(nqtor > 0) then
		write (*,360) nqtor
		write(*,361) ('state',i,i=1,nstates)
360		format (/,'No. of changing torsions = ',i5)
		allocate(iqtor(nqtor), &
			jqtor(nqtor), &
			kqtor(nqtor), &
			lqtor(nqtor), &
			qtorcod(nqtor,nstates))

		do i=1,nqtor
			if(.not. prm_get_line(line)) goto 1000
			read(line,*, err=1000) iqtor(i),jqtor(i),kqtor(i),lqtor(i),qtorcod(i,:)
			iqtor(i) = iqtor(i) + offset
			jqtor(i) = jqtor(i) + offset
			kqtor(i) = kqtor(i) + offset
			lqtor(i) = lqtor(i) + offset
			write (*,362) iqtor(i),jqtor(i),kqtor(i),lqtor(i), qtorcod(i,:)
			!check types
			do j=1,nstates
				if(qtorcod(i,j) >0) then
					if (allocated(type_read)) then
					if(.not. type_read(qtorcod(i,j))) then
						write(*,112) 'Q-torsion', qtorcod(i,j)
						qatom_load_fep = .false.
					end if
					else
					    write(*,112) 'Q-torsion', qtorcod(i,j)
					    qatom_load_fep = .false.
					end if
				end if
			end do !j
		end do !i
361	format('atom_i atom_j atom_k atom_l    torsion type in',4(1x,a5,i2))
362	format(4(i6,1x),t47,4i8)
	end if
	if (allocated(type_read)) then
	    deallocate(type_read)
	end if

	! --- Set new impropers
	section='improper_types'
	nqcod=prm_max_enum(section, type_count)
	if(nqcod > 0 ) then
		allocate(type_read(nqcod))
	type_read(:) = .false.
		allocate(qfkimp(nqcod), qimp0(nqcod))
		write (*,450)
450	format(/,'Q-improper types:',/,'type #       force-k     imp0')
		do i=1,type_count
			if(.not. prm_get_line(line)) goto 1000
			read(line,*, err=1000) j, qfkimp(j),qimp0(j)
			type_read(j) = .true.
			write(*,225) j,qfkimp(j),qimp0(j)
			qimp0(j) = deg2rad*qimp0(j)
		end do
		write (*,*)
	end if

	section = 'change_impropers'
	nqimp = prm_count(section)
	if(nqimp > 0) then
		write (*,460) nqimp
460		format (/,'No. of changing impropers= ',i5)
		write(*,461) ('state',i,i=1,nstates)
		allocate(iqimp(nqimp), &
			jqimp(nqimp), &
			kqimp(nqimp), &
			lqimp(nqimp), &
			qimpcod(nqimp,nstates))

		do i=1,nqimp
			if(.not. prm_get_line(line)) goto 1000
			read(line,*, err=1000) iqimp(i),jqimp(i),kqimp(i),lqimp(i),qimpcod(i,:)
			iqimp(i) = iqimp(i) + offset
			jqimp(i) = jqimp(i) + offset
			kqimp(i) = kqimp(i) + offset
			lqimp(i) = lqimp(i) + offset
			write (*,462) iqimp(i),jqimp(i),kqimp(i),lqimp(i), qimpcod(i,1:nstates)
			!check types
			do j=1,nstates
				if(qimpcod(i,j) >0) then
					if (allocated(type_read)) then
					if(.not. type_read(qimpcod(i,j))) then
						write(*,112) 'Q-impsion', qimpcod(i,j)
						qatom_load_fep = .false.
					end if
					else
					    write(*,112) 'Q-impsion', qimpcod(i,j)
					    qatom_load_fep = .false.
					end if
				end if
			end do !j
		end do !i
	end if
461	format('atom_i atom_j atom_k atom_l    improper type in',4(1x,a5,i2))
462	format(4(i6,1x),t48,4i8)
	if (allocated(type_read)) then
	    deallocate(type_read)
	end if

  ! --- Read angle,torsion and improper couplings to Morse bonds
	section = 'angle_couplings'
	nang_coupl = prm_count(section)
	if(nang_coupl > 0) then
		write (*,544) nang_coupl
544		format (/,'No. of angle-Morse couplings = ',i5)
		write(*,543)
		do i=1,nang_coupl
			if(.not. prm_get_int_int(j, k)) goto 1000
			iang_coupl(1,i) = j
			iang_coupl(2,i) = k
			write (*,545) iang_coupl(:,i)
			if(j < 1 .or. j > nqangle) then
				write(*,546) 'q-angle', j
				qatom_load_fep = .false.
			end if
			if(k < 1 .or. k > nqbond) then
				write(*,546) 'q-bond', j
				qatom_load_fep = .false.
			end if
			if(bond_harmonic_in_any_state(k)) then
				qatom_load_fep = .false.
			end if
		end do
	end if
543	format ('angle_i bond_j')
545	format(i7,1x,i6)
546	format('>>>>> ERROR: ',a,' number ',i2,' does not exist.')

	section='torsion_couplings'
	ntor_coupl=prm_count(section)
	if( ntor_coupl > 0) then
		write (*,548) ntor_coupl
548		format (/,'No. of torsion-Morse couplings = ',i5)
		write(*,547)
		do i=1,ntor_coupl
			if(.not. prm_get_int_int(j, k)) goto 1000
			itor_coupl(1,i) = j
			itor_coupl(2,i) = k
			write (*,549) itor_coupl(:,i)
			if(itor_coupl(1,i) < 1 .or. itor_coupl(1,i) > nqtor) then
				write(*,546) 'q-torsion', itor_coupl(1,i)
				qatom_load_fep = .false.
			end if
			if(itor_coupl(2,i) < 1 .or. itor_coupl(2,i) > nqbond) then
				write(*,546) 'q-bond', itor_coupl(2,i)
				qatom_load_fep = .false.
			end if
			if(bond_harmonic_in_any_state(k)) then
				qatom_load_fep = .false.
			end if
		end do
547	format ('torsion_i bond_j')
549	format (i9,1x,i6)
	end if

	section='improper_couplings'
	nimp_coupl=prm_count(section)

	if(nimp_coupl > 0 ) then
		write (*,552) nimp_coupl
552		format (/,'No. of improper-Morse couplings = ',i5)
		write(*,554)
		do i=1,nimp_coupl
		    if (.not. prm_get_line(line)) goto 1000
			read(line,*, iostat=stat) j,k,l
			if(stat .eq. 0) then
				iimp_coupl(1,i) = j
				iimp_coupl(2,i) = k
				iimp_coupl(3,i) = l
			else
				read(line,*, err=1000) j, k
				iimp_coupl(1,i) = j
				iimp_coupl(2,i) = k
				iimp_coupl(3,i) = 0
			end if

			write (*,553) iimp_coupl(:,i)
			if(j < 1 .or. j > nqimp) then
				write(*,546) 'q-improper', j
				qatom_load_fep = .false.
			end if
			if(k < 1 .or. k > nqbond) then
				write(*,546) 'q-bond', j
				qatom_load_fep = .false.
			end if
			if(bond_harmonic_in_any_state(k)) then
				qatom_load_fep = .false.
			end if
		end do
	end if
554	format ('improper_i bond_j making/breaking')
553	format (i10,1x,i6,1x,i15)

	! --- Read extra shake constraints
	section='shake_constraints'
	nqshake=prm_count(section)
	if(nqshake > 0) then
		write (*,560) nqshake, ('state',i,i=1,nstates)
560		format (/,'No. of fep/evb shake contraints = ',i5,/ &
				'atom_i atom_j    distance in',5(1x,a5,i2))
		do i=1,nqshake
			if(.not. prm_get_line(line)) goto 1000
			read(line,*, err=1000) iqshake(i),jqshake(i),qshake_dist(i,1:nstates)
			iqshake(i) = iqshake(i) + offset
			jqshake(i) = jqshake(i) + offset
			write (*,580) iqshake(i),jqshake(i),qshake_dist(i,1:nstates)
		end do
580		format (i6,1x,i6,t29,5f8.2)
	end if

  ! --- Read off-diagonal hamiltonial matrix functions
	section = 'off_diagonals'
	noffd = prm_count(section)
	!always allocate (needed for energy file writing)
	allocate(offd(noffd), offd2(noffd))
	if(noffd > 0) then
		write (*,600) noffd
600		format (/,'No. of offdiagonal (Hij) functions = ',i5)
		if(noffd > max_states) then
			write(*,'(a,i2,a)') &
				'>>> Error: Too many off-diagonal functions, (max_states is ', &
				max_states, ' Aborting.'
			qatom_load_fep = .false.
			call prm_close
			return
		end if
		write (*,'(a)') 'state_i state_j atom_k atom_l     Aij mu_ij'

		do i=1,noffd
			if(.not. prm_get_line(line)) goto 1000
			read(line,*, err=1000) offd(i)%i,offd(i)%j,offd2(i)%k, &
				offd2(i)%l,offd2(i)%A, offd2(i)%mu
			write (*,620) offd(i)%i,offd(i)%j,offd2(i)%k,&
				offd2(i)%l,offd2(i)%A, offd2(i)%mu
			!check states
			if(offd(i)%i < 1 .or. offd(i)%i > nstates .or. &
				offd(i)%j < 1 .or. offd(i)%j > nstates .or. &
				offd(i)%i == offd(i)%j) then
				write(*,622)
				qatom_load_fep = .false.
			end if
			!check q-atom numbers
			if(offd2(i)%k < 1 .or. offd2(i)%k > nqat .or. &
				offd2(i)%l < 1 .or. offd2(i)%l > nqat) then
				write(*,148) offd2(i)%k, offd2(i)%l
				qatom_load_fep = .false.
			else if(iqseq(offd2(i)%k) == 0 .or. iqseq(offd2(i)%l) == 0) then
				write(*,148) offd2(i)%k,offd2(i)%l
				qatom_load_fep = .false.
			end if
		end do
	end if
620	format (i7,1x,3i7,f8.2,f6.2)
622	format('>>>>> ERROR: Invalid combination of states: ',2i5)

	!************Softcore section************  MPA
	if (allocated(sc_lookup))    deallocate (sc_lookup)
	allocate (alpha_max(nqat,nstates),sc_lookup(nqat,natyps+nqat,nstates))

	!default is no softcore, i.e. alpha=zero (sc_lookup = 0)
	sc_lookup(:,:,:)=0.0

	section = 'softcore'
	if(.not. prm_open_section(section)) then
		write (*,1630)
	else
		write(*,'(a)') 'Reading softcore section'
		if(.not. qvdw_flag) then
			write(*,'(a)') '>>>>> ERROR: Q-atom types must be redefined in "change_atoms" section'
			qatom_load_fep = .false.
			return
		end if

		numsoftlines = prm_count(section)
		if(numsoftlines /= nqat) then
			write(*,'(a)') '>>>>> ERROR: Alpha must be given for every Q-atom!'
			qatom_load_fep = .false.
			return
		else
			write (*,1625)
1625			format (/,'Assigning softcore to all Q atoms:')
1626			format('Q atom    softcore alpha in',6(2x,a5,i2))
1627			format('Q atom      max potenial in',6(2x,a5,i2))

			if (softcore_use_max_potential) then
				write (*,'(a)') ('Using values in FEP file as desired maximum vdW potentials at r=0')
				write (*,1627) ('state',i,i=1,nstates)
			else
				write (*,'(a)') ('Using values in FEP file as explicit alphas')
				write (*,1626) ('state',i,i=1,nstates)
			end if

			do i=1,nqat   !read all the softcore max_alphas
				!reading should not fail at this stage - get_line was called by prm_count!
				if(.not. prm_get_line(line)) goto 1000
				read (line,*, err=1000) tmpindex, alpha_max(tmpindex,1:nstates)

				write (*,'(i6,t25,6(f8.2,1x))') tmpindex, alpha_max(tmpindex,1:nstates)

			end do

			do i=1,nqat  !make the sc_lookup table
				do i2=1,nstates
					do j=1,natyps !do q-surroundings first

						if (softcore_use_max_potential) then
							sc_aq = qavdw(qiac(i,i2),1)
							sc_bq = qbvdw(qiac(i,i2),1)
							sc_aj = iaclib(j)%avdw(1)
							sc_bj = iaclib(j)%bvdw(1)
							if (alpha_max(i,i2) /= 0) then
								if (ivdw_rule == 1) then !geometric vdw rule
									sc_lookup(i,j,i2) = (-sc_bq*sc_bj+sqrt(sc_bq*sc_bq*sc_bj* &
										sc_bj+4*alpha_max(i,i2)*sc_aq*sc_aj))/(2*alpha_max(i,i2))
								else !arithmetic vdw rule. OBS some epsilons (q atom epsilons, sc_bq)
										!	have not been square rooted yet. We'll take this into account
										!   when calculating the sc_lookup
									sc_lookup(i,j,i2) = (-2*sqrt(sc_bq)*sc_bj+2*sqrt(sc_bq*sc_bj**2+ &
									alpha_max(i,i2)*sqrt(sc_bq)*sc_bj))*(sc_aq+sc_aj)**6/(2*alpha_max(i,i2))
								end if
							end if
						else  !user has not requested alpha calculation, each q-atom has the same alpha for every atom type
							sc_lookup(i,j,i2) = alpha_max(i,i2)
						end if


					end do

					do j=1,nqat  !now do q-q
						if ((alpha_max(i,i2) .gt. 1E-6) .or. (alpha_max(j,i2) .gt. 1E-6)) then !if both alphas are 0 then no need to calculate alphas
							sc_aq = qavdw(qiac(i,i2),1)
							sc_bq = qbvdw(qiac(i,i2),1)
							sc_aj = qavdw(qiac(j,i2),1)
							sc_bj = qbvdw(qiac(j,i2),1)

							if (softcore_use_max_potential) then  !use the smallest alpha_max of the two q atoms
								alpha_max_tmp = min ( alpha_max(i,i2), alpha_max(j,i2) )

								if ((alpha_max(i,i2) == 0) .or. (alpha_max(j,i2) == 0)) then !unless one of them is zero
									alpha_max_tmp = max ( alpha_max(i,i2), alpha_max(j,i2) )
								end if



							else  !use the largest alpha_max if we're using plain alphas
								alpha_max_tmp = max ( alpha_max(i,i2), alpha_max(j,i2) )
							end if


							if (softcore_use_max_potential) then

								if (ivdw_rule == 1) then !geometric vdw rule
									sc_lookup(i,j+natyps,i2) = (-sc_bq*sc_bj+ &
										sqrt(sc_bq*sc_bq*sc_bj*sc_bj+ &
										4*alpha_max_tmp*sc_aq*sc_aj))/(2*alpha_max_tmp)
								else !arithmetic vdw rule   OBS some epsilons (q atom epsilons, sc_bq and sc_bj)
										!	have not been square-rooted yet. We'll take this into account
										!   when calculating the sc_lookup
									sc_lookup(i,j+natyps,i2) = (-2*sqrt(sc_bq*sc_bj)+ &
										2*sqrt(sc_bq*sc_bj+alpha_max_tmp*sqrt(sc_bq*sc_bj)))* &
										(sc_aq+sc_aj)**6/(2*alpha_max_tmp)
								end if

							else  !user has not requested alpha calculation, each q-atom has the same alpha for every atom type
								sc_lookup(i,j+natyps,i2) = alpha_max_tmp
							end if !softcore_use_max_potential

						end if


					end do





				end do  ! states


			end do !softcore lookup table
		end if
	end if  !prm_open_section(section)   softcore


1630		format('No softcore section found. Using normal LJ potentials.')
1631		format('>>>>> Erroneous softcore section.')

	!********END*****Softcore section****END********  MPA

	!load atom groups whose non-bonded interactions are to be monitored
	section='monitor_groups'
	monitor_groups=prm_count(section)
	allocate(monitor_atom_group(monitor_groups))
	if (monitor_groups>0) then
			write (*,650) monitor_groups
			do i=1,monitor_groups
					k=0
					do while(prm_get_field(line)) !get one topology atom number at a time
						k=k+1
						read(line,*) temp_atom(k)   !read line med fritt format
					end do
					allocate(monitor_atom_group(i)%atom(k))   !k är antal element i array nr i
					monitor_atom_group(i)%atom(:)=temp_atom(1:k) + offset! kopiera temp_atom arrayens element med index 1-k till
					monitor_atom_group(i)%n=k
					write (*,660,advance='no') i
					write (*,670) monitor_atom_group(i)%atom(:)
			end do
	end if
650	format (/,'No. atom groups to monitor   = ',i5)
660	format ('group' ,i2, ':')
670	format (7i8)

	section='monitor_group_pairs'
	monitor_group_pairs=prm_count(section)   !antal rader i avd.
	allocate(monitor_group_pair(monitor_group_pairs))
	if(monitor_group_pairs > 0) then
			write (*,*)
			write (*,630) monitor_group_pairs
			write (*,'(a)') 'group_i group_j'
			do i=1,monitor_group_pairs
				if (.not. prm_get_int_int(monitor_group_pair(i)%i, monitor_group_pair(i)%j)) goto 1000   !läser in i arrayen monitor_group_pair
				write (*,640)  monitor_group_pair(i)%i, monitor_group_pair(i)%j
			end do
	end if
630	format (/,'No. of group pairs to monitor= ',i5)
640	format (i7,1x,i7)


	call prm_close
	return

	!error handling code
1000	write(*,1900) i,section
	call prm_close
	qatom_load_fep = .false.
1900	format('>>>>> ERROR: Read error at line ',i2,' of section [',a,']')
  !.......................................................................
end function qatom_load_fep

!-------------------------------------------------------------------------

logical function bond_harmonic_in_any_state(k)
	integer, intent(in)			::	k
	integer						::	st

	bond_harmonic_in_any_state = .false.
	do st=1, nstates
		if(qbnd(k)%cod(st) /= 0) then
			if(qbondlib(qbnd(k)%cod(st))%fk /= 0.) then
				!the bond types is harmonic
				write(*,547) k, st
				bond_harmonic_in_any_state = .true.
			end if
		end if
	end do
547 format('>>>>> ERROR: Bond',i3,' is harmonic in state',i2)
end function bond_harmonic_in_any_state

end module QATOM
