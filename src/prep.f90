! (C) 2000 Uppsala Molekylmekaniska HB, Uppsala, Sweden
! prep.f90
! by Johan Åqvist & John Marelius
! topology preparation, solvation, validation and PDB I/O
!TODO: precision not fixed

MODULE PREP

	use TRJ
	use PARSE
	use PRMFILE
	use INDEXER
	use PREFS
	use MASKMANIP

	IMPLICIT none

!constants
	character(*), private, parameter	::	MODULE_VERSION = '5.06'
	character(*), private, parameter	::	MODULE_DATE    = '2014-01-01'

	!library
	!max number of library entries
	integer, parameter		::	max_entry = 1000
	!max atoms in a residue (only used when reading PDB file)
	integer, parameter		::	max_atlib = 500
	integer, parameter		::	max_cgplib = 100
	integer, parameter		::	max_atcgplib = 100


	!FF parameters
	integer, parameter			::	max_old_atyps = 100
	!Extra bonds(S-S bridges etc.)
	integer, parameter			::	max_extrabnd = 100

	!Topology
	integer, parameter			::	max_lib = 1000
	integer, parameter			::	max_long = 100000

	integer, parameter			:: max_conn = 10

!default values for user-settable variables !TODO: is it should be here?
	!minimum solvent - solute heavy atom distance for solvation
	real(kind=prec), target				::	solvent_pack = 2.4_prec
	!average heavy atom number density of proteins
	real(kind=prec), target				::  rho_solute = 0.05794_prec  ! atoms / A**3
	!maximum cross-linking bond distance
	real(kind=prec), target				::	max_xlink = 2.1_prec
	character(len=200), target	::	solvent_names = ''
	!Random numbers for H generation
	integer, target				::	random_seed_solute = 179857
	integer, target				::	random_seed_solvent = 758971
!variables
!     Library information
!-----------------------------------------------------------------------


	TYPE LIB_BOND_TYPE
		integer(AI)				::	i,j
	end TYPE LIB_BOND_TYPE

	TYPE LIB_IMP_TYPE !for explicit improper definitions +/- and 4-char. atom names
		character(len=1+4)		::	i, j, k, l
	end TYPE LIB_IMP_TYPE

	type LIB_RULE_TYPE
		integer					::	kind
		integer					::	atom(4)
		real(kind=prec)					::	value
	end type LIB_RULE_TYPE

	integer, parameter		::	BUILD_RULE_TORSION = 1

	TYPE LIB_ENTRY_TYPE
		integer					::	nat, nbnd, nimp, ncgp, nrules
		integer					::	head, tail
		logical					::	HETATM, solvent
		real(kind=prec)					::	density
		character(len=4)		::	nam
		character(len=8)		::	SYBYLTYPE
		character(len=4), pointer	::	atnam(:)
		character(len=KEYLENGTH), pointer :: tac_lib(:)
		real(kind=prec), pointer			::	crg_lib(:)
		type(LIB_RULE_TYPE), pointer::rules(:)
		type(LIB_BOND_TYPE), pointer::bnd(:)
		type(LIB_IMP_TYPE), pointer::imp(:)
		integer(AI), pointer	::	natcgp(:)
		integer(AI), pointer	::	switch(:)
		integer(AI), pointer	::	atcgp(:,:)
	end type LIB_ENTRY_TYPE

	integer						::	nlibres
	type(LIB_ENTRY_TYPE), target::	lib(max_entry)

! topology generation flags
!-----------------------------------------------------------------------
	LOGICAL						::	have_prm_file_name = .false.
	logical						::	have_solute_sphere = .false.
	logical						::	have_title = .false.
	logical						::	have_solvent_boundary = .false.
	logical						::  boundary_set = .false.

!     FF parameter information
! --- Bond parameters
	integer						::	nbnd_prm
	integer						::	nbnd_types
	TYPE BOND_TYPE_TYPE
		character(len=KEYLENGTH)::	taci, tacj
		integer					::	cod
	end type BOND_TYPE_TYPE
	TYPE BOND_PRM_TYPE
		type(BONDLIB_TYPE)		::	prm
		character(len=2)		::	SYBYLtype
	end type BOND_PRM_TYPE
	type(BOND_TYPE_TYPE), allocatable::	bnd_types(:)
	type(BOND_PRM_TYPE), allocatable ::	bnd_prm(:)
! --- Angle parameters
	TYPE ANGLE_TYPE_TYPE
		character(len=KEYLENGTH)::	taci, tacj, tack
		integer					::	cod
	end type ANGLE_TYPE_TYPE
	type(ANGLE_TYPE_TYPE), allocatable::	ang_types(:)
	type(ANGLIB_TYPE), allocatable::	ang_prm(:)
	integer						::	nang_types
	integer						::	nang_prm

! --- Torsion parameters
	TYPE TORSION_TYPE_TYPE
		character(len=KEYLENGTH)::	taci, tacj, tack, tacl
		integer					::	cod
	end type TORSION_TYPE_TYPE
	type(TORSION_TYPE_TYPE), allocatable :: tor_types(:)
	type(TORLIB_TYPE), allocatable :: tor_prm(:)
	integer						::	ntor_types
	integer						::	ntor_prm
	type TOR_CODES
		integer					::	ncod
		integer					::	cod(10)
	end type TOR_CODES

! --- Improper torsion parameters
	TYPE IMP_PRM_TYPE
		character(len=KEYLENGTH)::	taci, tacj, tack, tacl
		type(IMPLIB_TYPE)		::	prm
	end type IMP_PRM_TYPE
	type(IMP_PRM_TYPE), allocatable :: imp_prm(:)
	integer						::	nimp_prm
	logical						::	imp_explicit !explicit or automatic defs

! --- Co-ordinate information
	integer						::	natom, nat_wat, nwat
	character(len=100)			::	coord_source = ''
	character(len=80)			::	auto_name = ''
	integer, parameter			::	trj_unit = 17
	integer						::	trj_frame = 0
	character(len=180)			::	trj_filnam
	type(MASK_TYPE)				::	mask

! --- Extra bonds(S-S bridges etc.)
	integer						::	nextrabnd
	type(BOND_TYPE)				::	extrabnd(max_extrabnd)

!     Things needed for topology generation
!-----------------------------------------------------------------------
	logical, allocatable		::	makeH(:)
	integer, allocatable		::	nconn(:)

	integer, allocatable		::	iconn(:,:)

	real(kind=prec)						::	pi, deg2rad
	!temporary storage for solvent coordinates
	real(kind=prec), allocatable		::	xw(:,:,:) !solvent coordinates zyx,atom,molecule
	!flag for solvent molecules not clashing with any solute heavy atom
	logical, allocatable		::	keep(:)
	!residue code for solvent
	integer						::	irc_solvent
	!residue name for solvent
	character(len=4)			::	solvent_name


! --- Flag to keep track of missing parameters without bailing out
	LOGICAL						::	topo_ok = .false.
	logical						::	ff_ok = .false.

!	Memory management
	integer, private			::	alloc_status
	!private subroutine
	private						::	check_alloc

! --- New stuff to make connection lists right between molecules
	type RETTYPE
		logical				:: new
		integer				:: first
	end type RETTYPE

! --- Types for creation of larger solvent molecules
	TYPE MISSING_TYPE
		integer                                 ::      bonds(4), angles(12)
		logical                                 ::      bond_set(4), angle_set(12)
	end TYPE MISSING_TYPE
	TYPE MISSING_BOND_TYPE
		integer					::	i,j,cod
	end TYPE MISSING_BOND_TYPE
	TYPE MISSING_ANGLE_TYPE
		integer                                 ::      i,j,k,cod
	end TYPE MISSING_ANGLE_TYPE


contains

!-----------------------------------------------------------------------

subroutine prep_startup
	logical						::	l
	!initialise used modules
	call topo_startup			! empty
	call prmfile_startup		! empty
	call parse_startup			! some code
	call index_startup			! empty

	!initialize preference module & add user-settable variables
	call pref_initialize() !	! allocates mem for p(:)

	l = pref_add('solvent_pack', rval=solvent_pack)
	l = pref_add('solvent_names', sval=solvent_names)
	l = pref_add('solute_density', rval=rho_solute)
	l = pref_add('random_seed_solute', ival=random_seed_solute)
	l = pref_add('random_seed_solvent', ival=random_seed_solvent)
	pi = 4.0_prec * atan(1.0_prec)
	deg2rad = pi/180.0_prec

	!make sure all arrays are allocated so that solvation can be done
	!without reading solute pdb file
	call allocate_for_pdb(1,1,1)
	call cleartop		! clears pdb filename and prm name
	CALL clearpdb		!get rid of old PDB data
end subroutine prep_startup

!-----------------------------------------------------------------------

subroutine allocate_for_pdb(atoms, residues, molecules)
!arguments
	integer, intent(in)			::	atoms, residues, molecules
	allocate(xtop(3*atoms), heavy(atoms), makeH(atoms), &
		res(residues), istart_mol(molecules), &
		stat=alloc_status)

	if(alloc_status /= 0) then
		write(*,*) 'ERROR: Out of memory when allocating arrays for new atoms'
		stop 255
	end if

end subroutine allocate_for_pdb

!-----------------------------------------------------------------------

subroutine prep_shutdown
	call topo_deallocate

	call index_shutdown
	call mask_finalize(mask)
	call trj_shutdown
end subroutine prep_shutdown

!-----------------------------------------------------------------------

subroutine clearlib
	lib_files = ''
	nlibres = 0
end subroutine clearlib

!-----------------------------------------------------------------------

subroutine check_alloc(message)
!arguments
	character(*) message

	if(alloc_status .ne. 0) then
		write(*,*) &
			'>>> Out of memory trying to allocate ', message
		stop 255
	end if
end subroutine check_alloc

!-----------------------------------------------------------------------

subroutine addbond
! *** local variables
	integer						::	ia, ja
	real(kind=prec)						::	bond_dist, rij(3)
	character(len=80)			::	reply
	integer						::	readstat

	call get_string_arg(reply, &
		'-----> First atom (number or residue:atom_name): ')
	if(scan(reply, ':') > 0) then !got res:at
		ia=get_atom_from_descriptor(reply)
		if(ia== 0) then
			write(*,910) trim(reply)
910			format('>>>>> ERROR: Could not find atom ',a)
			return
		end if
	else !got atom number
		read(reply, *, iostat=readstat) ia
		if(readstat /= 0) then
			write(*,910) trim(reply)
			return
		end if
	end if
	if(ia > nat_pro .or. ia < 1) then
		write(*,'(a)') 'Out of range!'
		return
	end if
	call get_string_arg(reply, &
		'-----> Second atom (number or residue:atom_name): ')
	if(scan(reply, ':') > 0) then !got res:at
		ja=get_atom_from_descriptor(reply)
		if(ja== 0) then
			write(*,910) trim(reply)
			return
		end if
	else !got atom number
		read(reply, *, iostat=readstat) ja
		if(readstat /= 0) then
			write(*,910) trim(reply)
			return
		end if
	end if
	if(ja > nat_pro .or. ja < 1) then
		write(*,'(a)') 'Out of range!'
		return
	elseif(ja == ia) then
		write(*,'(a)') 'First and second atom must be different!'
		return
	end if
	rij(:) = xtop(3*ia-2:3*ia) - xtop(3*ja-2:3*ja)
	bond_dist = sqrt(dot_product(rij,rij))
	if(bond_dist >  max_xlink) then
		write(*,100) bond_dist
100		format('>>> WARNING: Bond distance is',f6.2,' A.')
		call get_string_arg(reply, '-----> Make bond anyway (y/n)? ')
		call upcase(reply)
		if(reply /= 'Y' .and. reply /= 'YES') then
			return
		end if
	end if
	if(nextrabnd < max_extrabnd) then
		nextrabnd = nextrabnd+1
		extrabnd(nextrabnd)%i = ia
		extrabnd(nextrabnd)%j = ja
	else
		write(*,900) max_extrabnd
900		format('>>>>> ERROR: Too many extra bonds. (Max ',i5,')')
	end if

end subroutine addbond

!-----------------------------------------------------------------------

subroutine clearbond
	!clear extra bonds
	nextrabnd = 0
end subroutine clearbond

!-----------------------------------------------------------------------

subroutine angle_ene(emax, nlarge, av_ene)
! *** local variables
	integer i, j, k, ia, ic, istart, iend, i3, j3, k3, nlarge
	real(kind=prec) rji(3), rjk(3), bji, bjk, scp, angle, da, ae, dv, f1, di(3),&
	dk(3)
	real(kind=prec) emax, av_ene


	istart = 1
	iend = nangles
	nlarge = 0
	av_ene = zero

	do ia = istart, iend

	i = ang(ia)%i
	j = ang(ia)%j
	k = ang(ia)%k
	ic = ang(ia)%cod
	if(ic == 0) then !missing parameter
		write(*, '(4i5,1x,a)') ia, i, j, k, 'MISSING PARAMETERS'
		cycle
	end if
	i3 = i * 3 - 3
	j3 = j * 3 - 3
	k3 = k * 3 - 3
	rji(1) = xtop(i3 + 1) - xtop(j3 + 1)
	rji(2) = xtop(i3 + 2) - xtop(j3 + 2)
	rji(3) = xtop(i3 + 3) - xtop(j3 + 3)
	rjk(1) = xtop(k3 + 1) - xtop(j3 + 1)
	rjk(2) = xtop(k3 + 2) - xtop(j3 + 2)
	rjk(3) = xtop(k3 + 3) - xtop(j3 + 3)
	bji = sqrt(rji(1) **2 + rji(2) **2 + rji(3) **2)
	bjk = sqrt(rjk(1) **2 + rjk(2) **2 + rjk(3) **2)
	scp =(rji(1) * rjk(1) + rji(2) * rjk(2) + rji(3) * rjk(3) )
	scp = scp /(bji * bjk)
	if(scp>one) scp = one
	if(scp< -one) scp = -one
	angle = acos(scp)
	da = angle - anglib(ic)%ang0 * pi / 180.0_prec
	ae = 0.5_prec * anglib(ic)%fk * da**2.0_prec
	av_ene = av_ene+ae

	if(ae>emax) then
		nlarge = nlarge+1
		write( * , '(5i5,4f8.2)') ia, i, j, k, ic, anglib(ic)%fk ,&
			anglib(ic)%ang0, angle * 180.0_prec / pi, ae
	endif

	enddo

	if(nangles/=0) av_ene = av_ene / real(nangles,kind=prec)

!.......................................................................
end subroutine angle_ene

!-----------------------------------------------------------------------

integer function anglecode(taci, tacj, tack)
!arguments
	character(*), intent(in)		::	taci, tacj, tack
! *** local variables
	integer i

	character(len=KEYLENGTH)		::	ti, tj, tk

	ti = wildcard_tac(taci)
	tj = wildcard_tac(tacj)
	tk = wildcard_tac(tack)

!.......................................................................
	do i = 1, nang_types
		if((ti == ang_types(i)%taci .and. tj == ang_types(i)%tacj &
			.and. tk == ang_types(i)%tack) .or. (ti==ang_types(i)%tack &
			.and. tj == ang_types(i)%tacj .and. tk == ang_types(i)%taci)) then
				anglecode = ang_types(i)%cod
				return
		endif
	enddo

	anglecode = 0
	write( * , '(a,3(1x,a8))') '>>> Missing angle type for atom types' &
		, ti, tj, tk
	topo_ok = .false.
!.......................................................................
end function anglecode

!-----------------------------------------------------------------------

subroutine bond_ene(emax, nlarge, av_ene)
! *** local variables
	integer i, j, ib, ic, istart, iend, i3, j3, nlarge
	real(kind=prec) rij(3), b, db, be, dv, emax, av_ene

!.......................................................................

	istart = 1
	iend = nbonds
	nlarge = 0
	av_ene = zero

	do ib = istart, iend

	i = bnd(ib)%i
	j = bnd(ib)%j
	ic = bnd(ib)%cod
	if(ic == 0) then !missing parameter
		write(*, '(3i5,1x,a)') ib, i, j,  'MISSING PARAMETERS'
		cycle
	end if

	i3 = i * 3 - 3
	j3 = j * 3 - 3
	rij(1) = xtop(j3 + 1) - xtop(i3 + 1)
	rij(2) = xtop(j3 + 2) - xtop(i3 + 2)
	rij(3) = xtop(j3 + 3) - xtop(i3 + 3)
	b = sqrt(rij(1) **2 + rij(2) **2 + rij(3) **2)
	db = b - bondlib(ic)%bnd0
	be = 0.5_prec * bondlib(ic)%fk * db**2.0_prec
	av_ene = av_ene+be

	if(be>emax) then
		nlarge = nlarge+1
		write( * , '(4i5,4f8.2)') ib, i, j, ic, bondlib(ic)%fk,&
			bondlib(ic)%bnd0, b, be
	endif

	enddo

	if(nbonds/=0) av_ene = av_ene / real(nbonds,kind=prec)

	return
!.......................................................................
end subroutine bond_ene

!-----------------------------------------------------------------------

character(KEYLENGTH) function wildcard_tac(taci)
!arguments
	character(*), intent(in)	::	taci
! *** local variables
	integer						::	hash_pos
	!find # sign and truncate there - will match anything after
	hash_pos = scan(taci, '#')
	if(hash_pos > 0) then
		wildcard_tac = taci(1:hash_pos-1)
	else
		wildcard_tac = taci
	endif
end function wildcard_tac

integer function bondcode(taci, tacj)
!arguments
	character(*), intent(in)	::	taci, tacj
! *** local variables
	character(len=KEYLENGTH)		::	ti, tj
	integer i
	integer						::	leni, lenj
	integer						::	score, max_score = 0, max_cod
!.......................................................................

	ti = wildcard_tac(taci)
	tj = wildcard_tac(tacj)

	do i = 1, nbnd_types
		if(ti==bnd_types(i)%taci .and. tj==bnd_types(i)%tacj &
			.or. tj==bnd_types(i)%taci .and. ti==bnd_types(i)%tacj) then

			bondcode = bnd_types(i)%cod
			return
		endif
	enddo

	write( * , '(a,2(1x,a8))') '>>> Missing bond type for atom types' &
		, ti, tj
	topo_ok = .false.
	bondcode = 0

!.......................................................................
end function bondcode

!-----------------------------------------------------------------------

subroutine changeimp
! *** local variables
	integer ii, ip, i, j, k, l, noptimp, nchange, nlarge
	real(kind=prec) emax, av_ene

!.......................................................................

	write( * , '(/,a)') '	   Two options:  1. change specified impropers'
	write( * , '(a)') '					 2. change all with E > Emax(k <---> l)'
	CALL prompt('-----> Which option [1/2]: ')
	READ( *, * ) noptimp

	if(noptimp==1) then
		CALL prompt('-----> No. of impropers to be changed = ')
		READ( *, * ) nchange
		do ii = 1, nchange
		CALL prompt('-----> Give: Impr-no., i, j, k, l : ')
		READ( *, * ) ip, i, j, k, l
		imp(ip)%i = i
		imp(ip)%j = j
		imp(ip)%k = k
		imp(ip)%l = l
		enddo
	endif

	if(noptimp==2) then
		CALL prompt('-----> Give an energy threshold [kcal/mol] = ')
		READ( *, * ) emax
		nchange = 0
		CALL impr_ene(emax, nchange, av_ene, 1)

		write( * , '(a,f8.2)') '< impr. tors energy > = ', av_ene
	write( * , '(a,i5)') 'No. of changed imps   = ', nchange

		CALL impr_ene(emax, nlarge, av_ene, 0)

	write( * , '(a,i5)') 'No. above Emax now	= ', nlarge
		write( * , '(a,f8.2)') '< impr. tors energy > = ', av_ene
	endif

	write( *, * )

	return
!.......................................................................
end subroutine changeimp

!-----------------------------------------------------------------------

subroutine checkangs
! *** local variables
	integer nlarge
	real(kind=prec) emax, av_ene

!.......................................................................

	write( *, * )
	emax = get_real_arg('-----> Give an energy threshold [kcal/mol] = ')

	write( * , '(/,a)') 'angle  atoms i--j--k code force-k   ang_0   angle  energy'

	CALL angle_ene(emax, nlarge, av_ene)

	write( * , '(a,i5)') 'No. above Emax   = ', nlarge
	write( * , '(a,f8.2)') '< angle energy > = ', av_ene

	write( *, * )

end subroutine checkangs

!-----------------------------------------------------------------------

subroutine checkbonds
! *** local variables
	integer nlarge
	real(kind=prec) emax, av_ene

!.......................................................................

	write( *, * )
	emax = get_real_arg('-----> Give an energy threshold [kcal/mol] = ')

	write( * , '(/,a)') ' bond atoms i-j code force-k  bond_0   dist.  energy'

	CALL bond_ene(emax, nlarge, av_ene)

	write( * , '(a,i5)') 'No. above Emax   = ', nlarge
	write( * , '(a,f8.2)') '< bond energy >  = ', av_ene

	write( *, * )

end subroutine checkbonds

!-----------------------------------------------------------------------

subroutine checkimps
! *** local variables
	integer nlarge
	real(kind=prec) emax, av_ene

!.......................................................................

	write( *, * )
	emax = get_real_arg('-----> Give an energy threshold [kcal/mol] = ')

	write( * , '(/,a)') &
		'impr atoms i--j--k--l    code force-k   imp_0     phi  energy'

	CALL impr_ene(emax, nlarge, av_ene, 0)

	write( * , '(a,i5)') 'No. above Emax		= ', nlarge
	write( * , '(a,f8.2)') '< impr. tors energy > = ', av_ene

	write( *, * )

	return
!.......................................................................
end subroutine checkimps

!-----------------------------------------------------------------------

subroutine checktors
! *** local variables
	integer nlarge
	real(kind=prec) emax, av_ene

!.......................................................................

	write( *, * )
	emax = get_real_arg('-----> Give an energy threshold [kcal/mol] = ')

	write( * , '(/,a,a)') ' tors	atoms i--j--k--l code force-k   mult.   delta', '	 phi  energy'

	CALL tors_ene(emax, nlarge, av_ene)

	write( * , '(a,i5)') 'No. above Emax	 = ', nlarge
	write( * , '(a,f8.2)') '< torsion energy > = ', av_ene

	write( *, * )

	return
!.......................................................................
end subroutine checktors

!-----------------------------------------------------------------------

function cross_product(a,b)
	real(kind=prec)						::	a(3), b(3), cross_product(3)

	cross_product(1) = a(2) * b(3) - a(3) * b(2)
	cross_product(2) = a(3) * b(1) - a(1) * b(3)
	cross_product(3) = a(1) * b(2) - a(2) * b(1)
end function cross_product

!-----------------------------------------------------------------------


subroutine xlink
!add cross-linking bonds like SS-bridges
!locals
	integer						::	ires, jres, iat, jat, i, j, b
	integer						::	irc, jrc
	real(kind=prec)						::	d2
	character(len=1)			::	reply

	if(.not. check_residues()) then
		return
	end if

	if(.not. topo_ok) then
		write(*,90)
90		format('Making solute bond list.')
		call make_solute_bonds
	end if


	write(*,92) max_xlink
92	format('Searching for non-bonded heavy atom pairs closer than', f6.2,' A.')

	!loop over all residue pairs
	do ires = 1, (nres_solute-1)
		irc = res(ires)%irc
		do jres = ires+1, nres_solute
		jrc = res(jres)%irc
			do iat=1, lib(irc)%nat
				i = iat + res(ires)%start-1
				if(.not. heavy(i)) cycle !skip hydrogens
jloop:			do jat=1, lib(jrc)%nat-1
					j = jat+res(jres)%start-1
					if(.not. heavy(j)) cycle !skip hydrogens
					!check distance
					d2=(xtop(3*i-2)-xtop(3*j-2))**2+ &
						(xtop(3*i-1)-xtop(3*j-1))**2+(xtop(3*i)-xtop(3*j))**2
					if(d2 < max_xlink**2) then
						!found atom pair within limit Å, loop over all bonds
						do b=1, nbonds
							if(bnd(b)%i == i .and. bnd(b)%j == j) cycle jloop
							if(bnd(b)%i == j .and. bnd(b)%j == i) cycle jloop
						end do
						!By now we know the atoms are not bonded
						write(*,100, advance='no') i, trim(lib(irc)%nam), &
							ires, trim(lib(irc)%atnam(iat)), &
							j, trim(lib(jrc)%nam), jres, &
							trim(lib(jrc)%atnam(jat)), sqrt(d2)
						!add
						read(*,*) reply
						if(reply=='N' .or. reply=='n') cycle
						if(nextrabnd < max_extrabnd) then
							nextrabnd = nextrabnd+1
							extrabnd(nextrabnd)%i = i
							extrabnd(nextrabnd)%j = j
						else
							write(*,900) max_extrabnd
900		format('>>>>> ERROR: Too many extra bonds. (Max ',i5,')')
						end if

					end if
				end do jloop
			end do
		end do
	end do

	!issue a warning if the topology was already made.
	if(nextrabnd > 0) then
		write(*,110) nextrabnd
110		format(i3,' bonds were added.')
		if(topo_ok) then
			write(*,120)
		end if
120		format('You need to regenerate the topology with the added bonds.')
	end if

100	format('Atoms',i6,1x,a,1x,i5,':',a,' and',i6,1x,a,1x,i5,':',a,' are',f5.2,' A apart. Connect [Y/n]?')
end subroutine xlink

!-----------------------------------------------------------------------

integer function impcode(taci, tacj, tack, tacl)
!arguments
	character(KEYLENGTH)			::	taci, tacj, tack, tacl
! *** local variables
	integer i
	character(KEYLENGTH)			::	ti, tj, tk, tl
	logical							::  found1, found2, found3, found4
	ti = wildcard_tac(taci)
	tj = wildcard_tac(tacj)
	tk = wildcard_tac(tack)
	tl = wildcard_tac(tacl)



!.......................................................................

	!find a match with four types (i and j and k and l)
!	do i = 1, nimp_prm
!		if( (ti == imp_prm(i)%taci .or. &
!			 ti == imp_prm(i)%tacj .or. &
!			 ti == imp_prm(i)%tack .or. &
!			 ti == imp_prm(i)%tacl) .and. &
!
!			(tj == imp_prm(i)%taci .or. &
!			 tj == imp_prm(i)%tacj .or. &
!			 tj == imp_prm(i)%tack .or. &
!			 tj == imp_prm(i)%tacl) .and. &
!
!			(tk == imp_prm(i)%taci .or. &
!			 tk == imp_prm(i)%tacj .or. &
!			 tk == imp_prm(i)%tack .or. &
!			 tk == imp_prm(i)%tacl) .and. &
!
!			(tl == imp_prm(i)%taci .or. &
!			 tl == imp_prm(i)%tacj .or. &
!			 tl == imp_prm(i)%tack .or. &
!			 tl == imp_prm(i)%tacl)) then
!
!			impcode = i
!			return
!		endif
!	enddo


	!find a match with four types (i and j and k and l)
	do i = 1, nimp_prm
			found1 = .false.
			found2 = .false.
			found3 = .false.
			found4 = .false.

		    if (ti == imp_prm(i)%taci) then
				found1 = .true.
			else if (ti == imp_prm(i)%tacj) then
				found2 = .true.
			else if (ti == imp_prm(i)%tack) then
				found3 = .true.
			else if (ti == imp_prm(i)%tacl) then
				found4 = .true.
			endif

			if ((.not. found1) .and. (tj == imp_prm(i)%taci)) then
				found1 = .true.
			else if ((.not. found2) .and. (tj == imp_prm(i)%tacj)) then
				found2 = .true.
			else if ((.not. found3) .and. (tj == imp_prm(i)%tack)) then
				found3 = .true.
			else if ((.not. found4) .and. (tj == imp_prm(i)%tacl)) then
				found4 = .true.
			endif

			if ((.not. found1) .and. (tk == imp_prm(i)%taci)) then
				found1 = .true.
			else if ((.not. found2) .and. (tk == imp_prm(i)%tacj)) then
				found2 = .true.
			else if ((.not. found3) .and. (tk == imp_prm(i)%tack)) then
				found3 = .true.
			else if ((.not. found4) .and. (tk == imp_prm(i)%tacl)) then
				found4 = .true.
			endif

			if ((.not. found1) .and. (tl == imp_prm(i)%taci)) then
				found1 = .true.
			else if ((.not. found2) .and. (tl == imp_prm(i)%tacj)) then
				found2 = .true.
			else if ((.not. found3) .and. (tl == imp_prm(i)%tack)) then
				found3 = .true.
			else if ((.not. found4) .and. (tl == imp_prm(i)%tacl)) then
				found4 = .true.
			endif

			if (found1 .and. found2 .and. found3 .and. found4) then
				impcode = i
				return
			endif
	enddo















	!now we have to handle force fields differently
	select case (ff_type)
		case (FF_GROMOS, FF_AMBER)
		!find a match with three types (j and k and l) or
		!(i and j and k)
		do i = 1, nimp_prm

			if(imp_prm(i)%taci == '' .and. imp_prm(i)%tacl /='') then
			! for ?  A*  B  C   ...

				if((tj == imp_prm(i)%tacj .and. &
					tk == imp_prm(i)%tack .and. &
					tl == imp_prm(i)%tacl) .or. &

					(tj == imp_prm(i)%tacj .and. &
					tk == imp_prm(i)%tacl .and. &
					tl == imp_prm(i)%tack) .or. &

                    (tj == imp_prm(i)%tacj .and. &
					ti == imp_prm(i)%tack .and. &
					tl == imp_prm(i)%tacl) .or. &

					(tj == imp_prm(i)%tacj .and. &
					ti == imp_prm(i)%tacl .and. &
					tl == imp_prm(i)%tack) .or. &

					(tj == imp_prm(i)%tacj .and. &
					ti == imp_prm(i)%tack .and. &
					tk == imp_prm(i)%tacl) .or. &

					(tj == imp_prm(i)%tacj .and. &
					ti == imp_prm(i)%tacl .and. &
					tk == imp_prm(i)%tack)) then

					impcode = i
					return
				endif
			elseif(imp_prm(i)%taci /= '' .and. imp_prm(i)%tacl =='') then
				! and for A*  B  C  ?  (2nd half) ...

				if((tj == imp_prm(i)%taci .and. &
					tk == imp_prm(i)%tacj .and. &
					tl == imp_prm(i)%tack) .or. &

					(tj == imp_prm(i)%taci .and. &
					tk == imp_prm(i)%tack .and. &
					tl == imp_prm(i)%tacj) .or. &

					(tj == imp_prm(i)%taci .and. &
					tj == imp_prm(i)%tacj .and. &
					ti == imp_prm(i)%tack) .or. &

					(tj == imp_prm(i)%taci .and. &
					tj == imp_prm(i)%tack .and. &
					ti == imp_prm(i)%tacj) .or. &

					(tj == imp_prm(i)%taci .and. &
					tk == imp_prm(i)%tacj .and. &
					ti == imp_prm(i)%tack) .or. &

					(tj == imp_prm(i)%taci .and. &
					tk == imp_prm(i)%tack .and. &
					ti == imp_prm(i)%tacj)) then

					impcode = i
					return
				endif
			endif
		enddo

		!find a match with two types (j and  (i or k or l))
		do i = 1, nimp_prm

			if(imp_prm(i)%taci == '' .and. imp_prm(i)%tacl =='') then

			if((tj == imp_prm(i)%tacj .and. tl == imp_prm(i)%tack) .or. &
			   (tj == imp_prm(i)%tack .and. tl == imp_prm(i)%tacj) .or. &

					   (tj == imp_prm(i)%tacj .and. tk == imp_prm(i)%tack) .or. &
			   (tj == imp_prm(i)%tack .and. tk == imp_prm(i)%tacj) .or. &

					   (tj == imp_prm(i)%tacj .and. ti == imp_prm(i)%tack) .or. &
			   (tj == imp_prm(i)%tack .and. ti == imp_prm(i)%tacj) .or. &

					   (tk == imp_prm(i)%tacj .and. tl == imp_prm(i)%tack) .or. &
			   (tk == imp_prm(i)%tack .and. tl == imp_prm(i)%tacj) .or. &

					   (tk == imp_prm(i)%tacj .and. tj == imp_prm(i)%tack) .or. &
			   (tk == imp_prm(i)%tack .and. tj == imp_prm(i)%tacj) .or. &

					   (tk == imp_prm(i)%tacj .and. ti == imp_prm(i)%tack) .or. &
			   (tk == imp_prm(i)%tack .and. ti == imp_prm(i)%tacj)) then



				impcode = i
				return
			endif
			endif
		enddo

		!find a match with only type of j atom
		do i = 1, nimp_prm
			if(imp_prm(i)%taci == '' .and. imp_prm(i)%tack =='' .and. &
				imp_prm(i)%tacl =='')then
			if(tj == imp_prm(i)%tacj) then
				impcode = i
				return
			endif
			endif
		enddo


		case (FF_CHARMM)

		!find a match with three types (j and k and l) or
		!(i and j and k)
		do i = 1, nimp_prm
			if(imp_prm(i)%taci /= '' .and. imp_prm(i)%tacl =='') then
				if((ti == imp_prm(i)%taci .and. &
					tj == imp_prm(i)%tacj .and. &
					tk == imp_prm(i)%tack) .or. &

					(ti == imp_prm(i)%taci .and. &
					tj == imp_prm(i)%tack .and. &
					tk == imp_prm(i)%tacj) .or. &

					(ti == imp_prm(i)%taci .and. &
					tj == imp_prm(i)%tacj .and. &
					tl == imp_prm(i)%tack) .or. &

					(ti == imp_prm(i)%taci .and. &
					tj == imp_prm(i)%tack .and. &
					tl == imp_prm(i)%tacj) .or. &

					(ti == imp_prm(i)%taci .and. &
					tk == imp_prm(i)%tacj .and. &
					tl == imp_prm(i)%tack) .or. &

					(ti == imp_prm(i)%taci .and. &
					tk == imp_prm(i)%tack .and. &
					tl == imp_prm(i)%tacj) .or. &

					(tl == imp_prm(i)%taci .and. &
					tj == imp_prm(i)%tacj .and. &
					tk == imp_prm(i)%tack) .or. &

					(tl == imp_prm(i)%taci .and. &
					tj == imp_prm(i)%tack .and. &
					tk == imp_prm(i)%tacj) .or. &

					(tl == imp_prm(i)%taci .and. &
					tj == imp_prm(i)%tacj .and. &
					ti == imp_prm(i)%tack) .or. &

					(tl == imp_prm(i)%taci .and. &
					tj == imp_prm(i)%tack .and. &
					ti == imp_prm(i)%tacj) .or. &

					(tl == imp_prm(i)%taci .and. &
					tk == imp_prm(i)%tacj .and. &
					ti == imp_prm(i)%tack) .or. &

					(tl == imp_prm(i)%taci .and. &
					tk == imp_prm(i)%tack .and. &
					ti == imp_prm(i)%tacj)) then

					impcode = i
					return
				endif
			endif
		enddo

		!find a match with two types (j and  (i or k or l))
		do i = 1, nimp_prm

			if(imp_prm(i)%taci == '' .and. imp_prm(i)%tacl =='') then

				if((ti == imp_prm(i)%tacj .and. tl == imp_prm(i)%tack) .or. &
				   (ti == imp_prm(i)%tack .and. tl == imp_prm(i)%tacj) .or. &

						   (ti == imp_prm(i)%tacj .and. tk == imp_prm(i)%tack) .or. &
				   (ti == imp_prm(i)%tack .and. tk == imp_prm(i)%tacj) .or. &

						   (ti == imp_prm(i)%tacj .and. tj == imp_prm(i)%tack) .or. &
				   (ti == imp_prm(i)%tack .and. tj == imp_prm(i)%tacj) .or. &

				   (tl == imp_prm(i)%tacj .and. ti == imp_prm(i)%tack) .or. &
				   (tl == imp_prm(i)%tack .and. ti == imp_prm(i)%tacj) .or. &

						   (tl == imp_prm(i)%tacj .and. tk == imp_prm(i)%tack) .or. &
				   (tl == imp_prm(i)%tack .and. tk == imp_prm(i)%tacj) .or. &

						   (tl == imp_prm(i)%tacj .and. tj == imp_prm(i)%tack) .or. &
				   (tl == imp_prm(i)%tack .and. tj == imp_prm(i)%tacj)) then

					impcode = i
					return
				endif
			endif
		enddo

		!find a match with only type of j atom
		do i = 1, nimp_prm
			if(imp_prm(i)%taci == '' .and. imp_prm(i)%tack =='' .and. &
				imp_prm(i)%tacl =='')then
				if(ti == imp_prm(i)%tacj) then
					impcode = i
					return
				endif
			endif
		enddo

	end select

	impcode = 0
	write( * , '(a,4(1x,a8))') '>>> Missing improper type for atom types', &
		 ti, tj, tk, tl
	topo_ok = .false.

!.......................................................................
end function impcode

!-----------------------------------------------------------------------

subroutine impr_ene(emax, nlarge, av_ene, mode)
! *** local variables
	integer i, j, k, l, ip, ic, i3, j3, k3, l3, nlarge, mode
	real(kind=prec) rji(3), rjk(3), rkl(3), rnj(3), rnk(3), bj, bk, scp,	&
	phi, sgn, pe, dv, rki(3), rlj(3), dp(12), arg, f1, di(3),	 &
	dl(3)
	real(kind=prec) emax, av_ene

!.......................................................................

	nlarge = 0
	av_ene = zero

	do ip = 1, nimps

		i = imp(ip)%i
		j = imp(ip)%j
		k = imp(ip)%k
		l = imp(ip)%l
		ic = imp(ip)%cod
		if(ic == 0) then !missing parameter
			write(*, '(5i5,1x,a)') ip, i, j, k, l, 'MISSING PARAMETERS'
			cycle
		end if
		i3 = i * 3 - 3
		j3 = j * 3 - 3
		k3 = k * 3 - 3
		l3 = l * 3 - 3
		rji(1) = xtop(i3 + 1) - xtop(j3 + 1)
		rji(2) = xtop(i3 + 2) - xtop(j3 + 2)
		rji(3) = xtop(i3 + 3) - xtop(j3 + 3)
		rjk(1) = xtop(k3 + 1) - xtop(j3 + 1)
		rjk(2) = xtop(k3 + 2) - xtop(j3 + 2)
		rjk(3) = xtop(k3 + 3) - xtop(j3 + 3)
		rkl(1) = xtop(l3 + 1) - xtop(k3 + 1)
		rkl(2) = xtop(l3 + 2) - xtop(k3 + 2)
		rkl(3) = xtop(l3 + 3) - xtop(k3 + 3)
		rnj(1) = rji(2) * rjk(3) - rji(3) * rjk(2)
		rnj(2) = rji(3) * rjk(1) - rji(1) * rjk(3)
		rnj(3) = rji(1) * rjk(2) - rji(2) * rjk(1)
		rnk(1) = - rjk(2) * rkl(3) + rjk(3) * rkl(2)
		rnk(2) = - rjk(3) * rkl(1) + rjk(1) * rkl(3)
		rnk(3) = - rjk(1) * rkl(2) + rjk(2) * rkl(1)
		bj = sqrt(rnj(1) **2.0_prec + rnj(2) **2.0_prec + rnj(3) **2.0_prec)
		bk = sqrt(rnk(1) **2.0_prec + rnk(2) **2.0_prec + rnk(3) **2.0_prec)
		scp =(rnj(1) * rnk(1) + rnj(2) * rnk(2) + rnj(3) * rnk(3) )/(bj * bk)
		if(scp>one) scp = one
		if(scp< -one) scp = -one
		phi = acos(scp)
		sgn = rjk(1) *(rnj(2) * rnk(3) - rnj(3) * rnk(2) ) + rjk(2)&
		*(rnj(3) * rnk(1) - rnj(1) * rnk(3) ) + rjk(3) *(rnj(1)   &
		* rnk(2) - rnj(2) * rnk(1) )
		if(sgn<zero) phi = - phi

! ---	energy
		if(imp_type == 1) then !harmonic
			arg = phi - implib(ic)%imp0 * pi / 180.0_prec
			arg = arg - 2.0_prec * pi * nint(arg /(2.0_prec * pi) )
			dv = implib(ic)%fk * arg
			pe = 0.5_prec * dv * arg
		else !periodic
			arg = 2.0_prec*phi - implib(ic)%imp0 * pi / 180.0_prec
			pe = implib(ic)%fk * (1 + cos(arg))
			dv = -2.0_prec*implib(ic)%fk * sin(arg)
		end if

		av_ene = av_ene+pe

		if(pe>emax) then
			nlarge = nlarge+1
			if(mode==0) write( * , '(6i5,5f8.2)') ip, i, j, k, l, ic, &
			implib(ic)%fk , implib(ic)%imp0 , phi * 180.0_prec / pi, pe
			if(mode==1) then
				imp(ip)%i = i
				imp(ip)%j = j
				imp(ip)%k = l
				imp(ip)%l = k
			endif
		endif

	enddo

	if(nimps/=0) av_ene = av_ene / real(nimps,kind=prec)

!.......................................................................
end subroutine impr_ene

!-----------------------------------------------------------------------

subroutine listres
! *** local variables
	integer i, iat, ires
	type(LIB_ENTRY_TYPE), pointer::	lp
!.......................................................................

	ires = get_int_arg('-----> Residue number: ')
	if(ires > nres .or. ires == 0) then
		write(*,*) 'Out of range!'
		return
	end if
	if(.not. check_residues()) return
	lp => lib(res(ires)%irc )
	write( * , '(a,i5,a,a)') 'Residue ', ires, ' is ', res(ires)%name

	write( * , '(a)') 'Atom name  tac     charge '
	do i = 1, lp%nat
		iat = res(ires)%start - 1 + i
		write( * , '(i5,1x,a4,a8,f8.4)') iat, lp%atnam(i),&
			lp%tac_lib(i), lp%crg_lib(i)
	enddo

	write( *, * )

end subroutine listres

!-----------------------------------------------------------------------

subroutine listseq
! *** local variables
	integer i

!.......................................................................

	write( * , '(a,/)') 'Sequence listing:'
	write( *, 100) (i, res(i)%name, i = 1, nres)
  100 format(9(i4,'-',a4))
	write( *, * )

end subroutine listseq

!-----------------------------------------------------------------------

subroutine make14list
! *** local variables
	integer						::	it, i, j, ioff, l




	n14nbrs = 0
	n14long = 0
	list14(:, :) = .false.

itloop:	do it = 1, ntors



		i = tor(it)%i
		j = tor(it)%l
		if(i<j) then
			ioff = j - i
			if(ioff<=max_nbr_range) then
				if(.not. list14(ioff, i)) then
					n14nbrs = n14nbrs + 1
					list14(ioff, i) = .true.
				end if
			ELSE
				do l = 1, n14long
					if(list14long(1,l) == i .and. list14long(2,l) == j) then
						!it's already in here
						cycle itloop
					end if
				end do
				n14long = n14long + 1
				list14long(1, n14long) = i
				list14long(2, n14long) = j
			endif
		ELSE
			ioff = i - j
			if(ioff<=max_nbr_range) then
				if(.not. list14(ioff, j)) then
					n14nbrs = n14nbrs + 1
					list14(ioff, j) = .true.
				end if
			ELSE
				do l = 1, n14long
					if(list14long(1,l) == j .and. list14long(2,l) == i) then
						!it's already in here
						cycle itloop
					end if
				end do
				n14long = n14long + 1
				list14long(1, n14long) = j
				list14long(2, n14long) = i
			endif
		endif
	enddo itloop

	write(*, 100) n14nbrs
	write(*, 110) n14long
100	format('Made    ',i6,' entries in the 1-4 neighbour list.')
110	format('Made    ',i6,' entries in the long-range 1-4 neighbour list.')

!.......................................................................
end subroutine make14list

!-----------------------------------------------------------------------

subroutine makeangles
! *** local variables
	integer ib, jb, ia, i, j, k, iaci, iacj, iack
	logical						::	used(nang_prm+max_extrabnd)
	integer						::	itrans(nang_prm+max_extrabnd)

!.......................................................................

	nangles = 0
	nangles_solute = 0

	do ib = 1, nbonds_solute - 1
		!skip special bonds used only to specify shake constraints
		if(bnd(ib)%cod == 0) cycle
		if(bondlib(bnd(ib)%cod)%fk == 0.) cycle
		do jb = ib + 1, nbonds_solute
			if(bnd(jb)%cod == 0) cycle
			if(bondlib(bnd(jb)%cod)%fk == 0.) cycle

			if(bnd(ib)%i == bnd(jb)%i ) then
				nangles = nangles + 1
				ang(nangles)%i = bnd(ib)%j
				ang(nangles)%j = bnd(ib)%i
				ang(nangles)%k = bnd(jb)%j
			elseif(bnd(ib)%i ==bnd(jb)%j ) then
				nangles = nangles + 1
				ang(nangles)%i = bnd(ib)%j
				ang(nangles)%j = bnd(ib)%i
				ang(nangles)%k = bnd(jb)%i
			elseif(bnd(ib)%j ==bnd(jb)%i ) then
				nangles = nangles + 1
				ang(nangles)%i = bnd(ib)%i
				ang(nangles)%j = bnd(ib)%j
				ang(nangles)%k = bnd(jb)%j
			elseif(bnd(ib)%j ==bnd(jb)%j ) then
				nangles = nangles + 1
				ang(nangles)%i = bnd(ib)%i
				ang(nangles)%j = bnd(ib)%j
				ang(nangles)%k = bnd(jb)%i
			endif

		enddo
	enddo
	nangles_solute = nangles
	write(*, 100) nangles_solute, 'solute'
100	format('Made    ',i6,1x,a,' angles.')

	do ib = nbonds_solute+1, nbonds - 1
		!skip special bonds used only to specify shake constraints
		if(bondlib(bnd(ib)%cod)%fk == 0.) cycle
		do jb = ib + 1, nbonds
			if(bondlib(bnd(jb)%cod)%fk == 0.) cycle

			if(bnd(ib)%i == bnd(jb)%i ) then
				nangles = nangles + 1
				ang(nangles)%i = bnd(ib)%j
				ang(nangles)%j = bnd(ib)%i
				ang(nangles)%k = bnd(jb)%j
			elseif(bnd(ib)%i ==bnd(jb)%j ) then
				nangles = nangles + 1
				ang(nangles)%i = bnd(ib)%j
				ang(nangles)%j = bnd(ib)%i
				ang(nangles)%k = bnd(jb)%i
			elseif(bnd(ib)%j ==bnd(jb)%i ) then
				nangles = nangles + 1
				ang(nangles)%i = bnd(ib)%i
				ang(nangles)%j = bnd(ib)%j
				ang(nangles)%k = bnd(jb)%j
			elseif(bnd(ib)%j ==bnd(jb)%j ) then
				nangles = nangles + 1
				ang(nangles)%i = bnd(ib)%i
				ang(nangles)%j = bnd(ib)%j
				ang(nangles)%k = bnd(jb)%i
			endif
		enddo
	enddo
	write(*, 100) nangles-nangles_solute, 'solvent'

	do ia = 1, nangles
		ang(ia)%cod = anglecode(tac(iac(ang(ia)%i)), tac(iac(ang(ia)%j)), tac(iac(ang(ia)%k)))
		if(ang(ia)%cod == 0) then
			write( * , * ) 'atoms: ', ang(ia)%i , ang(ia)%j , ang(ia)%k
		endif
	enddo

! --- Make a list of actual angle types to be used
	used(:) = .false.

	nangcod = 0
	do i = 1, nangles
	!	workaround to permit missing params
		if(ang(i)%cod /=0) then
			if(.not. used(ang(i)%cod)) then
				used(ang(i)%cod ) = .true.
				nangcod = nangcod+1
				itrans(ang(i)%cod ) = nangcod
				!copy all parameters for this angle type from parameter list to top
				anglib(nangcod) = ang_prm(ang(i)%cod)
				ang(i)%cod = nangcod
			ELSE
				ang(i)%cod = itrans(ang(i)%cod )
			endif
		endif
	enddo

end subroutine makeangles

!-----------------------------------------------------------------------

subroutine make_solute_bonds
	!make solute bond list
	!used only by xlink routine to make solute bond list to check new
	!bonds against

	!deallocate old bond list
	if(allocated(bnd)) deallocate(bnd)

	!allocate topology arrays (set max_atom etc
	call topo_set_max(nat_pro, max_lib, max_long) !make space for new topology
	allocate(bnd(max_bonds), stat=alloc_status)
	if(alloc_status /= 0) then
		write(*,*) 'ERROR: Out of memory when allocating bond list'
		stop 255
	end if

	nbonds = 0

	!make solute bonds
	nbonds_solute = makesomebonds(1,nres_solute)

end subroutine make_solute_bonds

!--------------------------------------------------------------------

subroutine makebonds
	!locals
	integer						::	nbonds_solvent, nbonds_extra
	integer						::	i

	!clear bonds
	nbonds = 0

	!make solute bonds
	nbonds_solute = makesomebonds(1,nres_solute)
	write(*, 100) nbonds_solute, 'solute'
100	format('Made    ',i6,1x,a,' bonds.')

	!make extra bonds
	nbonds_extra = makeextrabonds()
	!extra bonds count as solute bonds
	nbonds_solute = nbonds_solute + nbonds_extra
	write(*, 100) nbonds_extra, 'extra'


	!make solvent bonds
	nbonds_solvent = makesomebonds(nres_solute+1, nres)
	write(*, 100) nbonds_solvent, 'solvent'

end subroutine makebonds

!-----------------------------------------------------------------------

integer function makeextrabonds()
	!locals
	integer						::	ix, ib,iaci, iacj
	character(len=KEYLENGTH)	::	taci, tacj

! --- Add extra bonds like S-S bridges

	makeextrabonds = 0

extra:do ix = 1, nextrabnd
		do ib = 1, nbonds
			if((bnd(ib)%i==extrabnd(ix)%i .and. bnd(ib)%j==extrabnd(ix)%j) &
				.or. (bnd(ib)%j==extrabnd(ix)%i &
				.and. bnd(ib)%i==extrabnd(ix)%j)) then
				write(*,8) extrabnd(ix)%i, extrabnd(ix)%j
				cycle extra
			end if
		end do
		nbonds = nbonds + 1
		makeextrabonds=makeextrabonds+1
		bnd(nbonds)%i = extrabnd(ix)%i
		bnd(nbonds)%j = extrabnd(ix)%j
		taci = tac(iac(bnd(nbonds)%i ))
		tacj = tac(iac(bnd(nbonds)%j ))

		bnd(nbonds)%cod = bondcode(taci, tacj)
	end do extra

8	format('>>> WARNING: Atoms',i5,' and',i5,' are already bonded.')

end function makeextrabonds

!-----------------------------------------------------------------------

subroutine set_bondcodes
	!locals
	integer						::	ib, iaci, iacj
	logical						::	used(nbnd_prm)
	integer						::	itrans(nbnd_prm)
	character(len=keylength)	::	taci, tacj

	do ib = 1, nbonds
		iaci = iac(bnd(ib)%i)
		iacj = iac(bnd(ib)%j)
		taci = tac(iac(bnd(ib)%i))
		tacj = tac(iac(bnd(ib)%j))
		bnd(ib)%cod =  bondcode(taci, tacj)
		IF(bnd(ib)%cod == 0) then
			WRITE( * , * ) 'atoms: ', bnd(ib)%i , bnd(ib)%j
		ENDIF
	end do

	!Make a list of actual bond types to be used

	used(:) = .false.
	nbndcod = 0
	do ib = 1, nbonds
	!	workaround to permit missing param's
		if(bnd(ib)%cod /=0) then
			if(.not. used(bnd(ib)%cod)) then
				used(bnd(ib)%cod) = .true.
				nbndcod = nbndcod+1
				itrans(bnd(ib)%cod ) = nbndcod
				bondlib(nbndcod) = bnd_prm(bnd(ib)%cod)%prm
				SYBYL_bond_type(nbndcod) = bnd_prm(bnd(ib)%cod)%SYBYLtype
				bnd(ib)%cod = nbndcod
			ELSE
				bnd(ib)%cod = itrans(bnd(ib)%cod )
			endif
		endif
	enddo

end subroutine set_bondcodes

!-----------------------------------------------------------------------

integer function makesomebonds(startres, endres)
!arguments
	integer						::	startres, endres

! *** local variables
	integer						::	i, ib, ires, imol

	makesomebonds = 0

	if(startres > endres) return !if no solvent

	!determine starting molecule
	imol = 0
	do i = 1,nmol
		if(istart_mol(i) == res(startres)%start) exit
		imol = imol+1
	end do
	if(imol >= nmol) then
		write(*,'(a)') '>>>>> ERROR: Inconsistent molecule/residue start atoms'
		return
	end if

	do ires = startres, endres

! ---	start new molecule or attach to chain
		if(nmol>imol .and. res(ires)%start==istart_mol(imol + 1) )	then
			!it is a new molecule
			imol = imol + 1
		elseif(lib(res(ires-1)%irc)%tail > 0 .and. &
			lib(res(ires)%irc)%head > 0) then
			!make tail-to-head bond
			nbonds = nbonds + 1
			makesomebonds = makesomebonds + 1
			bnd(nbonds)%i = res(ires-1)%start - 1 + lib(res(ires-1)%irc)%tail
			bnd(nbonds)%j = res(ires)%start - 1 + lib(res(ires)%irc)%head
		endif

	! ---	Build bonds in the residue

		do i = 1, lib(res(ires)%irc)%nbnd

			nbonds = nbonds + 1
			makesomebonds = makesomebonds + 1
			bnd(nbonds)%i = res(ires)%start - 1 + lib(res(ires)%irc)%bnd(i)%i
			bnd(nbonds)%j = res(ires)%start - 1 + lib(res(ires)%irc)%bnd(i)%j
		enddo
	enddo
!.......................................................................
end function makesomebonds

!-----------------------------------------------------------------------

subroutine makeexlist
! *** local variables
	integer ib, ia, i, j, ioff

!.......................................................................

	nexnbrs = 0
	nexlong = 0
	listex(:, :) = .false.

	do ib = 1, nbonds

	i = bnd(ib)%i
	j = bnd(ib)%j
	if(i<j) then
		ioff = j - i
		if(ioff<=max_nbr_range) then
			nexnbrs = nexnbrs + 1
			listex(ioff, i) = .true.
		ELSE
			nexlong = nexlong + 1
			listexlong(1, nexlong) = i
			listexlong(2, nexlong) = j
		endif
	ELSE
		ioff = i - j
		if(ioff<=max_nbr_range) then
			nexnbrs = nexnbrs + 1
			listex(ioff, j) = .true.
		ELSE
			nexlong = nexlong + 1
			listexlong(1, nexlong) = j
			listexlong(2, nexlong) = i
		endif
	endif

	enddo

	do ia = 1, nangles

	i = ang(ia)%i
	j = ang(ia)%k
	if(i<j) then
		ioff = j - i
		if(ioff<=max_nbr_range) then
			nexnbrs = nexnbrs + 1
			listex(ioff, i) = .true.
		ELSE
			nexlong = nexlong + 1
			listexlong(1, nexlong) = i
			listexlong(2, nexlong) = j
		endif
	ELSE
		ioff = i - j
		if(ioff<=max_nbr_range) then
			nexnbrs = nexnbrs + 1
			listex(ioff, j) = .true.
		ELSE
			nexlong = nexlong + 1
			listexlong(1, nexlong) = j
			listexlong(2, nexlong) = i
		endif
	endif

	enddo


	write(*, 100) nexnbrs
	write(*, 110) nexlong
100	format('Made    ',i6,' entries in list of excluded neighbours.')
110	format('Made    ',i6,' entries in list of long-range excluded neighbours.')

	return
!.......................................................................
end subroutine makeexlist

!-----------------------------------------------------------------------

subroutine makehyds
!local variables
	integer						::	i, nH_solute, nH_solvent, nH_required
	real(kind=prec)						::	r
	integer						::	atom, residue

	nH_required = 0

	!reassign heavy atom flag based on masses from parameter file
	do i = 1, nat_pro
		heavy(i) = .false.
		if(iaclib(iac(i))%mass >=4.0_prec) heavy(i) = .true.
		if(makeH(i)) nH_required = nH_required+1
	enddo

	r=randm(random_seed_solute, seed_only=.true.)
	nH_solute = 0
	do residue = 1, nres_solute
		do atom = res(residue)%start, &
			res(residue)%start+lib(res(residue)%irc)%nat-1
			nH_solute = nH_solute + genH(atom, residue)
		end do
	end do
	write(*, 100) nH_solute, 'solute'
100	format('Made    ',i6,1x,a,' hydrogens.')

	!use a separate random number sequence for solvent
	r=randm(random_seed_solvent, seed_only=.true.)
	nH_solvent = 0
	do residue = nres_solute+1, nres
		do atom = res(residue)%start, &
			res(residue)%start+lib(res(residue)%irc)%nat - 1
			nH_solvent = nH_solvent + genH(atom, residue)
		end do
	end do
	write(*, 100) nH_solvent, 'solvent'
	if(nH_solute+nH_solvent/=nH_required) write( * , '(a,i5,/)') &
		'>>> Warning, no. of hydrogens required = ', nH_required
	makeH(:) = .false. !don't need to redo when making topology
end subroutine makehyds

!-----------------------------------------------------------------------

integer function genH(j, residue)
!arguments
	integer, intent(in)			::	j, residue
!locals
	real(kind=prec)						::	xj(3), xk(3)
	integer						::	ligand, H, kt, lt
	real(kind=prec)						::	old_xH(3), xH(3), V, Vtot, dV(3), dvLast(3), gamma, dVtot(3)
	real(kind=prec)						::	VtotLast
	real(kind=prec)						::	dx(3), dx_line, rms_dV
	real(kind=prec),parameter			::	convergence_criterum = 0.1_prec
	real(8),parameter			::	dV_scale = 0.025_prec
	real(8),parameter			::	max_dx = 1.0_prec !max_dx is max distance of line search step in first CG iteration (Å)
	real(kind=prec)						::	local_min = 30
	real(kind=prec)						::	tors_fk = 10.0_prec
	integer, parameter			::	max_cg_iterations = 100, max_line_iterations = 35
	integer						::	cgiter, lineiter
	real(kind=prec)						::	bnd0
	integer						::	nHang, Hang_atom(max_conn)
	integer						::	Hang_code(max_conn)
	integer						::	rule
	type(LIB_ENTRY_TYPE), pointer:: lp
	integer						::	a, b, axis
	real(kind=prec)						::	bond_length, db
	real(kind=prec)						::	rjH(3), rjk(3), bjHinv, bjkinv
	real(kind=prec)						::	scp, angle, angle_deg, dVangle, da, f1
	real(kind=prec)						::	xkt(3), xlt(3), rjkt(3), rktlt(3)
	real(kind=prec)						::	rnj(3), rnk(3), bj, bk
	real(kind=prec)						::	phi, phi_deg, sgn, dVtors, arg, dH(3)
	logical						::	flipped
	integer                     :: setH
	integer, parameter          :: nsetH = 5   !number of times to flip, if local min, and retry

	genH = 0
	xj(:) = xtop(3*j-2:3*j)

	!First find the bond, angle and torsion parameters to use

	!loop over connected atoms
	do ligand = 1, nconn(j)

		H = iconn(ligand, j)
		!is it an unmade H?
		if(makeH(H)) then
			!make this H
			!find the bond potential
			do b = 1, nbonds
				if((bnd(b)%i == j .and. bnd(b)%j == H) .or. &
				   (bnd(b)%j == j .and. bnd(b)%i == H)) then
					if(bnd(b)%cod > 0) then
						bnd0 = bondlib(bnd(b)%cod)%bnd0
					else
						bnd0 = 1.0_prec !default when missing parameters
					endif
					exit
				end if
			end do
			!a unique Hbond_code will always be found
			!because list of connections is based on bond list

			!find all the angles
			nHang = 0
			do a = 1, nangles
				if(ang(a)%j == j .and. ang(a)%cod > 0) then
					if(ang(a)%i == H .and. .not. makeH(ang(a)%k)) then
						nHang = nHang + 1
						Hang_atom(nHang) = ang(a)%k
						Hang_code(nHang) = ang(a)%cod
					elseif(ang(a)%k == H .and. .not. makeH(ang(a)%i)) then
						nHang = nHang + 1
						Hang_atom(nHang) = ang(a)%i
						Hang_code(nHang) = ang(a)%cod
					end if
				end if
			end do
			!look up build rule
			lp => lib(res(residue)%irc)
			do rule = lp%nrules, 1, -1
				if(lp%rules(rule)%kind == BUILD_RULE_TORSION) then
					if(lp%rules(rule)%atom(1) == H - res(residue)%start + 1 .and. &
						lp%rules(rule)%atom(2) == j - res(residue)%start + 1) then
						exit
					end if
				end if
			end do

			!find the torsion rule defined in the library, if any
			if(rule > 0) then
				kt = lp%rules(rule)%atom(3) + res(residue)%start - 1
				lt = lp%rules(rule)%atom(4) + res(residue)%start - 1
				xkt(:) = xtop(3*kt-2:3*kt)
				xlt(:) = xtop(3*lt-2:3*lt)
				rjkt(:) = xkt(:) - xj(:)
				rktlt(:) = xlt(:) - xkt(:)
				rnk(:) = -cross_product(rjkt, rktlt)
				bk = sqrt(dot_product(rnk, rnk))
			end if
			!generate initial co-ordinates for H
			!random vector
			do axis = 1,3
				xH(axis) = randm() - 0.5_prec
			end do
			!normalise to unit length and scale by bond length from lib.
			bond_length = sqrt(dot_product(xH,xH))
			xH(:) = xH(:) / bond_length * bnd0
			!place near atom j
			xH(:) = xH(:) + xj(:)
!			write(*,800) 0,H,xH(:)
			!conjugate gradient minimisation
			do setH = 1,nsetH
			do cgiter = 1, max_cg_iterations
!				write(*,9, advance='no') cgiter
9	format('cg step',i3,':')
				!do line search
				dx_line = max_dx*2.0_prec**(-cgiter)
				flipped = .false.
				do lineiter = 1, max_line_iterations
					Vtot = zero
					dVtot(:) = zero
					!calc. potential & gradient
					!angles
					rjH(:) = xH(:) - xj(:)
					do a = 1, nHang
						xk(:) = xtop(3*Hang_atom(a)-2:3*Hang_atom(a))
						rjk(:) = xk(:) - xj(:)
						bjHinv = one/sqrt(dot_product(rjH, rjH))
						bjkinv = one/sqrt(dot_product(rjk, rjk))
						! calculate scp and angv
						scp = dot_product(rjH, rjk) * bjHinv*bjkinv
						if(scp >  one) then
							scp =  one
						else if(scp < -one) then
							scp = -one
						end if
						angle = acos(scp)
						angle_deg = angle / deg2rad
						! calculate da and dv
						da = angle - anglib(Hang_code(a))%ang0*deg2rad
						V = 0.5_prec*anglib(Hang_code(a))%fk*da**2
						Vtot = Vtot + V
						dVangle = anglib(Hang_code(a))%fk*da
						! calculate f1
						f1 = sin ( angle )
						! avoid division by zero
						if ( abs(f1) < 1.e-12_prec ) then
							f1 = -1.e12_prec
						else
							f1 =  -one / f1
						end if
						dV(:) = dVangle &
							* (f1*(rjk(:)*bjHinv*bjkinv - scp*rjH(:)*bjHinv*bjHinv))
						dVtot(:) = dVtot(:) + dV(:)
					end do

					!torsion
					if(rule > 0) then
						rnj(1) = rjH(2) * rjkt(3) - rjH(3) * rjkt(2)
						rnj(2) = rjH(3) * rjkt(1) - rjH(1) * rjkt(3)
						rnj(3) = rjH(1) * rjkt(2) - rjH(2) * rjkt(1)
						bj = sqrt(dot_product(rnj, rnj))
						scp =dot_product(rnj, rnk)/(bj * bk)
						if(scp>one) scp = one
						if(scp< -one) scp = -one
						phi = acos(scp)
						phi_deg = phi / deg2rad
						sgn = rjkt(1) *(rnj(2) * rnk(3) - rnj(3) * rnk(2) ) &
							+ rjkt(2) *(rnj(3) * rnk(1) - rnj(1) * rnk(3) ) &
							+ rjkt(3) *(rnj(1) * rnk(2) - rnj(2) * rnk(1) )
						if(sgn<0) phi = - phi
						arg = phi-lp%rules(rule)%value*deg2rad
						V = tors_fk*(one+cos(arg))
						!note changed sign of dVtors to get min (not max) at rule value
						dVtors = +tors_fk*sin(arg)

						! ---       forces

						f1 = sin ( phi )
						if ( abs(f1) .lt. 1.e-12_prec ) f1 = 1.e-12_prec
						f1 =  -one/ f1
						dH(:) = f1*(rnk(:)/(bj*bk) - scp*rnj(:)/(bj*bj))
						dV(:) = dVtors*cross_product(rjk, dH)
						dVtot(:) = dVtot(:) + dV(:)
					end if
					if(cgiter == 1 .and. lineiter == 1) then
						!its the start of the search, use the gradient vector
						dvLast(:) = dVtot(:)
					elseif(VtotLast < Vtot) then
						dx_line = -dx_line !next step back
						flipped = .true.
					endif
					rms_dV = sqrt(abs(dot_product(dVtot, dVLast)))
					if(rms_dV < convergence_criterum / 10) then
						exit !reached a potential minimum along the search line
					endif

					!update position in line search
					dx(:) = -dVLast(:) / sqrt(dot_product(dVlast, dVLast)) * dx_line
					if(flipped) dx_line = 0.5_prec * dx_line !reduce step size only after first change of direction
					VtotLast = Vtot
					xH(:) = xH(:) + dx(:)
				end do !lineiter

!				write(*,*) lineiter

				rjH(:) = xH(:) - xj(:)
				bond_length = sqrt(dot_product(rjH, rjH))
				rjH(:) = rjH(:) / bond_length * bnd0      !adjust bond length
				xH(:) = xj(:) + rjH(:)
				rms_dV = sqrt(dot_product(dVtot,dVtot))
				if(rms_dV <  convergence_criterum) then
					if(Vtot < local_min) then
						!found global min
						exit
					else
						!it's a local min - flip 180 deg.
						dVlast(:) = -dVlast(:)
						xH(:) = xj(:) - rjH(:)
					end if
				end if

				!get new gradient search direction
				gamma = dot_product(dVtot, dVlast) / dot_product(dVlast, dVlast)
				!use conjugate gradient
				dVlast(:) = dVtot(:) - gamma*dVlast(:)

			end do  !cgiter

				!Check if local min -> restart iteration ; problem with conversion
					if(Vtot > local_min) then
						!it's a local min - flip 180 deg.
						dVlast(:) = -dVlast(:)
						xH(:) = xj(:) - rjH(:)
					else
					exit
					end if
			end do !setH
			!check if not converged
			if(cgiter >= max_cg_iterations .and. Vtot > 3.0_prec) then
				!display warning
				write(*,900) H, j, Vtot
900				format('>>> WARNING: Positioning of hydrogen',i6,&
					' bound to atom',i6,' didn''t converge.',/, &
					'Potential is ',f8.3,' kcal/mol.')
			end if
			!copy coordinates to topology
			xtop(3*H-2:3*H) = xH(:)
			!clear makeH flag
			makeH(H) = .false.
			genH = genH + 1
		end if
	end do !hydrogens to make
end function genh

integer function genHeavy(residue,missing_heavy,missing_bonds,missing_angles)
!arguments
	integer, intent(in)			::	residue
	type(MISSING_TYPE), intent(in)		::	missing_heavy(:)
	type(MISSING_BOND_TYPE),intent(in)	::	missing_bonds(:)
	type(MISSING_ANGLE_TYPE),intent(in)	::	missing_angles(:)
!locals
	real(kind=prec)						::	xj(3), xk(3)
	integer						::	ligand, H, kt, lt
	real(kind=prec)						::	old_xH(3), xH(3), V, Vtot, dV(3), dvLast(3), gamma, dVtot(3)
	real(kind=prec)						::	VtotLast
	real(kind=prec)						::	dx(3), dx_line, rms_dV
	real(kind=prec),parameter			::	convergence_criterum = 0.1_prec
	real(8),parameter			::	dV_scale = 0.025_prec
	real(8),parameter			::	max_dx = 1.0_prec !max_dx is max distance of line search step in first CG iteration (Å)
	real(kind=prec)						::	local_min = 30
	real(kind=prec)						::	tors_fk = 10.0_prec
	integer, parameter			::	max_cg_iterations = 100, max_line_iterations = 35
	integer						::	cgiter, lineiter
	real(kind=prec)						::	bnd0
	integer						::	nHang, Hang_atom(max_conn)
	integer						::	Hang_code(max_conn)
	integer						::	rule
	type(LIB_ENTRY_TYPE), pointer:: lp
	integer						::	a, b, axis
	real(kind=prec)						::	bond_length, db
	real(kind=prec)						::	rjH(3), rjk(3), bjHinv, bjkinv
	real(kind=prec)						::	scp, angle, angle_deg, dVangle, da, f1
	real(kind=prec)						::	xkt(3), xlt(3), rjkt(3), rktlt(3)
	real(kind=prec)						::	rnj(3), rnk(3), bj, bk
	real(kind=prec)						::	phi, phi_deg, sgn, dVtors, arg, dH(3)
	logical						::	flipped
	integer                     :: setH
	integer, parameter          :: nsetH = 5   !number of times to flip, if local min, and retry

	genHeavy = 0
	xj(:) = xtop(3*j-2:3*j)

! for the solvent all heavy atoms other than the center will be missing initially
! so generate all of them sequentially

	!loop over connected atoms for first atom
	do addatom = 1, lib((irc_solvent)%irc)%nat
		do mbond = 1, 4
		if (missing_heavy(addatom)%bond_notset(mbond) .eqv. .true.) then
			miss = missing_heavy(addatom)%bonds(mbond)
			if (missing_bonds(miss)%i .eq. addatom) then
				H = missing_bonds(miss)%j
			else
				H = missing_bonds(miss)%i
			end if
			if (missing_bonds(miss)%cod .gt. 0) then
				bnd0 = bondlib(missing_bonds(miss)%cod)%bnd0
			else
				bnd0 = 1.0_prec !default when missing parameters
			end if

! now check all the angles for the atoms, if we can add them already
			nHang = 0
			do mangle = 1 , 12
				if(missing_heavy(addatom)%angle_notset(mangle) .eqv. .true.) then
					miss = missing_heavy(addatom)%angles(mangle)
					if(missing_heavy(missing_angles(miss)%i)%atom_missing .or. missing_heavy(missing_angles(miss)%j)%atom_missing &
						.or. missing_heavy(missing_angles(miss)%k)%atom_missing) cycle
					if(
		!find all the angles
		nHang = 0
		do a = 1, nangles
			if(ang(a)%j == j .and. ang(a)%cod > 0) then
				if(ang(a)%i == H .and. .not. makeH(ang(a)%k)) then
					nHang = nHang + 1
					Hang_atom(nHang) = ang(a)%k
					Hang_code(nHang) = ang(a)%cod
				elseif(ang(a)%k == H .and. .not. makeH(ang(a)%i)) then
					nHang = nHang + 1
					Hang_atom(nHang) = ang(a)%i
					Hang_code(nHang) = ang(a)%cod
				end if
			end if
		end do
		!look up build rule
		lp => lib(res(residue)%irc)
		do rule = lp%nrules, 1, -1
			if(lp%rules(rule)%kind == BUILD_RULE_TORSION) then
				if(lp%rules(rule)%atom(1) == H - res(residue)%start + 1 .and. &
					lp%rules(rule)%atom(2) == j - res(residue)%start + 1) then
					exit
				end if
			end if
		end do

		!find the torsion rule defined in the library, if any
		if(rule > 0) then
			kt = lp%rules(rule)%atom(3) + res(residue)%start - 1
			lt = lp%rules(rule)%atom(4) + res(residue)%start - 1
			xkt(:) = xtop(3*kt-2:3*kt)
			xlt(:) = xtop(3*lt-2:3*lt)
			rjkt(:) = xkt(:) - xj(:)
			rktlt(:) = xlt(:) - xkt(:)
			rnk(:) = -cross_product(rjkt, rktlt)
			bk = sqrt(dot_product(rnk, rnk))
		end if
		!generate initial co-ordinates for H
		!random vector
		do axis = 1,3
			xH(axis) = randm() - 0.5_prec
		end do
		!normalise to unit length and scale by bond length from lib.
		bond_length = sqrt(dot_product(xH,xH))
		xH(:) = xH(:) / bond_length * bnd0
		!place near atom j
		xH(:) = xH(:) + xj(:)
!			write(*,800) 0,H,xH(:)
		!conjugate gradient minimisation
		do setH = 1,nsetH
		do cgiter = 1, max_cg_iterations
!				write(*,9, advance='no') cgiter
9	format('cg step',i3,':')
			!do line search
			dx_line = max_dx*2.0_prec**(-cgiter)
			flipped = .false.
			do lineiter = 1, max_line_iterations
				Vtot = zero
				dVtot(:) = zero
				!calc. potential & gradient
				!angles
				rjH(:) = xH(:) - xj(:)
				do a = 1, nHang
					xk(:) = xtop(3*Hang_atom(a)-2:3*Hang_atom(a))
					rjk(:) = xk(:) - xj(:)
					bjHinv = one/sqrt(dot_product(rjH, rjH))
					bjkinv = one/sqrt(dot_product(rjk, rjk))
					! calculate scp and angv
					scp = dot_product(rjH, rjk) * bjHinv*bjkinv
					if(scp >  one) then
						scp =  one
					else if(scp < -one) then
						scp = -one
					end if
					angle = acos(scp)
					angle_deg = angle / deg2rad
					! calculate da and dv
					da = angle - anglib(Hang_code(a))%ang0*deg2rad
					V = 0.5_prec*anglib(Hang_code(a))%fk*da**2
					Vtot = Vtot + V
					dVangle = anglib(Hang_code(a))%fk*da
					! calculate f1
					f1 = sin ( angle )
					! avoid division by zero
					if ( abs(f1) < 1.e-12_prec ) then
						f1 = -1.e12_prec
					else
						f1 =  -one / f1
					end if
					dV(:) = dVangle &
						* (f1*(rjk(:)*bjHinv*bjkinv - scp*rjH(:)*bjHinv*bjHinv))
					dVtot(:) = dVtot(:) + dV(:)
				end do

				!torsion
				if(rule > 0) then
					rnj(1) = rjH(2) * rjkt(3) - rjH(3) * rjkt(2)
					rnj(2) = rjH(3) * rjkt(1) - rjH(1) * rjkt(3)
					rnj(3) = rjH(1) * rjkt(2) - rjH(2) * rjkt(1)
					bj = sqrt(dot_product(rnj, rnj))
					scp =dot_product(rnj, rnk)/(bj * bk)
					if(scp>one) scp = one
					if(scp< -one) scp = -one
					phi = acos(scp)
					phi_deg = phi / deg2rad
					sgn = rjkt(1) *(rnj(2) * rnk(3) - rnj(3) * rnk(2) ) &
						+ rjkt(2) *(rnj(3) * rnk(1) - rnj(1) * rnk(3) ) &
						+ rjkt(3) *(rnj(1) * rnk(2) - rnj(2) * rnk(1) )
					if(sgn<0) phi = - phi
					arg = phi-lp%rules(rule)%value*deg2rad
					V = tors_fk*(one+cos(arg))
					!note changed sign of dVtors to get min (not max) at rule value
					dVtors = +tors_fk*sin(arg)

					! ---       forces

					f1 = sin ( phi )
					if ( abs(f1) .lt. 1.e-12_prec ) f1 = 1.e-12_prec
					f1 =  -one/ f1
					dH(:) = f1*(rnk(:)/(bj*bk) - scp*rnj(:)/(bj*bj))
					dV(:) = dVtors*cross_product(rjk, dH)
					dVtot(:) = dVtot(:) + dV(:)
				end if
				if(cgiter == 1 .and. lineiter == 1) then
					!its the start of the search, use the gradient vector
					dvLast(:) = dVtot(:)
				elseif(VtotLast < Vtot) then
					dx_line = -dx_line !next step back
					flipped = .true.
				endif
				rms_dV = sqrt(abs(dot_product(dVtot, dVLast)))
				if(rms_dV < convergence_criterum / 10) then
					exit !reached a potential minimum along the search line
				endif

				!update position in line search
				dx(:) = -dVLast(:) / sqrt(dot_product(dVlast, dVLast)) * dx_line
				if(flipped) dx_line = 0.5_prec * dx_line !reduce step size only after first change of direction
				VtotLast = Vtot
				xH(:) = xH(:) + dx(:)
			end do !lineiter

!				write(*,*) lineiter

			rjH(:) = xH(:) - xj(:)
			bond_length = sqrt(dot_product(rjH, rjH))
			rjH(:) = rjH(:) / bond_length * bnd0      !adjust bond length
			xH(:) = xj(:) + rjH(:)
			rms_dV = sqrt(dot_product(dVtot,dVtot))
			if(rms_dV <  convergence_criterum) then
				if(Vtot < local_min) then
					!found global min
					exit
				else
					!it's a local min - flip 180 deg.
					dVlast(:) = -dVlast(:)
					xH(:) = xj(:) - rjH(:)
				end if
			end if

			!get new gradient search direction
			gamma = dot_product(dVtot, dVlast) / dot_product(dVlast, dVlast)
			!use conjugate gradient
			dVlast(:) = dVtot(:) - gamma*dVlast(:)

		end do  !cgiter

			!Check if local min -> restart iteration ; problem with conversion
				if(Vtot > local_min) then
					!it's a local min - flip 180 deg.
					dVlast(:) = -dVlast(:)
					xH(:) = xj(:) - rjH(:)
				else
				exit
				end if
		end do !setH
		!check if not converged
		if(cgiter >= max_cg_iterations .and. Vtot > 3.0_prec) then
			!display warning
			write(*,900) H, j, Vtot
900				format('>>> WARNING: Positioning of heavy atom',i6,&
				' bound to atom',i6,' didn''t converge.',/, &
				'Potential is ',f8.3,' kcal/mol.')
		end if
		!copy coordinates to topology
		xtop(3*H-2:3*H) = xH(:)
		!clear makeH flag
		genHeavy = genHeavy + 1
	end do !heavy atoms to make
end function genheavy

!-----------------------------------------------------------------------

subroutine makeimps
! *** local variables
	integer						::	i

!.......................................................................

	nimps = 0
	nimps_solute = 0

	do i = 1, nat_pro
		if(nconn(i) ==3) then
			nimps = nimps + 1
			if(i <= nat_solute) nimps_solute = nimps_solute + 1
			imp(nimps)%i = iconn(1, i)
			imp(nimps)%j = i
			imp(nimps)%k = iconn(2, i)
			imp(nimps)%l = iconn(3,i)
		end if
	enddo

end subroutine makeimps

!-----------------------------------------------------------------------

subroutine imp_params
! *** local variables
	integer i, ii, iaci, iacj, iack, iacl
	integer iused(nimp_prm+max_extrabnd), itrans(nimp_prm+max_extrabnd)


	do i = 1, nimps
		iaci = iac(imp(i)%i )
		iacj = iac(imp(i)%j )
		iack = iac(imp(i)%k )
		iacl = iac(imp(i)%l )

		imp(i)%cod = impcode(tac(iaci), tac(iacj), tac(iack), tac(iacl))

		if(imp(i)%cod == 0) then
			write(*,20) imp(i)%i, imp(i)%j, imp(i)%k, imp(i)%l
		end if
	enddo
20	format('Atoms of missing improper:',4i5)

! --- Make a list of actual angle types to be used

	iused(:) = 0

	nimpcod = 0

	do i = 1, nimps
!	workaround to permit missing param's
		if(imp(i)%cod /=0) then
			if(iused(imp(i)%cod ) ==0) then
				iused(imp(i)%cod ) = 1
				nimpcod = nimpcod+1
				itrans(imp(i)%cod ) = nimpcod
				implib(nimpcod) = imp_prm(imp(i)%cod)%prm
				imp(i)%cod = nimpcod
			ELSE
				imp(i)%cod = itrans(imp(i)%cod )
			endif
		endif
	enddo

	write(*, 100) nimps_solute, 'solute'
	write(*, 100) nimps-nimps_solute, 'solvent'
100	format('Made    ',i6,1x,a,' impropers.')
!.......................................................................
end subroutine imp_params

!-----------------------------------------------------------------------

integer function find_atom(ires, atom)
!arguments
	integer						::	ires
	character(len=5)			::	atom
!locals
	integer						::	my_res, iat, irc
	character(len=4)			::	my_atom
	integer						::	mol1, mol2

	if(atom(1:1) == '-') then
		my_res = ires - 1
		my_atom = atom(2:5)
	elseif(atom(1:1) == '+') then
		my_res = ires + 1
		my_atom = atom(2:5)
	else
		my_res = ires
		my_atom = atom(1:4)
	end if

	find_atom = 0

	!first may not refer backwards and last may not refer ahead
	if(my_res < 1 .or. my_res > nres) then
		!can't move outside terminal residues - skip
		return
	end if

	!check if moving outside molecule
	if(ires /= my_res) then
		do mol1 = 1, nmol
			if(res(ires)%start < istart_mol(mol1+1)) exit
		end do
		do mol2 = 1, nmol
			if(res(my_res)%start < istart_mol(mol2+1)) exit
		end do
		if(mol1 /= mol2) return !crossed molecule boundary
	end if

	irc = res(my_res)%irc !get residue code

	do iat = 1, lib(irc)%nat
		if(lib(irc)%atnam(iat) == my_atom) then
			find_atom = res(my_res)%start + iat - 1
			exit
		end if
	end do

	if(ires == my_res .and. find_atom == 0) then
		!we have an error: an atom in _this_ residue was not found
		topo_ok = .false.
		write(*, 120) my_atom, lib(irc)%nam, ires
	end if

120	format('>>>>> ERROR: There is no atom named ',a4,' in residue ',a4,i5,'.')

end function find_atom

!-----------------------------------------------------------------------

subroutine makeimps_explicit
!generate impropers using definitions in library file
! *** local variables
	integer						::	ires, irc, iimp, i

	nimps = 0
	nimps_solute = 0

	do ires = 1, nres !loop over residues
		irc = res(ires)%irc
		do iimp = 1, lib(irc)%nimp !loop over improper defs
			nimps = nimps+1
			if(ires <= nres_solute) nimps_solute = nimps_solute + 1
			imp(nimps)%i = find_atom(ires, lib(irc)%imp(iimp)%i)
			imp(nimps)%j = find_atom(ires, lib(irc)%imp(iimp)%j)
			imp(nimps)%k = find_atom(ires, lib(irc)%imp(iimp)%k)
			imp(nimps)%l = find_atom(ires, lib(irc)%imp(iimp)%l)
			if(imp(nimps)%i == 0 .or. imp(nimps)%j == 0 .or. &
				imp(nimps)%k == 0 .or. imp(nimps)%l == 0) then
				!not found
				nimps = nimps - 1 !take a step back, discard this one
				if(ires <= nres_solute) nimps_solute = nimps_solute - 1
			end if
		end do
	end do

end subroutine makeimps_explicit

!-----------------------------------------------------------------------

subroutine maketop

	!check if library is loaded
	if(.not. check_residues()) return

!	bail out if readparm failed
	if(.not.ff_ok) then
		write(*,900)
900		format('>>>>> ERROR: Force field parameters not loaded.')
		return
	end if
	topo_ok = .false.

	!check if boundaries are set

    if(.not. boundary_set) then
		write(*,910)
910		format('>>>>> ERROR: Boundary has not been set.')
		return
	end if

	if( .not. have_solvent_boundary) then
        xwcent(:) = xpcent(:)
		write(*,911)
		write(*,912)
911		format('>>>>> WARNING: System has not been solvated.')
912		format('Using boundary center as solvent center.')
	end if

	if(.not. have_title) then
		CALL get_line_arg(title, '-----> Give this topology a title: ')
		have_title = .true.
	end if
	topo_ok = .true. !now hope for the best

	!deallocate old topology if any
	if(allocated(iac)) deallocate(iac)
	if(allocated(crg)) deallocate(crg)
	if(allocated(cgpatom)) deallocate(cgpatom)
	if(allocated(list14)) deallocate(list14)
	if(allocated(listex)) deallocate(listex)
	if(allocated(nconn)) deallocate(nconn)
	if(allocated(iconn)) deallocate(iconn)
	if(allocated(bnd)) deallocate(bnd)
	if(allocated(bondlib)) deallocate(bondlib)
	if(allocated(SYBYL_bond_type)) deallocate(SYBYL_bond_type)
	if(allocated(ang)) deallocate(ang)
	if(allocated(anglib)) deallocate(anglib)
	if(allocated(tor)) deallocate(tor)
	if(allocated(torlib)) deallocate(torlib)
	if(allocated(imp)) deallocate(imp)
	if(allocated(implib)) deallocate(implib)
	if(allocated(list14long)) deallocate(list14long)
	if(associated(listexlong)) deallocate(listexlong)
	if(allocated(cgp)) deallocate(cgp)
	if(allocated(nconn)) deallocate(nconn)
	if(allocated(iconn)) deallocate(iconn)
	if(allocated(excl)) deallocate(excl)
	if(allocated(shell)) deallocate(shell)

	!allocate topology arrays (set max_atom etc
	call topo_set_max(nat_pro, max_lib, max_long) !make space for new topology
	allocate(iac(max_atom), &
		crg(max_atom), &
		cgpatom(max_atom), &
		list14(max_nbr_range, max_atom), &
		listex(max_nbr_range, max_atom), &
		nconn(max_atom), iconn(max_conn, max_atom), &
		stat=alloc_status)
	if(alloc_status /= 0) then
		write(*,*) 'ERROR: Out of memory when allocating topology arrays'
		stop 255
	end if

	if( .not. use_PBC ) then
		allocate(excl(max_atom), shell(max_atom), stat=alloc_status)
		if(alloc_status /= 0) then
			write(*,*) 'ERROR: Out of memory when allocating topology arrays'
			stop 255
		end if
	end if

	allocate(bnd(max_bonds), &
		bondlib(max_bondlib), &
		SYBYL_bond_type(max_bondlib), &
		ang(max_angles), &
		anglib(max_anglib), &
		tor(max_tors), &
		torlib(max_torlib), &
		imp(max_imps), &
		implib(max_implib), &
		stat=alloc_status)
	if(alloc_status /= 0) then
		write(*,*) 'ERROR: Out of memory when allocating topology arrays'
		stop 255
	end if

	allocate(list14long(2, max_14long), &
		listexlong(2, max_exlong), &
		cgp(max_cgp), &
		stat=alloc_status)
	if(alloc_status /= 0) then
		write(*,*) 'ERROR: Out of memory when allocating topology arrays'
		stop 255
	end if

	!set date of creation
	call date_and_time(creation_date)

	!set atom types, charges, charge groups
	CALL set_iac
	if(.not. topo_ok) return

	call set_solvent_type
	CALL set_crg
	CALL makebonds
	call set_bondcodes
	call makeconn
	CALL makeangles
	CALL maketors
	if(imp_explicit) then
		call makeimps_explicit !explicit (from lib.) generation
	else
		CALL makeimps !automatic generation
	end if
	call imp_params !set codes & filter out unused types
	CALL make14list
	CALL makeexlist
	CALL makehyds   !generate hydrogens
	CALL set_cgp    !assigns all cgp() within exclusion radius

	call set_default_mask

!	Now check if we have missing parameters
	if(topo_ok) then
		write(*,'(a)') 'Topology successfully generated.'
	else
		write(*, '(/,a,/)') 'ERROR: There are missing parameters!'
		write(*, '(a)') 'You need to add parameters and do readprm + maketop again.'
		write(*, 100) trim(prm_file)
	endif
  100 format('parameter file: ', a)
end subroutine maketop

!----------------------------------------------------------------------

subroutine set_default_mask

	!make atom mask including all atoms
	!first get rid of old mask
	call mask_finalize(mask)
	!allocate new mask
	call mask_initialize(mask)
	!add all atoms (in_mask is not used)
	call mask_all(mask)
end subroutine set_default_mask
!-----------------------------------------------------------------------

subroutine makeconn
!locals
	integer						::	i, ib

	nconn(:) = 0
! --- Make a list of connections to each atom

	do i = 1, nat_pro
		nconn(i) = 0
		do ib = 1, nbonds
			if(bnd(ib)%i ==i) then
				nconn(i) = nconn(i) + 1
				iconn(nconn(i), i) = bnd(ib)%j
			ELSEif(bnd(ib)%j ==i) then
				nconn(i) = nconn(i) + 1
				iconn(nconn(i), i) = bnd(ib)%i
			endif
		enddo
	enddo
end subroutine makeconn

!-----------------------------------------------------------------------

subroutine maketors
! *** local variables
	integer i, j, ic, jc, ib, it, iaci, iacj, iack, iacl
	integer iused(ntor_prm+max_extrabnd), itrans(ntor_prm+max_extrabnd)
	type(TOR_CODES)				::	torcodes
	integer						::	icod
! --- Make the torsion list

	ntors = 0
	ntors_solute = 0

	do ib = 1, nbonds

	do ic = 1, nconn(bnd(ib)%i )
	i = iconn(ic, bnd(ib)%i )
	if(i/=bnd(ib)%j ) then
		do jc = 1, nconn(bnd(ib)%j )
		j = iconn(jc, bnd(ib)%j )
		if(j/=bnd(ib)%i ) then

			if(i/=j) then

				iaci = iac(i)
				iacj = iac(bnd(ib)%i)
				iack = iac(bnd(ib)%j)
				iacl = iac(j)

				torcodes = torcode(tac(iaci), tac(iacj), tac(iack), tac(iacl))
				if(torcodes%ncod == 0) then
					write( * , * ) 'atoms: ', i , bnd(ib)%i , bnd(ib)%j, j
				endif

				do icod = 1, torcodes%ncod
					ntors = ntors + 1
					if(i <= nat_solute) ntors_solute = ntors_solute + 1
					tor(ntors)%i = i
					tor(ntors)%j = bnd(ib)%i
					tor(ntors)%k = bnd(ib)%j
					tor(ntors)%l = j
					tor(ntors)%cod = torcodes%cod(icod)
				end do
!				tor(ntors)%cod = torcode(tac(iaci), tac(iacj), tac(iack), tac(iacl))
!						   workaround to allow missing param's
!				if(tor(ntors)%cod /=0) then
!					if(tor_prm(tor(ntors)%cod)%rmult <= - 1.0) then
!						more_terms = .true.
!						do while(more_terms)
!							ntors = ntors + 1
!							if(i <= nat_solute) ntors_solute = ntors_solute + 1
!							tor(ntors)%i = tor(ntors - 1)%i
!							tor(ntors)%j = tor(ntors - 1)%j
!							tor(ntors)%k = tor(ntors - 1)%k
!							tor(ntors)%l = tor(ntors - 1)%l
!							tor(ntors)%cod = tor(ntors - 1)%cod + 1
!							more_terms = .false.
!							if(tor_prm(tor(ntors)%cod)%rmult <= -1.0)  &
!								more_terms = .true.
!						enddo
!					endif
!				endif
			endif

		endif
		enddo
	endif
	enddo

	enddo

! --- Make a list of actual torsion types to be used

	iused(1:ntor_prm) = 0

	ntorcod = 0

	do i = 1, ntors
!	   workaround to permit missing param's
		if(tor(i)%cod /=0) then
			if(iused(tor(i)%cod ) ==0) then
				iused(tor(i)%cod ) = 1
				ntorcod = ntorcod+1
				itrans(tor(i)%cod ) = ntorcod
				torlib(ntorcod) = tor_prm(tor(i)%cod)
				tor(i)%cod = ntorcod
			ELSE
				tor(i)%cod = itrans(tor(i)%cod)
			endif
		endif
	enddo

	write(*, 100) ntors_solute, 'solute'
	write(*, 100) ntors-ntors_solute, 'solvent'
100	format('Made    ',i6,1x,a,' torsions.')
!.......................................................................
end subroutine maketors

!-----------------------------------------------------------------------

subroutine prompt(outtxt)
! *** local variables
	CHARACTER( * ) outtxt

!.......................................................................

	write( * , '(a,$)') outtxt

	return
!.......................................................................
end subroutine prompt

!-----------------------------------------------------------------------

real FUNCTION randm(seed, seed_only)
!arguments
	integer, intent(in), optional::	seed
	logical, intent(in), optional::	seed_only
! *** Local variables
	integer, parameter			::	m = 100000000, m1 = 10000, mult = 31415821
	integer, save				::	irand = 0
	integer						::	irandh, irandl, multh, multl
	real(kind=prec)						::	r

	if(present(seed)) then
		irand = mod(iabs(seed), m)
	end if
	if(present(seed_only)) return
!
! --- multiply irand by mult, but take into account that overflow must
! --- be discarded, and do not generate an error.
!
	irandh = irand / m1
	irandl = mod(irand, m1)
	multh = mult / m1
	multl = mod(mult, m1)
!
	irand = mod(irandh * multl + irandl * multh, m1) * m1 + irandl * &
	multl
	irand = mod(irand+1, m)
!
! --- convert irand to a real random number between 0 and 1.
!
	r = real(irand / 10_prec,kind=prec) * 10_prec / real(m,kind=prec)
	if((r<=0.e0_prec) .or.(r>1.e0_prec) ) r = 0.e0_prec
	randm = r

!.......................................................................
end FUNCTION randm

!-----------------------------------------------------------------------

subroutine oldreadlib(filnam)
!arguments
	character(*)				::	filnam
! *** local variables
	CHARACTER					::	line*80
	integer						::	irec, i, iat, ires, j, igp, ntot
	real(kind=prec)						::	qtot, qgrp, qtot_grp

	!some extra arrays used only for allocation. This circumvents
	!an unaligned access error on Digitail UNIX 4.0 / Digital FORTRAN
	!(Compiler bug?)

	character(len=KEYLENGTH), pointer::	tac_lib(:)
	character(len=4), pointer	::	atnam(:)
	real(kind=prec), pointer				::	crg_lib(:)
	integer(AI), pointer		::	natcgp(:), switch(:)
	integer(AI), pointer		::	atcgp(:,:)
	type(LIB_BOND_TYPE), pointer::bnd(:)
	integer						:: stat

1	format(a,t7,a,t14,a,t21,a,t28,a)
2	format(a4)
3	format(t9,i3)
4	format(t14,f5.2)
5	format(t23,i3)
6	format(t31,i3)
!.......................................................................

	if(openit(1, filnam, 'old', 'formatted', 'read') /= 0) then
		write(*, '(a,a)') '>>>>>ERROR: Could not open library file ', &
			trim(filnam)
		return
	end if

	write( * , '(/,a,a,/)') 'Reading molecular library ',trim(filnam)
	write(*,1) 'name', 'atoms','net Q', 'bonds', 'Q-grps'
	do irec = 1, 99999
		call skip_comments(1)
		READ(1, '(a80)') line
		if(line=='end') exit
		BACKSPACE(1)
		nlibres = nlibres + 1
		if(nlibres > max_entry) then
			write(*,*) 'ERROR: Too many library entries!'
			write(*,*) '(This problem may be solved by increasing max_entry.)'
			exit
		end if
		ires = nlibres

! ---	   Read entry name

		lib(nlibres)%SYBYLTYPE = '****'
		lib(nlibres)%HETATM =.false. !don't use new feature
		lib(nlibres)%solvent =.false. !don't use new feature
		lib(nlibres)%density =0 !don't use new feature
		read(1, '(a)') line
		READ(line, fmt=*, err=9, end=9) lib(nlibres)%nam, lib(nlibres)%SYBYLTYPE
9		write( *, 2, advance='no') lib(nlibres)%nam
! ---	   Read no. of atoms, at. no., name, iac, charge

		READ(1, * , iostat=stat) lib(nlibres)%nat
		if (stat /= 0) then
			write(*,*)
			write(*,*) 'Incorrect input!'
			return
		endif
		write(*,3, advance='no') lib(nlibres)%nat

		allocate(atnam(lib(nlibres)%nat))
		lib(nlibres)%atnam => atnam

		allocate(tac_lib(lib(nlibres)%nat))
		lib(nlibres)%tac_lib => tac_lib

		allocate(crg_lib(lib(nlibres)%nat))
		lib(nlibres)%crg_lib => crg_lib

		qtot = 0
		do i = 1, lib(nlibres)%nat
			READ(1, *) iat, lib(nlibres)%atnam(iat), lib(nlibres)%tac_lib(iat), &
				lib(nlibres)%crg_lib(iat)
			qtot = qtot + lib(nlibres)%crg_lib(iat)
		enddo
		write(*, 4, advance='no') qtot
! ---	   Read no. of bonds, bond list, tail/head connections

		READ(1, * ) lib(nlibres)%nbnd
		allocate(bnd(lib(nlibres)%nbnd), stat=alloc_status)
		lib(nlibres)%bnd => bnd
		write(*, 5, advance='no') lib(nlibres)%nbnd
		do i = 1, lib(nlibres)%nbnd
			READ(1, * ) lib(nlibres)%bnd(i)%i, lib(nlibres)%bnd(i)%j
		enddo

		READ(1, * ) lib(nlibres)%head, lib(nlibres)%tail

! ---	   Read charge group info

		READ(1, * ) lib(nlibres)%ncgp
		write(*, 6) lib(nlibres)%ncgp
		allocate(atcgp(max(max_atcgplib, lib(nlibres)%nat), lib(nlibres)%ncgp))
		lib(nlibres)%atcgp => atcgp
		allocate(natcgp(lib(nlibres)%ncgp), switch(lib(nlibres)%ncgp))

		lib(nlibres)%natcgp => natcgp
		lib(nlibres)%switch => switch

		qtot_grp = zero
		ntot = 0
		do i = 1, lib(nlibres)%ncgp
			READ(1, * ) lib(nlibres)%natcgp(i), lib(nlibres)%switch(i)
			ntot = ntot + lib(nlibres)%natcgp(i)
			READ(1, * )lib(nlibres)%atcgp(1:lib(nlibres)%natcgp(i), i)
			qgrp = 0
			do j = 1, lib(nlibres)%natcgp(i)
				!The Cray compiler doesn't like lib(nlibres)%crg_lib(lib(nlibres)%atcgp(j, i))
				igp = lib(nlibres)%atcgp(j, i)
				qgrp = qgrp + lib(nlibres)%crg_lib(igp)
			enddo
			!check fractional charges
			if(abs(qgrp - nint(qgrp) ) >0.000001_prec)  then
				write( *, 130) qgrp, i
			end if
			qtot_grp = qtot_grp + qgrp
		enddo
		!check total charge & atom count consistency
		if(ntot /= lib(nlibres)%nat) then
			write(*, 150)
		endif
		if(abs(qtot - qtot_grp)  >0.000001_prec) then
			write(*,160)
		end if
	enddo
	close(1)

	write( *, 110) nlibres

110 format(/,'Accumulated no. of library entries loaded =',i4,/)

130 format('>>> Warning: fractional charge = ',f8.5, ' on group no.',i3)
150	format('>>> ERROR: charge group atom count does not match total atom count!')
160	format('>>> ERROR: Sum of charge group charges not equal to sum of all charges!')

end subroutine oldreadlib

!-----------------------------------------------------------------------

subroutine readlib(file)
!arguments
	character(*), optional		::	file
! *** local variables
	CHARACTER					::	line*80, filnam*80
	integer						::	irec, i, iat, ires, j, igp, ntot
	real(kind=prec)						::	qtot, qgrp, qtot_grp
	character(len=80)			::	resnam
	integer						::	res_count
	logical						::	prm_res
	integer						::	cgp_read(max_atcgplib)

	!some extra arrays used only for allocation. This circumvents
	!an unaligned access error on Digital UNIX 4.0 / Digital FORTRAN
	!(Compiler bug?)

	character(len=KEYLENGTH), pointer::	tac_lib(:)
	character(len=4), pointer	::	atnam(:)
	real(kind=prec), pointer				::	crg_lib(:)
	integer(AI), pointer		::	natcgp(:), switch(:)
	integer(AI), pointer		::	atcgp(:,:)
	type(LIB_BOND_TYPE), pointer::	bnd(:)
	type(LIB_IMP_TYPE), pointer	::	imp(:)
	logical						::	lib_exists
	character(len=8)			::	atnam1, atnam2, atnam3, atnam4
	integer						::	atno1, atno2, atno3, atno4
	logical						::	yes
	integer						::	atoms_in_rule, readstat

1	format(10a7)
2	format(a4)
3	format(i7)
4	format(f7.2)
7	format(4x,a)
!.......................................................................

	res_count = 0
	if(present(file)) then
		filnam = file
	else
		call get_string_arg(filnam, '-----> Name of molecular library: ')
	end if
	inquire(file=filnam, exist=lib_exists)
	if(.not. lib_exists) then
		write(*,10) trim(filnam)
		return
	end if

	if(lib_files == '') then
		lib_files = filnam
	elseif(len_trim(lib_files) < 200-len_trim(filnam)) then
		lib_files = trim(lib_files) // ";" // trim(filnam)
	end if

	if(.not. prm_open(filnam)) then
		!if reading fails try old format
		write(*,'(a)') '>>>>>WARNING: Attempting to read old-style library file.'
		call oldreadlib(filnam)
		return
	end if

10	format('>>>>> ERROR: ',a,' does not exist.')

	write( * , '(/,a,a,/)') 'Reading molecular library ',trim(filnam)
	write(*,2, advance='no') 'name'
	write(*,1) 'atoms','net Q', 'bonds', 'rules', 'imps', 'Q-grps'
	!loop over all titles in file
	do while(prm_get_next_title(resnam))
		call check_overload(resnam)
		nlibres = nlibres + 1
		res_count = res_count + 1
		lib(nlibres)%nam = resnam(1:len(lib(nlibres)%nam)) !copy name
		write( *, 2, advance='no') lib(nlibres)%nam
		lib(nlibres)%SYBYLTYPE = '****'
		lib(nlibres)%HETATM = .false.
		lib(nlibres)%solvent = .false.
		if(prm_open_section('info')) then !read info
			if(prm_get_string_by_key('SYBYLTYPE', line)) then
				lib(nlibres)%SYBYLTYPE = line(1:len(lib(nlibres)%SYBYLTYPE))
			end if
			yes=prm_get_string_by_key('PDBTYPE', line, 'ATOM')
			call upcase(line)
			if(trim(line) == 'HETATM') then
				lib(nlibres)%HETATM = .true.
			elseif(trim(line) /= 'ATOM') then
				write(*,'(a,a)') '>>>>>WARNING: Unknown PBBtype ',trim(line)
			end if
			yes=prm_get_logical_by_key('solvent',lib(nlibres)%solvent,.false.)
			if(lib(nlibres)%solvent) then
				solvent_names = trim(solvent_names)//lib(nlibres)%nam//','
			end if
			yes=prm_get_real_by_key('density',lib(nlibres)%density,0.0_prec)
		end if

! ---	   Read no. of atoms, at. no., name, iac, charge
		lib(nlibres)%nat = prm_count('atoms')
		if(lib(nlibres)%nat == 0) then
			write(*,'(a,a)') '>>>>>ERROR: No atoms in residue ',trim(resnam)
			!undo this residue
			nlibres = nlibres - 1
			res_count = res_count - 1
			exit
		end if

		write(*,3, advance='no') lib(nlibres)%nat

		allocate(atnam(lib(nlibres)%nat))

		!initialise atom name -> number index
		call index_create(lib(nlibres)%nat)

		lib(nlibres)%atnam => atnam

		allocate(tac_lib(lib(nlibres)%nat))
		lib(nlibres)%tac_lib => tac_lib

		allocate(crg_lib(lib(nlibres)%nat))
		lib(nlibres)%crg_lib => crg_lib

		qtot = zero
		do i = 1, lib(nlibres)%nat
			prm_res = prm_get_line(line)
			READ(line, *) iat
			!Check numbering of atoms
			if (iat>lib(nlibres)%nat) then
			  write(*,119) iat
		    end if
			READ(line, *) iat, lib(nlibres)%atnam(iat), lib(nlibres)%tac_lib(iat), &
				lib(nlibres)%crg_lib(iat)
			if(.not. index_add(lib(nlibres)%atnam(iat), iat)) then
				!could not add - name duplication?
				write(*, 120) iat, lib(nlibres)%atnam(iat)
			end if
			!add also the number as an atom name
			write(atnam1,'(i4)') iat
			atnam1 = adjustl(atnam1)
			if(.not. index_add(atnam1, iat)) then
				!could not add - number duplication?
				write(*, 121) iat
			end if

			qtot = qtot + lib(nlibres)%crg_lib(iat)
		enddo
		write(*, 4, advance='no') qtot

! if this is a solvent molecule, set solv_atoms variable
                if(lib(nlibres)%solvent) then
                        solv_atom = lib(nlibres)%nat
                end if

! ---	   Read bond list
		lib(nlibres)%nbnd = prm_count('bonds')
		allocate(bnd(lib(nlibres)%nbnd), stat=alloc_status)
		lib(nlibres)%bnd => bnd

		write(*, 3, advance='no') lib(nlibres)%nbnd
		do i = 1, lib(nlibres)%nbnd
			prm_res = prm_get_line(line)
			read(line, *) atnam1, atnam2
			if(.not. index_get(atnam1, atno1)) then
				write(*,125) atnam1, i
				cycle
			end if
			if(.not. index_get(atnam2, atno2)) then
				write(*,126) atnam2, i
				cycle
			end if
			lib(nlibres)%bnd(i)%i = atno1
			lib(nlibres)%bnd(i)%j = atno2
		enddo

		!read head and tail connections, if present
		lib(nlibres)%head = 0
		lib(nlibres)%tail = 0
		if(prm_open_section('connections')) then
			!0 is returned if not found
			atno1 = 0
			if(prm_get_string_by_key('head', atnam1)) then
				if(atnam1 /= '0' .and. .not. index_get(atnam1, atno1)) then
					write(*,127) atnam1
				end if
			end if
			lib(nlibres)%head = atno1
			atno1 = 0
			if(prm_get_string_by_key('tail', atnam1)) then
				if(atnam1 /= '0' .and. .not. index_get(atnam1, atno1)) then
					write(*,128) atnam1
				end if
			end if
			lib(nlibres)%tail = atno1
		end if

! ---	   Read hydrogen build rules
		lib(nlibres)%nrules = prm_count('build_rules')
		allocate(lib(nlibres)%rules(lib(nlibres)%nrules))

		write(*, 3, advance='no') lib(nlibres)%nrules
ruleloop: do i = 1, lib(nlibres)%nrules
			!read kind
			prm_res = prm_get_field(line)
			call locase(line)
			if(line == 'torsion') then
				lib(nlibres)%rules(i)%kind = BUILD_RULE_TORSION
				atoms_in_rule = 4
			else
				lib(nlibres)%rules(i)%kind = 0
				write(*, 200) i, trim(line)
200				format('>>>>> ERROR: Unrecognised rule',i2,' kind ',a)
				!skip to end of line
				yes = prm_get_field(line, skip=.true.)
				!continue to next rule
				cycle ruleloop
			end if
			do iat = 1, atoms_in_rule
				if(.not. prm_get_field(line)) then
					lib(nlibres)%rules(i)%kind = 0
					write(*, 210) i
210					format('>>>>> ERROR: Read error at rule',i2)
					yes = prm_get_field(line, skip=.true.)
					!continue to next rule
					cycle ruleloop
				end if
				if(.not. index_get(line, lib(nlibres)%rules(i)%atom(iat))) then
					lib(nlibres)%rules(i)%kind = 0
					write(*, 220) trim(line), i
220					format('>>>>> ERROR: Atom ',a,' not found in rule',i2)
					yes = prm_get_field(line, skip=.true.)
					!continue to next rule
					cycle ruleloop
				end if
			end do
			if(.not. prm_get_field(line)) then
				lib(nlibres)%rules(i)%kind = 0
				write(*, 230) i
230				format('>>>>> ERROR: No value for rule ',i2)
				yes = prm_get_field(line, skip=.true.)
				cycle ruleloop
			end if
			read(line, *, iostat=readstat) lib(nlibres)%rules(i)%value
			if(readstat /= 0) then
				lib(nlibres)%rules(i)%kind = 0
				write(*, 240) trim(line), i
240				format('>>>>> ERROR: Read error at rule ',i2)
			end if
			yes = prm_get_field(line, skip=.true.)
		end do ruleloop


! ---	   Read improper list
		lib(nlibres)%nimp = prm_count('impropers')
		allocate(imp(lib(nlibres)%nimp), stat=alloc_status)
		lib(nlibres)%imp => imp

		write(*, 3, advance='no') lib(nlibres)%nimp
		do i = 1, lib(nlibres)%nimp
			prm_res = prm_get_line(line)
			!note the order j,i,l,k used to conform with amber and charmm
			!parameter file conventions
			READ(line, * ) lib(nlibres)%imp(i)%i, lib(nlibres)%imp(i)%j, &
				lib(nlibres)%imp(i)%k, lib(nlibres)%imp(i)%l
		enddo

! ---	   Read charge group info
		lib(nlibres)%ncgp = prm_count('charge_groups')
		!if no chargegroups construct one
		if(lib(nlibres)%ncgp == 0) then
			write(*,7) 'none, creating default'
			lib(nlibres)%ncgp = 1
			allocate(atcgp(max(max_atcgplib, lib(nlibres)%nat), 1))
			lib(nlibres)%atcgp => atcgp
			allocate(natcgp(1), switch(1))

			lib(nlibres)%natcgp => natcgp
			lib(nlibres)%switch => switch

			!create a charge group with all atoms
			lib(nlibres)%switch(1) = 1
			lib(nlibres)%natcgp(1) = lib(nlibres)%nat
			do i = 1, lib(nlibres)%nat
				lib(nlibres)%atcgp(i,1) = i
			end do
			!pass consistencty test below
			qtot_grp = qtot
			ntot = lib(nlibres)%nat
		else
			write(*, 3) lib(nlibres)%ncgp
			allocate(atcgp(max(max_atcgplib, lib(nlibres)%nat), lib(nlibres)%ncgp))
			lib(nlibres)%atcgp => atcgp
			allocate(natcgp(lib(nlibres)%ncgp), switch(lib(nlibres)%ncgp))

			lib(nlibres)%natcgp => natcgp
			lib(nlibres)%switch => switch

			qtot_grp = zero
			ntot = 0
			do i = 1, lib(nlibres)%ncgp
				qgrp = 0.
				lib(nlibres)%natcgp(i) = 0
				do while(prm_get_field(line)) !get one number at a time
					if(.not. index_get(line, igp)) then
						write(*,131) lib(nlibres)%natcgp(i) + 1, trim(line), i
						cycle
					end if
					!read(line, *) igp
					qgrp = qgrp + lib(nlibres)%crg_lib(igp)
					lib(nlibres)%natcgp(i) = lib(nlibres)%natcgp(i) + 1
					lib(nlibres)%atcgp(lib(nlibres)%natcgp(i),i) = igp
				end do
				lib(nlibres)%switch(i) = lib(nlibres)%atcgp(1,i)
				ntot = ntot + lib(nlibres)%natcgp(i)

				!check fractional charges
				if(abs(qgrp - nint(qgrp) ) >0.000001_prec)  then
					write( *, 130) qgrp, i
				end if
				qtot_grp = qtot_grp + qgrp
			enddo !charge groups
		end if !charge groups defined

		!check total charge & atom count consistency
		if(ntot /= lib(nlibres)%nat) then
			write(*, 150)
		endif
		if(abs(qtot - qtot_grp)  >0.000001_prec) then
			write(*,160)
		end if
	end do !entries

	call prm_close


	write( *, 110) nlibres

110 FORMAT(/,'Accumulated no. of library entries loaded =',i4,/)
119	format('>>> ERROR: Atom',i3,' has erraneos atom number.')
120	format('>>> ERROR: Atom',i3,' has duplicated or invalid name ',a)
121	format('>>> ERROR: Atom',i3,' has duplicated or invalid number.')
125	format('>>> ERROR: First atom named ',a,' in bond ',i3,' not found')
126	format('>>> ERROR: Second atom named ',a,' in bond ',i3,' not found')
127	format('>>> ERROR: Head atom named ',a,' not found')
128	format('>>> ERROR: Tail atom named ',a,' not found')
130 format('>>> Warning: fractional charge = ',f8.5, ' on group no.',i3)
131	format('>>>>> ERROR: Atom',i3,' named ',a,' of charge group',i3,' not found.')
150	format('>>> ERROR: charge group atom count does not match total atom count!')
160	format('>>> ERROR: Sum of charge group charges not equal to sum of all charges!')

end subroutine readlib

!---------------------------------------------------------------------------------

subroutine check_overload(resnam)
	! arguments
	character(len=80)			::	resnam
	!locals
	integer						::	i

	do i = 1, nlibres
		if(lib(i)%nam == resnam) then
			write(*,100) trim(resnam)
100			format('>>> WARNING: Overloading old definition of ',a)
			lib(i)%nam = 'gone'
			exit
		end if
	end do
end subroutine check_overload

!---------------------------------------------------------------------------------

subroutine oldreadparm(flag)
!	arguments
	LOGICAL flag
! *** local variables
	integer i, ityp, filestat
	character(200)				::	line
	character(len=KEYLENGTH)	::	taci, tacj, tack, tacl
	integer						::	iaci, iacj, iack, iacl
!.......................................................................

	!set default values for options which are not available in old param. file
	iuse_switch_atom = 1 !this is the default. It can be overridden
	!in old libraries by setting switch atoms to 0
	!new parameter files have the keyword switch_atoms (sec. options)
	!to control this
	coulomb_constant = 332.0
	imp_explicit = .false.
	imp_type = 1 !harmonic impropers
	ff_type = FF_GROMOS

	flag = .true.
	CALL skip_comments(2)

	i = 0
	READ(2, *, err = 1000) nbnd_prm
	CALL skip_comments(2)
	if(allocated(bnd_prm)) deallocate(bnd_prm)
	if(allocated(bnd_types)) deallocate(bnd_types)
	allocate(bnd_prm(nbnd_prm), stat=alloc_status)
	allocate(bnd_types(nbnd_prm), stat=alloc_status)
	call check_alloc('bond parameters')
	do i = 1, nbnd_prm
		READ(2, *, err = 1000) iaci, iacj, bnd_prm(i)%prm%fk, bnd_prm(i)%prm%bnd0
		bnd_types(i)%taci = tac(iaci)
		bnd_types(i)%tacj = tac(iacj)
		bnd_prm(i)%SYBYLtype = '1 ' !set to default
	enddo
	write( * , 110) nbnd_prm, 'bond types'
	CALL skip_comments(2)
	i = 0
	READ(2, *, err = 1010) nang_prm
	CALL skip_comments(2)
	if(allocated(ang_prm)) deallocate(ang_prm)
	allocate(ang_prm(nang_prm), stat=alloc_status)
	if(allocated(ang_types)) deallocate(ang_types)
	allocate(ang_types(nang_prm), stat=alloc_status)
	call check_alloc('angle parameters')
	!set optional parameters to 0.
	ang_prm(:)%ureyfk = 0.
	ang_prm(:)%ureyr0 = 0.
	do i = 1, nang_prm
		read(2, '(a80)', err=1010) line
		read(line, *, iostat=filestat) iaci, iacj, iack, ang_prm(i)
		ang_types(i)%taci = tac(iaci)
		ang_types(i)%tacj = tac(iacj)
		ang_types(i)%tack = tac(iack)
		!accept missing parameters but not read error
		if(filestat > 0) goto 1010
	enddo
	write( * , 110) nang_prm, 'angle types'

	CALL skip_comments(2)
	i = 0
	READ(2, *, err = 1020) ntor_prm
	CALL skip_comments(2)
	if(allocated(tor_prm)) deallocate(tor_prm)
	allocate(tor_prm(ntor_prm), stat=alloc_status)
	if(allocated(tor_types)) deallocate(tor_types)
	allocate(tor_types(ntor_prm), stat=alloc_status)
	call check_alloc('torsion parameters')
	do i = 1, ntor_prm
		READ(2, *, err = 1020) iaci, iacj, iack, iacl, tor_prm(i)
		tor_types(i)%taci = tac(iaci)
		tor_types(i)%tacj = tac(iacj)
		tor_types(i)%tack = tac(iack)
		tor_types(i)%tacl = tac(iacl)
	enddo
	write(*, 110) ntor_prm, 'torsion types'

	CALL skip_comments(2)
	i = 0
	READ(2, *, err = 1030) nimp_prm
	CALL skip_comments(2)
	if(allocated(imp_prm)) deallocate(imp_prm)
	allocate(imp_prm(nimp_prm), stat=alloc_status)
	call check_alloc('improper parameters')
	do i = 1, nimp_prm
		READ(2, *, err = 1030) iacj, iack, imp_prm(i)%prm
		imp_prm(i)%taci = '' !not used
		imp_prm(i)%tacj = tac(iacj)
		imp_prm(i)%tack = tac(iack)
		imp_prm(i)%tacl = '' !not used
	enddo
	write( * , 110) nimp_prm, 'improper types'

	CALL skip_comments(2)
	i = 0
	READ(2, *, err = 1040) natyps

	READ(2, *, err = 1040) ivdw_rule
	READ(2, *, err = 1040) el14_scale
	CALL skip_comments(2)

	!if no PDB file read (blank topology)
	if(allocated(iaclib)) deallocate(iaclib)
	max_atyps = max_old_atyps
	allocate(iaclib(max_atyps))
	if(allocated(SYBYL_atom_type)) deallocate(SYBYL_atom_type)
	allocate(SYBYL_atom_type(max_atyps))
	if(allocated(tac)) deallocate(tac)
	allocate(tac(max_atyps))
	!clear atom type parameters
	do i=1,max_atyps
		iaclib(i)%avdw(:) = 0.
		iaclib(i)%bvdw(:) = 0.
		iaclib(i)%mass = 0.
		SYBYL_atom_type(i) = '     '
		tac(i) = ''
	end do

	call index_create(natyps)
	do i = 1, natyps
		READ(2, *, err = 1040) ityp, iaclib(ityp)%avdw(1), iaclib(ityp)%avdw(2), &
			iaclib(ityp)%bvdw(1), iaclib(ityp)%avdw(3), iaclib(ityp)%bvdw(3), &
			iaclib(ityp)%mass
		iaclib(ityp)%bvdw(2) = iaclib(ityp)%bvdw(1)
		write(taci,'(i4)') ityp
		tac(ityp) = adjustl(taci)
		if(.not. index_add(tac(ityp), ityp)) then
			write(*,130) tac(i)
		end if
	enddo
	write( * , 110) natyps, 'atom types'

	CALL skip_comments(2)
	i = 0
	READ(2, *, err = 1050) nlj2
	CALL skip_comments(2)
	if(allocated(lj2)) deallocate(lj2)
	allocate(lj2(nlj2))
	do i = 1, nlj2
	READ(2, *, err = 1050) lj2(i)%i, lj2(i)%j
	enddo
	write( * , 110) nlj2, 'LJ type 2 pairs'

	write( * , '(a,/)') 'Force field parameters successfully read.'
	return

  110 format ('Read ',i4,' ', a)
  120 format ('>>> ERROR: Failed to read ',a,' number ',i4)
  130 format ('>>> ERROR: Could not enumerate atom type ',a, '. Duplicate name?')

 1000 write( * , 120) 'bond type', i
	flag = .false.
	return

 1010 write( * , 120) 'angle type', i
	flag = .false.
	return

 1020 write( * , 120) 'torsion type', i
	flag = .false.
	return

 1030 write( * , 120) 'improper type', i
	flag = .false.
	return

 1040 write( * , 120) 'atom type', i
	flag = .false.
	return

 1050 write( * , 120) 'LJ type 2 pair', i
	flag = .false.
	return
!.......................................................................
end subroutine oldreadparm

!-----------------------------------------------------------------------

subroutine readff
	!read parameters
	if(.not.have_prm_file_name) then
		call get_string_arg(prm_file, '-----> Name of FF parameter file: ')
		have_prm_file_name = .true.
	endif
	inquire(file=prm_file, exist=have_prm_file_name)
	if(.not. have_prm_file_name) then
		write(*,10) trim(prm_file)
		return
	end if
10	format('>>>>> ERROR: ',a,' does not exist.')
	if(.not. prm_open(prm_file)) then
		!if reading fails try old format
		write(*,'(a)') '>>>>>WARNING: Attempting to read old-style parameter file.'
		if(openit(2, prm_file, 'old', 'formatted', 'read') /= 0) return
		CALL oldreadparm(ff_ok)
		close(2)
	else
		call readparm(prm_file) !this will set global flag topo_ok
		call prm_close
	end if
end subroutine readff

!-----------------------------------------------------------------------

subroutine readparm(filnam)
!	arguments
	character(*)				::	filnam

! *** local variables
	integer						::	i, ityp, j
	character(len=200)			::	line, restofline
	real(kind=prec)						::	rdummy
	logical						::	ldummy
	character(len=KEYLENGTH)	::	taci, tacj, tack, tacl
	integer						::	iaci, iacj, iack, iacl
	character(len=80)			::	section
	logical						::	SYBYL_warn = .false.
	integer						::	naliases
	type(BOND_PRM_TYPE)			::	bnd_prm_tmp
	type(ANGLIB_TYPE)			::	ang_prm_tmp
	type(TORLIB_TYPE)			::	tor_prm_tmp
!.......................................................................
	ff_ok = .false.

	write(*,10)
10	format(/,'Reading force field parameter file.')
!parameter file is opened by maketop
	if(.not. prm_open_section('options')) then
		write(*,*) '>>> WARNING: Options section not found in parameter file.'
		write(*,*) '    Trying to read old style parameter file:'
		call prm_close()
		if(openit(2, filnam, 'old', 'formatted', 'read') /= 0) return
		CALL oldreadparm(ff_ok)
		close(2)
		return
	end if

	if(prm_get_string_by_key('name', forcefield)) then
		write(*, 20) trim(forcefield)
	end if
20	format('Forcefield name:',t32,a)

	if(.not. prm_get_string_by_key('type', line)) then
		!also accept the keyword force_field for the same purpose
		ldummy = prm_get_string_by_key('force_field', line, 'GROMOS')
	end if
	call upcase(line)
	write(*, 22) trim(line)
	select case (line)
	case('GROMOS')
		ff_type = FF_GROMOS
	case('AMBER')
		ff_type = FF_AMBER
	case('CHARMM')
		ff_type = FF_CHARMM
	case default
		write(*,24)
		return
	end select

22	format('Forcefield type:',t32,a)
24	format('>>>>> ERROR: forcefield type must be one of GROMOS, AMBER, CHARMM')

	if(.not. prm_get_string_by_key('vdw_rule', line)) then
		write(*,*) '>>> ERROR: vdw_rule in options section not found in parameter file.'
		return
	end if

	call upcase(line)
	write(*,30) trim(line)
30	format('Lennard-Jones combination rule:',t32,a)

	if(line(1:9) == 'GEOMETRIC') then
		ivdw_rule = 1
	else if(line(1:10) == 'ARITHMETIC') then
		ivdw_rule = 2
	else
		write(*,*) '>>> ERROR: vdw_rule must be geometric or arithmetic.'
		return
	end if

	if(.not. prm_get_real_by_key('scale_14', rdummy)) then
		write(*,*) '>>> ERROR: scale_14 in options section not found in parameter file.'
		return
	end if
	el14_scale = real(rdummy,8)
	write(*,40) el14_scale
40	format('Scaling of 1-4 electrostatics:',t31,f6.3)

	if(prm_get_logical_by_key('switch_atoms', ldummy)) then
		if(ldummy) then
			iuse_switch_atom = 1
			write(*,50) 'switching atoms'
		else
			iuse_switch_atom = 0
			write(*,50) 'any atoms'
		end if
	else
		iuse_switch_atom = 1 !default
		write(*,50) 'switching atoms'
	end if
50	format('Cut-off rule:',t32,a)

	!read improper potential type (harmonic (default) or periodic)
	ldummy = prm_get_string_by_key('improper_potential', line, 'HARMONIC')
	call upcase(line)
	write(*,60) trim(line)
60	format('Improper potential:',t32,a)

	if(line == 'HARMONIC') then
		imp_type = 1
	else if(line == 'PERIODIC') then
		imp_type = 2
	else
		write(*,*) '>>> ERROR: improper_potential must be harmonic or periodic.'
		return
	end if

	!read improper definition scheme (automatic (default) or explicit)
	ldummy = prm_get_string_by_key('improper_definition', line, 'AUTOMATIC')
	call upcase(line)
	write(*,70) trim(line)
70	format('Improper generation:',t32,a)
	if(line == 'AUTOMATIC') then
		imp_explicit = .false.
	else if(line == 'EXPLICIT') then
		imp_explicit = .true.
	else
		write(*,*) '>>> ERROR: improper_definition must be automatic or explicit.'
		return
	end if

	if(.not. prm_get_real_by_key('coulomb_constant', rdummy)) then
		write(*,*) 'Coulomb constant is set to the default value of 332.'
		coulomb_constant = 332.0
	else
		coulomb_constant = real(rdummy,8)
		write(*,80) coulomb_constant
	end if
80	format('Coulomb constant:',t32,f8.4)

	section = 'atom_types'
	natyps = prm_count(section)
	if(natyps == 0) then
		write(*,*) '>>> ERROR: No atom_types defined.'
		return
	end if
	!if no PDB file read (blank topology)
	if(allocated(iaclib)) deallocate(iaclib)
	allocate(iaclib(natyps))
	if(allocated(SYBYL_atom_type)) deallocate(SYBYL_atom_type)
	allocate(SYBYL_atom_type(natyps))
	if(allocated(tac)) deallocate(tac)
	allocate(tac(natyps))
	!clear atom type parameters
	do i=1,natyps
		iaclib(i)%avdw(:) = 0.
		iaclib(i)%bvdw(:) = 0.
		iaclib(i)%mass = 0.
		SYBYL_atom_type(i) = '     '
		tac(i) = ''
	end do
	!change number of atom types in topology
	max_atyps = natyps

	!initialize atom type lookup table
	call index_create(natyps)

	do i = 1, natyps
		if(.not. prm_get_line(line)) goto 1040
		READ(line, *, iostat=j) tac(i), iaclib(i)%avdw(1), iaclib(i)%avdw(2), &
			iaclib(i)%bvdw(1), iaclib(i)%avdw(3), &
			iaclib(i)%bvdw(3), iaclib(i)%mass, SYBYL_atom_type(i)
		iaclib(i)%bvdw(2) = iaclib(i)%bvdw(1)
		if(j > 0) goto 1040
		if(j < 0) SYBYL_warn = .true.
		if(.not. index_add(tac(i), i)) then
			write(*,130) tac(i)
		end if
	enddo
	write( * , 110) natyps, 'atom types'
	if(SYBYL_warn) write(*,'(a)') &
		'>>>>> WARNING: No SYBYL name for one or more atom types.'
	SYBYL_warn = .false. !reset

	section = 'atom_aliases'
	naliases = prm_count(section)
	if(naliases > 0) call index_resize(natyps+naliases)
	do i = 1, naliases
		if(.not. prm_get_string_string(taci, tacj)) then
			goto 1060
		!make tac1 an alias for tac2
		elseif(.not. index_alias(taci, tacj)) then
			write(*,140) i, trim(section), tacj
			goto 1060
		end if
	end do
	write( * , 110) naliases, 'atom type alias names'

	j = 0
	do i = 1, natyps
		if(.not. index_get(tac(i), ityp)) then
			write(*,'(a,a,a)') 'ERROR: Failed to retrieve atom type ' , tac(i), ' from index.'
		else
			j = j + 1
		end if
	end do

	if(j<natyps) then
		write(*,*) '>>>>> ERROR: Bug in indexer, not all atom types found.'
		call prm_close
		return
	end if

	section = 'LJ_type2_pairs'
	nlj2 = prm_count(section)
	if(allocated(lj2)) deallocate(lj2)
	allocate(lj2(nlj2))
	do i = 1, nlj2
		if(.not. prm_get_string_string(taci, tacj)) goto 1050
		if(.not. index_get(taci, iaci)) then
			write(*,140) i, trim(section), taci
			goto 1050
		end if
		lj2(i)%i = iaci
		if(.not. index_get(tacj, iacj)) then
			write(*,140) i, trim(section), tacj
			goto 1050
		end if
		lj2(i)%j = iacj
	enddo
	write( * , 110) nlj2, 'LJ type 2 pairs'


	section = 'bonds'
	nbnd_types = prm_count(section)
	nbnd_prm = 0 !will accumulate this while reading

	if(allocated(bnd_prm)) deallocate(bnd_prm)
	if(allocated(bnd_types)) deallocate(bnd_types)
	allocate(bnd_prm(nbnd_types), stat=alloc_status)
	allocate(bnd_types(nbnd_types), stat=alloc_status)
	call check_alloc('bond parameters')

	bnd_prm_tmp%SYBYLtype = '   '




bondloop: do i = 1, nbnd_types

		if(.not. prm_get_line(line)) goto 1000
		read(line, *, iostat=j) taci, tacj, bnd_prm_tmp
		if(j>0) goto 1000
		if(j<0) SYBYL_warn = .true.
		if(.not. index_get(taci, iaci)) then
			write(*,140) i, trim(section), taci
			goto 1060
		end if
		if(.not. index_get(tacj, iacj)) then
			write(*,140) i, trim(section), tacj
			goto 1060
		end if
		!make sure to expand any aliases
		bnd_types(i)%taci = tac(iaci)
		bnd_types(i)%tacj = tac(iacj)
		do j = 1, nbnd_prm
			if(bnd_prm(j)%prm%fk == bnd_prm_tmp%prm%fk .and. bnd_prm(j)%prm%bnd0 == bnd_prm_tmp%prm%bnd0) then
				bnd_types(i)%cod = j
				cycle bondloop !skip to next
			endif
		end do
		nbnd_prm = nbnd_prm + 1
		bnd_prm(nbnd_prm) = bnd_prm_tmp
		bnd_types(i)%cod = nbnd_prm
	end do bondloop
	write( * , 110) nbnd_types, 'bond types'
	write( * , 110) nbnd_prm, 'unique bond parameters'
	if(SYBYL_warn) write(*,'(a)') &
		'>>>>> WARNING: No SYBYL bond type for one or more bond types.'
	SYBYL_warn = .false. !reset

	section = 'angles'
	nang_types = prm_count(section)
	nang_prm = 0 !accumulate unique parameters while reading
	if(nang_types == 0) then
		write(*,*) '>>> ERROR: Angles section not found in parameter file.'
		return
	end if
	if(allocated(ang_prm)) deallocate(ang_prm)
	allocate(ang_prm(nang_types), stat=alloc_status)
	if(allocated(ang_types)) deallocate(ang_types)
	allocate(ang_types(nang_types), stat=alloc_status)
	call check_alloc('angle parameters')
	!set optional parameters to 0.
	ang_prm(:)%ureyfk = 0.
	ang_prm(:)%ureyr0 = 0.
	!need to set optional parameters to 0 in ang_prm_tmp, above has no effect (???)
	ang_prm_tmp%ureyfk = 0.
	ang_prm_tmp%ureyr0 = 0.

angloop:do i = 1, nang_types
		if(.not. prm_get_line(line)) goto 1010
		read(line, *, iostat=j) taci, tacj, tack, ang_prm_tmp
		!accept missing parameters but not read error
		if(j > 0) goto 1060
		if(.not. index_get(taci, iaci)) then
			write(*,140) i, trim(section), taci
			goto 1060
		end if
		if(.not. index_get(tacj, iacj)) then
			write(*,140) i, trim(section), tacj
			goto 1060
		end if
		if(.not. index_get(tack, iack)) then
			write(*,140) i, trim(section), tack
			goto 1060
		end if
		!make sure to expand any aliases
		ang_types(i)%taci = tac(iaci)
		ang_types(i)%tacj = tac(iacj)
		ang_types(i)%tack = tac(iack)
		do j = 1, nang_prm
			if(ang_prm(j)%fk == ang_prm_tmp%fk .and. &
				ang_prm(j)%ang0 == ang_prm_tmp%ang0 .and. &
				ang_prm(j)%ureyfk == ang_prm_tmp%ureyfk .and. &
				ang_prm(j)%ureyr0 == ang_prm_tmp%ureyr0) then
				ang_types(i)%cod = j
				cycle angloop !skip to next
			endif
		end do
		nang_prm = nang_prm + 1
		ang_prm(nang_prm) = ang_prm_tmp
		ang_types(i)%cod = nang_prm
	end do angloop
	write( * , 110) nang_types, 'angle types'
	write( * , 110) nang_prm, ' unique angle parameters'

	section = 'torsions'
	ntor_types = prm_count(section)
	ntor_prm = 0 !accumulate unique parameters while reading
	if(ntor_types == 0) then
		write(*,*) '>>> ERROR: Torsions section not found in parameter file.'
		return
	end if
	if(allocated(tor_prm)) deallocate(tor_prm)
	allocate(tor_prm(ntor_types), stat=alloc_status)
	if(allocated(tor_types)) deallocate(tor_types)
	allocate(tor_types(ntor_types), stat=alloc_status)
	call check_alloc('torsion parameters')
torloop: do i = 1, ntor_types
		if(.not. prm_get_line(line)) goto 1020
		read(line, *, err = 1020) taci, tacj, tack, tacl, tor_prm_tmp
		if(index_get(taci, iaci, allow_wildcard=.true.)) then
			if(iaci > 0) then
				tor_types(i)%taci = tac(iaci)
			else
				!wildcard found
				tor_types(i)%taci = ''
			endif
		else
			!atom type not found and not a wildcard => error
			write(*,140) i, trim(section), taci
			goto 1060
		end if
		if(index_get(tacj, iacj, allow_wildcard=.true.)) then
			if(iacj > 0) then
				tor_types(i)%tacj = tac(iacj)
			else
				tor_types(i)%tacj = ''
			endif
		else
			write(*,140) i, trim(section), tacj
			goto 1060
		end if
		if(index_get(tack, iack, allow_wildcard=.true.)) then
			if(iack > 0) then
				tor_types(i)%tack = tac(iack)
			else
				tor_types(i)%tack = ''
			endif
		else
			write(*,140) i, trim(section), tack
			goto 1060
		end if
		if(index_get(tacl, iacl, allow_wildcard=.true.)) then
			if(iacl > 0 ) then
				tor_types(i)%tacl = tac(iacl)
			else
				tor_types(i)%tacl = ''
			endif
		else
			!this field may be 0, no it may not be!!!!!
			!we fail topology generation on any missing parameter
			!Paul Bauer 17092015
			write(*,140) i, trim(section), tacl
			goto 1060
		end if
		do j = 1, ntor_prm
			if(tor_prm(j)%fk == tor_prm_tmp%fk .and. &
				tor_prm(j)%rmult == tor_prm_tmp%rmult .and. &
				tor_prm(j)%deltor == tor_prm_tmp%deltor .and. &
				tor_prm(j)%paths == tor_prm_tmp%paths) then
				tor_types(i)%cod = j
				cycle torloop !skip to next
			endif
		end do
		ntor_prm = ntor_prm + 1
		tor_prm(ntor_prm) = tor_prm_tmp
		tor_types(i)%cod = ntor_prm

	end do torloop
	write(*, 110) ntor_types, 'torsion types'
	write(*, 110) ntor_prm, 'unique torsion parameters'

	section = 'impropers'
	nimp_prm = prm_count(section)
	if(ntor_prm == 0) then
		write(*,*) '>>> ERROR: Impropers section not found in parameter file.'
		return
	end if
	if(allocated(imp_prm)) deallocate(imp_prm)
	allocate(imp_prm(nimp_prm), stat=alloc_status)
	call check_alloc('improper parameters')
	do i = 1, nimp_prm
		if(.not. prm_get_line(line)) goto 1030
		read(line, *, err = 1030) taci, tacj, tack, tacl, imp_prm(i)%prm
		if(index_get(taci, iaci)) then
			imp_prm(i)%taci = tac(iaci)
		else
			!this field may be wild-card
			imp_prm(i)%taci = ''
		end if
		if(index_get(tacj, iacj)) then
			imp_prm(i)%tacj = tac(iacj)
		else
			!this field may be wild-card
			imp_prm(i)%tacj = ''
		end if
		if(index_get(tack, iack)) then
			imp_prm(i)%tack = tac(iack)
		else
			!this field may be wild-card
			imp_prm(i)%tack = ''
		end if
		if(index_get(tacl, iacl)) then
			imp_prm(i)%tacl = tac(iacl)
		else
			!this field may be wild-card
			imp_prm(i)%tacl = ''
		end if
	end do
	write( * , 110) nimp_prm, 'improper types'


	write( * , '(a,/)') 'Force field parameters successfully read.'
	ff_ok = .true.
	call prm_close()
	return

  110 format ('Read ',i4,' ', a)
  120 format ('>>> ERROR: Failed to read ',a,' number ',i4)
  130 format ('>>> ERROR: Could not enumerate atom type ',a, '. Duplicate name?')
  140 format ('>>> ERROR at line ',i3,' of section ',a,': Atom type ',a,' not found.')
1000 write( * , 120) 'bond type', i
	call prm_close()
	return

1010 write( * , 120) 'angle type', i
	call prm_close()
	return

1020 write( * , 120) 'torsion type', i
	call prm_close()
	return

1030 write( * , 120) 'improper type', i
	call prm_close()
	return

1040 write( * , 120) 'atom type', i
	call prm_close()
	return

1050 write( * , 120) 'LJ type 2 pair', i
	call prm_close()
	return

1060 write( * , 120) trim(section), i
	call prm_close()
	return
!.......................................................................
end subroutine readparm

!-----------------------------------------------------------------------

subroutine clearpdb

	nres = 0
	nres_solute = 0
	nmol = 0
	nat_pro = 0
	nat_solute = 0
	topo_ok = .false.
	if(allocated(makeH)) deallocate(makeH)
	if(allocated(heavy)) deallocate(heavy)
	if(allocated(xtop)) deallocate(xtop)
	if(allocated(res)) deallocate(res)
	if(allocated(istart_mol)) deallocate(istart_mol)

end subroutine clearpdb

!-----------------------------------------------------------------------

subroutine cleartop
	!forget the name of the parameter file
	have_prm_file_name = .false.
	have_solute_sphere = .false.
	have_title = .false.
	topo_ok = .false.
end subroutine cleartop

!-----------------------------------------------------------------------

logical function countpdb(pdb_fileno, atoms, residues, molecules)
!count atoms, residues and molecules in a pdb file

!arguments
	integer						::	pdb_fileno, atoms, residues, molecules

!locals
	integer						::	resno, oldno
	character(len=4)			::	resnam, atnam, oldresnam
	character(len=80)			::	line
	real(kind=prec)						::	xtmp(3)
	integer						::	atoms_in_res, atoms_in_file
	integer						::	rescode, oldrescode

! old 10	format(13x,a4,a4,i5,4x,3f8.3)
10	format(12x,a5,a4,1x,i4,4x,3f8.3)

	countpdb = .true.

	atoms = 0
	residues = 0
	molecules = 1
	atoms_in_file = 0
	rescode = 0
	oldrescode = 0
	oldno = 0
	do
		read(pdb_fileno,'(a)', end=100) line
		if(adjustl(line) == 'GAP' .or. line(1:6) == 'TER   ') then
			molecules = molecules + 1
			oldrescode = 0 !avoid implicit gap check
		else if(line(1:6) /= 'HETATM' .and. line(1:6) /= 'ATOM  ') then
			write(*,'(a,/,a)') '>>>WARNING: ignoring unrecognised line in PDB file', trim(line)
		else
			READ(line, 10, end = 100, err = 200) atnam, resnam, resno, xtmp(1:3)
			if(resno /= oldno) then
				if(rescode /= 0) then
					!check it
					if(atoms_in_res > lib(rescode)%nat) then
						write(*,12) oldresnam, oldno
						countpdb = .false.
					end if
					atoms = atoms + lib(rescode)%nat
				end if
				!lookup the new one
				do rescode = nlibres,1,-1
					if(lib(rescode)%nam == resnam) exit
				end do
				if(rescode == 0) then
					write(*,11) resno, resnam
					countpdb = .false.
				end if
				!check for implicit GAP
				if(oldrescode /= 0 .and. rescode /=0) then
					if(lib(oldrescode)%tail == 0 .or. &
						lib(rescode)%head == 0) then
						molecules = molecules+1
					end if
				end if
				!remember previous rescode
				oldrescode = rescode

				oldno = resno
				oldresnam = resnam
				atoms_in_res = 0
				residues = residues + 1
			else if(oldrescode /= 0 .and. resnam /= oldresnam) then
				!numbering problem?
				write(*, 13) oldresnam, resno, resnam
				countpdb = .false.
			end if
			atoms_in_file = atoms_in_file + 1
			atoms_in_res = atoms_in_res + 1
		end if
	!rewind but leave open
	end do

	!branch here at EOF
	!take care of last residue

100	if(rescode /= 0) then
		!check it
		if(atoms_in_res > lib(rescode)%nat) then
			write(*,12) oldresnam, oldno
			countpdb = .false.
		end if
		atoms = atoms + lib(rescode)%nat
	end if

	write(*,15) atoms_in_file, atoms
	rewind(pdb_fileno)
	if(atoms_in_file == 0) countpdb = .false.
	return

200	rewind(pdb_fileno) !error exit point
	write(*,16) line
	countpdb = .false.
16	format('>>>>> ERROR found in line: ', a80)

11	format('>>>>> ERROR: Residue number ',i5,' is of unknown type ',a4)
12	format('>>>>> ERROR: Too many atoms in residue ',a4,1x,i5)
13	format('>>>>> ERROR: Two residues with same number: ',a,i5,1x,a)

15	format( 'PDB file contains',i6,' atoms.', /, &
			'The number of atoms in the topology will be',i6)

end function countpdb

!-----------------------------------------------------------------------

subroutine readpdb()

! *** local variables
	character(len=256)		:: pdb_file
	CHARACTER(len=4)		:: atnam_tmp , resnam_tmp , resnam_tmp2
	character(len=80)			::	line
	integer resnum_tmp, oldnum, irec, i, atom_id(max_atlib), j , oldnum2 , resnum_tmp2
	real(kind=prec)				:: xtmp(3)
	LOGICAL res_found, at_found
	integer						::	first_res_of_mol
	logical						::	last_line_was_gap
	integer						::	atoms, residues, molecules
	type(RETTYPE)					:: glob
!.......................................................................

	write( *, * )
	call get_string_arg(pdb_file, '-----> Name of PDB file: ')
	if(openit(3, pdb_file, 'old', 'formatted', 'read') /= 0) return
 	REWIND(3)
!
!	PDB format(we need only atom name, res. name, number and coords):
!   The format is
!  1. |    1 -  6    |   A6    | Record ID (eg ATOM, HETATM)
!   2. |    7 - 11    |   I5    | Atom serial number
!   -  |   12 - 12    |   1X    | Blank
!   3. |   13 - 16    |   A4    | Atom name (eg " CA " , " ND1")
!   4. |   17 - 17    |   A1    | Alternative location code (if any)
!   5. |   18 - 20    |   A3    | Standard 3-letter amino acid code for residue
!   -  |   21 - 21    |   1X    | Blank
!   6. |   22 - 22    |   A1    | Chain identifier code
!   7. |   23 - 26    |   I4    | Residue sequence number
!   8. |   27 - 27    |   A1    | Insertion code (if any)
!   -  |   28 - 30    |   3X    | Blank
!   9. |   31 - 38    |  F8.3   | Atom's x-coordinate
!  10. |   39 - 46    |  F8.3   | Atom's y-coordinate
!  11. |   47 - 54    |  F8.3   | Atom's z-coordinate
!  12. |   55 - 60    |  F6.2   | Occupancy value for atom
!  13. |   61 - 66    |  F6.2   | B-value (thermal factor)
!

!Old format line  10	FORMAT(13x,a4,a4,i5,4x,3f8.3)
!Format if all is read  10	format(a6,i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3)
! Changed format so that residues for 4 letter codes can also be read P. Bauer
10	format(13x,a4,a4,1x,i4,4x,3f8.3)

!	progress output formats
20	format('molecule ',i4,': ',a4,i5)
21	format(' - ',a4,i5)

22	format('>>>>> ERROR: The check of the PDB file failed.')
23	format('>>> WARNING: Multiple GAP lines at line',i5)
	oldnum = 0

	!get rid of old topology but not FF params
	call topo_deallocate(keep_ff=.true.)

	if(.not. countpdb(3, atoms, residues, molecules)) then
		write(*,22)
		close(3)
		return
	end if

	CALL clearpdb !get rid of old PDB data.
	call allocate_for_pdb(atoms, residues, molecules) !make space for new topology
	!clear hydrogen make flags
	makeH(:) = .false.
	!reset molecule counter
	nmol = 1
	istart_mol(1) = 1
	glob%new = .false.
	last_line_was_gap = .false.
	irec = 0
	do
		read(3,'(a)', end=100) line
		irec = irec + 1
		if(adjustl(line) == 'GAP' .or. line(1:6) == 'TER   ') then
			if(last_line_was_gap) then
!			if(glob%new) then
				write(*, 23) irec
			else
				!set gap flag - new molecule will be recognised later
!				glob%new = .true.
				last_line_was_gap = .true.
			end if
		else if(line(1:6) /= 'HETATM' .and. line(1:6) /= 'ATOM  ') then
			!do nothing
		else
			resnum_tmp2 = resnum_tmp
			resnam_tmp2 = resnam_tmp
			READ(line, 10, end = 100, err = 200) atnam_tmp, resnam_tmp, &
				resnum_tmp, xtmp(1:3)

		! ---	New residue ?
			if(resnum_tmp/=oldnum) then
		! ---	   Check if old residue was OK...
				if(nres/=0) then
					do i = 1, lib(res(nres)%irc)%nat
						if(lib(res(nres)%irc)%atnam(i)(1:1) /='H') then !heavy
							heavy(res(nres)%start - 1 + i) = .true.
							makeH(res(nres)%start - 1 + i) = .false.
							if(atom_id(i) == 0) then !not found
								write(*, '(/,a,a,a,i5,/)') '>>> Heavy atom ',&
									lib(res(nres)%irc )%atnam(i), &
									' missing in residue ', oldnum
								goto 210
							end if
						else !hydrogen
							heavy(res(nres)%start - 1 + i) = .false.
							!flag hydrogens to be generated
							if(atom_id(i) == 0) then
								!it was not in the file and needs to be made
								makeH(res(nres)%start - 1 + i) = .true.
							end if
						end if
					enddo
					if (nres>1) then
					glob=is_new_mol(nres,nres-1,glob%new,oldnum2,resnam_tmp2,resnum_tmp2,glob%first)
					if(last_line_was_gap) then
						glob%new=.true.
						last_line_was_gap = .false.
					end if
					end if
				endif

				!look up new residue in library
				res_found = .false.
				nres = nres + 1
				res(nres)%start = nat_pro + 1
				do i = 1, nlibres
					if(resnam_tmp==lib(i)%nam ) then
						res_found = .true.
						res(nres)%name = lib(i)%nam
						res(nres)%irc = i
						exit
					endif
				enddo
				if(.not.res_found) then
					write( * , '(/,a,a,/)') '>>> Residue not found in library: ', &
						resnam_tmp
					goto 210
				endif
				!clear atom read flags
				do i = 1, lib(res(nres)%irc )%nat
					atom_id(i) = 0
				enddo

                                if (nres .eq. 1 ) then
                                        if(nat_pro + 1 == istart_mol(nmol)) then
                                                write(*,20,advance='no') nmol, resnam_tmp, resnum_tmp
                                                first_res_of_mol = resnum_tmp
                                        end if
                                end if


				nat_pro = nat_pro + lib(res(nres)%irc )%nat
				!set nat_solute = nat_pro unless residue is water
				if(index(solvent_names, trim(resnam_tmp)) == 0) then
!					nat_solute = nat_pro
					nat_solute = nat_solute + lib(res(nres)%irc )%nat
					nres_solute = nres_solute + 1
					if (nat_solute .ne. nat_pro) then
!Means we are having a molecule of solvent between solute molecules
						write(*,*) '>>>> ERROR: Solvent molecule between solute molecules'
						write(*,'(a,i4,a)') 'Please remove molecule ',nres-1,' and try again'
						stop 'Abnormal termination of Qprep'
					end if
				end if
!Check if we have previously set a solvent molecule and are now setting
!again a solute molecule
!If yes, abort topology generation to avoid mixing of solvent and solute descriptors
				oldnum2 = oldnum
				oldnum = resnum_tmp

			end if !new residue, first part, make check for new molecule after all information
				!has been read in

				!check for implicit GAP, i.e. this residue and previous
				!have no tail or head connections, respectively.
				!Do this unless this is the first residue
				!Don't do it if previous line was gap (already done)


			at_found = .false.
			do i = 1, lib(res(nres)%irc )%nat
				if(atnam_tmp==lib(res(nres)%irc)%atnam(i)) then
					at_found = .true.
					atom_id(i) = 1
					j = res(nres)%start - 1 + i
					xtop(j * 3 - 2) = xtmp(1)
					xtop(j * 3 - 1) = xtmp(2)
					xtop(j * 3) = xtmp(3)
					exit
				endif
			enddo !i

			if(.not.at_found) then
				write( * , '(/,a,a,a,i5,a,a,/)') '>>> Atom ', atnam_tmp, &
				' in residue no. ', resnum_tmp, ' not found in library entry for ', &
				lib(res(nres)%irc)%nam
				goto 210
			endif !.not.at_found
		end if !line type
	enddo

!branch here at EOF
100	close(3)
	if(nres>1) then
		glob=is_new_mol(nres,nres-1,glob%new,oldnum2,resnam_tmp2,resnum_tmp2,glob%first)
	end if
	if(last_line_was_gap) then
		write(*,'(a)') '>>> Warning: PDB file ends with GAP line.'
		nmol = nmol - 1 !correct molecule count
	else
		!print last residue of last molecule if more than one
		if(glob%first /= oldnum) then
			write(*,21) resnam_tmp, oldnum
		else
			write(*,*)
		end if
	end if

! --- Check if last residue was OK...
	if(nres > 0) then !avoid crashing when reading empty file
		do i = 1, lib(res(nres)%irc )%nat
			if(lib(res(nres)%irc )%atnam(i)(1:1) /='H') then
				!it is a heavy (non-H) atom
				heavy(res(nres)%start - 1 + i) = .true.
				makeH(res(nres)%start - 1 + i) = .false.
				if(atom_id(i) == 0) then !it was not found
					write( * , '(/,a,a,a,i5,/)') '>>> Heavy atom ', &
					lib(res(nres)%irc )%atnam(i),' missing in residue ',oldnum
					goto 210
				end if
			else !it is a hydrogen
				heavy(res(nres)%start - 1 + i) = .false.
				if(atom_id(i) == 0) then
					!it was not read from the file and should be gerenated
					makeH(res(nres)%start - 1 + i) = .true.
				end if
			end if
		enddo
	end if

	nwat = (nat_pro - nat_solute) / solv_atom
	write(*,110) nmol,nres, nres_solute, nat_pro, nat_solute
110	format(/,'Successfully read PDB file with', i5,' molecules,',/,&
		i5,' residues totally (',i5,' in solute).',/,&
		i5,' atoms totally (',i5,' in solute).')

	write( * , '(/,a,/)') 'Sequence listing:'
	write( * , '(16(a4,1x))') res(1:nres)%name
	write( *, * )

	return

!	error handling
  200 write( * , '(/,a,i5,/)') '>>> Read error on line ', irec
  210 write( * , '(a)') 'Correct the PDB file and try readpdb again!'
	write( * , '(a)') 'If the problem is in the library you need to correct it and do',&
		'clearlib and readlib before trying readpdb again.'
	close(3)
	CALL clearpdb
	pdb_file = ''
end subroutine readpdb

!move checking of new molecule to separate subroutine
type(RETTYPE) function is_new_mol(thisres,prevres,gap,old,tnam,tnum,firstres)
	! *** input variables
	integer					:: thisres,prevres,old,tnum,firstres
	logical					:: gap
	character*4				:: tnam
	! *** local variables
	real*8					:: dist,xdist,ydist,zdist
	integer					:: dhead,dtail
	integer,allocatable			:: start_temp(:)
	type(RETTYPE)				:: local

20      format('molecule ',i4,': ',a4,i5)
21      format(' - ',a4,i5)


	!check for implicit GAP, i.e. this residue and previous
	!have no tail or head connections, respectively.
	!Do this unless this is the first residue
	!Don't do it if previous line was gap (already done)

	if(.not.gap) then
		if (lib(res(prevres)%irc)%tail .ne. 0 .and. &
			lib(res(thisres)%irc)%head .ne. 0 ) then
			dtail = lib(res(prevres)%irc)%tail
			dhead = lib(res(thisres)%irc)%head
			dtail = res(prevres)%start - 1 + dtail
			dhead = res(thisres)%start - 1 + dhead
			dtail = (dtail*3)-3
			dhead = (dhead*3)-3
			xdist = (xtop(dtail+1)-xtop(dhead+1))
			ydist = (xtop(dtail+2)-xtop(dhead+2))
			zdist = (xtop(dtail+3)-xtop(dhead+3))
			dist  = sqrt(xdist**2 + ydist**2 + zdist**2)
		end if
	end if
        if(gap .or. &
		(lib(res(prevres)%irc)%tail == 0) .or. &
		(lib(res(thisres)%irc)%head == 0) .or. &
		dist .gt. 5 ) then 
		nmol = nmol + 1
		if (nmol .gt. size(istart_mol)) then
			allocate(start_temp(nmol))
			start_temp(1:nmol-1)=istart_mol(1:nmol-1)
			deallocate(istart_mol)
			allocate(istart_mol(nmol))
			istart_mol(1:nmol-1)=start_temp(1:nmol-1)
			deallocate(start_temp)
		end if
		istart_mol(nmol) = res(thisres)%start 
		if(firstres /= old) then
			write(*,21) res(prevres)%name, old
		else
			write(*,*)
		end if
		!if new molecule write output
		write(*,20,advance='no') nmol, tnam, tnum 
		firstres = tnum
	end if
	local%new = .false.
	local%first = firstres
	is_new_mol = local
end function  is_new_mol

!-----------------------------------------------------------------------

subroutine readtop
	! *** local variables
	character(len=80)			::	filnam
	character(len=200)			::	files_to_load
	integer						::	 i, j, ires
	logical						::	loaded, res_found
	integer						::	fn_start, fn_end, totlen
	character(len=4)            ::  resnam_tmp

	CALL get_string_arg(filnam, '-----> Give name of(old) topology file: ')
	if(openit(10, filnam, 'old', 'formatted', 'read') /= 0) return
	!tell topo_read to allocate some extra space for bonds
	topo_ok = topo_read(u=10, require_version=2.0, extrabonds=max_extrabnd)
	close(10)
	if(.not. topo_ok) return

	have_title = .true.
	call set_default_mask

	!sphere is defined if topology is version 4 or later
	if(version >= 4) then
		have_solute_sphere = .true.
		have_solvent_boundary = .true.
	end if

	!boundary condition is defined if topology is version 4.30 or later
	if(version >= 4.30) then
	boundary_set = .true.
	end if

	coord_source = 'topology'
	auto_name = filnam !for automatic naming of pdb and mol2 files
	nwat = (nat_pro - nat_solute) / solv_atom
	!reset makeH - don't need to make any hydrogens
	deallocate(makeH, stat=alloc_status)
	allocate(makeH(nat_pro))
	makeH(:) = .false.

	if(nlibres == 0) then
		!auto-load libraries
		write(*,100)
100		format('Auto-loading libraries:')
		files_to_load = lib_files
		!clear lib file list
		lib_files = ''
		fn_start = 1
		totlen = len_trim(files_to_load)
		do while(fn_start < totlen)
			fn_end = string_part(files_to_load, ';', fn_start)
			filnam = files_to_load(fn_start:fn_end)
			call readlib(trim(filnam))
			fn_start = fn_end + 2
		end do
	end if

	!Annotate irc to residues to find atoms in residues.
	do ires = 1, nres
		res_found = .false.
		resnam_tmp = res(ires)%name
		do i = 1, nlibres
			if(resnam_tmp==lib(i)%nam ) then
				res_found = .true.
				res(ires)%irc = i
				exit
			end if
		enddo
		if(.not.res_found) then
			write( * , '(/,a,a,/)') '>>> ERROR Residue not found in loaded libraries: ', &
				resnam_tmp
		end if
	end do

end subroutine readtop

!-----------------------------------------------------------------------

subroutine readx
! fix for handling new restart files
! needs to be done for all functions
! *** local variables
	CHARACTER(len=80)			::	filnam
	integer						::	i
	integer(4)					::	nat3,canary
	integer						::	u
        logical                                         :: old_restart = .false.


	CALL get_string_arg(filnam, '-----> Give name of Q-coordinate file: ')
	u = freefile()
	if(openit(u, filnam, 'old', 'unformatted', 'read') /= 0) return

        !first try to read canary
        read(u) canary
        if ((canary .ne. -137).and.(canary.ne.-1337).and.(canary.ne.-13337)) then
        ! old restart
                rewind (u)
                old_restart = .true.
        end if

	!read number of co-ordinates (=3*natom)
	READ(u) nat3
	rewind(u)
	natom = nat3 / 3
	if(natom/=nat_pro) then
		write(*, 20) natom
20		format('>>> Error:',i6,' atoms in restart file does not match topology.')
		close(u)
		return
	endif
        if (.not.old_restart) read(u) canary
	READ(u) nat3,xtop(1:nat3)

	write( * , '(a,/)') 'Coordinate file successfully read.'
	close(u)
	write(coord_source, 100) trim(filnam)
	auto_name = filnam

100	format('restart file ', a)
end subroutine readx

!-----------------------------------------------------------------------

subroutine readnext
	integer						::	filestat, nat3
	logical						::	isopen
	character(len=10)			::	buf

	!trj_frame is the number of the record about to be read
	if(.not. trj_read(xtop)) then
		call trj_close
		trj_frame = 0
		write(*,900)
900		format('>>>>> ERROR: Trajectory read failure')
	else
		write(buf,'(i8)') trj_frame
		write( * , '(a,a)') 'Loaded trajectory frame ',trim(buf)
		write(coord_source, 100) trim(trj_filnam), trj_frame
		write(auto_name, 110) trim(trj_filnam),trim(adjustl(buf))
		trj_frame = trj_frame + 1
	end if

100	format('trajectrory file ', a,' frame ',i6)
110	format(a,'.',a) !for auto file naming trj_filnam.frame
end subroutine readnext

!-----------------------------------------------------------------------

subroutine readframe
! *** local variables
	integer						::	frame

	frame = get_int_arg('-----> Give frame number: ')

	if(trj_seek(frame)) then
		trj_frame = frame
		call readnext
	end if

end subroutine readframe

!-----------------------------------------------------------------------

subroutine trajectory
!locals
	character(len=80)			::	reply

	if(.not. topo_ok) then
		!can't do this without a topology
		return
	end if

	CALL get_string_arg(trj_filnam, '-----> Give name of trajectory file: ')

	if(trj_open(trj_filnam)) then
		call get_string_arg(reply, '-----> Use the trajectory atom mask (y/n)? ')
		if(reply == 'y' .or. reply == 'Y') then
			call trj_clone_mask(mask)
			write(*,100) mask%included
		else
			write(*,110) trj_count(), trj_intersection(mask)
		end if
	end if
	trj_frame = 1
100	format('There are',i6,' atoms in the trajectory.')
110	format('There are',i6,' atoms in the trajectory.',i6,&
		' of these are in the mask')
end subroutine trajectory

!-----------------------------------------------------------------------

subroutine modify_mask
!locals
	character(len=80)			::	mask_def
	integer						::	atoms_added

	if(.not. topo_ok) then
		write(*,900)
900		format('>>>>> ERROR: Cannot define mask without topology.')
	else
		write(*,120) mask%included
120		format('The current mask contains',i6,' atoms.')
		call get_line_arg(mask_def, &
			'-----> Enter atom mask (''none'' to clear):')
		call upcase(mask_def)
		if(mask_def == 'NONE' .or. mask_def == '0') then
			call mask_clear(mask)
			write(*,110)
110			format('Mask cleared.')
		else
			atoms_added = mask_add(mask, mask_def)
			write(*,100) atoms_added, mask%included
100			format(i6,' atoms added to mask which now contains',i6,' atoms.')
		end if
	end if
end subroutine modify_mask

!-----------------------------------------------------------------------

subroutine set_cgp
! *** local variables
	integer ires, igp, i, ntot, ntot_solute, i3
	integer						::	switchatom, ia, ncgp_skipped = 0
	real(kind=prec)						::	r2
	real(kind=prec)						::	cgp_cent(3)
	integer						::	nheavy

	ntot = 0
	ntot_solute = 0
	ncgp = 0
	ncgp_solute = 0

	!do exclusions here as well
	nexats = 0
	nexwat = 0
	if (.not. use_PBC) excl(:) = .false.

	!set chargegroups for all atoms!
	do ires = 1, nres
		do igp = 1, lib(res(ires)%irc)%ncgp
			if(use_PBC ) then
				r2 = -1 !do not exclude
			elseif(iuse_switch_atom == 1) then
				!exclude on basis of switching atom
				switchatom = res(ires)%start - 1 + lib(res(ires)%irc)%switch(igp)
				if(.not. heavy(switchatom)) then
					topo_ok = .false.
					write(*,900) igp, lib(res(ires)%irc)%nam
900	format('>>>>> ERROR: Switching atom of charged group',i2,' of residue type ',a4,' is not a heavy atom.')
					cycle !continue looping, do not exclude
				endif
				i3 = 3*switchatom-3
				r2 = ( xtop(i3+1) - xpcent(1) )**2 &
					+( xtop(i3+2) - xpcent(2) )**2 &
					+( xtop(i3+3) - xpcent(3) )**2
			else
				!exclude by charge group centre
          cgp_cent(:) = 0
				  do i = 1, lib(res(ires)%irc)%natcgp(igp)
			    	ia = res(ires)%start - 1 + lib(res(ires)%irc)%atcgp(i, igp)
					  i3 = ia*3
					  cgp_cent(:) = cgp_cent(:) + xtop(i3-2:i3)
				  end do
          cgp_cent(:) = cgp_cent(:)/real(lib(res(ires)%irc)%natcgp(igp))
          r2 = dot_product(cgp_cent(:)-xpcent(:),cgp_cent(:)-xpcent(:))
			endif

			if( (r2 < rexcl_o**2) .or. &
				(ires > nres_solute .and. r2 < (rexcl_o + 2.)**2)) then
				!is it inside? For solvent include up to 2Å outside rexcl_o too
				ncgp = ncgp + 1
				if(ires <= nres_solute) then
					ncgp_solute = ncgp_solute + 1
					ntot_solute = ntot_solute + lib(res(ires)%irc)%natcgp(igp)
				endif
				cgp(ncgp)%first = ntot + 1
				cgp(ncgp)%last = ntot + lib(res(ires)%irc)%natcgp(igp)
				cgp(ncgp)% iswitch = res(ires)%start - 1 + lib(res(ires)%irc)%switch(igp)
				do i = 1, cgp(ncgp)%last - cgp(ncgp)%first + 1
					ntot = ntot + 1
					cgpatom(ntot) = res(ires)%start - 1 + lib(res(ires)%irc)%atcgp(i, igp)
				enddo
			else
				nexats = nexats + lib(res(ires)%irc)%natcgp(igp)
				if(ires > nres_solute) then
					nexwat = nexwat + 1
				endif
				ncgp_skipped = ncgp_skipped + 1
				do i = 1, lib(res(ires)%irc)%natcgp(igp)
					excl(res(ires)%start - 1 + lib(res(ires)%irc)%atcgp(i, igp)) = .true.
				enddo
			endif
		enddo
	enddo

	if(ntot + nexats /= nat_pro) then
	write( * , '(/,a,/)') '>>> Error in set_cgp: ntot+nexats /= nat_pro.'
		STOP
	endif

	write(*, 100) 'Made', ncgp_solute, ' solute', ntot_solute
	write(*, 100) 'Made', ncgp-ncgp_solute, ' solvent', ntot - ntot_solute
	if(nexats > 0) then
		write(*, 100) 'Skipped', ncgp_skipped, 'excluded', nexats
	endif
	if(nexwat > 0) then
		write(*, 110) nexwat
	endif

	if (.not. use_PBC) then
		rwat = rwat_eff()
		write(*,140) rwat
	endif

100 format(a,t8,i6,a,t23,'charge groups comprising a total of ', i6,' atoms.')
110 format('Excluded ',i6,' water molecules.')
140	format('Effective water radius                     :',f6.2)
end subroutine set_cgp

!-----------------------------------------------------------------------

subroutine set_crg
! *** local variables
	integer ires, i, ntot

	ntot = 0

	do ires = 1, nres
	do i = 1, lib(res(ires)%irc )%nat

	ntot = ntot + 1
	crg(ntot) = lib(res(ires)%irc )%crg_lib(i)

	enddo
	enddo

	if(ntot/=nat_pro) then
	write( * , '(/,a,/)') '>>> Error in set_crg: ntot /= nat_pro.'
		STOP
	endif

	write(*, 100) ntot
100 format('Assigned',i6,' atomic charges.')
end subroutine set_crg

!-----------------------------------------------------------------------

subroutine set_iac
! *** local variables
	integer						::	ires, i, ntot, iaci
	logical						::	used(max_atyps)

	used(:) = .false.
	ntot = 0
	!set atom types
	do ires = 1, nres
		do i = 1, lib(res(ires)%irc )%nat
			ntot = ntot + 1
			if(.not. index_get(lib(res(ires)%irc)%tac_lib(i), iaci)) then
				write(*,110) lib(res(ires)%irc)%tac_lib(i), i, ires
				topo_ok = .false.
			else
				used(iaci) = .true.
				iac(ntot) = iaci
			end if
		enddo
	enddo

	!this should never happen
	if(ntot/=nat_pro) then
		write( * , '(/,a,/)') '>>>>> ERROR in set_iac: ntot /= nat_pro.'
		topo_ok = .false.
	endif

	do i=1, natyps
		if(used(i) .and. iaclib(i)%mass <= 0) then
			write(*, 120) tac(i), iaclib(i)%mass
			topo_ok = .false.
		end if
	end do

	if(topo_ok) write(*, 100) nat_pro
100	format('Assigned',i6,' atom types.')

110	format('>>>>> ERROR: Atom type ',a8,' not defined (atom',i3,' of residue',i5,'.)')
120	format('>>>>> ERROR: Atom type ',a8,' has mass ',f6.2)
end subroutine set_iac

!-----------------------------------------------------------------------

logical function set_irc_solvent()

	set_irc_solvent = .false.

	do irc_solvent = 1, nlibres
		if(lib(irc_solvent)%nam == solvent_name) exit
	end do
	if(irc_solvent > nlibres) then
		write(*,90) solvent_name
90		format('>>>>> ERROR: No solvent library entry named ',a,&
			' has been loaded.')
		return
	end if
	write(*,95) solvent_name
95	format('Library entry used to generate solvent  : ',a10)

	set_irc_solvent = .true.

end function set_irc_solvent

!-----------------------------------------------------------------------

subroutine set_solvent_type
	!set solvent type (SPC-like, TIP3P-like or general)
	!also set the number of atoms per solvent molecule
	!has to be done here to prevent the program from crashing 
	!if it was not set before, better save than sorry
	!locals
	integer						::	irc_solvent
	integer						::	i
	integer,allocatable				::	iac_H(:)
	logical						::	dummy, fail = .false.

	solvent_type = SOLVENT_GENERAL

	!any solvent at all?
	if(nwat == 0) return
	
	!look up residue code of 1st solvent molecule
	irc_solvent = res(nres_solute + 1)%irc
	!get number of atoms in each solvent molecule
        solv_atom = lib(irc_solvent)%nat

	!check if all solvent molecules are of same type
	do i=nres_solute + 2, nres
		if(res(i)%irc /= irc_solvent) return
	end do

	!now we know all solvent molecules are the same type of n-atomic molecules
	solvent_type = SOLVENT_ALLATOM

	!further checks for SPC water, for which special optimisations can be used
	!check if atoms 2 and 3 are hydrogens
iloop:	do i = 2,solv_atom
	if(lib(irc_solvent)%atnam(i)(1:1) /= 'H') then
	fail = .true.
	end if
	end do iloop
	if (fail) return
	allocate(iac_H(solv_atom-1))

	do i = 2,solv_atom
	dummy = index_get(lib(irc_solvent)%tac_lib(i), iac_H(i-1))
	end do
	do i = 2,solv_atom-1
	if ((iac_H(i-1) /= iac_H(i)).or.(iaclib(iac_H(i))%avdw(1) .ne. 0.0_prec) .or. (iaclib(iac_H(i))%bvdw(1) .ne. 0.0_prec)) then
		fail = .true.
	end if
	if (fail) then
		deallocate(iac_H)
		return
	end if
	end do
		solvent_type = SOLVENT_SPC
end subroutine set_solvent_type

!-----------------------------------------------------------------------

integer function get_atom_from_descriptor(aid)
	!arguments
	character(*), intent(in)	::	aid	!string=residue:atom

	!locals
	integer						::	separator_pos
	character(len=20)			::	res_str
	character(len=5)			::	atom_str
	integer						::	filestat
	integer						::	resnum, atnum

	get_atom_from_descriptor = 0

	separator_pos = scan(aid, ':')
	if(separator_pos < 2 .or. separator_pos == len_trim(aid)) return !no valid colon found
	res_str = aid(1:separator_pos-1)
	atom_str = aid(separator_pos+1:len_trim(aid))
	read(res_str, *, iostat=filestat) resnum
	if(filestat > 0) return
	get_atom_from_descriptor = find_atom(resnum, atom_str)

end function get_atom_from_descriptor

!-----------------------------------------------------------------------

!Routine to call from other modules
subroutine define_boundary_condition
   if (.not. set_boundary_condition()) then
	  call parse_reset ! clear command line
	  return
   end if
end subroutine define_boundary_condition

!-----------------------------------------------------------------------

!Set the boundary condition & read parameterts connected to the boundary.
logical function set_boundary_condition()

	!locals
	character(len=80)		::	kind_of_boundary

	!ask if to use simulation sphere or periodic box
	call get_string_arg(kind_of_boundary, &
			'-----> Select kind of boundary: simulation sphere(1) or periodic box(2):')
	call upcase(kind_of_boundary)
	select case(kind_of_boundary)

	!sphere selected
	case('SPHERE', '1')
		use_PBC = .false. !set the flag
		!set the center and radius of the sphere
		if(.not. set_simulation_sphere()) then
			call parse_reset !clear command line
			return
		end if
		boundary_set = .true.

	!box selected
	case('BOX', '2')
		use_PBC = .true.
		!set the center and radius of the box
		if(.not. set_solvent_box()) then
			call parse_reset !clear command line
			return
		end if
		boundary_set = .true.

	case default
		write(*,'(a)') '>>>>> ERROR: Unknown boundary'
		set_boundary_condition = .false.
		return
	end select

	set_boundary_condition = .true.

end function set_boundary_condition

!-----------------------------------------------------------------------

logical function set_simulation_sphere()

	!get centre co-ordinates
	!as residue:atom or x y z

	!locals
	character(len=80)			::	line
	integer						::	filestat
	integer						::	centre_atom
	real(kind=prec)						::	rwat_in, xwat_in

	set_simulation_sphere = .false.

	call get_string_arg(line, '-----> Sphere centre (<x y z> or <residue:atom_name> or <"mass">): ')
	if(scan(line, ':') > 0) then !got res:at
		centre_atom=get_atom_from_descriptor(line)
		if(centre_atom == 0) then
			write(*,900) trim(line)
900			format('>>>>> ERROR: Could not find centre atom ',a)
			return
		end if
		xpcent(:) = xtop(3*centre_atom-2:3*centre_atom)
	elseif(scan(line, 'mass') > 0) then   !define center by center of mass
		if (.not. get_centre_by_mass(xpcent(:))) then
			write(*,*) ('>>>>> ERROR: Could not create centre ')
			return
		end if
	else !got x
		read(line, *, iostat=filestat) xpcent(1)
		if(filestat > 0) then !invalid x coordinate
			return
		end if
		xwat_in=get_real_arg('-----> Sphere centre y: ')
		xpcent(2) = xwat_in
		xwat_in=get_real_arg('-----> Sphere centre z: ')
		xpcent(3) = xwat_in
	end if

	rwat_in = get_real_arg('-----> Simulation sphere radius: ')
	rexcl_o = rwat_in
	if(rexcl_o <= 0.) then
		write(*,800) rexcl_o
800		format('>>>>> ERROR: Invalid radius ',f8.2)
		return
	end if

    write(*,*)
	write(*,100) xpcent(:)
    write(*,110) rexcl_o
100	format('Simulation sphere centre                   :   ',3f8.3)
110 format('Simulation radius                          :   ',f8.3)
	set_simulation_sphere = .true.



end function set_simulation_sphere

!-----------------------------------------------------------------------

logical function set_solvent_box()

	!locals
	integer					::	i !loop index
	character(len=80)		::	line
	integer					::	centre_atom
	integer					::	filestat
	real(kind=prec)					::	xwat_in, coord_in

	have_solvent_boundary = .false.
	set_solvent_box = .false.


	!the centre of the PCB-box
	call get_string_arg(line, '-----> Periodic box centre (<x y z> or <residue:atom_name> or <"mass">): ')
	call upcase(line)
	if(scan(line, ':') > 0) then !got res:at
		centre_atom=get_atom_from_descriptor(line)
		if(centre_atom == 0) then
			write(*,900) trim(line)
900			format('>>>>> ERROR: Could not find centre atom ',a)
			return
		end if
		boxcentre(:) = xtop(3*centre_atom-2:3*centre_atom)
	elseif(line == 'MASS') then   !define center by center of mass
		if (.not. get_centre_by_mass(boxcentre(:))) then
			write(*,*) ('>>>>> ERROR: Could not create centre ')
			return
		end if
	else !got x
		read(line, *, iostat=filestat) boxcentre(1)
		if(filestat > 0) then !invalid x coordinate
			return
		end if
		xwat_in=get_real_arg('-----> Box centre y: ')
		boxcentre(2) = xwat_in
		xwat_in=get_real_arg('-----> Box centre z: ')
		boxcentre(3) = xwat_in
	end if


	!read the size of the box
	coord_in = get_real_arg('-----> Boxlength x-direction: ')
	boxlength(1) = coord_in
	coord_in = get_real_arg('-----> Boxlength y-direction: ')
	boxlength(2) = coord_in
	coord_in = get_real_arg('-----> Boxlength z-direction: ')
	boxlength(3) = coord_in

	if( any(boxlength(:) == zero) ) then
		inv_boxl(:) = zero
	else
		inv_boxl(:) = one/boxlength(:)
	end if

	do i=1,3
		if( boxlength(i) < zero ) then
			write(*, '(a)') '>>>>> WARNING: Boxlength with negative sign. Converting to positive.'
			boxlength(i) = -boxlength(i)
		end if
	end do

	have_solvent_boundary = .true.
	set_solvent_box = .true.

end function set_solvent_box

!--------------------------------------------------------------------------------

subroutine solvate
	! Make sure boundary condition is set.
	if (.not. boundary_set) then
		write(*,'(a)') '>>>>> ERROR Boundary unknown'
		write(*,'(a)') "Use <boundary> to define boundary condition"
		return
	end if

	if (use_PBC) then
		call solvate_box
	else
		call solvate_sphere
	end if

end subroutine solvate

!-------------------------------------------------------------------------

subroutine solvate_box
	character(len=80)			::	solvate_mode

	write(*,'(a)') 'Using predefined boxcentre and boxlengths.'
	write(*,101) boxcentre
	write(*,102) boxlength
101 format('Boxcentre (x,y,z)                :  ',3f8.3)
102 format('Boxlengths (x,y,z)               :  ',3f8.3)

	!chose solvation mode
	call get_string_arg(solvate_mode, &
		'-----> Select solvation mode: grid(1), file(2), restart(3): ')
	call upcase(solvate_mode)

	select case(solvate_mode)
	case ('GRID', '1')
		call solvate_box_grid

	case ('FILE', '2')
		call solvate_box_file

	case('RESTART', '3')
		call solvate_restart

	case default
		write(*,'(a)') '>>>>> ERROR: Unknown mode of solvation.'
		call parse_reset
		return
	end select

end subroutine solvate_box

!----------------------------------------------------------------------------

subroutine solvate_box_grid

!solvate sphere using grid

	!locals

	real(kind=prec)						::	xmin, xmax, ymin, ymax, zmin, zmax
	real(kind=prec)						::	xgrid, ygrid, zgrid
	integer						::	max_wat !max number of molecules
	integer						::	waters_in_box
	real(kind=prec)						::	radius2, solvent_grid
	character(len=200)			::	solvent

	!set water residue name
	call get_string_arg(solvent, &
		'-----> Enter library entry name for water molecule: ')
	solvent_name = solvent(1:4)
	!get residue code for solvent
	if(.not. set_irc_solvent()) return

	!calc cubic grid spacing from density
	if(lib(irc_solvent)%density > 0.) then
		solvent_grid = lib(irc_solvent)%density**(-1/3.)
	else
		write(*,900) lib(irc_solvent)%nam
900		format('>>>>> ERROR: Density not set in library entry ',a)
		return
	end if


 	!Make sure boxsize is consistent with the periodic boundary condition
 	!changed to nearest integer (nint) from truncation (int)  /M.A.
    boxlength(1) = nint(boxlength(1)/solvent_grid)*solvent_grid
    boxlength(2) = nint(boxlength(2)/solvent_grid)*solvent_grid
    boxlength(3) = nint(boxlength(3)/solvent_grid)*solvent_grid

!Calculate max number of waters possible
	max_wat = lib(irc_solvent)%density * boxlength(1) * boxlength(2) * boxlength(3)

	allocate(xw(3,lib(irc_solvent)%nat,max_wat), keep(max_wat), stat=alloc_status)
	call check_alloc('water sphere co-ordinate array')


	xmin = boxcentre(1) - boxlength(1)/2.0_prec + solvent_grid/2.0_prec
	xmax = boxcentre(1) + boxlength(1)/2.0_prec - solvent_grid/2.0_prec
	ymin = boxcentre(2) - boxlength(2)/2.0_prec + solvent_grid/2.0_prec
	ymax = boxcentre(2) + boxlength(2)/2.0_prec - solvent_grid/2.0_prec
	zmin = boxcentre(3) - boxlength(3)/2.0_prec + solvent_grid/2.0_prec
	zmax = boxcentre(3) + boxlength(3)/2.0_prec - solvent_grid/2.0_prec

	write(*,100) boxlength(1), boxlength(2), boxlength(3), solvent_grid
100	format('New boxlength                     = ',3f8.2,' A',/ &
		   'Grid spacing                      = ',f10.2,' A ')

	!Fill box with water
 	!xmax+0.1 is needed for intel/windows. Otherwise the last loop step is skipped
	waters_in_box = 0
	xgrid = xmin
	do while (xgrid <= xmax + 0.1)
		ygrid = ymin
		do while (ygrid <= ymax + 0.1)
			zgrid = zmin
			do while (zgrid <= zmax + 0.1)
				waters_in_box = waters_in_box + 1
				xw(1,1,waters_in_box) = xgrid
				xw(2,1,waters_in_box) = ygrid
				xw(3,1,waters_in_box) = zgrid
				!all the molecules inside are inside
				keep(waters_in_box) = .true.
				zgrid = zgrid + solvent_grid
			end do
			ygrid = ygrid + solvent_grid
		end do
		xgrid = xgrid + solvent_grid
	end do

	call add_solvent_to_topology(waters_in_sphere=waters_in_box, &
		max_waters=waters_in_box, make_hydrogens=.true., pack=solvent_pack)
	deallocate(xw,keep)

end subroutine solvate_box_grid

!-----------------------------------------------------------------------
subroutine solvate_box_file

!local variables
	character(len=80)		::	xwat_file
	integer 				::	fstat
	character(len=80)		::	line
	real(kind=prec) 				::	boxl, waterbox_v, waterbox(1:3)
	character(len=6)		::	sphere
	logical					::	replicate
	integer					::	extension(1:3)
	real(kind=prec)					::	extensionbox_v, ext(3,3)
	integer					::	nw, nnw !water molecule counters
	integer					::	i, j, k !loop indecis
	integer					::	filestat !error variable
	character(len=3)		::	atomnames
	character(len=4)		::	resnam(3)
	integer					::	resno(3)
	integer					::	nbox !number of replicated boxes
	integer					::	nwat_allocate
	real(kind=prec)					::	xcm(3) !center of the waterbox
	real(kind=prec)					::	wshift(3) !distanse to move waters
	integer					::	nwat_keep !how many waters to keep
	real(kind=prec)					::	temp(3) !temporary coordinate

!get the name of the file and open the file in unit 13
	call get_string_arg(xwat_file, '-----> Solvent file name: ')
	open (unit=13, file=xwat_file, status='old', form='formatted', action='read', iostat=fstat)

	if( fstat /= 0 ) then
		write(*,'(a)') '>>>>> ERROR: Could not open water coordinate file.'
		call parse_reset
		return
	end if

!read the size of the waterbox
	read(13, '(a80)') line
	read(line, *, iostat=fstat) boxl
	if( fstat /= 0 ) then
		write(*, '(a)') '>>>>> ERROR: Size not specified in water file.'
		close(13)
		call parse_reset
		return
	end if
	if( boxl < 0 ) then
		write(*, '(a)' ) '>>>>> WARNING: Size with negative sign. Converting to positive.'
		boxl = - boxl
	end if

!check if the file contains a sphere of water instead of a box
	read(line, *, iostat=fstat) boxl, sphere
	call upcase(sphere)
	if(sphere == 'SPHERE' ) then
		write(*, '(a)') '>>>>> ERROR: This file containts a sphere of water. Use a file with a box of water instead.'
		close(13)
		call parse_reset
		return
	end if

!Compute waterbox volume and display the waterbox sidelength
	waterbox(:) = boxl
	waterbox_v = boxl**3
	write(*, '(a, f10.3)') 'Boxlength of solvent file 		=', waterbox(1)

!Determine residue name to use for solvent molecule
	read(13,1) atomnames(1:1), solvent_name
	backspace(13)
	if( .not. set_irc_solvent() ) then
		close(13)
		call parse_reset
		return
	else if( lib(irc_solvent)%nat /= 3 ) then
		write(*,'(a)') '>>>>> ERROR: Solvate only works for 3-atom solvents (in this version).'
		close(13)
		call parse_reset
		return
	else if( lib(irc_solvent)%density <= 0. ) then
		write(*, '(a, a)' ) '>>>>> ERROR: Density not set in library entry ', lib(irc_solvent)%nam
		close(13)
		call parse_reset
		return
	end if

!Estimate amount of memory to allocate for temporary waters
	if( all(boxlength(:)<waterbox(1)) ) then ! don't need to replicate. 5% margin
		replicate = .false.
		nwat_allocate = int( lib(irc_solvent)%density*1.05_prec*waterbox_v )

	else !the waterbox is not big enough
		replicate = .true.
		!find out in wich direction replication is needed
		extension(:) = ceiling( boxlength(:)/boxl )
		extensionbox_v = extension(1)*extension(2)*extension(3)*waterbox_v
		nwat_allocate = int( lib(irc_solvent)%density*1.05_prec*extensionbox_v )
	end if

	allocate( xw(3, lib(irc_solvent)%nat, nwat_allocate), keep(nwat_allocate), stat=alloc_status )
	call check_alloc('temporary solvent coord. arrays')

!The reading of the coordinates
1	format(13x, a1, 3x, a4, i5, 4x, 3f8.3)

	nw = 0

	do i = 1, nwat_allocate-1
		read(13, 1, iostat=filestat, end=10) &
			atomnames(1:1), resnam(1), resno(1), xw(:, 1, nw+1), &
			& atomnames(2:2), resnam(2), resno(2), xw(:, 2, nw+1), &
			& atomnames(3:3),	resnam(3), resno(3), xw(:, 3, nw+1)
		!checking the read info
		if(filestat > 0) then
			write(*, 7) nw
			close(13)
			call parse_reset
			deallocate(xw, keep)
			return
		else if ( any( resnam(:)/=solvent_name ) ) then
			write(*, 8) solvent_name, nw
			close(13)
			call parse_reset
			deallocate(xw, keep)
			return
		else if ( any( resno(:)/=resno(1) ) ) then
			write(*, 9) nw
			close(13)
			call parse_reset
			deallocate(xw, keep)
			return
		else
			nw = nw + 1
		end if
	end do

7	format('>>>>> ERROR: Read failure at molecule', i6)
8	format('>>>>> ERROR: Residue name other than ', a4, ' found at molecule ', i6)
9	format('>>>>> ERROR: Inconsistent residue numbering at molecule ', i6)

10	write(*, '(a, i10)') 'No. of molecules in solvent file ', nw

	nnw = nw
!Replicate if necessary
	if(replicate) then
		write(*, '(a)') 'Replicating box of water...'
		ext(:,1) = (/boxl, 0.0_prec, 0.0_prec/)
		ext(:,2) = (/0.0_prec, boxl, 0.0_prec/)
		ext(:,3) = (/0.0_prec, 0.0_prec, boxl/)

		nbox = 0
		do k=1,3 !the three coordinates
			if( .not. extension(k) > 1) cycle

			do i = 1,extension(k)-1
				nbox = nbox + 1
				do j = 1,nw
					nnw = nnw + 1
					xw(:,1,nnw) = xw(:,1,j) + i*ext(k,:)
					xw(:,2,nnw) = xw(:,2,j) + i*ext(k,:)
					xw(:,3,nnw) = xw(:,3,j) + i*ext(k,:)
				end do
				waterbox(k) = waterbox(k) + boxl
				end do
            nw = nnw
		end do

		nw = nnw
		write(*, '(a, i5, a)') 'Added', nbox, ' boxes of water'
		write(*,'(a, i10)') 'No of molecules after replication		= ', nw

	end if


!Compute the centre coordinates for the waterbox
	xcm(:) = sum( sum(xw, dim=3), dim=2) / (nw*lib(irc_solvent)%nat)
	wshift(:) = boxcentre(:) - xcm(:)

!Remove water molecules outside of periodic box
	nwat_keep = 0

	do i = 1,nw
		temp = xw(:,1,i)

		if(all(temp(:)<xcm(:)+0.5_prec*boxlength(:)) &
			.and. all(temp(:)>xcm(:)-0.5_prec*boxlength(:))) then
			!keep this molecule
			nwat_keep = nwat_keep + 1
			keep(i) = .true.
			!shift to periodic box centre
			do j = 1, lib(irc_solvent)%nat
				xw(:, j, i) = xw(:, j, i) + wshift(:)
			end do

		else !do not keep this molecule
			keep(i) = .false.

		end if

	end do

!Add solvent to topology
	call add_solvent_to_topology(waters_in_sphere=nwat_keep, max_waters=nw, &
		make_hydrogens=.false., pack=solvent_pack)

	close(13)
	deallocate(xw, keep)

end subroutine solvate_box_file

!---------------------------------------------------------------------------------------------

subroutine solvate_sphere

	character(len=80)			::	solvate_mode

	if(.not. set_solvent_sphere()) then
		call parse_reset !clear command line
		return
	end if

	!chose solvation mode
	call get_string_arg(solvate_mode, &
		'-----> Select solvation mode: grid(1), file(2), restart(3), DWF(4): ')
	call upcase(solvate_mode)

	select case(solvate_mode)
	case ('GRID','1')
		call solvate_sphere_grid

	case ('FILE','2')
		call solvate_sphere_file(.true.)

	case ('RESTART','3')
		call solvate_restart

	case ('DWF','4')
		call solvate_sphere_file(.false.)

	case default
		write(*,'(a)') '>>>>> ERROR: Unknown mode of solvation.'
		call parse_reset
		return
	end select

end subroutine solvate_sphere

!------------------------------------------------------------------------------------------------

logical function set_solvent_sphere()

	!locals
	character(len=80)			::	line
	integer						::	filestat
	integer						::	centre_atom
	real(kind=prec)						::	rwat_in, xwat_in

	set_solvent_sphere = .false.
	have_solvent_boundary = .false.

	if (have_title) then !Print old centre if available
		write(*,'(a,3f8.3)') 'Previous solvent centre:', xwcent(1), xwcent(2), xwcent(3)
	end if
	call get_string_arg(line, '-----> Sphere centre (<x y z> or <residue:atom_name> or <"mass"> or <"boundary">): ')
	call upcase(line)

	if(scan(line, ':') > 0) then !got res:at
		centre_atom=get_atom_from_descriptor(line)
		if(centre_atom == 0) then
			write(*,900) trim(line)
900			format('>>>>> ERROR: Could not find centre atom ',a)
			return
		end if
		xwcent(:) = xtop(3*centre_atom-2:3*centre_atom)
	elseif(line == 'MASS') then   !define center by center of mass
		if (.not. get_centre_by_mass(xwcent(:))) then
			write(*,*) ('>>>>> ERROR: Could not create centre ')
			return
		end if
	elseif(line == 'BOUNDARY') then   !define center same as boundary center
		xwcent = xpcent
	else !got x
		read(line, *, iostat=filestat) xwcent(1)
		if(filestat > 0) then !invalid x coordinate
			return
		end if
		xwat_in=get_real_arg('-----> Sphere centre y: ')
		xwcent(2) = xwat_in
		xwat_in=get_real_arg('-----> Sphere centre z: ')
		xwcent(3) = xwat_in
	end if

	if (have_title) then !Print old radius if available
		write(*,'(a,f8.3)') 'Previous effective solvent radius:', rwat_eff()
		rwat_in = get_real_arg('-----> Solvation sphere radius: ')
	else
		rwat_in = get_real_arg('-----> Solvation sphere radius: ')
	end if
	rwat = rwat_in
	if(rwat == 0.) then
		return
	end if

    write(*,*)
	write(*,100) xwcent(:)
100	format('Solvation sphere centre                   :   ',3f8.3)

	set_solvent_sphere = .true.
	have_solvent_boundary = .true.

end function set_solvent_sphere

!-----------------------------------------------------------------------

subroutine solvate_sphere_grid

!solvate sphere using grid

	!locals

	real(kind=prec)						::	xmin, xmax, ymin, ymax, zmin, zmax
	real(kind=prec)						::	xgrid, ygrid, zgrid
	integer						::	max_wat !max number of molecules
	integer						::	waters_in_sphere
	real(kind=prec)						::	radius2, solvent_grid
	character(len=200)			::	solvent

	!set water residue name
	call get_string_arg(solvent, &
		'-----> Enter library entry name for water molecule: ')
	solvent_name = solvent(1:4)
	!get residue code for solvent
	if(.not. set_irc_solvent()) return

	!calc cubic grid spacing from density
	if(lib(irc_solvent)%density > 0.) then
		solvent_grid = lib(irc_solvent)%density**(-1/3.)
	else
		write(*,900) lib(irc_solvent)%nam
900		format('>>>>> ERROR: Density not set in library entry ',a)
		return
	end if

	max_wat = (2.0_prec*rwat+solvent_grid)**3 / solvent_grid**3
	radius2 = rwat**2

	allocate(xw(3,lib(irc_solvent)%nat,max_wat), keep(max_wat), stat=alloc_status)
	xw = zero
	call check_alloc('water sphere co-ordinate array')

	xmin = xwcent(1) - int(rwat/solvent_grid)*solvent_grid
	xmax = xwcent(1) + int(rwat/solvent_grid)*solvent_grid
	ymin = xwcent(2) - int(rwat/solvent_grid)*solvent_grid
	ymax = xwcent(2) + int(rwat/solvent_grid)*solvent_grid
	zmin = xwcent(3) - int(rwat/solvent_grid)*solvent_grid
	zmax = xwcent(3) + int(rwat/solvent_grid)*solvent_grid

	write(*,100) rwat, solvent_grid
100	format('Radius of solvent sphere                = ',f10.2,' A',/ &
		   'Grid spacing                            = ',f10.2,' A ')

	!constuct water-only sphere
 	!xmax+0.1 is needed for intel/windows. Otherwise the last loop step is skipped
	waters_in_sphere = 0
	xgrid = xmin
	do while (xgrid <= xmax + 0.1)
		ygrid = ymin
		do while (ygrid <= ymax + 0.1)
			zgrid = zmin
			do while (zgrid <= zmax + 0.1)
				if( .not. ((xgrid-xwcent(1))**2 + (ygrid-xwcent(2))**2 &
					+ (zgrid-xwcent(3))**2 > radius2) ) then
				waters_in_sphere = waters_in_sphere + 1
				!if not outside keep these coordinates
				xw(1,1,waters_in_sphere) = xgrid
				xw(2,1,waters_in_sphere) = ygrid
				xw(3,1,waters_in_sphere) = zgrid
				!all the molecules inside are inside
				keep(waters_in_sphere) = .true.
				endif
				zgrid = zgrid + solvent_grid
			end do
			ygrid = ygrid + solvent_grid
		end do
		xgrid = xgrid + solvent_grid
	end do

	call add_solvent_to_topology(waters_in_sphere=waters_in_sphere, &
		max_waters=waters_in_sphere, make_hydrogens=.true., pack=solvent_pack)
	deallocate(xw,keep)

end subroutine solvate_sphere_grid

!-----------------------------------------------------------------------

subroutine solvate_sphere_file(shift)
! parameters
	logical, intent(in), optional::	shift

! local variables
	integer						::	i,j,nw,nnw
	integer						::	nwat_allocate, nwat_keep
	real(kind=prec)						::	rmax2,dx2,boxl, newboxl
	real(kind=prec), save				::	xcm(3),wshift(3)
	integer						::	fstat
	logical						::	is_box
	character(len=6)			::	sphere
	character(len=80)			::	line
	real(kind=prec)						::	volume
	real						::	r4dum
	character(len=80)			::	xwat_file
	real(kind=prec)						::	xwshift(3,7)
	integer						::	box
	character(len=3)			::	atomnames
	integer						::	filestat
	character(len=10)			::	filepos
	integer						::	resno(3)
	character(len=4)			::	resnam(3)

	real(kind=prec)						::	xwtmp(3, max_atlib)
	integer						::	at_id(max_atlib)

	call get_string_arg(xwat_file, '-----> Solvent file name: ')
	open (unit=13, file=xwat_file, status='old', form='formatted', &
		action='read', iostat=fstat)
	if(fstat /= 0) then
		write(*,'(a)') '>>>>> ERROR: Could not open water coordinate file.'
		call parse_reset
		return
	end if

2	format('Boxlength of solvent file               = ',f10.3)
3	format('Radius of sphere in solvent file        = ',f10.3)

  ! --> water coordinate file for initial generation (13)

	read (13,'(a80)') line
	read(line,*, iostat=fstat) boxl
	if(fstat /=0 ) then
		write(*,'(a)') 'ERROR >>>>>: Size not specified in water file.'
		close(13)
		call parse_reset
		return
	end if
	read(line, *, iostat=fstat) boxl, sphere
	call upcase(sphere)
	if(sphere == 'SPHERE') then
		volume = boxl**3*4.0_prec*pi/3.0_prec
		is_box = .false.
		write(*,3) boxl
	else
		volume = boxl**3
		is_box = .true.
		write(*,2) boxl
	end if

	!determine residue name to use for solvent molecules
	read(13,1) atomnames(1:1), solvent_name
	backspace(13)
	if(.not. set_irc_solvent()) then
		goto 999
	else if(lib(irc_solvent)%nat /= 3) then
		write(*,900)
900		format('>>>>> ERROR: Solvate only works for 3-atom solvents (in this version).')
		goto 999
	else if(lib(irc_solvent)%density <= 0.) then
		write(*,910) lib(irc_solvent)%nam
910		format('>>>>> ERROR: Density not set in library entry ',a)
		goto 999
	end if

  ! estimate amount of memory to allocate for temporary water coordinates
	if((is_box .and. abs(rwat) <= boxl/2.)) then
		!don't need to replicate
		!5% safety margin
		nwat_allocate = int(lib(irc_solvent)%density*1.05*volume)
	elseif(is_box) then
		newboxl = boxl*2**(int(log(abs(rwat)*2/boxl)/log(2.)+0.9999))
		nwat_allocate = int(lib(irc_solvent)%density*1.05*newboxl**3)
	elseif(abs(rwat) > boxl) then !it's a sphere and it's too small
		write(*,'(a)') '>>>>> ERROR: Water sphere in the file is too small.'
		call parse_reset
		close(13)
		return
	else !its a sufficiently big sphere
		nwat_allocate = int(lib(irc_solvent)%density*1.05*volume)
	endif
	allocate(xw(3,lib(irc_solvent)%nat,nwat_allocate), keep(nwat_allocate), &
       stat=alloc_status)
	call check_alloc('temporary solvent coord. arrays')

	rmax2 = rwat**2
	nw = 0
1	format(13x,a1,3x,a4,i5,4x,3f8.3)
	do i=1,nwat_allocate-1
		read (13,1,iostat=filestat, end=10) &
			atomnames(1:1), resnam(1), resno(1), xw(:, 1,nw+1), &
			atomnames(2:2), resnam(2), resno(2), xw(:, 2,nw+1), &
			atomnames(3:3), resnam(3), resno(3), xw(:, 3,nw+1)
		if(filestat > 0) then
			write(*,8) nw
			goto 999
		elseif(any(resnam(:) /= solvent_name)) then
			write(*,9) solvent_name, nw
			goto 999
		elseif(any(resno(:) /= resno(1))) then
			write(*,7) nw
			goto 999
		else
			nw = nw+1
		end if
	end do
7	format('>>>>> ERROR: Inconsistent residue numbering at molecule',i6)
8	format('>>>>> ERROR: Read failure at molecule',i6)
9	format('>>>>> ERROR: Residue name other than ',a4,&
		   ' found at molecule ',i6)
10	write (*,20) nw
20	format ('No. of molecules in solvent file           = ',i10)

  ! --- Replicate box if necessary

  nnw = nw
  if (is_box .and.  abs(rwat) .gt. boxl/2. ) then
     do while ( abs(rwat) .gt. boxl/2. )
        write (*,'(a)') 'Replicating periodic box...'

		xwshift(:,1) = (/boxl, 0.0_prec, 0.0_prec/)
		xwshift(:,2) = (/0.0_prec, boxl, 0.0_prec/)
		xwshift(:,3) = (/0.0_prec, 0.0_prec, boxl/)
		xwshift(:,4) = (/boxl, boxl, 0.0_prec/)
		xwshift(:,5) = (/boxl, 0.0_prec, boxl/)
		xwshift(:,6) = (/0.0_prec, boxl, boxl/)
		xwshift(:,7) = (/boxl, boxl, boxl/)

        do i=1,nw
			do box = 1,7
				nnw = nnw+1
				xw(:,1,nnw) = xw(:,1,i) + xwshift(:,box)
				xw(:,2,nnw) = xw(:,2,i) + xwshift(:,box)
				xw(:,3,nnw) = xw(:,3,i) + xwshift(:,box)
			end do
		end do
        boxl = 2.*boxl
        nw = nnw
        write (*,28) nw
28      format ('Number of molecules after replication      = ',i10)
     end do
	end if

	!no water coord shift
!	shift(:) = 0.
!	cm(:) = xwcent(:)

	! shift water coords from water center of mass to sphere center
	! do not shift for "dedicated water file"
	if(shift) then ! shift == TRUE => normal water file equilibrated without solute
		!xcm = center of mass for solvent atoms
		xcm(:) = sum(sum(xw, dim=3), dim=2)/(nw*lib(irc_solvent)%nat)
		wshift(:) = xwcent(:) -xcm(:) + epsilon(r4dum)

		nnw = 0
		nwat_keep = 0
		do i=1,nw
			dx2 = dot_product(xw(:,1,i)-xcm(:),xw(:,1,i)-xcm(:))
			if(dx2 <= rmax2) then
				nwat_keep = nwat_keep+1
				!shift to xwcent
				do j = 1, lib(irc_solvent)%nat
					xw(:,j,i) = xw(:,j,i) + wshift(:)
				end do
				keep(i) = .true. ! water i is in sphere
			else
				keep(i) = .false.
			end if
		end do
	else ! shift == FALSE => water file equilibrated with solute
	    do i=1,nw
			dx2 = dot_product(xw(:,1,i)-xwcent(:),xw(:,1,i)-xwcent(:))
			if(dx2 <= rmax2) then
				nwat_keep = nwat_keep+1
				keep(i) = .true. ! water i is in sphere
			else
				keep(i) = .false.
			end if
		end do
	endif



	call add_solvent_to_topology(waters_in_sphere=nwat_keep, &
		max_waters=nw, make_hydrogens=.false., pack=solvent_pack)
	!branch here on error
999	close(13)
	deallocate(xw,keep)

end subroutine solvate_sphere_file

!-------------------------------------------------------------------------

subroutine solvate_restart

!locals
	integer(4)					::	natom, nat3, waters_added
	character(len=80)			::	xfile
	integer						::	u, fstat
	real(kind=prec),allocatable			::	xtmp(:)
	character(len=200)			::	solvent

	u=freefile()

	call get_string_arg(xfile, '-----> Restart file name: ')
	open (unit=u, file=xfile, status='old', form='unformatted', &
		action='read', iostat=fstat)
	if(fstat /= 0) then
		write(*,'(a)') '>>>>> ERROR: Could not open restart file.'
		call parse_reset
		return
	end if

	!set water residue name
	call get_string_arg(solvent, &
		'-----> Enter library entry name for water molecule: ')
	solvent_name = solvent(1:4)
	!get residue code for solvent
	if(.not. set_irc_solvent()) return

	read(u) nat3
	natom = nat3/3
	if(natom == nat_pro) then
		write(*,10)
10		format('>>>>> ERROR: No. of atoms in restart file not greater than in topology.')
		close(u)
		return
	elseif(mod(natom-nat_pro,solv_atom) /= 0) then
		write(*,20) natom,solv_atom
20		format('>>>>> ERROR:',i6,' new atoms in restart file is not a multiple of',i6,'.')
		close(u)
		return
	end if


	waters_added = (natom-nat_pro)/solv_atom
	call grow_arrays_for_solvent(waters_added, solv_atom)
	allocate(xtmp(nat3))
	backspace(u)
	read(u) nat3, xtmp(1:nat3)
	close(u)
	allocate(xw(3,solv_atom,waters_added), keep(waters_added))
	xw(:,:,:) = reshape(xtmp(3*nat_pro+1:nat3),(/3,solv_atom,waters_added/))
	keep(:) = .true.

	!allow very tight packing of waters to solute - this is a restart file!
	call add_solvent_to_topology(waters_in_sphere=waters_added, &
		max_waters=waters_added, make_hydrogens=.false., pack=0.0_prec)
	deallocate(xtmp,xw,keep)

end subroutine solvate_restart

!-------------------------------------------------------------------------
subroutine add_solvent_to_topology(waters_in_sphere, max_waters, make_hydrogens, pack)
!arguments
	integer						::	waters_in_sphere
	integer						::	max_waters
	logical						::	make_hydrogens
	real(kind=prec)						::	pack
!locals
	integer						::	waters_added
	real(kind=prec)						::	rpack2
	integer						::	w_at, w_mol, p_atom
	logical						::	wheavy(max_atlib)
	real(kind=prec)						::	dx, dy, dz, r2
	integer						::	next_wat, next_atom, num_heavy, added_heavy
	integer						::	heavy_bonds, heavy_angles, added, i, j
	type(MISSING_TYPE),allocatable			::	missing_heavy(:)
	type(MISSING_BOND_TYPE),allocatable		::	missing_bonds(:)
	type(MISSING_ANGLE_TYPE),allocatable		::	missing_angles(:)
	type

	if(use_PBC) then
		write(*,111) waters_in_sphere
	else
		write(*,110) waters_in_sphere
	end if
110	format('No. of solvent molecules in the sphere:', i6)
111 format('No. of solvent molecules in the box:', i6)

	!exclude water molecues close to heavy solute atoms
	write(*,120) pack
120	format('Minimum distance to solute heavy atoms:',f6.2,' A')

	waters_added = 0
	rpack2 = pack**2
	num_heavy = 0
	do w_at = 1, lib(irc_solvent)%nat
		if(lib(irc_solvent)%atnam(w_at)(1:1) == 'H') then
			wheavy(w_at) = .false.
		else
			wheavy(w_at) = .true.
			num_heavy = num_heavy + 1
		end if
	end do

! before checking clashes we will now have to add all the missing heavy atoms
! so we do one loop over all missing heavy atoms
	if (num_heavy .gt. 1) then
! make the connectivity for the solvent atoms and temporary storage of variables
		heavy_bonds = 0
		allocate(missing_bonds(lib(irc_solvent)%nbnd))
		do i = 1, lib(irc_solvent)%nbnd
			if ((.not.wheavy(lib(irc_solvent)%bnd(i)%i)).or.(.not.wheavy(lib(irc_solvent)%bnd(i)%j))) cycle
			missing_bonds(i)%i=lib(irc_solvent)%bnd(i)%i
			missing_bonds(i)%j=lib(irc_solvent)%bnd(i)%j
			do j = 1, nbnd_types
				if((wildcard_tac(lib(irc_solvent)%atnam(missing_bonds(i)%i))==bnd_types(j)%taci .and. wildcard_tac(lib(irc_solvent)%atnam(missing_bonds(i)%j))==bnd_types(j)%tacj) &
					.or. (wildcard_tac(lib(irc_solvent)%atnam(missing_bonds(i)%j))==bnd_types(j)%taci .and. wildcard_tac(lib(irc_solvent)%atnam(missing_bonds(i)%i))==bnd_types(j)%tacj)) then
					missing_bonds(i)%cod = bnd_types(j)%cod
					return
				endif
			enddo
			heavy_bonds = heavy_bonds + 1
		end do
		allocate(missing_angles(heavy_bonds*4))
		heavy_angles = 0
		do i = 1, heavy_bonds - 1
			do j = i + 1, heavy_bonds
				if(missing_bonds(i)%i == missing_bonds(j)%i ) then
					heavy_angles = heavy_angles + 1
					missing_angles(heavy_angles)%i = missing_bonds(i)%j
					missing_angles(heavy_angles)%j = missing_bonds(i)%i
					missing_angles(heavy_angles)%k = missing_bonds(j)%j
				elseif(missing_bonds(i)%i == missing_bonds(j)%j ) then
					heavy_angles = heavy_angles + 1
					missing_angles(heavy_angles)%i = missing_bonds(i)%j
					missing_angles(heavy_angles)%j = missing_bonds(i)%i
					missing_angles(heavy_angles)%k = missing_bonds(j)%i
				elseif(missing_bonds(i)%j == missing_bonds(j)%i ) then
					heavy_angles = heavy_angles + 1
					missing_angles(heavy_angles)%i = missing_bonds(i)%i
					missing_angles(heavy_angles)%j = missing_bonds(i)%j
					missing_angles(heavy_angles)%k = missing_bonds(j)%j
				elseif(missing_bonds(i)%j == missing_bonds(j)%j ) then
					heavy_angles = heavy_angles + 1
					missing_angles(heavy_angles)%i = missing_bonds(i)%i
					missing_angles(heavy_angles)%j = missing_bonds(i)%j
					missing_angles(heavy_angles)%k = missing_bonds(j)%i
				endif
				do j = 1, nang_types
					if((wildcard_tac(lib(irc_solvent)%atnam(missing_angles(heavy_angles)%i))==ang_types(j)%taci .and. wildcard_tac(lib(irc_solvent)%atnam(missing_angles(heavy_angles)%j))==ang_types(j)%tacj &
						.and. wildcard_tac(lib(irc_solvent)%atnam(missing_angles(heavy_angles)%k))==bnd_types(j)%tack).or.(wildcard_tac(lib(irc_solvent)%atnam(missing_angles(heavy_angles)%i))==ang_types(j)%tack &
						.and. wildcard_tac(lib(irc_solvent)%atnam(missing_angles(i)%j))==ang_types(j)%tacj .and. wildcard_tac(lib(irc_solvent)%atnam(missing_angles(i)%k))==bnd_types(j)%taci)) then
						missing_angles(heavy_angles)%cod = ang_types(j)%cod
						return
					end if
				end do
			end do
		end do
! now make data structure that has the information dependent on the individual atoms that are missing
! remember, first atom is the one generated before
		allocate(missing_heavy(lib(irc_solvent)%nat))
		do i = 1, lib(irc_solvent)%nat
			missing_heavy(i)%atom_missing = .true.
			missing_heavy(i)%bond_notset(:)=.false.
			missing_heavy(i)%angle_notset(:)=.false.
			do j = 1, heavy_bonds
				added = 0
				if((i.eq.missing_bonds(j)%i).or.(i.eq.missing_bonds(j)%j)) then
					added = added + 1
					missing_heavy(i)%bonds(added) = j
					missing_heavy(i)%bond_notset(added) = .true.
				end if
			end do
			do j = 1, heavy_angles
				added = 0
				if((i.eq.missing_angles(j)%i).or.(i.eq.missing_angles(j)%j).or.(i.eq.missing_angles(j)%k)) then
					added = added + 1
					missing_heavy(i)%angles(added) = j
					missing_heavy(i)%angle_notset(added) = .true.
				end if
			end do
		end do
		missing_heavy(1)%atom_missing = .false.

		added_heavy = 0
		do w_mol = 1, max_waters
! Q always generates the first solvent atom on a grid, so we start generating missing heavy atoms after that
! means first atom has to be heavy center atom
! we stole the genH function for that
			added_heavy = added_heavy + genHeavy(w_mol,missing_heavy,missing_bonds,missing_angles)
		end do
166	format('Number of heavy solvent atoms added:',i6)
		write(*,166) added_heavy
		deallocate(missing_bonds,missing_angles,missing_heavy)
	end if


	!loop over solvent molecules
wloop:do w_mol = 1, max_waters
		!skip if not in sphere or box
		if(.not. keep(w_mol)) cycle wloop

		!loop over atoms in solvent molecule
		do w_at = 1, lib(irc_solvent)%nat
			!skip hydrogens
			if(.not. wheavy(w_at)) cycle
			!check clash with all other heavy atoms
ploop:		do p_atom = 1, nat_pro !for each water check all other atoms
				if(heavy(p_atom)) then
				!clash with a heavy atom?
				!*****PWchanged way of computing packing
					dx = xtop(3*p_atom-2) - xw(1,w_at,w_mol)
					dy = xtop(3*p_atom-1) - xw(2,w_at,w_mol)
					dz = xtop(3*p_atom  ) - xw(3,w_at,w_mol)
					if( use_PBC ) then
						dx = dx - boxlength(1)*nint( dx*inv_boxl(1) )
						dy = dy - boxlength(2)*nint( dy*inv_boxl(2) )
						dz = dz - boxlength(3)*nint( dz*inv_boxl(3) )
					end if
					r2 = dx**2 + dy**2 + dz**2
					if( r2<rpack2 ) then

					!if(dot_product(xtop(3*p_atom-2:3*p_atom) &
					!	-xw(:,w_at,w_mol), &
					!	xtop(3*p_atom-2:3*p_atom)-xw(:,w_at,w_mol)) &
					!	< rpack2) then
						keep(w_mol) = .false.
						cycle wloop
					end if
				end if
			end do ploop

		!*****PWadded new loop checking water-water distance
			if(use_PBC) then
				do next_wat = w_mol+1, max_waters
					if(.not. keep(next_wat) ) cycle
					do next_atom = 1,lib(irc_solvent)%nat
						if(.not. wheavy(next_atom) ) cycle
						dx = xw(1, w_at, w_mol) - xw(1, next_atom, next_wat)
						dy = xw(2, w_at, w_mol) - xw(2, next_atom, next_wat)
						dz = xw(3, w_at, w_mol) - xw(3, next_atom, next_wat)
						dx = dx - boxlength(1)*nint( dx*inv_boxl(1) )
						dy = dy - boxlength(2)*nint( dy*inv_boxl(2) )
						dz = dz - boxlength(3)*nint( dz*inv_boxl(3) )
						r2 = dx**2 + dy**2 + dz**2
						if( r2<rpack2 ) then
							keep(w_mol) = .false.
							cycle wloop
						end if
					end do
				end do
			end if

		end do !w_at
		if(keep(w_mol)) then
			!this water molecules does not clash with any other heavy atom.
			waters_added = waters_added + 1
		end if
	end do wloop
	write(*,130) waters_added
130	format('No. of solvent molecules to add to the system:',i6)

	!increase size of xtop, imakeh
	call grow_arrays_for_solvent(waters_added, lib(irc_solvent)%nat)

	do w_mol = 1, max_waters
		if(.not. keep(w_mol)) cycle
		!move this vater to topology
		nwat=nwat+1
		nmol = nmol + 1
		istart_mol(nmol) = nat_pro + 1
		nres = nres + 1
		res(nres)%start = nat_pro + 1
		res(nres)%irc = irc_solvent
		res(nres)%name = lib(irc_solvent)%nam
		do w_at = 1, lib(irc_solvent)%nat
			nat_pro = nat_pro + 1
			!copy coordinates for each atom
			xtop(3*nat_pro-2:3*nat_pro) = xw(:,w_at,w_mol)
			heavy(nat_pro) = wheavy(w_at)
			if(wheavy(w_at)) then
				makeH(nat_pro) = .false.
			else
				makeH(nat_pro) = make_hydrogens
			end if
		end do
	end do
end subroutine add_solvent_to_topology

!-----------------------------------------------------------------------

subroutine grow_arrays_for_solvent(nmore, atoms_per_molecule)
!make space for nmore more water molecules
!arguments
	integer, intent(in)			::	nmore, atoms_per_molecule
!locals
	real(kind=prec), allocatable		::	r8temp(:)
	logical, allocatable		::	ltemp(:)
	integer, allocatable		::	itemp(:)
	integer						::	new_nat, nat3old
	type(residue_type), allocatable	:: restemp(:)

	if(nmore == 0) return

	new_nat = nat_pro + nmore*atoms_per_molecule
	nat3old = nat_pro*3
	if(allocated(xtop)) then
		allocate(r8temp(nat3old), stat=alloc_status)
		call check_alloc('reallocating topology atom array')
		r8temp(1:nat3old) = xtop(1:nat3old)
		deallocate(xtop)
		allocate(xtop(new_nat*3),  stat=alloc_status)
		call check_alloc('reallocating topology atom array')
		xtop(1:nat3old) = r8temp(1:nat3old)
		deallocate(r8temp)
	else
		allocate(xtop(new_nat*3),  stat=alloc_status)
		call check_alloc('reallocating topology atom array')
	end if

	if(allocated(makeH)) then
		allocate(ltemp(nat_pro), stat=alloc_status)
		call check_alloc('reallocating topology atom array')
		ltemp(1:nat_pro) = makeH(1:nat_pro)
		deallocate(makeH)
		allocate(makeH(new_nat), stat=alloc_status)
		call check_alloc('reallocating topology atom array')
		makeH(1:nat_pro) = ltemp(1:nat_pro)
		makeH(nat_pro+1:new_nat) = .false.
		ltemp(1:nat_pro) = heavy(1:nat_pro)
		deallocate(heavy)
		allocate(heavy(new_nat), stat=alloc_status)
		call check_alloc('reallocating topology atom array')
		heavy(1:nat_pro) = ltemp(1:nat_pro)
		heavy(nat_pro+1:new_nat) = .false.
		deallocate(ltemp)
	else
		allocate(makeH(new_nat), heavy(new_nat), stat=alloc_status)
		call check_alloc('reallocating topology atom array')
		makeH(1:new_nat) = .false.
		heavy(1:new_nat) = .false.
	end if

	if(allocated(res)) then
		allocate(restemp(nres), stat=alloc_status)
		call check_alloc('reallocating topology atom array')
		restemp(1:nres) = res(1:nres)
		deallocate(res)
		allocate(res(nres + nmore), stat=alloc_status)
		call check_alloc('reallocating topology atom array')
		res(1:nres) = restemp(1:nres)
		deallocate(restemp)
	else
		allocate(res(nres + nmore), stat=alloc_status)
		call check_alloc('reallocating topology atom array')
	end if

	if(allocated(istart_mol)) then
		allocate(itemp(nmol), stat=alloc_status)
		call check_alloc('reallocating topology atom array')
		itemp(1:nmol) = istart_mol(1:nmol)
		deallocate(istart_mol)
		allocate(istart_mol(nmol+nmore), stat=alloc_status)
		call check_alloc('reallocating topology atom array')
		istart_mol(1:nmol) = itemp(1:nmol)
		deallocate(itemp)
	else
		allocate(istart_mol(nmol+nmore), stat=alloc_status)
		call check_alloc('reallocating topology atom array')
	end if
end subroutine grow_arrays_for_solvent

!-----------------------------------------------------------------------

real(kind=prec) function rwat_eff()
  ! local variables
  integer						:: i,kr,isort,bins
  real(kind=prec)						:: rc,rnwat
  integer, allocatable			:: npro_of_r(:)
	real(kind=prec)						::	rho_ratio, solvent_volume, rho_solvent

	rwat_eff = zero
	if ( nwat .eq. 0 ) return

	!make sure to include all waters in topology
  bins = int(max(abs(rwat), rexcl_o)*100)+1000
  allocate(npro_of_r(bins), stat=alloc_status)
  call check_alloc('protein radial distribution array')


	npro_of_r(:) = 0
	do i = 1, nat_solute
		rc    = sqrt(dot_product(xtop(3*i-2:3*i)-xwcent(:),&
			xtop(3*i-2:3*i)-xwcent(:)))
		isort = int ( 100. * rc )
		if ( isort .eq. 0 ) isort = 1
		if ( isort .le. bins ) then
			if(heavy(i) .and. (.not. excl(i)))  &
				npro_of_r(isort) = npro_of_r(isort) + 1
		end if
	end do

	!calc solvent density
	solvent_volume = 0.
	do i = nres_solute + 1, nres
		if(lib(res(i)%irc)%density > 0) then
			solvent_volume = solvent_volume + 1./lib(res(i)%irc)%density
		end if
	end do
	rho_solvent = nwat/solvent_volume
	! --- rho_ratio (0.577) comes from dividing the average volume of
	!     protein atoms (17.26 A**3) by the volume of a water molecule (29.9 A**3).
	rho_ratio = rho_solvent/rho_solute
	rnwat = zero
	do kr = 1, bins
		rc = real(kr,kind=prec)*0.01_prec
		rnwat = rnwat + 4.0_prec*pi*rc*rc*rho_solvent &
			*0.01_prec-rho_ratio*npro_of_r(kr)
		if ( int(rnwat) >= nwat - nexwat ) exit
	end do
	rwat_eff = rc

	deallocate(npro_of_r)
end function rwat_eff

!-----------------------------------------------------------------------

type(TOR_CODES) function torcode(taci, tacj, tack, tacl)
!arguments
	character(KEYLENGTH)			::	taci, tacj, tack, tacl
! *** local variables
	integer i
	character(KEYLENGTH)			::	ti, tj, tk, tl

	ti = wildcard_tac(taci)
	tj = wildcard_tac(tacj)
	tk = wildcard_tac(tack)
	tl = wildcard_tac(tacl)

	torcode%ncod = 0
	torcode%cod(:) = 0

	!search for exact match first
	do i = 1, ntor_types
		if((ti == tor_types(i)%taci .and. tj == tor_types(i)%tacj &
			.and. tk == tor_types(i)%tack .and. tl == tor_types(i)%tacl) &
			.or. &
			(ti == tor_types(i)%tacl .and. tj == tor_types(i)%tack .and. tk == tor_types(i)%tacj &
			.and. tl == tor_types(i)%taci)) then

			torcode%ncod = torcode%ncod + 1
			torcode%cod(torcode%ncod) = tor_types(i)%cod
		endif
	enddo
	if(torcode%ncod > 0) return !stop looking if something was found

	!then search with wildcard at i
	do i = 1, ntor_types
		if(tor_types(i)%taci == '') then
		    if((tj == tor_types(i)%tacj .and. tk == tor_types(i)%tack .and. tl == tor_types(i)%tacl) .or. &
			    (ti == tor_types(i)%tacl .and. tj == tor_types(i)%tack .and. tk == tor_types(i)%tacj)) then

				torcode%ncod = torcode%ncod + 1
				torcode%cod(torcode%ncod) = tor_types(i)%cod
		    endif
		endif
	enddo
	if(torcode%ncod > 0) return !stop looking if something was found

	!then search with wildcard at l
	do i = 1, ntor_types
		if(tor_types(i)%tacl == '') then
		    if((ti == tor_types(i)%taci .and. tj == tor_types(i)%tacj .and. tk == tor_types(i)%tack) .or. &
			    (tj == tor_types(i)%tack .and. tk == tor_types(i)%tacj .and. tl == tor_types(i)%taci)) then

				torcode%ncod = torcode%ncod + 1
				torcode%cod(torcode%ncod) = tor_types(i)%cod
		    endif
		endif
	enddo
	if(torcode%ncod > 0) return !stop looking if something was found

	!then search with wildcards at i & l
	do i = 1, ntor_types
		if(tor_types(i)%taci == '' .and. tor_types(i)%tacl == '') then
		    if((tj == tor_types(i)%tacj .and. tk == tor_types(i)%tack) .or. &
			    (tj == tor_types(i)%tack .and. tk == tor_types(i)%tacj)) then

				torcode%ncod = torcode%ncod + 1
				torcode%cod(torcode%ncod) = tor_types(i)%cod
		    endif
		endif
	enddo
	if(torcode%ncod > 0) return !stop looking if something was found

	!then search with wildcards at i & j or k & l
	do i = 1, ntor_types
		if(tor_types(i)%taci == '' .and. tor_types(i)%tacj == '') then
		    if((tk == tor_types(i)%tack .and. tl == tor_types(i)%tacl) .or. &
			    (ti == tor_types(i)%tacl .and. tj == tor_types(i)%tack)) then

				torcode%ncod = torcode%ncod + 1
				torcode%cod(torcode%ncod) = tor_types(i)%cod
		    endif
		ELSEif(tor_types(i)%tack == '' .and. tor_types(i)%tacl == '') then
		    if((ti == tor_types(i)%taci .and. tj == tor_types(i)%tacj) .or. &
			    (tl == tor_types(i)%taci .and. tk == tor_types(i)%tacj)) then

				torcode%ncod = torcode%ncod + 1
				torcode%cod(torcode%ncod) = tor_types(i)%cod
		    endif
		endif
	enddo
	if(torcode%ncod > 0) return !stop looking if something was found

	write( * , '(a,4(1x,a8))') '>>> Missing torsion type for atom types', &
		 ti, tj, tk, tl
	topo_ok = .false.

!.......................................................................
end function torcode

!-----------------------------------------------------------------------

subroutine tors_ene(emax, nlarge, av_ene)
! *** local variables
	integer i, j, k, l, ip, ic, i3, j3, k3, l3, nlarge
	real(kind=prec) rji(3), rjk(3), rkl(3), rnj(3), rnk(3), bj, bk, scp,	&
	phi, sgn, pe, dv, rki(3), rlj(3), dp(12), arg, f1, di(3),	 &
	dl(3)
	real(kind=prec) emax, av_ene

!.......................................................................

	nlarge = 0
	av_ene = zero

	do ip = 1, ntors

	i = tor(ip)%i
	j = tor(ip)%j
	k = tor(ip)%k
	l = tor(ip)%l
	ic = tor(ip)%cod
	if(ic == 0) then !missing parameter
		write(*, '(5i5,1x,a)') ip, i, j, k, l, 'MISSING PARAMETERS'
		cycle
	end if
	i3 = i * 3 - 3
	j3 = j * 3 - 3
	k3 = k * 3 - 3
	l3 = l * 3 - 3
	rji(1) = xtop(i3 + 1) - xtop(j3 + 1)
	rji(2) = xtop(i3 + 2) - xtop(j3 + 2)
	rji(3) = xtop(i3 + 3) - xtop(j3 + 3)
	rjk(1) = xtop(k3 + 1) - xtop(j3 + 1)
	rjk(2) = xtop(k3 + 2) - xtop(j3 + 2)
	rjk(3) = xtop(k3 + 3) - xtop(j3 + 3)
	rkl(1) = xtop(l3 + 1) - xtop(k3 + 1)
	rkl(2) = xtop(l3 + 2) - xtop(k3 + 2)
	rkl(3) = xtop(l3 + 3) - xtop(k3 + 3)
	rnj(1) = rji(2) * rjk(3) - rji(3) * rjk(2)
	rnj(2) = rji(3) * rjk(1) - rji(1) * rjk(3)
	rnj(3) = rji(1) * rjk(2) - rji(2) * rjk(1)
	rnk(1) = - rjk(2) * rkl(3) + rjk(3) * rkl(2)
	rnk(2) = - rjk(3) * rkl(1) + rjk(1) * rkl(3)
	rnk(3) = - rjk(1) * rkl(2) + rjk(2) * rkl(1)
	bj = sqrt(rnj(1) **2 + rnj(2) **2 + rnj(3) **2)
	bk = sqrt(rnk(1) **2 + rnk(2) **2 + rnk(3) **2)
	scp =(rnj(1) * rnk(1) + rnj(2) * rnk(2) + rnj(3) * rnk(3) )&
	/(bj * bk)
	if(scp>one) scp = one
	if(scp< -one) scp = -one
	phi = acos(scp)
	sgn = rjk(1) *(rnj(2) * rnk(3) - rnj(3) * rnk(2) ) + rjk(2)&
	*(rnj(3) * rnk(1) - rnj(1) * rnk(3) ) + rjk(3) *(rnj(1)   &
	* rnk(2) - rnj(2) * rnk(1) )
	if(sgn<zero) phi = - phi

! ---	energy

	arg = torlib(ic)%rmult * phi - torlib(ic)%deltor * pi / 180.0_prec
	pe = torlib(ic)%fk *(one + cos(arg) ) / real(torlib(ic)%paths )

	av_ene = av_ene+pe

	if(pe>emax) then
		nlarge = nlarge+1
		write( * , '(6i5,5f8.2)') ip, i, j, k, l, ic, torlib(ic)%fk ,   &
		torlib(ic)%rmult , torlib(ic)%deltor , phi * 180. / pi, pe
	endif

	enddo

	if(ntors/=0) av_ene = av_ene / real(ntors,kind=prec)

end subroutine tors_ene

!-----------------------------------------------------------------------

logical function check_residues()
!locals
	integer						::	i,j

	!stop immediately if no library
	if(nlibres == 0 .and. nat_pro > 0) then
		check_residues = .false.
		write(*, '(a)') '>>>>>ERROR: No library loaded.'
		return
	end if

	check_residues = .true.
	!check that all residues are in the library
resloop:do i = 1, nres
libloop:do j = 1, nlibres
			if(res(i)%name ==lib(j)%nam ) then
				res(i)%irc = j
				cycle resloop
			endif
		enddo libloop
		!if we reach here then the library entry was not found
		write(*,'(a,i3,a)') '>>> Residue ',i, ' not found in library: ',	&
			res(i)%name
		check_residues = .false.
	enddo resloop
end function check_residues

!-----------------------------------------------------------------------

subroutine writepdb
! *** local variables
	CHARACTER(len=80)			::	filnam
	character(len=1)			::	reply
	character(*), parameter	::	gap = &
		'                 GAP                                      '
	character(*), parameter	::	ter = 'TER   '
	integer						::	i, j, k, l, imol
	integer						::	iat
	integer						::	iwrite_g
	logical						::	wrote_atom_in_molecule
	character(len=6)			::	PDBtype
	integer						::	lig(4) !atoms connected to atom i

	if(.not. check_residues()) then
		call parse_reset()
		return
	end if
	if(mask%included == 0) then
		write(*,900)
900		format('>>>>> ERROR: The mask is empty - no atoms to write!')
		return
	end if

	write( *, * )
	call get_string_arg(filnam, '-----> Name of PBD file: ')
	if(openit(3, filnam, 'unknown', 'formatted', 'write') /= 0) then
		call parse_reset
		return
	end if
	REWIND(3)


	CALL get_string_arg(reply, '-----> Write out TER cards [y/n] ? ')
	select case(reply)
	case('y', 'Y')
		iwrite_g = 1
	case default
		iwrite_g = 0
	end select

!PDB standard minus B-factors 10	format(a6,i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3)
!See readpdb for specification of the format
!Old format  10	format(a6,i5,2x,a4,a4,i5,4x,3f8.3)
!TODO: Fix writing of chainID, requires chain info in topology.
   !(PDBtype,atomNr,atomName,resName,resNr,coords)
10	format(a6,i5,1x,a5,a4,1x,i4,4x,3f8.3)
11	format(a6,11x,a4,1x,i4) !For TER cards
	iat = 0
	imol = 1
	wrote_atom_in_molecule = .false.
	do i = 1, nres
		if(lib(res(i)%irc )%HETATM) then
			PDBtype = 'HETATM'
		else
			PDBtype = 'ATOM  '
		end if
		do j = 1, lib(res(i)%irc)%nat
			iat = iat + 1
			!write only atoms in mask
			if(mask%mask(iat)) then
				write(3, 10) PDBtype, iat, lib(res(i)%irc )%atnam(j), res(i)%name,&
					i,xtop(iat*3-2:iat*3)
				wrote_atom_in_molecule = .true.
			end if
		enddo
		if(i < nres_solute .and. imol < nmol) then
			if(istart_mol(imol + 1) ==iat + 1) then
				imol = imol + 1
				if(iwrite_g>0 .and. wrote_atom_in_molecule) then
					write(3, 11)  ter, res(i)%name, i
				end if
				wrote_atom_in_molecule = .false.
			endif
		end if
	enddo

	!now write connect records for HETATM groups
	do i = 1, nres_solute !loop over all solute residues
		if(lib(res(i)%irc)%HETATM) then !only for HETATM groups
			j = 1
			!work through all bonds
			do while(j <= lib(res(i)%irc)%nbnd)
				!this is the first atom
				iat = res(i)%start - 1 + lib(res(i)%irc)%bnd(j)%i
				!find a set of up to four bonds from this atom
				l = 0
				do k = 1,4
					lig(k) = res(i)%start - 1 + lib(res(i)%irc)%bnd(j)%j
					j=j+1
					!include only ligand atoms which are in the mask
					if(mask%mask(lig(k))) l = l + 1
					if(j>lib(res(i)%irc)%nbnd) exit
					!stop if new bond has different i atom
					if(res(i)%start - 1 + lib(res(i)%irc)%bnd(j)%i /= iat) exit
				end do
				if(mask%mask(iat) .and. l>0) write(3,20) iat, lig(1:l)
			end do
		end if
	end do
20	format('CONECT',5i5)

	write( * , '(/,a,/)') 'PDB file successfully written.'

	close(3)

end subroutine writepdb

!-----------------------------------------------------------------------

subroutine writemol2
! *** local variables
	CHARACTER filnam*80, reply*1
	character					::	mol_name*16, ti*1, tj*1
	integer						::	i, cnt, j, iat
	integer						::	u
	integer						::	mol, mol_start, mol_end, nat_mol, nbnd_mol
	integer						::	at_start, at_end, res_start, res_end
	integer						::	bnd_start, bnd_end, new_last, new_res
	logical						::	iwrite_h, iwrite_w
	integer						::	nres_mol,res_atoms(nres), new_num(nat_pro)
	integer						::	new_resnum(nres)
	integer						::	dict_type

	if(.not. check_residues()) then
		call parse_reset()
		return
	end if

	if(mask%included == 0) then
		write(*,900)
900		format('>>>>> ERROR: The mask is empty - no atoms to write!')
		return
	end if

	call get_string_arg(filnam, '-----> Name of Mol2 file (or auto): ')

5	format('#       Mol2 file written by Qprep version 5',/,&
		   '#       Title of topology: ',a,/,&
		   '#       Co-ordinate source: ',a)
! atom records
10	format(i6,1x,a4,4x,3f10.3,2x,a,t55,i6,2x,a4,2x,f6.3) !atom record
!              #    nam    xyz     atyp  res# rnam    q
!bond record
20	format(i6,1x,i6,1x,i6,4x,a3)
	!          #     i      j    typ
!molecule record : name,atoms, bonds, substr., features, sets, mol_type, charge_type
30	format(/,'@<TRIPOS>MOLECULE',/,'MOLECULE',/,5(i6,1x),/,a,/,a,/)
!substructure record
40	format(i6,1x,a3,a,t18,i6,1x,a,t34,i1,1x,a1,1x,a4) !substructure record
!           #    name#    start kind  dict chain, restyp

	!get new atom numbers
	new_num(:) = 0 !reset
	j = 0
	do i = 1, nat_pro
		if(mask%mask(i)) then
			j = j + 1
			new_num(i) = j
		end if
	end do

	!count atoms in residues
	do i = 1, nres-1
		res_atoms(i) = count(mask%mask(res(i)%start:res(i+1)%start-1))
	end do
	res_atoms(nres) = count(mask%mask(res(nres)%start:nat_pro))

	!renumber residues
	nres_mol = 0
	do i = 1, nres
		if(res_atoms(i) > 0) then
			nres_mol = nres_mol + 1
			new_resnum(i) = nres_mol
		end if
	end do

	!count bonds
	nbnd_mol = 0
	do i = 1, nbonds
		if(mask%mask(bnd(i)%i) .and. mask%mask(bnd(i)%j)) then
			!both atoms are in mask
			nbnd_mol = nbnd_mol + 1
		end if
	end do


	!name & open file
	if(filnam == 'auto' .or. filnam == 'AUTO') then !generate name
		write(filnam, 110) trim(auto_name), trim(adjustl(mol_name))
	end if
110 format(a,a,'.mol2')
	u = freefile()
	if(openit(u, filnam, 'unknown', 'formatted', 'write') /= 0) then
		call parse_reset()
		return
	end if

	!write header
	write(u,5) title, coord_source

	!write molecule record
	write(u,30) mask%included, nbnd_mol, nres_mol, 0, 0,'SMALL', 'USER_CHARGES'

	iat = 0
	new_res = 0
	write(u,'(a)') '@<TRIPOS>ATOM'
	do i = 1, nres
		do j = 1, lib(res(i)%irc)%nat
			iat = iat + 1
			if(mask%mask(iat)) then
				write(u, 10) new_num(iat), lib(res(i)%irc )%atnam(j), &
					xtop(iat*3-2:iat*3), SYBYL_atom_type(iac(iat)), &
					new_resnum(i), res(i)%name, crg(iat)
			end if
		enddo
	enddo

	!Write bond records
	write(u,'(a)') '@<TRIPOS>BOND'
	iat = 0
	do i = 1, nbonds
		if(mask%mask(bnd(i)%i) .and. mask%mask(bnd(i)%j)) then
			!both atoms are in mask
			iat = iat + 1
			if(SYBYL_bond_type(bnd(i)%cod) == '') then
				write(u,20) iat, new_num(bnd(i)%i), new_num(bnd(i)%j),"1"
			else
				write(u,20) iat, new_num(bnd(i)%i), new_num(bnd(i)%j), &
				 SYBYL_bond_type(bnd(i)%cod)
			end if
		end if
	end do

	!write substructure records
	write(u,'(a)') '@<TRIPOS>SUBSTRUCTURE'
	new_res = 0
	do i = 1, nres
		if(res_atoms(i) > 0) then
			new_res = new_res + 1
			at_start = res(i)%start
			!if no hydrogens then skip forward to first heavy atom
			do while(.not. mask%mask(at_start))
				at_start = at_start + 1
			end do
			dict_type = 0
			if(lib(res(i)%irc)%SYBYLTYPE == 'RESIDUE') dict_type = 1
			write(mol_name, '(i6)') i
			write(u,40) new_resnum(i), res(i)%name,adjustl(mol_name), new_num(at_start),&
				lib(res(i)%irc)%SYBYLTYPE, dict_type, 'A', res(i)%name
		end if
	end do

	close(u)

	write( * , '(/,a,/)') 'Mol2 file successfully written.'


end subroutine writemol2

!-----------------------------------------------------------------------

subroutine writetop
! *** local variables
	CHARACTER filnam * 80
	integer i, j, ig
	CHARACTER answer * 10

! --- Warn if missing parameters were found by maketop
	if(.not. topo_ok) then
		write( * ,  * ) 'WARNING: The topology is incomplete due to missing parameters!'
		write( * ,  * ) 'Do you realLY want to write this erronenous topology?'
		write( * ,  * ) 'Enter yes to proceed, anything else to cancel.'
		call parse_reset
		CALL get_string_arg(answer, '-----> Write _ERRONENOUS_ topology [yes/NO] ? ')
		if(answer/='yes') then
			return
		endif
	endif

	write( *, * )

	CALL get_string_arg(filnam, '-----> Give name of new topology file: ')

	call topo_save(filnam)

end subroutine writetop

!-----------------------------------------------------------------------

subroutine listprefs
	call pref_list('Preferences (use set command to change):')
end subroutine listprefs

!-----------------------------------------------------------------------

subroutine set
	!locals
	character(len=PREF_LEN)		::	name, value
	logical						::	l

	call get_string_arg(name, '-----> Enter name (or number): ')
	call get_string_arg(value, '-----> Enter value: ')

	l = pref_set(name, value)
end subroutine set

!*************************************************************************
!*Sort out atoms in restrained shell. Use protein center to calculate distance.
!*Use coordinates from topology unless 'implicit_rstr_from_file' is specified
!*Swiped from md.f90, needed when using atom masks "not restrained"
!*************************************************************************
subroutine make_shell2
! *** Local variables
	integer						::	i,ig,i3,k
	real(kind=prec)						::	rout2,rin2,r2
	real(kind=prec), allocatable		::	cgp_cent(:,:)
	nshellats = 0

	rexcl_i = get_real_arg('-----> Inner radius of restrained shell ')

	rin2  = rexcl_i**2

	shell(:) = .false.

	allocate(cgp_cent(3,ncgp+nwat))

	cgp_cent(:,:) = zero

	do ig=1,ncgp_solute
    if (.not. excl(cgp(ig)%iswitch)) then
                do i = cgp(ig)%first,cgp(ig)%last
			i3 = cgpatom(i)*3
			cgp_cent(:,ig) = cgp_cent(:,ig) + xtop(i3-2:i3)
		end do
        cgp_cent(:,ig) = cgp_cent(:,ig)/real(cgp(ig)%last - cgp(ig)%first +1)
		r2 = dot_product(cgp_cent(:,ig)-xpcent(:),cgp_cent(:,ig)-xpcent(:))

		if ( r2 .gt. rin2 ) then
			do i=cgp(ig)%first, cgp(ig)%last
				nshellats = nshellats + 1
				shell(cgpatom(i)) = .true.
			end do
		end if
   end if
	end do

	deallocate(cgp_cent)
	write(*,105) nshellats, rexcl_i, rexcl_o
105	format('Found   ',i6,' solute atoms in the restrained shell region (',f6.2,' to ',f6.2,')')
end subroutine make_shell2

!*************************************************************************
!* Returns true if centre of mass can be assigned
!*  and returns centre of mass for a mask of atoms in the vector 'centre'
!*************************************************************************
logical function get_centre_by_mass(centre)
real(kind=prec), intent(OUT)   :: centre(3)
type(MASK_TYPE)        :: mask
integer                :: ats, imaskat, iat
real(kind=prec)                :: totmass, mass
get_centre_by_mass = .false.
centre = zero
totmass=zero
call mask_initialize(mask)
ats = maskmanip_make_pretop(mask)
if (.not. allocated(iac))allocate(iac(nat_pro))
call set_iac        !Set info about masses
do iat = 1,nat_pro
  if (mask%mask(iat)) then
    if(heavy(iat)) then !skip hydrogens
      mass=iaclib(iac(iat))%mass
      totmass=totmass+mass
	  centre = centre + mass*xtop(iat*3-2:iat*3)
    end if
  end if
end do
centre = centre/totmass
get_centre_by_mass = .true.
end function get_centre_by_mass
!*************************************************************************

END module PREP
