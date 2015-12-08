! (C) 2014 Uppsala Molekylmekaniska HB, Uppsala, Sweden
! topo.f90
! by John Marelius & Johan Åqvist
! molecular topology data and I/O
!TODO: fix in accordance with best practice

module TOPO
	use SIZES
	use MISC

	implicit none

! constants
        real, private, parameter                ::      MODULE_VERSION = 5.06
        character(*), private, parameter        ::      MODULE_DATE    = '2014-04-21'

        integer, parameter                      ::      nljtyp = 3	!TINY

	integer, parameter			::	max_nbr_range	= 25

! types
	type BOND_TYPE
		integer(AI)				::	i, j
		integer(SHRT)			::	cod
	end type BOND_TYPE

	type ANG_TYPE
		integer(AI)				::	i, j, k
		integer(SHRT)			::	cod
	end type ANG_TYPE

	type BONDLIB_TYPE
		real(kind=prec)					::	fk, bnd0
	end type BONDLIB_TYPE

	type ANGLIB_TYPE
		real(kind=prec)					::	fk, ang0
		real(kind=prec)					::	ureyfk, ureyr0
	end type ANGLIB_TYPE

	type TOR_TYPE
		integer(AI)				::	i, j, k, l
		integer(SHRT)			::	cod
	end type TOR_TYPE

	type TORLIB_TYPE
		real					::	fk, rmult, deltor
		real					::	paths
	end type TORLIB_TYPE

	type IMPLIB_TYPE
		real					::	fk, imp0
	end type IMPLIB_TYPE

	type CGP_TYPE
		integer(AI)				::	iswitch	!switching atom
		integer(AI)				::	first	!index in cgpatom array for 1st member
		integer(AI)				::	last	!index in cgpatom array for last member
	end type CGP_TYPE

	type IAC_TYPE
		real(kind=prec)					::	mass
		real(kind=prec)					::	avdw(nljtyp)
		real(kind=prec)					::	bvdw(nljtyp)
	end type IAC_TYPE

	type LJ2_TYPE
		integer(TINY)			::	i, j
	end type LJ2_TYPE

	type RESIDUE_TYPE
		integer					::	irc, start
		character(4)			::	name
	end type RESIDUE_TYPE

!	Memory management
integer, private			::	alloc_status
!private subroutine
private						::	topo_check_alloc

!topology information
real						::	version
character(len=80)			::	title = ''
character(len=80)			::	forcefield = ''
character(len=8)			::	creation_date
character(len=256)			::	pdb_file = ''
character(len=256)			::	lib_files = ''
character(len=256)			::	prm_file

!atom information
	!nat_pro = total # of atoms in topology, nat_solute = # solute atoms (no water)
	integer											::	nat_pro, nat_solute, max_atom
	real(kind=prec), allocatable				::	xtop(:)			! topology/coordinates
	integer(TINY), allocatable	::	iac(:)			! integer atom codes
	logical, allocatable				::	heavy(:)		! boolean flag, true if atom >= He
	real(kind=prec), allocatable						::	crg(:)			! charges
	real(kind=prec), allocatable						::	chg_solv(:)		! solvent charges
	real(kind=prec), allocatable						::	aLJ_solv(:,:),bLJ_solv(:,:)	! solvent LJ types
        integer :: solv_atom
! number of atoms/solvent molecule
	integer(AI), allocatable		::	cgpatom(:)		! charge groups
	ENUM, bind(c)
        	ENUMERATOR      :: SOLVENT_SPC, SOLVENT_ALLATOM, SOLVENT_GENERAL 
	END ENUM

!	integer, parameter					::	SOLVENT_SPC=0, SOLVENT_3ATOM=1, SOLVENT_GENERAL=2
	integer											::	solvent_type
	integer,allocatable					::	glb_cofactor(:)	! 0 = protein or ligand atom, 1,2,3...= cofactor 1, 2, 3...

	integer, public	:: ligand_offset ! topology index offset for ligand. xtop(1:ligand_offset) = protein

!sphere information
	!!sim. sphere & water sphere centres
	real(kind=prec)						::	xpcent(3), xwcent(3)
	real(kind=prec)						::	rwat !solvation radius
	real(kind=prec)						::	rexcl_o,rexcl_i
	integer						::	nexats, nshellats, nexwat
	logical, allocatable		::	shell(:)
	logical, allocatable		::	excl(:)

!box information
	real(kind=prec)						::	boxlength(3) !length of the boxedges
	real(kind=prec)						::	boxcentre(3) !center coordinates of the box
	real(kind=prec)						::	inv_boxl(3)  !inverse of the boxedges

!flag indication if simulation sphere (.false.) or periodic box (.true.)is used.
	logical						::	use_PBC = .false.

!bond information
	integer												::	nbonds, nbonds_solute, max_bonds
	type(BOND_TYPE), allocatable	::	bnd(:)

	integer						::	nbndcod, max_bondlib
	type(BONDLIB_TYPE), allocatable:: bondlib(:)
	character(len=3), allocatable::	SYBYL_bond_type(:)

!angle information
	integer						::	nangles, nangles_solute, max_angles
	type(ANG_TYPE), allocatable	::	ang(:)

	integer						::	nangcod, max_anglib
	type(ANGLIB_TYPE), allocatable:: anglib(:)

!torsion information
	integer						::	ntors, ntors_solute, ntorcod
	integer						::	max_tors, max_torlib
	integer						::	nimps, nimps_solute, nimpcod
	integer						::	max_imps, max_implib, imp_type
	type(TOR_TYPE),target,allocatable:: tor(:)
	type(TOR_TYPE), target,allocatable:: imp(:)
	type(TORLIB_TYPE), target,allocatable:: torlib(:)
	type(IMPLIB_TYPE), target, allocatable:: implib(:)

!charge group information
	integer						::	ncgp, ncgp_solute, max_cgp
	type(CGP_TYPE), allocatable	::	cgp(:)

!atom type information
	integer												::	natyps, max_atyps
	integer												::	nlj2, max_lj2
	type(IAC_TYPE),allocatable		::	iaclib(:)
	character(len=5),allocatable	::	SYBYL_atom_type(:)
	character(len=8), allocatable ::	tac(:) !text atom codes
	type(LJ2_TYPE), allocatable		::	lj2(:)

!force field options
	integer						::	ivdw_rule !combination rule
	integer, parameter			::	VDW_GEOMETRIC=1, VDW_ARITHMETIC=2
	real(kind=prec)						::	el14_scale !scaling of 1-4 electrostatics
	integer						::	iuse_switch_atom !switching atoms in charge group
	real(kind=prec)						::	coulomb_constant !constant in Coulombs law
	!type of forcefield, one of FF_GROMOS, FF_AMBER, FF_CHARMM
	!this only (right now) affects how improper parameters are assigned
	integer						::	ff_type
	integer, parameter			::	FF_GROMOS=1, FF_AMBER=2, FF_CHARMM=3

!neighbour lists
	integer						::	n14nbrs
	integer						::	nexnbrs, max_exbnrs
	integer						::	n14long, max_14long
	integer						::	nexlong, max_exlong
	logical, allocatable		::	list14(:,:)
	logical, allocatable		::	listex(:,:)
	integer(AI), allocatable	::	list14long(:,:)
	!listexlong may be reallocated to include special exclusions
	integer(AI), pointer		::	listexlong(:,:)

!Residue and molecule bookkeeping information
	INTEGER						::	nres, nres_solute, max_res
	type(RESIDUE_TYPE), allocatable	::	res(:)

	INTEGER						::	nmol, max_mol
	INTEGER, allocatable		::	istart_mol(:)

!Solvent dielectric information
!saved as integer (num * 1000) to make reading easier
        integer                         :: dielectric        

contains

!----------------------------------------------------------------------

subroutine topo_startup
end subroutine topo_startup

!----------------------------------------------------------------------

subroutine topo_set_max(max, max_lib, max_nbr)
!arguments
	integer						::	max, max_lib, max_nbr

	max_atom	= max
	max_bonds	= 2*max
	max_angles	= 3*max
	max_tors	= 12*max
	max_imps	= max

	if(max_lib > MAX_SHRT) then
		write(*,*) 'ERROR: Cannot have larger libraries than ',MAX_SHRT,&
			' with current word length.'
		stop 255
	end if

	max_bondlib = max_lib
	max_anglib	= max_lib
	max_torlib	= max_lib
	max_implib	= max_lib
	max_lj2		= max_lib

	max_14long	= max_nbr
	max_exlong	= max_nbr
	max_cgp		= max

	max_res		= max
	max_mol		= max

end subroutine topo_set_max

!----------------------------------------------------------------------

subroutine topo_allocate_potential(stat_out)
!arguments
	integer, optional, intent(out)::	stat_out

!!$                bnd(max_bonds), &
!!$		bondlib(max_bondlib), &
!!$		SYBYL_bond_type(max_bondlib), &
!!$		ang(max_angles), &
!!$		anglib(max_anglib), &
!!$		tor(max_tors), &
!!$		torlib(max_torlib), &
!!$		imp(max_imps), &
!!$		implib(max_implib), &
!!$		lj2(max_lj2), &
	allocate(iaclib(max_atyps), &
		tac(max_atyps), &
		SYBYL_atom_type(max_atyps), &
		list14long(2, max_14long), &
		listexlong(2, max_exlong), &
		cgp(max_cgp), &
	stat = alloc_status)

	if(alloc_status /= 0) then
		write(*,*) 'ERROR: Out of memory when allocating topology arrays'
		if(.not. present(stat_out)) then
			stop 255
		else
			stat_out = alloc_status
			return
		end if
	end if
end subroutine topo_allocate_potential

!----------------------------------------------------------------------

subroutine topo_allocate_atom(stat_out)
!arguments
	integer, optional, intent(out)::	stat_out

	allocate(glb_cofactor(max_atom), &
						iac(max_atom), &
						crg(max_atom), &
						xtop(3*max_atom), &
						heavy(max_atom), &
						cgpatom(max_atom), &
						list14(max_nbr_range, max_atom), &
						listex(max_nbr_range, max_atom), &
						excl(max_atom), shell(max_atom), &
						stat=alloc_status)

	if(alloc_status /= 0) then
		write(*,*) 'ERROR: Out of memory when allocating topology arrays'
		if(.not. present(stat_out)) then
			stop 255
		else
			stat_out = alloc_status
			return
		end if
	end if

end subroutine topo_allocate_atom

!----------------------------------------------------------------------

subroutine topo_deallocate (keep_ff)
!arguments
	logical, optional, intent(in)::	keep_ff
!locals
	logical						::	keep

	keep = .false.

	if(present(keep_ff)) then
		if(keep_ff) keep = .true.
	end if

	if(.not. keep .and. allocated(iaclib)) deallocate(iaclib)
	if(.not. keep .and. allocated(SYBYL_atom_type)) deallocate(SYBYL_atom_type)
	if(.not. keep .and. allocated(lj2)) deallocate(lj2)
	if(.not. keep .and. allocated(tac)) deallocate(tac)
	if(allocated(glb_cofactor)) deallocate(glb_cofactor)
	if(allocated(bnd)) deallocate(bnd)
	if(allocated(bondlib)) deallocate(bondlib)
	if(allocated(SYBYL_bond_type)) deallocate(SYBYL_bond_type)
	if(allocated(ang)) deallocate(ang)
	if(allocated(anglib)) deallocate(anglib)
	if(allocated(tor)) deallocate(tor)
	if(allocated(torlib)) deallocate(torlib)
	if(allocated(imp)) deallocate(imp)
	if(allocated(implib)) deallocate(implib)
	if(allocated(iac)) deallocate(iac)
	if(allocated(crg)) deallocate(crg)
	if(allocated(xtop)) deallocate(xtop)
	if(allocated(heavy)) deallocate(heavy)
	if(allocated(cgpatom)) deallocate(cgpatom)
	if(allocated(list14)) deallocate(list14)
	if(allocated(listex)) deallocate(listex)
	if(allocated(res)) deallocate(res)
	if(allocated(istart_mol)) deallocate(istart_mol)
	if(allocated(cgp)) deallocate(cgp)
	if(allocated(list14long)) deallocate(list14long)
	if(allocated(excl)) deallocate(excl)
	if(allocated(shell)) deallocate(shell)
	if(allocated(lj2)) deallocate(lj2)
	deallocate(listexlong, stat=alloc_status)
end subroutine topo_deallocate

!----------------------------------------------------------------------

subroutine topo_reallocate_xtop(atoms)
!arguments
	integer						::	atoms

	real(kind=prec), allocatable		::	r8temp(:)
	integer						::	nat3old_array(1), nat3old

	nat3old_array = ubound(xtop)
	nat3old = nat3old_array(1)
	allocate(r8temp(nat3old), stat=alloc_status)
	call topo_check_alloc('reallocating topology atom array')
	r8temp(1:nat3old) = xtop(1:nat3old)
	deallocate(xtop)
	allocate(xtop(atoms*3),  stat=alloc_status)
	call topo_check_alloc('reallocating topology atom array')
	xtop(1:nat3old) = r8temp(1:nat3old)
	deallocate(r8temp)

end subroutine topo_reallocate_xtop

!----------------------------------------------------------------------

subroutine topo_reallocate(oldatoms, atoms, waters)
!arguments
	integer						::	oldatoms, atoms, waters
!locals
	integer						::	oldbonds, bonds
	integer						::	oldangles, angles
	real(kind=prec), allocatable			::	r4temp(:)
	integer(1), allocatable		::	i1temp(:)
	type(BOND_TYPE), allocatable::	bndtemp(:)
	type(ANG_TYPE), allocatable	::	angtemp(:)

	oldbonds = nbonds
	bonds = oldbonds + 2*waters
	oldangles = nangles
	angles = oldangles + 2*waters

	!realloc bonds
	allocate(bndtemp(oldbonds),  stat=alloc_status)
	call topo_check_alloc('reallocating atom arrays')
	bndtemp(1:oldbonds) = bnd(1:oldbonds)
	deallocate(bnd)
	allocate(bnd(bonds),  stat=alloc_status)
	call topo_check_alloc('reallocating atom arrays')
	bnd(1:oldbonds) = bndtemp(1:oldbonds)
	deallocate(bndtemp)

	!realloc angles
	allocate(angtemp(oldangles),  stat=alloc_status)
	call topo_check_alloc('reallocating atom arrays')
	angtemp(1:oldangles) = ang(1:oldangles)
	deallocate(ang)
	allocate(ang(angles),  stat=alloc_status)
	call topo_check_alloc('reallocating atom arrays')
	ang(1:oldangles) = angtemp(1:oldangles)
	deallocate(angtemp)

	allocate(r4temp(oldatoms), stat=alloc_status)
	call topo_check_alloc('reallocating atom arrays')

	!realloc. crg
	r4temp(1:oldatoms) = crg(1:oldatoms)
	deallocate(crg)
	allocate(crg(atoms), stat=alloc_status)
	call topo_check_alloc('reallocating atom arrays')
	crg(1:oldatoms) = r4temp(1:oldatoms)

	!done with real's
	deallocate(r4temp)
	allocate(i1temp(oldatoms), stat=alloc_status)
	call topo_check_alloc('reallocating atom arrays')

	!realloc. iac
	i1temp(1:oldatoms) = iac(1:oldatoms)
	deallocate(iac)
	allocate(iac(atoms), stat=alloc_status)
	call topo_check_alloc('reallocating atom arrays')
	iac(1:oldatoms) = i1temp(1:oldatoms)

	!done with i1's
	deallocate(i1temp)

end subroutine topo_reallocate

!----------------------------------------------------------------------

subroutine topo_check_alloc(message)
	implicit none

	character(*) message

	if(alloc_status .ne. 0) then
		write(*,*) &
			'>>> Out of memory trying to allocate ', message
		stop 255
	end if
end subroutine topo_check_alloc

!----------------------------------------------------------------------

logical function topo_load(filename, require_version)
!arguments
	character(*)			::	filename
	real, intent(in)		::	require_version

!locals
	integer					::	u, readflag

	topo_load = .false.

	! try opening the file
	u = topo_open(filename)
	if(u <= 0) then
		write(*,'(a,a)') 'ERROR: Failed to open topology file ', &
			trim(filename)
		return
	end if

	! topo_read
	topo_load = topo_read(u, require_version)

! close file
	close(u)
end function topo_load

!----------------------------------------------------------------------

integer function topo_open(filename)
!arguments
	character(*)			::	filename

!locals
	integer					::	u, stat_out

	u = freefile()

	open(unit = u, file=filename, status='old',form='formatted', &
		   action='read', iostat = stat_out)
	if(stat_out /= 0) then
		write(*,10) trim(filename)
		topo_open = 0
		return
	end if

	topo_open = u

10	format('ERROR: Failed to open topology file ', a)

end function topo_open

!----------------------------------------------------------------------
logical function topo_read(u, require_version, extrabonds)
	! arguments
	integer						::	u
	real, intent(in)			::	require_version
	integer, optional,intent(in)::	extrabonds

	! local variables
	integer                     :: rd
	integer						:: nat3, nwat
	integer						:: i,j,k,n, si
	integer						:: nhyds
	integer						:: paths
	character(len=256)			:: line, restofline
	character(len=10)			:: key, boundary_type !PWadded
	integer						:: filestat
	integer(1), allocatable		:: temp_list(:,:)
	integer						:: extra
	real(kind=prec)						:: deprecated

	topo_read = .false.
	! get rid of old topology
	call topo_deallocate

	!	This is the default / cutoffs based on switching atoms
	iuse_switch_atom = 1

	!	read topology file
	call centered_heading('Reading topology file', '-')

	!handle extra bonds etc
	if(present(extrabonds)) then
		extra = extrabonds
	else
		extra = 0
	end if

  !	read topology file

  ! --> 1. header or title row
  ! ==========================

	!clear header information
	pdb_file = ''
	prm_file = ''
	lib_files = ''
	forcefield = ''
	ff_type = FF_GROMOS !this is the default

	read (u, '(a80)', err=1000) title
	if(title == 'Q topology file') then !it's a file with header
		do
			read(u,'(a256)') line
			read(line,*) key
			!workaround to allow / in strings
			restofline = adjustl(line(len_trim(key)+1:len(line)))
			call upcase(key)
			select case (key)
			case ('TITLE')
				title = restofline
			case ('DATE')
				creation_date = restofline
			case ('VERSION')
				read(restofline, *) version
			case ('PDB_FILE')
				pdb_file = restofline
			case ('LIB_FILES')
				lib_files = restofline
			case ('FORCEFIELD')
				forcefield = restofline
			case ('FF_TYPE')
				read(restofline, *) ff_type
			case ('PRM_FILE')
				prm_file = restofline
			case ('END')
				exit
			case default
				write(*,'(a,a)') '>>>WARNING: Unrecognised key in header: ', key
			end select
		end do
	else
		version = 2
	end if
	write (*,'(a80,/)') title

	if(version < require_version) goto 1100

  ! --> 2. nat_pro (no. of atoms in topology)
  ! =========================================
	read (u, '(a)', err=1000) line
	read (line,*,err=1000) nat_pro !first read what must be there

	read (line,*,iostat=filestat) nat_pro, nat_solute, solv_atom, dielectric !then try to read optional things too
	if (filestat .ne. 0) then
! maybe just does not have dielectric and solv atom, try without first
          read (line,*,iostat=filestat) nat_pro, nat_solute, solv_atom
          if (filestat .ne. 0) then
! ok, maybe even solv_atom is not there??
            read (line,*,iostat=filestat) nat_pro, nat_solute
	    if (filestat .ne. 0) then
	      nat_solute = nat_pro !default is nat_solute = nat_pro
            end if
            solv_atom  = 3        !default for solv atom to keep rest from crashing
            dielectric = int(80.1*1000) !default for dielectric (of water) to keep rest working
          end if
          dielectric = int(80.1*1000) !default for dielectric (of water) to keep rest working
        end if

	write (*,20) nat_solute, nat_pro-nat_solute
20	format ('No. of solute atoms     = ',i10,/,&
			'No. of solvent atoms    = ',i10)
	max_atom = nat_pro
	nwat = (max_atom - nat_solute) / solv_atom
        write (*,'(a,i10)') 'Solvent dielectric*1000 =',dielectric
	call topo_allocate_atom

	if ( nat_pro .eq. 0 ) write (*,'(a)')  &
	   '>>> WARNING: Zero topology atoms.'


  ! --> 3. topology coordinates    ---->     xtop
  ! ===========================    ==============
   if(nat_pro > 0) read (unit=u, fmt=*, err=1000) (xtop(rd),rd=1,3*nat_pro)
   write (*,40) 3*nat_pro
   40 format ('No. of coordinates      = ',i10)


  ! --> 4. integer atom codes      ---->     iac
  ! =========================
  read (unit=u, fmt=*, err=1000)                      ! skip
  if(nat_pro > 0) read (unit=u, fmt=*, err=1000) (iac(rd),rd=1,nat_pro)
  write (*,50) nat_pro
50 format ('No. of atom type codes  = ',i10)


  ! --> 5. bond list and params    ---->     bnd, bondlib
  ! ===========================
	read(u,'(a)',err=1000) line
	read(line,*,err=1000) nbonds !first read what must be there
	read (line, fmt=*, iostat=filestat) nbonds, nbonds_solute
	if(filestat /= 0) nbonds_solute = nbonds
	write (*,60) nbonds_solute
	write (*,61) nbonds-nbonds_solute

60 format ('No. of solute bonds     = ',i10)
61 format ('No. of solvent bonds    = ',i10)
	max_bonds = nbonds + 2*nwat + extra
  allocate(bnd(max_bonds), stat=alloc_status)
  call topo_check_alloc('bond list')

  if(nbonds > 0) read (unit=u, fmt=*, err=1000) (bnd(rd),rd=1,nbonds)

	read(unit=u, fmt=*, err=1000) nbndcod
  !allocate one extra for water bond
  max_bondlib = nbndcod + 1 + extra
  allocate(bondlib(max_bondlib), SYBYL_bond_type(max_bondlib), &
	stat=alloc_status)
  call topo_check_alloc('bond library')

  do i=1,nbndcod
    ! initialize sybyl bond type to nothing
    SYBYL_bond_type(i) = '   '

	read (unit=u, fmt='(a)', err=1000) line
	read (line, fmt=*, end=700, err=700) j, bondlib(i)%fk, bondlib(i)%bnd0, SYBYL_bond_type(i)

	!restofline
!	call sscan(restofline, SYBYL_bond_type(i)) ! get WS-separated string
700	end do


	! --> 6. angle list and params    ---->     ang, anglib
	! ============================
	read(u,'(a)',err=1000) line
	read(line,*,err=1000) nangles
	read(line, *, iostat=filestat) nangles, nangles_solute
	if(filestat /= 0) nangles_solute = nangles
	write (*,80) nangles_solute
	write (*,81) nangles-nangles_solute
80  format ('No. of solute angles    = ',i10)
81  format ('No. of solvent angles   = ',i10)
	max_angles = nangles + nwat + 2*extra
	allocate(ang(max_angles), stat=alloc_status)
	call topo_check_alloc('angle list')

	if(nangles > 0) read (unit=u, fmt=*, err=1000) ( ang(si), si = 1,nangles )

	read(unit=u, fmt=*, err=1000) nangcod

	!allocate one extra for water angle
	max_anglib = nangcod + 1 + extra
	allocate(anglib(max_anglib), stat=alloc_status)
	call topo_check_alloc('angle library')
	!set optional parameters to zero
	anglib(:)%ureyfk=0.0_prec
	anglib(:)%ureyr0=0.0_prec

	do i=1,nangcod
		!handle optional Urey-Bradley params
		read(unit=u, fmt='(a)', err=1000) line
		!try to read all things
		read(unit=line, fmt=*, iostat=filestat) j, anglib(i)
		!accept end-of-data but not read errors
		if(filestat > 0) goto 1000
	end do


	! --> 7. torsion list and params    ---->     tor, torlib
	! ==============================
	read(u,'(a)',err=1000) line
	read(line,*,err=1000) ntors
	read(line, *, iostat=filestat) ntors, ntors_solute
	if(filestat /= 0) ntors_solute = ntors
	write (*,100) ntors_solute
	write (*,101) ntors-ntors_solute
100 format ('No. of solute torsions  = ',i10)
101 format ('No. of solvent torsions = ',i10)
	max_tors = ntors + 4*extra
	allocate(tor(max_tors), stat=alloc_status)
	call topo_check_alloc('torsion list')
	if(ntors > 0) read (unit=u, fmt=*, err=1000) (tor(si), si = 1,ntors)

	read (unit=u, fmt=*, err=1000) ntorcod
	max_torlib = ntorcod + extra
	allocate(torlib(max_torlib), stat=alloc_status)
	call topo_check_alloc('torsion library')

	do i=1,ntorcod
		read (unit=u, fmt=*, err=1000) j,torlib(i)%fk,torlib(i)%rmult, &
			torlib(i)%deltor,torlib(i)%paths
	end do

	! --> 8. improper torsion list and params    ---->     imp, implib
	! =======================================
	read(u,'(a)',err=1000) line
	read(line,*,err=1000) nimps
	read(line, *, iostat=filestat) nimps, nimps_solute
	if(filestat /= 0) nimps_solute = nimps
	write (*,120) nimps_solute
	write (*,121) nimps-nimps_solute
120 format ('No. of solute impropers = ',i10)
121 format ('No. of solvent impropers= ',i10)
	max_imps = nimps + 4*extra
	allocate(imp(max_imps), stat=alloc_status)
	call topo_check_alloc('impoper torsion list')
	if(nimps > 0) read (unit=u, fmt=*, err=1000) imp(1:nimps)

	read (unit=u, fmt=*, iostat=filestat) nimpcod, imp_type
	if(filestat /= 0) then !water atom type data present
		imp_type = 1 !default to harmonic
	end if
	max_implib = nimpcod + extra
	allocate(implib(max_implib), stat=alloc_status)
	call topo_check_alloc('improper torsion library')

	do i=1,nimpcod
		read (unit=u, fmt=*, err=1000) j,implib(i)
	end do


	! --> 9. charges charges charges		    ---->     crg
	! ==============================
	read (unit=u, fmt='(a)', err=1000) line
	read (unit=line, fmt=*, iostat=filestat) nat_pro

	if(nat_pro > 0) read (unit=u, fmt=*, err=1000) crg(1:nat_pro)
	write (*,130) nat_pro
130 format ('No. of atomic charges   = ',i10)


	! --> 10. charge groups					    ---->     cgp
	! =====================
	read(u,'(a)',err=1000) line
	read(line,*,err=1000) ncgp !first read what must be there
	read (line, fmt=*, iostat=filestat) i, ncgp_solute
	if(filestat /= 0) ncgp_solute = ncgp
	read (line, fmt=*, iostat=filestat) i, i, iuse_switch_atom
	if(filestat /= 0) iuse_switch_atom = 1 !default to using sw. atoms
	write (*,140) ncgp_solute
	write (*,141) ncgp-ncgp_solute
140 format ('No. of solute chargegrps= ',i10)
141 format ('No. of solvent chargegps= ',i10)
	max_cgp = ncgp
	allocate(cgp(max_cgp), stat=alloc_status)
	call topo_check_alloc('charge group list')

	k=1
	do i=1,ncgp
		read (unit=u, fmt=*, err=1000) n, cgp(i)%iswitch
		cgp(i)%first = k
		cgp(i)%last = k + n - 1
		read(unit=u, fmt=*, err=1000) (cgpatom(j), j=cgp(i)%first,cgp(i)%last)
		k = k + n
	end do

	if(iuse_switch_atom == 0 ) then
		write (*,142) 'Any atoms'
	else
		write (*,142) 'Switching atoms'
	end if
142	format ('Charge group cut-off    : ',a)



	! --> 11. Masses and van der Waals params
	! =======================================

	!read water atom types to temp. variables to avoid clearing default values
	read (unit=u, fmt=*, err=1000) natyps
	max_atyps = natyps
	allocate(iaclib(max_atyps), SYBYL_atom_type(max_atyps), tac(max_atyps), &
		stat=alloc_status)
	call topo_check_alloc('atom type library')

	read (unit=u, fmt=*, err=1000) ivdw_rule
	write (*,162) ivdw_rule
162 format ('vdW  rule [ 1=G / 2=A ] = ',i10)

	read (u, '(a)', err=1000) line
	read (line,*, err=1000) el14_scale !first read what must be there
	!then try to read optional things
	read (line,*, iostat=filestat) el14_scale, coulomb_constant
	!set coulomb_constant to default if not present
	if(filestat /= 0. .or. coulomb_constant <= 0) coulomb_constant = 332.
	write (*,164) el14_scale
	write (*,165) coulomb_constant
164 format ('El-static 1-4 damping   = ',f10.3)
165 format ('Coulomb constant        = ' ,f10.4)

	read (unit=u, fmt=*, err=1000)                      ! skip
	read (unit=u, fmt=*, err=1000) (iaclib(i)%mass, i=1, natyps)

	do j=1,nljtyp
		read (unit=u, fmt=*, err=1000)                   ! skip
		read (unit=u, fmt=*, err=1000) (iaclib(i)%avdw(j),i=1,natyps)
		read (unit=u, fmt=*, err=1000)                   ! skip
		read (unit=u, fmt=*, err=1000) (iaclib(i)%bvdw(j),i=1,natyps)
	end do

	write (*,160) natyps
160 format ('No. of atom types       = ',i10)

	read (unit=u, fmt=*, err=1000) nlj2
	max_lj2 = nlj2
	allocate(lj2(nlj2), stat=alloc_status)
	call topo_check_alloc('vdw type 2 pair list')

	do i=1,nlj2
		read (unit=u, fmt=*, err=1000) lj2(i)%i,lj2(i)%j
	end do

	write (*,161) nlj2
161 format ('No. of LJ type 2 pairs  = ',i10)

	!
	!	flag heavy atoms for various reasons
	!
	nhyds = 0
	do i=1,nat_pro
		if(iaclib(iac(i))%mass < 4.0_prec) then     ! lighter than He = H
			heavy(i) = .false.
			nhyds = nhyds+1
		else
			heavy(i) = .true.
		end if
	end do
	write (*,170) nat_pro-nhyds
170 format ('No. of heavy atoms      = ',i10)



	! --> 12. 1-4 neighbour and exclusion lists
	! =========================================
	read (unit=u, fmt=*, err=1000) n14nbrs
	write (*,180) n14nbrs
180 format ('No. of 1-4 neighbours   = ',i10)

	allocate(temp_list(max_nbr_range,nat_solute), stat=alloc_status)
	if(nat_solute > 0) then
	! work-around to preserve file format
	! read (unit=u,fmt='(80i1)', err=1000) ((list14(i,j),i=1,max_nbr_range),j=1,nat_solute)
	  call topo_check_alloc('temporary neighbor list')
	  read (unit=u,fmt='(80i1)', err=1000) ((temp_list(i,j),i=1,max_nbr_range),j=1,nat_solute)
	  do j=1,nat_solute
		do i=1,max_nbr_range
		  list14(i,j) = (temp_list(i,j) .eq. 1)
		end do
	  end do
	end if

	read (unit=u, fmt=*, err=1000) n14long
	write (*,200) n14long
200 format ('No. long-range 1-4 nbrs = ',i10)
	max_14long = n14long
	allocate(list14long(2, max_14long), stat=alloc_status)
	call topo_check_alloc('long-range neighbour list')
	do i = 1, n14long
		read (unit=u, fmt=*, err=1000) (list14long(j,i),j=1,2)
	end do

	read (unit=u, fmt=*, err=1000) nexnbrs
	write (*,220) nexnbrs
220 format ('No. of nbor exclusions  = ',i10)

	if(nat_solute > 0) then
	  read (unit=u,fmt='(80i1)', err=1000) ((temp_list(i,j),i=1,max_nbr_range),j=1,nat_solute)
	  do j=1,nat_solute
		do i=1,max_nbr_range
		  listex(i,j) = (temp_list(i,j) .eq. 1)
		end do
	  end do
	end if

	deallocate(temp_list)

	read (unit=u, fmt=*, err=1000) nexlong
	write (*,240) nexlong
240 format ('No. of long-range excls = ',i10)
	max_exlong = nexlong
	allocate(listexlong(2, max_exlong), stat=alloc_status)
	call topo_check_alloc('long-range exclusion list')

	do i = 1, nexlong
		read (unit=u, fmt=*, err=1000) (listexlong(j,i),j=1,2)
	end do

! --- RESIDUE/MOLECULE BOOKKEEPING
	read(u,'(a)', err=1000) line
	!Fix for inconsistent write/read of empty lists on SGI
	if(line == '') read(u,'(a)', err=1000) line
	READ(line, fmt=*, err=1000) nres
	!attempt to read #solute residues
	read(line, *, iostat=filestat) nres, nres_solute
	if (filestat .ne. 0) then
	  nres_solute = nres ! set # nolute residues to total #
	end if

	allocate(res(nres), stat=alloc_status)
	call topo_check_alloc('residue list')

	if(nres > 0) READ(unit=u, fmt=*, err=1000) res(1:nres)%start
	READ(unit=u, fmt=*, err=1000)
	if(nres > 0) then
	    do i = 0, int((nres+15)/16)-1
		READ(unit=u, fmt='(16(a4,1x))', err=1000) (res(i*16+j)%name, j = 1, min(16, nres-i*16))
	    end do
	endif
	WRITE(*, '(a,i10)') 'No. of residues         = ', nres
	if(nres_solute < nres) then
		WRITE(*, '(a,i10)') 'No of solute residues   = ', nres_solute
	end if

	READ(unit=u, fmt=*, err=1000) nmol
	allocate(istart_mol(nmol+1), stat=alloc_status)
	call topo_check_alloc('molecule list')
	if(nmol > 0) READ(unit=u, fmt=*, err=1000)(istart_mol(i), i = 1, nmol)
	WRITE(*, '(a,i10)') 'No. of molecules        = ', nmol

	read(unit=u, fmt=*, iostat=filestat)
	if(natyps > 0) then
	    do i = 0, int((natyps+7)/8)-1
		read(unit=u, fmt='(8(a8,1x))', err=900, end=900) (tac(i*8+j), j = 1, min(8, natyps-i*8))
	    end do
	endif
	if(filestat /= 0) then
		if(version > 3.5) then
			goto 1000 !its an error
		elseif(require_version > 3.5) then
			goto 1100
		else
			goto 900
		end if
	end if
	WRITE(*, '(a,i10)') 'Atom type names         = ', natyps

	read(unit=u, fmt=*, iostat=filestat)
	if(natyps > 0) then
	    do i = 0, int((natyps+12)/13)-1
		read(unit=u, fmt='(13(a5,1x))', iostat=filestat) &
		    (SYBYL_atom_type(i*13+j), j = 1, min(13, natyps-i*13))
	    end do
	endif
	if(filestat /= 0) then
		if(version > 3.5) then
			goto 1000 !its an error
		elseif(require_version > 3.5) then
			goto 1100
		else
			goto 900
		end if
	end if

	WRITE(*, '(a,i10)') 'SYBYL atom types        = ', natyps

	if(version < 4) then
		if(require_version >= 4) then
			goto 1100
		else
			goto 900
		end if
	end if

	read(unit=u, fmt=*, err=1000) solvent_type

	!******PWadded input reading 2001-10-10
	read (u, '(a)', err=1000) line
	read(line, *, err=1000) boundary_type
	if(boundary_type == 'PBC') then
		use_PBC = .true.
		write(*,'(a)')'Boundary: periodic box'
		read(unit=u, fmt=*, err=1000) boxlength(:)
		write(*,'(a,3f10.3)') 'Box size                = ', boxlength(:)
		read(unit=u, fmt=*, err=1000) boxcentre(:)
		write(*,'(a,3f10.3)') 'Box centre              = ', boxcentre(:)
                excl(:) = .false.
                shell(:) = .false.

	else
		use_PBC = .false.
		write(*,'(a)')'Boundary: sphere'

    !For backward compatibility the shell radius from the topology can be used.
    if (version < 5.01) then
      print*,'Using shell radius from topology < v5.01! '
		  read(line, *, err=1000) rexcl_o, rexcl_i, rwat
		  write(*, '(a,f10.3)') 'Exclusion radius        = ', rexcl_o
      write(*, '(a,f10.3)') 'Restrained shell radius = ', rexcl_i
		  write(*, '(a,f10.3)') 'Eff. solvent radius     = ', rwat
    else
		  read(line, *, err=1000) rexcl_o, rwat
		  write(*, '(a,f10.3)') 'Exclusion radius        = ', rexcl_o
		  write(*, '(a,f10.3)') 'Eff. solvent radius     = ', rwat
    end if
		read(unit=u, fmt=*, err=1000) xpcent(:)
		write(*, '(a,3f10.3)') 'Solute centre           = ', xpcent(:)
		read(unit=u, fmt=*, err=1000) xwcent(:)
		write(*, '(a,3f10.3)') 'Solvent centre          = ', xwcent(:)

		excl(:) = .false.
		shell(:) = .false.
		read(unit=u, fmt=*, err=1000, end=1000) nexats, nexwat
		read(unit=u, fmt='(80l1)') excl(1:nat_pro)
		WRITE(*, '(a,i10)') 'No. of excluded atoms   = ', nexats

    endif
900	WRITE(*, '(a,/)') 'Molecular topology read successfully.'

	topo_read = .true.
	return

1000 write(*,1001)
1001 format('>>>>> ERROR: Could not read topology file.')
	return

1100 write(*,1101) version, require_version
1101 format('>>>>> ERROR: Incompatible topology version ',f5.2, &
				' found. Version >= ',f5.2,' required.')
	return
end function topo_read

!----------------------------------------------------------------------
subroutine topo_save(name)
	!arguments
	character(*)				:: name

	!locals
	integer						:: i, j, u, ig, si
	integer(1), allocatable		:: temp_list(:,:)
    real                        :: crgtot

10	format(a, t29,': ')
20	format(i6)
30	format(f6.2)
	u = freefile() !get an unused unit number
	open(unit=u, file=name, status='unknown', form='formatted', &
		action='write', err = 100)
	call centered_heading('Writing topology file', '-')
!---  HEADER
2	format(a,t12,a)
4	format(a,t12,i1)

	write(u, '(a)') 'Q topology file'
	if(title > '') write(u, 2) 'TITLE', trim(title)
	write(u, 2) 'DATE', trim(creation_date)
	write(u, '(a,t12,f5.2)') 'VERSION', MODULE_VERSION
	if(pdb_file > '') write(u, 2) 'PDB_FILE', trim(pdb_file)
	if(lib_files > '') write(u, 2) 'LIB_FILES', trim(lib_files)
	if(forcefield > '') write(u,2) 'FORCEFIELD', trim(forcefield)
	write(u, 4) 'FF_TYPE', ff_type
	write(u, 2) 'PRM_FILE', trim(prm_file)
	write(u, 2) 'END', 'of header'

	write(*, 10, advance='no') 'solute atoms'
	write(*, 20) nat_solute
	write(*, 10, advance='no') 'solvent atoms'
	write(*, 20) nat_pro-nat_solute
        write(*, 10, advance='no') 'solvent mol atoms'
        write(*, 20) solv_atom
        write(*, 10, advance='no') 'solvent dielectric as num*1000'
        write(*, 20) dielectric

! --- NAT_PRO / COORDINATES
	write(*, 10, advance='no') 'co-ordinates'
	write(u, '(4i8,a)') nat_pro, nat_solute, solv_atom, dielectric, &
		' = Total no. of atoms, no. of solute atoms, atoms per solvent molecule, solvent dielectric as eps*1000. Coordinates: (2*3 per line)'
	if(nat_pro > 0) write(u, '(2(3(f9.3,1x),1x))') ( xtop(si), si = 1,3*nat_pro )
	write(*, 20) 3*nat_pro

! --- INTEGER ATOM CODES
	write(*, 10, advance='no') 'integer atom codes'
	write(u, '(i9,a)') nat_pro, ' = No. of integer atom codes. iac''s: '
	if(nat_pro > 0)	write(u, '(16(i4,1x))') iac(1:nat_pro)
	write(*, 20) nat_pro

! --- BONDS
	write(*, 10, advance='no') 'solute bonds'
	write(u, '(2i9,a)') nbonds, nbonds_solute, &
		 ' = No. of bonds, no. of solute bonds. i - j - icode: (5 per line)'
	if(nbonds>0) write(u, '(5(i7,1x,i7,1x,i3,1x))') ( bnd(si), si = 1,nbonds )
	write(*, 20) nbonds_solute
	write(*, 10, advance='no') 'solvent bonds'
	write(*, 20) nbonds-nbonds_solute

	write(*, 10, advance='no') 'bond parameters'
	write(u, '(i9,a)') nbndcod, ' = No. of bond codes. Parameters: '
	DO i = 1, nbndcod
		write(u, '(i7,f9.3,1x,f10.4,1x,a2)') i, bondlib(i), SYBYL_bond_type(i)
	enddo
	write(*, 20) nbndcod
! --- ANGLES
	write(*, 10, advance='no') 'solute angles'
	write(u, '(2i9,a)') nangles, nangles_solute, &
		' = No. of angles, no. of solute angles. i - j - k - icode: (3 per line)'
	if(nangles>0) write(u, '(3(i7,1x,i7,1x,i7,1x,i3,1x))') ( ang(si), si = 1,nangles )
	write(*, 20) nangles_solute
	write(*, 10, advance='no') 'solvent angles'
	write(*, 20) nangles-nangles_solute

	write(*, 10, advance='no') 'angle parameters'
	write(u, '(i8,a)') nangcod, ' = No. of angle codes. Parameters:'
	DO i = 1, nangcod
		write(u, '(i7,f9.3,1x,f9.3,1x,f9.3,1x,f11.5)') i, anglib(i)
	enddo
	write(*, 20) nangcod

! --- TORSIONS
	write(*, 10, advance='no') 'torsions'
	write(u, '(2i9,a)') ntors, ntors_solute, &
		' = No. of torsions, solute torsions. i - j - k - l - icode: (2 per line)'
	if(ntors>0) write(u, '(2(4(i7,1x),i3,1x))') (tor(si), si = 1,ntors)



	write(*, 20) ntors

	write(*, 10, advance='no') 'torsion parameters'
	write(u, '(i9,a)') ntorcod, ' = No. of torsion codes. Parameters:'
	DO i = 1, ntorcod
		!%paths is a real but we write it as an integer for compatibility
		write(u, '(i7,3(f9.3,1x),i5 )') i, &
			torlib(i)%fk , torlib(i)%rmult, torlib(i)%deltor, int(torlib(i)%paths)
	enddo
	write(*, 20) ntorcod

! --- IMPROPERS
	write(*, 10, advance='no') 'impropers'
	write(u, '(2i9,a)') nimps, nimps_solute, &
		' = No. of impropers, solute impr. i - j - k - l - icode: (2 per line)'
	if(nimps>0) write(u, '(2(4(i7,1x),i3,1x))') imp(1:nimps)
	write(*, 20) nimps

	write(*, 10, advance='no') 'improper parameters'
	write(u, '(2i9,a)') nimpcod, imp_type, &
		' = No. of improper codes, type (1=harmonic,2=periodic). Parameters:'
	DO i = 1, nimpcod
		write(u, '(i7,f10.3,1x,f10.3)') i, implib(i)
	enddo
	write(*, 20) nimpcod

! --- CHARGES
	write(*, 10, advance='no') 'charges'
	write(u, '(i9,a)') nat_pro, ' = No. of atomic charges'
	if(nat_pro >0) write(u, '(10(f7.4,1x))') crg(1:nat_pro)
	write(*, 20) nat_pro

    crgtot = 0.0
    if (use_PBC) then
      do i = 1, nat_solute
        crgtot = crgtot + crg(i)
      end do
	write(*, 10, advance='no') 'total charge of system'
    else
      do i = 1, nat_solute
        if (.not. excl(i)) crgtot = crgtot + crg(i)
      end do
	  write(*, 10, advance='no') 'total charge of not excluded'
	end if
    write (*, 30) crgtot

! --- CHARGE GROUPS
	write(*, 10, advance='no') 'charge groups'
	write(u, '(3i9,a)') ncgp, ncgp_solute, iuse_switch_atom, &
		' = No. of charge groups, no of solvent cgps, switch atoms flag. nat_cgp, iswitch / atom list: '
	DO ig = 1, ncgp
		write(u, '(i7,1x,i7)') cgp(ig)%last - cgp(ig)%first + 1, cgp(ig)%iswitch
		write(u, '(13(i7,1x))') cgpatom(cgp(ig)%first : cgp(ig)%last)
	enddo
	write(*, 20) ncgp

! --- ATOM TYPES
	write(*, 10, advance='no') 'atom type parameters'
	write(u, '(i9,a)') max_atyps, ' = No. of atom types'
	write(u, '(i9,a)') ivdw_rule, &
		' = vdW combination rule (1 = Geom. / 2 = Arit.)'
	write(u, '(f8.5,f9.4,a)') el14_scale, coulomb_constant, &
		' = Electrostatic 1-4 scaling factor and  Coulomb constant'

	write(u, '(a)') 'Masses: '
	write(u, '(10(f7.3,1x))') (iaclib(i)%mass, i=1, max_atyps)
	IF(ivdw_rule==VDW_GEOMETRIC) then
		write(u, '(a)') 'sqrt (Aii) normal:'
		write(u, '(8(f8.2,1x))') (iaclib(i)%avdw(1), i=1, max_atyps)
		write(u, '(a)') 'sqrt (Bii) normal:'
		write(u, '(8(f8.2,1x))') (iaclib(i)%bvdw(1), i=1, max_atyps)
		write(u, '(a)') 'sqrt (Aii) polar:'
		write(u, '(8(f8.2,1x))') (iaclib(i)%avdw(2), i=1, max_atyps)
		write(u, '(a)') 'sqrt (Bii) polar:'
		write(u, '(8(f8.2,1x))') (iaclib(i)%bvdw(2), i=1, max_atyps)
		write(u, '(a)') 'sqrt (Aii) 1-4:'
		write(u, '(8(f8.2,1x))') (iaclib(i)%avdw(3), i=1, max_atyps)
		write(u, '(a)') 'sqrt (Bii) 1-4:'
		write(u, '(8(f8.2,1x))') (iaclib(i)%bvdw(3), i=1, max_atyps)
	elseif(ivdw_rule==VDW_ARITHMETIC) then
		write(u, '(a)') 'R* normal:'
		write(u, '(8(f7.4,1x))') (iaclib(i)%avdw(1), i=1, max_atyps)
		write(u, '(a)') 'epsilon normal:'
		write(u, '(8(f10.7,1x))') (iaclib(i)%bvdw(1), i=1, max_atyps)
		write(u, '(a)') 'R* polar:'
		write(u, '(8(f8.4,1x))') (iaclib(i)%avdw(2), i=1, max_atyps)
		write(u, '(a)') 'epsilon polar:'
		write(u, '(8(f10.7,1x))') (iaclib(i)%bvdw(2), i=1, max_atyps)
		write(u, '(a)') 'R* 1-4:'
		write(u, '(8(f8.4,1x))') (iaclib(i)%avdw(3), i=1, max_atyps)
		write(u, '(a)') 'epsilon 1-4:'
		write(u, '(8(f10.7,1x))') (iaclib(i)%bvdw(3), i=1, max_atyps)
	ENDIF
	write(*, 20) max_atyps

	write(*, 10, advance='no') 'polar LJ pairs'
	write(u, '(i9,a)') nlj2, ' = No. of type-2 vdW interactions. pairs: '
	DO i = 1, nlj2
		write(u, '(i4,1x,i4)') lj2(i)
	enddo
	write(*, 20) nlj2

! --- 1-4 NEIGHBOUR LIST

	write(*, 10, advance='no') 'neighbour list'
	write(u, '(i9,a)') n14nbrs, ' = No. of 1-4 neighbours. nborlist (range=max_nbr_range): '

	! work-around to preserve file format
!	write(u, '(80i1)') list14(1:max_nbr_range, 1:nat_solute)
	allocate(temp_list(max_nbr_range,nat_solute), stat=alloc_status)
	call topo_check_alloc('temporary neighbor list')
	temp_list(:,:) = 0
	do j=1,nat_solute
	  do i=1,max_nbr_range
		if (list14(i,j)) temp_list(i,j) = 1
	  end do
	end do
	if(nat_solute > 0) write(u, '(80i1)') temp_list(1:max_nbr_range, 1:nat_solute)

	write(*, 20) n14nbrs
	write(*, 10, advance='no') 'long-range neighbour list'
	write(u, '(i9,a)') n14long, &
		' = No. of long 1-4 nbrs (>max_nbr_range). pairlist: '
	DO i = 1, n14long
		write(u, '(i7,1x,i7)') list14long(1:2, i)
	enddo
	write(*, 20) n14long

! --- EXCLUSION NEIGHBOUR LIST
	write(*, 10, advance='no') 'neighbour exclusion list'
	write(u, '(i9,a)') nexnbrs, &
		' = No. of exclusions. exclusion list (range=max_nbr_range): '

	! work-around to preserve file format
!	write(u, '(80i1)') listex(1:max_nbr_range, 1:nat_solute)
	temp_list(:,:) = 0
	do j=1,nat_solute
	  do i=1,max_nbr_range
		if (listex(i,j)) temp_list(i,j) = 1
	  end do
	end do
	if(nat_solute > 0) write(u, '(80i1)') temp_list(1:max_nbr_range, 1:nat_solute)
	deallocate(temp_list)

	write(*, 20) nexnbrs
	write(*, 10, advance='no') 'long-range exclusions'
	write(u, '(i9,a)') nexlong, &
		' = No. of long exclusions (>max_nbr_range). pairlist: '
	if(nexlong > 0) write(u, '(i7,1x,i7)') listexlong(1:2, 1:nexlong)
	write(*, 20) nexlong

! --- RESIDUE/MOLECULE BOOKKEEPING
	write(*, 10, advance='no') 'residues'
	if(nres_solute < nres) then
		write(u, '(2i9,a)') nres, nres_solute, ' = No. of residues, No of solute residues. start atoms: '
	else
		write(u, '(i9,a)') nres, ' = No. of residues. start atoms: '
	end if
	if(nres > 0) write(u, '(13(i7,1x))') res(1:nres)%start

	write(u, '(a)') 'Sequence: '
	if(nres > 0) write(u, '(16(a4,1x))') res(1:nres)%name
	write(*, 20) nres

	write(*, 10, advance='no') 'molecules'
	write(u, '(i8,a)') nmol, &
		' = No. of separate molecules. start atoms:'
	if(nmol > 0) write(u, '(13(i7,1x))') istart_mol(1:nmol)
	write(*, 20) nmol

	write(*, 10, advance='no') 'Atom type names'
	write(u, '(i9,a)') max_atyps, &
		' = No. of atom types:'
	if(max_atyps > 0) write(u, '(8(a8,1x))') tac(1:max_atyps)
	write(*, 20) max_atyps

	write(*, 10, advance='no') 'SYBYL type information'
	write(u, '(i9,a)') max_atyps, &
		' = No. of SYBYL atom types:'
	if(max_atyps > 0) write(u, '(13(a5,1x))') SYBYL_atom_type(1:max_atyps)
	write(*, 20) max_atyps

	!solvent type & atom types
	write(u, '(i9,a)') solvent_type, &
		' = solvent type (0=SPC,1=3-atom,2=general)'

	!boundary
	!******PWedited this part
	if( use_PBC ) then !use periodic box
		write(u,'(a8,a)') 'PBC', ' = kind of boundary'
		!PBC-parameters
		write(u,'(3f8.3, a)') boxlength(:), ' = Size of box, x y z'
		write(u,'(3f8.3, a)') boxcentre(:), ' = Centre coordinate of box'
	else !use simulation sphere
		!radii & centres
		write(u, '(2f8.3,a)') rexcl_o, rwat, &
			' = Exclusion, solvent radii'
		write(u, '(3f8.3,a)') xpcent(:), ' = Solute centre'
		write(u, '(3f8.3,a)') xwcent(:), ' = Solvent centre'
		!exluded atoms
		write(*, 10, advance='no') 'excluded atom list'
		write(u, '(2i9,a)') nexats, nexwat, &
			' = No. of excluded atoms (incl. water), no. of excluded waters'
		write(u, fmt='(80l1)') excl(1:nat_pro)
		write(*, 20) nexats
	end if


	close(u)
	write(*,*)
	return

100	write(*,*) 'ERROR: Failed to open topology file ', trim(name)
end subroutine topo_save

!----------------------------------------------------------------------

end module TOPO
