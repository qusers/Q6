! (C) 2014 Uppsala Molekylmekaniska HB, Uppsala, Sweden
! calc_chemscore.f90
! Kajsa Ljungberg 1998-06-11
! Implementation of a scoring function by Eldridge et al.
! JCAMD 11 (1997) 425-445.
! Integration into QCalc5 by Peter Hanspers
!TODO: precision not fixed

module CALC_CHEMSCORE

	use CALC_BASE
	use MASKMANIP
	use TRJ
!	use PRMFILE
	use INDEXER
!	use QATOM
	implicit none

	real(kind=prec)					:: score
	integer						:: frame
	integer						:: i
	character(len=4)	:: trj_type

	type SCORE_PRECALC_TYPE
		character(len=80)	:: chFilename		! name of re file if one is given
		integer						:: iType				! 1 = top calculation, 2 = restart calculation
	end type SCORE_PRECALC_TYPE

	type SCORE_TYPE											! structure tp keep track of scoring data
		character(len=80)		:: chFilename	! either name of restart file, 'top' or empty if trajectory frame
		integer							:: frame			! empty if chReFilename is present
		real(kind=prec)							:: score
		real(kind=prec)							:: h_bonds
		real(kind=prec)							:: metal
		real(kind=prec)							:: lipophil
		real(kind=prec)							:: rot_bond
	end type SCORE_TYPE

!constants
	integer, parameter						::	MAX_MASKS = 10
!module variables
	type(SCORE_PRECALC_TYPE),pointer,private	:: aPrecalc(:)	! array of pre-calculations
	integer, private							:: iPrecalc		! current number of pre-calcs added
	integer, private							:: maxPrecalc	! number of pre-calc elements allocated

	type(SCORE_TYPE),pointer,private		::	aScore(:)				! framewise scoring info
	integer									::	nScores, maxScores		! number of elements in aScore, number of allocated elements of aScore

	type(MASK_TYPE), private, target	::	masks(MAX_MASKS)
	integer, private					::	Nmasks = 0

	type SCORE_COORD_TYPE
		real(kind=prec), pointer				::	xr(:)
		real(kind=prec)							::	xrcm(3)
	end type SCORE_COORD_TYPE
	type(SCORE_COORD_TYPE), private		::	coords(MAX_MASKS)

	type wat
		integer(AI)					::	O, H1, H2		!topology numbers
		real						::	score			!H-bond score with receptor
	end type wat
	
	type donor
		integer(AI)					::	heavy,hydr		!topology numbers for atoms in H-bond donor
	end type donor

	type lipos											!return type for lip_atom
		integer						::	n_lip			!number of lipophilic atoms
		integer						::	n_nonlip		!number of heavy, non-lipophilic atoms
	end type lipos
	
	type bond_pointer									!pointer to the datatype q_bond
		type(q_bond),pointer		::	qb				!needed for array of q_bonds in type q_atom
	end type bond_pointer
	
	type q_atom			
		integer						::	top_nr									! atom number in topology		
		integer						::	n												!number of bonds								
		type(bond_pointer),dimension(1:6)		::	bd		!pointers to the bonds
		logical						::	contact				!if	closer than r to any receptor atom 
		integer						::	at_type				!if lipophilic 2, if other heavy 1, if hydrogen 0 
		integer						::	hybrid				!if included sp3 then 3, if included sp2 then 2, others 0
		logical						::	been_there		!flag used in recursive routines, to avoid eternal loops
		logical,pointer		::	cyclic(:)			!array of logicals, one for each ring, true if part of that ring
		logical						::	active				!same as been_there
	end type q_atom										
	
	type q_bond
		type(q_atom),pointer		::	a,b					!pointers to two bonded atoms
		logical						::	rotatable			!if not part of ring and sp3-sp3 or sp3-sp2 bond
		logical						::	acontact,bcontact	!if a and b in contact with receptor, only used if rotatable
		logical						::	been_there
		integer						::	a_nonlip,b_nonlip	!number of non-lipophilic heavy atoms on a-side and b-side of bond
		integer						::	a_lip,b_lip			!number of lipophilic heavy atoms on a-side and b-side of bond
		logical,pointer				::	cyclic(:)			
		logical						::	active
	end type q_bond
	
	integer(AI), allocatable	::	lph_r(:)	!topology number for all lipophilic atoms in the receptor
	integer(AI), allocatable	::	lph_l(:)	!				''						     the ligand
	type(donor), allocatable	::	hbd_r(:)	!H-bond donor in receptor
	type(donor), allocatable	::	hbd_l(:)	!  ''            ligand
	integer(AI), allocatable	::	hba_r(:)	!H-bond acceptor atoms in receptor 
	integer(AI), allocatable	::	hba_l(:)	!      ''                 ligand
	integer(AI), allocatable	::	met_r(:)	!metal atoms in the receptor
	type(wat), allocatable		::	waters(:)	!water molecules
	
	integer, allocatable		::	pol_con(:)	!connections for potentially polar atoms
												
	
	type(q_atom),private,allocatable,target	::	q_atoms(:)		! atoms in ligand
	type(q_bond),private,allocatable,target	::	q_bonds(:)		! bonds in ligand
	
	!numbers of each of the atom types
	integer						::	nlph_r, nlph_l, nhbd_r, nhbd_l
	integer						::	nhba_r, nhba_l, nmet_r, nwaters,nqbonds
	integer						::	nhbd_prot, nhba_prot
	integer						::	nrings		!number of rings in the ligand
	
	!lists of the atom type categories
	integer, parameter				::	LIPO=1, METAL=2, SULPHUR=3, ACCEPTOR=4
	integer, parameter				::	SP3=5, SP2=6, SP3_N=7, PL_N=8
	integer, parameter				::	HYDROGEN=9, OXYGEN=10, CARBON_LIPO=11
	integer, parameter				::	CARBONYL_CARBON=12, CARBONYL_OXYGEN=13
	integer, parameter				::	NPROPS=13
	character*(*), parameter		::	file_heading(NPROPS) = [ &
		'lipophilic     ', &
		'metal          ', &
		'sulphur        ', &
		'acceptor       ', &
		'sp3            ', &
		'sp2            ', &
		'sp3_nitrogen   ', &
		'planar_nitrogen', &
		'hydrogen       ', &
		'oxygen         ', &
		'carbon_lipo    ', &
		'carbonyl_carbon', &
		'carbonyl_oxygen']
	
	type ATOM_DATA_TYPE
		real						::	radius
		logical						::	prop(NPROPS)
	end type ATOM_DATA_TYPE
	type(ATOM_DATA_TYPE), allocatable ::atom_data(:)
	
	!No record is kept of polar(non H-bonding) atoms
	
!	integer, allocatable		::	iqatom(:)	!one element per atom, 0 if not Q-atom, 
												!else number in iqseq (ligand)
	real,allocatable			::	vdwr(:)		!van der Waals radius of each atom 	

	character*80				::	atom_data_file
	character*80				::  coord_file	
	
	real(kind=prec)				::	hbond_term,metal_term,lipo_term,rot_term	!result from score
	real(kind=prec), parameter		::	dGconst=-5.48_prec, dGhbond=-3.34_prec, dGmetal=-6.03_prec
	real(kind=prec), parameter		::	dGlipo=-0.117_prec, dGrot=2.56_prec
	
	integer, private			::	bDoTopcalc		! boolean flag
	integer, private			::	iRestartCalc=-1	! flag to indicate if doing calcs on restart file or not

	logical, private			::	bUseXIN	= .false.		! boolean flag to indicate use of xin instead of xtop

	integer, private			:: warn						! flag to indicate if any warnings were displayed

contains
	
subroutine score_initialize
	nScores = 0
	maxScores = 16
	allocate(aScore(maxScores))
        calc_xtop(1:3*nat_pro-2:3) = xtop(1:nat_pro)%x
        calc_xtop(2:3*nat_pro-1:3) = xtop(1:nat_pro)%y
        calc_xtop(3:3*nat_pro  :3) = xtop(1:nat_pro)%z
end subroutine score_initialize

subroutine chemscore_finalize
	if(warn.ne.0) then
		write(*,'(a,i3,a)') 'WARNING: There were ', warn, ' warnings.'
	else
		write(*,'(a)') 'ChemScore finished succesfully.'
	end if
end subroutine chemscore_finalize

subroutine log_frame(iFrame)					! logs scoring results for frame
	! args
	integer, intent(in)		:: iFrame			! current frame
	
	! locals
	type(SCORE_TYPE),pointer:: new_aScore(:)	! tmp pointer
	
	nScores = nScores +1
	if(nScores>maxScores) then
		maxScores = maxScores + 16		! allocate 16 elements a time
		allocate(new_aScore(maxScores))	! allocate new mem
		new_aScore(1:nScores-1) = aScore(1:nScores-1)	! copy old
		
		deallocate(aScore)				! return old mem to OS
		aScore => new_aScore				! update pointer
	end if

	! finally store values
	aScore(nScores)%frame		= iFrame	! iFrame may be redundant
	aScore(nScores)%score		= score
	aScore(nScores)%h_bonds		= hbond_term
	aScore(nScores)%metal		= metal_term
	aScore(nScores)%lipophil	= lipo_term
	aScore(nScores)%rot_bond	= rot_term
end subroutine log_frame
	
subroutine score_heading(i)
	integer		::	i

	write(*,100)  'score', 'h-bonds', 'metal', 'lipophilic', 'rot_bond'
100	format(        t29,a5,  t38,a7,    t51,a5,  t58,a10,       t71,a8)

!	write(*,'(a)', advance='no') 'SCORE'
end subroutine score_heading

subroutine calc_hbonds		  !in subroutine calc_scores
	integer		:: i,j
	real(kind=prec)	:: score1,score2

	hbond_term = 0

	do i = 1,nhbd_l																			!ligand h-bond donors scored against
		do j = 1,nhba_r																		!receptor acceptors
			score1 = g1_hb(hbd_l(i)%hydr,hba_r(j))					!distance term
			score2 = g2_hb(hbd_l(i)%heavy,hbd_l(i)%hydr,hba_r(j))   !angle term
			hbond_term = hbond_term+score1*score2
		end do
	end do
	
	do i = 1,nhba_l																			!ligand acceptors scored against
		do j = 1,nhbd_r																		!receptor donors
			score1 = g1_hb(hba_l(i),hbd_r(j)%hydr)
			score2 = g2_hb(hba_l(i),hbd_r(j)%hydr,hbd_r(j)%heavy)
			hbond_term = hbond_term+score1*score2
		end do
	end do

end subroutine calc_hbonds

subroutine calc_lipo				!in subroutine calc_scores
	integer					::i,j	
	real(kind=prec)				::fij

	lipo_term = 0

	do i = 1,nlph_l								!lipophilic atoms in ligand
		do j = 1,nlph_r							!    ''              receptor
			fij = fr_lip(lph_l(i),lph_r(j))
			lipo_term = lipo_term+fij
		end do
	end do
	
end subroutine calc_lipo

subroutine calc_metals				!in subroutine calc_scores
	integer			:: i,j
	real(kind=prec)		:: score

	metal_term = 0
	do i = 1,nmet_r												!metal atoms in receptor
		do j = 1,nhba_l											!H-bond acceptors in ligand		
			score = fr_met(met_r(i),hba_l(j))
			metal_term = metal_term+score
		end do
	end do

end subroutine calc_metals

subroutine calc_rot					!in subroutine calc_scores
	integer			:: i,nfrozen
	real(kind=prec)		:: pa,pb,hrot,sum
	
	nfrozen = 0
	sum = 0

	do i = 1,nqbonds
		if (q_bonds(i)%acontact .and. q_bonds(i)%bcontact) then  !if the bond is rotatable and frozen
				nfrozen = nfrozen + 1
				pa = real(q_bonds(i)%a_nonlip)/(q_bonds(i)%a_nonlip + q_bonds(i)%a_lip) ! % non-lipophilic heavy atoms on a-side
				pb = real(q_bonds(i)%b_nonlip)/(q_bonds(i)%b_nonlip + q_bonds(i)%b_lip)
				sum = sum + ((pa + pb)/2)		!add to sum with weight according to pa and pb
				write(*,*) 'qbond ', i, ', atoms ', q_bonds(i)%a%top_nr, q_bonds(i)%b%top_nr, ', contrib = ', ((pa + pb)/2)
		end if
	end do
	if(nfrozen == 0) then 
		hrot = 0					
	else
		hrot = 1 + (1-1.0/nfrozen)*sum
	end if

	rot_term = hrot
end subroutine calc_rot

subroutine calc_scores
	call set_waters				! does calcs on waters
	call calc_hbonds			! calc contribution from hbond interaction between ligand and receptor
	
	call calc_metals			! calc contrib from metals in rec interacting with h-bond acceptors in ligand
	
	call calc_lipo				! calc lipophilic contrib

	call q_contacts				! set contact parameter =.true. for qatoms in contact with receptor
	call frozen						! checks if rotatable bonds are frozen
	call calc_rot					! calcs contrib from frozen rotatable bonds
end subroutine calc_scores

subroutine count_qlipo		
	!count (non)lipophilic heavy atoms on each side of every bond
	integer					::	i
	type(lipos)				::	lipo

	q_bonds(:)%a_lip = 0		!reset
	q_bonds(:)%b_lip = 0
	q_bonds(:)%a_nonlip = 0
	q_bonds(:)%b_nonlip = 0

	do i = 1,nqbonds 
		if (q_bonds(i)%rotatable) then			!go through all rotatable bonds
			q_bonds(:)%been_there = .false.
			q_atoms(:)%been_there = .false.		!reset

			q_bonds(i)%been_there = .true.
			lipo = lip_atom(q_bonds(i)%a)		!lip_atom recursively counts heavy atoms on a-side
			q_bonds(i)%a_lip = lipo%n_lip
			q_bonds(i)%a_nonlip = lipo%n_nonlip
			lipo = lip_atom(q_bonds(i)%b)		!on b-side
			q_bonds(i)%b_lip = lipo%n_lip
			q_bonds(i)%b_nonlip = lipo%n_nonlip
		end if
	end do
end subroutine count_qlipo


recursive function lip_atom(atom) result(ans)
!arguments
	type(q_atom),pointer			::	atom
		
!local
	type(lipos)						::  lipo
	type(lipos)						::	temp
	integer							::	i
	type(q_bond),pointer			::	bond
	type(lipos)						::  ans

	lipo%n_lip = 0
	lipo%n_nonlip = 0
	atom%been_there = .true.					!flag to avoid cycling in rings
	do i = 1,atom%n								!go through every bond to this atom
		bond => atom%bd(i)%qb							
		if ( .not. bond%been_there ) then		!if the bond has not already been passed
			temp%n_lip = 0
			temp%n_nonlip = 0
			bond%been_there = .true.
												!count heavy atoms on the opposite side
			if (.not. bond%a%been_there ) then	!i.e. the side with an atom which has not been passed
				temp = lip_atom(bond%a)			
			end if
			if (.not. bond%b%been_there ) then 
				temp = lip_atom(bond%b)
			end if	
	
			lipo%n_lip = lipo%n_lip+temp%n_lip  !add result from each bond
			lipo%n_nonlip = lipo%n_nonlip+temp%n_nonlip
		end if
	end do
	
	if (atom%at_type == 2) then					!add atom itself to sum
		lipo%n_lip = lipo%n_lip + 1
	end if
	if (atom%at_type == 1) then
		lipo%n_nonlip = lipo%n_nonlip + 1
	end if

	ans = lipo									!return result to atom which made the call
end function lip_atom
		

subroutine frozen
	! check each rotatable bond to see if it is frozen
	integer					::	i,j
	type(q_bond),pointer	::	bond
	type(q_atom),pointer	::	atom

	q_bonds(:)%acontact = .false.	!reset
	q_bonds(:)%bcontact = .false.
				
	do i = 1,nqbonds										! go through all bonds
		bond => q_bonds(i)								! pointer to bond i
		if (bond%rotatable) then					! only have to check rotatable bonds
			q_bonds(:)%been_there = .false.	! reset
			q_atoms(:)%been_there = .false.
			
			bond%been_there = .true.
			bond%acontact = find_contact(bond%a)	!if any contacts on a-side of rotatable bond
			bond%bcontact = find_contact(bond%b)	!if any contacts on b-side
		end if
	end do
end subroutine frozen

recursive function find_contact(atom) result(ans)	
!argument
	type(q_atom),pointer	:: atom
!local
	integer					::	i
	logical					::	cont
	type(q_bond),pointer	::	bond
	logical					::	ans

	cont = .false.						!no contact found yet
	atom%been_there = .true.

	if (atom%contact) then				!if this atom is in contact there is no need to keep looking
		cont = .true.
	else								!else go through all the bonds to this atom which have not been passed yet
		do i = 1,atom%n	
			bond => atom%bd(i)%qb	
			if ( .not. bond%been_there ) then
				bond%been_there = .true.

				if ( .not. bond%a%been_there ) then		!if not passed atom a (call came from atom b)
					cont = find_contact(bond%a)			!check if contact on this side
				else if ( .not. bond%b%been_there ) then		!if call came from atom a
					cont = find_contact(bond%b)			!check b-side
				end if

				if (cont) then							!if contact found on other side of bond
					exit								!no need to keep looking
				end if
			end if
		end do
	end if

	ans = cont											!return result to atom which made the call
end function find_contact


subroutine get_atom_data
	!---------------- Local variables -----------------!
	integer						::	i, j
	real(kind=prec)						::	rad(natyps)
	real(kind=prec)						::	radius
	logical						::	stat
	logical						::	novdw(natyps)
	character(len=8)			::	tac_in
	!--------------------------------------------------!

	write(*,'(a,a)') 'Reading atom parameters from ', trim(adjustl(atom_data_file))

	allocate(atom_data(natyps))
	do i = 1, natyps
		atom_data(i)%radius = 0.
		atom_data(i)%prop(:) = .false.
	end do

	if(.not. prm_open_section('vdw_radii', atom_data_file)) then
		write(*,*) 'get_atom_data:: Failed to open vdw_radii section of parameter file ', trim(atom_data_file)
		stop 'Parameter file read failure'
	end if

	novdw(:) = .false.
	allocate(vdwr(1:nat_solute))

	do while(prm_get_string_real(tac_in, radius))		! tac_in and radius are returned
		if(.not. index_get(tac_in, i)) then
!			write(*,10) trim(tac_in)
10			format('>>> WARNING: Atom type ',a,' ignored: not in topology.')
		else
			atom_data(i)%radius = radius
		end if
	end do
	
	do i = 1,nat_solute			!!! Go through all atoms and give them a vdW radius.
		vdwr(i) = atom_data(iac(i))%radius	!!! I represents atom number. 
		if (vdwr(i) == zero) novdw(iac(i)) = .true.		
	end do

	do i = 1, natyps
		if (novdw(i)) then
			write(*,20) tac(i)
			warn = warn +1
20		format('WARNING: No vdW radius designated to atom type ', a)
		end if
	end do

	!read atom type categories
	do j = 1, NPROPS
		if(.not. prm_open_section(file_heading(j))) then
			write(*,30) file_heading(j)
30			format('Failed to open ',a,' section of parameter file.')
			stop 'Parameter file read failure'
		end if
		do while(prm_get_line(tac_in))
			if(.not. index_get(tac_in, i)) then
!				write(*,10) trim(tac_in)
			else
				atom_data(i)%prop(j) = .true.
			end if
		end do
	end do
	
end subroutine get_atom_data

subroutine make_tac_index
	! creates index of atom types	
	integer		::	i

	call index_create(natyps)		! natyps is read from file
	do i = 1, natyps
		if(.not. index_add(tac(i), i)) then
			stop '>>>>> ERROR: Could not make atom type index'
		end if
	end do
end subroutine make_tac_index

subroutine q_contacts		
	!check which q-atoms are in contact with receptor, in subroutine calc_scores
	integer				::	i,j,tnr
	real					::	distance
	
	q_atoms(:)%contact = .false.

	do i = 1,nqat																							! loop over heavy q-atoms
		tnr = q_atoms(i)%top_nr
		if (heavy(tnr)) then																		! if heavy 
			do j = 1,nat_solute																		! all protein atoms, not the waters
				if ((iqatom(j) == 0) .and. (heavy(j)) ) then				! if heavy receptor atom (not q-atom)
					distance = dist(tnr,j)														! measure distance
					if (distance <= ( vdwr(tnr)+vdwr(j)+0.5 )) then		! if distance < vdw1+vdw2+0.5
						q_atoms(i)%contact = .true.											! then contact
						exit		!no need to keep looking once one contact has been found
					end if
				end if
			end do
		end if
	end do
end subroutine q_contacts

subroutine q_types			
	! assign q-atom types and hybrids, used in entropy
	! calculation, in subroutine set_ligand
	integer						::	i,nr, a_type


	!set atom type 0:hydrogen, 2:heavy lipophilic, 1:non-lipophilic heavy	
	do i = 1,nqat
		nr = iqseq(i)										! topology number
		q_atoms(i)%top_nr = nr					! used in other subroutines
		a_type = iac(nr)
		if(atom_data(a_type)%prop(LIPO)) then !all C, P, F
			if (pol_con(nr) < 2)	then	!if carbons have >=2 polar connections then it is polar, 62 always have 0
				q_atoms(i)%at_type = 2		!lipophilic
			else
				q_atoms(i)%at_type = 1		!non-lipophilic
			end if
		elseif(atom_data(a_type)%prop(SULPHUR)) then	!sulphur
			if (pol_con(nr) < 1) then		!which is not attached to any polar atom
				q_atoms(i)%at_type = 2		!is lipophilic
			else
				q_atoms(i)%at_type = 1
			end if
		elseif(iaclib(a_type)%mass < 4.) then !hydrogen
			q_atoms(i)%at_type = 0			!not heavy
		else						!all others are non-lipophilic heavy
				q_atoms(i)%at_type = 1
		end if


		!set hybridization
		if(atom_data(a_type)%prop(SP3)) then !sp3
			q_atoms(i)%hybrid = 3
		elseif(atom_data(a_type)%prop(SP3_N)) then !sp3 nitrogens
			if (pol_con(nr) == 21 .or. pol_con(nr) == 31) then	!if bonded to two H and one other, or three H and one other
				q_atoms(i)%hybrid = 0							!then terminal NH2 or NH3 and should not be included	
			else
				q_atoms(i)%hybrid = 3	
			end if 
		elseif(atom_data(a_type)%prop(SP2)) then !sp2
			q_atoms(i)%hybrid = 2
		elseif(atom_data(a_type)%prop(PL_N)) then !planar nitrogens
			if (pol_con(nr) == 21 ) then	!if bonded to two H and one other then terminal NH2 and not included
				q_atoms(i)%hybrid = 0
			else
				q_atoms(i)%hybrid = 2	
			end if 
		else
			q_atoms(i)%hybrid = 0
		end if
	end do

end subroutine q_types
              
subroutine report_ligand
	integer					::i,top,cont

	write(*,100) '#',  'atom', 'type', 'hybrid', 'type', 'bonds', 'contact'
	100	format(   t8,a, t14,a,  t23,a,  t35,a,    t48,a,  t58,a,   t67,a)

	do i = 1,nqat
		top = q_atoms(i)%top_nr
		if (q_atoms(i)%contact) then
			cont = 1
		else
			cont = 0
		end if
		110	format(t1,i8, t10,i8, t23,a8, t33,i8, t45,i8, t55,i8, t66,i8)
!		110	format((1x,i8), t9,(1x,i8), t18,(1x,a8), t27,(1x,i8), t36,(1x,i8), t45,(1x,i8), t54,(1x,i8))
!		110 format(2(1x,i8),1x,a8,4(1x,i8))
!		110 format(7i8)
	
		write(*, 110) i, top, tac(iac(top)), q_atoms(i)%hybrid, q_atoms(i)%at_type, q_atoms(i)%n, cont
!		write(*, 110) i, top, iac(top), q_atoms(i)%hybrid, q_atoms(i)%at_type, q_atoms(i)%n, cont
	end do

	write(*,120) 'lipophilic'
	120 format(1a11)
	130 format(1i7)
	140 format(2i7)
	do i = 1,nlph_l	
		write(*,130) lph_l(i)
	end do
	write(*,120) 'Hb acc'
	do i = 1,nhba_l
		write(*,130) hba_l(i)
	end do
	write(*,120) 'Hb don'
	do i = 1,nhbd_l
		write(*,140) hbd_l(i)%heavy, hbd_l(i)%hydr
	end do
end subroutine report_ligand	

subroutine report_protein

100	format('# Receptor atom summary:')
110	format('# ',a,t20,i6)
	write(*,100)
	write(*,110) 'total', nat_pro-nqat
	write(*,110) 'lipophilic', nlph_r
	write(*,110) 'H-bond donors', nhbd_r
	write(*,110) 'H-bond acceptors', nhba_r
	write(*,110) 'metal ions', nmet_r
	write(*,110) 'water molecules', nwaters
	write(*,*)
end subroutine report_protein		

subroutine report_rings
	integer					::i,j

	write (*,'(/,a,I3)') 'number of rings = ',nrings

	do i = 1,nrings
		write (*,'(a,I3,a)') 'ring ',i,': '
			do j = 1,nqbonds
				if (q_bonds(j)%cyclic(i)) then
					write (*,'(I5,a,I5)') q_bonds(j)%a%top_nr,' ',q_bonds(j)%b%top_nr
				end if
			end do
	end do

end subroutine report_rings

subroutine reset_waters	!in subroutine set_waters
		nhbd_r = nhbd_prot
		nhba_r = nhba_prot
end subroutine reset_waters
		
			


subroutine score_waters
	!calculate hydrogen bonds between each water molecule and receptor !in set_waters
	integer				::i,j,k
	real(kind=prec)			::score1,score2, sum_score


	do i = 1,nwaters
		sum_score = zero
		
		do j = 1,nhba_r								!water hydrogens as donors, scored against receptor acceptors
			score1 = g1_hb(waters(i)%H1,hba_r(j))				!distance term
			score2 = g2_hb(waters(i)%O,waters(i)%H1,hba_r(j))	!angle term
			sum_score = sum_score+score1*score2

			score1 = g1_hb(waters(i)%H2,hba_r(j))
			score2 = g2_hb(waters(i)%O,waters(i)%H2,hba_r(j))	
			sum_score = sum_score+score1*score2
		end do

		do j = 1,nhbd_r		!water oxygen as acceptor, scored against receptor donors
			score1 = g1_hb(waters(i)%O,hbd_r(j)%hydr)
			score2 = g2_hb(waters(i)%O,hbd_r(j)%hydr,hbd_r(j)%heavy)	
			sum_score = sum_score+score1*score2
		end do
		waters(i)%score = sum_score
	end do
end subroutine score_waters

subroutine set_ligand
	call q_types				! assign vdw radii to qatoms
	call set_rings			! detect rings
	call set_rotatable	! set parameter 'rotatable' to apropriate value
	call count_qlipo		! count (non)lipophilic heavy atoms on each side of every bond
end subroutine set_ligand

subroutine set_rings		!in set_ligand
	
	integer					::	ir
	type(q_atom),pointer	::	start_atom
	integer					::	i,j

	nrings = nqbonds-nqat+1						!number of rings in ligand

	if (nrings > 0) then							!if at least one ring
		do i = 1,nqbonds 
			allocate (q_bonds(i)%cyclic(nrings))
			do j = 1,nrings
				q_bonds(i)%cyclic(j) = .false.		!reset	
			end do
		end do
		do i = 1,nqat
			allocate (q_atoms(i)%cyclic(nrings))
			do j = 1,nrings
				q_atoms(i)%cyclic(j) = .false.
			end do
		end do
		
		q_bonds(:)%been_there = .false.		!reset
		q_atoms(:)%been_there = .false.

		start_atom => q_atoms(1)			!does not matter which q-atom search starts from
		ir = 0								
		
		start_atom%been_there = .true.
		do i = 1,start_atom%n				!go through all bonds of the start atom
			if (.not. start_atom%bd(i)%qb%been_there ) then
				call search_ring(start_atom%bd(i)%qb,ir)
			end if
		end do
		
	end if
end subroutine set_rings

recursive subroutine search_ring(bond,ir)	
!arguments
	type(q_bond),pointer			::	bond
	integer							::	ir
!local
	type(q_atom),pointer			::	atom
	integer							::	i

	bond%been_there = .true.
	if ( bond%a%been_there .and. bond%b%been_there ) then   !if both atoms of a bond have been passed, but not the bond 
		call mark_ring(bond,ir)								!then a ring has been found. call mark_ring to find all the atoms of this ring
	end if
	
	atom => bond%a											!pekartilldelning
	if (.not. atom%been_there ) then						!if atom a has not been passed			
		atom%been_there = .true.
		do i = 1,atom%n										!look for rings in all directions which have not been passed 
			if (.not. atom%bd(i)%qb%been_there ) then
				call search_ring(atom%bd(i)%qb,ir)
			end if
		end do	
	else 
		atom => bond%b				
		if (.not. atom%been_there ) then						!if it is atom b which has not been passed	
			atom%been_there = .true.
			do i = 1,atom%n
				if (.not. atom%bd(i)%qb%been_there ) then
					call search_ring(atom%bd(i)%qb,ir)
				end if
			end do	
		end if
	end if
end subroutine search_ring

subroutine mark_ring(bond,ir)	!identify atoms involved in rings
!arguments
		type(q_bond),pointer		::	bond		!starting bond
		integer						::	ir
!local
		logical						::	found		!not used for anything, but the function trace_ring returns a value

		ir = ir+1									!one more ring has been found
		q_bonds(:)%active = .false.					!reset
		q_atoms(:)%active = .false.	
		bond%cyclic(ir) = .true.					!this bond is part of the ring
		bond%active = .true.			
		bond%a%cyclic(ir) = .true.					!atom on a-side also part of the ring
		found = trace_ring(bond%b,ir)				!try to close ring starting at atom on b-side 
end subroutine mark_ring

recursive function trace_ring(atom,ir) result(ans)
!arguments
	type(q_atom),pointer			::	atom
	integer							::	ir
!local	
	integer							::	i
	logical							::	found
	type(q_bond),pointer			::	bond
	logical							::	ans

	found = .false.
	atom%active = .true.							!this atom has been passed trying to close the ring
	do i = 1,atom%n	
		bond => atom%bd(i)%qb						!go through all bonds which could be part of the ring recently found (been_there)
		if (bond%been_there .and. (.not. bond%active)) then	!but have not been passed trying to close the ring (not active)
			bond%active = .true.
			if ( bond%a%been_there .and. (.not. bond%a%active ))then	!a could be part of the ring, and have not checked it yet
				if ( bond%a%cyclic(ir) ) then								!if a is part of the right ring
					found = .true.											!then the ring is closed
				else
					found = trace_ring(bond%a,ir)							!else look in a direction
				end if
			else if ( bond%b%been_there .and. (.not. bond%b%active ))then	!b could be part of the ring, and have not checked it yet
				if ( bond%b%cyclic(ir) ) then								!if b part of the right ring
					found = .true.											!ring closed
				else
					found = trace_ring(bond%b,ir)							!else look in b direction
				end if
			end if

			if (found) then					!if possible to close ring using this bond
				bond%cyclic(ir) = .true.	!this bond is part of the ring
											!if possible to close the ring using bond i
				atom%cyclic(ir) = .true.	!then this atom is also part of the ring
				exit						!no need to keep looking when ring closed
			end if
		end if
	end do
	
	ans = found													!if possible to close ring using this atom
end function trace_ring 

subroutine set_rotatable
	!check which q-bonds are rotatable

	integer					:: i,j

	q_bonds(:)%rotatable = .false.

	do i = 1,nqbonds		!go through all qbonds, rotatable if sp2-sp3 or sp3-sp3
		if  ( ( (q_bonds(i)%a%hybrid == 3) .and. (q_bonds(i)%b%hybrid > 0) ) &
			.or. ( (q_bonds(i)%b%hybrid == 3) .and. (q_bonds(i)%a%hybrid > 0) ) ) then
			q_bonds(i)%rotatable = .true.
			
			do j = 1,nrings			!check if bond is member of any ring
				if (q_bonds(i)%cyclic(j)) then		!then not rotatable
					q_bonds(i)%rotatable = .false.
					exit
				end if
			end do	

		end if
	end do

end subroutine set_rotatable


subroutine set_waters		!in calc_scores
	call reset_waters			! does some reseting of two ints
	call score_waters			! calculate hydrogen bonds between each water molecule and receptor
	call sort_waters			! does some calculation on the waters
end subroutine set_waters


subroutine sort_atoms		!must be done after sort_bonds ,in subroutine start

	integer					::i, a_type
	
	allocate(hba_l(nqat),lph_l(nqat),lph_r(nat_solute),hba_r(nat_solute),met_r(64))
	nlph_l = 0
	nlph_r = 0
	nmet_r = 0
	nhba_l = 0
	nhba_r = 0
	
	
	do i = 1,nat_solute				!go through all atoms, put in right category	
		a_type = iac(i)
		if(atom_data(a_type)%prop(LIPO)) then
		!atoms which may be lipophilic: all carbons, Cl, Br, I which are not ions
			if (pol_con(i) < 2)	then	!If carbon has >=2 polar connections then it is polar (not lipophilic), CL always included
				if (iqatom(i) /= 0) then
					nlph_l = nlph_l+1
					lph_l(nlph_l) = i
				else
					nlph_r = nlph_r+1
					lph_r(nlph_r) = i
				end if
			end if
		elseif(atom_data(a_type)%prop(SULPHUR)) then
			!sulphurs with more than one connection
			if (pol_con(i) < 1) then	!which are not attached to a polar atom are lipophilic
				if (iqatom(i) /= 0) then
					nlph_l = nlph_l+1
					lph_l(nlph_l) = i
				else
					nlph_r = nlph_r+1
					lph_r(nlph_r) = i
				end if
			end if
		elseif(atom_data(a_type)%prop(METAL)) then
			!all metal atoms
				nmet_r = nmet_r+1			!only metals in receptor
				met_r(nmet_r) = i
		elseif(atom_data(a_type)%prop(ACCEPTOR)) then
			! all N, all O, all halogen ions, sulphurs with only one 
			! connection can be Hbond acceptors		
			if (pol_con(i) .le. 2) then	!If N has >2 connections then it is polar, all other atoms have value 0
				if (iqatom(i) /= 0) then
					nhba_l = nhba_l+1
					hba_l(nhba_l) = i
				else
					nhba_r = nhba_r+1
					hba_r(nhba_r) = i
				end if			
			end if
		end if
	end do

	nhba_prot = nhba_r

end subroutine sort_atoms

subroutine sort_bonds		!identify pairs of atoms involved in q-bonds and/or hydrogen bonding !in start 
	integer						::	atom1,atom2,hydrogen_atom
	integer						::	heavy_at,a_type,type1,type2
	integer						::	i,j,k
	logical						::	polar

	!allokera listor

	allocate (q_atoms(nqat))
	allocate (q_bonds(nqat*2))		! worst case?
	allocate (hbd_r(nat_solute),hbd_l(nqat))
	allocate (pol_con(nat_solute))

	nhbd_l = 0
	nhbd_r = 0
	q_atoms(:)%n = 0
	nqbonds = 0
	pol_con(:) = 0

	do i = 1,nbonds_solute		! loop over all bonds between solute atoms
		polar = .false. 
		atom1 = bnd(i)%i		
		atom2 = bnd(i)%j
		
		if(atom_data(iac(atom1))%prop(HYDROGEN)) then		!all polar H	
			polar = .true.				!if one of the atoms is a polar hydrogen, then the bond is polar
			hydrogen_atom = atom1
			heavy_at = atom2
		elseif(atom_data(iac(atom2))%prop(HYDROGEN)) then
			polar = .true.
			hydrogen_atom = atom2
			heavy_at = atom1
		end if

		if (polar) then
			a_type = iac(heavy_at)					!check if the bond matches any of the h-bond categories
			if(atom_data(a_type)%prop(OXYGEN)) then	!all oxygens (not water)	
				if (iqatom(heavy_at) /= 0) then		!if q-atom
					nhbd_l = nhbd_l+1
					hbd_l(nhbd_l)%heavy = heavy_at	
					hbd_l(nhbd_l)%hydr = hydrogen_atom
				else								!else receptor
					nhbd_r = nhbd_r+1
					hbd_r(nhbd_r)%heavy = heavy_at		
					hbd_r(nhbd_r)%hydr = hydrogen_atom
				end if
			elseif(atom_data(a_type)%prop(SP3_N) .or. &
				atom_data(a_type)%prop(PL_N)) then		!all nitrogens			
				if (iqatom(heavy_at) /= 0) then
					nhbd_l = nhbd_l+1
					hbd_l(nhbd_l)%heavy = heavy_at		
					hbd_l(nhbd_l)%hydr = hydrogen_atom
				else
					nhbd_r = nhbd_r+1
					hbd_r(nhbd_r)%heavy = heavy_at		
					hbd_r(nhbd_r)%hydr = hydrogen_atom
				end if
				pol_con(heavy_at) = pol_con(heavy_at)+10		!N bonded to H should not be included in hbond acceptors
			end if
		else										!if not polar hydrogen, count connections to N,S and C
			type1 = iac(atom1)
			type2 = iac(atom2)
			if(atom_data(type1)%prop(SP3_N) .or. &
				atom_data(type1)%prop(PL_N)) then !all nitrogens
				pol_con(atom1) = pol_con(atom1)+1		!count all connections to N
			elseif(atom_data(type1)%prop(CARBON_LIPO)) then
				!all carbons (except CH3 and carbons which might be in 
				!carbonyls or nitriles)
				if(atom_data(type2)%prop(PL_N) .or. &
					atom_data(type2)%prop(ACCEPTOR)) then
					!	case (30:37,40,41,43,44,49,51,52,61)	
					! all polar and H-bonding atoms (but not H, not C, not F) 
					!i.e. all N, all O, all halogen ions, P, F, S with only one connection
					!count polar connections to C
					pol_con(atom1) = pol_con(atom1)+1	
				end if
			elseif(atom_data(type1)%prop(SULPHUR)) then	!sulphur						
				!count polar connections to S with more than one connection 
				if(atom_data(type2)%prop(PL_N) .or. &
					atom_data(type2)%prop(ACCEPTOR)) then
					!and not always bonded to a polar atom (49, sulphate, is always polar and no need to count)
					!case (30:37,40,41,43,44,49,51,52,61)	! all polar and H-bonding atoms (but not H, not C, not F)
						pol_con(atom1) = pol_con(atom1)+1
				end if
			elseif(atom_data(type1)%prop(CARBONYL_CARBON)) then
				!case (21) 
				!sp2 carbons which might be carbonyl, 
				!(carbons which could be part of nitriles should be included)
				if(atom_data(type2)%prop(CARBONYL_OXYGEN)) then
					!case (40,43)	
					!carbonyl O, (nitrile N should come here)
					pol_con(atom1) = pol_con(atom1)+2 !carbonyl always polar
				elseif(atom_data(type2)%prop(PL_N) .or. &
					atom_data(type2)%prop(ACCEPTOR)) then
					!case (30:37,41,44,49,51,52,61)	
					!all other polar atoms, see list above
						pol_con(atom1) = pol_con(atom1)+1 !else only count polar connections
				end if
			end if

			!same for atom 2
			if(atom_data(type2)%prop(SP3_N) .or. &
				atom_data(type2)%prop(PL_N)) then	
				!case (30:37)
				pol_con(atom2) = pol_con(atom2)+1
			elseif(atom_data(type2)%prop(CARBON_LIPO)) then
			!case (20,22,23,25:29)
				!select case (type1)
				if(atom_data(type1)%prop(PL_N) .or. &
					atom_data(type1)%prop(ACCEPTOR)) then
						!all polar atoms
						!case (30:37,40,41,43,44,49,51,52,61)	
						pol_con(atom2) = pol_con(atom2)+1
				end if
			elseif(atom_data(type2)%prop(SULPHUR)) then
				!select case (type1)
				if(atom_data(type1)%prop(PL_N) .or. &
					atom_data(type1)%prop(ACCEPTOR)) then
					!case (30:37,40,41,43,44,49,51,52,61)	! polar atoms
						pol_con(atom2) = pol_con(atom2)+1
				end if
			elseif(atom_data(type2)%prop(CARBONYL_CARBON)) then
				!case (21) !sp2 carbon
				!select case (type1)
				if(atom_data(type1)%prop(CARBONYL_OXYGEN)) then
					!case (40,43)
					pol_con(atom2) = pol_con(atom2)+2 !carbonyl carbons are always polar
				elseif(atom_data(type1)%prop(PL_N) .or. &
					atom_data(type1)%prop(ACCEPTOR)) then
						!case (30:37,41,44,49,51,52,61)
						pol_con(atom2) = pol_con(atom2)+1
				end if
			end if
		end if

		!store all bonds between q-atoms, pointers relating arrays q_atoms and q_bonds 

		if ((iqatom(atom1) /= 0) .and. (iqatom(atom2) /= 0)) then
			j = iqatom(atom1)
			k = iqatom(atom2)  

			nqbonds = nqbonds+1
			q_bonds(nqbonds)%a => q_atoms(j)		
			q_bonds(nqbonds)%b => q_atoms(k)	!pointer to atom in q_atoms

			q_atoms(j)%n = q_atoms(j)%n +1		!number of bonds to atom j
			q_atoms(j)%bd(q_atoms(j)%n)%qb => q_bonds(nqbonds)	!pointer to bond in q_bonds
				
			q_atoms(k)%n = q_atoms(k)%n +1
			q_atoms(k)%bd(q_atoms(k)%n)%qb => q_bonds(nqbonds)
		end if
			
	end do

	!go through bond list once more, count bonds to polar sulphurs

	do i = 1,nbonds_solute
		atom1 = bnd(i)%i		
		atom2 = bnd(i)%j
		type1 = iac(atom1)
		type2 = iac(atom2)

		if(atom_data(type1)%prop(SULPHUR) .and. &
			pol_con(atom1) >= 1)  then 
			!if ((type1 == 50) .and. (pol_con(atom1) .ge. 1))  then  
		!	if atom1 is a polar sulphur
			!select case (type2)
			if(atom_data(type2)%prop(CARBON_LIPO) .or. &
				atom_data(type2)%prop(CARBONYL_CARBON)) then
				!case (20:23,25:29)								!all carbons
				pol_con(atom2) = pol_con(atom2)+1			!then atom2 has one more connection to a polar atom
			end if
		end if
			
		if(atom_data(type2)%prop(SULPHUR) .and. &
			pol_con(atom2) >= 1) then
			!if ((type2 == 50) .and. (pol_con(atom2) .ge. 1))  then 
			!if atom2 is a polar sulphur
			!select case (type1)
			if(atom_data(type1)%prop(CARBON_LIPO) .or.  &
				atom_data(type1)%prop(CARBONYL_CARBON)) then
				!case (20:23,25:29)							!all carbons
				pol_con(atom1) = pol_con(atom1)+1		!then atom1 has one more connection to a polar atom
			end if
		end if
	end do

	nhbd_prot = nhbd_r		!used later to keep track of waters in receptor list
	
end subroutine sort_bonds

subroutine sort_waters	!gors ev bara en gang		!in set_waters
	integer					::i

	do i = 1,nwaters						!all waters
		if (waters(i)%score > 1) then		!if the water is attached to receptor
			nhbd_r = nhbd_r+1
			hbd_r(nhbd_r)%heavy = waters(i)%O
			hbd_r(nhbd_r)%hydr = waters(i)%H1
			
			nhbd_r = nhbd_r+1
			hbd_r(nhbd_r)%heavy = waters(i)%O
			hbd_r(nhbd_r)%hydr = waters(i)%H2

			nhba_r = nhba_r+1
			hba_r(nhba_r) = waters(i)%O	
		end if
	end do

end subroutine sort_waters


subroutine start		
!	character*80				::	top_file, fep_file
	integer(4)					::	u,nat3
	integer						::	i

	call centered_heading('Reading input','-')

	! --- topology file is now read by qcalc
	!write(*,'(a)', advance='no') 'Topology file: '
	!read (*,'(a)') top_file			
	!top_file = 'hapc.top'

	! --- fep file is read by score_add
	!write(*,'(/,a)', advance='no') 'Q-atom (FEP) file: '
	!read (*,'(a)') fep_file
	!fep_file = 'lie.fep'
	
	! --- prm file is read by score_add
	!write(*,'(/,a)', advance='no') 'Parameter file: '
	!read (*,'(a)') atom_data_file
	!atom_data_file = 'eldro.prm'

	! --- topology is loaded by qcalc
	!if(.not. topo_load(top_file, 4.00)) then
	!	write(*,'(a)') '>>>>> ERROR: Failed to load topology.'
	!	stop
	!end if  
	if (solv_atom .ne. 0 ) then
	nwaters = (nat_pro-nat_solute)/solv_atom
	else
	nwaters = 0
	end if
	call make_tac_index			! creates text atom code index
	call get_atom_data			!!! Takes vdw_radii from vdw_file and puts into vector vdwr.
													!!! Also assigns radii to every solute atom
	if(.not. qatom_load_atoms(fep_file)) then
		stop 'Failed to read Q-atom list from FEP file.'
	end if
	
	allocate(iqatom(nat_solute))! Allocate mem for every solute atom
	iqatom(:) = 0
	do i=1,nqat
		iqatom(iqseq(i)) = i	!Set references so that iqatom = 0 if protein and {number} of ligand
	end do

	call sort_bonds
	call sort_atoms
	call store_waters
		
end subroutine start



subroutine store_waters !order of atoms in coordinate list must be O,H,H		!in subroutine start
	
	integer			::i,j
	j = 0
	allocate(waters(1:nwaters))

	do i = nat_solute+1,nat_pro,3
		j = j+1
		waters(j)%O = i
		waters(j)%H1 = i+1
		waters(j)%H2 = i+2
	end do

end subroutine store_waters

real function angle(a,b,c)
!determines the angle a-b-c in degrees, using the cosine theorem

	integer					:: a,b,c	!parameters

	real					::ab_sq,bc_sq,ca_sq, scp !local variables
	real, parameter :: pi = 4.0*atan(1.0)

	ab_sq = distsq(a,b)
	bc_sq = distsq(b,c)
	ca_sq = distsq(c,a)

	scp = (ab_sq+bc_sq-ca_sq)/(2*sqrt(ab_sq*bc_sq))
	!it may happen that scp is very slightly outside [-1;1]
    if ( scp >  1.0 ) scp =  1.0
    if ( scp < -1.0 ) scp = -1.0

	angle = acos(scp)*180.0/pi
end function angle

real function dist (a,b)
	!!! Gives the distance between two atoms.
	
	integer				::	a,b			!topology numbers
	real,dimension(3)	::	delta
	
	if(bUseXIN) then
		delta = xin(3*a-2:3*a)-xin(3*b-2:3*b)
	else
		delta = calc_xtop(3*a-2:3*a)-calc_xtop(3*b-2:3*b)			
	end if
	dist = sqrt(dot_product(delta,delta))

end function dist

real function distsq (a,b)
	!!! Gives the squared distance between two atoms.

	
	integer				::	a,b			
	real,dimension(3)	::	delta

	if(bUseXIN) then
		delta = xin(3*a-2:3*a)-xin(3*b-2:3*b)
	else
		delta = calc_xtop(3*a-2:3*a)-calc_xtop(3*b-2:3*b)
	end if
	distsq = dot_product(delta,delta)

end function distsq

real(kind=prec) function fr_lip(a,b)
	!!! Gives f(rlL) i.e. contribution from one lipophilic - lipophilic pair
	
	integer					::	a,b	!topology numbers of the pair
	
	real(kind=prec)				::	dr	!distance between the two atoms
	real(kind=prec)				::  R1

	dr = dist(a,b)
	R1 = vdwr(a)+vdwr(b)+0.5

	if(dr < (R1+3)) then				!if the atoms interact at all
		if(dr <= R1) then			
			fr_lip = one					!if they are closer than R1 the interaction is given full weight
		else
			fr_lip = one-(dr-R1)/3.0_prec		!else the interaction is given reduced weight
		end if
	else
		fr_lip = zero
	end if
end function fr_lip


real(kind=prec) function fr_met(a,b)
	! Gives f(raM) i.e. contribution from one H-bond acceptor
	! (or H-bond acceptor/donor) - metal pair
	integer				::	a,b	! topology numbers
	real(kind=prec)			::	dr	! distance between the two atoms
	
	dr = dist(a,b)
	if(dr < 2.6_prec) then
		if(dr <= 2.2_prec) then
			fr_met = one
		else
			fr_met = one-(dr-2.2_prec)/0.4_prec
		end if
	else
		fr_met = zero
	end if
end function fr_met

real(kind=prec) function g1_hb(a,b)	
	! Gives g1(dev_dr) i.e. contribution depending on the distance between
	! the H of an H-bond donor and an H-bond acceptor
	
	integer			::	a,b					!	top. nr
	real(kind=prec)		::	dr, dev_dr	! distance between the two atoms, deviation from the
															! ideal value 1.85
	dr = dist(a,b)
	dev_dr = q_abs(1.85_prec-dr)

	if(dev_dr < 0.65_prec) then		!if atoms are interacting at all
		if(dev_dr <= 0.25_prec) then
			g1_hb = one
		else
			g1_hb = one-(dev_dr-0.25_prec)/0.4_prec
		end if
	else
		g1_hb = zero
	end if
end function g1_hb

real(kind=prec) function g2_hb(a,b,c)
	!!! Gives g2(dev_ang) i.e. contribution depending on the H-bond angle at the hydrogen atom	
	integer		:: a,b,c			!top.nr, the hydrogen must be atom b
	real(kind=prec) :: ang, dev_ang	!H-bond angle, deviation from ideal value of 180 degrees
	
	ang = angle(a,b,c)
	dev_ang = q_abs(180.0_prec-ang)

	if(dev_ang < 80.0_prec) then
		if(dev_ang <= 30.0_prec) then
			g2_hb = one
		else
			g2_hb = one-(dev_ang-30.0_prec)/50.0_prec
		end if
	else
		g2_hb = zero
	end if
end function g2_hb

subroutine score_precalc
	! 1. Prepares ligand and other stuff for calculations
	! 2. Does calculations on topology or restart file before frame by frame calculations start
	
	integer							:: i
	character(len=200)	:: chBuf

	call start						! Prepares coordinates etc
	call set_ligand				! Prepare ligand
	call report_rings			! output results
	call report_ligand		! output results
	call report_protein		! output results

	if(bDoTopcalc.eq.0) return	! abort if no top-calc is to be output

	! continue scoring topology
	call calc_scores			! does the calcs
	score = dGconst+hbond_term*dGhbond+metal_term*dGmetal+lipo_term*dGlipo+rot_term*dGrot

	write(*,*)					!two empty lines for gnuplot as a data set separator
	write(*,*)
	call centered_heading('Scoring results','=')
	write(*,100) 'file', 'frame', 'score', 'h-bonds', 'metal', 'lipophilic', 'rot_bond'
	write(*,110) 'topology', frame, score, hbond_term, metal_term, lipo_term, rot_term

100	format(a,     t18,a5, t29,a5,  t38,a7,   t51,a5,   t58,a10,  t71,a8)
110	format(a,     t19,i4, t25,f9.2,t36,f9.2, t47,f9.2, t59,f9.2, t70,f9.2)
end subroutine score_precalc

integer function score_add(desc)
	!arguments
	character(*)				:: desc
	integer							:: ats, iRe
	! locals	
	character(len=400)	:: chInFile		! long buffers to meet long paths
	character(len=400)	:: chBuf

	if(Nmasks == MAX_MASKS) then
		write(*,10) MAX_MASKS
		return
	end if
10	format('Sorry, the maximum number of score calculations is ',i2)

	!add a new score mask
	Nmasks = Nmasks + 1
	call mask_initialize(masks(Nmasks))
	ats =  maskmanip_make(masks(Nmasks))

	!discard if no atoms in mask
	if(ats == 0) then
		call mask_finalize(masks(Nmasks))
		Nmasks = Nmasks - 1
		score_add = 0
		return
	end if
	
	chBuf = ''
	do while((trim(chBuf).ne.'yes').and.(trim(chBuf).ne.'no'))
		call getlin(chBuf,'Score initial topology? (yes/no)')
		call locase(chBuf)
	end do

	if(chBuf.eq.'yes') then
		bDoTopcalc = 1
	else
		bDoTopCalc = 0
	end if
		
	call getlin(fep_file,'Q-atom (FEP) file: ')
	write(chBuf,'(a,a,a)') 'Parameter file (blank to use ', trim(prm_file), '): '
	call getlin(atom_data_file, trim(chBuf))

	if(len(trim(atom_data_file)).eq.0) then
		atom_data_file = prm_file
	end if

	score_add = Nmasks
	write(desc, 20) nat_pro
20	format('Score calculation using ',i6,' atoms')

end function score_add

subroutine score_mean
	! calculate and print mean score for all contributions
	! locals
	real(kind=prec) :: m_score, m_hbond, m_metal, m_lipo, m_rot

	do i = 1,nScores
		m_score = m_score + aScore(i)%score
		m_hbond = m_hbond + aScore(i)%h_bonds
		m_metal = m_metal + aScore(i)%metal
		m_lipo  = m_lipo  + aScore(i)%lipophil
		m_rot   = m_rot   + aScore(i)%rot_bond
	end do

	m_score = m_score / nScores
	m_hbond = m_hbond / nScores
	m_metal = m_metal / nScores
	m_lipo  = m_lipo  / nScores
	m_rot   = m_rot   / nScores
	
	! output mean values
	write(*,100) 'processed frames: ', nScores,  'score', 'h-bonds', 'metal', 'lipophilic', 'rot_bond'
100	format(       t1,a,                t21,i4,    t29,a5,  t38,a7,    t51,a5,   t58,a10,     t71,a8)

	write(*,110) 'mean value', m_score,   m_hbond,    m_metal,   m_lipo,    m_rot
110	format(       t1,a,        t25,f9.2,  t36, f9.2,  t47,f9.2,  t59,f9.2,  t70,f9.2)
	write(*,111) '------------------------------------------------------------------------------'                                                                     
111	format(t1,a)	
	m_score = zero
	m_hbond = zero
	m_metal = zero
	m_lipo  = zero
	m_rot   = zero
	nScores = 0

end subroutine score_mean

subroutine score_calc(iCalc, iFrame)		! calc topmost routine
	!arguments
	integer, intent(in)	::	iCalc			! calculation index, not used
	integer, intent(in)	::	iFrame			! frame index

	! <begin fulhack>
	! added flag to have calc routines use xin instead of xtop. not a pretty solution but better than xtop = xin <end fulhack>
	bUseXIN = .true.
	!call start
	!call set_ligand
	call calc_scores

	score = dGconst+hbond_term*dGhbond+metal_term*dGmetal+lipo_term*dGlipo+rot_term*dGrot

!	if((iFrame.eq.0).and.(iRestartCalc.eq.-1)) then		! first restart calc, no calcs have been stored
!		iRestartCalc = 1
!		write(*,*) ' begin fresh restart'
!	elseif((iFrame.eq.0).and.(iRestartCalc.eq.0)) then	! begin restart calc, last calc was trajectory calc so calc mean of that
!		call score_calc_mean							! calc mean of traj scores
!		nScores = 0										! reset frame count
!		write(*,*) ' begin restart'
!		iRestartCalc = 1
!		call log_frame(iFrame)
!	elseif((iFrame.eq.0).and.(iRestartCalc.eq.0)) then	! continued restart calc
!		write(*,*) ' cont restart'
!		call log_frame(iFrame)
!	elseif((iFrame.ne.0).and.(iRestartCalc.eq.1)) then	! begin traj calc, last was restart so calc mean of that
!		write(*,*) ' begin traj'
!		call score_calc_mean
!		nScores = 0
!		iRestartCalc = 0
!		call log_frame(iFrame)
!	elseif((iFrame.ne.0).and.(iRestartCalc.eq.0)) then	! continued traj calc
!		write(*,*) ' cont traj'
!		call log_frame(iFrame)
!	else
!		write(*,*) ' >>>> ERROR: Unrecognized calculation type in score_calc'
!	end if

	call log_frame(iFrame)

	!write(110,*) score, hbond_term, metal_term, lipo_term, rot_term
!	write(100,*) score, hbond_term, metal_term, lipo_term, rot_term
 	!write(*,100, advance='no') score
100	format(f10.3)

	write(*, 110) score,    hbond_term, metal_term, lipo_term, rot_term
110	format(       t25,f9.2, t36,f9.2,   t47,f9.2,   t59,f9.2,  t70,f9.2)

end subroutine score_calc

end module CALC_CHEMSCORE
