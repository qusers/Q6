! (C) 2014 Uppsala Molekylmekaniska HB, Uppsala, Sweden
! calc_xscore.f90
! by Peter Hanspers & Martin Nervall
! Qcalc trajectory analysis main program
!TODO: precision not fixed

module CALC_XSCORE
	!	v.2003.11, Structure-Based Drug Design Toolkits, developed by Dr. Renxiao Wang
	! F90 implementation and QCalc5 integration by Peter Hanspers
	
	! Most of the structure from the original xscore source code is ported to f90.
	! Before scoring, protein and ligand data structures are translated to xscore format.
  
	use CALC_BASE
	use MASKMANIP
	use TRJ
	use	TOPO
!	use PRMFILE
	use INDEXER
!	use QATOM
	use MISC

	implicit none

	type tXMolecule
		integer ::  xtool_format
		integer ::  id
		integer ::  valid
		character(len=256) :: name
		character(len=256) :: formula
		real(kind=prec) ::  weight
		integer ::  num_hb_atom
		integer(2) ::  num_rotor
		real(kind=prec) ::  logp
		real(kind=prec) ::  surface,nsur,psur
		real(kind=prec) ::  volume

		integer ::  num_atom
		type(tXAtom),pointer :: atom(:)

		integer ::  num_bond
		type(tXBond),pointer ::  bond(:)

		integer ::  num_subst
		integer ::  num_feature
		integer ::  num_set
		character(len=256) :: mol_type
		character(len=256) :: charge_type

		integer ::  num_ring
		type(tXRing),pointer		 :: ring(:)			
		integer									 :: ring_count	! actual count

		type(tXDot),pointer		:: vol_dot(:)
		integer								:: num_vol_dot
		type(tXDot),pointer		:: sur_dot(:)
		integer								:: num_sur_dot
	end type

	type tXProtein	! full version
		integer :: xtool_format
		character(len=256)	:: name
		real(kind=prec) ::  surface,bnsur,bpsur

		integer	::  num_atom
		type(tXAtom),pointer :: atom(:)
		integer ::  num_bond
		type(tXBond),pointer ::  bond(:)

		integer :: num_chain

		integer :: num_ring
		type(tXRing), pointer	:: ring(:)
	end type tXProtein

	! ATOMIC BINDING SCORE
	type tABS    
		 real(kind=prec) ::  pkd1,pkd2,pkd3
		 real(kind=prec) ::  vdw,hb,hm,hp,hs,rt
		 real(kind=prec) ::  score
	end type tABS

	type tXLigand
		real(kind=prec) ::  bind_score,chem_score
		real(kind=prec) ::  vdw,sb,hb,hp,hm,hs,ar,rt,pmf,uhb,bnsur,bpsur
		real(kind=prec) ::  pkd1,pkd2,pkd3
	
		integer	:: cofactor_offset	! index of first cofactor atom in protein
	
		type(tABS), pointer ::  abs_inf(:)

		type(tXMolecule)	:: mol		! member inherited from molecule class
	end type tXLigand

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
		logical,pointer		::	cyclic(:)			!array of logicals, one for each ring, 1 if part of that ring
		logical						::	active				!same as been_there
	end type q_atom										

	type q_bond
		type(q_atom),pointer		::	a,b					!pointers to two bonded atoms
		logical						::	rotatable			!if not part of ring and sp3-sp3 or sp3-sp2 bond
		logical						::	acontact,bcontact	!if a and b in contact with receptor, only used if rotatable
		logical						::	been_there
		integer						::	a_nonlip,b_nonlip	!number of non-lipophilic heavy atoms on a-side and b-side of bond
		integer						::	a_lip,b_lip			!number of lipophilic heavy atoms on a-side and b-side of bond
		logical,pointer		::	cyclic(:)			
		logical						::	active

		integer						:: cod
	end type q_bond

	type tAtomStruct
		character(len=10)::  name
		character(len=10)::  ttype
		character(len=10)::  xtype
		real(kind=prec) ::  r
		real(kind=prec) ::  eps
		real(kind=prec) ::  q
		character(len=3)::  hb
		real(kind=prec) ::  logp
		real(kind=prec) ::  solv
		integer ::  ring
		character(len=3)::  pmftype
	end type tAtomStruct

	type tBondStruct
		 character(len=10)::  atom1
		 character(len=10)::  atom2
		 character(len=3)::  ttype
	end type tBondStruct 
	
	type tAtom_Def
		 character(len=10)::  ttype
		 real(kind=prec) ::  weight
		 real(kind=prec) ::  r
		 real(kind=prec) ::  eps
		 real(kind=prec) ::  q
		 character(len=3)::  hb
	end type tAtom_Def
	
	type tXatom_Def
		 character(len=20) ::  ttype
		 real(kind=prec) ::  logp
	end type tXatom_Def
	
	type tBond_Def
		 character(len=10)::  atom_1
		 character(len=10)::  atom_2
		 character(len=3)::  ttype
		 real(kind=prec) ::  length
	end type tBond_Def
	
	type tTors_Def
		 character(len=10)::  atom_1
		 character(len=10)::  atom_2
		 character(len=10)::  atom_3
		 character(len=10)::  atom_4
		 character(len=3)::  ttype
		 real(kind=prec) ::  V         ! twisting force constant
		 integer ::  n     ! periodicity
		 integer ::  S     ! sign of torsion angle type
	end type tTors_Def
	
	type tPMF_Def
		 character(len=3)::  ltype
		 character(len=3)::  ptype
		 real(kind=prec),dimension(0:59) ::  d		! distance
		 real(kind=prec),dimension(0:59) ::  p		! PMF potential
	end type tPMF_Def
	
	type tHB_Def
		 character(len=10)::  ttype
		 real(kind=prec) ::  low_cutoff_d, high_cutoff_d, step_d
		 real(kind=prec) ::  low_cutoff_a1, high_cutoff_a1, step_a1
		 real(kind=prec) ::  low_cutoff_a2, high_cutoff_a2, step_a2
		 integer ::  num_bin_d, num_bin_a1, num_bin_a2, num_bin_total
		 real(kind=prec), pointer :: pmf(:)
	end type
	
	type tWPMF_Def
		 character(len=10)::  ttype
		 real(kind=prec) ::  low_cutoff, high_cutoff, step
		 integer ::  num_bin
		 real(kind=prec),pointer :: pmf(:)
	end type
	
	type tResidue_Def
		character(len=10)::  name
		integer ::  num_atom, num_bond
		type(tAtomStruct),dimension(0:49) :: atom
				
		type(tBondStruct),dimension(0:49) :: bond
	end type tResidue_Def

	type tXFF	! full version
		integer ::  num_restype
 		type(tResidue_Def),pointer :: residue(:)

		integer ::  num_atomtype
		type(tAtom_Def),pointer :: atom(:)

		integer ::  num_xatomtype
		type(tXatom_Def),pointer :: xatom(:)

		integer ::  num_bondtype
		type(tBond_Def),pointer :: bond(:)

		integer ::  num_torstype
		type(tTors_Def),pointer :: torsion(:)

		integer ::  num_pmftype
		type(tPMF_Def),pointer :: pmf(:)

		integer ::  num_sdot_type		! pre-calculated surface dots
		type(tXDotSet),pointer :: sdot(:)

		integer ::  num_vdot_type		! pre-calculated volume dots

		integer ::  num_hbtype
		type(tHB_Def),pointer :: hb_pmf(:)

		integer ::  num_wpmftype
		type(tWPMF_Def),pointer :: wpmf(:)
	end type tXFF

	type tXInput
		integer						::	num_method
		character(len=10) ::  apply_hpscore
		character(len=10) ::  apply_hmscore
		character(len=10) ::  apply_hsscore
		character(len=10) ::  apply_pmfscore
		character(len=20) ::  show_abs
		character(len=10)	::	show_total

		character(len=3)	::	show_ligand
		character(len=3)	::	show_protein
		character(len=3)	::	show_cofactor
		character(len=3)	::  show_bonds

		character(len=256)	:: residue_def
		character(len=256)	:: atom_def
		character(len=256)	:: logp_def
		character(len=256)	:: surface_def
		character(len=256)	:: atom_translations

		integer						::	qatom

		real(kind=prec) ::  hpscore_cvdw
		real(kind=prec) ::  hpscore_chb
		real(kind=prec) ::  hpscore_chp
		real(kind=prec) ::  hpscore_crt
		real(kind=prec) ::  hpscore_c0

		real(kind=prec) ::  hmscore_cvdw
		real(kind=prec) ::  hmscore_chb
		real(kind=prec) ::  hmscore_chm
		real(kind=prec) ::  hmscore_crt
		real(kind=prec) ::  hmscore_c0

		real(kind=prec) ::  hsscore_cvdw
		real(kind=prec) ::  hsscore_chb
		real(kind=prec) ::  hsscore_chs
		real(kind=prec) ::  hsscore_crt
		real(kind=prec) ::  hsscore_c0
	end type tXInput

	type tLogPFactor
		character(len=30) :: symbol
		real(kind=prec) ::  num 
		real(kind=prec) ::  coeff 
	end type tLogPFactor

	type tXAtom
		integer							:: topindex			!  topology index
		integer							:: id  	        !* atom id 
		integer							:: valid				!* valid indicator
		integer							:: mask					!  mask flag
		character(len=10)		:: name         !  atom name
		character(len=10)		:: ttype				!* basic atom type
		character(len=10)		:: xtype        !* xtool atom type
		character(len=10)		:: type2        !  another atom type when necessary
		character(len=10)		:: residue      !* residue name
		character(len=10)		:: old_residue	!  old residue name, helps debugging force field translations
		character(len=10)		:: res_id       !  residue id: must be a string for PDB
		integer							:: atom_res_id	!	 atom number within residue
		integer							:: iscofactor		!  flag indicating if atom is part of a cofactor (1) or not (0)
		character						:: chain        !  chain label
		real(kind=prec),dimension(0:2)	:: coor					!* coordinates
		real(kind=prec),dimension(0:2)	:: root					!* HB root's coordinates 
		real(kind=prec)								:: weight       !* atomic weight
		real(kind=prec)								:: r            !* vdw radius
		real(kind=prec)								:: eps          !* vdw epsilon value
		real(kind=prec)								:: q						!* partial atomic charge
		real(kind=prec)								:: XR						!  X radius
		real(kind=prec)								:: logp         !* atomic hydrophobic scale
		real(kind=prec)								:: solv					!* atomic solvation parameter
		character(len=3)		:: hb           !* HB property
		real(kind=prec)								:: occupancy		!  occupancy probability, for protein atoms
		real(kind=prec)								:: bfactor			!  B-factor, for protein atoms
		real(kind=prec)								:: score        !  atomic binding score
		integer							:: ring					!  ring indicator: 1=normal 2=aromatic
		integer							:: origin				!  atom origin indicator: 1=ligand 2=protein 
		integer							:: part					!  component indicator: 
																				!	 for protein atoms: 1=ATOM 2=HETATM
		integer							:: num_neib			!  number of neighboring atoms
		integer							:: num_nonh			!  number of non-H neighboring atoms
		integer,dimension(0:6) :: neib			!  ID of neighboring atoms
		integer,dimension(0:6) :: bond			!  ID of neighboring bonds
		integer							:: temp					!  for misc uses, e.g. cofactor merge
	end type tXAtom

	type tXWater
		integer							:: id      	!  atom id
		integer							:: valid  !  valid indicator
		character(len=10)		:: name    !  atom name
		character(len=10)		:: type				!	atom type
		character(len=10)		:: xtype	!  another atom type
		character(len=10)		:: residue       !  residue name
		real(kind=prec),dimension(0:3) :: coor    !  coordinates
		real(kind=prec)								:: r    !  vdw radius
		real(kind=prec)								:: eps        !  vdw epsilon value
		real(kind=prec) 								:: q    !  partial atomic charge
		real(kind=prec) 								:: logp       !  atomic hydrophobic scale
		character(len=3)		:: hb       !  HB property
		real(kind=prec) 								:: depth		!  buried depth
		real(kind=prec) 								:: score		!  score value
	end type

	type tXBond	!  full version
		integer :: id
		integer :: valid  !  valid indicator
		integer :: atom_1       !  ID of the atom_1
		integer :: atom_2       !  ID of the atom_2
		character(len=3) :: ttype     !  bond type
		integer :: part   !  ID of the component
		integer :: ring   !  ring indicator
		real(kind=prec) :: length     !  bond length
		integer :: num_neib     !  number of neighboring bonds
		integer,dimension(0:5) :: neib     !  ID of neighboring bonds
	end type

	type tXTorsion	!  full version
		type(tXAtom) atom_1
		type(tXAtom) atom_2
		type(tXAtom) atom_3
		type(tXAtom) atom_4
		character(len=3) :: ttype     !  type of the torsion between 2-3
		integer :: angle  !  torsion angle, in degree
		real(kind=prec) :: V    !  potential barrier
		integer :: n      !  periodicity
		integer :: S      !  sign
		real(kind=prec) :: e    !  torsion energy
	end type

	type tXGroup	!  full version
		integer :: valid
		integer :: num_neib    	!  number of neighboring atoms
		integer :: num_nonh    	!  number of neighboring non-h atoms
		integer :: num_h       	!  number of neighboring hydrogen atoms
		integer :: num_hetero  	!  number of neighboring heteroatoms
		integer :: num_pi      	!  number of neighboring pi atoms
		integer :: num_car     	!  number of neighboring C.ar
		integer :: num_nar     	!  number of neighboring N.ar
		integer :: num_db	!  number of double bonds
		integer :: num_tb	!  number of triple bonds
		integer :: db_type      !  double bond type
		integer :: amide	!  indicator for amide group
		type(tXAtom) center      !  center atom of the group
		type(tXAtom),dimension(0:6) :: neib   !  neighboring atoms
		type(tXBond),dimension(0:6) :: bond   !  bonds
	end type

	type tXHBond
		integer :: valid
		integer :: ttype	!  1: donor=latom, acceptor=patom
											!  2: donor=patom, acceptor=latom
											!  3: donor=metal, acceptor=latom/patom
											!  4: donor=latom, acceptor=latom
											!  5: donor=patom, acceptor=patom
											!  0: invalid, no H-bond
		integer :: d_type	!  donor type: 1=straight 2=angled
		integer :: a_type	!  acceptor type: 1=straight 2=angled
		integer :: sb			!  flag for neutral HB or SB

		integer :: latom	!  id of the ligand atom
		integer :: patom	!  id of the protein atom

		type(tXAtom) H,D,A		! not pointers to atoms but real(kind=prec) atoms

		real(kind=prec) :: dd	!  D-A distance
		real(kind=prec) :: a0	!  D-H-A angle
		real(kind=prec) :: a1	!  DR-D-A angle
		real(kind=prec) :: a2	!  D-A-AR angle

		real(kind=prec) :: score	!  strength of this HBond
	end type

	type tXRing
		integer :: valid	!  0 = invalid 1 = valid
		integer :: ttype	!  1 = normal 2 = aromatic

		integer :: num_member
		integer, pointer					:: atom_id(:)			! worst case
		integer										:: atom_count = 0	! actual count
		integer										:: atom_max = 0		! size of allocated mem
		integer, pointer					:: bond_id(:)			! worst case
		integer										:: bond_count = 0	! actual count
		integer										:: bond_max = 0		! size of allocated mem

		real(kind=prec),dimension(0:3) :: centroid
	end type

	type tXDot
		integer							:: valid		!  status indicator
		character(len=10)		:: ttype		!  type
		real(kind=prec),dimension(0:3) :: coor			!  coordinates
		real(kind=prec)								:: unit			!  contribution, in either A^2 or A^3
		real(kind=prec)								:: score		!  score on this dot
	end type

	type tXDotSet 
		integer							:: num_dot = 0	! current number of dots
		integer							:: max_dot = 0	! number of dots for which mem is allocated
		type(tXDot),pointer	:: dot(:)				! array of dots
		real(kind=prec)								:: r
		character(len=10)		:: ttype
		real(kind=prec)								:: unit					!  default contribution of each dot to total
		real(kind=prec)								:: total				!  total volume or surface
	end type

	type tXResidue
		integer valid
		character(len=10) :: name		!  for PDB files it is a 3-letter string
		character :: chain
		character(len=10) :: id

		integer :: num_atom
	end type

	type tXChain
		integer valid
		character :: label

		integer :: length
	end type

	type tXScore
		integer	:: frame
		real(kind=prec)		:: total
		real(kind=prec)		:: vdw
		real(kind=prec)		:: hb
		real(kind=prec)		:: rt
		real(kind=prec)		:: hp
		real(kind=prec)		:: hs
		real(kind=prec)		:: hm
		real(kind=prec)		:: score
	end type tXScore

	type donor
		integer(AI)			::	heavy,hydr		!topology numbers for atoms in H-bond donor
	end type donor

	type wat
		integer(AI)			::	O, H1, H2		!topology numbers
		real(kind=prec)						::	score			!H-bond score with receptor
	end type wat

	type ATOM_DATA_TYPE
		real(kind=prec)						::	radius
	end type ATOM_DATA_TYPE
	type(ATOM_DATA_TYPE), allocatable ::atom_data(:)

	type ATOM_TYPE_CONVERSION
		character(len=20)	:: xscore, other
	end type

	type RESIDUE_NAME_CONVERSION
		character(len=8)	:: xscore, other
	end type

	type PATOM_NAME_CONVERSION
		character(len=20)					:: res									! residue name
		character(len=20),pointer	:: qname(:), xname(:)		! atom name table
		integer										:: natom								! #atoms in residue
	end type

	type heavybond
		integer						:: a,b					! topology indecies of bond atoms
	end type heavybond

	character*80		:: atom_data_file
	character*80		:: coord_file	

	type(q_atom),private,allocatable,target	::	q_atoms(:)		! atoms in ligand
	type(q_bond),private,allocatable,target	::	q_bonds(:)		! bonds in ligand
!	integer, allocatable					::	iqatom(:)		!one element per atom, 0 if not Q-atom, 
	
	type(q_atom),private,allocatable,target	::	aHQ(:)			! heavy atoms in ligand
	type(q_bond),private,allocatable,target	::	aHB(:)			! heavy bonds in ligand
																!else number in iqseq (ligand)
	integer	:: numHQ		! number of heavy q-atoms
	integer	:: nHB		! number of bonds between heavy q-atoms

	integer(AI), allocatable	::	lph_r(:)	!topology number for all lipophilic atoms in the receptor
	integer(AI), allocatable	::	lph_l(:)	!				''						     the ligand
	type(donor), allocatable	::	hbd_r(:)	!H-bond donor in receptor
	type(donor), allocatable	::	hbd_l(:)	!  ''            ligand
	integer(AI), allocatable	::	hba_r(:)	!H-bond acceptor atoms in receptor 
	integer(AI), allocatable	::	hba_l(:)	!      ''                 ligand
	integer(AI), allocatable	::	met_r(:)	!metal atoms in the receptor
	type(wat), allocatable		::	waters(:)	!water molecules
	
	integer, allocatable			::	pol_con(:)!connections for potentially polar atoms

	!numbers of each of the atom types
	integer						::	nlph_r, nlph_l, nhbd_r, nhbd_l
	integer						::	nhba_r, nhba_l, nmet_r, nwaters,nqbonds
	integer						::	nhbd_prot, nhba_prot
	integer						::	nrings		!number of rings in the ligand

	! System
	type(tXProtein),public							:: protein
	type(tXLigand), public							:: ligand
	type(tXLigand), public,allocatable	:: cofactor(:)
	type(tXFF),public										:: ff
	integer,pointer											:: top2prot(:)				! topology atom index -> protein atom index translation matrix, used in bond translations

	! Mask
	integer, public,parameter						::	MAX_MASKS = 10
	type(MASK_TYPE), public, target			::	masks(MAX_MASKS)
	integer, public											::	Nmasks = 0

	! Results
	type(tXScore), pointer			:: xscores(:)		! xscoring stats for each frame
	integer											:: nXScores, maxXScores

	! Input 
	integer						:: bDoTopCalc		! flag to indicate whether to score initial topology or not
	type(tXInput)			:: input
	character(len=20)	:: chTranslationKey
	character(len=80)	:: chInput
	character(len=80),allocatable	:: cofactor_def(:)
	integer												:: num_cofactor

	! Conversion
	type(ATOM_TYPE_CONVERSION), allocatable	:: atomtypetable(:)
	integer																	:: num_atomtype
	type(PATOM_NAME_CONVERSION), pointer		:: patomnametable(:)
	integer																	:: num_patomrestype
	type(RESIDUE_NAME_CONVERSION), pointer	:: residuetable(:)
	integer																	:: num_residue

	! Consts, params and misc
	integer :: MAX_BOND_NEIB = 10
	integer :: MAX_ATOM_NEIB = 6

	real(kind=prec) ::  LOGP_HYDROPHOBIC_CARBON = 0.211_prec
	real(kind=prec) ::  LOGP_INTERNAL_HBOND = 0.429_prec
	real(kind=prec) ::  LOGP_HALOGEN_PAIR = 0.137_prec
	real(kind=prec) ::  LOGP_NAR_PAIR = 0.485_prec
	real(kind=prec) ::  LOGP_O3_PAIR = -0.268_prec
	real(kind=prec) ::  LOGP_ACCEPTOR_PAIR = 0.580_prec
	real(kind=prec) ::  LOGP_AMINO_ACID = -2.166_prec
	real(kind=prec) ::  LOGP_SALICYLIC_ACID = 0.554_prec
	real(kind=prec) ::  LOGP_SULFONIC_ACID = -0.501_prec
	type(tLogPFactor), dimension(0:9),public	::	logp_factor

	real(kind=prec) ::  LARGE = 1.0e+6_prec
	real(kind=prec) ::  DIST_CUTOFF = 8.0_prec
	real(kind=prec) ::  WATER_R = 1.4_prec
	real(kind=prec) ::  POCKET_DEPTH = 4.0_prec
	real(kind=prec) ::  LAYER_DEPTH = 3.0_prec

	integer, parameter		::	max_atlib = 500

	integer,public				:: err								! public allocation error indicator
	integer, private			:: warn								! flag to inidacte if warnings were displayd

contains

subroutine xscore_initialize
	nXScores = 0
	maxXScores = 16
	allocate(xscores(maxXScores))
        calc_xtop(1:3*nat_pro-2:3) = xtop(1:nat_pro)%x
        calc_xtop(2:3*nat_pro-1:3) = xtop(1:nat_pro)%y
        calc_xtop(3:3*nat_pro  :3) = xtop(1:nat_pro)%z
end subroutine xscore_initialize

subroutine xscore_finalize
	if(warn.ne.0) then
		write(*,'(a,i3,a)') 'WARNING: There were ', warn, ' warnings.'
	else
		write(*,'(a)') 'X-Score finished succesfully.'
	end if
end subroutine xscore_finalize

integer function xscore_add(desc)
	!arguments
	character(*),intent(out)		:: desc
	! locals
	character(len=400)			:: chBuf
	character(len=80)				:: chSmallBuf
	integer									:: i,ats

	if(Nmasks == MAX_MASKS) then
		write(*,1) MAX_MASKS
		return
	end if
1 format('Sorry, the maximum number of score calculations is ',i2)
	
	!add a new score mask
	Nmasks = Nmasks + 1
	call mask_initialize(masks(Nmasks))
	ats =  maskmanip_make(masks(Nmasks))

	!discard if no atoms in mask
	if(ats == 0) then
		call mask_finalize(masks(Nmasks))
		Nmasks = Nmasks - 1
		xscore_add = 0
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

	allocate(cofactor_def(1:16))			! initially allow up to 16 cofactors
	num_cofactor = 0
	do
		call getlin(chSmallBuf,'Cofactor (. or EOL terminates): ')
		if(chSmallBuf.eq.'.'.or.chSmallBuf.eq.'') exit
		num_cofactor = num_cofactor +1
		cofactor_def(num_cofactor) = chSmallBuf
	end do

	call getlin(chTranslationKey, 'FF translation key: ')
	call getlin(chInput,'Scroring parameters: ')

	xscore_add = Nmasks
	write(desc, 2) nat_pro
2	format('X-Score calculation using ',i6,' atoms')
end function xscore_add

subroutine xsort_bonds		!identify pairs of atoms involved in q-bonds and/or hydrogen bonding !in start 
	integer						::	atom1,atom2
	integer						::	i,j,k
	logical						::	polar

	!allokera listor
	allocate (q_atoms(nqat))
	allocate (q_bonds(nqat*5))		! worst case?
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
		
		!store all bonds between q-atoms, pointers relating arrays q_atoms and q_bonds 
		if ((iqatom(atom1) /= 0) .and. (iqatom(atom2) /= 0)) then
			j = iqatom(atom1)
			k = iqatom(atom2)  

			nqbonds = nqbonds+1
			q_bonds(nqbonds)%a => q_atoms(j)		
			q_bonds(nqbonds)%b => q_atoms(k)	! pointer to atom in q_atoms
			q_bonds(nqbonds)%cod = bnd(i)%cod	! some sort of bond type

			q_atoms(j)%n = q_atoms(j)%n +1		! number of bonds to atom j
			q_atoms(j)%bd(q_atoms(j)%n)%qb => q_bonds(nqbonds)	!pointer to bond in q_bonds
				
			q_atoms(k)%n = q_atoms(k)%n +1
			q_atoms(k)%bd(q_atoms(k)%n)%qb => q_bonds(nqbonds)
		end if
			
	end do
end subroutine xsort_bonds

subroutine XReadAtomTypeConversions(filename, chType)
	character(len=256)	:: filename
	character(len=20)		:: chType
	integer							:: fp, p,p2, num_heading = 0, column, count = 0,filestat,mark,i,j
	character(len=256)	:: buf,buf2
	character(len=12),allocatable		:: tmp(:)		

	fp = freefile()
	open(unit=fp,file=filename,err=999,action='READ')
		do
			! First, find the right section
			read(fp,'(a256)',iostat=filestat) buf
			if(filestat.eq.0) then
				call locase(buf)
				if(trim(buf).eq.'') then
					goto 11					!	skip blank line
				elseif(trim(buf(1:1)).eq.'#') then
					goto 11		! skip comment
				elseif(trim(buf).eq.'{ligand atom type translations}') then
					exit
				end if
			else
				write(*,'(a,a)') 'Fatal: Unable to read ligand atom type translations from ', trim(filename)
				stop
			end if
11	end do

		do
			read(fp,'(a256)',iostat=filestat) buf
			if(filestat.eq.0) then
				call locase(buf)
				if(trim(buf).eq.'') then
					goto 1					!	skip blank line
				elseif(trim(buf(1:1)).eq.'#') then
					goto 1		! skip comment
				elseif(buf(1:5).eq.'sybyl') then
					p = index(buf,'#')
					if(p.ne.0) buf = buf(1:p-1)							! remove traling comment

					! First find out how many headings (including xscore) are present
					p = 7
					num_heading = 1
					do
						read(buf(p:256),*,iostat=filestat) buf2
						if(filestat.eq.0) then
							num_heading = num_heading +1
							p2 = index(buf(p:255),buf2(1:len(trim(buf2))))
							if(buf2(1:len(trim(buf2))).eq.chType) then
								column = num_heading
							end if
							p = p + p2 +len(trim(buf2))
						else
							exit
						end if
					end do
				elseif(buf(1:1).eq.'{') then
					exit
				else
					count = count +1
				end if
			else
				exit
			end if
1		end do

		rewind(fp) 

		if(column.eq.0) then
			write(*,'(a,a,a,a)') 'Fatal error: Column ', adjustl(trim(chType)), ' not found in ', adjustl(trim(filename))
			stop
		end if

		num_atomtype = count
		allocate(atomtypetable(1:num_atomtype))
		allocate(tmp(1:num_heading))

		count = 0
		do
			read(fp,'(a256)',iostat=filestat) buf
			if(filestat.ne.0) goto 10
			if(trim(buf).eq.'') goto 2		! skip blank line
			if(buf(1:1).eq.'#') goto 2		! skip REM
		
			call locase(buf)
			if(trim(buf).eq.'{ligand atom type translations}') then
				do
					read(fp,'(a256)',iostat=filestat) buf
					if(filestat.ne.0) goto 10

					buf = adjustl(trim(buf))
					if(buf.eq.'') goto 3		! skip blank line
					if(buf(1:1).eq.'#') goto 3		! skip REM

					if(buf(1:5).eq.'SYBYL') then
						mark = 1
						goto 3
					end if

					p = index(buf,'#')
					buf = adjustl(buf(1:p-1))			! remove traling comment

					if(mark.ne.0) then									! read translation table
						tmp(1) = ''
						tmp(column) = ''
						read(buf,*,iostat=filestat) tmp(1:column)
						if(filestat.ne.0) goto 10

						if(buf(1:1).eq.'{') then
						 goto 10
						elseif(tmp(column)(1:12).ne.'') then	
							count = count +1
							atomtypetable(count)%xscore = tmp(1)
							atomtypetable(count)%other  = tmp(column)
						end if
					end if
3				end do
			end if
2		end do
10 close(fp)

	do i = 1, size(atomtypetable)
		if(atomtypetable(i)%other.eq.'.') cycle
		do j = i,size(atomtypetable)
			if(atomtypetable(j)%other.eq.'.') cycle
			if(trim(atomtypetable(i)%other).eq.trim(atomtypetable(j)%other) &
			.and.trim(atomtypetable(i)%xscore).ne.trim(atomtypetable(j)%xscore)) then
				write(*,'(a,a,a,a,a,a,a,a)') &
				'Fatal: Ambiguous entries in type translation matrix (file ', &
				trim(adjustl(filename)), '). ', &
				trim(adjustl(atomtypetable(j)%other)), &
				' maps to both ', &
				adjustl(trim(atomtypetable(i)%xscore)), &
				' and ', &
				adjustl(trim(atomtypetable(j)%xscore))
				stop
			end if
		end do
	end do


	write(*,'(a,i4,a,a,a)') 'Read ', count, ' atom translations for force field ', adjustl(trim(chType)), '.'
	return

999 write(*,'(a,a,a)') '>>>>> ERROR: Unable to read ', adjustl(trim(filename)), '. Exiting.'
		stop
end subroutine XReadAtomTypeConversions

subroutine XReadResidueConversions(filename, chType)
	character(len=256)	:: filename
	character(len=20)		:: chType
	integer							:: fp, p, filestat
	character(len=256)	:: buf
	type(RESIDUE_NAME_CONVERSION),pointer		:: tmp(:)		

	fp = freefile()
	open(unit=fp,file=filename,err=999,action='READ')
		do
			! First, find the right section
			read(fp,'(a256)',iostat=filestat) buf
			if(filestat.eq.0) then
				call locase(buf)
				if(trim(buf).eq.'') then
					goto 11					!	skip blank line
				elseif(trim(buf(1:1)).eq.'#') then
					goto 11		! skip comment
				elseif(trim(buf).eq.'{residue type translations}') then
					exit
				end if
			else
				write(*,'(a,a)') 'Fatal: Unable to read residue translations from ', trim(filename)
				stop
			end if
11	end do

		num_residue = 0
		do
			read(fp,'(a256)',iostat=filestat) buf
			if(filestat.eq.0) then
				if(trim(buf).eq.'') then						!	skip blank line
					goto 12					
				elseif(trim(buf(1:1)).eq.'#') then	! skip comment
					goto 12					
				elseif(buf(1:1).eq.'{') then
					exit
				else
					p = index(buf,'#')
					if(p.ne.0) buf = buf(1:p-1)				! remove traling comment

					num_residue = num_residue +1
					! Make sure storage has been allocated
					if(.not. associated(residuetable)) allocate(residuetable(9))
					
					! Expand table if it's too small
					if(num_residue>size(residuetable)) then			 
						allocate(tmp(size(residuetable)+4))				
						tmp(1:size(residuetable)) = residuetable(1:size(residuetable))
						deallocate(residuetable)
						residuetable => tmp
					end if

					! Read table
					read(buf,*) residuetable(num_residue)%xscore,residuetable(num_residue)%other
				end if
			else
				exit
			end if
12	end do

	close(fp)
	return

999 write(*,'(a,a,a)') '>>>>> ERROR: Unable to read ', adjustl(trim(filename)), '. Exiting.'
		stop
end subroutine XReadResidueConversions

subroutine XPAtomNameConversions(filename, chType)
	character(len=256)	:: filename
	character(len=20)		:: chType,name1,name2
	integer							:: fp, p, i,  filestat, mark
	character(len=256)	:: buf
	character(len=10)		:: head
	character(len=10)		:: section(1:5)
	integer							:: stackpointer,current_res
	!type(RESIDUE_NAME_CONVERSION),pointer		:: tmp(:)		

	call locase(chType)

	fp = freefile()
	open(unit=fp,file=filename,err=999,action='READ')
		do
			! First, find the right section
			read(fp,'(a256)',iostat=filestat) buf
			if(filestat.eq.0) then
				call locase(buf)
				if(trim(buf).eq.'') then						!	skip blank line
					goto 1					
				elseif(trim(buf(1:1)).eq.'#') then	! skip comment
					goto 1		
				elseif(trim(buf).eq.'{protein atom name translations}') then

					! Read and make corrections
					mark = 0
					stackpointer = 0
					do
						read(fp,'(a256)',iostat=filestat) buf
						if(filestat.eq.0) then
							if(trim(buf).eq.'') then						!	skip blank line
								goto 2					
							elseif(trim(buf(1:1)).eq.'#') then	! skip comment
								goto 2					
							elseif(buf(1:1).eq.'{'.and.mark.eq.1) then
								exit
							else
								p = index(buf,'#')
								if(p.ne.0) buf = buf(1:p-1)				! remove traling comment

								! Prepare buffer
								call locase(buf)
								read(buf,*), head

								if(head(1:1).eq.'<') then
									p = index(head,'>')
									head = head(2:p-1)							! extract head
									
									if(head(1:1).eq.'.') then				! end of section
										head = head(2:len(head))
										if(head.ne.section(stackpointer)) then
											write(*,*) 'Syntax error in ', filename
											stop
										else
											section(stackpointer) = ''
											stackpointer = stackpointer -1

											if(stackpointer.eq.0.and.mark.ne.0) goto 3
										end if
									else														! beginning of section
										stackpointer = stackpointer +1
										section(stackpointer) = head

										if(stackpointer.eq.2) then		! residue name
											do i = 1,num_patomrestype
												call upcase(head)
												if(patomnametable(i)%res.eq.head) then
													current_res = i
													exit
												end if
											end do
										elseif(stackpointer.eq.1) then	! forcefield name
											if(head.eq.chType) mark = 1
										end if
									end if
								else
									if((mark.ne.0).and.(current_res.ne.0)) then
										read(buf,*) name1, name2
										call upcase(name1)
										call upcase(name2)
										do i = 1,patomnametable(current_res)%natom
											if(name1.eq.patomnametable(current_res)%xname(i)) then
												patomnametable(current_res)%xname(i) = name2
												exit
											end if
										end do
									end if
								end if
							end if
						else
							exit
						end if
2					end do				
					exit
				end if
			else
				write(*,'(a,a)') 'Fatal: Unable to read residue translations from ', trim(filename)
				stop
			end if
1		end do

3	close(fp)
	return

999 write(*,'(a,a,a)') '>>>>> ERROR: Unable to read ', adjustl(trim(filename)), '. Exiting.'
		stop
end subroutine XPAtomNameConversions

character(len=12) function XTranslateAtomType(othertype)
	character	:: othertype, buf  !removed (len=12) to compile on linux
	integer						:: i

	read(othertype, *, iostat=i) buf
	
	do i = 1,num_atomtype
		if(atomtypetable(i)%other.eq.buf) then
			XTranslateAtomType = atomtypetable(i)%xscore
			return
		end if
	end do

	XTranslateAtomType = 'Un'
end function XTranslateAtomType

subroutine XReadLib(filename)
	character(len=256)	:: filename
	integer							:: i,fp, filestat,p,count,b_atomsection
	character(len=6)		:: slask
	character(len=256)	:: buf
	integer							:: atom_count(1:200)						! worst case
	type(PATOM_NAME_CONVERSION), pointer :: tmp_table(:)

	fp = freefile()
	count = 0
	open(unit=fp,file=filename,err=999,action='READ')
		do
			read(fp,'(a256)',iostat=filestat) buf
			if(filestat.eq.0) then
				if(trim(buf).eq.'') then
					goto 1															!	skip blank line
				elseif(trim(buf(1:1)).eq.'#') then
					goto 1															! skip comment
				elseif(trim(buf(1:1)).eq.'*') then
					goto 1															! skip comment
				elseif(trim(buf(1:1)).eq.'!') then
					goto 1															! skip comment
				else
					p = index(buf,'!')
					if(p.ne.0) buf = buf(1:p-1)										! remove traling comment

					if(buf(1:1).eq.'{') then
						count = count +1

						! now count the number of atoms in this residue
						b_atomsection = 0
						do
							read(fp,'(a256)',iostat=filestat) buf
							if(filestat.ne.0) then
								exit
							elseif(trim(buf).eq.'') then
								goto 2															!	skip blank line
							elseif(trim(buf(1:1)).eq.'#') then
								goto 2															! skip comment
							elseif(trim(buf(1:1)).eq.'*') then
								goto 2															! skip comment
							elseif(trim(buf(1:1)).eq.'!') then
								goto 2															! skip comment
							elseif(index(buf,'[atoms]').ne.0) then
								b_atomsection =	1										! set flag 
							elseif(index(buf,'[').ne.0.and.b_atomsection.ne.0) then			! reach end of [atoms]-section
								b_atomsection = 0
								exit
							else
								if(b_atomsection.ne.0) atom_count(count) = atom_count(count) +1	! if no rules triggered, this is an atom line
							end if
	2					end do
					end if
				end if
			else
				exit
			end if
1		end do

		rewind(fp) 

		if(.not. associated(patomnametable)) allocate(patomnametable(1:count))
		
		if(size(patomnametable)<num_patomrestype+count) then									! increase table size
			allocate(tmp_table(1:num_patomrestype+count))
			tmp_table(1:num_patomrestype) = patomnametable(1:num_patomrestype)
			deallocate(patomnametable)
			patomnametable => tmp_table
		end if

		do i = 1,count
			patomnametable(num_patomrestype + i)%natom = atom_count(i)
			allocate(patomnametable(num_patomrestype + i)%qname(1:atom_count(i)))	! allocate name tables
			allocate(patomnametable(num_patomrestype + i)%xname(1:atom_count(i)))	!  - " -
		end do

		! now read the residue and atom name information
		count = 0
		do
			read(fp,'(a256)',iostat=filestat) buf
			if(filestat.eq.0) then
				if(trim(buf).eq.'') then
					goto 3															!	skip blank line
				elseif(trim(buf(1:1)).eq.'#') then
					goto 3															! skip comment
				elseif(trim(buf(1:1)).eq.'*') then
					goto 3															! skip comment
				elseif(trim(buf(1:1)).eq.'!') then
					goto 3															! skip comment
				else
					p = index(buf,'!')
					if(p.ne.0) buf = buf(1:p-1)										! remove traling comment

					if(buf(1:1).eq.'{') then						! begin residue information
						count = count +1
						p = index(buf,'}')
						patomnametable(num_patomrestype + count)%res = buf(2:p-1)	! read residue name

						! now read atom name table
						i = 0
						b_atomsection = 0
						do
							read(fp,'(a256)',iostat=filestat) buf
							if(filestat.ne.0) then
								exit
							elseif(trim(buf).eq.'') then
								goto 4															!	skip blank line
							elseif(trim(buf(1:1)).eq.'#') then
								goto 4															! skip comment
							elseif(trim(buf(1:1)).eq.'*') then
								goto 4															! skip comment
							elseif(trim(buf(1:1)).eq.'!') then
								goto 4															! skip comment
							elseif(index(buf,'[atoms]').ne.0) then
								b_atomsection =	1										! set flag
								i = 0										
							elseif(index(buf,'[').ne.0.and.b_atomsection.ne.0) then			! reached end of [atoms]-section
								b_atomsection = 0
								exit
							else
								if(b_atomsection.ne.0) then
									i = i +1														! read name information
									read(buf(1:256),*,iostat=filestat) slask, &
																										 patomnametable(num_patomrestype+count)%xname(i),&
																										 patomnametable(num_patomrestype+count)%qname(i)
								end if
							end if
4						end do
					end if
				end if
			else
				exit
			end if
3		end do
	close(fp)

	num_patomrestype = num_patomrestype + count
	return

999 write(*,'(a,a,a)') '>>>>> ERROR: Unable to read ', adjustl(trim(filename)), '. Exiting.'
		stop
end subroutine XReadLib

character(len=10) function XTranslatePAtomName(atom_res_id, res)
	! returns the correct aton name as found in residue table
	! dependant on the correct ordering of atoms within residue!
	integer						:: atom_res_id
	character(len=10)	:: res
	integer						:: i


	do i = 1,num_patomrestype
		if(patomnametable(i)%res.eq.res) then
			XTranslatePatomName = patomnametable(i)%xname(atom_res_id)
			return
		end if
	end do

	XTranslatePAtomName = 'Un'
end function XTranslatePAtomName

character(len=10) function XTranslateResidue(res)
	! Translates some residue names used by Q to XScore equivalents
	character(len=8)	:: res
	integer						:: i

	do i = 1,num_residue
		if(trim(res).eq.trim(residuetable(i)%other)) then
			XTranslateResidue = residuetable(i)%xscore
			return
		end if
	end do

	XTranslateResidue = res
end function XTranslateResidue

subroutine XReadInput()
	integer							:: i,fp,filestat,p
	character(len=256)	:: buf, head, chVal
	character(len=256)			:: line, restofline
	integer							:: res

	chInput = adjustl(trim(chInput))

	if(chInput.ne.'default') then
		! Open file
		fp = freefile()
		open(unit=fp,file=chInput,err=999,action='READ')
			do
				read(fp,'(a256)',iostat=filestat) buf

				if(filestat.eq.0) then
					if(trim(buf).eq.'') then
						goto 1															!	skip blank line
					elseif(trim(buf(1:1)).eq.'#') then
						goto 1															! skip comment
					else
						p = index(buf,'#')
						if(p.ne.0) buf = buf(1:p-1)										! remove traling comment
						
						read(buf, *,iostat=filestat) head
						call upcase(head)
						head = adjustl(trim(head))

						read(buf,*,iostat=filestat) head, chVal
						
						select case(head)			
							case('SHOW_ATOM_BIND_SCORE')
								read(buf,*,iostat=filestat) head, input%show_abs
							case('SHOW_TOTAL_SCORE')
								read(buf,*,iostat=filestat) head, input%show_total

							case('APPLY_HPSCORE')
								read(buf,*,iostat=filestat) head, input%apply_hpscore
							case('HPSCORE_CVDW')
								read(buf,*,iostat=filestat) head, input%hpscore_cvdw
							case('HPSCORE_CHB')
								read(buf,*,iostat=filestat) head, input%hpscore_chb
							case('HPSCORE_CHP')
								read(buf,*,iostat=filestat) head, input%hpscore_chp
							case('HPSCORE_CRT')
								read(buf,*,iostat=filestat) head, input%hpscore_crt
							case('HPSCORE_C0')
								read(buf,*,iostat=filestat) head, input%hpscore_c0

							case('APPLY_HMSCORE')
								read(buf,*,iostat=filestat) head, input%apply_hmscore
							case('HMSCORE_CVDW')
								read(buf,*,iostat=filestat) head, input%hmscore_cvdw
							case('HMSCORE_CHB')
								read(buf,*,iostat=filestat) head, input%hmscore_chb
							case('HMSCORE_CHM')
								read(buf,*,iostat=filestat) head, input%hmscore_chm
							case('HMSCORE_CRT')
								read(buf,*,iostat=filestat) head, input%hmscore_crt
							case('HMSCORE_C0')
								read(buf,*,iostat=filestat) head, input%hmscore_c0

							case('APPLY_HSSCORE')
								read(buf,*,iostat=filestat) head, input%apply_hsscore
							case('HSSCORE_CVDW')
								read(buf,*,iostat=filestat) head, input%hsscore_cvdw
							case('HSSCORE_CHB')
								read(buf,*,iostat=filestat) head, input%hsscore_chb
							case('HSSCORE_CHS')
								read(buf,*,iostat=filestat) head, input%hsscore_chs
							case('HSSCORE_CRT')
								read(buf,*,iostat=filestat) head, input%hsscore_crt
							case('HSSCORE_C0')
								read(buf,*,iostat=filestat) head, input%hsscore_c0

							case('SHOW_LIGAND')
								read(buf,*,iostat=filestat) head, input%show_ligand
							case('SHOW_PROTEIN')
								read(buf,*,iostat=filestat) head, input%show_protein
							case('SHOW_COFACTOR')
								read(buf,*,iostat=filestat) head, input%show_cofactor
							case('SHOW_BONDS')
								read(buf,*,iostat=filestat) head, input%show_bonds

							case('RESIDUE_DEFINITIONS')
								read(buf,'(a256)') line
								read(line,*) head
								!workaround to allow / in strings
								input%residue_def = adjustl(line(len_trim(head)+1:len(line)))
!								read(buf,*,iostat=filestat) head, input%residue_def
							case('ATOM_DEFINITIONS')
								read(buf,'(a256)') line
								read(line,*) head
								!workaround to allow / in strings
								input%atom_def = adjustl(line(len_trim(head)+1:len(line)))
!								read(buf,*,iostat=filestat) head, input%atom_def
							case('LOGP_DEFINITIONS')
								read(buf,'(a256)') line
								read(line,*) head
								!workaround to allow / in strings
								input%logp_def = adjustl(line(len_trim(head)+1:len(line)))
!								read(buf,*,iostat=filestat) head, input%logp_def
							case('SURFACE_DEFINITIONS')
								read(buf,'(a256)') line
								read(line,*) head
								!workaround to allow / in strings
								input%surface_def = adjustl(line(len_trim(head)+1:len(line)))
!								read(buf,*,iostat=filestat) head, input%surface_def
							case('ATOM_TRANSLATIONS')
								read(buf,'(a256)') line
								read(line,*) head
								!workaround to allow / in strings
								input%atom_translations = adjustl(line(len_trim(head)+1:len(line)))
!								read(buf,*,iostat=filestat) head, input%atom_translations

						end select	
					end if
				else
					exit
				end if
1			end do
		close(fp)
		
		input%num_method = 0
		if(input%apply_hpscore.eq.'YES') input%num_method = input%num_method +1
		if(input%apply_hmscore.eq.'YES') input%num_method = input%num_method +1
		if(input%apply_hsscore.eq.'YES') input%num_method = input%num_method +1

		input%qatom = -1
	else
		input%show_abs			= 'NO'
		input%show_total		= 'YES'
		input%show_ligand		= 'YES'
		input%show_protein	= 'NO'
		input%show_cofactor	= 'YES'
		input%show_cofactor	= 'NO'

		input%num_method		=  3
		input%apply_hpscore = 'YES'
		input%apply_hmscore = 'YES'
		input%apply_hsscore = 'YES'

		input%hpscore_cvdw  =  0.004 
		input%hpscore_chb   =  0.054
		input%hpscore_chp   =  0.009
		input%hpscore_crt		= -0.061
		input%hpscore_c0    =  3.441

		input%hmscore_cvdw  =  0.004
		input%hmscore_chb   =  0.101
		input%hmscore_chm   =  0.387
		input%hmscore_crt		= -0.097
		input%hmscore_c0    =  3.567

		input%hsscore_cvdw  =  0.004
		input%hsscore_chb   =  0.073
		input%hsscore_chs   =  0.004
		input%hsscore_crt		= -0.090
		input%hsscore_c0    =  3.328

		input%residue_def					= 'residue_def_xtool.dat'
		input%atom_def						= 'atom_def_xtool.dat'
		input%logp_def						= 'atom_def_xlogp.dat'
		input%surface_def					= 'surface_def_xtool.dat'
		input%atom_translations		= 'atom_type_translations.dat'

	end if

	return

999 write(*,'(a,a,a)') '>>>>> ERROR: Unable to read ', adjustl(trim(chInput)), '. Exiting.'
		stop
end subroutine XReadInput

subroutine XReadInputQatom
	integer							:: i,fp,filestat,p
	character(len=256)	:: buf, head, chVal
	integer							:: res

	if(input%show_abs(1:6).eq.'QATOM=') then
		chVal = input%show_abs(7:20)

		p = index(chVal,':')
		if(p.ne.0) then

			head = adjustl(trim(chVal(1:p-1)))
			buf  = adjustl(trim(chVal(p+1:256)))

			read(buf,*,iostat=filestat) res

			do i = 0,ligand%mol%num_atom -1
				if(adjustl(trim(ligand%mol%atom(i)%res_id)).eq.adjustl(trim(head))) then
					input%qatom = i + res
					exit
				end if
			end do 
		else
			read(chVal,*,iostat=filestat) input%qatom
		end if
	
		input%show_abs = 'QATOM'
		if(input%qatom<0.or.input%qatom>ligand%mol%num_atom -1) then
			write(*,'(a,a,a)') '>>> WARNING: Q-atom ', trim(chVal), ' does not exist. ABS will not be displayed.'
			warn = warn +1
			input%show_abs = 'NO'
		end if

	end if
end subroutine XReadInputQatom

subroutine xscore_precalc
	integer							:: i,j, p
	character(len=256)	:: chBuf
	character(len=256)	:: filename
	real(kind=prec)								:: score

	! ======================================
	! 0. READ PARAMETERS AND DO PREPARATIONS
	! ======================================


			call XReadInput()
			
			filename = input%residue_def								!'RESIDUE_DEF_XTOOL'
			call ForceField_Read_RESIDUE_DEF(filename)

			filename = input%atom_def										!'ATOM_DEF_XTOOL'
			call ForceField_Read_ATOM_DEF(filename)

			filename = input%logp_def										!'ATOM_DEF_XLOGP'
			call ForceField_Read_XATOM_DEF(filename)

			filename = input%surface_def								!'SURFACE_DEF_XTOOL'
			call ForceField_Read_Surface_Def(filename)

			filename = input%atom_translations					!'XSCORE_ATOM_TYPE_CONVERSIONS.dat'		! atom type conversions are only used on ligand atoms
			call XReadAtomTypeConversions(filename,chTranslationKey)		! read atom type conversions
			call XReadResidueConversions(filename,chTranslationKey)			! read residue conversions

			! Read atom names from residue tables of lib-files
			! This is done prior to residue renameing and atom renameing
			chBuf = lib_files
			do
				p = index(chBuf,';')
				if(p>0) then
					filename = chBuf(1:p-1)
					call XReadLib(filename)
					chBuf(1:256) = chBuf(p+1:256)
				else	
					call XReadLib(chBuf)
					exit
				end if
			end do					
			write(*,'(a,i4,a)') 'Read ', num_patomrestype, ' residues from library.'

			! Make corrections to what was read from lib-files
			filename = input%atom_translations					!'XSCORE_ATOM_TYPE_CONVERSIONS.dat'
			call XPAtomNameConversions(filename,chTranslationKey)

	! ============================
	! 1. READ AND TRANSLATE LIGAND
	! ============================

			if(.not. qatom_load_atoms(fep_file)) then
				stop 'Failed to read Q-atom list from FEP file.'
			end if

			allocate(iqatom(nat_solute))	! Allocate mem for every solute atom
			iqatom(:) = 0
			do i=1,nqat
				iqatom(iqseq(i)) = i		! Set references so that iqatom = 0 if protein and {number} of ligand
			end do

			! Extract bonds between q-atoms
			call xsort_bonds						

			! Set topology numbers for q-atoms
			q_atoms(:)%top_nr = iqseq(:)

			! Translate q-atom data to ligand structure
			call xligand_translate(ligand,nqat,nqbonds,q_atoms,q_bonds)
			ligand%mol%name = adjustl(trim(title))
			do i = 0,ligand%mol%num_atom -1
				ligand%mol%atom(i)%type2 = ligand%mol%atom(i)%ttype
				ligand%mol%atom(i)%ttype = XTranslateAtomType(ligand%mol%atom(i)%ttype)
			end do
			write(*,'(a,i4,a,i4,a)') 'Translated ', ligand%mol%num_atom, ' atoms and ', ligand%mol%num_bond, ' bonds to ligand'

			call Ligand_Value_Atom(ligand)

			! Find monitored q-atom if specified. This has to be done after ligand atom names are assigned.
			call XReadInputQatom()

			j = 0
			do i = 0,ligand%mol%num_atom -1
				if(ligand%mol%atom(i)%valid.eq.0) then
!					write(*,'(a,i6,a,a,a,a,a)') 'WARNING: Atom ', mol%atom(i)%id, ' (', trim(adjustl(mol%atom(i)%qtype)),') in residue ', trim(adjustl(mol%atom(i)%residue)), ' does not match any type rule and will be ignored.'
					j = j +1
				end if
			end do
			if(j>0) then
				write(*,'(a,i5,a,i5,a,a,a)') &
				'WARNING: ', j, ' atoms (out of ', &
				ligand%mol%num_atom, ') in ', &
				adjustl(trim(ligand%mol%name)), &
				' (ligand) will be ignored. Note that nonpolar hydrogens are automatically ignored.'
				warn = warn +1
			end if
			if(adjustl(trim(input%show_ligand)).eq.'YES') call Molecule_Show_Contents(ligand%mol)

	! =============================
	! 2. READ AND TRANSLATE PROTEIN
	! =============================

			! First count protein bonds
			j = 0
			do i = 1,nbonds_solute
!				if((bnd(i)%i<=ligand_offset).and.(bnd(i)%j<=ligand_offset)) j = j +1
				if((iqatom(bnd(i)%i).eq.0).and.(iqatom(bnd(i)%j).eq.0)) j = j +1
			end do

			allocate(top2prot(1:nat_solute))						! top -> protein index translation matrix

!			call xprotein_translate(protein, offset, xtop)
			call xprotein_translate(protein, nat_solute, calc_xtop,nbonds_solute,bnd)
			protein%name = adjustl(trim(title))

			! Need to value cofactor-protein atoms using Molecule_Value_Atom instead
			! of Protein_Value_Atom. Cofactor atoms are invalidated and added separately later.
			! Need to flag cofactors before residue names are translated.
			call xflag_cofactors

			do i = 0,protein%num_atom -1
				! Rename atoms according to lib-files
				protein%atom(i)%name		= XTranslatePAtomName(protein%atom(i)%atom_res_id, protein%atom(i)%residue)
				
				! Rename residues according to translation table
				protein%atom(i)%old_residue = protein%atom(i)%residue
				protein%atom(i)%residue = XTranslateResidue(protein%atom(i)%residue)
			end do
			
			write(*,'(a,i5,a,i5,a)') 'Translated ', protein%num_atom, ' atoms and ', protein%num_bond, ' bonds to protein'

			call Protein_Value_Atom(protein,flag=1)

			allocate(cofactor(1:num_cofactor))
			do i = 1, num_cofactor
!				call xextract_cofactor(cofactor(i),cofactor_def(i),i,offset, nbonds,xtop,bnd)
				call xextract_cofactor(cofactor(i),cofactor_def(i),i,nat_solute, nbonds_solute,calc_xtop,bnd)

				if(cofactor(i)%mol%num_atom.eq.0) then
					write(*,'(a,a,a)') '>>>>ERROR: Cofactor definition ', adjustl(trim(cofactor_def(i))), ' is empty.'
					stop
				end if

				do j = 0,cofactor(i)%mol%num_atom -1
					cofactor(i)%mol%atom(j)%type2 = cofactor(i)%mol%atom(j)%ttype
					cofactor(i)%mol%atom(j)%ttype = XTranslateAtomType(cofactor(i)%mol%atom(j)%ttype)
				end do
				write(*,'(a,i4,a,i4,a,a)') 'Translated ', &
				cofactor(i)%mol%num_atom, ' atoms and ', &
				cofactor(i)%mol%num_bond, ' bonds to cofactor ', cofactor_def(i)

				call Ligand_Value_Atom(cofactor(i))
				if(adjustl(trim(input%show_cofactor)).eq.'YES') call Molecule_Show_Contents(cofactor(i)%mol)
				
				call Protein_Merge_Cofactor(protein,cofactor(i))
			end do
			
			! Update neighbour count for protein atoms (some waters may mess this up)
			do i = 0,protein%num_atom -1
				do j = 0,6
					if(protein%atom(i)%neib(j).eq.0) then
						protein%atom(i)%num_neib = j
						exit
					end if
				end do
			end do

			call Protein_Define_Pocket(protein, ff,ligand, 10.0_prec)

			j = 0
			do i = 0,protein%num_atom -1
				if(protein%atom(i)%valid.eq.0) then
					if(protein%atom(i)%xtype.ne.'H' .and. &
					protein%atom(i)%xtype.ne.'Un') then
					    write(*,'(a,i6,a,a,a,a,a)') &
					    'WARNING: Atom ', &
					    protein%atom(i)%id, ' (', &
					    trim(adjustl(protein%atom(i)%xtype)), &
					    ') in residue ', &
					    trim(adjustl(protein%atom(i)%residue)), &
					    ' will be ignored.'
					endif
					j = j +1
				end if
			end do
			if(j>0) then
				write(*,'(a,i4,a,i4,a,a,a)') &
				'WARNING: ', j, ' atoms (out of ', &
				protein%num_atom, ') in ', &
				adjustl(trim(protein%name)), &
				' (protein) will be ignored. Note that nonpolar hydrogens are automatically ignored.'
				warn = warn +1
			end if
			if(adjustl(trim(input%show_protein)).eq.'YES') call Protein_Show_Contents(protein)


	! ==============================
	! 3. SCORE INITIAL CONFIGURATION
	! ==============================

		if(bDoTopcalc.eq.1) then 
			write(*,'(a)') 'Scoring initial configuration'
			score = Ligand_Calculate_Binding_Score(ligand, input,protein)
			call Ligand_Scoring_Stats(ligand)
			write(*,'(a, t28,f6.1,  t36,f5.2, t43,f5.1, t51,f5.2, t57,f5.1, t65,f4.1, t73,f6.3)') &
			'Total score of topology', ligand%vdw,ligand%hb, &
			ligand%hp,ligand%hm,ligand%hs,ligand%rt,ligand%bind_score
		end if

end subroutine xscore_precalc

subroutine xflag_cofactors
	integer							:: i,j,count
	character(len=10)		:: restype

	! valid cofactor definitions: restype=HEM			-> all atoms in residues named HEM

	protein%atom(:)%iscofactor = 0
	do i = 1, num_cofactor
		if(cofactor_def(i)(1:7).eq.'restype'.or.cofactor_def(i)(1:7).eq.'RESTYPE') then
			restype = adjustl(trim(cofactor_def(i)(9:80)))
			do j = 0,protein%num_atom -1
				if(protein%atom(j)%residue.eq.restype) protein%atom(j)%iscofactor = i
			end do
		end if
	end do
end subroutine xflag_cofactors

subroutine xextract_cofactor(cf,cf_def,conumber,nAtoms,nBonds,coordinates,bonds)
	type(tXLigand)			:: cf
	character(len=80)		:: cf_def
	integer							:: conumber
	integer							:: nAtoms,nBonds
	real(kind=prec)							:: coordinates(:)
	type(BOND_TYPE)			:: bonds(:)

	integer							:: i,j,count,mark,floating_offset
	character(len=10)		:: restype
	integer							:: iRes					! current residue
	integer							:: iNextRes			! starting point of next residue

	integer,pointer			:: mask(:),mask_res(:)
	type(tXBond)				:: tmp
	integer,pointer			:: prt2cof(:)

	! valid cofactor definitions: restype=HEM			-> all atoms in residues named HEM

	! Create atom mask for this cofactor
	if(cf_def(1:7).eq.'restype'.or.cf_def(1:7).eq.'RESTYPE') then
		restype = trim(adjustl(cf_def(9:80)))
		count = 0

		iRes = 0
		iNextRes = 1		 
		do i = 1,nAtoms
			if(i.eq.iNextRes) then
				iRes = iRes +1
				if(iRes<size(res)) iNextRes = res(iRes +1)%start
			end if
			if(res(iRes)%name.eq.restype) then
				call Int_Pushback(mask,i)					! push back this atom
				call Int_Pushback(mask_res,iRes)	! push back corresponding residue id
				glb_cofactor(i) = conumber				! set global cofactor flag
			end if
		end do
		
		count = 0
		do i = 0,protein%num_atom -1
			if(protein%atom(i)%residue.eq.restype) count = count +1
		end do

		cf%mol%num_atom = count
		allocate(cf%mol%atom(0:cf%mol%num_atom-1))
		
		allocate(prt2cof(0:protein%num_atom-1))

		! Copy protein atom data to cofactor
		cf%mol%name = cf_def
		cf%mol%num_ring = 0
		nullify(cf%mol%ring)
		j = -1
		floating_offset = 0
		do i = 0,protein%num_atom -1
			if(protein%atom(i)%residue.eq.restype) then
				j = j +1
				cf%mol%atom(j)%iscofactor		= conumber				! set cofactor number
				cf%mol%atom(j)%coor(0)			= protein%atom(i)%coor(0)
				cf%mol%atom(j)%coor(1)			= protein%atom(i)%coor(1)
				cf%mol%atom(j)%coor(2)			= protein%atom(i)%coor(2)
				cf%mol%atom(j)%id						= j +1						! need to set this temporarily
				cf%mol%atom(j)%atom_res_id	= j +1
				cf%mol%atom(j)%name					= ''							! not used. names are assigned using lib file and index within residue
				cf%mol%atom(j)%residue			= 'COF'						! res(iRes)%name
				cf%mol%atom(j)%res_id				= protein%atom(i)%res_id
				cf%mol%atom(j)%ttype				= protein%atom(i)%ttype			! ttype and xtype are set by Value_Atom
				cf%mol%atom(j)%xtype				= 'Un'						! cf%mol%atom(i -1)%ttype
				cf%mol%atom(j)%type2				= 'Un' 
				cf%mol%atom(j)%hb						= 'N'
				cf%mol%atom(j)%valid				= 1								! default
				cf%mol%atom(j)%origin				= 2								! This is a cf%mol atom
				cf%mol%atom(j)%part					= 1								! regular atom
				cf%mol%atom(j)%occupancy		= 0.0							! default
				cf%mol%atom(j)%solv					= 0
				cf%mol%atom(j)%bfactor			= 0
				cf%mol%atom(j)%logp					= 0
				cf%mol%atom(j)%num_neib     = protein%atom(i)%num_neib
				where(protein%atom(i)%neib.ne.0)
					cf%mol%atom(j)%neib				= protein%atom(i)%neib -floating_offset
				elsewhere
					cf%mol%atom(j)%neib       = 0
				end where
				cf%mol%atom(j)%bond(:)			= protein%atom(i)%bond(:)
				cf%mol%atom(j)%root(:)			= 0
				cf%mol%atom(j)%temp					= i								! save insertion index
				prt2cof(i)									= j								! save protein index -> cof index
			else
				floating_offset = floating_offset +1
			end if
		end do
	
		! Translate bonds
		count = 0
		do i = 0,protein%num_bond -1
!			write(*,*) i, protein%bond(i)%atom_1-1, protein%bond(i)%atom_2-1, size(protein%atom)
			if(protein%atom(protein%bond(i)%atom_1-1)%iscofactor.eq.conumber.and. &
			protein%atom(protein%bond(i)%atom_2-1)%iscofactor.eq.conumber) then
				count = count +1
			end if
		end do

		cf%mol%num_bond = count
		allocate(cf%mol%bond(0:cf%mol%num_bond -1))
		count = 0

		! No need to copy bonds later since they are never removed from protein
		count = 0
		do i = 0,protein%num_bond -1
			if(protein%atom(protein%bond(i)%atom_1-1)%iscofactor.eq.conumber.and. &
			protein%atom(protein%bond(i)%atom_2-1)%iscofactor.eq.conumber) then
				cf%mol%bond(count)%id				= count
				cf%mol%bond(count)%valid		= 1
				cf%mol%bond(count)%atom_1		= prt2cof(protein%bond(i)%atom_1-1) +1
				cf%mol%bond(count)%atom_2		= prt2cof(protein%bond(i)%atom_2-1) +1
				cf%mol%bond(count)%ttype		= protein%bond(i)%ttype	 !SYBYL_bond_type(bonds(i)%cod)
				cf%mol%bond(count)%neib(:)	= 0
				count = count +1
			end if
		end do

		do
			mark = 0
			do i = 0, cf%mol%num_bond -2
				if((cf%mol%bond(i)%atom_1.eq.cf%mol%bond(i+1)%atom_1).and. &
					(cf%mol%bond(i)%atom_2.eq.cf%mol%bond(i+1)%atom_2)) goto 2
				if(cf%mol%bond(i)%atom_1 <  cf%mol%bond(i+1)%atom_1) goto 2
				if((cf%mol%bond(i)%atom_1.eq.cf%mol%bond(i+1)%atom_1) .and. &
					(cf%mol%bond(i)%atom_2<=cf%mol%bond(i+1)%atom_2)) goto 2

				! If no rule triggered, the bonds are in the wrong order
				mark = 1
				call Bond_Swap(cf%mol%bond(i),cf%mol%bond(i+1))
	2		end do
			if(mark.eq.0) exit
		end do

		! Reset bond id
		do i = 0,cf%mol%num_bond -1
			cf%mol%bond(i)%id = i +1
		end do
	end if

	deallocate(prt2cof)
	return

	if(.not. associated(mask).or.size(mask).eq.0) then
		write(*,'(a,a)') 'Fatal: No atoms match cofactor description ', trim(adjustl(cf_def))
		stop
	end if
	
!	! Extract atom data for this cofactor
!	cf%mol%num_atom = size(mask)
!	allocate(cf%mol%atom(0:cf%mol%num_atom-1))
!	cf%cofactor_offset = mask(0) -1
!	cf%mol%name = cf_def
!	cf%mol%num_ring = 0
!	nullify(cf%mol%ring)
!	
!	do i = 0,size(mask) -1
!		j			= mask(i) 
!		iRes	= mask_res(i)	
!
!		cf%mol%atom(i)%coor(0)			= coordinates(3*(j-1) +1)
!		cf%mol%atom(i)%coor(1)			= coordinates(3*(j-1) +2)
!		cf%mol%atom(i)%coor(2)			= coordinates(3*(j-1) +3)
!		cf%mol%atom(i)%id						= i +1
!		cf%mol%atom(i)%atom_res_id	= i +1
!		cf%mol%atom(i)%name					= ''							! not used. names are assigned using lib file and index within residue
!		cf%mol%atom(i)%residue			= 'COF'						! res(iRes)%name
!		write(cf%mol%atom(i)%res_id, '(I8)') iRes			! use internal read to convert integer -> character
!		cf%mol%atom(i)%ttype				= tac(iac(j))			! ttype and xtype are set by Value_Atom
!		cf%mol%atom(i)%xtype				= 'Un'						! cf%mol%atom(i -1)%ttype
!		cf%mol%atom(i)%type2				= 'Un' 
!		cf%mol%atom(i)%hb						= 'N'
!		cf%mol%atom(i)%valid				= 1								! default
!		cf%mol%atom(i)%origin				= 2								! This is a cf%mol atom
!		cf%mol%atom(i)%part					= 1								! regular atom
!		cf%mol%atom(i)%occupancy		= 0.0							! default
!		cf%mol%atom(i)%solv					= 0
!		cf%mol%atom(i)%bfactor			= 0
!		cf%mol%atom(i)%logp					= 0
!		cf%mol%atom(i)%neib(:)			= 0
!		cf%mol%atom(i)%bond(:)			= 0
!		cf%mol%atom(i)%root(:)			= 0
!	end do
!	
!	j = -1
!	do i = 0,protein%num_atom -1
!		if(protein%atom(i)%residue.eq.restype) then	
!			j = j +1
!			cf%mol%atom(j)%temp = i											! save insertion index, this works
!		end if
!	end do

		
	! Extract bond data for this cofactor. Note: covalent bonds to non-cofactor protein atoms are broken (but are restored later)
	! Find all bonds between atoms flagged as cofactor "conumber".
	count = 0
	do i = 1,nBonds
		if(glb_cofactor(bnd(i)%i).eq.conumber.and.glb_cofactor(bnd(i)%j).eq.conumber) then
			count = count +1
		end if
	end do

	cf%mol%num_bond = count
	allocate(cf%mol%bond(0:cf%mol%num_bond -1))
	count = 0
	do i = 1,nBonds
		if(glb_cofactor(bnd(i)%i).eq.conumber.and.glb_cofactor(bnd(i)%j).eq.conumber) then
			cf%mol%bond(count)%id			= count
			cf%mol%bond(count)%valid	= 1
			cf%mol%bond(count)%atom_1	= bnd(i)%i - cf%cofactor_offset
			cf%mol%bond(count)%atom_2	= bnd(i)%j - cf%cofactor_offset
			cf%mol%bond(count)%ttype	= SYBYL_bond_type(bonds(i)%cod)
			cf%mol%bond(count)%neib(:)	= 0
			count = count +1
		end if
	end do

	! Bubble sort bonds in increasing atom order
	do
		mark = 0
		do i = 0, cf%mol%num_bond -2
			if((cf%mol%bond(i)%atom_1.eq.cf%mol%bond(i+1)%atom_1).and. &
				(cf%mol%bond(i)%atom_2.eq.cf%mol%bond(i+1)%atom_2)) cycle
			if(cf%mol%bond(i)%atom_1 <  cf%mol%bond(i+1)%atom_1) cycle
			if((cf%mol%bond(i)%atom_1.eq.cf%mol%bond(i+1)%atom_1) .and. &
				(cf%mol%bond(i)%atom_2<=cf%mol%bond(i+1)%atom_2)) cycle

			! If no rule triggered, the bonds are in the wrong order
			mark = 1
			call Bond_Swap(cf%mol%bond(i),cf%mol%bond(i+1))
		end do
		if(mark.eq.0) exit
	end do

	! Reset bond id
	do i = 0,cf%mol%num_bond -1
		cf%mol%bond(i)%id = i +1
	end do
	
	deallocate(mask)
	deallocate(mask_res)
end subroutine xextract_cofactor

subroutine xscore_calc(iCalc, iFrame)	! calc topmost routine
	!arguments
	integer, intent(in)	:: iCalc				! calculation index, not used
	integer, intent(in)	:: iFrame				! frame index
	real(kind=prec)				:: score				! scoring results
	integer			:: i

	! Update coordinates
	do i = 0,nqat -1
		ligand%mol%atom(i)%coor(0) = xin(3*(iqseq(i+1)-1) +1)
		ligand%mol%atom(i)%coor(1) = xin(3*(iqseq(i+1)-1) +2)
		ligand%mol%atom(i)%coor(2) = xin(3*(iqseq(i+1)-1) +3)
	end do

	do i = 0,protein%num_atom -1
		protein%atom(i)%coor(0) = xin(3*(protein%atom(i)%topindex-1) +1)
		protein%atom(i)%coor(1) = xin(3*(protein%atom(i)%topindex-1) +2)
		protein%atom(i)%coor(2) = xin(3*(protein%atom(i)%topindex-1) +3)
	end do

	call Ligand_Value_Atom(ligand,1)
	call Protein_Value_Atom(protein, 1,1)		! do lite run
	call Protein_Define_Pocket(protein, ff,ligand, 10.0_prec)

	score = Ligand_Calculate_Binding_Score(ligand,input,protein)	

	call xlog_frame(iFrame,ligand)

	if(input%show_abs.eq.'YES')		call Ligand_Scoring_Stats(ligand)
	if(input%show_abs.eq.'QATOM') call Ligand_Scoring_Stats(ligand,input%qatom-1)

	if(input%show_total.eq.'YES') then
	    write(*,'(t28,f6.1,  t36,f5.2, t43,f5.1, t51,f5.2, t57,f5.1, t65,f4.1, t73,f6.3)') &
		ligand%vdw,ligand%hb,ligand%hp,ligand%hm,ligand%hs,ligand%rt, &
		ligand%bind_score
	endif
	if(input%show_abs.eq.'YES')		write(*,'(a)') '------------------------------------------------------------------------------'
end subroutine xscore_calc

subroutine xlog_frame(iFrame,ligand)				! logs scoring results for frame
	integer, intent(in)		:: iFrame						! current frame
	type(tXLigand)				:: ligand
	type(tXScore),pointer	:: new_aXScore(:)		! tmp pointer
	
	nXScores = nXScores +1
	if(nXScores>maxXScores) then
		maxXScores = maxXScores + 16										! allocate 16 elements a time
		allocate(new_aXScore(maxXScores))								! allocate new mem
		new_aXScore(1:nXScores-1) = xscores(1:nXScores-1)	! copy old
		
		deallocate(xscores)				! return old mem to OS
		xscores => new_aXScore			! update pointer
	end if

	! finally store values
	xscores(nXScores)%frame	= iFrame	! iFrame may be redundant
	xscores(nXScores)%vdw		= ligand%vdw
	xscores(nXScores)%hb		= ligand%hb
	xscores(nXScores)%rt		= ligand%rt
	xscores(nXScores)%hp		=	ligand%hp
	xscores(nXScores)%hs		= ligand%hs
	xscores(nXScores)%hm		= ligand%hm
	xscores(nXScores)%score	= ligand%bind_score
end subroutine xlog_frame

subroutine xscore_mean
	real(kind=prec)					:: vdw,hb,rt,hp,hs,hm,score
	integer					:: i
	
	do i = 1,nXScores
		vdw		= vdw		+ xscores(i)%vdw
		hb		= hb		+ xscores(i)%hb
		rt		= rt	  + xscores(i)%rt
		hp		= hp	  + xscores(i)%hp
		hs		= hs		+ xscores(i)%hs
		hm		= hm		+ xscores(i)%hm
		score = score + xscores(i)%score
	end do

	vdw		= vdw		/ nXScores
	hb		= hb		/ nXScores
	rt		= rt		/ nXScores
	hp		= hp		/ nXScores
	hs		= hs		/ nXScores
	hm		= hm		/ nXScores
	score = score / nXScores


	! output mean values
	write(*,1)  'processed frames: ', nXScores, 'VDW',   'HB',    'HP', &
	    'HM',      'HS',       'RT',			'SCORE'
1	format(      t1,a,                t19,i4, t31,a3,  t38,a2,  t46,a2,  t52,a2,   t60,a2, t67,a2,    t74,a5)
	write(*,'(a, t28,f6.1,  t36,f4.1, t43,f5.1, t50,f4.1, t57,f5.1, t65,f4.1, t74,f5.2)') 'mean score', vdw,hb,hp,hm,hs,rt,score

end subroutine xscore_mean

subroutine xscore_heading(i)
	integer, intent(in)		:: i

	write(*,1)  'VDW',   'HB',    'HP',    'HM',      'HS',       'RT',			'SCORE'
1	format(      t31,a3,  t38,a2,  t46,a2,  t52,a2,   t60,a2, t67,a2,   t74,a5)
end subroutine xscore_heading

! ==========================================================================================================================================
! ==========================================================================================================================================

! XMOLECULE MODULE

! ==========================================================================================================================================
! ==========================================================================================================================================

integer function Molecule_Get_Num_Heavy_Atom(molecule)
	type(tXMolecule)	:: molecule
	integer						:: i,num

	num = 0
	do i = 0,molecule%num_atom -1 
		 if(molecule%atom(i)%valid.eq.0) goto 1
		 if(molecule%atom(i)%ttype.eq.'H') goto 1
		 num = num +1 
1	end do

	Molecule_Get_Num_Heavy_Atom = num 
end function Molecule_Get_Num_Heavy_Atom

subroutine	Molecule_SurDot_Clear(molecule)
	type(tXMolecule)	:: molecule

	molecule%num_sur_dot = 0
	if(associated(molecule%sur_dot)) deallocate(molecule%sur_dot)
end subroutine Molecule_SurDot_Clear

subroutine Molecule_Generate_Surface_Dots(molecule,probe_r)
	type(tXMolecule)	:: molecule
	real(kind=prec)							:: probe_r
	integer						:: i,j,k,cc 
	integer						:: mark 
	real(kind=prec)							:: d,dd,dmin 
	type(tXDotSet)		:: tmp_set 
	type(tXDot)				:: tmp_dot 
	integer,allocatable ::  ligand_check_list(:)

	! molecule%sur_dot is a vector of type tXDot
	call Molecule_SurDot_Clear(molecule)

	allocate(ligand_check_list(0:molecule%num_atom -1),stat=err) 
	if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'

	do i = 0,molecule%num_atom -1 
		if(molecule%atom(i)%valid<=0) goto 1
		if(molecule%atom(i)%xtype.eq.'H') goto 1

		do j = 0,molecule%num_atom -1 
			if((i.eq.j).or.molecule%atom(j)%valid<=0) then
				ligand_check_list(j)=0
				goto 2
			elseif(molecule%atom(j)%xtype.eq.'H') then
				ligand_check_list(j)=0
				goto 2
			end if

			d=Distance(molecule%atom(i)%coor,molecule%atom(j)%coor) 

			if(d>(molecule%atom(i)%R+molecule%atom(j)%R+2*probe_r)) then
				ligand_check_list(j)=0 
			else
				ligand_check_list(j)=1 
			end if
2		end do
			
		if(associated(tmp_set%dot)) tmp_set%dot(:)%valid = 0
		tmp_set = ForceField_Get_Surface_Dot(molecule%atom(i),probe_r) 

		do j = 0,tmp_set%num_dot -1			! check whether this dot is on the surface or not 
			mark=1
			dmin=1000.0

			do k = 0,molecule%num_atom -1 
				if(ligand_check_list(k).eq.0) goto 3

				d=Distance(tmp_set%dot(j)%coor(0:2),molecule%atom(k)%coor(0:2)) 
				dd=d/(molecule%atom(k)%R+probe_r)-1.00 
				if(dd<0.000) then
					mark=0
					exit
				elseif(dd<dmin) then
					dmin=dd
					goto 3
				else
					goto 3
				end if
3			end do

			if(mark.eq.0) then	
				tmp_set%dot(j)%valid=0 
			elseif(dmin>=0.10) then 
				tmp_set%dot(j)%valid=1  ! regular 
			else 
				tmp_set%dot(j)%valid=2 	! dots at edge
			end if
		end do

		! now record the surface dots of the current atom
		cc = 0
		do j = 0,tmp_set%num_dot -1
			if(tmp_set%dot(j)%valid.eq.0) goto 4

			tmp_dot=tmp_set%dot(j)
			tmp_dot%valid=i+1

			if(tmp_set%dot(j)%valid.eq.1) then
				tmp_dot%unit=tmp_set%unit 
			else  ! correct the overlapping of edging dots
				tmp_dot%unit=tmp_set%unit*0.500 
			end if

			call Dot_PushBack(molecule%sur_dot, molecule%num_sur_dot, tmp_dot)
			cc = cc +1
4		end do

		deallocate(tmp_set%dot)
1	end do

	deallocate(ligand_check_list) 
end subroutine Molecule_Generate_Surface_Dots

subroutine Molecule_Assign_Apparent_Charge(molecule)
	type(tXMolecule)	:: molecule
	integer						:: i 

	molecule%charge_type = 'FORMAL_CHARGES'

	do i = 0,molecule%num_atom -1 
		if(molecule%atom(i)%xtype.eq.'O.co2') then
			molecule%atom(i)%q=-0.500 
		elseif(molecule%atom(i)%xtype.eq.'N.4') then
			molecule%atom(i)%q=1.000 
		else 
			molecule%atom(i)%q=0.000 
		end if
	end do
	
	do i = 0,molecule%num_atom -1 
		 if(molecule%atom(i)%xtype.ne.'N.pl3.h') goto 1
		 if(molecule%atom(i)%num_nonh.ne.1) goto 1
		 if(molecule%atom(molecule%atom(i)%neib(0)-1)%ttype.ne.'C.cat') goto 1
		 molecule%atom(i)%q=0.500 
1	end do
end subroutine Molecule_Assign_Apparent_Charge

subroutine Molecule_Show_Contents(mol)
	type(tXMolecule) :: mol

	write(*,'(a,a,a)')  'Molecule: ', trim(mol%name)
	write(*,'(a,f9.3)') 'MW = ', mol%weight 

	write(*,'(a,i5)') 'Number of HB atoms = ', mol%num_hb_atom 
	write(*,'(a,i5)') 'Number of rotors    = ', mol%num_rotor 

	write(*,'(a,f6.2)') 'LogP    = ', mol%logp 

	call Molecule_Show_Atoms(mol) 
	call Molecule_Show_Bonds(mol) 
	call Molecule_Show_Rings(mol) 
end subroutine Molecule_Show_Contents

subroutine Molecule_Show_Atoms(mol)
	type(tXMolecule) :: mol
	integer ::  i 

	write(*,'(a,i5)') 'Total number of atoms in this molecule = ',mol%num_atom
	write(*,'(t1,a)') ' ID  topoID  vld. type   xtype   type2   res r_id  &
	    &name   weight  r   eps   q   XR  logp solv hb  ring origin part &
	    &#neib #nonh  neib[0..6].................  bond[0..6]................  hb root[0..2]'

	write(*,'()')
	do i = 0,mol%num_atom -1
		write(*,'(i4, 2x, i5, 3x, i1, t20, a, t27, a, t36, a, t42, a, t47, a, t55, a, t59, f5.2, 6(1x,f4.1), &
			& 2x, a, t96, 5(4x,i2), 1x, 14(1x,i3), 3(1x,f5.1))'), &

			mol%atom(i)%id, &
			mol%atom(i)%topindex, &
			mol%atom(i)%valid, &
			mol%atom(i)%ttype, &
			mol%atom(i)%xtype, &
			mol%atom(i)%type2, &
			mol%atom(i)%residue, &
			adjustl(trim(mol%atom(i)%res_id)), &
			mol%atom(i)%name, &
			mol%atom(i)%weight, &
			mol%atom(i)%r, &
			mol%atom(i)%eps, &
			mol%atom(i)%q, &
			mol%atom(i)%XR, &
			mol%atom(i)%logp, &
			mol%atom(i)%solv, &
			mol%atom(i)%hb, &
!			mol%atom(i)%occupancy,&
!			mol%atom(i)%bfactor,&
			mol%atom(i)%ring,&
			mol%atom(i)%origin,&
			mol%atom(i)%part,&
			mol%atom(i)%num_neib,&
			mol%atom(i)%num_nonh,&
			mol%atom(i)%neib(0),mol%atom(i)%neib(1),&
			mol%atom(i)%neib(2),mol%atom(i)%neib(3),&
			mol%atom(i)%neib(4),mol%atom(i)%neib(5),&
			mol%atom(i)%neib(6), &
			mol%atom(i)%bond(0),mol%atom(i)%bond(1),&
			mol%atom(i)%bond(2),mol%atom(i)%bond(3),&
			mol%atom(i)%bond(4),mol%atom(i)%bond(5),&
			mol%atom(i)%bond(6), &
			mol%atom(i)%root(0),mol%atom(i)%root(1),mol%atom(i)%root(2)
	end do
end subroutine Molecule_Show_Atoms

subroutine Molecule_Show_Bonds(mol)
	type(tXMolecule) :: mol
	integer ::  i 

	write(*,'(a,i6)') 'Total number of bonds in this molecule = ',mol%num_bond
	write(*,'(a)') '  ID valid atom_1 atom_2 topo_1 topo_2 type ring part  num_neib  neib[0:5]'

	do i = 0,mol%num_bond -1 
		write(*,'(i4, 1x, i3, 2x,i4, 3x, i4, 4x, i5, 2x i5, t42 ,a, t46, i2, 1x, i3, 4x, i4, 2x, 6(i5,1x))')  &
		 			mol%bond(i)%id, mol%bond(i)%valid, mol%bond(i)%atom_1, mol%bond(i)%atom_2, &
					mol%atom(mol%bond(i)%atom_1-1)%topindex, mol%atom(mol%bond(i)%atom_2-1)%topindex, &
					mol%bond(i)%ttype, mol%bond(i)%ring, mol%bond(i)%part, mol%bond(i)%num_neib, mol%bond(i)%neib(0), mol%bond(i)%neib(1), &
					mol%bond(i)%neib(2), mol%bond(i)%neib(3), mol%bond(i)%neib(4), mol%bond(i)%neib(5)
	end do
end subroutine Molecule_Show_Bonds

subroutine Molecule_Show_Rings(mol)
	type(tXMolecule) :: mol
	integer ::  i 

	write(*,'(a,i4)') 'Total number of aromatic rings in this molecule = ', size(mol%ring)

	if(size(mol%ring)>0 .and. associated(mol%ring)) then
		write(*,'(t1,a)') 'id  vld. type size    atoms (first-last)  centroid'
		do i = 0,size(mol%ring) -1 
!			write(*,*) i, associated(mol%ring), size(mol%ring)
			if(mol%ring(i)%valid.ne.0) then
				write(*,'(t1,i3, t6,i1, t10,i1, t15,i2, t23,i5,t29,i5, t43,f5.1,t49,f5.1,t56,f5.1)') &
						i+1, mol%ring(i)%valid,mol%ring(i)%ttype,mol%ring(i)%num_member, &
						mol%ring(i)%atom_id(0),mol%ring(i)%atom_id(mol%ring(i)%num_member-1), &
						mol%ring(i)%centroid(0),mol%ring(i)%centroid(1),mol%ring(i)%centroid(2)
			else
				write(*,'(t1,i3, t6,a)') i+1, 'invalid'
			end if
		end do
	end if
end subroutine Molecule_Show_Rings

subroutine Molecule_DetectConnections(molecule)
	type(tXMolecule),intent(inout)	:: molecule
	integer ::  i,j,k,mark,tmp1,tmp2,old_part,new_part,id1,id2 
	integer,allocatable ::  part_list(:)

	! task 1: detect all the fragments and determine num_subst 
	! first, detect how many fragments are in the molecule
	! initialize the molecule as the assembly of isolated atoms
	do i = 0,molecule%num_atom -1
		molecule%atom(i)%part=molecule%atom(i)%id 
	end do
	
	!  then merge the atoms into fragments 
	do i = 0,molecule%num_bond -1 
		if(molecule%bond(i)%valid<=0) goto 1 

		tmp1=molecule%bond(i)%atom_1
		tmp2=molecule%bond(i)%atom_2 
!		write(*,*) tmp1, tmp2, size(molecule%atom)
		if(molecule%atom(tmp1-1)%part>=molecule%atom(tmp2-1)%part) then
			 old_part=molecule%atom(tmp1-1)%part  
			 new_part=molecule%atom(tmp2-1)%part 
		else  
			 old_part=molecule%atom(tmp2-1)%part  
			 new_part=molecule%atom(tmp1-1)%part 
		end if

		do j = 0,molecule%num_atom -1 
			 if(molecule%atom(j)%part.ne.old_part) then
				 goto 2 
			 else
				molecule%atom(j)%part=new_part 
			 end if 
2		end do
1	end do

	! then count the number of fragments 
	! also re-arrange the id of fragments to make them continuous
	allocate(part_list(0:molecule%num_atom -1),stat=err)
	if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'

	part_list(:)=0

	molecule%num_subst=0 
	do i = 0,molecule%num_atom -1 
		 mark=0 

		 do j = 0,molecule%num_subst -1 
			 if(molecule%atom(i)%part.ne.part_list(j)) then
				goto 3 
			 else 
				mark=1  
				exit
			 end if
3		 end do

		 if(mark==1) then
			molecule%atom(i)%part=j+1 	! an existing fragment 
		 else					! a new fragment 
			 part_list(molecule%num_subst)=molecule%atom(i)%part 
			 molecule%num_subst = molecule%num_subst +1 
			 molecule%atom(i)%part=molecule%num_subst 	! 1,2,3 ...
		 end if
	 end do

	deallocate(part_list)

	! also define all the bonds
	do i = 0,molecule%num_bond -1  
		 if(molecule%bond(i)%valid<=0) then
			molecule%bond(i)%part=0 
		 else
			molecule%bond(i)%part=molecule%atom(molecule%bond(i)%atom_1-1)%part 
		 end if
	end do

	! task 2: set up the connection tables for atoms and bonds
	! assign atomic weight to atoms, which is necessary to rank 
	! the neighboring atoms for a central atom
	do i = 0,molecule%num_atom -1
		 if(molecule%atom(i)%ttype.eq.'F') then
				molecule%atom(i)%weight=19.00 
     elseif(molecule%atom(i)%ttype.eq.'Cl')  then
				molecule%atom(i)%weight=35.45 
     elseif(molecule%atom(i)%ttype.eq.'Br')  then
				molecule%atom(i)%weight=79.90 
		 elseif(molecule%atom(i)%ttype.eq.'I')  then
				molecule%atom(i)%weight=126.90 
		 elseif(molecule%atom(i)%ttype.eq.'Si')  then
				molecule%atom(i)%weight=28.09 
		 elseif(molecule%atom(i)%ttype(1:1).eq.'C')  then
				molecule%atom(i)%weight=12.01 
		 elseif(molecule%atom(i)%ttype(1:1).eq.'H')  then
				molecule%atom(i)%weight=1.00 
		 elseif(molecule%atom(i)%ttype(1:1).eq.'N')  then
				molecule%atom(i)%weight=14.01 
		 elseif(molecule%atom(i)%ttype(1:1).eq.'O')  then
				molecule%atom(i)%weight=16.00 
		 elseif(molecule%atom(i)%ttype(1:1).eq.'P')  then
				molecule%atom(i)%weight=30.97 
		 elseif(molecule%atom(i)%ttype(1:1).eq.'S') then
				molecule%atom(i)%weight=32.06 
		 else 
				molecule%atom(i)%weight=0.00 
		 end if
  end do

	! now detect the neighboring atoms for each atom
	do i = 0,molecule%num_atom -1
		molecule%atom(i)%num_neib=0
		molecule%atom(i)%num_nonh=0 
  end do
	do i = 0,molecule%num_bond -1 
		 if(molecule%bond(i)%valid<=0) goto 4 
		 id1=molecule%bond(i)%atom_1
		 id2=molecule%bond(i)%atom_2 
		 tmp1=molecule%atom(id1-1)%num_neib  
		 if(tmp1<MAX_ATOM_NEIB) then
		 	 molecule%atom(id1-1)%neib(tmp1)=id2  
			 molecule%atom(id1-1)%bond(tmp1)=i+1  
		 	 molecule%atom(id1-1)%num_neib = molecule%atom(id1-1)%num_neib +1
		 end if

		 tmp2=molecule%atom(id2-1)%num_neib  

		 if(tmp2<MAX_ATOM_NEIB) then
		 	 molecule%atom(id2-1)%neib(tmp2)=id1  
			 molecule%atom(id2-1)%bond(tmp2)=i+1 
		 	 molecule%atom(id2-1)%num_neib = molecule%atom(id2-1)%num_neib +1 
		 end if
4	 end do

	! now arrange the neighboring atoms in a decreasing order
	! according to atomic weights  bonds are also re-arranged
	do i = 0,molecule%num_atom -1 
		 if(molecule%atom(i)%valid==0) then
			goto 5 
		 elseif(molecule%atom(i)%num_neib<=1) then
			goto 5 
		 end if
	
		do j = 0,molecule%atom(i)%num_neib-1 -1 
		 do k=j+1,molecule%atom(i)%num_neib -1
			 tmp1=molecule%atom(i)%neib(j)
			 tmp2=molecule%atom(i)%neib(k) 
			 if(molecule%atom(tmp1-1)%weight>=molecule%atom(tmp2-1)%weight) then
				 goto 6 
			 else
				 mark=molecule%atom(i)%neib(j) 
				 molecule%atom(i)%neib(j)=molecule%atom(i)%neib(k) 
				 molecule%atom(i)%neib(k)=mark 
				 mark=molecule%atom(i)%bond(j) 
         molecule%atom(i)%bond(j)=molecule%atom(i)%bond(k) 
         molecule%atom(i)%bond(k)=mark 
				end if
6			end do
		end do
5	end do
	do i = 0,molecule%num_atom -1 
		do j = 0,molecule%atom(i)%num_neib -1 
			 if(molecule%atom(molecule%atom(i)%neib(j)-1)%ttype(1:1)=='H') then
				 goto 7 
		 	 else
				 molecule%atom(i)%num_nonh = molecule%atom(i)%num_nonh +1 
			 end if
7		 end do
	end do

	! now detect the neighboring bonds
	molecule%bond(:)%num_neib=0 
	do i = 0,molecule%num_bond-1
		molecule%bond(i)%neib(:)=0
	end do
	 	
	do i = 0,molecule%num_bond-1 -1
		do j=i+1,molecule%num_bond  -1
			 if((molecule%bond(i)%valid<=0).or.(molecule%bond(j)%valid<=0)) goto 8 

			 if(Molecule_Two_Bonds_Connection_Check(molecule%bond(i),molecule%bond(j)).ne.0) then
				 tmp1=molecule%bond(i)%num_neib 
				 if(tmp1<MAX_BOND_NEIB) then
			 		 molecule%bond(i)%neib(tmp1)=molecule%bond(j)%id  
					 molecule%bond(i)%num_neib = molecule%bond(i)%num_neib +1 
				 end if

				 tmp2=molecule%bond(j)%num_neib 
				 if(tmp2<MAX_BOND_NEIB) then
			 		 molecule%bond(j)%neib(tmp2)=molecule%bond(i)%id  
					 molecule%bond(j)%num_neib = molecule%bond(j)%num_neib +1
				 end if
			 else
				goto 8
			 end if
8		 end do
	 end do
end subroutine Molecule_DetectConnections

integer function Molecule_Two_Bonds_Connection_Check(bond1, bond2)
	! this function returns the ID of the joint atom
  type(tXBond)		:: bond1, bond2
	integer					:: id 

  id=bond1%atom_1 

  if(id.eq.bond2%atom_1) then
		Molecule_Two_Bonds_Connection_Check = id
		return 
  elseif(id.eq.bond2%atom_2) then
		Molecule_Two_Bonds_Connection_Check = id 
		return
	end if

  id=bond1%atom_2 

  if(id==bond2%atom_1) then
		Molecule_Two_Bonds_Connection_Check = id 
		return
  elseif(id==bond2%atom_2) then
		Molecule_Two_Bonds_Connection_Check = id 
		return
	end if

  Molecule_Two_Bonds_Connection_Check = 0        ! two bonds are not connected
end function Molecule_Two_Bonds_Connection_Check

integer function Molecule_ValueAtom(molecule,lite_run)
	type(tXMolecule),intent(inout)	:: molecule
	integer, optional								:: lite_run
	integer :: i
	integer :: mark,success
	character(len=20) :: xtype

	success = 1
	call Molecule_DetectConnections(molecule)

	! check if this molecule is valid or not
	!success = Molecule_CheckAtomType(molecule)
	success = 1

	! now detect the ring systems in this molecule
	call Molecule_Detect_Rings(molecule)
	call Molecule_Detect_Aromatic_Rings(molecule)

	! now assign XTOOL atom type and corresponding parameters
	do i=0,molecule%num_atom-1
		 mark = Molecule_Get_XTOOL_Type(molecule,molecule%atom(i)%id,xtype)

	 	 if(mark.eq.1) then
			 molecule%atom(i)%xtype = xtype
			 molecule%atom(i)%valid=1
	 	 else 
			 molecule%atom(i)%xtype = xtype
			 molecule%atom(i)%valid=0
			 success=0
		 end if
	end do
	do i = 0,molecule%num_atom-1
	 	 if(molecule%atom(i)%valid<=0) then
			goto 1
		 end if

	 	 mark = ForceField_Assign_Atom_Parameters(molecule%atom(i))
						
	 	 if(mark==0) then 
			molecule%atom(i)%valid=0
			success=0
		 else
			goto 1
		 end if
1	end do
	
	molecule%weight = 0
	do i = 0,molecule%num_atom -1
		molecule%weight = molecule%weight + molecule%atom(i)%weight
	end do

	! now return value
	if(success==1) then
		Molecule_ValueAtom = 1
	else
		Molecule_ValueAtom = 0
	end if
end function Molecule_ValueAtom

integer function Molecule_Get_XTOOL_Type(molecule,atom_id, xtype)
	type(tXMolecule) :: molecule
	integer	:: atom_id
	character(len=20) :: xtype
	type(tXGroup) :: group 
	type(tXAtom) :: atom 

	atom = molecule%atom(atom_id-1)			! check this
	group = Molecule_Find_A_Group(molecule,atom%id)

	xtype = 'Un'

	if((atom%ttype.eq.'H').or.(atom%ttype.eq.'H.spc')) then
	 	if(group%neib(0)%ttype(1:1).eq.'O') then
			xtype = 'H.hb'
	 	elseif(group%neib(0)%ttype(1:1).eq.'N') then
			xtype = 'H.hb'
	 	else 
			xtype = 'H'  
		end if
	end if
	if((atom%ttype.eq.'C.3')) then
		 if(group%num_hetero.eq.0) then
			xtype = 'C.3'
		 elseif(group%num_hetero>0) then
			xtype = 'C.3.x'
		 else 
			xtype = 'C.3.un'
		end if
	end if
	if((atom%ttype.eq.'C.2').and.(atom%ring.ne.2)) then
		 if(group%num_hetero.eq.0) then
			xtype = 'C.2'
		 elseif(group%num_hetero>0) then
			xtype = 'C.2.x'
		 else 
			xtype = 'C.2.un'
		 end if
	end if

	if((atom%ttype.eq.'C.ar').or.((atom%ttype.eq.'C.2').and.atom%ring.eq.2)) then
		 if(group%num_hetero.eq.0) then
			xtype = 'C.ar'
		 elseif(group%num_hetero>0) then
			xtype = 'C.ar.x'
		 else 
			xtype = 'C.ar.un'
		 end if
	end if
	
	if((atom%ttype.eq.'C.1')) then
		 if(group%num_hetero.eq.0) then
			xtype = 'C.1'
		 elseif(group%num_hetero>0) then
			xtype = 'C.1.x'
		 else
			xtype = 'C.1.un'
		 end if
	end if

	if((atom%ttype.eq.'C.cat')) xtype = 'C.cat'

	if((atom%ttype.eq.'N.4').or.(atom%ttype.eq.'N.3')) then
		 if(group%num_nonh<=2) then
			xtype = 'N.4'
		 elseif(group%num_nonh.eq.3) then
			xtype = 'N.3'
		 else
			xtype = 'N.3.un'
		 end if
	end if

	if(((atom%ttype.eq.'N.am').and.atom%ring.ne.2).or.((atom%ttype.eq.'N.pl3').and.atom%ring.ne.2)) then
		 if(group%num_nonh.eq.1) then
			xtype = 'N.pl3.h'
		 elseif(group%num_nonh.eq.2) then
			xtype = 'N.pl3.h'
		 elseif(group%num_nonh.eq.3) then
			xtype = 'N.pl3' 
		 else
			xtype = 'N.pl3.un'
		 end if
	end if

	if((atom%ttype.eq.'N.2').and.atom%ring.ne.2) then
		 if(group%num_nonh.eq.1) then	
			xtype = 'N.2.h'
		 elseif(group%num_nonh.eq.2) then
			xtype = 'N.2'
		 else 
			xtype = 'N.2.un'
		 end if
	end if
	if((atom%ttype.eq.'N.ar').or.((atom%ttype.eq.'N.2').and.atom%ring.eq.2).or. &
	  ((atom%ttype.eq.'N.pl3').and.atom%ring.eq.2).or.((atom%ttype.eq.'N.am').and.atom%ring.eq.2)) then
		 if(group%num_h.eq.1) then
			xtype = 'N.ar.h'
		 elseif(group%num_h.eq.0) then
			xtype = 'N.ar'
		 else 
			xtype = 'N.ar.un'
		 end if
	end if

	if((atom%ttype.eq.'N.1')) then
		 if(group%num_nonh.eq.1) then
			xtype = 'N.1'
		 else 
			xtype = 'N.1.un'
		 end if
	end if

	if((atom%ttype.eq.'O.3')) then
		 if(group%num_nonh.eq.1) then
			xtype = 'O.3.h'
		 elseif(group%num_nonh.eq.2) then
			xtype = 'O.3'
		 else 
			xtype = 'O.3.un'
		 end if
	end if

	if((atom%ttype.eq.'O.2')) xtype = 'O.2'
	if((atom%ttype.eq.'O.co2')) xtype = 'O.co2'

	if((atom%ttype.eq.'S.3')) then
		 if(group%num_nonh.eq.1) then
			xtype = 'S.3.h'
		 elseif(group%num_nonh.eq.2) then
			xtype = 'S.3'
		 else 
			xtype = 'S.3.un'
		 end if
	end if

	if((atom%ttype.eq.'S.2'))		xtype = 'S.2'
	if((atom%ttype.eq.'S.o'))		xtype = 'S.o'
	if((atom%ttype.eq.'S.o2'))	xtype = 'S.o'
	if((atom%ttype.eq.'P.3'))		xtype = 'P.3'
	if((atom%ttype.eq.'F'))			xtype = 'F'
	if((atom%ttype.eq.'Cl'))		xtype = 'Cl'
	if((atom%ttype.eq.'Br'))		xtype = 'Br'
	if((atom%ttype.eq.'I'))			xtype = 'I'
	if((atom%ttype.eq.'Si'))		xtype = 'Si'
	
	if((atom%ttype.eq.'Na'))		xtype =	'M+'

	if((atom%ttype.eq.'Na'))		xtype =	'M+'
	if((atom%ttype.eq.'K'))			xtype =	'M+'
	if((atom%ttype.eq.'Ca'))		xtype =	'M+'
	if((atom%ttype.eq.'Li'))		xtype =	'M+'
	if((atom%ttype.eq.'Al'))		xtype =	'M+'
	if((atom%ttype.eq.'Fe'))		xtype =	'M+'

	if((atom%ttype.eq.'Un')) then
		Molecule_Get_XTOOL_Type = 0
	else
	  Molecule_Get_XTOOL_Type = 1
	end if
end function Molecule_Get_XTOOL_Type

type(tXGroup) function Molecule_Find_A_Group(molecule,atom_id)
	type(tXMolecule) :: molecule
	integer	:: atom_id
  integer ::  i,j,id,num 
	integer ::  mark 
  type(tXGroup) :: group 

	! define the center
  group%center = molecule%atom(atom_id-1) 

  ! find the center's neighbours
  num=group%center%num_neib 

	do i = 0,num -1 
		group%neib(i)=molecule%atom(group%center%neib(i)-1) 
		group%bond(i)=molecule%bond(group%center%bond(i)-1) 
  end do
    
	! count the necessary parameters
  group%num_neib=group%center%num_neib 
  group%num_h=0
	group%num_nonh=0

	do i = 0,num -1 
		if(group%neib(i)%ttype(1:1).eq.'H') then
			group%num_h = group%num_h +1 
		else
			group%num_nonh = group%num_nonh +1 
		end if
  end do

  group%num_hetero=0 

	do i = 0,group%num_nonh -1 
		if(group%neib(i)%ttype.eq.'F') then
			group%num_hetero = group%num_hetero +1 
		elseif(group%neib(i)%ttype.eq.'Cl') then
			group%num_hetero = group%num_hetero +1 
		elseif(group%neib(i)%ttype.eq.'Br') then
			group%num_hetero = group%num_hetero +1  
		elseif(group%neib(i)%ttype.eq.'I') then
			group%num_hetero = group%num_hetero +1 
		elseif(group%neib(i)%ttype.eq.'Si') then
			goto 1  
		elseif(group%neib(i)%ttype(1:1).eq.'N') then
			group%num_hetero = group%num_hetero +1 
		elseif(group%neib(i)%ttype(1:1).eq.'O') then
			group%num_hetero = group%num_hetero +1 
		elseif(group%neib(i)%ttype(1:1).eq.'P') then
			group%num_hetero = group%num_hetero +1 
		elseif(group%neib(i)%ttype(1:1).eq.'S') then
			group%num_hetero = group%num_hetero +1 
		else
			goto 1 
		end if
1	end do

  group%db_type=0
	group%num_db=0
	group%num_tb=0

	do i = 0,group%num_nonh -1 
		if(group%bond(i)%ttype.eq.'2') then
			group%num_db = group%num_db +1 
			if(group%neib(i)%ttype(1:1).eq.'C') then
				group%db_type=1 
			else if(group%neib(i)%ttype(1:1).eq.'N') then
				group%db_type=2 
			else if(group%neib(i)%ttype(1:1).eq.'O') then
				group%db_type=3 
			else if(group%neib(i)%ttype(1:1).eq.'S') then 
				group%db_type=4 
			else if(group%neib(i)%ttype(1:1).eq.'P') then
				group%db_type=5 
			else
				goto 2 
			end if
		elseif((group%bond(i)%ttype.eq.'1').and.(group%neib(i)%ttype.eq.'O.co2.')) then
			group%db_type=3
			group%num_db = group%num_db +1
		elseif((group%bond(i)%ttype.eq.'ar').and.(group%neib(i)%ttype.eq.'O.co2')) then
	    group%db_type=3
		  group%num_db = group%num_db +1
    elseif(group%bond(i)%ttype.eq.'3') then
	    group%num_tb = group%num_tb +1 
    else
			goto 2 
		end if
2 end do
  
	group%num_pi=0 

	do i = 0,group%num_nonh -1 
		if(group%bond(i)%ttype.eq.'2') then 
			goto 3 
		elseif(group%bond(i)%ttype.eq.'3') then 
			goto 3 
		elseif(group%neib(i)%ttype.eq.'C.ar') then 
			group%num_pi = group%num_pi +1 
		elseif(group%neib(i)%ttype.eq.'C.2') then 
			group%num_pi = group%num_pi +1 
		elseif(group%neib(i)%ttype.eq.'C.1') then 
			group%num_pi = group%num_pi +1 
		elseif(group%neib(i)%ttype.eq.'C.cat') then 
			group%num_pi = group%num_pi +1 
		elseif(group%neib(i)%ttype.eq.'N.2') then 
			group%num_pi = group%num_pi +1 
		elseif(group%neib(i)%ttype.eq.'N.1') then 
			group%num_pi = group%num_pi +1 
		elseif(group%neib(i)%ttype.eq.'N.ar') then 
			group%num_pi = group%num_pi +1 
		else 
			goto 3 
		end if
3	end do

	! check if the central atom is adjacent to any -SO-, -PO-, or -CO- 
	mark=0 
	do i = 0,group%num_nonh -1 
		if((group%neib(i)%ttype.eq.'P.3').or.(group%neib(i)%ttype.eq.'S.o').or. &
		    (group%neib(i)%ttype.eq.'S.o2').or.(group%neib(i)%ttype.eq.'C.2')) then
			 
			num=group%neib(i)%num_nonh  

			do j = 0,num -1 
				 id=molecule%atom(group%neib(i)%id-1)%neib(j) 
				 if(id.eq.group%center%id) then
					goto 4 
				 elseif(molecule%atom(id-1)%ttype.eq.'O.2') then
					 mark=1
					 exit
				 elseif(molecule%atom(id-1)%ttype.eq.'O.co2') then
					 mark=1
					 exit
				 else
					goto 4 
				 end if
4			 end do

			 if(mark.eq.1) then
				 exit
			 else
				goto 5 
			 end if
    else
			goto 5 
		end if
5	end do

	if(mark.eq.0) then
		group%amide=0 
	else 
		group%amide=j+1    ! assign the value of the amide bond
	end if
  
	group%valid=1
	Molecule_Find_A_Group = group

end function Molecule_Find_A_Group

subroutine Molecule_Detect_Aromatic_Rings(molecule)
	type(tXMolecule)	:: molecule
	integer ::  i,j ,slask
	type(tXRing) :: tmp_ring 
	integer,pointer ::  atom_path(:)
	integer,pointer::  bond_path(:)
	
	allocate(atom_path(0:6),stat=err)
	if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'
	allocate(bond_path(0:6),stat=err)
	if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'

	allocate(tmp_ring%atom_id(0:6),stat=err)
	if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'
	allocate(tmp_ring%bond_id(0:6),stat=err)
	if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'

	! Clear old ring data
	if(molecule%num_ring>0) then  !if(associated(molecule%ring)) then
		do i = 0,size(molecule%ring) -1
			call Ring_Clear(molecule%ring(i))
		end do
		deallocate(molecule%ring)
	end if
	if(associated(molecule%ring)) deallocate(molecule%ring)
	molecule%ring_count = 0
	molecule%num_ring = 0
	
! check all the 6-membered rings first 
	do i = 0,molecule%num_bond -1 
		if(molecule%bond(i)%ring.eq.1) then
			atom_path(:) = 0
			bond_path(:) = 0 

			if(Molecule_Look_For_A_Ring(molecule, molecule%bond(i)%id,atom_path,bond_path,6).eq.0) goto 1
			if(Molecule_Aromatic_Ring_Check_6(molecule,atom_path,bond_path).eq.0) goto 1 

			call Ring_Clear(tmp_ring)
		
			do j = 0,6 -1 
				molecule%atom(atom_path(j)-1)%ring=2 
				molecule%bond(bond_path(j)-1)%ring=2 

				!tmp_ring.atom_id.push_back(atom_path(j)) 
				call Int_PushBack(tmp_ring%atom_id,atom_path(j))

				!tmp_ring.bond_id.push_back(bond_path(j)) 
				call Int_PushBack(tmp_ring%bond_id,bond_path(j))

				tmp_ring%centroid(0)=tmp_ring%centroid(0)+molecule%atom(atom_path(j)-1)%coor(0) 
				tmp_ring%centroid(1)=tmp_ring%centroid(1)+molecule%atom(atom_path(j)-1)%coor(1) 
				tmp_ring%centroid(2)=tmp_ring%centroid(2)+molecule%atom(atom_path(j)-1)%coor(2) 
			end do

			tmp_ring%num_member=size(tmp_ring%atom_id)

			if(tmp_ring%num_member>0) then
				tmp_ring%centroid(0)=tmp_ring%centroid(0)/tmp_ring%num_member 
				tmp_ring%centroid(1)=tmp_ring%centroid(0)/tmp_ring%num_member 
				tmp_ring%centroid(2)=tmp_ring%centroid(0)/tmp_ring%num_member 

				tmp_ring%valid=1
				tmp_ring%ttype=2 

				Call Molecule_Add_Ring(molecule, tmp_ring)
				molecule%ring_count = molecule%ring_count +1
			end if
		end if
1	end do

	! then check all the 5-membered rings
	do i = 0,molecule%num_bond -1 
 		if(molecule%bond(i)%ring.ne.1) goto 2 

		atom_path(:)=0
		bond_path(:)=0 
	
		if(Molecule_Look_For_A_Ring(molecule,molecule%bond(i)%id,atom_path,bond_path,5).eq.0) goto 2 
		if(Molecule_Aromatic_Ring_Check_5(molecule,atom_path,bond_path).eq.0) goto 2 

		call Ring_Clear(tmp_ring)

		do j = 0,5 -1 
			molecule%atom(atom_path(j)-1)%ring=2 
			molecule%bond(bond_path(j)-1)%ring=2 

			!tmp_ring%atom_id%push_back(atom_path(j)) 
			call Int_PushBack(tmp_ring%atom_id,atom_path(j))

			!tmp_ring%bond_id%push_back(bond_path(j)) 
			call Int_PushBack(tmp_ring%bond_id,bond_path(j))

			tmp_ring%centroid(0)=tmp_ring%centroid(0)+molecule%atom(atom_path(j)-1)%coor(0) 
			tmp_ring%centroid(1)=tmp_ring%centroid(1)+molecule%atom(atom_path(j)-1)%coor(1) 
			tmp_ring%centroid(2)=tmp_ring%centroid(2)+molecule%atom(atom_path(j)-1)%coor(2) 
		end do

		tmp_ring%num_member=size(tmp_ring%atom_id)
		
		if(tmp_ring%num_member>0) then
			tmp_ring%centroid(0)=tmp_ring%centroid(0)/tmp_ring%num_member 
			tmp_ring%centroid(1)=tmp_ring%centroid(1)/tmp_ring%num_member 
			tmp_ring%centroid(2)=tmp_ring%centroid(2)/tmp_ring%num_member 

			tmp_ring%valid=1
			tmp_ring%ttype=2 
	
			Call Molecule_Add_Ring(molecule, tmp_ring)
			molecule%ring_count = molecule%ring_count +1
		end if
2	end do

	!this->num_ring=this->ring%size() 
	molecule%num_ring=size(molecule%ring)

	deallocate(atom_path)
	deallocate(bond_path)

	call Ring_Clear(tmp_ring)
end subroutine Molecule_Detect_Aromatic_Rings

subroutine Molecule_Add_Ring(m,r)
	type(tXMolecule)			:: m
	type(tXRing)					:: r
	integer								:: i
	type(tXRing),pointer	:: tmp(:)

	if(.not. associated(m%ring)) then
		! Adding first ring to molecule, need to init pointer
		allocate(m%ring(0:0),stat=err)
		if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'

		Call Ring_Copy(m%ring(0),r)
	else
		! Adding a second (or more) ring to molecule
	
		i = size(m%ring)
				
		! Increase size of ring-array
		allocate(tmp(0:i-1),stat=err)
		if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'
		do i  = 0,i -1
			call Ring_Copy(tmp(i),m%ring(i))
		end do

		deallocate(m%ring)
		allocate(m%ring(0:i),stat=err)
		if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'
		do i = 0, i-1
			call Ring_Copy(m%ring(i),tmp(i))
			call Ring_Clear(tmp(i))
		end do

		deallocate(tmp)

		! Add ring
		call Ring_Copy(m%ring(i), r)
	end if
end subroutine Molecule_Add_Ring

integer function Molecule_Aromatic_Ring_Check_6(molecule,atom_path,bond_path)
	type(tXMolecule)	:: molecule
	integer,pointer	::	atom_path(:)
	integer,pointer	::	bond_path(:)
	integer ::  i,count 
	integer ::  mark,mark1,mark2 

	! check the number of pi electrons 

	mark=1
	count=0

	do i = 0,6 -1 
		 if(molecule%atom(atom_path(i)-1)%ttype.eq.'C.ar') then
			 count = count +1  
		 elseif(molecule%atom(atom_path(i)-1)%ttype.eq.'N.ar') then
			 count = count +1  
		 elseif(molecule%atom(atom_path(i)-1)%ttype.eq.'C.2') then
			 count = count +1  
		 elseif(molecule%atom(atom_path(i)-1)%ttype.eq.'N.2') then
			 count = count +1  
		 else
			 mark=0
			 exit
		 end if
	end do

	if((mark.eq.0).or.(count.ne.6)) then
		Molecule_Aromatic_Ring_Check_6 = 0
		return
	end if 

	! check if there are two continuous single bonds on the path
	! this algorithm is superior to check -=-=-=
	mark=0 
	do i = 0,5 -1 
		if((molecule%bond(bond_path(i)-1)%ttype.eq.'1').or.(molecule%bond(bond_path(i)-1)%ttype.eq.'am')) then
			mark1=1 
		else
			mark1=0 
		end if
		if((molecule%bond(bond_path(i+1)-1)%ttype.eq.'1').or.(molecule%bond(bond_path(i+1)-1)%ttype.eq.'am')) then
			mark2=1 
		else
			mark2=0 
		end if
		if(mark1.ne.0.and.mark2.ne.0) then
			mark=1
			exit
		else
			goto 1 
		end if
1	end do
	if(mark.ne.0) then
		Molecule_Aromatic_Ring_Check_6 = 0
	else
		Molecule_Aromatic_Ring_Check_6  = 1
	end if
end function Molecule_Aromatic_Ring_Check_6

integer function Molecule_Aromatic_Ring_Check_5(molecule,atom_path, bond_path)
	type(tXMolecule)	:: molecule
	integer	::	atom_path(0:6)
	integer	::	bond_path(0:6)
	integer ::  i,count 
	integer ::  mark,mark1,mark2 
	integer ::  pi_path(0:4) 

	! check the number of pi electrons 
	mark=1
  count=0

	do i = 0,5 -1 
		 if(molecule%atom(atom_path(i)-1)%ttype.eq.'C.ar')  then
			 count = count +1
			 pi_path(i)=1
		 elseif(molecule%atom(atom_path(i)-1)%ttype.eq.'N.ar') then
			 count = count +1
			 pi_path(i)=1
		 elseif(molecule%atom(atom_path(i)-1)%ttype.eq.'C.2') then
			 count = count +1
			 pi_path(i)=1
		 elseif(molecule%atom(atom_path(i)-1)%ttype.eq.'N.2') then
			 count = count +1
			 pi_path(i)=1
		 elseif(molecule%atom(atom_path(i)-1)%ttype.eq.'N.pl3') then
			 count = count +2
			 pi_path(i)=2
		 elseif(molecule%atom(atom_path(i)-1)%ttype.eq.'N.am') then
			 count = count +2
			 pi_path(i)=2
		 elseif(molecule%atom(atom_path(i)-1)%ttype.eq.'O.3') then
			 count = count +2
			 pi_path(i)=2
		 elseif(molecule%atom(atom_path(i)-1)%ttype.eq.'S.3') then
			 count = count +2
			 pi_path(i)=2
		 else 
			 mark=0
			 exit
		 end if
	end do

	if((mark.eq.0).or.(count.ne.6)) then
		Molecule_Aromatic_Ring_Check_5 = 0
		return
	end if

	! check if there are two continuous single bonds on the path
	! but it is okay if these two single bonds explain the special atom
	mark=0 
	do i = 0,4 -1 
!		write(*,*) i, pi_path(i), '  ', adjustl(trim(molecule%bond(bond_path(i)-1)%ttype)),'  ',adjustl(trim(molecule%bond(bond_path(i)-1)%ttype))
		if((adjustl(trim(molecule%bond(bond_path(i)-1)%ttype)).eq.'1') &
		.or.(adjustl(trim(molecule%bond(bond_path(i)-1)%ttype)).eq.'am')) then
			mark1=1 
		else
			mark1=0
		end if

!		write(*,*) i+1, pi_path(i+1),'  ',adjustl(trim(molecule%bond(bond_path(i+1)-1)%ttype)),'  ',adjustl(trim(molecule%bond(bond_path(i+1)-1)%ttype))
		if((adjustl(trim(molecule%bond(bond_path(i+1)-1)%ttype)).eq.'1') &
		.or.(adjustl(trim(molecule%bond(bond_path(i+1)-1)%ttype)).eq.'am')) then
			mark2=1 
		else
			mark2=0 
		end if

!		write(*,*) 'mark1, mark2, pi_path = ', mark1, mark2, pi_path(i)
	
		if((mark1.ne.1).or.(mark2.ne.1)) then
			goto 1 
		elseif(pi_path(i+1).eq.2) then
			goto 1 
		else
!			write(*,*) 'not good: ', i, '  ',adjustl(trim(molecule%bond(bond_path(i)-1)%ttype)),'  ',adjustl(trim(molecule%bond(bond_path(i)-1)%ttype))
!			write(*,*) 'not good: ', i, '  ',adjustl(trim(molecule%bond(bond_path(i)-1)%ttype)),'  ',adjustl(trim(molecule%bond(bond_path(i)-1)%ttype))
			mark=1
			exit
		end if
1	end do

	if(mark.eq.1) then
		Molecule_Aromatic_Ring_Check_5 = 0
	else
		Molecule_Aromatic_Ring_Check_5 = 1
!		write(*,*) 'success'
	end if
end function Molecule_Aromatic_Ring_Check_5

subroutine Molecule_Detect_Rings(molecule)
! find the ring systems in the given molecule
! define all the atoms in a ring as atom%ring=1, otherwise atom%ring=0
	type(tXMolecule),intent(inout)	:: molecule
	integer													::  i,j,mark,count 
	integer,pointer									::  atom_path(:)
	integer,pointer									::  bond_path(:)
	
	! first of all, check whether there is any ring in the molecule
	if((molecule%num_bond-molecule%num_atom+molecule%num_subst).eq.0) then
		molecule%atom(:)%ring=0 
		molecule%bond(:)%ring=0 
		return 
	end if

	! detect terminal atoms
	do i = 0,molecule%num_atom -1  	
		if(molecule%atom(i)%num_neib<=1) then
			molecule%atom(i)%ring=0 
		else
			molecule%atom(i)%ring=-1 	! still ambiguous here 
		end if
	end do
	
	! collapse the structure to eliminate branches
	do		
		mark=0 
		do i = 0,molecule%num_atom -1 
			if(molecule%atom(i)%ring.ne.-1) goto 11
				
			count=0 	! count possible ring neighbors	
			do j = 0,molecule%atom(i)%num_neib -1 
				if(molecule%atom(molecule%atom(i)%neib(j)-1)%ring.eq.0) then
					goto 1 
				else
					count = count +1 
				end if
1			end do

			if(count<=1) then
				molecule%atom(i)%ring=0
				mark = mark +1
			end if
11	end do
		if(mark.eq.0) exit
	end do

	! detect branching bonds
	mark=0 
	do i = 0,molecule%num_bond -1  
		if(molecule%bond(i)%valid<=0) then
			goto 2  
		elseif(molecule%atom(molecule%bond(i)%atom_1-1)%ring.eq.0) then
			molecule%bond(i)%ring=0 
		elseif(molecule%atom(molecule%bond(i)%atom_2-1)%ring.eq.0) then
			molecule%bond(i)%ring=0 
		else
			molecule%bond(i)%ring=-1
			mark = mark +1	! ambiguous
		end if
2	end do

	if(mark.eq.0) return  	! finished, no ring detected

	! now comes the ring perception algorithm
	allocate(atom_path(0:molecule%num_bond),stat=err)
	if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'
	allocate(bond_path(0:molecule%num_bond),stat=err)
	if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'

	do i = 0,molecule%num_bond -1 	
		if(molecule%bond(i)%ring.ne.-1) goto 3 
		atom_path(:) = 0
		bond_path(:) = 0 

		! =============================== CHECK LAST ARGUMENT
		count = Molecule_Look_For_A_Ring(molecule, molecule%bond(i)%id,atom_path,bond_path,0) 

		if(count.eq.0) then
			molecule%bond(i)%ring=0
			goto 3
		end if

		do j = 0,count -1  
			 molecule%bond(bond_path(j)-1)%ring=1 
			 molecule%atom(atom_path(j)-1)%ring=1 
		end do
3	end do

	! define all the left atoms which are not in ring
	do i = 0,molecule%num_atom -1 
		if(molecule%atom(i)%ring.ne.-1) then
			goto 4 
		else 
			molecule%atom(i)%ring=0 
		end if

		if(associated(atom_path)) deallocate(atom_path)
		if(associated(bond_path)) deallocate(bond_path)
4	end do
end subroutine Molecule_Detect_Rings

integer function Molecule_Look_For_A_Ring(molecule, bond_id, atom_path, bond_path, required_size)
! this is a bond-based algorithm to detect rings in a molecule%
! it returns the size of the ring (if zero, the given bond is not in a ring)
! also the corresponding atom path and bond path
	type(tXMolecule),intent(inout)	:: molecule
	integer :: bond_id
	integer,pointer :: atom_path(:)
	integer,pointer :: bond_path(:)
	integer :: required_size

	integer ::  i,ring_size=0 
	integer ::  p,next 
	integer,allocatable ::  choice(:,:)

	! initialize the variables first
	allocate(choice(0:molecule%num_bond-1,0:MAX_BOND_NEIB-1),stat=err)
	if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'
	choice(:,:) = 0

	! note that atom_path() and bond_path() is not cleaned here
	! they are supposed to be cleaned at the upper level 
	do i = 0,molecule%num_bond -1  
		 if(molecule%bond(i)%valid<=0) goto 1 
		 call Molecule_Reset_Choices(molecule,molecule%bond(i)%id,choice(i,:)) 
1	end do

	! put the given bond at the beginning of the searching chain 
	p=0  
	bond_path(0)=bond_id
	atom_path(0)=molecule%bond(bond_id-1)%atom_1		! bond_id -1

	do
		! search the next possible step 
		do
			next = Molecule_Get_Choices(molecule, choice(bond_path(p)-1,:),atom_path(p),p+1,bond_path) 
			if(next.eq.0) then 	! next step is not available 
				if(p.eq.0) then	! the given bond is not in a ring%
					ring_size=0 
					goto 10
				end if

				! trace back to the previous step 
				call Molecule_Reset_Choices(molecule,bond_path(p),choice(bond_path(p)-1,:)) 
				bond_path(p)=0 
				atom_path(p)=0 
				p = p -1 
				goto 2 
			else
				exit
			end if
2		end do

		! define the next step 

		p = p +1 
		bond_path(p)=next 
		atom_path(p)=Molecule_Two_Bonds_Connection_Check(molecule%bond(bond_path(p-1)-1), molecule%bond(bond_path(p)-1))

		! record the searching history, prevent repeating  
		call Molecule_Clean_Choices(molecule,bond_path(p),choice(bond_path(p-1)-1,:)) 

		! now check if a ring has formed
		if(Molecule_Two_Bonds_Connection_Check(molecule%bond(bond_path(0)-1),molecule%bond(bond_path(p)-1)).eq.atom_path(0)) then
			if(required_size<=0) then ! a ring is found
				ring_size=p+1
				goto 10
      elseif((p+1).eq.required_size) then  ! a ring with required size
				ring_size=p+1
				goto 10	 
			else   ! pre-mature ring, trace back
				call Molecule_Reset_Choices(molecule,bond_path(p), choice(bond_path(p)-1,:)) 
				bond_path(p)=0 
				atom_path(p)=0 
				p = p -1 
				goto 3 
			end if
		else   ! no ring is formed at this step
			if(required_size<=0) then
				goto 3 
			elseif((p+1).eq.required_size) then ! no ring on this path, trace back
				call Molecule_Reset_Choices(molecule,bond_path(p), choice(bond_path(p)-1,:)) 
				bond_path(p)=0 
				atom_path(p)=0 
				p = p -1 
				goto 3 
			else
				goto 3
			end if
		end if
3	end do

10	deallocate(choice) 
	Molecule_Look_For_A_Ring = ring_size 
end function Molecule_Look_For_A_Ring

subroutine Molecule_Reset_Choices(molecule, bond_id,choice)
! find all the neighboring bonds for the given bond 
! it is exclusively used by the ring perception algorithm
	type(tXMolecule) :: molecule
	integer	:: bond_id
	integer	:: choice(0:MAX_BOND_NEIB -1)

	integer ::  i,tmp 

	choice(:) = 0
	
	do i = 0,molecule%bond(bond_id -1)%num_neib -1		! bond_id -1
		tmp=molecule%bond(bond_id-1)%neib(i) 
		if(molecule%bond(tmp-1)%ring.eq.0) then
			choice(i)=0
		elseif(tmp>molecule%num_bond) then
			choice(i)=0 
		else
			choice(i)=tmp  
		end if
	end do
end subroutine Molecule_Reset_Choices

integer function Molecule_Get_Choices(molecule, choice,wrong_end,current_path_length,current_path)
! notice that wrong_end helps to define the direction of the searching path 
	type(tXMolecule)	:: molecule
	integer	:: choice(0:MAX_BOND_NEIB -1)
	integer	:: wrong_end,current_path_length
	integer,pointer	:: current_path(:)

	integer ::  i,j,tresult 
	integer ::  mark 

	tresult=0 

	do i = 0,MAX_BOND_NEIB -1 
		if(choice(i).eq.0) then
			goto 1 
		elseif(choice(i)>molecule%num_bond) then
			goto 1 
		elseif(molecule%bond(choice(i)-1)%ring.eq.0) then
			goto 1
		elseif(molecule%bond(choice(i)-1)%atom_1.eq.wrong_end) then
			goto 1 
    elseif(molecule%bond(choice(i)-1)%atom_2.eq.wrong_end) then
			goto 1
		end if
	
		! Is this bond already included in the current path?
		mark=0 
		do j = 0,current_path_length -1 
			if(current_path(j).ne.choice(i)) then
				goto 2 
			else
				mark=1
				exit
			end if
2		end do
		
		if(mark.eq.1) then
			goto 1 
		else
			tresult = choice(i)
			exit
		end if
1	end do

	Molecule_Get_Choices =  tresult 
end function Molecule_Get_Choices

subroutine Molecule_Clean_Choices(molecule,tried_choice,choice)
	type(tXMOlecule)	:: molecule
	integer :: tried_choice
	integer :: choice(0:MAX_BOND_NEIB-1)
	integer ::  i 

		do i = 0,MAX_BOND_NEIB -1 
		if(choice(i).eq.0) then
			goto 1 
		elseif(choice(i).eq.tried_choice) then
			choice(i)=0 
		else
			goto 1 
		end if
1	end do

end subroutine Molecule_Clean_Choices

integer function Molecule_Get_Num_HB_Atom(molecule)
	type(tXMolecule)	:: molecule
	integer ::  i,num=0 

	do i = 0,molecule%num_atom -1 
		if(molecule%atom(i)%valid.ne.0) then
			if(molecule%atom(i)%hb.eq.'D') then
				num = num +1 
			elseif(molecule%atom(i)%hb.eq.'DA') then
				num = num +1 
			elseif(molecule%atom(i)%hb.eq.'A') then
				num = num +1 
			end if
		end if
	end do

	Molecule_Get_Num_HB_Atom = num 
end function Molecule_Get_Num_HB_Atom

real(kind=prec) function Molecule_Count_Rotor(molecule)
! if a single bond is normal, bond%valid=1 
! if a single bond is a rotor, bond%valid=2 
	type(tXMolecule)	::	molecule
	integer						::  i,mark,tmp 
	integer						::  id1,id2 
	real(kind=prec)							::  sum 

	! clean the variables
	molecule%atom(:)%score=0.000 

	! eliminate all the bonds in_ring and all the non-single bonds
	do i = 0,molecule%num_bond -1 
		 if(molecule%bond(i)%ring<=0) then
			if(molecule%bond(i)%ttype.ne.'1') then
				molecule%bond(i)%valid=2 
			end if
		end if
	end do

	! eliminate all the R-H, R-X, R-OH, R-NH2, R-CH3 bonds 
	do i = 0,molecule%num_bond -1 
		 if(molecule%bond(i)%valid.ne.2) goto 2

		 id1=molecule%bond(i)%atom_1
		 id2=molecule%bond(i)%atom_2

		 if(molecule%atom(id1-1)%num_nonh.eq.1.or.molecule%atom(id2-1)%num_nonh.eq.1) then
			molecule%bond(i)%valid=1 
		 else
			goto 2
		 end if
2	end do

	! sp2-sp2 rotors
	do i = 0,molecule%num_bond -1 
		 if(molecule%bond(i)%valid.ne.2) goto 3

		 id1=molecule%bond(i)%atom_1
		 id2=molecule%bond(i)%atom_2

		 mark=0  
		 tmp = Molecule_Get_Atom_Hybridizing_Type(molecule,molecule%atom(id1-1)%ttype) 
		 if(tmp.eq.1.or.tmp.eq.2) mark = mark +1 

		 tmp=Molecule_Get_Atom_Hybridizing_Type(molecule,molecule%atom(id2-1)%ttype) 
     if(tmp.eq.1.or.tmp.eq.2) mark = mark +1 
		 if(mark.eq.2) molecule%bond(i)%valid=1
3	end do

	! eliminate terminal rotors, e.g. -PO3, -CF3, -CMe3, -NMe3
	do i = 0,molecule%num_bond -1 
     if(molecule%bond(i)%valid.ne.2) goto 4

		 id1=molecule%bond(i)%atom_1
		 id2=molecule%bond(i)%atom_2

		 if(Molecule_Judge_Terminal_Atom(molecule,molecule%atom(id1-1),id2).ne.0) then
			 molecule%bond(i)%valid=1  
		 elseif(Molecule_Judge_Terminal_Atom(molecule,molecule%atom(id2-1),id1).ne.0) then
			 molecule%bond(i)%valid=1 
		 end if
4	 end do

	! eliminate abnormal rotors
	do i = 0,molecule%num_bond -1 
		if(molecule%bond(i)%valid.ne.2) goto 5

		id1=molecule%bond(i)%atom_1
	  id2=molecule%bond(i)%atom_2

		if(molecule%atom(id1-1)%valid<=0) then
			molecule%bond(i)%valid=1 
		elseif(molecule%atom(id2-1)%valid<=0) then
			molecule%bond(i)%valid=1 
		end if
5	end do

	! now count the frozen rotors, all the rotors have been labeled as 2
	sum=0.000 
	do i = 0,molecule%num_bond -1 
		if(molecule%bond(i)%valid.ne.2) then
			goto 6
		else 
			sum = sum + 1.00 
		end if
6	end do
	
	Molecule_Count_Rotor = sum
end function Molecule_Count_Rotor

integer function Molecule_Judge_Terminal_Atom(molecule,atm, partner)
	type(tXMolecule) :: molecule
	type(tXAtom)		 :: atm
	integer					 ::	partner
	integer					 ::  i,j 

	if(Molecule_Get_Atom_Hybridizing_Type(molecule,atm%ttype).ne.3) then
		Molecule_Judge_Terminal_Atom = 0
		return
	end if

	if(atm%num_nonh.ne.4) then
		Molecule_Judge_Terminal_Atom = 0
		return
	end if

	do i = 0,atm%num_nonh -1 
		if(atm%neib(i).eq.partner) then
			goto 1
		elseif(molecule%atom(atm%neib(i)-1)%num_nonh.ne.1) then
			Molecule_Judge_Terminal_Atom = 0
			return
		end if
1	end do

	do i = 0,atm%num_nonh-1 -1 
		if(atm%neib(i).eq.partner) goto 2

		do j = i,atm%num_nonh -1 
			 if(atm%neib(j).eq.partner) then
				goto 3
			 elseif(molecule%atom(atm%neib(j)-1)%xtype.eq.molecule%atom(atm%neib(i)-1)%xtype) then
				goto 3
			 else 
				Molecule_Judge_Terminal_Atom = 0
			 end if
3		end do
2	end do

	Molecule_Judge_Terminal_Atom = 1
end function Molecule_Judge_Terminal_Atom

integer function Molecule_Get_Atom_Hybridizing_Type(molecule,ttype)	
	type(tXMolecule)	:: molecule
	character(len=10) :: ttype
	integer ::  mark  	! sp=1, sp2=2  sp3=3, none=4

	if(ttype.eq.'C.3') then  
		mark=3 
	elseif(ttype.eq.'C.2') then    
		mark=2 
	elseif(ttype.eq.'C.1') then    
		mark=1  
	elseif(ttype.eq.'C.cat') then  
		mark=2 
	elseif(ttype.eq.'C.ar') then   
		mark=2 
	elseif(ttype.eq.'H') then      
		mark=4 
	elseif(ttype.eq.'H.spc') then  
		mark=4 
	elseif(ttype.eq.'N.4') then    
		mark=3 
	elseif(ttype.eq.'N.3') then    
		mark=3
	elseif(ttype.eq.'N.2') then    
		mark=2 
	elseif(ttype.eq.'N.1') then    
		mark=1  
	elseif(ttype.eq.'N.ar') then   
		mark=2 
	elseif(ttype.eq.'N.pl3') then	
		mark=2 
	elseif(ttype.eq.'N.am') then   
		mark=2 
	elseif(ttype.eq.'O.3') then    
		mark=3 
	elseif(ttype.eq.'O.2') then    
		mark=2 
	elseif(ttype.eq.'O.co2') then  
		mark=2 
	elseif(ttype.eq.'P.3') then    
		mark=3 
	elseif(ttype.eq.'S.3') then    
		mark=3 
	elseif(ttype.eq.'S.2') then    
		mark=2 
	elseif(ttype.eq.'S.o') then    
		mark=3 
	elseif(ttype.eq.'S.o2') then   
		mark=3 
	elseif(ttype.eq.'F') then      
		mark=4 
	elseif(ttype.eq.'Cl') then     
		mark=4 
	elseif(ttype.eq.'Br') then     
		mark=4 
	elseif(ttype.eq.'I') then      
		mark=4 
	elseif(ttype.eq.'Si') then     
		mark=3 
	else 
		mark=4 
	end if

	Molecule_Get_Atom_Hybridizing_Type = mark 
end function Molecule_Get_Atom_Hybridizing_Type

real(kind=prec) function Molecule_Calculate_LogP(molecule)
	type(tXMolecule)	:: molecule
	integer						:: i,flag 
	real(kind=prec)							:: xlogp 
	character(len=20) :: ttype
	
	! coordinates are not used in this routine

	! logp_factor(0:9) is declared in calc_xscore
	logp_factor(0)%symbol = 'Hydrophobic carbon' 
	logp_factor(1)%symbol = 'Internal H-bond' 
	logp_factor(2)%symbol = 'Halogen 1-3 pair' 
	logp_factor(3)%symbol = 'Aromatic nitrogen 1-4 pair' 
	logp_factor(4)%symbol = 'Ortho sp3 oxygen pair' 
	logp_factor(5)%symbol = 'Acceptor 1-5 pair' 
	logp_factor(6)%symbol = 'Paralleled donor pair' 
	logp_factor(7)%symbol = 'Alpha amino acid' 
	logp_factor(8)%symbol = 'Salicylic acid' 
	logp_factor(9)%symbol = 'P-amino sulfonic acid' 

	logp_factor(0)%coeff = 0.211 
	logp_factor(1)%coeff = 0.429 
	logp_factor(2)%coeff = 0.137 
	logp_factor(3)%coeff = 0.485 
	logp_factor(4)%coeff =-0.268 
	logp_factor(5)%coeff = 0.580 
	logp_factor(6)%coeff =-0.423 
	logp_factor(7)%coeff =-2.166 
	logp_factor(8)%coeff = 0.554 
	logp_factor(9)%coeff =-0.501 

	xlogp=0.000 

	do i = 0,molecule%num_atom -1 
		if(molecule%atom(i)%valid.ne.0) then
			call Molecule_Get_XLOGP_Type(molecule,molecule%atom(i)%id,ttype)
			molecule%atom(i)%logp = FF_Get_Atom_LogP(ttype) 
      xlogp= xlogp + molecule%atom(i)%logp 
		end if
	end do

  logp_factor(0:10 -1)%num=0 

	! this flag determines if correction factors are also distributed
	! onto each atoms or not
	flag=0 
	logp_factor(0)%num=Molecule_Count_Hydrophobic_Carbon(molecule,flag) ! no coords used
	logp_factor(1)%num=Molecule_Count_Internal_HBond(molecule,flag)			! no coords used
	logp_factor(2)%num=Molecule_Count_Halogen_1_3_Pair(molecule,flag)		! no coords used
	logp_factor(3)%num=Molecule_Count_Nar_1_4_Pair(molecule,flag)				! no coords used
	logp_factor(4)%num=Molecule_Count_O3_1_4_Pair(molecule,flag)				! no coords used
	logp_factor(5)%num=Molecule_Count_Acceptor_1_5_Pair(molecule,flag)	! no coords used
	logp_factor(6)%num=0.000																						! no coords used
	logp_factor(7)%num=Molecule_Count_Amino_Acid(molecule,flag)					! no coords used
	logp_factor(8)%num=Molecule_Count_Salicylic_Acid(molecule,flag)			! no coords used
	logp_factor(9)%num=Molecule_Count_Sulfonic_Acid(molecule,flag)			! no coords used

	do i = 0,10 -1  
		xlogp=xlogp + (logp_factor(i)%num*logp_factor(i)%coeff) 
	end do

	Molecule_Calculate_LogP = xlogp 
end function Molecule_Calculate_LogP

real(kind=prec) function Molecule_Count_Sulfonic_Acid(molecule,flag)
	type(tXMolecule)	:: molecule
	integer						:: flag
	integer ::  i,j,num,tmp1,tmp2 

	num=0 
	do i = 0,molecule%num_atom -1 
		if((molecule%atom(i)%ttype.ne.'S.o2')) goto 1  
		if(molecule%atom(i)%ring.ne.0) goto 1 

		tmp1=Molecule_Adjacent_Aromatic_Check(molecule,(molecule%atom(i)%id)) 
		if(tmp1.eq.0) goto 1 

		do j = 0,molecule%num_atom -1 
			if((molecule%atom(j)%xtype.ne.'N.3.h2.pi=1')) goto 2 

			tmp2=Molecule_Adjacent_Aromatic_Check(molecule,(molecule%atom(j)%id)) 
			if(tmp2.eq.0) goto 2 

			if(Molecule_Connection_1_6_Check(molecule,molecule%atom(i)%id,molecule%atom(j)%id).eq.0) goto 2 
			if(Molecule_Connection_1_4_Check(molecule,tmp1,tmp2).eq.0) goto 2 

			num = num +1  

			if(flag.ne.0) then
				molecule%atom(i)%logp=molecule%atom(i)%logp+(LOGP_SULFONIC_ACID/2.0) 
				molecule%atom(j)%logp=molecule%atom(j)%logp+(LOGP_SULFONIC_ACID/2.0) 
			end if

			exit 
2		end do
1	end do

	Molecule_Count_Sulfonic_Acid = num 
end function Molecule_Count_Sulfonic_Acid

real(kind=prec) function Molecule_Count_Salicylic_Acid(molecule,flag)
	type(tXMolecule)	:: molecule
	integer						:: flag
	integer ::  i,j,num,tmp,mark 

	num=0 

	do i = 0,molecule%num_atom -1 
		if((molecule%atom(i)%ttype.ne.'O.2')) goto 1 

		tmp=molecule%atom(i)%neib(0)-1 
		if(molecule%atom(tmp)%ttype(1:1).ne.'C') goto 1 
		if(molecule%atom(tmp)%ring.ne.0) goto 1 
		if(Molecule_Adjacent_Aromatic_Check(molecule,(molecule%atom(tmp)%id)).eq.0) goto 1 
		
		mark=0 

		do j = 0,molecule%num_atom -1 
			if(i.eq.j) goto 2 
			if((molecule%atom(j)%ttype.ne.'O.3')) goto 2 
			if(molecule%atom(j)%ring.ne.0) goto 2 
			if(Molecule_Connection_1_3_Check(molecule,(molecule%atom(i)%id), (molecule%atom(j)%id)).eq.0) goto 2 
			mark=0
			exit
2		end do
	

		if(mark.eq.0) goto 1 

		mark=0 
		do j = 0,molecule%num_atom -1 
			if(i.eq.j) goto 3 
			if((molecule%atom(j)%ttype.ne.'O.3')) goto 3 
			if((molecule%atom(j)%hb.ne.'DA')) goto 3 
			if(Molecule_Adjacent_Aromatic_Check(molecule,(molecule%atom(j)%id)).eq.0) goto 3 
			if(Molecule_Connection_1_5_Check(molecule,(molecule%atom(i)%id), (molecule%atom(j)%id)).eq.0) goto 3 

			num = num +1 
			if(flag.ne.0) then
				molecule%atom(i)%logp=molecule%atom(i)%logp+((LOGP_SALICYLIC_ACID)/2.0) 
				molecule%atom(j)%logp=molecule%atom(j)%logp+((LOGP_SALICYLIC_ACID)/2.0) 
			end if
			exit 
3		end do
		
		if(num.ne.0) exit 
1	end do

	Molecule_Count_Salicylic_Acid =  num 
end function Molecule_Count_Salicylic_Acid

real(kind=prec) function Molecule_Count_Amino_Acid(molecule,flag)
	type(tXMolecule)	:: molecule
	integer						:: flag
  integer ::  i,j,tmp,num,mark 

	num=0 

	do i = 0,molecule%num_atom -1 
		if((molecule%atom(i)%ttype.ne.'O2')) goto 1 

		tmp=molecule%atom(i)%neib(0)-1 
		if(molecule%atom(tmp)%ttype(1:1).ne.'C') goto 1 
		if(molecule%atom(tmp)%ring.ne.0) goto 1 

		mark=0 

		do j = 0,molecule%num_atom -1 
			if(i.eq.j) goto 2
			if((molecule%atom(j)%ttype.ne.'O.3')) goto 2 
			if((molecule%atom(j)%hb.ne.'DA')) goto 2 
			if(molecule%atom(j)%ring.ne.0) goto 2 
			if(Molecule_Connection_1_3_Check(molecule,(molecule%atom(i)%id), (molecule%atom(j)%id)).eq.0) goto 2 
			
			mark=1
			exit
2		end do

		if(mark.eq.0) goto 1 

		do j = 0,molecule%num_atom -1 
			if((molecule%atom(j)%xtype.ne.'N.3.h2.pi=0')) goto 3 
			if(Molecule_Connection_1_4_Check(molecule,(molecule%atom(i)%id), (molecule%atom(j)%id)).eq.0) goto 3 

			num = num +1 
			if(flag.ne.0) then
				molecule%atom(i)%logp=molecule%atom(i)%logp+((LOGP_AMINO_ACID)/2.0) 
				molecule%atom(j)%logp=molecule%atom(j)%logp+((LOGP_AMINO_ACID)/2.0)		
			end if

			exit 
3		end do

		if(num.ne.0) exit 

		do j = 0,molecule%num_atom -1 
			if((molecule%atom(j)%ttype.ne.'N.ar')) goto 4 
			if(Molecule_Connection_1_4_Check(molecule,(molecule%atom(i)%id), (molecule%atom(j)%id)).eq.0) goto 4 

			num = num +1 

			if(flag.ne.0) then
				molecule%atom(i)%logp=molecule%atom(i)%logp+((LOGP_AMINO_ACID)/2.0) 
				molecule%atom(j)%logp=molecule%atom(j)%logp+((LOGP_AMINO_ACID)/2.0) 
			end if

			exit 
4		end do

		if(num.ne.0) exit 
1	end do	

	Molecule_Count_Amino_Acid = num 
end function Molecule_Count_Amino_Acid

real(kind=prec) function Molecule_Count_Acceptor_1_5_Pair(molecule,flag)
	type(tXMolecule)	:: molecule
	integer						:: flag
  integer ::  i,j,num,tmp1,tmp2 

  num=0 

	do i = 0,molecule%num_atom-1 -1 
		if((molecule%atom(i)%hb.ne.'A')) goto 1 
		if((molecule%atom(i)%ttype.eq.'O.3')) goto 1 
		if((molecule%atom(i)%ttype.eq.'N.2')) goto 1 
		if((molecule%atom(i)%ttype.eq.'N.ar')) goto 1 

		tmp1=molecule%atom(i)%neib(0)-1 
		if(molecule%atom(tmp1)%ttype(1:1).eq.'S') goto 1 
		if(molecule%atom(tmp1)%ttype(1:1).eq.'P') goto 1 
		if(molecule%atom(tmp1)%ring.ne.0) goto 1 

		do j = i,molecule%num_atom -1 
			if((molecule%atom(j)%hb.ne.'A')) goto 2 
			if((molecule%atom(j)%ttype.eq.'O.')) goto 2 
			if((molecule%atom(j)%ttype.eq.'N.2')) goto 2 
			if((molecule%atom(j)%ttype.eq.'N.ar')) goto 2 

			tmp2=molecule%atom(j)%neib(0)-1 
			if(molecule%atom(tmp2)%ttype(1:1).eq.'S') goto 2 
			if(molecule%atom(tmp2)%ttype(1:1).eq.'P') goto 2 
			if(molecule%atom(tmp2)%ring.ne.0) goto 2 

			if(Molecule_Connection_1_5_Check(molecule,(molecule%atom(i)%id), (molecule%atom(j)%id)).eq.0) goto 2 

			num = num +1 

			if(flag.ne.0) then
				molecule%atom(i)%logp=molecule%atom(i)%logp+((LOGP_ACCEPTOR_PAIR)/2.0) 
				molecule%atom(j)%logp=molecule%atom(j)%logp+((LOGP_ACCEPTOR_PAIR)/2.0) 
			else 
				continue 
			end if
2   end do
1 end do      

  Molecule_Count_Acceptor_1_5_Pair = num 
end function Molecule_Count_Acceptor_1_5_Pair

real(kind=prec) function Molecule_Count_O3_1_4_Pair(molecule,flag)
	type(tXMolecule) :: molecule
	integer					 :: flag
	integer ::  i,j 
	real(kind=prec) ::  num 
	integer ,allocatable::  trecord(:) 

	allocate(trecord(0:molecule%num_atom -1),stat=err)
	if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'

	trecord(:)=0 

	num=0 

	do i = 0,molecule%num_atom-1 -1 
		if((molecule%atom(i)%ttype.ne.'O.3')) goto 1 
		if((molecule%atom(i)%hb.eq.'DA')) goto 1 
		if(molecule%atom(i)%ring.ne.0) goto 1 
		if(Molecule_Adjacent_Aromatic_Check(molecule,(molecule%atom(i)%id)).eq.0) goto 1 

		do j = i,molecule%num_atom -1 
			if((molecule%atom(j)%ttype.ne.'O.3')) goto 2 
			if((molecule%atom(j)%hb.eq.'DA')) goto 2 
			if(molecule%atom(j)%ring.ne.0) goto 2 
      if(Molecule_Adjacent_Aromatic_Check(molecule,(molecule%atom(j)%id)).eq.0) goto 2 
			if(Molecule_Connection_1_4_Check(molecule,(molecule%atom(i)%id),(molecule%atom(j)%id)).eq.0) goto 2 

			trecord(i) = trecord(i) +1
			trecord(j) = trecord(j) +1 
2		end do
1	end do

	do i = 0,molecule%num_atom -1 
		if(trecord(i).eq.0) goto 3 

		num=num+0.500 

		if(flag.ne.0) then
			molecule%atom(i)%logp=molecule%atom(i)%logp+(LOGP_O3_PAIR/2.0) 
		else 
			goto 3 
		end if
3	end do

	deallocate(trecord) 

	Molecule_Count_O3_1_4_Pair = num 
end function Molecule_Count_O3_1_4_Pair

integer function Molecule_Adjacent_Aromatic_Check(molecule,id_in)
	type(tXMolecule)	:: molecule
	integer						::id, id_in
  integer ::  i,num,tmp,mark 

	id = id_in
	id = id -1 

	mark=0
	num=molecule%atom(id)%num_neib

	do i = 0,num -1 
		tmp=molecule%atom(id)%neib(i) 
		if((molecule%atom(tmp-1)%ttype.eq.'C.ar')) then
			mark=tmp
			exit
		end if
	end do
  
	Molecule_Adjacent_Aromatic_Check = mark 
end function Molecule_Adjacent_Aromatic_Check

real(kind=prec) function Molecule_Count_Nar_1_4_Pair(molecule,flag)
	type(tXMolecule) :: molecule
	integer					 :: flag
	integer ::  i,j,num,tmp1,tmp2,tmp3,tmp4 

	num=0 

	do i = 0,molecule%num_atom-1 -1 
		if((molecule%atom(i)%ttype.ne.'N.ar')) goto 1 

		tmp1=molecule%atom(i)%neib(0)
		tmp2=molecule%atom(i)%neib(1)

		do j = i,molecule%num_atom -1 
			if((molecule%atom(j)%ttype.ne.'N.ar')) then
				goto 2 
			elseif(Molecule_Connection_1_4_Check(molecule,(molecule%atom(i)%id), (molecule%atom(j)%id)).eq.0) then
				goto 2 
			else
				tmp3=molecule%atom(j)%neib(0)
				tmp4=molecule%atom(j)%neib(1)

				if(Molecule_Connection_1_2_Check(molecule,(tmp1),(tmp3)).eq.1) then
					if(Molecule_Connection_1_2_Check(molecule,(tmp2),(tmp4)).eq.1) then
						num = num +1 
						if(flag.ne.0) then
							molecule%atom(i)%logp=molecule%atom(i)%logp+((LOGP_NAR_PAIR)/2.0) 
							molecule%atom(j)%logp=molecule%atom(j)%logp+((LOGP_NAR_PAIR)/2.0) 
						end if
					else 
						goto 2  
					end if
				elseif(Molecule_Connection_1_2_Check(molecule,tmp1,tmp4).eq.1) then
					if(Molecule_Connection_1_2_Check(molecule,tmp2,tmp3).eq.1) then
						num = num +1 
						if(flag.ne.0) then
							molecule%atom(i)%logp=molecule%atom(i)%logp+((LOGP_NAR_PAIR)/2.0) 
							molecule%atom(j)%logp=molecule%atom(j)%logp+((LOGP_NAR_PAIR)/2.0) 
						end if
					end if
				else 
					goto 2  
				end if
			end if
2		end do
1	end do

	Molecule_Count_Nar_1_4_Pair = num 
end function Molecule_Count_Nar_1_4_Pair

real(kind=prec) function Molecule_Count_Halogen_1_3_Pair(molecule,flag)
	type(tXMolecule)	:: molecule
	integer						:: flag
	integer ::  i,j 
	integer ::  num1,num2 

	num1=0
	num2=0 

	do i = 0,molecule%num_atom-1 -1 
		if((molecule%atom(i)%ttype.ne.'F')) goto 1 

		do j = i,molecule%num_atom -1 
			if((molecule%atom(j)%ttype.ne.'F')) then
				goto 2 
			elseif(Molecule_Connection_1_3_Check(molecule,molecule%atom(i)%id,molecule%atom(j)%id).eq.0) then
				goto 2 
			end if

			num1 = num1 +1 

			if(flag.ne.0) then
				molecule%atom(i)%logp=molecule%atom(i)%logp+((LOGP_HALOGEN_PAIR)/2.0) 
				molecule%atom(j)%logp=molecule%atom(j)%logp+((LOGP_HALOGEN_PAIR)/2.0) 
			else 
				goto 2 
			end if
2		end do
1	end do

	do i = 0,molecule%num_atom-1 -1 
		if((molecule%atom(i)%ttype.ne.'Cl').and.(molecule%atom(i)%ttype.ne.'Br').and.(molecule%atom(i)%ttype.ne.'I')) goto 3 

		do j = i,molecule%num_atom -1 
			if((molecule%atom(j)%ttype.ne.'Cl').and.(molecule%atom(j)%ttype.ne.'Br').and.(molecule%atom(j)%ttype.ne.'I')) then
				goto 4 
			elseif(Molecule_Connection_1_3_Check(molecule,molecule%atom(i)%id,molecule%atom(j)%id).eq.0) then
				goto 4 
			end if

			num2 = num2 +1 

			if(flag.ne.0) then
				molecule%atom(i)%logp=molecule%atom(i)%logp+((LOGP_HALOGEN_PAIR)/2.0) 
				molecule%atom(j)%logp=molecule%atom(j)%logp+((LOGP_HALOGEN_PAIR)/2.0) 
			else 
				goto 4
			end if
4		end do
3	end do

	Molecule_Count_Halogen_1_3_Pair = num1+num2
end function Molecule_Count_Halogen_1_3_Pair

integer function Molecule_Adjacent_Ring_Check(molecule,id_in)
	type(tXMolecule) :: molecule
  integer					:: id, id_in
	integer ::  i,num,tmp,mark 
	
	id = id_in
	id = id -1 

  mark=0
	num=molecule%atom(id)%num_neib

	do i = 0,num -1 
		tmp=molecule%atom(id)%neib(i) 
		if(molecule%atom(tmp-1)%ring.ne.0) then
			mark=tmp
			exit
		end if
	end do

  Molecule_Adjacent_Ring_Check = mark 
end function Molecule_Adjacent_Ring_Check

real(kind=prec) function Molecule_Count_Internal_HBond(molecule,flag)
	type(tXMolecule) :: molecule
	integer	:: flag
	integer ::  i,j,mark1,mark2 
	real(kind=prec) ::  num 
	integer,allocatable ::  trecord(:)

	allocate(trecord(0:molecule%num_atom -1),stat=err)
	if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'

	trecord(:)=0 

	num=0 

	do i = 0,molecule%num_atom -1 
		if((molecule%atom(i)%hb.ne.'D').and.(molecule%atom(i)%hb.ne.'DA')) then
			goto 1
		elseif(molecule%atom(i)%ring.ne.0) then
			goto 1  ! not allowed in ring
		end if

		if(Molecule_Adjacent_Ring_Check(molecule,molecule%atom(i)%id).eq.0) then
			mark1=0 
		else 
			mark1=1 
		end if

		do j = 0,molecule%num_atom -1 
			if(i.eq.j) then
				goto 2  
			elseif((molecule%atom(j)%hb.ne.'A').and.(molecule%atom(j)%hb.ne.'DA')) then
				goto 2 
			elseif((molecule%atom(j)%ttype.eq.'O.3').and.(molecule%atom(j)%hb.eq.'A')) then
				goto 2 
			elseif((molecule%atom(j)%ttype.eq.'N.2')) then
				goto 2 
			elseif((molecule%atom(j)%ttype.eq.'N.ar')) then	
				goto 2 
			elseif(molecule%atom(j)%ring.ne.0) then	
				goto 2  ! not in ring 
			end if

			if(Molecule_Adjacent_Ring_Check(molecule, molecule%atom(j)%id).eq.0) then
				mark2=0 
			else 
				mark2=1 
			end if

			if(mark1.eq.1.and.mark2==1) then
				if(Molecule_Connection_1_4_Check(molecule,molecule%atom(i)%id,molecule%atom(j)%id).eq.0) then
					goto 2
				else 
					trecord(i) = trecord(i) +1
					trecord(j) = trecord(j) +1
				end if
			elseif(mark1.eq.1.and.mark2==0) then
				if(Molecule_Connection_1_5_Check(molecule,molecule%atom(i)%id,molecule%atom(j)%id).eq.0) then
					goto 2
				else
          trecord(i) =  +1
				end if
				trecord(j) = trecord(j) +1
			elseif(mark2.eq.0.and.mark2==1) then
				if(Molecule_Connection_1_5_Check(molecule,molecule%atom(i)%id,molecule%atom(j)%id).eq.0) then
					goto 2
				else
          trecord(i) = trecord(i) +1
					trecord(j) = trecord(j) +1
				end if
			end if
2		end do
1	end do

	do i = 0,molecule%num_atom -1 
		if(trecord(i).eq.0) then
			goto 3 
		end if

		num=num+0.500 

		if(flag.ne.0) then
			molecule%atom(i)%logp=molecule%atom(i)%logp+(LOGP_INTERNAL_HBOND/2.0) 
		else 
			goto 3 
		end if
3	end do

	deallocate(trecord) 

	Molecule_Count_Internal_HBond = num 
end function Molecule_Count_Internal_HBond

real(kind=prec) function Molecule_Count_Hydrophobic_Carbon(molecule,flag)
	type(tXMolecule) :: molecule
	integer	:: flag
	integer :: i,num 

  num=0 

	do i = 0,molecule%num_atom -1 
		if((molecule%atom(i)%ttype.ne.'C.3').and.(molecule%atom(i)%ttype.ne.'C.2'))	goto 1
		if(Molecule_Hydrophobic_Neighbor_Check(molecule,molecule%atom(i)%id).eq.0) goto 1

		num = num +1  
		if(flag.ne.0) then
			molecule%atom(i)%logp=molecule%atom(i)%logp+(LOGP_HYDROPHOBIC_CARBON) 
		end if
1	end do
 
	if(num>=10) num = num/2 

	Molecule_Count_Hydrophobic_Carbon = num 
end function Molecule_Count_Hydrophobic_Carbon

integer function Molecule_Hydrophobic_Neighbor_Check(molecule,id_in)
	type(tXMolecule) :: molecule
	integer					 :: id,id_in
  integer					 ::  i,mark 

	id = id_in
	id = id -1 
	mark=1 

	molecule%atom(2)%id = 3				! most cryptic bug ever...

	do i = 0,molecule%num_atom -1 
		if(i.eq.id) then
			goto 1
		elseif((molecule%atom(i)%ttype.eq.'F').or.   &
					 (molecule%atom(i)%ttype.eq.'Cl').or.  &
					 (molecule%atom(i)%ttype.eq.'Br').or.  &
					 (molecule%atom(i)%ttype.eq.'I').or.   &
					 (molecule%atom(i)%ttype.eq.'Si').or.  & 
						molecule%atom(i)%ttype(1:1).eq.'N'.or. &
						molecule%atom(i)%ttype(1:1).eq.'O'.or. &
						molecule%atom(i)%ttype(1:1).eq.'S'.or. &
						molecule%atom(i)%ttype(1:1).eq.'P') then
			 if(Molecule_Connection_1_2_Check(molecule,(molecule%atom(id)%id),(molecule%atom(i)%id)).eq.1) then
				mark=0
				exit
			 elseif(Molecule_Connection_1_3_Check(molecule,(molecule%atom(id)%id),(molecule%atom(i)%id)).eq.1) then
				mark=0
				exit
			 elseif(Molecule_Connection_1_4_Check(molecule,(molecule%atom(id)%id),(molecule%atom(i)%id)).eq.1) then
				mark=0
				exit
			 end if
		 end if
1		end do
 
    Molecule_Hydrophobic_Neighbor_Check = mark 
end function Molecule_Hydrophobic_Neighbor_Check

integer function Molecule_Connection_1_2_Check(molecule,id1_in, id2_in)
	type(tXMolecule) :: molecule
	integer :: id1, id2,id1_in,id2_in
	integer :: i 

	id1 = id1_in
	id2 = id2_in
	if(id1.eq.id2) then
		Molecule_Connection_1_2_Check = 0
		return
	end if

	do i = 0,molecule%atom(id2-1)%num_neib -1 
		if(id1.eq.molecule%atom(id2-1)%neib(i)) then
			Molecule_Connection_1_2_Check= 1
			return
		end if
	end do
	
	Molecule_Connection_1_2_Check = 0
end function Molecule_Connection_1_2_Check

integer function Molecule_Connection_1_3_Check(molecule,id1_in, id2_in)
	type(tXMolecule)	:: molecule
	integer						:: id1, id2,id1_in,id2_in
	integer ::  i,j 

	id1 = id1_in
	id2 = id2_in
	if(id1.eq.id2) then
		Molecule_Connection_1_3_Check = 0
		return
	end if

	if(Molecule_Connection_1_2_Check(molecule,id1,id2).eq.1) then
		Molecule_Connection_1_3_Check = 0
		return
	end if

	id1 = id1 -1
	id2 = id2 -1

	do i = 0,molecule%atom(id1)%num_neib -1
		do j = 0,molecule%atom(id2)%num_neib -1 
			if(molecule%atom(id1)%neib(i).eq.molecule%atom(id2)%neib(j)) then
				Molecule_Connection_1_3_Check = 1
				return
			end if
		end do
	end do

	Molecule_Connection_1_3_Check = 0
end function Molecule_Connection_1_3_Check

integer function Molecule_Connection_1_4_Check(molecule, id1_in, id2_in)	
	type(tXMolecule) :: molecule
	integer	:: id1, id2,id1_in,id2_in
	integer ::  i,j 

	id1 = id1_in
	id2 = id2_in

	if(id1.eq.id2) then
		Molecule_Connection_1_4_Check = 0
		return
	end if
	if(Molecule_Connection_1_2_Check(molecule,id1,id2).eq.1) then
		Molecule_Connection_1_4_Check = 0
		return
	end if
	if(Molecule_Connection_1_3_Check(molecule,id1,id2).eq.1) then
		Molecule_Connection_1_4_Check = 0
		return 
	end if

	id1 = id1 -1
	id2 = id2 -1

	do i = 0,molecule%atom(id1)%num_neib -1 
		do j = 0,molecule%atom(id2)%num_neib -1 
			if(Molecule_Connection_1_2_Check(molecule,molecule%atom(id1)%neib(i),molecule%atom(id2)%neib(j)).eq.1) then
				Molecule_Connection_1_4_Check = 1
				return
			end if
		end do
	end do

	Molecule_Connection_1_4_Check = 0
end function Molecule_Connection_1_4_Check

integer function Molecule_Connection_1_5_Check(molecule,id1_in, id2_in)
	type(tXMolecule) :: molecule
	integer					 :: id1,id2,id1_in,id2_in
	integer					 ::  i,j 

	id1 = id1_in
	id2 = id2_in

  Molecule_Connection_1_5_Check = 0
	if(id1.eq.id2) return 
	if(Molecule_Connection_1_2_Check(molecule,id1,id2).eq.1) return  
	if(Molecule_Connection_1_3_Check(molecule,id1,id2).eq.1) return 
	if(Molecule_Connection_1_4_Check(molecule,id1,id2).eq.1) return 

	id1 = id1 -1
	id2 = id2 -1

	do i = 0,molecule%atom(id1)%num_neib -1 
		do j = 0,molecule%atom(id2)%num_neib -1 
			if(Molecule_Connection_1_3_Check(molecule,molecule%atom(id1)%neib(i),molecule%atom(id2)%neib(j)).eq.1) then
				Molecule_Connection_1_5_Check = 1
				return
			end if
		end do
	end do

	Molecule_Connection_1_5_Check = 0
end function Molecule_Connection_1_5_Check

integer function Molecule_Connection_1_6_Check(molecule, id1_in, id2_in)
	type(tXMolecule)	:: molecule
	integer						:: id1, id2,id1_in,id2_in
  integer ::  i,j 

	id1 = id1_in
	id2 = id2_in

	Molecule_Connection_1_6_Check = 0
	if(id1.eq.id2) return  
	if(Molecule_Connection_1_2_Check(molecule,id1,id2).eq.1) return  
	if(Molecule_Connection_1_3_Check(molecule,id1,id2).eq.1) return  
	if(Molecule_Connection_1_4_Check(molecule,id1,id2).eq.1) return  
	if(Molecule_Connection_1_5_Check(molecule,id1,id2).eq.1) return  

	id1 = id1 -1
	id2 = id2 -1

	do i = 0,molecule%atom(id1)%num_neib -1 
		do j = 0,molecule%atom(id2)%num_neib -1 
			if(Molecule_Connection_1_4_Check(molecule,molecule%atom(id1)%neib(i),molecule%atom(id2)%neib(j)).eq.1) then
				Molecule_Connection_1_6_Check = 1
				return 
			end if
		end do
	end do

	Molecule_Connection_1_6_Check = 0
end function Molecule_Connection_1_6_Check

type(tXGroup) function Molecule_Find_X_Group(molecule,atom_id)
	type(tXMolecule)	:: molecule
	integer						:: atom_id
  integer						:: i,num 
  type(tXGroup)			:: group 

  ! define the center
  group%center=molecule%atom(atom_id-1) 

  ! find the center's neighbours
  num=group%center%num_neib 

	do i = 0,num -1 
		group%neib(i)=molecule%atom(group%center%neib(i)-1) 
		group%bond(i)=molecule%bond(group%center%bond(i)-1) 
	end do
        
	! count the necessary parameters
	group%num_neib=num 

	group%num_h=0
	group%num_nonh=0

	do i = 0,num -1 
		if(group%neib(i)%ttype(1:1).eq.'H') then
			group%num_h = group%num_h +1 
		else 
			group%num_nonh = group%num_nonh +1 
		end if
	end do

  group%num_hetero=0 

	do i = 0,group%num_nonh -1 
		if(group%bond(i)%ttype.eq.'2') goto 1 
		if(group%bond(i)%ttype.eq.'3') goto 1 
		if(group%bond(i)%ttype.eq.'ar') goto 1 
		if(group%neib(i)%ttype(1:1).eq.'N') then
			group%num_hetero = group%num_hetero +1 
		elseif(group%neib(i)%ttype(1:1).eq.'O') then
			group%num_hetero = group%num_hetero +1 
		end if
1	end do

  group%num_pi=0 

	do i = 0,group%num_nonh -1 
		if(group%bond(i)%ttype.eq.'2') then
			goto 2
		elseif(group%bond(i)%ttype.eq.'3') then
			goto 2
		elseif(group%neib(i)%ttype.eq.'C.ar') then
			group%num_pi = group%num_pi +1 
		elseif(group%neib(i)%ttype.eq.'C.2') then
			group%num_pi = group%num_pi +1 
		elseif(group%neib(i)%ttype.eq.'C.1') then
			group%num_pi = group%num_pi +1 
		elseif(group%neib(i)%ttype.eq.'C.cat') then
			group%num_pi = group%num_pi +1 
		elseif(group%neib(i)%ttype.eq.'N.2') then
			group%num_pi = group%num_pi +1 
		elseif(group%neib(i)%ttype.eq.'N.1') then
			group%num_pi = group%num_pi +1 
		elseif(group%neib(i)%ttype.eq.'N.ar') then
			group%num_pi = group%num_pi +1 
		end if
2	end do

  group%num_nar=0
	group%num_car=0 

	do i = 0,group%num_nonh -1 
		if((group%bond(i)%ttype.ne.'ar')) then
			goto 3
		elseif(group%neib(i)%ttype.eq.'N.ar') then
			group%num_nar = group%num_nar +1 
		elseif(group%neib(i)%ttype.eq.'C.ar') then
			group%num_car = group%num_car +1 
		elseif(group%neib(i)%ttype(1:1).eq.'C') then
			group%num_car = group%num_car +1 
		else 
			group%num_nar = group%num_nar +1 
		end if
3	end do

  group%db_type=0 

	do i = 0,group%num_nonh -1 
		if((group%center%ttype.eq.'O.co2').or. &
			 (group%center%ttype.eq.'O.2').or. &
			 (group%center%ttype.eq.'S.2')) then
				 if(group%neib(i)%ttype(1:1).eq.'C') then
					group%db_type=1 
				 elseif(group%neib(i)%ttype(1:1).eq.'N') then
					group%db_type=2 
				 elseif(group%neib(i)%ttype(1:1).eq.'O') then
					group%db_type=3 
				 elseif(group%neib(i)%ttype(1:1).eq.'S') then
					group%db_type=4 
				 elseif(group%neib(i)%ttype(1:1).eq.'P') then
					group%db_type=5 
				 else 
					goto 4
				 end if
		elseif((group%bond(i)%ttype.eq.'2').or.	&
					 (group%neib(i)%ttype.eq.'O.co2').or.	&
				   (group%neib(i)%ttype.eq.'O.2').or. &
				   (group%neib(i)%ttype.eq.'S.2')) then
						 if(group%neib(i)%ttype(1:1).eq.'C') then
							group%db_type=1 
						 elseif(group%neib(i)%ttype(1:1).eq.'N') then
							group%db_type=2 
						 elseif(group%neib(i)%ttype(1:1).eq.'O') then
							group%db_type=3 
						 elseif(group%neib(i)%ttype(1:1).eq.'S') then
							group%db_type=4 
						 elseif(group%neib(i)%ttype(1:1).eq.'P') then
							group%db_type=5 
						 else 
							goto 4
						 end if
    end if
4	end do
  
	group%valid=1
	
	Molecule_Find_X_Group = group
end function Molecule_Find_X_Group


subroutine Molecule_Get_XLOGP_Type(molecule,atom_id, ttype)
	type(tXMolecule)	:: molecule
	integer						:: atom_id
	character(len=20)	:: ttype
	type(tXGroup)			:: group 
	type(tXAtom)			:: atom 

	atom=molecule%atom(atom_id-1)
	group=Molecule_Find_X_Group(molecule,atom%id)

	ttype = 'Un'  

	if(atom%ttype.eq.'H'.or.atom%ttype.eq.'H.spc') then
		if(group%neib(0)%ttype(1:1).eq.'O') then
			ttype = 'H.hb' 
		elseif(group%neib(0)%ttype(1:1).eq.'N') then
			ttype = 'H.hb' 
		else
			ttype = 'H'
		end if
	end if
	
	if(atom%ttype.eq.'C.3') then
		 if(group%num_nonh.eq.1) then
			 if(group%num_hetero.eq.0) then
				 if(group%num_pi.eq.0) then
					 ttype = 'C.3.h3.pi=0' 
				 else
					 ttype = 'C.3.h3.pi=1'
				 end if
			 else
	       ttype = 'C.3.h3.x' 
			 end if
		 elseif(group%num_nonh.eq.2) then
			 if(group%num_hetero.eq.0) then
				 if(group%num_pi.eq.0) then
					 ttype = 'C.3.h2.pi=0' 
				 elseif(group%num_pi.eq.1) then
					 ttype = 'C.3.h2.pi=1'  
				 else
					 ttype = 'C.3.h2.pi=2' 
				 end if
			 else 
				 if(group%num_pi.eq.0) then
           ttype = 'C.3.h2.x.pi=0' 
		     elseif(group%num_pi.eq.1) then
           ttype = 'C.3.h2.x.pi=1' 
				 else
           ttype = 'C.3.h2.x.pi=2' 
				 end if
			 end if
		 elseif(group%num_nonh.eq.3) then
			 if(group%num_hetero.eq.0) then
				 if(group%num_pi.eq.0) then
					 ttype = 'C.3.h.pi=0' 
				 elseif(group%num_pi.eq.1) then
					 ttype = 'C.3.h.pi=1' 
				 else
					 ttype = 'C.3.h.pi>1' 
				 end if
			 else 
         if(group%num_pi.eq.0) then
           ttype = 'C.3.h.x.pi=0' 
         elseif(group%num_pi.eq.1) then
           ttype = 'C.3.h.x.pi=1' 
         else
           ttype = 'C.3.h.x.pi>1' 
         end if
			 end if
		 elseif(group%num_nonh.eq.4) then
			 if(group%num_hetero.eq.0) then
				 if(group%num_pi.eq.0) then
           ttype = 'C.3.pi=0' 
			   elseif(group%num_pi.eq.1) then
					 ttype = 'C.3.pi=1' 
				 else
           ttype = 'C.3.pi>1' 
         end if
		   else 
				 if(group%num_pi.eq.0) then 
				   ttype = 'C.3.x.pi=0' 
         else 
				   ttype = 'C.3.x.pi>0' 
			   end if
			 end if
		 else	
		   ttype = 'C.3.unknown' 
		 end if
	 end if

	if(atom%ttype.eq.'C.2') then
		 if(group%num_nonh.eq.1) then
			 ttype = 'C.2.h2' 	
		 elseif(group%num_nonh.eq.2) then
			 if(group%num_hetero.eq.0) then
				 if(group%num_pi.eq.0) then
					ttype = 'C.2.h.pi=0' 
				 else 
				  ttype = 'C.2.h.pi=1' 
				 end if
			 else
				 if(group%num_pi.eq.0) then
					 ttype = 'C.2.h.x.pi=0' 
         else 
					ttype = 'C.2.h.x.pi=1' 
				 end if
			 end if
		 elseif(group%num_nonh.eq.3) then
			 if(group%num_hetero.eq.0) then
				 if(group%num_pi.eq.0) then
					ttype = 'C.2.pi=0'  
				 else 
					ttype = 'C.2.pi>0' 
				 end if
			 elseif(group%num_hetero.eq.1) then
				 if(group%num_pi.eq.0) then
					 ttype = 'C.2.x.pi=0'  
         else 
					ttype = 'C.2.x.pi>0' 
				end if
			 else
				 if(group%num_pi.eq.0) then
					ttype = 'C.2.x2.pi=0' 
         else 
				  ttype = 'C.2.x2.pi>0' 
				 end if
				end if
		 else
			 ttype = 'C.2.unknown' 
		 end if
	end if

	 if(atom%ttype.eq.'C.cat') then
		 ttype = 'C.2.x2.pi>0' 
	 end if
	 if(atom%ttype.eq.'C.ar') then 
		 if(group%num_nonh.eq.2) then
			 if(group%num_nar.eq.0) then
			  ttype = 'C.ar.h' 
			 else 
				ttype = 'C.ar.h.(X)' 
			 end if
		 elseif(group%num_nonh.eq.3) then
			 if(group%num_nar.eq.0) then
				 if(group%num_hetero.eq.0) then
					ttype = 'C.ar' 
				 else 
					ttype = 'C.ar.x' 
				 endif
			 else 
				 if(group%num_hetero.eq.0) then
	         ttype = 'C.ar.(X)' 
				 else 
           ttype = 'C.ar.(X).x' 
         end if
			end if
		else	
			ttype = 'C.ar.unknown' 
		end if
	end if

	if(atom%ttype.eq.'C.1') then
		 if(group%db_type.ne.0) then
			ttype = 'C.1..eq.' 
		 elseif(group%num_nonh.eq.1) then
			ttype = 'C.1.h' 
		 elseif(group%num_nonh.eq.2) then
			ttype = 'C.1' 
		 else 
			ttype = 'C.1.unknown' 
		 end if
	end if

	if(atom%ttype.eq.'N.4'.or.atom%ttype.eq.'N.3'.or.atom%ttype.eq.'N.pl3') then
		 if(group%num_nonh.eq.1) then
			 if(group%num_hetero.eq.0) then
				 if(group%num_pi.eq.0) then
					 ttype = 'N.3.h2.pi=0' 
				 else 
					 ttype = 'N.3.h2.pi=1' 
					end if
			 else 
				 ttype = 'N.3.h2.x' 
				end if
		 elseif(group%num_nonh.eq.2) then
			 if(group%num_hetero.eq.0) then
				 if(group%num_pi.eq.0) then
					 ttype = 'N.3.h.pi=0' 
				 elseif(group%num_pi.eq.1) then
					 ttype = 'N.3.h.pi>0' 
				 elseif(atom%ring.eq.0) then
					 ttype = 'N.3.h.pi>0' 
				 else 
					ttype = 'N.3.h.ring' 
				 end if
       else 
				 if(group%num_pi<=1) then
           ttype = 'N.3.h.x' 
			   elseif(atom%ring.eq.0) then
           ttype = 'N.3.h.x' 
				 else
           ttype = 'N.3.h.x.ring' 
         end if
			end if
		 elseif(group%num_nonh.eq.3) then
			 if(group%num_hetero.eq.0) then
				 if(group%num_pi.eq.0) then
					 ttype = 'N.3.pi=0' 
				 elseif(group%num_pi.eq.1) then
					 ttype = 'N.3.pi>0' 
				 elseif(atom%ring.eq.0) then
					 ttype = 'N.3.pi>0' 
				 else
					 ttype = 'N.3.ring' 
				 end if
       else 
         if(group%num_pi<=1) then
					ttype = 'N.3.x' 
         elseif(atom%ring.eq.0) then
          ttype = 'N.3.x' 
         else
          ttype = 'N.3.x.ring' 
         end if
			 end if
		 else
			 ttype = 'N.3.unknown' 
		 end if
	end if

	if(atom%ttype.eq.'N.am') then
		 if(group%num_nonh.eq.1) then
			ttype = 'N.am.h2' 
		 elseif(group%num_nonh.eq.2) then
			 if(group%num_hetero.eq.0) then
				ttype = 'N.am.h' 
			 else 
				ttype = 'N.am.h.x' 
			 end if
		 elseif(group%num_nonh.eq.3) then
			 if(group%num_hetero.eq.0) then
				ttype = 'N.am' 
			 else 
				ttype = 'N.am.x' 
			 end if
		 else 
			ttype = 'N.am.unknown' 
		 end if
	end if

	if(atom%ttype.eq.'N.2') then
		 ! N=C, N=S, N=P
		 if(group%db_type.eq.1.or.group%db_type.eq.4.or.group%db_type.eq.5) then
			 if(group%num_hetero.eq.0) then
				 if(group%num_pi.eq.0) then
					 ttype = 'N.2.(=C).pi=0' 
				 else
					 ttype = 'N.2.(=C).pi=1' 
					end if
			 else
				 if(group%num_pi.eq.0) then
					 ttype = 'N.2.(=C).x.pi=0' 
			   else
					 ttype = 'N.2.(=C).x.pi=1' 
        end if
			end if
		 elseif(group%db_type.eq.2) then
			 if(group%num_hetero.eq.0) then
				ttype = 'N.2.(=N)' 
			 else 
				ttype = 'N.2.(=N).x' 
			 endif
		 elseif(group%db_type.eq.3) then
			 if(group%num_nonh.eq.2) then
				ttype = 'N.2.o' 
			 elseif(group%num_nonh.eq.3) then
				ttype = 'N.2.o2' 
			 else 
				ttype = 'N.2.o' 
			 end if
		 else
			 ttype = 'N.2.unknown' 
		 end if
	end if

	if(atom%ttype.eq.'N.ar') ttype = 'N.ar' 
	if(atom%ttype.eq.'N.1') ttype = 'N.1' 

	if(atom%ttype.eq.'O.3') then
		 if(group%num_nonh.eq.1) then
			 if(group%num_hetero.eq.0) then 
				 if(group%num_pi.eq.0) then
					ttype = 'O.3.h.pi=0' 
				 else 
					ttype = 'O.3.h.pi=1' 
				 end if
			 else 
        ttype = 'O.3.h.x' 
			 end if
		 elseif(group%num_nonh.eq.2) then
			 if(group%num_hetero.eq.0) then
				 if(group%num_pi.eq.0) then
					ttype = 'O.3.pi=0' 
				 else 
					ttype = 'O.3.pi>0' 
				 end if
			 else
				 ttype = 'O.3.x' 
			 end if
		 else
			 ttype = 'O.3.unknown' 
		 end if
	end if

	if(atom%ttype.eq.'O.2') ttype = 'O.2' 
	if(atom%ttype.eq.'O.co2') ttype = 'O.co2' 

	if(atom%ttype.eq.'S.3') then
		 if(group%num_nonh.eq.1) then
			ttype = 'S.3.h' 
		 elseif(group%num_nonh.eq.2) then
			ttype = 'S.3' 
		 else 
			ttype = 'S.3.unknown' 
		 end if
	end if

	if(atom%ttype.eq.'S.2') ttype = 'S.2' 
	if(atom%ttype.eq.'S.o') ttype = 'S.o' 
	if(atom%ttype.eq.'S.o2') ttype = 'S.o2' 

	if(atom%ttype.eq.'P.3') then
		 if(group%db_type.eq.3) then
			ttype = 'P.3.(=O)' 
		 elseif(group%db_type.eq.4) then
			ttype = 'P.3.(=S)' 
		 else 
			ttype = 'P.3.unknown' 
		end if
	end if

	if(atom%ttype.eq.'F') then
   if(group%num_pi.eq.0) then
		ttype = 'F.pi=0' 
   elseif(group%num_pi.eq.1) then
		 ttype = 'F.pi=1' 
		 else 
			ttype = 'F.unknown' 
		end if
	end if

	if(atom%ttype.eq.'Cl') then
	 if(group%num_pi.eq.0) then
		ttype = 'Cl.pi=0' 
	 elseif(group%num_pi.eq.1) then
		ttype = 'Cl.pi=1' 
   else 
		ttype = 'Cl.unknown' 
	 end if
	end if

	if(atom%ttype.eq.'Br') then
   if(group%num_pi.eq.0) then
		 ttype = 'Br.pi=0' 
   elseif(group%num_pi.eq.1) then
		ttype = 'Br.pi=1' 
	 else 
		ttype = 'Br.unknown' 
	 end if
	end if

	if(atom%ttype.eq.'I') then
		 if(group%num_pi.eq.0) then
			ttype = 'I.pi=0' 
		 elseif(group%num_pi.eq.1) then
			 ttype = 'I.pi=1' 
		 else 
			ttype = 'I.unknown' 
		 end if
	end if

	if(atom%ttype.eq.'Si') ttype = 'Si' 
end subroutine Molecule_Get_XLOGP_Type 


! ==========================================================================================================================================
! ==========================================================================================================================================

! XLIGAND MODULE

! ==========================================================================================================================================
! ==========================================================================================================================================
	
subroutine Ligand_Scoring_Stats(ligand,iAtom)
	type(tXLigand)	:: ligand
	integer,optional:: iAtom
	integer					:: i
	real(kind=prec)						:: vdw, hb, hp, hm, hs, rt, score

	vdw = 0
	hb = 0
	hp = 0
	hm = 0
	hs = 0
	rt = 0
	score = 0


	if(present(iAtom)) then
		i = iAtom
		write(*,'(t81,i3, t85,a,   t28,f6.1,  t36,f5.2, t43,f5.2, t51,f5.2, t57,f5.2, t65,f4.1, t73,f6.3)') &
							ligand%mol%atom(i)%id, ligand%mol%atom(i)%ttype,  ligand%abs_inf(i)%vdw,	ligand%abs_inf(i)%hb,		ligand%abs_inf(i)%hp, &
							ligand%abs_inf(i)%hm,  ligand%abs_inf(i)%hs,		  ligand%abs_inf(i)%rt,		ligand%abs_inf(i)%score
	else
		write(*,'(t31,a3,  t38,a2,  t46,a2,  t52,a2,   t60,a2,&
		    t67,a2,   t74,a5)')  'VDW',   'HB',    'HP',    'HM', &
		    'HS',       'RT',			'SCORE'
		do i = 0,ligand%mol%num_atom -1
		!	write(*,'(t28,f6.1,  t36,f4.1, t43,f5.1, t51,f3.1, t57,f5.1, t65,f4.1, t74,f5.2)')
			write(*,'(t1,i3, t6,a,   t28,f6.1,  t36,f5.2, t43,f5.2, t51,f5.2, t57,f5.2, t65,f4.1, t73,f6.3)') &
			    ligand%mol%atom(i)%id, ligand%mol%atom(i)%ttype, &
			    ligand%abs_inf(i)%vdw, ligand%abs_inf(i)%hb, &
			    ligand%abs_inf(i)%hp, &
			    ligand%abs_inf(i)%hm, ligand%abs_inf(i)%hs,	&
			    ligand%abs_inf(i)%rt, ligand%abs_inf(i)%score
			
			vdw		= vdw		+ ligand%abs_inf(i)%vdw
			hb		= hb		+ ligand%abs_inf(i)%hb
			hp		= hp		+ ligand%abs_inf(i)%hp
			hm		= hm		+ ligand%abs_inf(i)%hm
			hs		= hs		+	ligand%abs_inf(i)%hs
			rt		= rt		+ ligand%abs_inf(i)%rt
			score = score + ligand%mol%atom(i)%score

		end do

		write(*,'(a)') '                            --------------------------------------------------'
	end if

end subroutine Ligand_Scoring_Stats

real(kind=prec) function Ligand_Calculate_Binding_Score(ligand,input,protein)
	type(tXLigand)		:: ligand
	type(tXInput)			:: input
	type(tXProtein)		:: protein
	integer						:: i 
	type(tXLigand)		:: bak 
	integer						:: num_nonh=0
	real(kind=prec)							:: sum=0.000 

	bak = ligand   ! backup the atomic charge information

	! these scoring functions may need formal charges
	call Molecule_Assign_Apparent_Charge(ligand%mol) 


	! prepare for atomic binding score
	if(.not. associated(ligand%abs_inf)) then
		allocate(ligand%abs_inf(0:ligand%mol%num_atom -1),stat=err)	! create array of atomic binding score, one element for each atom in ligand
		if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'
	end if

	! initialize binding scores to 0%0
	ligand%abs_inf(:)%pkd1	= 0.000 
	ligand%abs_inf(:)%pkd2	= 0.000 
	ligand%abs_inf(:)%pkd3	= 0.000 
	ligand%abs_inf(:)%vdw		= 0.000 
	ligand%abs_inf(:)%hb		= 0.000 
	ligand%abs_inf(:)%hm		= 0.000 
	ligand%abs_inf(:)%hp		= 0.000 
	ligand%abs_inf(:)%hs		= 0.000 
	ligand%abs_inf(:)%rt		= 0.000 
	ligand%abs_inf(:)%score	= 0.000 

	! clean other variables
	ligand%pkd1	= 0 
	ligand%pkd2	= 0
	ligand%pkd3	= 0
	ligand%vdw	= 0
	ligand%hb		= 0
	ligand%hm		= 0
	ligand%hp		= 0
	ligand%hs		= 0
	ligand%rt		= 0 

	! first, generate dot surface
	call Molecule_Generate_Surface_Dots(ligand%mol, WATER_R) 

	! now calculate the van der Waals term
	ligand%vdw = Ligand_Calculate_VDW(ligand,protein) 
	ligand%abs_inf(:)%vdw = ligand%mol%atom(:)%score

	! now calculate the H-bond term
	ligand%hb=Ligand_Calculate_HB(ligand,protein) 
  ligand%abs_inf(:)%hb=ligand%mol%atom(:)%score 

	! now calculate the rotor term
	ligand%rt=Ligand_Calculate_RT(ligand,protein) 
	ligand%abs_inf(:)%rt=ligand%mol%atom(:)%score 

	! now calculate the hydrophobic term
	if(input%apply_hpscore.eq.'YES') then
		ligand%hp=Ligand_Calculate_HP(ligand,protein) 
		ligand%abs_inf(:)%hp=ligand%Mol%atom(:)%score 
	  ligand%pkd1=input%hpscore_c0 +				&
                input%hpscore_cvdw*ligand%vdw	+ &
                input%hpscore_chb*ligand%hb	+		&
                input%hpscore_chp*ligand%hp	+		&
                input%hpscore_crt*ligand%rt 
	end if
	
 	if(input%apply_hmscore.eq.'YES') then
	  ligand%hm=Ligand_Calculate_HM(ligand,protein) 
    ligand%abs_inf(:)%hm=ligand%mol%atom(:)%score 
		ligand%pkd2=input%hmscore_c0+	&
                    input%hmscore_cvdw*ligand%vdw+	&
                    input%hmscore_chb*ligand%hb+		&
                    input%hmscore_chm*ligand%hm+		&
              	    input%hmscore_crt*ligand%rt 

	end if

 	if(input%apply_hsscore.eq.'YES') then
	  ligand%hs=Ligand_Calculate_HS(ligand,protein) 
		ligand%abs_inf(:)%hs=ligand%mol%atom(:)%score 
	  ligand%pkd3=input%hsscore_c0+	&
                input%hsscore_cvdw* ligand%vdw+	&
                input%hsscore_chb * ligand%hb+		&
                input%hsscore_chs * ligand%hs+		&
                input%hsscore_crt * ligand%rt 
	end if

	ligand%bind_score=(ligand%pkd1+ligand%pkd2+ligand%pkd3)/(input%num_method) 

	! now compute atomic binding scores
	num_nonh = Molecule_Get_Num_Heavy_Atom(ligand%mol)

	if(input%apply_hpscore.eq.'YES') then
		do i = 0,ligand%mol%num_atom -1 
			 if(ligand%mol%atom(i)%valid<=0) goto 1
			 if(ligand%mol%atom(i)%ttype.eq.'H') goto 1

	 	 	 ligand%abs_inf(i)%pkd1 =  input%hpscore_c0/num_nonh + &
																 input%hpscore_cvdw * ligand%abs_inf(i)%vdw + &
																 input%hpscore_chb  * ligand%abs_inf(i)%hb +	&
																 input%hpscore_chp  * ligand%abs_inf(i)%hp +	&
																 input%hpscore_crt  * ligand%abs_inf(i)%rt 
1		end do
	end if

	if(input%apply_hmscore.eq.'YES') then
		do i = 0,ligand%mol%num_atom -1 
			if(ligand%mol%atom(i)%valid<=0) goto 2
			if(ligand%mol%atom(i)%ttype.eq.'H') goto 2

	 	 	ligand%abs_inf(i)%pkd2 = input%hmscore_c0/num_nonh + &
															 input%hmscore_cvdw * ligand%abs_inf(i)%vdw + &
															 input%hmscore_chb  * ligand%abs_inf(i)%hb + &
															 input%hmscore_chm	* ligand%abs_inf(i)%hm + &
															 input%hmscore_crt	* ligand%abs_inf(i)%rt 
2		end do
	end if

	if(input%apply_hsscore.eq.'YES') then
		do i = 0,ligand%mol%num_atom -1 
			if(ligand%mol%atom(i)%valid<=0) goto 3
			if(ligand%mol%atom(i)%ttype.eq.'H') goto 3

	 	 	ligand%abs_inf(i)%pkd3 = input%hsscore_c0/num_nonh + &
															 input%hsscore_cvdw * ligand%abs_inf(i)%vdw + &
															 input%hsscore_chb	* ligand%abs_inf(i)%hb + &
															 input%hsscore_chs	* ligand%abs_inf(i)%hs + &
															 input%hsscore_crt	* ligand%abs_inf(i)%rt 
3		end do
	end if

	sum = 0
	do i = 0,ligand%mol%num_atom -1 
		 ligand%abs_inf(i)%score=0
		 ligand%abs_inf(i)%score=ligand%abs_inf(i)%score+ligand%abs_inf(i)%pkd1 
		 ligand%abs_inf(i)%score=ligand%abs_inf(i)%score+ligand%abs_inf(i)%pkd2 
		 ligand%abs_inf(i)%score=ligand%abs_inf(i)%score+ligand%abs_inf(i)%pkd3 
		 ligand%abs_inf(i)%score=ligand%abs_inf(i)%score/input%num_method 
		 ligand%mol%atom(i)%score=ligand%abs_inf(i)%score 
		 sum= sum + ligand%mol%atom(i)%score 
	end do

	if(input%show_abs.eq.'YES') then ! use abs as charges
	 	ligand%mol%charge_type = 'USER_CHARGES'
		ligand%mol%atom(:)%q = ligand%abs_inf(:)%score
	else		! restore original charges
		ligand%mol%charge_type = bak%mol%charge_type 
 		ligand%mol%atom(:)%q=bak%mol%atom(:)%q 
	end if


	Ligand_Calculate_Binding_Score =  ligand%bind_score 
end function Ligand_Calculate_Binding_Score

subroutine Ligand_Sum_HBonds(ligand, protein, num_candidate, candidate)
	type(tXLigand)	:: ligand
	type(tXProtein)	:: protein
	integer					:: num_candidate
	type(tXHBond)		:: candidate(0:999)
	type(tXHBond)		:: tmp_hbond		! used in swap
	integer					:: i,j,k,count,limit,latom,patom 
	real(kind=prec)						:: angle,v1(0:2),v2(0:2) 
	type(tXGroup)		:: tmp_group 

	! Step 1: rank candidates according their strength in decreasing order
	! note 'fabs' is applied because SB strength could be negative
	do i = 0,num_candidate-1 -1 
		do j = i+1,num_candidate -1 
			if(abs(candidate(i)%score)>=abs(candidate(j)%score)) goto 1
!			write(*,*) 'swapping ', i, j
			!call SWAP(candidate(i),candidate(j))
			call HBond_operator_copy(tmp_hbond,candidate(i))
			call HBond_operator_copy(candidate(i),candidate(j))
			call HBond_operator_copy(candidate(j),tmp_hbond)
1   end do
	end do

	! Step 2: check the angular limit: the angle between any two H-bonds
	! on the same atom must be larger than 45 degrees
	! note that this filter could be applied to both ligand and protein
	do i = 0,num_candidate-1 -1 
		do j = i+1,num_candidate -1 
      if(candidate(i)%latom.ne.candidate(j)%latom) goto 2
			do k = 0,3 -1  
				latom = candidate(i)%latom
				patom = candidate(i)%patom
				v1(k) = protein%atom(patom-1)%coor(k)-ligand%mol%atom(latom-1)%coor(k) 
				latom = candidate(j)%latom
				patom = candidate(j)%patom
				v2(k) = protein%atom(patom-1)%coor(k)-ligand%mol%atom(latom-1)%coor(k) 
			end do

			angle=abs(Angle_of_Two_Vectors(v1,v2)) 

			if(angle<45.0) then
				candidate(j)%score=0.000 
			else 
				goto 2
			end if
2		end do
	end do

! THIS SECTION NOT USED IN ORIGINAL SOURCE
!	do i = 0,num_candidate-1 -1 
!		do j = i,num_candidate -1 
!			if(candidate(i)%patom.ne.candidate(j)%patom) goto 3
!			do k = 0,3 -1 
!				latom = candidate(i)%latom
!				patom = candidate(i)%patom
!				v1(k) = ligand%mol%atom(latom-1)%coor(k)-protein%atom(patom-1)%coor(k) 
!				latom = candidate(j)%latom
!				patom = candidate(j)%patom
!				v2(k) = ligand%mol%atom(latom-1)%coor(k)-protein%atom(patom-1)%coor(k) 
!      end do
!
!      angle=abs(Angle_of_Two_Vectors(v1,v2)) 
!      if(angle<45.0) then
!				candidate(j)%score=0.000 
!      else 
!				goto 3
!			end if
!3		end do
!	end do
 
	! Step 3: a donor atom shall not form more H-bonds than the H atoms it has
	! note that this filter is applied only to the ligand side
	do i = 0,num_candidate-1 -1 
		if(candidate(i)%ttype.ne.1) goto 4

		count=1
		latom=candidate(i)%latom

		if(ligand%mol%atom(latom-1)%hb.eq.'DA') then
			limit=1  
		else
			tmp_group=Molecule_Find_A_Group(ligand%mol,latom) 
			limit=tmp_group%num_h 
		end if

		do j = i+1,num_candidate -1 
			if(candidate(j)%ttype.ne.1) goto 5
			if(candidate(i)%latom.ne.candidate(j)%latom) goto 5

			count = count +1 

			if(count<=limit) then
				goto 5 
			else 
				candidate(j)%score=0.000 
			end if
5		end do
4	end do

	! Step 4: an acceptor atom shall not form more H-bonds than its LPs 
	! note that this filter is applied only to the ligand side
	do i = 0,num_candidate-1 -1 
		if(candidate(i)%ttype.ne.2.and.candidate(i)%ttype.ne.3) goto 6

		count=1
		latom=candidate(i)%latom

		if(ligand%mol%atom(latom-1)%ttype(1:1).eq.'O') then
			limit=2 
		elseif(ligand%mol%atom(latom-1)%ttype(1:1).eq.'N') then
			limit=1 
		elseif(ligand%mol%atom(latom-1)%ttype(1:1).eq.'S') then
			limit=2 
		end if

		do j = i+1,num_candidate -1 
			if(candidate(j)%ttype.ne.2.and.candidate(j)%ttype.ne.3) goto 7
			if(candidate(i)%latom.ne.candidate(j)%latom) goto 7

			count = count +1 

			if(count<=limit) then
				goto 7
			else 
				candidate(j)%score=0.000 
			end if
7		end do
6	end do
end subroutine Ligand_Sum_HBonds

real(kind=prec) function Ligand_Calculate_RT(ligand,protein)
	type(tXLigand) :: ligand
	type(tXProtein) :: protein
	integer ::  i,j,mark,tmp 
	integer ::  id1,id2 
	real(kind=prec) ::  sum 

	! if a single bond is normal, bond%valid=1 
	! if a single bond is a rotor, bond%valid=2 

	! clean the variables
	ligand%mol%atom(:)%score=0.000 

	!call Molecule_Show_Contents(ligand%mol)
	!chSlask = getcharqq()

	! eliminate all the bonds in_ring and all the non-single bonds
	do i = 0,ligand%mol%num_bond -1 
		 if(ligand%mol%bond(i)%ring>0) then
!			write(*,*) 'in ring: ', i, ligand%mol%bond(i)%atom_1, ligand%mol%bond(i)%atom_2
			goto 1
		 elseif(ligand%mol%bond(i)%ttype.ne.'1') then
!			write(*,*) 'non-single: ', i, ligand%mol%bond(i)%atom_1, ligand%mol%bond(i)%atom_2
			goto 1
		 end if
		 ligand%mol%bond(i)%valid=2 
1	end do

	! eliminate all the R-H, R-X, R-OH, R-NH2, R-CH3 bonds
	do i = 0,ligand%mol%num_bond -1 
		 if(ligand%mol%bond(i)%valid.ne.2) then
!			write(*,*) 'R- bond: ', i, ligand%mol%bond(i)%atom_1, ligand%mol%bond(i)%atom_2
			goto 2
		 end if
		 id1=ligand%mol%bond(i)%atom_1
		 id2=ligand%mol%bond(i)%atom_2
		 if(ligand%mol%atom(id1-1)%num_nonh.eq.1.or.ligand%mol%atom(id2-1)%num_nonh.eq.1) ligand%mol%bond(i)%valid=1 
2	end do

	! sp2-sp2 rotors
	do i = 0,ligand%mol%num_bond -1 
		 if(ligand%mol%bond(i)%valid.ne.2) goto 3
		 id1=ligand%mol%bond(i)%atom_1
		 id2=ligand%mol%bond(i)%atom_2
		 mark=0 

		 
		 tmp=Molecule_Get_Atom_Hybridizing_Type(ligand%mol,ligand%mol%atom(id1-1)%ttype) 
		 if((tmp.eq.1).or.(tmp.eq.2)) mark = mark +1 

		 tmp=Molecule_Get_Atom_Hybridizing_Type(ligand%mol,ligand%mol%atom(id2-1)%ttype) 
     if((tmp.eq.1).or.(tmp.eq.2)) mark = mark +1 
 
		 if(mark.eq.2) then
			ligand%mol%bond(i)%valid=1
			goto 3
		 end if
3	end do

	! eliminate terminal rotors, e%g% -PO3, -CF3, -CMe3, -NMe3
	do i = 0,ligand%mol%num_bond -1 
     if(ligand%mol%bond(i)%valid.ne.2) goto 4

		 id1=ligand%mol%bond(i)%atom_1
		 id2=ligand%mol%bond(i)%atom_2

		 if(Molecule_Judge_Terminal_Atom(ligand%mol,ligand%mol%atom(id1-1),id2).eq.1) then
			 ligand%mol%bond(i)%valid=1 
		 elseif(Molecule_Judge_Terminal_Atom(ligand%mol,ligand%mol%atom(id2-1),id1).eq.1) then
			 ligand%mol%bond(i)%valid=1 
		 else
			goto 4
		 end if
4	end do

	! eliminate abnormal rotors
	do i = 0,ligand%mol%num_bond -1
    if(ligand%mol%bond(i)%valid.ne.2) goto 5
    id1=ligand%mol%bond(i)%atom_1
	  id2=ligand%mol%bond(i)%atom_2

		if(ligand%mol%atom(id1-1)%valid<=0) then
			ligand%mol%bond(i)%valid=1 
		elseif(ligand%mol%atom(id2-1)%valid<=0) then	
			ligand%mol%bond(i)%valid=1 
		else 
			goto 5
		end if
5	end do



	! now count the frozen rotors, all the rotors have been labeled as 2
	! notice: the following part is different from subroutine Molecule_Count_Rotor()
	sum=0.000 
	do i = 0,ligand%mol%num_atom -1
		 if(ligand%mol%atom(i)%valid<=0) goto 6

		 mark=0 

		 do j = 0,ligand%mol%atom(i)%num_neib -1
			 tmp=ligand%mol%atom(i)%bond(j)
			 if(tmp.eq.0) then
				goto 7
			 elseif(ligand%mol%bond(tmp-1)%valid.ne.2) then
				goto 7
			 else 
				mark = mark +1 
			 end if
7		 end do

		 if(mark.eq.1) then
			ligand%mol%atom(i)%score=ligand%mol%atom(i)%score+0.50 
		 elseif(mark.eq.2) then
			ligand%mol%atom(i)%score=ligand%mol%atom(i)%score+1.00 
		 elseif(mark>=3) then	
			ligand%mol%atom(i)%score=ligand%mol%atom(i)%score+0.50 
		 end if
		 
		 sum=sum+ligand%mol%atom(i)%score 
6	end do

	Ligand_Calculate_RT = sum 
end function Ligand_Calculate_RT

real(kind=prec) function Ligand_Calculate_HP(ligand,protein)
	type(tXLigand)	:: ligand
	type(tXProtein)	:: protein 
	integer					:: i,j 
	real(kind=prec)						:: d,d1,d2,cutoff 
	real(kind=prec)						:: tmp,asum,sum 

	sum=0.000 
	ligand%mol%atom(:)%score=0.000 

	cutoff=DIST_CUTOFF 
	
	do i = 0,ligand%mol%num_atom -1 
		if(ligand%mol%atom(i)%valid<=0) goto 1
		if(ligand%mol%atom(i)%hb.ne.'H') goto 1

		asum=0.000 

		do j = 0,protein%num_atom -1 
			if(protein%atom(j)%valid.ne.2) then
				goto 2
			elseif(protein%atom(j)%ttype.eq.'H') then
				goto 2
			elseif(protein%atom(j)%ttype.eq.'O.w') then
				goto 2
			elseif(protein%atom(j)%hb.ne.'H') then
				goto 2
			end if

			d=Distance(ligand%mol%atom(i)%coor, protein%atom(j)%coor) 
			if(d>=cutoff) goto 2

			d1=ligand%mol%atom(i)%R+protein%atom(j)%R+0.500 
			d2=ligand%mol%atom(i)%R+protein%atom(j)%R+2.200 

      if(d<d1) then
				tmp=1.000 
      elseif(d<d2) then
				tmp=(1/(d1-d2))*(d-d2) 
      else
				tmp=0.000 
			end if

      asum=asum+tmp 
2		end do

		ligand%mol%atom(i)%score=asum
		sum=sum+asum
1	end do

	Ligand_Calculate_HP = sum 
end function Ligand_Calculate_HP

real(kind=prec) function Ligand_Calculate_HM(ligand,protein)
	type(tXLigand)	:: ligand
	type(tXProtein)	:: protein
	integer ::  i,j 
	real(kind=prec) ::  asum,sum,tmp,d,d1,d2,total,cutoff 

	ligand%mol%atom(:)%score=0.000 

	sum=0.000
	cutoff=DIST_CUTOFF

	do i = 0,ligand%mol%num_atom -1 
		if(ligand%mol%atom(i)%valid<=0) goto 1
		if(ligand%mol%atom(i)%ttype.eq.'H') goto 1 
		if(ligand%mol%atom(i)%hb.ne.'H') goto 1 
		if(ligand%mol%atom(i)%logp<=0.00) goto 1 

		total=0.000
		asum=0.000

		do j = 0,protein%num_atom -1 
			if(protein%atom(j)%valid.ne.2) goto 2
			if(protein%atom(j)%ttype.eq.'H') goto 2
			if(protein%atom(j)%ttype.eq.'O.w') goto 2

			d=Distance(ligand%mol%atom(i)%coor, protein%atom(j)%coor) 
			if(d>cutoff) goto 2

			d1=ligand%mol%atom(i)%R+protein%atom(j)%R+0.50 
			d2=ligand%mol%atom(i)%R+protein%atom(j)%R+2.20 

			if(d<d1) then	
				tmp=1.000  
			elseif(d<d2) then
				tmp=(1/(d1-d2))*(d-d2) 
			else 
				tmp=0.000 
			end if

			total=total + (protein%atom(j)%logp*tmp) 
2		end do

		if(ligand%mol%atom(i)%logp>=0.50) then
			asum=ligand%mol%atom(i)%logp 
 		elseif(total>-0.50) then
			asum=ligand%mol%atom(i)%logp 
		else 
			asum=0.000 
		end if

		ligand%mol%atom(i)%score=asum
		sum=sum+asum
1	end do

 Ligand_Calculate_HM = sum 
end function Ligand_Calculate_HM

real(kind=prec) function Ligand_Calculate_HS(ligand,protein)
	type(tXLigand)	:: ligand
	type(tXProtein)	:: protein
	integer ::  i 
	real(kind=prec) ::  sum,total,buried
	real(kind=prec) :: slask

	! clear the variables first
	sum=0.000 
	ligand%mol%atom(:)%score=0.000 

	! then get the buried surface atom by atom
	do i = 0,ligand%mol%num_atom -1 
	 if(ligand%mol%atom(i)%valid<=0) goto 1 
	 if(ligand%mol%atom(i)%ttype.eq.'H') goto 1 
	 if(ligand%mol%atom(i)%hb.eq.'DH') goto 1 
	 if(ligand%mol%atom(i)%hb.eq.'D') goto 1 
	 if(ligand%mol%atom(i)%hb.eq.'A') goto 1 
	 if(ligand%mol%atom(i)%hb.eq.'DA') goto 1 
	 if(ligand%mol%atom(i)%hb.eq.'P') goto 1 
	 if(ligand%mol%atom(i)%hb.eq.'N') goto 1 

	 slask = Ligand_Atom_Buried_Surface(ligand,protein,i+1,total,buried) 
				
	 sum=sum+buried
	 ligand%mol%atom(i)%score=buried

	 ! sum+=(buried*this->atom(i)%logp) 
	 ! this option works *slightly* better
1	end do

  Ligand_Calculate_HS = sum 
end function Ligand_Calculate_HS

real(kind=prec) function Ligand_Calculate_VDW(ligand,protein)
	type(tXLigand)	:: ligand
	type(tXProtein)	:: protein
	integer ::  i,j 
	real(kind=prec) ::  cutoff,d0,d,tmp,tmp1,tmp2,sum,asum 

	! clear the variables
  cutoff=DIST_CUTOFF
	sum=0.000
	
	ligand%mol%atom(:)%score=0.000 

	! now calculate the P-L vdw interaction
	do i = 0,ligand%mol%num_atom -1 
		if(ligand%mol%atom(i)%valid<=0) goto 1
		 if(ligand%mol%atom(i)%ttype.eq.'H') goto 1

		asum=0.000 

		do j = 0,protein%num_atom -1 
			if(protein%atom(j)%valid.ne.2) goto 2
			if(protein%atom(j)%ttype.eq.'H') goto 2
			if(protein%atom(j)%ttype.eq.'O.w') goto 2

			d0=ligand%mol%atom(i)%R+protein%atom(j)%R 
			d=Distance(ligand%mol%atom(i)%coor,protein%atom(j)%coor) 

			if(d>cutoff) goto 2

			tmp1=d0/d 
			tmp1=tmp1*tmp1*tmp1*tmp1
			tmp2=tmp1*tmp1
			tmp=tmp2-2.0*tmp1 

			asum = asum + tmp 
2		end do

		! revert the sign
 		asum=asum*(-1.00) 

		! if ligand atom has unfavorable contribution then neglect it		=============== FRÅGA NERVALL OM DETTA ===============
 		if(asum<0.00) then
			goto 1
		else 
			sum=sum+asum
			ligand%mol%atom(i)%score=asum
		end if
1	end do

	Ligand_Calculate_VDW = sum 
end function Ligand_Calculate_VDW

real(kind=prec) function Ligand_Calculate_HB(ligand,protein)
	type(tXLigand)	:: ligand
	type(tXProtein) :: protein
	integer		::  i
	real(kind=prec)			::  sum 
	integer		::  num_candidate 
	type(tXHBond),allocatable :: candidate(:)

	! clear the variables
	sum=0.000
	num_candidate=0
	ligand%mol%atom(:)%score=0.000 

	allocate(candidate(0:1000-1),stat=err)					! worst case
	if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'
	
	! first, get the HB pairs between protein and ligand
	num_candidate=Ligand_Get_HBond_Pair_PL(ligand,protein,candidate) 
!	write(*,*) 'num_candidate = ', num_candidate
!	do i = 0,num_candidate -1
!		write(*,*) candidate(i)%valid,candidate(i)%ttype,candidate(i)%d_type,candidate(i)%a_type,candidate(i)%sb,candidate(i)%latom,candidate(i)%patom
!	end do
!	chSlask = getcharqq()

	call Ligand_Sum_HBonds(ligand,protein,num_candidate,candidate) 

	! then sum their contributions
	do i = 0,num_candidate -1 
		if(abs(candidate(i)%score)<0.01) then
			goto 1
		else
	 	 sum=sum+candidate(i)%score 
		 ligand%mol%atom(candidate(i)%latom-1)%score= ligand%mol%atom(candidate(i)%latom-1)%score + candidate(i)%score 
		end if
1	end do

 	deallocate(candidate)
	Ligand_Calculate_HB = sum 
end function Ligand_Calculate_HB

integer function Ligand_Get_HBond_Pair_PL(ligand,protein,candidate)
	type(tXLigand)	:: ligand
	type(tXPRotein)	:: protein
	type(tXHBond)		:: candidate(0:999)
	integer ::  i,j,num,ttype ,slask
	integer ::  sb 
	real(kind=prec) ::  d,cutoff 
	type(tXHbond)	:: tmp_candidate 

	num=0
	cutoff=5.00

	! note that the following H-bond algorithm is not based on any
	! explicit hydrogen atom at all, and this is the right thing to do%
	do i = 0,ligand%mol%num_atom -1 
		if(ligand%mol%atom(i)%valid<=0) goto 1
		if(ligand%mol%atom(i)%ttype.eq.'H') goto 1
		if(ligand%mol%atom(i)%hb.eq.'N') goto 1
		if(ligand%mol%atom(i)%hb.eq.'H') goto 1
		if(ligand%mol%atom(i)%hb.eq.'P') goto 1
		if(ligand%mol%atom(i)%hb.eq.'DH') goto 1

		do j = 0,protein%num_atom -1 
			if(protein%atom(j)%valid.ne.2) goto 2
			if(protein%atom(j)%ttype.eq.'H') goto 2
			if(protein%atom(j)%ttype.eq.'O.w') goto 2
			if(protein%atom(j)%hb.eq.'H') goto 2
			if(protein%atom(j)%hb.eq.'P') goto 2
			if(protein%atom(j)%hb.eq.'N') goto 2

			! determine the ttype of ligand H-bond first
			! ttype=0, no H-bond
			! ttype=1, latom is the donor.eq. patom is the acceptor 
			! ttype=2, latom is the acceptor.eq. patom is the donor 
			! ttype=3, latom bound with metal ion  
			if(ligand%mol%atom(i)%hb.eq.'D') then
				if(protein%atom(j)%hb.eq.'D') then 
					ttype=0 
				elseif(protein%atom(j)%hb.eq.'A') then 
					ttype=1 
				elseif(protein%atom(j)%hb.eq.'DA') then 
					ttype=1 
				elseif(protein%atom(j)%hb.eq.'M') then 
					ttype=0 
				else 
					ttype=0 
				end if
			elseif(ligand%mol%atom(i)%hb.eq.'A') then
		 	 if(protein%atom(j)%hb.eq.'D') then 
				ttype=2 
		 	 elseif(protein%atom(j)%hb.eq.'A') then 
				ttype=0 
     	 elseif(protein%atom(j)%hb.eq.'DA') then 
				ttype=2 
		 	 elseif(protein%atom(j)%hb.eq.'M') then 
				ttype=3 
     	 else 
				ttype=0 
			 end if
		 elseif(ligand%mol%atom(i)%hb.eq.'DA') then
		 	 if(protein%atom(j)%hb.eq.'A') then 
				ttype=1 
		 	 elseif(protein%atom(j)%hb.eq.'D') then 
				ttype=2 
		 	 elseif(protein%atom(j)%hb.eq.'DA') then 
				ttype=1 
		 	 elseif(protein%atom(j)%hb.eq.'M') then 
				ttype=3 
		 	 else 
				ttype=0 
			 end if
		 else 
			 ttype=0 
		 end if

		 if(ttype.eq.0) goto 2		! no H-bond

		 ! a crude distance check
		 d=Distance(ligand%mol%atom(i)%coor,protein%atom(j)%coor) 
		 if(d>cutoff) goto 2

		 ! this section is used to differentiate HB and SB
		 if((ligand%mol%atom(i)%ttype.eq.'O.co2').and.(protein%atom(j)%ttype.eq.'O.co2')) then
			 sb=0
	 	 elseif(abs(ligand%mol%atom(i)%q)>0.01.and.abs(protein%atom(j)%q)>0.01)  then
			 sb=1
		 else 
			 sb=0
		 end if

		 ! now handle ligand h-bond
		 call HBond_Clear(tmp_candidate)
		 tmp_candidate%ttype=ttype 
		 tmp_candidate%sb=sb 
		 tmp_candidate%latom=i+1 
		 tmp_candidate%patom=j+1 

		 if(ttype.eq.1) then
			 tmp_candidate%D=ligand%mol%atom(i) 
			 tmp_candidate%A=protein%atom(j) 
		 elseif(ttype.eq.2) then
			 tmp_candidate%A=ligand%mol%atom(i) 
			 tmp_candidate%D=protein%atom(j) 
		 elseif(ttype.eq.3) then
			 tmp_candidate%A=ligand%mol%atom(i) 
			 tmp_candidate%D=protein%atom(j) 
		 end if

		 slask = HBond_Value_HBond_2(tmp_candidate) 

		 if(abs(tmp_candidate%score)>0.000) then
			candidate(num)=tmp_candidate
			num = num +1
		 else
			goto 2
		 end if
2		end do
1	end do

 Ligand_Get_HBond_Pair_PL =num 
end function Ligand_Get_HBond_Pair_PL

subroutine Ligand_Value_Atom(ligand,lite_run)
	type(tXLigand),intent(inout)			:: ligand
	integer,optional									:: lite_run
	integer														:: mark

	if(.not. present(lite_run)) mark = Molecule_ValueAtom(ligand%mol)

	! now compute the H-bond root for each HB atom
	call Ligand_Calculate_HB_Root(ligand)

	! assign atomic logp values, which is needed by binding score
	! atomic logp is not read from ATOM_DEF_XTOOL!
	if(.not. present(lite_run)) ligand%mol%logp = Molecule_Calculate_LogP(ligand%mol)

	! some other molecular properties
	if(.not. present(lite_run)) ligand%mol%num_hb_atom  = Molecule_Get_Num_HB_Atom(ligand%mol)
	if(.not. present(lite_run)) ligand%mol%num_rotor		= Molecule_Count_Rotor(ligand%mol)
end subroutine Ligand_Value_Atom

subroutine Ligand_Calculate_HB_Root(lig)
	type(tXLigand) :: lig
	integer :: i,j,tmp,num_nonh
	real(kind=prec) :: tmpx,tmpy,tmpz

	do i = 0,lig%mol%num_atom -1
		if(lig%mol%atom(i)%valid==0) then
			goto 1
		else if(lig%mol%atom(i)%hb.eq.'N') then
			goto 1
		else if(lig%mol%atom(i)%hb.eq.'H') then
			goto 1
		else if(lig%mol%atom(i)%hb.eq.'P') then
			goto 1
		end if

		tmpx=0
		tmpy=0
		tmpz=0.000
		num_nonh=0
		
		do j = 0,lig%mol%atom(i)%num_neib -1
			tmp=lig%mol%atom(i)%neib(j)-1
			if(lig%mol%atom(tmp)%ttype(1:1)=='H') then
				goto 2
			else
				tmpx=tmpx+lig%mol%atom(tmp)%coor(0)
				tmpy=tmpy+lig%mol%atom(tmp)%coor(1)
				tmpz=tmpz+lig%mol%atom(tmp)%coor(2)
				num_nonh = num_nonh + 1
			end if
2		end do

		if(num_nonh==0) then
			lig%mol%atom(i)%hb = 'P'
		else
		 	 tmpx=tmpx/num_nonh
			 tmpy=tmpy/num_nonh
			 tmpz=tmpz/num_nonh
			 lig%mol%atom(i)%root(0)=tmpx
			 lig%mol%atom(i)%root(1)=tmpy
			 lig%mol%atom(i)%root(2)=tmpz
		end if
1	end do
end subroutine Ligand_Calculate_HB_Root


!################################################################
!# Translate list of q-atoms and bonds to xscore-type ligand.
!# If offset input eq 0 then 
!################################################################
subroutine xligand_translate(ligand,nAtoms,nBonds,aQ,aB)
	type(tXLigand),intent(out)				:: ligand
	integer,intent(in)						:: nAtoms, nBonds
	type(q_atom),dimension(:),intent(in)	:: aQ
	type(q_bond),dimension(:),intent(in)	:: aB
!	integer																:: offset
	integer																:: i, mark,iRes,atom_res_id
	type(tXBond)													:: tmp

	! set number of qatoms
	ligand%mol%num_atom = nAtoms
	allocate(ligand%mol%atom(0:nAtoms -1),stat=err)		! old c++ code req. 0:...
	if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'
	ligand%mol%num_bond = nBonds
	allocate(ligand%mol%bond(0:nBonds -1),stat=err)		! old c++ code req. 0:...
	if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'

	!Set residue offset, provided that ligand is after protein 
	if (offset .eq. 0) then
		do i=1,nres_solute
			if (iqseq(1) <= res(i)%start) then
				iRes = i
				exit
			end if
		end do
!		offset=res(iRes)%start-1
	end if

atom_res_id=0
	do i = 0,nAtoms-1
		ligand%mol%atom(i)%id				 = i +1
		ligand%mol%atom(i)%topindex		= iqseq(i+1)
		ligand%mol%atom(i)%name			 = ''
		ligand%mol%atom(i)%coor(0)	 = calc_xtop(3*(iqseq(i+1)-1) +1)
		ligand%mol%atom(i)%coor(1)	 = calc_xtop(3*(iqseq(i+1)-1) +2)
		ligand%mol%atom(i)%coor(2)	 = calc_xtop(3*(iqseq(i+1)-1) +3)
		ligand%mol%atom(i)%origin		 = 1	! This is a ligand atom
		ligand%mol%atom(i)%ttype		 = tac(iac(iqseq(i+1)))
		ligand%mol%atom(i)%weight		 = iaclib(iac(iqseq(i+1)))%mass
		ligand%mol%atom(i)%solv			 = 0
		ligand%mol%atom(i)%bfactor	 = 0
		ligand%mol%atom(i)%logp			 = 0
		ligand%mol%atom(i)%valid		 = 1
!		if(masks(1)%mask(i+1)) ligand%mol%atom(i)%mask = 1 !masks(1)%mask(i+1)
		ligand%mol%atom(i)%occupancy = 1.0			! default
		ligand%mol%atom(i)%neib(:)	 = 0
		ligand%mol%atom(i)%bond(:)	 = 0
		ligand%mol%atom(i)%root(:)	 = 0

			if(iRes<nres) then
				if(atom_res_id + res(iRes)%start >= res(iRes+1)%start) then
					iRes = iRes +1
					atom_res_id = 0
				end if
			end if
			
			atom_res_id = atom_res_id +1
			ligand%mol%atom(i)%atom_res_id = atom_res_id
			ligand%mol%atom(i)%residue = res(iRes)%name      !What is function of this?
			write(ligand%mol%atom(i)%res_id, '(I8)') iRes		! use internal read to convert integer -> character
	end do	

	do i = 0,nBonds-1
		ligand%mol%bond(i)%id = i +1			! meanwhile
		ligand%mol%bond(i)%valid  = 1
		ligand%mol%bond(i)%atom_1 = aB(i +1)%a%top_nr - aQ(1)%top_nr +1		! Remove offset so that first q-atom has index 1
		ligand%mol%bond(i)%atom_2 = aB(i +1)%b%top_nr	- aQ(1)%top_nr +1
		ligand%mol%bond(i)%ttype = SYBYL_bond_type(aB(i +1)%cod)
		
		if(adjustl(trim(ligand%mol%bond(i)%ttype)).eq.'') then
			write(*,'(a,i4,a)') '>>> WARNING: Ligand bond ', i+1, ' is of unknown type. Rings may not be correctly identified.'
			warn = warn +1
		end if
	end do

	! Bubble sort bonds in increasing atom order
	do
		mark = 0
		do i = 0, nBonds -2
			if((ligand%mol%bond(i)%atom_1.eq.ligand%mol%bond(i+1)%atom_1).and. &
				(ligand%mol%bond(i)%atom_2.eq.ligand%mol%bond(i+1)%atom_2)) goto 2
			if(ligand%mol%bond(i)%atom_1 <  ligand%mol%bond(i+1)%atom_1) goto 2
			if((ligand%mol%bond(i)%atom_1.eq.ligand%mol%bond(i+1)%atom_1) .and. &
				(ligand%mol%bond(i)%atom_2<=ligand%mol%bond(i+1)%atom_2)) goto 2

			! If no rule triggered, the bonds are in the wrong order
			mark = 1
			tmp = ligand%mol%bond(i)
			ligand%mol%bond(i) = ligand%mol%bond(i+1)
			ligand%mol%bond(i+1) = tmp
2		end do
		if(mark.eq.0) exit
	end do

	! Reset bond id
	do i = 0,nBonds -1
		ligand%mol%bond(i)%id = i +1
	end do

end subroutine xligand_translate

real(kind=prec) function Ligand_Atom_Buried_Surface(ligand,protein,id_in, total, buried)
	type(tXLigand)	:: ligand
	type(tXProtein)	:: protein
	integer	:: id,id_in
	real(kind=prec)		:: total
	real(kind=prec)		:: buried

	integer ::  j,k,num 
	integer ::  mark 
	real(kind=prec) ::  d,ratio 
	integer,allocatable :: atom_check_list(:)

	id = id_in

	allocate(atom_check_list(0:protein%num_atom-1),stat=err)
	if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'

	total=0
	buried=0.000 

	do j = 0,protein%num_atom -1 
	  if((protein%atom(j)%valid.ne.2)  .or. &
			 (protein%atom(j)%xtype.eq.'H').or. &
			 (protein%atom(j)%xtype.eq.'O.w')) then

			atom_check_list(j)=0
			goto 1
		end if

		d=Distance(ligand%mol%atom(id-1)%coor,protein%atom(j)%coor) 

		if(d>(ligand%mol%atom(id-1)%R+protein%atom(j)%R+2*WATER_R)) then
			atom_check_list(j)=0 
		else
			atom_check_list(j)=1 
		end if
1	end do

	num=ligand%mol%num_sur_dot
	do j = 1,num										! === NOTE: CHANGED BOUNDARIES FROM 0,num-1 TO 1,num
		if(ligand%mol%sur_dot(j)%valid.ne.id) goto 2

		! check if this dot is buried 
		mark=0 
		do k = 0,protein%num_atom -1 
	    if(atom_check_list(k).eq.0) goto 3

      d=Distance(ligand%mol%sur_dot(j)%coor(0:2), protein%atom(k)%coor) 
			
			if(d>(protein%atom(k)%R+WATER_R)) then
				goto 3
      else
				mark=1
				exit
			end if
3   end do

		total=total + ligand%mol%sur_dot(j)%unit 

 		if(mark.eq.0) then
			goto 2
		else 
			buried=buried + ligand%mol%sur_dot(j)%unit 
    end if
2	end do

	deallocate(atom_check_list) 

	if(total<=0.00) then
		ratio=0.00 
	else 
		ratio=buried/total  
	end if

 Ligand_Atom_Buried_Surface = ratio 
end function Ligand_Atom_Buried_Surface

! ==========================================================================================================================================
! ==========================================================================================================================================

! XFORECFIELD MODULE

! ==========================================================================================================================================
! ==========================================================================================================================================

subroutine ForceField_Read_SURFACE_DEF(filename)
	character(len=80)	:: filename
	integer						:: fp
	integer						:: i,count, filestat
	character(len=256):: buf,head 

	fp = freefile()
	open(unit=fp,file=filename,err=999,action='READ')
		! first, determine num_sdot_type 
    count=0 
		do
			read(fp,'(a256)',iostat=filestat) buf
			if(filestat.eq.0) then
				if(buf(1:1).eq.'#') goto 1		! skip REM
				if(trim(buf).eq.'') goto 1		! skip blank line

				if(buf(1:6).eq.'DOTSET') count = count +1 
			else
				exit
			end if
1		end do

    rewind(fp) 

    ff%num_sdot_type=count 

    allocate(ff%sdot(0:ff%num_sdot_type),stat=err)
		if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'
		ff%sdot(:)%max_dot = 0
		ff%sdot(:)%num_dot = 0

		! now read pre-calculated volume dots 
		count=0 
		do
			read(fp,'(a256)',iostat=filestat) buf
			if(filestat.eq.0) then
				if(buf(1:1).eq.'#') goto 2		! skip REM
				if(trim(buf).eq.'') goto 2		! skip blank line

				if(buf(1:3).eq.'END') exit 
				
				head = buf(1:6)
				if(head.eq.'DOTSET') then
					read(buf(7:256),*,iostat=filestat) ff%sdot(count)%r,ff%sdot(count)%num_dot,ff%sdot(count)%unit,ff%sdot(count)%total
					ff%sdot(count)%ttype = 'Un'

!					write(*,'(t1,a,t35,f3.1,t40,i4,t46,f5.2,t54,f6.2)') 'count, r, num_dot, unit, total = ', ff%sdot(count)%r,ff%sdot(count)%num_dot,ff%sdot(count)%unit,ff%sdot(count)%total
					allocate(ff%sdot(count)%dot(0:ff%sdot(count)%num_dot -1),stat=err)
					if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'

					do i = 0,ff%sdot(count)%num_dot -1 
						read(fp,'(a256)',iostat=filestat) buf
						read(buf(1:256),*,iostat=filestat) &
							ff%sdot(count)%dot(i)%coor(0), &
							ff%sdot(count)%dot(i)%coor(1), &
							ff%sdot(count)%dot(i)%coor(2)

							ff%sdot(count)%dot%valid=1 
							ff%sdot(count)%dot%unit=ff%sdot(count)%unit 
							ff%sdot(count)%dot%ttype = 'Un'
					end do

					count = count +1 
				end if
			end if
2		end do

	close(fp) 
	return

999 write(*,'(a,a,a)') '>>>>> ERROR: Unable to read ', adjustl(trim(filename)), '. Exiting.'
		stop
end subroutine ForceField_Read_SURFACE_DEF

type(tXDotSet) function ForceField_Get_Surface_Dot(atom, r)
	type(tXAtom)		:: atom
	real(kind=prec)						:: r
	real(kind=prec)						:: bigR 
	type(tXDotSet)	:: tmp_set 

	bigR=atom%R+r

	tmp_set=ForceField_Get_Surface_Dot_Rxyz(bigR,atom%coor(0),atom%coor(1),atom%coor(2)) 

	tmp_set%ttype = atom%ttype

	tmp_set%dot(0:tmp_set%num_dot -1)%valid = atom%id 
	tmp_set%dot(0:tmp_set%num_dot -1)%ttype = atom%ttype

	ForceField_Get_Surface_Dot = tmp_set
end function ForceField_Get_Surface_Dot

type(tXDotSet) function ForceField_Get_Surface_Dot_Rxyz(bigR, x, y, z)
	real(kind=prec)						:: bigR,x,y,z
	integer					:: i,j,num 
	integer					:: mark 
	real(kind=prec)						:: tmp,theta,phi,theta_step,phi_step 
	real(kind=prec)						:: r,d,total 
	type(tXDotSet)	:: tmp_set 
	type(tXDot)			:: tmp_dot 

	! check the pre-calculated surface dot sets
  mark=0

  do i = 0,ff%num_sdot_type -1
    if(abs(bigR-ff%sdot(i)%r)>0.025) then
			goto 1
    else
			!tmp_set=ff%sdot(i)
			call DotSet_operator_Copy(tmp_set,ff%sdot(i))

!		write(*,*) 'selected sdot ', i, tmp_set%num_dot;
			mark=1
			exit
		end if
1 end do

	if(mark.eq.1) then
		do i = 0,tmp_set%num_dot -1
			tmp_set%dot(i)%coor(0)=tmp_set%dot(i)%coor(0)+x
			tmp_set%dot(i)%coor(1)=tmp_set%dot(i)%coor(1)+y
			tmp_set%dot(i)%coor(2)=tmp_set%dot(i)%coor(2)+z
		end do
	 
		ForceField_Get_Surface_Dot_Rxyz = tmp_set
		return
	end if

	! if it is not pre-calculated, calculate it now
	! === NOTE: THERE IS A BUG IN THIS SECTION, PROGRAM WILL NOT EXECUTE. 
	! ===       SAME BUG FOUND IN ORIGINAL c++ SOURCE
	total=(4.000*PI*bigR*bigR) 

	! d=0.500 		// spacing between two dots
	d=sqrt(total/300) 
	if(d<0.500) d=0.500 

	tmp=int(PI*bigR/d+0.500)
	theta_step=PI/tmp

	num=0
	tmp_set%num_dot = 0

	!do theta = 0,PI -1, theta_step
	do i = 0,int(PI/theta_step)
		theta = theta_step*i

		r=bigR*sin(theta) 
		tmp=int(2*PI*r/d+0.500)
		phi_step=2*PI/tmp

!		do phi = 0,(2*PI) -1, phi_step
		do j = 0,int((2*PI)/phi_step)
			phi = phi_step*j
      tmp_dot%coor(0)=r*cos(phi) 
      tmp_dot%coor(1)=r*sin(phi) 
      tmp_dot%coor(2)=bigR*cos(theta) 
			tmp_dot%valid=1 
			tmp_dot%ttype = 'Un'
			call DotSet_PushBack(tmp_set,tmp_dot) 
      num = num +1 
    end do
  end do

	tmp_set%ttype = 'Un'
	tmp_set%r=bigR 
  tmp_set%num_dot=num
	tmp_set%unit=total/num
	tmp_set%total=total 

	do i = 0,tmp_set%num_dot -1  
		 tmp_set%dot(i)%unit=tmp_set%unit 
		 tmp_set%dot(i)%coor(0)=tmp_set%dot(i)%coor(0)+x 
		 tmp_set%dot(i)%coor(1)=tmp_set%dot(i)%coor(1)+y 
		 tmp_set%dot(i)%coor(2)=tmp_set%dot(i)%coor(2)+z 
	end do

	ForceField_Get_Surface_Dot_Rxyz = tmp_set 
end function ForceField_Get_Surface_Dot_Rxyz

subroutine ForceField_Read_RESIDUE_DEF(filename)
	character(len=256)		:: filename

	integer ::	fp						! file number
	integer ::  i,num,count 
	character(len=256) :: buf,head
	integer	:: filestat
	
	fp = freefile()
	open(unit=fp,file=filename,err=999,action='READ')

    ! first, determine num_restype
    count=0 
    do
			read(fp,'(a256)',iostat=filestat) buf
			if(filestat.eq.0) then
				if(buf(1:1).eq.'#') goto 1		! skip REM
				if(trim(buf).eq.'') goto 1		! skip blank line

				if(buf(1:4).eq.'RESI') count = count +1 
			else
				exit
			end if
1		end do      

    rewind(fp) 

		ff%num_restype=count 

		allocate(ff%residue(0:ff%num_restype -1),stat=err)
		if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'

		! second, read the residue templates
		count=0 

		do
			read(fp,'(a256)',iostat=filestat) buf
			if(filestat.eq.0) then
				if(buf(1:1).eq.'#') goto 2		! skip REM
				if(trim(buf).eq.'') goto 2		! skip blank line

				if(buf(1:3).eq.'END') exit 
				
				head = buf(1:4)
				if(head.eq.'RESI') then
					num=0 
					read(buf(5:256),*,iostat=filestat) ff%residue(count)%name,ff%residue(count)%num_atom,ff%residue(count)%num_bond

					! read the ATOM information first
					num = ff%residue(count)%num_atom 

					do i = 0,num -1 
						read(fp,'(a256)',iostat=filestat) buf
						read(buf(5:256),*,iostat=filestat) &
							ff%residue(count)%atom(i)%name,		&
							ff%residue(count)%atom(i)%ttype,	&
							ff%residue(count)%atom(i)%xtype,	&
							ff%residue(count)%atom(i)%r,			&
							ff%residue(count)%atom(i)%eps,		&
							ff%residue(count)%atom(i)%q,			&
							ff%residue(count)%atom(i)%hb,			&
							ff%residue(count)%atom(i)%logp,		&
							ff%residue(count)%atom(i)%solv,		&
							ff%residue(count)%atom(i)%ring,		&
							ff%residue(count)%atom(i)%pmftype 
						end do

					! then read the BOND information
					num=ff%residue(count)%num_bond 

					if(num>0) then
						do i = 0,num -1 
							read(fp,'(a256)',iostat=filestat) buf
							read(buf(5:256),*,iostat=filestat) &
							ff%residue(count)%bond(i)%atom1, &
							ff%residue(count)%bond(i)%atom2, &
							ff%residue(count)%bond(i)%ttype 
            end do
					end if
		
					count = count +1 
				end if
			end if
2		end do

	close(fp) 

	write(*,'(a,i4,a,a)') 'Read ', ff%num_restype, ' residue types from ', trim(filename)
  if(count.ne.ff%num_restype) write(*,*) 'Error occured when reading ', filename
	return

999 write(*,'(a,a,a)') '>>>>> ERROR: Unable to read ', adjustl(trim(filename)), '. Exiting.'
		stop
end subroutine ForceField_Read_RESIDUE_DEF

subroutine ForceField_Read_ATOM_DEF(filename)
	character(len=256)	:: filename
	integer							:: fp, filestat
	integer							:: count 
	character(len=256)	:: buf
	integer							:: slask
	
	fp = freefile()
	open(unit=fp,file=filename,err=999,action='READ')
	
		! first, determine num_atomtype
    count=0 
		do
			read(fp,'(a256)',iostat=filestat) buf
			if(filestat.eq.0) then
				if(buf(1:1).eq.'#') goto 1		! skip REM
				if(trim(buf).eq.'') goto 1		! skip blank line
				count = count +1 
			else
				exit
			end if
1		end do

		rewind(fp) 

		ff%num_atomtype=count 

		allocate(ff%atom(0:ff%num_atomtype -1),stat=err)
		if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'

		! second, read parameters for each atom type
    count=0 
    do
			read(fp,'(a256)',iostat=filestat) buf
			if(filestat.eq.0) then
				if(buf(1:1).eq.'#') goto 2		! skip REM
				if(trim(buf).eq.'') goto 2		! skip blank line

				read(buf,*,iostat=filestat) & 
					slask,									&
					ff%atom(count)%ttype,		&
					ff%atom(count)%weight,	&
					ff%atom(count)%r,				&
					ff%atom(count)%eps,			&
					ff%atom(count)%q,				&
					ff%atom(count)%hb				
				count = count +1 
			else
				exit
			end if
2		end do
	close(fp) 

	write(*,'(a,i4,a,a)') 'Read ', ff%num_atomtype, ' atom types from ', trim(filename)
	if(ff%num_atomtype.ne.count) write(*,*) 'Error reading ', filename
	return

999 write(*,'(a,a,a)') '>>>>> ERROR: Unable to read ', adjustl(trim(filename)), '. Exiting.'
		stop
end subroutine ForceField_Read_ATOM_DEF

subroutine ForceField_Read_XATOM_DEF(filename)
	character(len=256)	:: filename
	integer							:: fp, filestat
	integer							:: count 
	character(len=256)	:: buf
	integer							:: slask
	character(len=4)		:: slask2

	fp = freefile()
	open(unit=fp,file=filename,err=999,action='READ')
	
		! first, determine num_xatomtype
    count=0 
		do
			read(fp,'(a256)',iostat=filestat) buf
			if(filestat.eq.0) then
				if(buf(1:1).eq.'#') goto 1		! skip REM
				if(trim(buf).eq.'') goto 1		! skip blank line
				count = count +1 
			else
				exit
			end if
1		end do

		rewind(fp) 

		ff%num_xatomtype=count 

		allocate(ff%xatom(0:ff%num_xatomtype),stat=err)
		if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'

		! second, read parameters for each atom type
    count=0 
    do
			read(fp,'(a256)',iostat=filestat) buf
			if(filestat.eq.0) then
				if(buf(1:1).eq.'#') goto 2		! skip REM
				if(trim(buf).eq.'') goto 2		! skip blank line

				read(buf,*,iostat=filestat) & 
					slask,									&
					ff%xatom(count)%ttype,			&
					slask2,									&
			    ff%xatom(count)%logp 
				count = count +1 
			else
				exit
			end if
2		end do
	close(fp) 

	write(*,'(a,i4,a,a)') 'Read ', ff%num_xatomtype, ' x-atom types from ', trim(filename)
	if(ff%num_xatomtype.ne.count) write(*,*) 'Error reading ', filename
	return

999 write(*,'(a,a,a)') '>>>>> ERROR: Unable to read ', adjustl(trim(filename)), '. Exiting.'
		stop
end subroutine ForceField_Read_XATOM_DEF

real(kind=prec) function FF_Get_Atom_LogP(ttype)
	character(len=20):: ttype
	integer :: i
	integer ::mark
	real(kind=prec) ::logp

	mark=0
	do i = 0,ff%num_xatomtype -1
	 if(ttype.eq.ff%xatom(i)%ttype) then
	 	 logp=ff%xatom(i)%logp
		 mark=1
		 exit
	 end if
	end do

	if(mark.eq.1) then
		FF_Get_Atom_LogP = logp
	else 	
	  write(*,'(a,a,a)') '>>> WARNING: No LogP parameter for atom type ', ttype, '  zero value assigned'
	  warn = warn +1
		FF_Get_Atom_LogP = 0
	end if
end function FF_Get_Atom_LogP

integer function ForceField_Assign_Atom_Parameters(atm)
	type(tXAtom)	 :: atm
	integer ::  i 

	if(ff%num_atomtype.eq.0) write(*,*) 'Empty force field'
	do i = 0,ff%num_atomtype -1
    if(atm%xtype.ne.ff%atom(i)%ttype) then
			goto 1
    else
			atm%weight	= ff%atom(i)%weight 
			atm%r				= ff%atom(i)%r 
			atm%xr			= ff%atom(i)%r
			atm%eps			= ff%atom(i)%eps 
			atm%q				= ff%atom(i)%q
			atm%hb			= ff%atom(i)%hb
			atm%valid		= 1 
			ForceField_Assign_Atom_Parameters = 1
			return
    end if
1 end do

	write(*,*) '>>> WARNING: No parameter for atom ', atm%id, atm%ttype, atm%xtype
	warn = warn +1

	atm%valid=0
	ForceField_Assign_Atom_Parameters = 0
end function ForceField_Assign_Atom_Parameters

subroutine ForceField_translate

end subroutine ForceField_translate

integer function ForceField_Assign_Patom_Parameters(atm,suppress_warnings)
	type(tXAtom) :: atm
	integer,optional	:: suppress_warnings
  integer ::  i,j 

	! First, check the standard residue templates to map the atom 
	do i = 0,ff%num_restype -1 
		if(atm%residue.ne.ff%residue(i)%name) goto 1
		do j = 0,ff%residue(i)%num_atom -1  ! residue template found
      if(ff%residue(i)%atom(j)%name.ne.atm%name) goto 2

			atm%ttype = ff%residue(i)%atom(j)%ttype 
			atm%xtype = ff%residue(i)%atom(j)%xtype 
			atm%type2 = ff%residue(i)%atom(j)%pmftype 
			atm%r			= ff%residue(i)%atom(j)%r 
			atm%xr		= ff%residue(i)%atom(j)%r 
			atm%eps		= ff%residue(i)%atom(j)%eps 
			atm%q			= ff%residue(i)%atom(j)%q 
			atm%R			= ff%residue(i)%atom(j)%r 
			atm%hb		= ff%residue(i)%atom(j)%hb 
			atm%logp	= ff%residue(i)%atom(j)%logp 
			atm%solv	= ff%residue(i)%atom(j)%solv 
			atm%ring	= ff%residue(i)%atom(j)%ring 
			atm%valid	= 1
			ForceField_Assign_Patom_Parameters = 1
			return
2		end do
1	end do

	! check if they are terminal atoms
	do i = 0,ff%num_restype -1 
		if(ff%residue(i)%name.ne.'TER') goto 3
		do j = 0,ff%residue(i)%num_atom -1   ! residue template found
			if(ff%residue(i)%atom(j)%name.ne.atm%name) goto 4

			atm%ttype = ff%residue(i)%atom(j)%ttype 
			atm%xtype = ff%residue(i)%atom(j)%xtype 
			atm%type2 = ff%residue(i)%atom(j)%pmftype 
			atm%r			= ff%residue(i)%atom(j)%r 
			atm%eps		= ff%residue(i)%atom(j)%eps 
			atm%q			= ff%residue(i)%atom(j)%q 
			atm%R			= ff%residue(i)%atom(j)%r 
			atm%hb		= ff%residue(i)%atom(j)%hb 
			atm%logp	= ff%residue(i)%atom(j)%logp 
			atm%solv	= ff%residue(i)%atom(j)%solv 
			atm%ring	= ff%residue(i)%atom(j)%ring 
			atm%valid	= 1
			ForceField_Assign_Patom_Parameters = 1
			return
4		end do
3	end do

	! check if they are valid hetero atoms 
	do i = 0,ff%num_restype -1 
		if(ff%residue(i)%name.ne.'HET') goto 5
		do j = 0,ff%residue(i)%num_atom -1   ! residue template found
			if(ff%residue(i)%atom(j)%name.ne.atm%name) goto 6

			atm%ttype = ff%residue(i)%atom(j)%ttype 
			atm%xtype = ff%residue(i)%atom(j)%xtype 
			atm%type2 = ff%residue(i)%atom(j)%pmftype 
			atm%r			= ff%residue(i)%atom(j)%r 
			atm%eps		= ff%residue(i)%atom(j)%eps 
			atm%q			= ff%residue(i)%atom(j)%q 
			atm%R			= ff%residue(i)%atom(j)%r 
			atm%hb		= ff%residue(i)%atom(j)%hb 
			atm%logp	= ff%residue(i)%atom(j)%logp 
			atm%solv	= ff%residue(i)%atom(j)%solv 
			atm%ring	= ff%residue(i)%atom(j)%ring 
			atm%valid	= 1
			ForceField_Assign_Patom_Parameters = 1
			return
6		end do
5 end do

	if(.not. present(suppress_warnings) .and. atm%iscofactor .ne. 0)then
		 write(*,'(a,a,a,a,a,a,a,a,a)') &
		    'WARNING: No parameters for protein atom ', &
		    trim(atm%name), ' in residue nr. ', &
		    trim(adjustl(atm%res_id)), ' (', trim(atm%residue),'/', &
		    trim(atm%old_residue),'). This atom is ignored.'
		 warn = warn +1
	end if
	atm%valid=0
	ForceField_Assign_Patom_Parameters = 0
end function ForceField_Assign_Patom_Parameters

integer function ForceField_Patom_Connection_Test(atm1, atm2)
	type(tXAtom) :: atm1,atm2
	integer ::  i,j 

	if(atm1%chain.ne.atm2%chain) then
		ForceField_Patom_Connection_Test = 0
		return
	end if
	if((atm1%residue.ne.atm2%residue)) then
		ForceField_Patom_Connection_Test = 0
		return
	end if
	if((atm1%res_id.ne.atm2%res_id)) then
		ForceField_Patom_Connection_Test = 0
		return
	end if
	if(atm1%id.eq.atm2%id) then
		ForceField_Patom_Connection_Test = 0
		return
	end if


	do i = 0,ff%num_restype -1 
		if((ff%residue(i)%name.ne.atm1%residue)) then 
			goto 1
		elseif((ff%residue(i)%name.ne.atm2%residue)) then 
			goto 1
		end if

		if(ff%residue(i)%num_bond.eq.0) exit 

		do j = 0,ff%residue(i)%num_bond -1 
		 if((ff%residue(i)%bond(j)%atom1.eq.atm1%name).and.(ff%residue(i)%bond(j)%atom2.eq.atm2%name)) then
			ForceField_Patom_Connection_Test = 1
			return
		 elseif((ff%residue(i)%bond(j)%atom1.eq.atm2%name).and.(ff%residue(i)%bond(j)%atom2.eq.atm1%name)) then
			ForceField_Patom_Connection_Test = 1
			return
     else 
			 goto 2
		 end if
2		end do
	 exit 
1 end do

	ForceField_Patom_Connection_Test = 0
end function ForceField_Patom_Connection_Test

! ==========================================================================================================================================
! ==========================================================================================================================================

! XMISC MODULE

! ==========================================================================================================================================
! ==========================================================================================================================================

subroutine Ring_Copy(ring,r)
	type(tXRing)	:: ring
	type(tXRing)	:: r

	integer				:: na,nb

	call Ring_Clear(ring)

	na = size(r%atom_id)
	nb = size(r%bond_id)

	ring%valid				= r%valid
	ring%ttype				= r%ttype
	ring%num_member		= r%num_member

	if(na>0) then	
		allocate(ring%atom_id(0:na-1),stat=err)
		if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'
		ring%atom_id(:) = r%atom_id(:)
	end if
	ring%atom_count	= r%atom_count
	ring%atom_max		= r%atom_max

	if(nb>0) then
		allocate(ring%bond_id(0:nb-1),stat=err)
		if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'
		ring%bond_id(:) = r%bond_id(:)
	end if
	ring%bond_count	= r%bond_count
	ring%bond_max		= r%bond_max

	ring%centroid(0:3) = r%centroid(0:3)
end subroutine Ring_Copy

subroutine Bond_Swap(b1,b2)
	type(tXBond)	:: b1,b2,tmp


	call Bond_Copy(tmp,b1)
	call Bond_Copy(b1,b2)
	call Bond_Copy(b2,tmp)
end subroutine 

subroutine Bond_Copy(b1,b2)
	type(tXBond)	:: b1,b2

	b1%id			= b2%id
	b1%valid	= b2%valid
	b1%atom_1 = b2%atom_1
	b1%atom_2 = b2%atom_2
	b1%ttype	= b2%ttype
	b1%part		= b2%part
	b1%ring		= b2%ring
	b1%length	= b2%length
	b1%num_neib		= b2%num_neib
	b1%neib(0:5)	= b2%neib(0:5)
end subroutine 

real(kind=prec) function Angle(a,b,c)	
	! the a-b-c angle in degree 0-180
	real(kind=prec),dimension(0:2) :: a,b,c
	integer :: i
	real(kind=prec),dimension(0:2)	:: v1,v2

	do i=0,2
		v1(i)=b(i)-a(i)
		v2(i)=b(i)-c(i)
	end do

	Angle=Angle_of_Two_Vectors(v1,v2)
end function Angle

!integer function HBond_Value_HBond(hbond)
!	type(tXHBond)	:: hbond
!	real(kind=prec) ::  d,d0,d1,d2,d3,d4 
!	real(kind=prec) ::  a0,a1,a2,angle1,angle2,angle3,angle4 
!	real(kind=prec) ::  tmp1,tmp2,tmp3,tmp4 
!	integer ::  mark0,mark1,mark2 
!	real(kind=prec) ::  a,b,c,bigA,bigC 
!
!	! first, calculate the necessary parameters
!	d=Distance(hbond%D%coor,hbond%A%coor) 
!
!	! angle D-H-bigA 
!	if((hbond%D%hb.eq.'M')) then
!		 a0=0.1
!		 mark0=0
!	elseif((hbond%D%ttype.eq.'O.3')) then
!	 	 ! align H to the optimal position and then recalculate <D-H-bigA 
!	 	 b=0.98  ! D-H length 
!	 	 c=d 	 ! D-bigA length 
!	 	 bigA=abs(Angle(hbond%D%root,hbond%D%coor,hbond%A%coor)-109.0)  
!	 	 a=sqrt(b*b+c*c-2*b*c*cos(bigA/180.0)) 	! H-bigA length
!		 bigC=acos((a*a+b*b-c*c)/(2*a*b))/PI*180.0 
!	 	 a0=180.0-bigC 
!		 mark0=1 
!	else 
!		 a0=180.0-Angle(hbond%D%coor,hbond%H%coor,hbond%A%coor) 
!		 mark0=1 
!	end if
!
!	! angle DR-D-bigA
!	if((hbond%D%hb.eq.'M').or.(hbond%D%ttype.eq.'O.w')) then
!		 a1=0.1
!		 mark1=0 
!	else 
!		 a1=180.0-Angle(hbond%D%root,hbond%D%coor,hbond%A%coor) 
!		 mark1=1 
!	end if
!
!	! D-bigA-bigAR
!	if((hbond%A%hb.eq.'M').or.(hbond%A%ttype.eq.'O.w'))  then
!		 a2=0.1
!		 mark2=0 
!	else 
!		 a2=180.0-Angle(hbond%D%coor,hbond%A%coor,hbond%A%root) 
!		 mark2=1 
!	end if
!
!	! now determine the geometry type of this h-bond
!	call HBond_Determine_DbigA_Type(hbond) 
!
!	! now compute the strength of this h-bond
!	d0=hbond%D%R+hbond%A%R  
!	d1=0.00
!	d2=1.00  
!	d3=d0-0.40
!	d4=d0+0.10
!
!	if(d<d1) then
!		tmp1=0.000 
!	elseif(d<d2) then
!		tmp1=(d-d1)/(d2-d1) 
!	elseif(d<d3) then
!		tmp1=1.000 
!	elseif(d<d4) then	
!		tmp1=(d4-d)/(d4-d3) 
!	else 
!		tmp1=0.000 
!	end if
!
!	if(mark0.eq.1) then
!		 angle1=0.0  
!		 angle2=0.001  
!		 angle3=60.0
!		 angle4=90.0
!		 if(a0<angle1) then
!			tmp2=0.000 
!     elseif(a0<angle2) then
!			tmp2=(a0-angle1)/(angle2-angle1) 
!     elseif(a0<angle3) then
!			tmp2=1.000 
!     elseif(a0<angle4) then
!			tmp2=(angle4-a0)/(angle4-angle3) 
!     else 
!			tmp2=0.000 
!		 end if
!	else 
!		tmp2=1.000 
!	end if
!
!	if(mark1.eq.1)  then
!		 if(d_type.eq.1) then
!			angle1=0.0  
!			angle2=0.001  
!			angle3=30.0
!			angle4=60.0
!		 else 
!			angle1=20.0  
!			angle2=40.0  
!			angle3=70.0
!			angle4=90.0
!		 end if
!
!     if(a1<angle1) then
!			tmp3=0.000 
!     elseif(a1<angle2) then	
!			tmp3=(a1-angle1)/(angle2-angle1) 
!     elseif(a1<angle3) then	
!			tmp3=1.000 
!     elseif(a1<angle4) then
!			tmp3=(angle4-a1)/(angle4-angle3) 
!     else 
!			tmp3=0.000
!		 end if
!	else 
!		tmp3=1.000 
!
!	if(mark2.eq.1) then
!		 if(a_type.eq.1) then
!			angle1=0.0
!			angle2=0.001  
!			angle3=30.0
!			angle4=60.0
!     else 
!			angle1=0.0  
!			angle2=5.0  
!			angle3=60.0
!			angle4=90.0
!		 end if
!
!     if(a2<angle1) then
!			tmp4=0.000 
!     elseif(a2<angle2) then
!			tmp4=(a2-angle1)/(angle2-angle1) 
!     elseif(a2<angle3) then
!			tmp4=1.000 
!     elseif(a2<angle4) then
!			tmp4=(angle4-a2)/(angle4-angle3) 
!     else 
!			tmp4=0.000 
!		 end if
!	else 
!		tmp4=1.000 
!	end if
!
!	hbond%d=d  
!	hbond%a0=a0  
!	hbond%a1=a1
!	hbond%a2=a2
!
!	hbond%score=tmp1*tmp2*tmp3*tmp4 
!
!	if(hbond%score<0.001) then
!		HBond_Value_HBond = 0
!	else 
!		HBond_Value_HBond = e
!	end if
!end function HBond_Value_HBond

integer function HBond_Value_HBond_2(hbond)
	type(tXHBond)	:: hbond
	real(kind=prec) ::  d,d0,d1,d2,d3,d4 
	real(kind=prec) ::  a1,a2,angle1,angle2,angle3,angle4 
	real(kind=prec) ::  tmp1,tmp3,tmp4 
	integer ::  mark1,mark2 

	! first, calculate the necessary parameters
	d=Distance(hbond%D%coor,hbond%A%coor) 

	! angle DR-D-A
	if((hbond%D%hb.eq.'M').or.(hbond%D%ttype.eq.'O.w')) then
		 mark1=0
		 a1=0.0 
	else 
		 a1=180.0-Angle(hbond%D%root,hbond%D%coor,hbond%A%coor) 
		 mark1=1 
	end if

	! D-A-AR
	if((hbond%A%hb.eq.'M').or.(hbond%A%ttype.eq.'O.w')) then
		 mark2=0
		 a2=0.0 
	else 
		 a2=180.0-Angle(hbond%D%coor,hbond%A%coor,hbond%A%root) 
		 mark2=1 
	end if

	! now determine the geometry type of this h-bond
	call HBond_Determine_DA_Type(hbond) 

	! now compute the strength of this h-bond

	d0=hbond%D%R+hbond%A%R  
	d1=0.00  
	d2=1.00  
	d3=d0-0.40
	d4=d0+0.20

	if(d<d1) then
		tmp1=0.000 
	elseif(d<d2) then	
		tmp1=(d-d1)/(d2-d1) 
	elseif(d<d3) then	
		tmp1=1.000 
	elseif(d<d4) then	
		tmp1=(d4-d)/(d4-d3) 
	else 
		tmp1=0.000 
	end if

	if(mark1.eq.1) then
		 if(hbond%d_type.eq.1) then
			angle1=0.0  
			angle2=0.001  
			angle3=25.0
			angle4=50.0
		 else 
			angle1=25.0  
			angle2=50.0  
			angle3=75.0
			angle4=100.0
		 end if
     if(a1<angle1) then
			tmp3=0.000 
     elseif(a1<angle2) then
			tmp3=(a1-angle1)/(angle2-angle1) 
     elseif(a1<angle3) then
			tmp3=1.000 
     elseif(a1<angle4) then
			 tmp3=(angle4-a1)/(angle4-angle3) 
     else 
			tmp3=0.000 
		 end if	
	else 
		tmp3=1.000 
	end if

	if(mark2.eq.1) then
		 if(hbond%a_type.eq.1) then
       angle1=0.0  
			 angle2=0.001  
			 angle3=30.0
			 angle4=55.0
     else 
			 angle1=0.0  
			 angle2=5.0  
			 angle3=70.0
			 angle4=95.0
		 end if

     if(a2<angle1) then
			tmp4=0.000 
     elseif(a2<angle2) then
			tmp4=(a2-angle1)/(angle2-angle1) 
     elseif(a2<angle3) then
			tmp4=1.000 
     elseif(a2<angle4) then
			tmp4=(angle4-a2)/(angle4-angle3) 
     else 
			tmp4=0.000 
		 end if
	else 
		tmp4=1.000 
	endif

	hbond%dd=d  
	hbond%a1=a1
	hbond%a2=a2

	if(tmp3>=tmp4) then
		hbond%score=tmp1*tmp4 
	else 
		hbond%score=tmp1*tmp3 
	end if

	if(hbond%score<0.001) then
		HBond_Value_HBond_2 = 0 
	else 
		HBond_Value_HBond_2 = 1 
	end if
end function HBond_Value_HBond_2

subroutine HBond_Determine_DA_Type(hbond)
	type(tXHBond)	:: hbond
	
	hbond%d_type=Atom_Get_Donor_Type(hbond%D) 
	hbond%a_type=Atom_Get_Acceptor_Type(hbond%A) 
end subroutine HBond_Determine_DA_Type

integer function Atom_Get_Donor_Type(atom)
	type(tXAtom)	:: atom
	integer ::  ttype  

	if((atom%hb.ne.'D').and.(atom%hb.ne.'DA').and.(atom%hb.ne.'M')) then
		Atom_Get_Donor_Type= 0 
		return
	end if

 if(atom%origin.eq.2)  then ! protein atoms
 	 if((atom%hb.eq.'M')) then
		ttype=1 
 	 elseif((atom%ttype.eq.'O.3').or.(atom%ttype.eq.'O.co2').or.(atom%ttype.eq.'O.w')) then
			ttype=2 
 	 elseif((index(atom%xtype,'N.4').ne.0).or.(index(atom%xtype,'N.3').ne.0)) then
		ttype=2 
 	 elseif((index(atom%xtype,'N.pl3').ne.0).or.(index(atom%xtype,'N.2').ne.0).or.(index(atom%xtype,'N.ar').ne.0)) then
	 	 if((index(atom%name,'NH').ne.0).and.(index(atom%residue,'ARG').ne.0)) then
			ttype=2 
	 	 elseif((index(atom%name,'ND').ne.0).and.(index(atom%residue,'ASN').ne.0)) then
			ttype=2 
	 	 elseif((index(atom%name,'NE').ne.0).and.(index(atom%residue,'GLN').ne.0)) then
			ttype=2 
     else 
			ttype=1 
		 end if
 	 elseif((atom%ttype.eq.'S.3')) then
			ttype=2 
 	 else 
		ttype=2 
	 end if
 else  ! ligand atoms
 	 if((atom%ttype.eq.'O.3').or.(atom%ttype.eq.'O.co2')) then
		ttype=2 
   elseif((index(atom%xtype,'N.4').ne.0).or.(index(atom%xtype,'N.3').ne.0)) then
		 if(atom%num_nonh<=2) then
			ttype=2 
		 else 
			ttype=1 
		 end if
 	 elseif((index(atom%xtype,'N.pl3').ne.0).or.(index(atom%xtype,'N.2').ne.0).or.(index(atom%xtype,'N.ar').ne.0)) then
    if(atom%num_nonh<=1) then
			ttype=2 
    else 
			ttype=1 
		end if
 	 elseif((atom%ttype.eq.'S.3')) then
		ttype=2 
 	 else 
		ttype=2 
	 end if
 end if

 Atom_Get_Donor_Type =ttype 
end function Atom_Get_Donor_Type

integer function Atom_Get_Acceptor_Type(atom)
	type(tXAtom)	:: atom
  integer ::  ttype 

 if((atom%hb.ne.'A').and.(atom%hb.ne.'DA')) then
	Atom_Get_Acceptor_Type = 0 
	return
 end if


 if(atom%origin.eq.2)  then ! protein atoms
 	 if((atom%ttype.eq.'O.3').or.(atom%ttype.eq.'O.2').or.(atom%ttype.eq.'O.co2').or.(atom%ttype.eq.'O.w')) then
  	 ttype=2 
 	 elseif((index(atom%xtype,'N.pl3').ne.0).or.(index(atom%xtype,'N.2').ne.0).or. (index(atom%xtype,'N.ar').ne.0)) then
	   ttype=1 
 	 else 
		ttype=2 
	 endif
 else
 	 if((atom%ttype.eq.'O.3').or.(atom%ttype.eq.'O.2').or.(atom%ttype.eq.'O.co2')) then
   	 ttype=2 
 	 elseif((index(atom%xtype,'N.pl3').ne.0).or.(index(atom%xtype,'N.2').ne.0).or.(index(atom%xtype,'N.ar').ne.0)) then
		 if(atom%num_nonh<=1) then
			ttype=2 
     else 
			ttype=1 
		 end if
 	 else 
		ttype=2 
	 end if
	end if

 Atom_Get_Acceptor_Type = ttype 
end function Atom_Get_Acceptor_Type


real(kind=prec) function Angle_of_Two_Vectors(v1, v2)
	real(kind=prec), dimension(0:2)	:: v1,v2
  real(kind=prec) angle
  real(kind=prec) l1,l2,tmp1,tmp2

	if((v1(0).eq.v2(0)).and.(v1(1).eq.v2(1)).and.(v1(2).eq.v2(2))) then					! need to treat this case
		Angle_of_Two_Vectors = 0
		return
	end if

	l1=sqrt(v1(0)*v1(0)+v1(1)*v1(1)+v1(2)*v1(2))
	l2=sqrt(v2(0)*v2(0)+v2(1)*v2(1)+v2(2)*v2(2))

	tmp1=v1(0)*v2(0)+v1(1)*v2(1)+v1(2)*v2(2)
	tmp2=l1*l2

	angle=acos(tmp1/tmp2)
	angle=angle/(PI)*180.0

	Angle_of_Two_Vectors = angle  ! return angle in degree, 0-180
end function Angle_of_Two_Vectors

subroutine HBond_Clear(hbond)
	type(tXHBond) :: hbond

	hbond%valid= 0
	hbond%latom= 0
	hbond%patom= 0
	hbond%ttype= 0
	hbond%d_type= 0
	hbond%a_type=0 
	hbond%score=0.000
	hbond%sb=0
	hbond%dd=0
	hbond%a0=0
	hbond%a1=0
	hbond%a2=0.000
	call taClear(hbond%D)
	call taClear(hbond%H)
	call taClear(hbond%A)
end subroutine HBond_Clear

subroutine Dot_PushBack(dots,numdot,newdot)
	type(tXDot),pointer	:: dots(:)
	integer			:: numdot
	type(tXDot)	:: newdot
	type(tXDot),allocatable	:: tmpdot(:)

	if(.not. associated(dots)) then			! need to init dot array
		allocate(dots(16),stat=err)
		if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'
	end if		

	numdot = numdot +1
	if(numdot>size(dots,1)) then								! need to increase the size of dots
		allocate(tmpdot(numdot-1),stat=err)				! allocate tmp array of old size
		if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'
		tmpdot(:) = dots(:)												! copy dots to tmp

		deallocate(dots)													! deallocate dots
		allocate(dots(numdot +16),stat=err)				! allocate new and bigger mem
		if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'

		dots(1:numdot -1) = tmpdot(1:numdot -1)		! copy old part from tmp
	end if

	dots(numdot) = newdot
end subroutine Dot_PushBack

subroutine Int_PushBack(array, int)
	integer, pointer	:: array(:)
	integer						:: int
	integer, pointer	:: tmp(:)
	integer						:: s

	if(associated(array)) then
		s = size(array)
		allocate(tmp(0:s),stat=err)
		if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'
		tmp(0:s -1) = array(0:s-1)

		deallocate(array)
		tmp(s) = int

		array => tmp
	else
		allocate(array(0:0),stat=err)
		if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'
		array(0) = int
	end if
end subroutine Int_PushBack		

subroutine Ring_PushBack(ring,atom_id)
	type(tXRing)		:: ring
	integer					:: atom_id
	integer,allocatable	:: tmp(:)

	if(ring%atom_max.eq.0) then				! need to init dot array
		allocate(ring%atom_id(0:15),stat=err)
		if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'
		ring%atom_max = 16
	end if		
	
	ring%atom_count = ring%atom_count +1
	if(ring%atom_count>ring%atom_max) then
		allocate(tmp(1:ring%atom_count-1),stat=err)
		if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'
		tmp(:) = ring%atom_id(:)
		deallocate(ring%atom_id)
		ring%atom_max = ring%atom_max +16						! increase by 16 elements
		allocate(ring%atom_id(1:ring%atom_max),stat=err)
		if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'
		ring%atom_id(1:ring%atom_count -1) = tmp(1:ring%atom_count -1)
		deallocate(tmp)
	end if

	ring%atom_id(ring%atom_count -1) = atom_id
end subroutine Ring_PushBack

subroutine DotSet_PushBack(dotset,dot)
	type(tXDotSet)	:: dotset
	type(tXDot)			:: dot
	type(tXDot),allocatable	:: tmpdot(:)

	if(dotset%max_dot.eq.0) then			! need to init dot array
		allocate(dotset%dot(0:15),stat=err)
		if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'
		dotset%max_dot = 16
	end if		
	
	dotset%num_dot = dotset%num_dot +1
	if(dotset%num_dot>dotset%max_dot) then
		allocate(tmpdot(1:dotset%num_dot),stat=err)
		if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'
		tmpdot(1:dotset%num_dot -1) = dotset%dot(1:dotset%num_dot -1)
		deallocate(dotset%dot)
		dotset%max_dot = dotset%max_dot +16						! increase by 16 elements
		allocate(dotset%dot(1:dotset%max_dot),stat=err)
		if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'
		dotset%dot(1:dotset%num_dot -1) = tmpdot(1:dotset%num_dot -1)
	end if

	dotset%dot(dotset%num_dot -1) = dot
end subroutine DotSet_PushBack

subroutine DotSet_operator_Copy(s1,s2)
	type(tXDotSet)	:: s1,s2
	integer					:: i

	s1%num_dot = s2%num_dot
	s1%max_dot = s2%max_dot
	s1%r			 = s2%r
	s1%ttype(:)= s2%ttype(:)
	s1%unit		 = s2%unit
	s1%total	 = s2%total

	if(associated(s2%dot)) then
		allocate(s1%dot(0:size(s2%dot)-1),stat=err)
		if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'
		do i = 0,size(s2%dot) -1
			s1%dot(i)%valid		= s2%dot(i)%valid
			s1%dot(i)%ttype		= s2%dot(i)%ttype
			s1%dot(i)%coor(:)	= s2%dot(i)%coor(:)
			s1%dot(i)%unit		= s2%dot(i)%unit
			s1%dot(i)%score		= s2%dot(i)%score
		end do
	end if
end subroutine DotSet_operator_Copy

subroutine HBond_operator_Copy(b1,b2)
	type(tXHBond)	:: b1,b2

	b1%valid	= b2%valid
	b1%ttype	= b2%ttype
	b1%d_type	=	b2%d_type
	b1%a_type	=	b2%a_type
	b1%sb			=	b2%sb
	b1%latom	=	b2%latom
	b1%patom	=	b2%patom

	call Atom_operator_Copy(b1%H,b2%H)
	call Atom_operator_Copy(b1%D,b2%D)
	call Atom_operator_Copy(b1%A,b2%A)

	b1%dd			=	b2%dd
	b1%a0			=	b2%a0
	b1%a1			=	b2%a1
	b1%a2			= b2%a2

	b1%score	= b2%score
end subroutine HBond_operator_copy

subroutine Atom_operator_copy(a1,a2)
	type(tXAtom)	:: a1,a2

	a1%id						= a2%id
	a1%valid				= a2%valid
	a1%name(:)			= a2%name(:)
	a1%ttype(:)			= a2%ttype(:)
	a1%xtype(:)			= a2%xtype(:)
	a1%type2(:)			= a2%type2(:)
	a1%residue(:)		= a2%residue(:)
	a1%res_id(:)		= a2%res_id(:)
	a1%atom_res_id	= a2%atom_res_id
	a1%iscofactor		= a2%iscofactor
	a1%chain(:)			= a2%chain(:)
	a1%coor(0:2)		= a2%coor(0:2)
	a1%root(0:2)		= a2%root(0:2)
	a1%weight				= a2%weight
	a1%r						= a2%r
	a1%eps					= a2%eps
	a1%q						= a2%q
	a1%XR						= a2%XR
	a1%logp					= a2%logp
	a1%solv					= a2%solv
	a1%hb(:)				= a2%hb(:)
	a1%occupancy		= a2%occupancy
	a1%bfactor			= a2%bfactor
	a1%score				= a2%score
	a1%ring					= a2%ring
	a1%origin				= a2%origin
	a1%part					= a2%part
	a1%num_neib			= a2%num_neib
	a1%num_nonh			= a2%num_nonh
	a1%neib(0:6)		= a2%neib(0:6)
	a1%bond(0:6)		= a2%bond(0:6)
	a1%temp					= a2%temp
end subroutine Atom_operator_copy

subroutine Group_operator_Copy(g1,g2)
	type(tXGroup) :: g1,g2

	g1%valid =			g2%valid
	g1%num_neib =		g2%num_neib    	
	g1%num_nonh =		g2%num_nonh    	
	g1%num_h =			g2%num_h       	
	g1%num_hetero = g2%num_hetero  	
	g1%num_pi =			g2%num_pi      	
	g1%num_car =		g2%num_car     	
	g1%num_nar =		g2%num_nar     	
	g1%num_db =			g2%num_db	
	g1%num_tb =			g2%num_tb	
	g1%db_type =		g2%db_type 
	g1%amide =			g2%amide	
	g1%center  =		g2%center
	g1%neib(:) =		g2%neib(:)
	g1%bond(:) =		g2%bond(:)
end subroutine Group_operator_Copy

subroutine DotSet_Show_Contents(dotset)
	type(tXDotSet)	:: dotset
	integer					:: i

	write(*,*) 'DOTSET: ', dotset%r, dotset%num_dot, dotset%unit, dotset%total, dotset%ttype
	do i = 0, dotset%num_dot -1
		write(*,*) dotset%dot(i)%coor(0),dotset%dot(i)%coor(1),dotset%dot(i)%coor(2),dotset%dot(i)%valid
	end do
end subroutine DotSet_Show_Contents

subroutine taClear(this)
	type(tXAtom),intent(out)	:: this
		
	this%id=0
	this%valid=0
	this%part=1
	this%origin=1

	this%name = 'Un' 
	this%ttype = 'Un' 
	this%xtype = 'Un' 
	this%type2 = 'Un' 
	this%residue = 'Un' 
	this%res_id = '0' 
	this%chain=' '

	this%coor(0)=0
	this%coor(1)=0
	this%coor(2)=0.000
	this%root(0)=0
	this%root(1)=0
	this%root(2)=0.000

	this%r=0.000 
	this%eps=0.000 
	this%q=0.000 
	this%weight=0.000 
	this%R=0.000 
	this%hb = 'N' 
	this%logp=0.000 
	this%solv=0.000 
	this%score=0.000
	this%ring=0

	this%occupancy=1.000 
	this%bfactor=0.000

	this%num_neib=0
	this%num_nonh=0 

	this%neib(:) = 0
	this%bond(:) = 0
	
	this%temp=0
end subroutine taClear

subroutine Ring_Clear(ring,size)
	type(tXRing) :: ring
	integer,optional :: size  ! determines size of atom_id and bond_id
	
	if(associated(ring%atom_id)) deallocate(ring%atom_id)
	if(associated(ring%bond_id)) deallocate(ring%bond_id)

	if(present(size)) then
		allocate(ring%atom_id(0:size),stat=err)
		if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'
		allocate(ring%bond_id(0:size),stat=err)
		if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'
	end if

	ring%valid=0
	ring%ttype=0
	ring%num_member=0

	ring%centroid(0)=0
	ring%centroid(1)=0
	ring%centroid(2)=0.000
end subroutine Ring_Clear

real(kind=prec) function Distance(a, b)
	real(kind=prec),dimension(0:2)	:: a,b
	real(kind=prec)								:: d,tmpx,tmpy,tmpz

	tmpx = (a(0)-b(0))*(a(0)-b(0))
	tmpy = (a(1)-b(1))*(a(1)-b(1))
	tmpz = (a(2)-b(2))*(a(2)-b(2))

	d = sqrt(tmpx+tmpy+tmpz)

	Distance = d
end function Distance

! ================================================================================================================================================
! ================================================================================================================================================

!	XPROTEIN MODULE

! ================================================================================================================================================
! ================================================================================================================================================

subroutine Protein_Merge_Cofactor(protein,cf)
	type(tXProtein)	:: protein
	type(tXLigand)	:: cf
	integer					:: i,j
	integer,dimension(0:6)		:: tmp

	do i = 0,cf%mol%num_atom -1
		tmp(:) = protein%atom(cf%mol%atom(i)%temp)%neib(:)
		call Atom_operator_copy(protein%atom(cf%mol%atom(i)%temp), cf%mol%atom(i))
		protein%atom(cf%mol%atom(i)%temp)%id = cf%mol%atom(i)%temp 
		protein%atom(cf%mol%atom(i)%temp)%neib(:) = tmp(:)
		protein%atom(cf%mol%atom(i)%temp)%valid = 1
		protein%atom(cf%mol%atom(i)%temp)%bond(:) = 0
	end do
end subroutine Protein_Merge_Cofactor

subroutine xprotein_translate(protein,nAtoms,coordinates,num_bond,bond)
	type(tXProtein)		:: protein
	integer						:: nAtoms				! nat_solute being sent in
	real(kind=prec)						:: coordinates(:)
	integer						:: num_bond
	type(BOND_TYPE)		:: bond(:)
	integer						:: i,j,mark,offset
	type(tXBond)			:: tmp
	integer						:: iRes					! current residue
	integer						:: iNextRes			! starting point of next residue
	integer						:: atom_res_id	! atom number within residue
	integer						:: last_valid_atom

	! First find out how many atoms to translate
	iRes = 0
	iNextRes = 1		 
	j = 0
	do i = 1,nAtoms
		if(i.eq.iNextRes) then
			iRes = iRes +1
			if (iRes<nres) 	iNextRes = res(iRes +1)%start
			atom_res_id = 0
		end if
		
		if(iqatom(i).ne.0) cycle																	! skip ligand atoms
		if(res(iRes)%name.eq.'DUM') cycle													! skip dummy atoms
		if((tac(iac(i)).eq.'Du').or.(tac(iac(i)).eq.'DU')) cycle	! skip dummy
		j = j +1
	end do

	protein%num_atom = j

	! Translate atom coordinates and other stuff
	allocate(protein%atom(0:protein%num_atom -1),stat=err) 
	if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'
	protein%atom(:)%valid = 0
	iRes = 0
	iNextRes = 1		 
	j = -1
	do i = 1,nAtoms
		if(i.eq.iNextRes) then
			iRes = iRes +1
			if (iRes<nres) iNextRes = res(iRes +1)%start
			atom_res_id = 0
		end if
		
		if(iqatom(i).ne.0) cycle																	! skip ligand atoms
		if(res(iRes)%name.eq.'DUM') cycle													! skip dummy atoms
		if((tac(iac(i)).eq.'Du').or.(tac(iac(i)).eq.'DU')) cycle	! skip dummy

		last_valid_atom = i
		atom_res_id = atom_res_id +1
		j = j +1

		top2prot(i) = j																						! set topology -> protein translation

		protein%atom(j)%topindex		= i								! save topology index of this atom
		protein%atom(j)%coor(0) = coordinates(3*(i-1) +1)
		protein%atom(j)%coor(1) = coordinates(3*(i-1) +2)
		protein%atom(j)%coor(2) = coordinates(3*(i-1) +3)
		protein%atom(j)%id      = j
		protein%atom(j)%atom_res_id = atom_res_id
		protein%atom(j)%name		 = ''		! not used. names are assigned using lib file and index within residue
		protein%atom(j)%residue = res(iRes)%name
		write(protein%atom(j)%res_id, '(I8)') iRes		! use internal read to convert integer -> character
		protein%atom(j)%ttype	 = tac(iac(i))									! ttype and xtype are set by Value_Atom
		protein%atom(j)%xtype	 = 'Un' !protein%atom(i -1)%ttype
		protein%atom(j)%type2	 = 'Un' 
		protein%atom(j)%hb			 = 'N'
		protein%atom(j)%valid	 = 1	! default
!		if(masks(1)%mask(i+1)) protein%atom(i)%mask = 1 !masks(1)%mask(i+1)
		protein%atom(j)%origin	 = 2		! This is a protein atom
		protein%atom(j)%part		 = 1		! regular atom
		protein%atom(j)%occupancy = 0.0			! default
	end do


	! Count protein bonds
	protein%num_bond = 0
	do i = 1,num_bond
		if(((bond(i)%i<=nAtoms).and.(bond(i)%j<=nAtoms)).and. &
		(iqatom(bond(i)%i).eq.0.and.iqatom(bond(i)%j).eq.0)) protein%num_bond = protein%num_bond +1
	end do
	
	allocate(protein%bond(0:protein%num_bond-1))
	j = 0
	offset = 0
	do i = 1, num_bond
		if(((bond(i)%i<=nAtoms).and.(bond(i)%j<=nAtoms)).and.(iqatom(bond(i)%i).eq.0.and.iqatom(bond(i)%j).eq.0)) then
			protein%bond(j)%id = j +1
			protein%bond(j)%valid  = 1
			protein%bond(j)%atom_1 = top2prot(bond(i)%i) +1
			protein%bond(j)%atom_2 = top2prot(bond(i)%j) +1
			protein%bond(j)%ttype = SYBYL_bond_type(bond(i)%cod)
!			write(*,*) i, '  ', protein%atom(protein%bond(j)%atom_1-1)%ttype, '  ', protein%atom(protein%bond(j)%atom_2-1)%ttype, '  ', protein%bond(j)%ttype
			j = j +1
		else
			offset = offset +1
		end if
	end do

	! Bubble sort bonds in increasing atom order
	do
		mark = 0
		do i = 0, protein%num_bond -2
			if((protein%bond(i)%atom_1.eq.protein%bond(i+1)%atom_1).and. &
				(protein%bond(i)%atom_2.eq.protein%bond(i+1)%atom_2)) goto 2
			if(protein%bond(i)%atom_1 <  protein%bond(i+1)%atom_1) goto 2
			if((protein%bond(i)%atom_1.eq.protein%bond(i+1)%atom_1) .and. &
				(protein%bond(i)%atom_2<=protein%bond(i+1)%atom_2)) goto 2

			! If no rule triggered, the bonds are in the wrong order
			mark = 1
			tmp = protein%bond(i)
			protein%bond(i) = protein%bond(i+1)
			protein%bond(i+1) = tmp
2		end do
		if(mark.eq.0) exit
	end do

	! Reset bond id
	do i = 0,protein%num_bond -1
		protein%bond(i)%id = i +1
	end do

end subroutine xprotein_translate

subroutine Protein_Show_Contents(protein)
	type(tXProtein) :: protein

	write(*,*) 'Protein molecule: ', adjustl(trim(protein%name))

	call Protein_Show_Atoms(protein) 
	call Protein_Show_Rings(protein)
	call Protein_Show_Bonds(protein) 
end subroutine Protein_Show_Contents

subroutine Protein_Show_Atoms(protein)
	type(tXProtein) :: protein
	integer ::  i 

	write(*,*) 'Total number of atoms in this molecule = ',protein%num_atom
	write(*,'(t1,a)') 'id  vld.  type   xtype   type2   res  r_id  name   &
	    &weight  r   eps  q  XR logp solv hb  occ. bfct  ring origin &
	    &part   #neib #nonh  neib[0..6]                  bond[0..6]                 hb root[0..2]'

	write(*,'()')
	do i = 0,protein%num_atom -1 
		write(*,'(t1,i4, t6,i2, t11,a, t18,a, t27,a, t34,a, t39,a, &
		    &t45,a, t52,f5.2, t59,f4.1, t63,f4.1, t67,f4.1, t71,f4.1, &
		    &t75,f4.1, t79,f4.1, t85,a, t89,f4.1, t93,f4.1, t100,i1, &
		    &t105,i1, t112,i1, t117,i1, t123,i1, t130,i4,t134,i4,t138, &
		    &i4,t142,i4,t146,i4,t150,i4,t154,i4, t158,i4,t162,i4,t166, &
		    &i4,t170,i4,t174,i4,t178,i4,t182,i4, t186,f5.1,t191,f5.1, &
		    &t196,f5.1)') &

			protein%atom(i)%id+1, &
			protein%atom(i)%valid, &
			protein%atom(i)%ttype, &
			protein%atom(i)%xtype, &
			protein%atom(i)%type2, &
			trim(protein%atom(i)%residue), &
			trim(protein%atom(i)%res_id(5:10)), &
			trim(protein%atom(i)%name), &
			protein%atom(i)%weight, &
			protein%atom(i)%r, &
			protein%atom(i)%eps, &
			protein%atom(i)%q, &
			protein%atom(i)%XR, &
			protein%atom(i)%logp, &
			protein%atom(i)%solv, &
			protein%atom(i)%hb, &
			protein%atom(i)%occupancy,&
			protein%atom(i)%bfactor,&
			protein%atom(i)%ring,&
			protein%atom(i)%origin,&
			protein%atom(i)%part,&
			protein%atom(i)%num_neib,&
			protein%atom(i)%num_nonh,&
			protein%atom(i)%neib(0),protein%atom(i)%neib(1), &
			protein%atom(i)%neib(2),protein%atom(i)%neib(3), &
			protein%atom(i)%neib(4),protein%atom(i)%neib(5), &
			protein%atom(i)%neib(6), &
			protein%atom(i)%bond(0),protein%atom(i)%bond(1), &
			protein%atom(i)%bond(2),protein%atom(i)%bond(3), &
			protein%atom(i)%bond(4),protein%atom(i)%bond(5), &
			protein%atom(i)%bond(6), &
			protein%atom(i)%root(0),protein%atom(i)%root(1),protein%atom(i)%root(2)
	end do
end subroutine Protein_Show_Atoms

subroutine Protein_Show_Rings(mol)
	type(tXProtein) :: mol
	integer ::  i 

	write(*,'(a,i4)') 'Total number of rings in this molecule = ', size(mol%ring)
	write(*,'(t1,a)') 'id  vld. type size    atoms (first-last)  centroid'

	do i = 0,size(mol%ring) -1 
		if(mol%ring(i)%valid.ne.0) then
			write(*,'(t1,i3, t6,i1, t10,i1, t15,i2, t23,i5,t29,i5, t43,f5.1,t49,f5.1,t56,f5.1)') &
					i+1, mol%ring(i)%valid,mol%ring(i)%ttype,mol%ring(i)%num_member, &
					mol%ring(i)%atom_id(0),mol%ring(i)%atom_id(mol%ring(i)%num_member-1), &
					mol%ring(i)%centroid(0),mol%ring(i)%centroid(1),mol%ring(i)%centroid(2)
		else
			write(*,'(t1,i3, t6,a)') i+1, 'invalid'
		end if
	end do
end subroutine Protein_Show_Rings

subroutine Protein_Show_Bonds(protein)
	type(tXPRotein)		:: protein
	integer							:: i,j,filestat
	character(len=4)		:: buf
	character(len=16)		:: ring,slask,ar,pol,planar,charge,hb,bondtype,mask
	character(len=64)   :: ringpath,neibs,bonds


	write(*,'(a,i6)') 'Total number of bonds in this molecule = ',protein%num_bond
	if(input%show_bonds.eq.'YES') then
		write(*,'(a)') '  ID TYPE    ATOM-ATOM  RING  NEIGHBOURING BONDS'

		do i = 0,protein%num_bond -1
!			if(protein%atom(protein%bond(i)%atom_1-1)%mask.eq.0 .OR. protein%atom(protein%bond(i)%atom_1-1)%mask.eq.0) cycle
			ring = ''
			if(protein%bond(i)%ring.eq.1)			ring		 = 'ring'
			if(protein%bond(i)%ttype.eq.'1')  bondtype = 'single'
			if(protein%bond(i)%ttype.eq.'2')  bondtype = 'double'
			if(protein%bond(i)%ttype.eq.'3')  bondtype = 'triple'
			if(protein%bond(i)%ttype.eq.'ar') bondtype = 'ar'
			if(protein%bond(i)%ttype.eq.'am') bondtype = 'amide'

			neibs = ''
			do j = 0, protein%bond(i)%num_neib-1
				write(buf,'(i4)') protein%bond(i)%neib(j)
				neibs((j)*5+1:(j+1)*5) = buf
			end do

			write(*,'(t1,i4,t6,a,  t14,i4,t19,i4,  t25,a, t31,a)') &
						i+1, bondtype, protein%bond(i)%atom_1,protein%bond(i)%atom_2, ring, neibs
		end do
	end if
end subroutine

subroutine Protein_Define_Pocket(protein,ff, ligand, cutoff)
	type(tXProtein)	:: protein
	type(tXLigand)	:: ligand
	type(tXFF)			:: ff
	real(kind=prec)						:: cutoff
	integer					:: i,j,mark,count 
	real(kind=prec)						:: d 
	integer					:: num_pocket_res  
  type(tXResidue),pointer	:: pocket_res(:)
	type(tXRing)		:: tmp_ring 

	num_pocket_res = 0 

	allocate(pocket_res(0:300 -1),stat=err)			! worst case? not my idea...
	if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'
	! find pocket atoms on the protein and label their valid as '2',
	! including the waters close to the ligand
	do i = 0,protein%num_atom -1 
		if(protein%atom(i)%valid<=0) then
			goto 1
		end if

		! filter out all of the non-polar hydrogen atoms
		! since X-Score uses unit-atom model
		if((protein%atom(i)%xtype.eq.'H').and.(protein%atom(i)%residue.ne.'COF')) then
			protein%atom(i)%valid=0
			goto 1
		end if


		! check if this atom is close to the ligand
		mark=0 
		do j = 0,ligand%mol%num_atom -1 
			if(ligand%mol%atom(j)%valid<=0) then
				goto 2
			elseif(ligand%mol%atom(j)%ttype.eq.'H') then
				goto 2
			else
				d=Distance(protein%atom(i)%coor,ligand%mol%atom(j)%coor) 
				if(d>cutoff) then
					goto 2
				else
					mark=1
					exit
				end if
			end if
2		end do

		
		if(mark.eq.0) goto 1 ! not a pocket atom

		protein%atom(i)%valid=2 

		! if this atom is a water or a metal ion, do not check its residue
		if(protein%atom(i)%xtype.eq.'O.w') goto 1
		if(protein%atom(i)%xtype.eq.'M+') goto 1

		! for regular atoms, check if it has been included in a 
		! known pocket residue  if not, add that newly found pocket residue	
		mark=0 

		do j = 0,num_pocket_res -1
			if((protein%atom(i)%residue.eq.pocket_res(j)%name).and. &
				 (protein%atom(i)%res_id.eq.pocket_res(j)%id)		.and. &
		     (protein%atom(i)%chain.eq.pocket_res(j)%chain)) then
				mark=1
				exit
			else
				goto 3
			end if
3		end do

		if(mark.eq.0)  then ! new pocket residue found
			pocket_res(num_pocket_res)%name		= protein%atom(i)%residue 
			pocket_res(num_pocket_res)%id			= protein%atom(i)%res_id 
			pocket_res(num_pocket_res)%chain	= protein%atom(i)%chain 
			pocket_res(num_pocket_res)%valid	= 1  
			num_pocket_res = num_pocket_res +1 
		end if
1	end do

	! now check if the binding pocket is well defined
	count=0 
	do i = 0,num_pocket_res -1 
		if(pocket_res(i)%name.eq.'HET') goto 4 
		if(pocket_res(i)%name.eq.'COF') goto 4 
		if(pocket_res(i)%name.eq.'WAT') goto 4 
		if(pocket_res(i)%name.eq.'HOH') goto 4 
		count = count +1 
4	end do

	if(count<3) then ! there should be at least 3 residues as pocket!
		write(*,'(a,i1)') 'Error: not enough binding pocket residues on the protein. Residues found = ', count
		write(*,*) 'Probably the ligand has not been docked with the protein.' 
		stop
	end if

	! define all the left atoms in pocket residues as pocket atoms 
	do i = 0,protein%num_atom -1 
		if(protein%atom(i)%valid<=0)				goto 5 
		if(protein%atom(i)%valid.eq.2)			goto 5 
		if(protein%atom(i)%ttype.eq.'O.w') goto 5 
		if(protein%atom(i)%hb.eq.'M')			goto 5 

		do j = 0,num_pocket_res -1 
			if(protein%atom(i)%residue.ne.pocket_res(j)%name)	goto 6 
			if(protein%atom(i)%res_id.ne.pocket_res(j)%id)			goto 6 
			if(protein%atom(i)%chain.ne.pocket_res(j)%chain)		goto 6 
			protein%atom(i)%valid=2
			exit
6		end do
5	end do

	! now detect all of the aromatic rings within pocket residues
	! this information is need for later scoring purposes
	if(associated(protein%ring)) then
		do i = 0,size(protein%ring) -1
			call Ring_Clear(protein%ring(i))
		end do
		deallocate(protein%ring)
	end if
	protein%num_ring=0

	do i = 0,num_pocket_res -1
		if(pocket_res(i)%name.ne.'PHE' .and. &
			 pocket_res(i)%name.ne.'TYR' .and. &
			 pocket_res(i)%name.ne.'HIS' .and. &
       pocket_res(i)%name.ne.'TRP') goto 7

		call Ring_Clear(tmp_ring)

		do j = 0,protein%num_atom -1 
			if(protein%atom(j)%valid.ne.2) goto 8 
			if(protein%atom(j)%part.ne.1) goto 8 
			if(protein%atom(j)%ring.ne.2) goto 8 
			if(protein%atom(j)%chain.ne.pocket_res(i)%chain) goto 8 
			if(protein%atom(j)%residue.ne.pocket_res(i)%name) goto 8 
			if(protein%atom(j)%res_id.ne.pocket_res(i)%id) goto 8 

			!tmp_ring%atom_id%push_back(j+1)    
			call Int_PushBack(tmp_ring%atom_id,j+1)

			tmp_ring%centroid(0) = tmp_ring%centroid(0) + protein%atom(j)%coor(0) 
			tmp_ring%centroid(1) = tmp_ring%centroid(1) + protein%atom(j)%coor(1) 
			tmp_ring%centroid(2) = tmp_ring%centroid(2) + protein%atom(j)%coor(2) 
8		end do

		tmp_ring%num_member = size(tmp_ring%atom_id)

		if(tmp_ring%num_member>0) then
			tmp_ring%centroid(0)=tmp_ring%centroid(0)/tmp_ring%num_member 
	 		tmp_ring%centroid(1)=tmp_ring%centroid(1)/tmp_ring%num_member 
	 		tmp_ring%centroid(2)=tmp_ring%centroid(2)/tmp_ring%num_member 
	 		tmp_ring%valid=1
			tmp_ring%ttype=2 
	 		
			!this->ring%push_back(tmp_ring) 
			Call Protein_Add_Ring(protein, tmp_ring)
			protein%num_ring = protein%num_ring +1
		end if
7	end do

	call Ring_Clear(tmp_ring)

	deallocate(pocket_res)
end subroutine Protein_Define_Pocket

subroutine Protein_Add_Ring(m,r)
	type(tXProtein)				:: m
	type(tXRing)					:: r
	integer								:: i,j
	type(tXRing),pointer	:: tmp(:)

	if(.not. associated(m%ring)) then
		! Adding first ring to molecule, need to init pointer
		allocate(m%ring(0:0),stat=err)
		if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'

		Call Ring_Copy(m%ring(0),r)
	else
		! Adding a second (or more) ring to molecule
	
		i = size(m%ring)
				
		! Increase size of ring-array
		allocate(tmp(0:i-1),stat=err)
		if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'
		do j  = 0,i -1
			call Ring_Copy(tmp(j),m%ring(j))
		end do

		deallocate(m%ring)
		allocate(m%ring(0:i),stat=err)
		if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'
		do j = 0, i-1
			call Ring_Copy(m%ring(j),tmp(j))
		end do

		do j = 0,i-1
			call Ring_Clear(tmp(j))
		end do
		deallocate(tmp)

		! Add ring
		call Ring_Copy(m%ring(i), r)
	end if
end subroutine Protein_Add_Ring


subroutine Protein_Value_Atom(protein,flag,lite_run)
	type(tXProtein),intent(inout)			:: protein
	integer,optional									:: lite_run
	integer														:: flag
	integer														:: slask,i

	! assign the parameters for all atoms, only done once
	if(.not. present(lite_run)) then
		do i=0, protein%num_atom -1
			if(protein%atom(i)%iscofactor.eq.0) slask = ForceField_Assign_Patom_Parameters(protein%atom(i))
		end do
	end if

	do i = 0,protein%num_atom -1
		if(protein%atom(i)%iscofactor.ne.0) protein%atom(i)%valid = 1
	end do

	! build connection tables, must assign parameters first
	call Protein_Detect_Connections(protein,lite_run)

	! calculate H-bond root
	if(flag.ne.0)	call Protein_Calculate_HB_Root(protein)
end subroutine Protein_Value_Atom

subroutine Protein_Calculate_HB_Root(protein)
	type(tXProtein)	:: protein
	integer ::  i,j,id,count 
	real(kind=prec) ::  tmpx,tmpy,tmpz 

	! (1) we do calculate root for polar hydrogens 
	! (2) we do NOT calculate root for non-HB atoms
	! (3) we do NOT calculate root for metal ions
	! (4) we do NOT calculate root for waters   
	! (5) we do calculate root for the oxygen atoms on SO4 and PO4

	do i = 0,protein%num_atom -1 
		if(protein%atom(i)%valid<=0) then
			goto 1
		elseif(protein%atom(i)%xtype.eq.'H') then 
			goto 1 
		elseif(protein%atom(i)%xtype.eq.'O.w') then 
			goto 1 
		elseif(protein%atom(i)%hb.eq.'N') then 
			goto 1 
		elseif(protein%atom(i)%hb.eq.'H') then 
			goto 1 
		elseif(protein%atom(i)%hb.eq.'P') then 
			goto 1 
		elseif(protein%atom(i)%hb.eq.'M') then 
			goto 1 
		end if

		if(protein%atom(i)%num_neib.eq.0) goto 1

		tmpx=0
		tmpy=0
		tmpz=0
		count=0

		do j = 0,protein%atom(i)%num_neib -1 
		 id=protein%atom(i)%neib(j)-1 
!		 write(*,*) i,j, id
		 if(protein%atom(id)%ttype.eq.'H') then
			goto 2
		 else
		 	 tmpx=tmpx + (protein%atom(id)%coor(0)) 
		 	 tmpy=tmpy + (protein%atom(id)%coor(1)) 
		 	 tmpz=tmpz + (protein%atom(id)%coor(2)) 
			 count = count +1 
		 end if
2		end do

		if(count.eq.0) goto 1

		protein%atom(i)%root(0)=tmpx/count 
		protein%atom(i)%root(1)=tmpy/count 
		protein%atom(i)%root(2)=tmpz/count 
1	end do
end subroutine Protein_Calculate_HB_Root

subroutine Protein_Detect_Connections(protein,lite_run)
	type(tXProtein)		:: protein
	integer,optional	:: lite_run
	integer ::  i,j,k,count,id 
	integer ::  mark 
	real(kind=prec) ::  cutoff,d 

	! build the connections for regular atoms first
	cutoff=2.00   ! covalent bond distance cutoff

	do i = 0,protein%num_atom -1 
	 if((protein%atom(i)%valid.eq.0)) then		! allow invalid cofactor atoms
		goto 1
	 elseif(protein%atom(i)%part.ne.1) then
		goto 1
	 elseif(protein%atom(i)%num_neib>0) then
		goto 1  ! already done 
	 end if
	 count=0 

	 if(present(lite_run)) goto 100

	 ! find internal connections inside the same residue first
	 do j = 0,protein%num_atom -1 
		 if(i.eq.j) then	
			goto 2 
		 elseif(protein%atom(j)%valid.eq.0) then 
			goto 2 
		 elseif(protein%atom(j)%part.ne.1) then 
			goto 2 
		 elseif(protein%atom(i)%chain.ne.protein%atom(j)%chain) then 
			goto 2 
		 elseif(protein%atom(i)%residue.ne.protein%atom(j)%residue) then 
			goto 2 
		 elseif(protein%atom(i)%res_id.ne.protein%atom(j)%res_id) then 
			goto 2 
		 end if
		 
		 ! there is no bond between hydrogen atoms
		 if((protein%atom(i)%ttype.eq.'H').and.(protein%atom(j)%ttype.eq.'H')) then
			goto 2
		end if

		! now check whether atom i and atom j are covalently bound
		if(ForceField_Patom_Connection_Test(protein%atom(i),protein%atom(j)).eq.1) then
			 if(count>=MAX_ATOM_NEIB) then
				goto 2   ! already full
			 else
			 	 protein%atom(i)%neib(count)=j+1
				 count = count +1
			 end if
		 else 
			goto 2
		end if
2	 end do

! JUMP HERE IF LITE RUN
100 continue
	 ! now detect peptide amide bonds
	 if((protein%atom(i)%name.eq.'N').and.(protein%atom(i)%ttype.eq.'N.am')) then
		do j = 0,protein%num_atom -1 
			 if(protein%atom(j)%valid.eq.0) then 
				goto 3 
			 elseif(protein%atom(j)%part.ne.1) then 
				goto 3 
			 elseif(protein%atom(i)%chain.ne.protein%atom(j)%chain) then 
				goto 3 
			 elseif(protein%atom(i)%res_id.eq.protein%atom(j)%res_id) then 
				goto 3 
			 elseif(protein%atom(j)%ttype.ne.'C.2') then 
				goto 3 
			 elseif(protein%atom(j)%name.ne.'C') then 
				goto 3 
			 end if

			 d=Distance(protein%atom(i)%coor,protein%atom(j)%coor) 

			 if(d>cutoff) then 
				goto 3 
			 elseif(count>=MAX_ATOM_NEIB) then 
				goto 3  
			 else
			 	 protein%atom(i)%neib(count)=j+1  
				 count = count +1 
			 end if
3		end do
	elseif((protein%atom(i)%name.eq.'C').and.(protein%atom(i)%ttype.eq.'C.2')) then
		do j = 0,protein%num_atom -1 
			if(protein%atom(j)%valid.eq.0) then 
				goto 4 
			elseif(protein%atom(j)%part.ne.1) then 
				goto 4 
			elseif(protein%atom(i)%chain.ne.protein%atom(j)%chain) then 
				goto 4 
			elseif(protein%atom(i)%res_id.eq.protein%atom(j)%res_id) then 
				goto 4 
			elseif(protein%atom(j)%ttype.ne.'N.am') then 
				goto 4 
			elseif(protein%atom(j)%name.ne.'N') then 
				goto 4 
			end if

			d=Distance(protein%atom(i)%coor,protein%atom(j)%coor) 

			if(d>cutoff) then 
				goto 4 
			elseif(count>=MAX_ATOM_NEIB) then 
				goto 4  
			else
				 protein%atom(i)%neib(count)=j+1  
				 count = count +1 
			end if
4		end do
	 end if

	 protein%atom(i)%num_neib=count 
1 end do

	! now find the patoms bound to the metal ions
	cutoff=3.00   ! M-bond distance cutoff

	do i = 0,protein%num_atom -1 
	 if(protein%atom(i)%valid.eq.0) then 
		goto 5 
	 elseif(protein%atom(i)%part.eq.1) then 
		goto 5 
	 elseif(protein%atom(i)%xtype.ne.'M+') then 
		goto 5 
	 elseif(protein%atom(i)%num_neib>0) then 
		goto 5   ! already done
	 end if

	 count=0 

	 do j = 0,protein%num_atom -1 
		 if(protein%atom(j)%valid.eq.0) then 
			goto 6 
		 elseif(protein%atom(j)%part.ne.1) then 
			goto 6 
		 elseif((protein%atom(j)%hb.ne.'A').and.(protein%atom(j)%hb.ne.'DA')) then 
			goto 6 
		 end if

		 d=Distance(protein%atom(i)%coor,protein%atom(j)%coor) 

		 if(d>cutoff) then 
			goto 6 
		 elseif(count>=MAX_ATOM_NEIB) then 
			goto 6 
		 else 
			 protein%atom(i)%neib(count)=j+1   ! note this
			 count = count +1 
			 ! here we do not add connections to atom(j)
		 end if
6		end do
		protein%atom(i)%num_neib=count 
5	end do

	! Return if lite run
	if(present(lite_run)) return

	! now build the connections for SO4s and PO4s
	do i = 0,protein%num_atom -1 
	 if(protein%atom(i)%valid.eq.0) then 
		goto 7 
	 elseif(protein%atom(i)%part.eq.1) then 
		goto 7 
	 elseif((protein%atom(i)%residue.ne.'SO4').and. &		! check this
					(protein%atom(i)%residue.ne.'PO4')) then 
		goto 7 
	 elseif((protein%atom(i)%name.ne.'S').and.	&				! check this
					(protein%atom(i)%name.ne.'P')) then 
		goto 7 
	 end if

	 if(protein%atom(i)%num_neib>0) goto 7   ! done 

	 ! now check the satellite atoms for this P or S atom
	 do j = 0,protein%num_atom -1 
		 if(protein%atom(j)%valid.eq.0) then 
			goto 8 
		 elseif(protein%atom(j)%part.eq.1) then 
			goto 8 
		 elseif((protein%atom(j)%residue.ne.'SO4').and.(protein%atom(j)%residue.ne.'PO4')) then 
			goto 8 
		 elseif((protein%atom(j)%res_id.ne.protein%atom(i)%res_id)) then 
			goto 8 
		 elseif(index(protein%atom(j)%name,'S').ne.0) then 
			goto 8 
		 elseif(index(protein%atom(j)%name,'P').ne.0) then 
			goto 8 
		 end if

		 count=protein%atom(i)%num_neib
		 mark=0
		
		 do k = 0,count -1 
			 if(protein%atom(i)%neib(k).eq.(j+1)) then
				 mark=1
				 exit
			 else
				goto 9
			 end if
9		 end do
	
		 if(mark.eq.0) then
			 protein%atom(i)%neib(count)=j+1   ! note this
			 protein%atom(i)%num_neib = protein%atom(i)%num_neib +1 
		 end if

		 count=protein%atom(j)%num_neib
		 mark=0

		 do k = 0,count -1 
			 if(protein%atom(j)%neib(k).eq.(i+1)) then
				 mark=1
				 exit
			 else
				goto 10
			 end if
10	 end do
		 
	 if(mark.eq.0) then
			 protein%atom(j)%neib(count)=i+1   ! note this
			 protein%atom(j)%num_neib = protein%atom(j)%num_neib +1 
	  end if
8	 end do
7 end do

	! now determine number of heavy atoms for each valid atom
	do i = 0,protein%num_atom -1 
		 if(protein%atom(i)%valid<=0) goto 12

		 count=0 
		 do j = 0,protein%atom(i)%num_neib -1 
			 id=protein%atom(i)%neib(j) 
			 if(protein%atom(id-1)%ttype.eq.'H') then
				goto 13
			 else 
				count = count +1 
			 end if
13	 end do

		 protein%atom(i)%num_nonh=count 
12 end do
end subroutine Protein_Detect_Connections

end module CALC_XSCORE
