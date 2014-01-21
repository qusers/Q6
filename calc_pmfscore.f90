!  (C) 2000 Uppsala Molekylmekaniska HB, Uppsala, Sweden
!
!  calc_pmfscore.f90
!
! Implementation by Peter Hanspers

module CALC_PMF
	use CALC_BASE
	use MASKMANIP
	use TRJ
	use	TOPO
	use PRMFILE
	use INDEXER
	use QATOM
	use MISC

	implicit none

	type tPMFEdge
		integer							:: valid
		integer,pointer			:: p(:)																					! connectivity path, 0..n-1
	end type

	type tPMFVertex
		integer							:: valid
		integer							:: id																						! ids are 0..n-1
		integer							:: num_conn																			! number of connections
	end type

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

	type tABSDetail
		real										:: score,r
		integer									:: atom_id,m
		character(len=3)				:: k,l
	end type

	type tPMFAtom
		integer									:: topindex			! topology index of this atom
		character(len=3)				:: ttype				! PMF type, two letter code
		character(len=8)				:: qtype				! Q type atom type short version
		character(len=8)				:: qtype_long	  ! Q type atom type long version
		character(len=8)				:: mol2type			! MOL2 type atom type
		real(8),dimension(0:2)	:: coor					! coordinates

		integer									:: polar				! 1 if atom is polar, 0 otherwise
		character(len=1)				:: hb						! D if atom is h-bond donor, A if acceptor
		integer									:: aromatic			! 1 if atom is aromatic, 0 if aliphatic
		integer									:: charge				! -1/0/1 if atom is negatively/un/positively charged
		character(len=3)				:: hybr					! hybridizing state
		integer									:: planar				! 1 if in planar ring, otherwise 0

		integer									:: ring					! ring indicator: 1=normal 2=aromatic
		integer,dimension(3)		:: ring_id			! id's of what rings this atom is part of
		integer									:: ring_count		! number of rings this atom is part of

		integer									:: bSecondPass	! flag inidicating if atom needs a second pass in assigning types

		real										:: abs					! atomic binding score
		type(tABSDetail),pointer:: absdetail(:) ! array of binding scores for each counter-atom considered
		integer									:: abscount			! count for absdetail(:)
		integer									:: nonzero_score	! flag indicating if atom ever got a non zero score
		
		! XScore members
		integer									:: id  	        !* atom id 
		integer									:: valid				!* valid indicator
		integer									:: mask					!  mask flag
		character(len=10)				:: residue      !* residue name
		character(len=10)				:: res_id       !  residue id: must be a string for PDB
		integer									:: atom_res_id	!	 atom number within residue
		integer									:: iscofactor		!  flag indicating if atom is part of a cofactor (1) or not (0)
		character								:: chain        !  chain label
		real										:: q						!* partial atomic charge
		real										:: weight
		integer									:: origin				!  atom origin indicator: 1=ligand 2=protein 
		integer									:: part					!  component indicator: for protein atoms: 1=ATOM 2=HETATM
		integer									:: num_neib			!  number of neighboring atoms
		integer									:: num_nonh			!  number of non-H neighboring atoms
		integer,dimension(0:6)	:: neib					!  ID of neighboring atoms
		integer,dimension(0:6)	:: bond					!  ID of bonds
		integer									:: num_bond			!	 number of bonds to/from this atom
	end type

	type tPMFBond
		character(len=3)				:: ttype				!  bond type
		integer									:: id
		integer									:: valid				!  valid indicator
		integer									:: atom_1       !  ID of the atom_1
		integer									:: atom_2       !  ID of the atom_2
		integer									:: part					!  ID of the component
		integer									:: ring					!  ring indicator
		integer									:: num_neib     !  number of neighboring bonds
		integer,dimension(0:5)	:: neib					!	 ID of neighboring bonds
	end type

	type tPMFRing
		integer,pointer					:: path(:)
		integer									:: aromatic
		integer									:: planar
	end type

	type tPMFMolecule
		character(len=64)				:: name
		real										:: score

		integer									:: num_subst

		integer									:: num_atom
		type(tPMFAtom),pointer	:: atom(:)

		integer									:: num_bond
		type(tPMFBond),pointer	:: bond(:)


		type(tPMFRing),pointer	:: ring(:)		
		integer									:: num_ring
	end type tPMFMolecule

	type tPMFInput
		character(len=256)			:: ff						! force field used
		character(len=20)				:: show_abs
		character(len=20)				:: show_detailed_abs
		character(len=10)				:: show_total
		character(len=3)				:: show_ligand
		character(len=3)				:: show_protein
		character(len=3)				:: show_bonds
		character(len=3)				:: show_hydrogens
		character(len=3)				:: show_masked

		character(len=3)				:: ignore_waters
		integer									:: maximum_ring_size
		character(len=3)				:: aromatic_require_ring

		integer									:: qatom

		character(len=256)			:: pmfs
		character(len=256)			:: pmf_names
		character(len=256)			:: rules
		character(len=256)			:: atom_translations
	end type

	type tPMFRule
		character(len=2)									:: ttype		! atom type that this rule applies to
		character(len=48),dimension(1:24) :: rule
		integer														:: count
	end type

	type donor
		integer(AI)							::	heavy,hydr				! topology numbers for atoms in H-bond donor
	end type donor

	type tPMFScore
		integer	:: frame
		real		:: score
	end type tPMFScore

	type(tPMFRule), pointer							:: protRule(:)
	integer															:: prot_num_rule
	type(tPMFRule), pointer							:: ligRule(:)
	integer															:: lig_num_rule

	type ATOM_TYPE_CONVERSION
		character(len=20)	:: mol2, other
	end type

	integer, public						:: SINGLE		= 1
	integer, public						:: DOUBLE		= 2
	integer, public						:: AROMATIC = 3
	integer, public						:: AMIDE		= 4

	character*80		:: top_file, fep_file
	character*80		:: atom_data_file
	character*80		:: coord_file	

	type(q_atom),private,allocatable,target	::	q_atoms(:)		! atoms in ligand
	type(q_bond),private,allocatable,target	::	q_bonds(:)		! bonds in ligand
	integer, allocatable					::	iqatom(:)		!one element per atom, 0 if not Q-atom, 
	
	type(q_atom),private,allocatable,target	::	aHQ(:)			! heavy atoms in ligand
	type(q_bond),private,allocatable,target	::	aHB(:)			! heavy bonds in ligand
																!else number in iqseq (ligand)
	integer	:: nHQ		! number of heavy q-atoms
	integer	:: nHB		! number of bonds between heavy q-atoms

	integer(AI), allocatable	::	lph_r(:)	!topology number for all lipophilic atoms in the receptor
	integer(AI), allocatable	::	lph_l(:)	!				''						     the ligand
	type(donor), allocatable	::	hbd_r(:)	!H-bond donor in receptor
	type(donor), allocatable	::	hbd_l(:)	!  ''            ligand
	integer(AI), allocatable	::	hba_r(:)	!H-bond acceptor atoms in receptor 
	integer(AI), allocatable	::	hba_l(:)	!      ''                 ligand
	integer(AI), allocatable	::	met_r(:)	!metal atoms in the receptor
!	type(wat), allocatable		::	waters(:)	!water molecules

	!numbers of each of the atom types
	integer						::	nlph_r, nlph_l, nhbd_r, nhbd_l
	integer						::	nhba_r, nhba_l, nmet_r, nwaters,nqbonds
	integer						::	nhbd_prot, nhba_prot
	integer						::	nrings		!number of rings in the ligand
	
	integer, allocatable							::	pol_con(:)	!connections for potentially polar atoms

	integer, parameter								::	MAX_MASKS = 10
	type(MASK_TYPE), private, target	::	masks(MAX_MASKS)
	integer, private									::	Nmasks = 0


	integer														:: bDoTopCalc
	integer														:: bSecondPass

	! System
	type(tPMFMolecule), public							:: protein
	type(tPMFMolecule), public							:: ligand
	type(tPMFMolecule), public,allocatable	:: cofactor(:)
	integer,pointer													:: top2prot(:)				! topology atom index -> protein atom index translation matrix, used in bond translations

	! Input
	type(tPMFInput),public									:: input
	character(len=256)											:: chInput
	type tPMFNames
		character(len=8)											:: pmfname, qname
	end type
	type(tPMFNames),pointer									:: pmf_names(:)
	type(ATOM_TYPE_CONVERSION), allocatable	:: atomtypetable(:)
	integer																	:: num_atomtype

	! Data
	real, dimension(34,16,60)								:: pmfdata																! (ligand type, protein type, bin)
	character(len=2),dimension(2,34)				:: pmftype

	integer, parameter :: nProt = 16, nLig = 34

  data pmftype(1,:) / 'CF','CP','cF','cP','C3','CW','CO','CN','NC','NP', &					! ligand atom types
											'NA','ND','NR','N0','NS','OC','OA','OE','OR','OS', &
											'OD','P ','SA','SD','HL','Zn','CL','Mn','Mg','F ', &
											'Fe','Br','V ','C0'/

	data pmftype(2,:) / 'CF','CP','cF','cP','CO','CN','NC','ND','NR','OC', &					! protein atom types
											'OA','OD','OW','SA','SD','HH','--','--','--','--', &
											'--','--','--','--','--','--','--','--','--','--', &
											'--','--','--','--'/

	! Results
	type(tPMFScore), pointer		:: pmfscores(:)		! pmfscoring stats for each frame
	integer											:: nPMFScores, maxPMFScores


	integer :: MAX_BOND_NEIB = 10
	integer :: MAX_ATOM_NEIB = 6
	integer,public				:: err									! public allocation error indicator
	integer,public				:: warn,waterwarn

contains

subroutine PMFReadAtomTypeConversions(filename, chType)
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
							atomtypetable(count)%mol2  = tmp(1)
							atomtypetable(count)%other = tmp(column)
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
			if(trim(atomtypetable(i)%other).eq.trim(atomtypetable(j)%other).and. &
			trim(atomtypetable(i)%mol2).ne.trim(atomtypetable(j)%mol2)) then
				write(*,'(a,a,a,a,a,a,a,a)') &
				'Fatal: Ambiguous entries in type translation matrix (file ', trim(adjustl(filename)), '). ', &
					 trim(adjustl(atomtypetable(j)%other)), ' maps to both ', &
					 adjustl(trim(atomtypetable(i)%mol2)), ' and ', &
					 adjustl(trim(atomtypetable(j)%mol2))
				stop
			end if
		end do
	end do

	write(*,'(a,i4,a,a,a)') 'Read ', count, ' atom translations for force field ', adjustl(trim(chType)), '.'
	return

999 write(*,'(a,a,a)') '>>>>> ERROR: Unable to read ', adjustl(trim(filename)), '. Exiting.'
		stop
end subroutine PMFReadAtomTypeConversions


subroutine pmfRead_Input
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
						cycle																!	skip blank line
					elseif(trim(buf(1:1)).eq.'#') then
						cycle																! skip comment
					else
						p = index(buf,'#')
						if(p.ne.0) buf = buf(1:p-1)					! remove traling comment
						
						read(buf, *,iostat=filestat) head
						call upcase(head)
						head = adjustl(trim(head))

						read(buf,*,iostat=filestat) head, chVal
						
						select case(head)			
							case('SHOW_ATOM_BIND_SCORE')
								read(buf,*,iostat=filestat) head, input%show_abs
							case('SHOW_DETAILED_ABS')
								read(buf,*,iostat=filestat) head, input%show_detailed_abs
							case('SHOW_TOTAL_SCORE')
								read(buf,*,iostat=filestat) head, input%show_total
							case('SHOW_LIGAND')
								read(buf,*,iostat=filestat) head, input%show_ligand
							case('SHOW_PROTEIN')
								read(buf,*,iostat=filestat) head, input%show_protein
							case('SHOW_BONDS')
								read(buf,*,iostat=filestat) head, input%show_bonds
							case('SHOW_HYDROGENS')
								read(buf,*,iostat=filestat) head, input%show_hydrogens
							case('SHOW_MASKED')
								read(buf,*,iostat=filestat) head, input%show_masked

							case('IGNORE_WATERS')
								read(buf,*,iostat=filestat) head, input%ignore_waters
							case('MAXIMUM_RING_SIZE')
								read(buf,*,iostat=filestat) head, input%maximum_ring_size
 
							case('POTENTIALS_OF_MEAN_FORCE')
								read(buf,'(a256)') line
								read(line,*) head
								!workaround to allow / in strings
								input%pmfs = adjustl(line(len_trim(head)+1:len(line)))
!								read(buf,*,iostat=filestat) head, input%pmfs
							case('NAME_TRANSLATIONS')
								read(buf,'(a256)') line
								read(line,*) head
								!workaround to allow / in strings
								input%pmf_names = adjustl(line(len_trim(head)+1:len(line)))
!								read(buf,*,iostat=filestat) head, input%pmf_names
							case('ATOM_TYPE_RULES')
								read(buf,'(a256)') line
								read(line,*) head
								!workaround to allow / in strings
								input%rules = adjustl(line(len_trim(head)+1:len(line)))
!								read(buf,*,iostat=filestat) head, input%rules
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
			end do
		close(fp)
		
		input%qatom = -1
	else
		input%show_abs							= 'NO'
		input%show_detailed_abs			= 'NO'
		input%show_total						= 'YES'
		input%show_ligand						= 'YES'
		input%show_protein					= 'NO'
		input%show_hydrogens				= 'YES'
		input%show_masked						= 'YES'
		input%aromatic_require_ring = 'YES'
		input%maximum_ring_size			=  7
		input%pmfs									=	'pmfs.dat'
		input%pmf_names							= 'pmf_names.dat'
		input%rules									= 'atdt_pmf.dat'
		input%atom_translations			= 'atom_type_translations.dat'
		input%qatom = -1
	end if

	return
999 write(*,'(a,a,a)') 'ERROR: Unable to read ', adjustl(trim(chInput)), '. Exiting.'
		stop
end subroutine pmfRead_Input

subroutine pmfRead_InputQatom
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

			do i = 0,ligand%num_atom -1
				if(adjustl(trim(ligand%atom(i)%res_id)).eq.adjustl(trim(head))) then
					input%qatom = i + res
					input%show_abs = 'QATOM'
					exit
				end if
			end do 
		else
			read(chVal,*,iostat=filestat) input%qatom
		end if
	
		input%show_abs = 'QATOM'
		if(input%qatom<0.or.input%qatom>ligand%num_atom -1) then
			write(*,'(a,a,a)') 'WARNING: Q-atom ', trim(chVal), ' does not exist. ABS will not be displayed.'
			warn = warn +1
			input%show_abs = 'NO'
		end if

	end if
end subroutine pmfRead_InputQatom

subroutine pmf_extract_rules(filename)
	character(len=256)			:: filename, buf
	integer									:: i,fp,num,section,filestat,p,count
	integer,dimension(1:2)	:: cnt

	character(len=48)									:: word
	character(len=48),dimension(1:24)	:: w
	
	type(tPMFRule)					:: rule

	fp = freefile()
	open(unit=fp,file=filename,err=999,action='READ')

    cnt(:) = 0																							! reset counters
    do
			read(fp,'(a256)',iostat=filestat) buf
			if(filestat.eq.0) then
				if(index(buf,'!')>0) buf = buf(1:index(buf,'!')-1)		! remove trailing comment
				if(buf(1:1).eq.'!') cycle															! skip REM
				if(trim(buf).eq.'') cycle															! skip blank line

				if(buf.eq.'Ligand atom types') then
					section = 1
					cycle
				elseif(buf.eq.'Protein atom types') then
					section = 2
					cycle
				elseif(section.eq.1.or.section.eq.2) then
					cnt(section) = cnt(section) +1
				end if
			else
				exit
			end if
		end do      

		rewind(fp)

		write(*,'(a,i3,a,i3,a,a)') 'Reading ', cnt(1), &
		' ligand atom type rules and ', cnt(2), &
		' protein atom type rules from ', trim(adjustl(filename))

		allocate( ligRule(1:cnt(1)))
		allocate(protRule(1:cnt(2)))
		

    do
			read(fp,'(a256)',iostat=filestat) buf
			if(filestat.eq.0) then
				if(index(buf,'!')>0) buf = buf(1:index(buf,'!')-1)		! remove trailing comment
				if(buf(1:1).eq.'!') cycle															! skip REM
				if(trim(buf).eq.'') cycle															! skip blank line

				if(buf.eq.'Ligand atom types') then
					section = 1
					cycle
				elseif(buf.eq.'Protein atom types') then
					section = 2
					cycle
				elseif(section.eq.1.or.section.eq.2) then
					! Analyze this
					w(:) = ''
					count = 1
					p = index(buf,':')
					word = adjustl(buf(1:p-1))
					w(count) = word
					buf = adjustl(buf(p+1:256))

					do
						p = index(buf,',')
						if(p>1) then	
							count = count +1
							w(count) = adjustl(buf(1:p-1))
							buf	 = adjustl(buf(p+1:256))
						elseif(p.eq.0.and.trim(adjustl(buf)).ne.'') then
							count = count +1
							w(count) = trim(adjustl(buf))
							exit
						end if
					end do
				end if
			else
				exit
			end if
			
			if(section.eq.1) then
				lig_num_rule = lig_num_rule +1
				ligRule(lig_num_rule)%count = count -1
				ligRule(lig_num_rule)%ttype = adjustl(trim(w(1)))
				ligRule(lig_num_rule)%rule(1:23) = w(2:24)
			elseif(section.eq.2) then
				prot_num_rule = prot_num_rule +1
				protRule(prot_num_rule)%count = count -1
				protRule(prot_num_rule)%ttype = adjustl(trim(w(1)))
				protRule(prot_num_rule)%rule(1:23) = w(2:24)
			end if
		end do      

	close(fp)

	return

999 write(*,*) 'FATAL: read error'
		stop

end subroutine

subroutine pmf_show_rule(rule)
	type(tPMFRule)	:: rule
	integer					:: i

	write(*,'(a,a)') rule%ttype , ' ->'
	do i = 1,rule%count
		write(*,'(a,a)') '      ', trim(rule%rule(i))
	end do
end subroutine pmf_show_rule

subroutine pmfRead_Data(filename)
	character(len=256)	:: filename,buf
	integer							:: i,j,ligtype,prottype,m,fp,p,filestat
	real								:: slask

	i = 1
	j = 1
	fp = freefile()
	write(*,'(a,a)') 'Reading pmf''s from ', adjustl(trim(filename))
	open(unit=fp,file=filename,err=999,action='READ')
		do i = 1,nProt
			do j = 1,nLig
				read(fp,'(a256)',iostat=filestat) buf
				buf = adjustl(trim(buf))
				p = index(buf,'--')
				if(p.ne.3) stop 'FATAL: error in parameter file'

				ligtype = -1
				do m = 1,nLig
					if(pmftype(1,m).eq.buf(5:6)) then
						ligtype = m																						
					end if
				end do
				if(ligtype.eq.-1) then
					write(*,'(a,a,a)') 'Fatal: Ligand type ', buf(5:6), ' is not a recognized atom type.'
					stop
				end if
				
				prottype = -1
				do m = 1,nProt
					if(pmftype(2,m).eq.buf(1:2)) then
						prottype = m																					
					end if
				end do
				if(prottype.eq.0) then
					write(*,'(a,a,a)') 'Fatal: Protein type ', buf(1:2), ' is not a recognized atom type.'
					stop
				end if
	
				do m = 1,60
					read(fp,*,iostat=filestat) slask, pmfdata(ligtype,prottype,m)
				end do
			
			end do
		end do
	close(fp)

	return
999 write(*,*) 'FATAL: read error'
		stop
end subroutine pmfRead_Data	

subroutine pmfRead_Names(filename)
	character(len=256)	:: filename,buf
	integer							:: i,j,fp,p,filestat,count
	character(len=8)		:: pmfname,qname

	i = 1
	j = 1
	fp = freefile()
	write(*,'(a,a)') 'Reading name translations from ', adjustl(trim(filename))
	open(unit=fp,file=filename,err=999,action='READ')
		count = 0
		do
			read(fp,'(a256)',iostat=filestat) buf

			if(filestat.eq.0) then
				if(trim(buf).eq.'') then
					cycle																!	skip blank line
				elseif(trim(buf(1:1)).eq.'#') then
					cycle																! skip comment
				else
					p = index(buf,'#')
					if(p.ne.0) buf = buf(1:p-1)					! remove traling comment
					
					read(buf,*,iostat=filestat) pmfname, qname
					count = count +1
				end if
			else
				exit
			end if
		end do

		rewind(fp)
		allocate(pmf_names(count))
	
		count = 0
		do
			read(fp,'(a256)',iostat=filestat) buf

			if(filestat.eq.0) then
				if(trim(buf).eq.'') then
					cycle																!	skip blank line
				elseif(trim(buf(1:1)).eq.'#') then
					cycle																! skip comment
				else
					p = index(buf,'#')
					if(p.ne.0) buf = buf(1:p-1)					! remove traling comment
					
					count = count +1
					read(buf,*,iostat=filestat) pmf_names(count)%pmfname, pmf_names(count)%qname
				end if
			else
				exit
			end if
		end do
	close(fp)

	return
999 write(*,*) 'FATAL: read error'
		stop
end subroutine

subroutine pmf_precalc
	character(len=256)	:: buf
	integer							:: i,j,iRes,iNextRes
	real								:: score


	! ============================
	! 0. READ RULES AND PARAMTERS
	! ============================

			call pmfRead_Input

			write(*,'(a)')   'Input parameters and switches:'
			write(*,'(t3,a,t31,a)')			'FORCE FIELD',						  trim(adjustl(input%ff))
			if(input%qatom>-1) then
				write(*,'(t3,a,t31,i6)')	'QATOM FOCUS',						  input%qatom
			else
				write(*,'(t3,a,t31,a)')		'QATOM FOCUS',						  '(no qatom focused)'
			end if
			write(*,'(t3,a,t31,i1)')		'MAXIMUM RING_SIZE',				input%maximum_ring_size
			write(*,'(t3,a,t31,a)')			'POTENTIALS OF MEAN FORCE', trim(adjustl(input%pmfs))
			write(*,'(t3,a,t31,a)')			'ATOM TYPE RULES',				    trim(adjustl(input%rules))
			write(*,'(t3,a,t31,a)')			'NAME TRANSLATIONS',				trim(adjustl(input%pmf_names))
			write(*,'(t3,a,t31,a)')			'TYPE TRANSLATIONS',				trim(adjustl(input%atom_translations))

			call pmf_extract_rules(input%rules)
			call pmfRead_Data(input%pmfs)
			call pmfRead_Names(input%pmf_names)
			call PMFReadAtomTypeConversions(input%atom_translations,input%ff)		! read atom type conversions

	! ============================
	! 1. READ AND TRANSLATE LIGAND
	! ============================
 
			if(.not. qatom_load_atoms(fep_file)) stop 'Failed to read Q-atom list from FEP file.'

			allocate(iqatom(nat_solute))	! Allocate mem for every solute atom
			iqatom(:) = 0
			do i=1,nqat
				iqatom(iqseq(i)) = i				! Set references so that iqatom = 0 if protein and {number} of ligand
			end do

			! Extract bonds between q-atoms
			call pmfsort_bonds						

			! Set topology numbers for q-atoms
			q_atoms(:)%top_nr = iqseq(:)

			! Translate q-atom data to ligand structure
			call pmfLigand_Translate(ligand,nqat,nqbonds,q_atoms,q_bonds,offset)
			ligand%name = adjustl(trim(title))
			write(*,'(a,i4,a,i4,a)') 'Translated ', ligand%num_atom, ' atoms and ', ligand%num_bond, ' bonds to ligand'

			! Assign atom types using ligand rules
			call pmfValue_Molecule(ligand,ligRule,lig_num_rule)

			! Find monitored q-atom if specified. This has to be done after ligand atom names are assigned.
			call pmfRead_InputQatom

			if(input%show_ligand.eq.'YES') call pmfShow_Molecule(ligand)
			ligand%atom(:)%nonzero_score = 0

	! =============================
	! 2. READ AND TRANSLATE PROTEIN
	! =============================

			allocate(top2prot(1:nat_solute))						! top -> protein index translation matrix
!			call pmfProtein_Translate(protein, offset, xtop, bnd,nbonds_solute)
			call pmfProtein_Translate(protein, nat_solute, xtop, bnd,nbonds_solute)

			protein%name = adjustl(trim(title))
			write(*,'(a,i5,a,i5,a)') 'Translated ', protein%num_atom, ' atoms and ', protein%num_bond, ' bonds to protein'

			! Assign atom types using protein rules
			call pmfValue_Molecule(protein,protRule,prot_num_rule,prot=1)
			if(input%show_protein.eq.'YES') call pmfShow_Molecule(protein)
			protein%atom(:)%nonzero_score = 0

	! ==============================
	! 3. SCORE INITIAL CONFIGURATION
	! ==============================

		if(bDoTopcalc.eq.0) return

		write(*,'(a)') 'Scoring initial configuration'
		score = PMF_Score(ligand, protein,input)
		call pmfScoring_Stats(ligand)

end subroutine pmf_precalc

subroutine pmfScoring_Stats(mol,iAtom)
	type(tPMFMolecule)	:: mol
	integer,optional		:: iAtom
	integer							:: i,j
	real								:: sum

	sum = 0
	if(present(iAtom)) then
		i = iAtom
		write(*,'(t40,a6, t47,i4, t52,a1,t53,a2,t55,a2, t58,f7.2)') &
		'qatom ', mol%atom(i)%id, '(', mol%atom(i)%ttype, '):', mol%atom(i)%abs
	else
		write(*,'(t30,a4, t38,a4, t50,a5)')  'ID', 'TYPE', 'SCORE'
		do i = 0,mol%num_atom -1
			write(*,'(t30,i4, t38,a, t49,f6.2)') mol%atom(i)%id, mol%atom(i)%ttype, mol%atom(i)%abs
			
			if(input%show_detailed_abs.eq.'YES') then
				do j = 1, mol%atom(i)%abscount
					write(*,'(t1,i4,t8,f6.2,t16,f6.2, t24,a,t29,a,t34,i3)') &
					mol%atom(i)%absdetail(j)%atom_id, &
					mol%atom(i)%absdetail(j)%score, &
					mol%atom(i)%absdetail(j)%r, &
					mol%atom(i)%absdetail(j)%k, &
					mol%atom(i)%absdetail(j)%l, &
					mol%atom(i)%absdetail(j)%m
				end do
			end if

			sum = sum + mol%atom(i)%abs
		end do
		write(*,'(t30,a)') '-------------------------'
		write(*,'(a,t48,f7.2)') 'Total score of topology',sum
	end if
end subroutine pmfScoring_Stats

subroutine pmfShow_ABS(mol)
	type(tPMFMolecule)	:: mol
	integer							:: i

	write(*,'(a)') 'id   type qtype  abs'
	do i = 0, mol%num_atom -1
		write(*,'(t1,i4, t6,a2, t11,a, t18,f7.3)') mol%atom(i)%id, mol%atom(i)%ttype, mol%atom(i)%qtype, mol%atom(i)%abs
	end do
end subroutine pmfShow_ABS

real function PMF_Score(lig,prot,inp)
	type(tPMFMolecule)	:: lig,prot
	type(tPMFInput)			:: inp
	integer							:: i,j,k,l,m
	real								:: r, cutoff, resolution, score, sum

	resolution = 0.2

	sum = 0
	do i = 0, lig%num_atom -1
		if(lig%atom(i)%valid.eq.0) cycle
		if(lig%atom(i)%mask.eq.0)  cycle
		do j = 0,prot%num_atom -1
			if(prot%atom(j)%valid.eq.0) cycle
			if(prot%atom(j)%mask.eq.0)  cycle

			if((lig%atom(i)%ttype(1:1).eq.'C' .or. &
			lig%atom(i)%ttype(1:1).eq.'c') .and. &
			(prot%atom(j)%ttype(1:1).eq.'C' .or. &
			prot%atom(j)%ttype(1:1).eq.'c')) then
				cutoff = 6.0
			else
				cutoff = 9.0
			end if

			r = pmf_distance(lig%atom(i)%coor,prot%atom(j)%coor)
			if(r>cutoff) cycle
			
			k = pmfGetIndex( lig%atom(i))													! ligand atom type
			l = pmfGetIndex(prot%atom(j))													! protein atom type
			m = int(r/resolution)																	! bin
						
			if(m<1) then
				write(*,'(a,i6,a,a,a,a,a,i6,a,a,a,a,a,f6.3,a)') &
				'WARNING: Ligand atom ', i+1, ' (', &
				trim(adjustl(lig%atom(i)%qtype_long)),',', &
				trim(adjustl(lig%atom(i)%residue)), &
				') and protein atom ', j+1, ' (', &
				trim(adjustl(prot%atom(j)%qtype_long)),',', &
				trim(adjustl(prot%atom(j)%residue)), &
				') are too close (r = ', r,' A). Using r = 0.2 A.'
				warn = warn +1
				m = 1
				cycle
			end if
			
			score = pmfdata(k,l,m)
			lig%atom(i)%abs  = lig%atom(i)%abs + score
			if(abs(score)>0.000001) lig%atom(i)%nonzero_score = 1

			if(input%show_detailed_abs.eq.'YES') then
				lig%atom(i)%abscount = lig%atom(i)%abscount +1						! save detailed score info
				lig%atom(i)%absdetail(lig%atom(i)%abscount)%score		= score
				lig%atom(i)%absdetail(lig%atom(i)%abscount)%atom_id = prot%atom(j)%id
				lig%atom(i)%absdetail(lig%atom(i)%abscount)%r				= r

				lig%atom(i)%absdetail(lig%atom(i)%abscount)%k				= pmftype(1,k)
				lig%atom(i)%absdetail(lig%atom(i)%abscount)%l				= pmftype(2,l)
				lig%atom(i)%absdetail(lig%atom(i)%abscount)%m				= m
			end if

			prot%atom(j)%abs = prot%atom(j)%abs + score
			sum = sum + score
		end do
	end do

	lig%score = sum
	prot%score = sum
	PMF_Score = sum
end function PMF_Score

integer function pmfGetIndex(atom)
	type(tPMFAtom)	:: atom
	integer					:: i

	do i = 1,34
		if(atom%ttype.eq.pmftype(atom%origin,i)) then
			pmfGetIndex = i
			return
		end if
	end do

	write(*,'(a,i6,a,a)') 'FATAL: atom ', atom%id, ' has an invalid type specifier: ', atom%ttype
	stop
	pmfGetIndex = 1
end function pmfGetIndex

real function pmf_distance(a, b)
	real(8),dimension(0:2)	:: a,b
	real									:: d,tmpx,tmpy,tmpz

	tmpx = (a(0)-b(0))*(a(0)-b(0))
	tmpy = (a(1)-b(1))*(a(1)-b(1))
	tmpz = (a(2)-b(2))*(a(2)-b(2))

	d = sqrt(tmpx+tmpy+tmpz)

	pmf_distance = d
end function pmf_distance

subroutine pmfprotein_translate(protein,nAtoms,coordinates,bond,num_bond,transfer_waters)
	type(tPMFMolecule)	:: protein
	integer							:: nAtoms
	real(8)							:: coordinates(:)
	type(BOND_TYPE)			:: bond(:)
	integer,optional		:: transfer_waters
	type(tPMFBond)			:: tmp
	integer							:: num_bond
	integer							:: i,j,mark
	integer							:: iRes					! current residue
	integer							:: iNextRes			! starting point of next residue
	integer							:: atom_res_id	! atom number within residue

	! First find out how many atoms to translate
	iRes = 0
	iNextRes = 1		 
	j = 0
	do i = 1,nAtoms
		if(i.eq.iNextRes) then
			iRes = iRes +1
			iNextRes = res(iRes +1)%start
			atom_res_id = 0
		end if
		
		if(iqatom(i).ne.0) cycle																	! skip ligand atoms
		if(res(iRes)%name.eq.'DUM') cycle													! skip dummy atoms
		if((tac(iac(i)).eq.'Du').or.(tac(iac(i)).eq.'DU')) cycle	! skip dummy
		j = j +1
	end do

	protein%num_atom = j
	allocate(protein%atom(0:nAtoms -1),stat=err) 
	if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'
	
	! Transfer all solute atom
	iRes = 0
	iNextRes = 1		 
	j = -1
	do i = 1,nAtoms
		if(i.eq.iNextRes) then
			iRes = iRes +1
			iNextRes = res(iRes +1)%start
			atom_res_id = 0
		end if
		atom_res_id = atom_res_id +1

		if(iqatom(i).ne.0) cycle																	! skip ligand atoms
		if(res(iRes)%name.eq.'DUM') cycle													! skip dummy atoms
		if((tac(iac(i)).eq.'Du').or.(tac(iac(i)).eq.'DU')) cycle	! skip dummy
		j = j +1
		top2prot(i) = j																						! set topology -> protein translation
		
		protein%atom(j)%topindex		= i								! save topology index of this atom
		write(protein%atom(j)%res_id, '(I8)') iRes		! use internal write to convert integer -> character
		protein%atom(j)%coor(0)			= coordinates(3*(i-1) +1)
		protein%atom(j)%coor(1)			= coordinates(3*(i-1) +2)
		protein%atom(j)%coor(2)			= coordinates(3*(i-1) +3)
		protein%atom(j)%id					= j +1
		protein%atom(j)%atom_res_id = atom_res_id
		protein%atom(j)%residue			= res(iRes)%name
		protein%atom(j)%qtype				= tac(iac(i))
		protein%atom(j)%ttype				= ' .'								
		protein%atom(j)%hb					= ''
		protein%atom(j)%weight		  = iaclib(iac(i))%mass
		protein%atom(j)%origin			= 2		! This is a protein atom
		protein%atom(j)%part				= 1		! regular atom
		protein%atom(j)%valid				= 1
		protein%atom(j)%mask				= 0
		if(masks(1)%mask(i)) protein%atom(j)%mask = 1 !masks(1)%mask(i+1)
		protein%atom(j)%neib(:)			= 0
		protein%atom(j)%bond(:)			= 0
		protein%atom(j)%num_bond		= 0
	end do

	! Count protein bonds
	protein%num_bond = 0
	do i = 1,num_bond
		if(((bond(i)%i<=nAtoms).and.(bond(i)%j<=nAtoms)).and. &
		(iqatom(bond(i)%i).eq.0.and.iqatom(bond(i)%j).eq.0)) protein%num_bond = protein%num_bond +1
	end do
	
	allocate(protein%bond(0:protein%num_bond-1))
	j = 0
	do i = 1, num_bond
		if(((bond(i)%i<=nAtoms).and.(bond(i)%j<=nAtoms)).and.(iqatom(bond(i)%i).eq.0.and.iqatom(bond(i)%j).eq.0)) then
			protein%bond(j)%id = j +1
			protein%bond(j)%valid  = 1
			protein%bond(j)%atom_1 = top2prot(bond(i)%i) +1
			protein%bond(j)%atom_2 = top2prot(bond(i)%j) +1
			protein%bond(j)%ttype = SYBYL_bond_type(bond(i)%cod)
			j = j +1
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

!	call PMFShow_Molecule(protein)

end subroutine pmfprotein_translate

subroutine pmfLigand_Translate(ligand,nAtoms,nBonds,aQ,aB,offset)
	! translate list of q-atoms and bonds to xscore-type ligand
	type(tPMFMolecule),intent(out)				:: ligand
	integer,intent(in)										:: nAtoms, nBonds
	type(q_atom),dimension(:),intent(in)	:: aQ
	type(q_bond),dimension(:),intent(in)	:: aB
	integer																:: offset
	integer																:: i, mark,iRes,atom_res_id
	type(tPMFBond)												:: tmp

	! set number of qatoms
	ligand%num_atom = nAtoms
	allocate(ligand%atom(0:nAtoms -1),stat=err)		! old c++ code req. 0:...
	if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'
	ligand%num_bond = nBonds
	allocate(ligand%bond(0:nBonds -1),stat=err)		! old c++ code req. 0:...
	if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'

	!Set residue offset 
	if (offset .eq. 0) then
		do i=1,nres_solute
			if (iqseq(1) <= res(i)%start) then
				offset = i
				exit
			end if
		end do
	end if
	
	iRes=offset	
	atom_res_id=0

	do i = 0,nAtoms-1
		allocate(ligand%atom(i)%absdetail(1:256))			! worst case
		ligand%atom(i)%topindex		= iqseq(i+1)
		ligand%atom(i)%absdetail(:)%score   = 0
		ligand%atom(i)%absdetail(:)%atom_id = 0
		ligand%atom(i)%abscount							= 0

		ligand%atom(i)%id				 = i +1
!		ligand%atom(i)%name			 = ''
		ligand%atom(i)%coor(0)	 = xtop(3*(iqseq(i+1)-1) +1)
		ligand%atom(i)%coor(1)	 = xtop(3*(iqseq(i+1)-1) +2)
		ligand%atom(i)%coor(2)	 = xtop(3*(iqseq(i+1)-1) +3)
		ligand%atom(i)%origin		 = 1	! This is a ligand atom
		ligand%atom(i)%ttype		 = ' .'
		ligand%atom(i)%qtype		 = tac(iac(iqseq(i+1)))
		ligand%atom(i)%weight		 = iaclib(iac(iqseq(i+1)))%mass
		ligand%atom(i)%valid		 = 1
		ligand%atom(i)%mask			 = 0
		if(masks(1)%mask(i+1)) ligand%atom(i)%mask = 1 !masks(1)%mask(i+1)
		ligand%atom(i)%neib(:)	 = 0
		ligand%atom(i)%bond(:)	 = 0
		ligand%atom(i)%num_bond	 = 0

			if(iRes<nres) then
				if(atom_res_id + res(iRes)%start >= res(iRes+1)%start) then
					iRes = iRes +1
					atom_res_id = 0
				end if
			end if
			
			atom_res_id = atom_res_id +1
			ligand%atom(i)%atom_res_id = atom_res_id
			ligand%atom(i)%residue = res(iRes)%name
			write(ligand%atom(i)%res_id, '(I8)') iRes		! use internal write to convert integer -> character
	end do	

	do i = 0,nBonds-1
		ligand%bond(i)%id = i +1			! meanwhile
		ligand%bond(i)%valid  = 1
		ligand%bond(i)%atom_1 = aB(i +1)%a%top_nr - aQ(1)%top_nr +1		! Remove offset so that first q-atom has index 1
		ligand%bond(i)%atom_2 = aB(i +1)%b%top_nr	- aQ(1)%top_nr +1
		ligand%bond(i)%ttype = SYBYL_bond_type(aB(i +1)%cod)
		
		if(adjustl(trim(ligand%bond(i)%ttype)).eq.'') then
			write(*,'(a,i4,a)') 'WARNING: Ligand bond ', i+1, ' is of unknown type. Rings may not be correctly identified.'
			warn = warn +1
		end if
	end do

	! Bubble sort bonds in increasing atom order
	do
		mark = 0
		do i = 0, nBonds -2
			if((ligand%bond(i)%atom_1.eq.ligand%bond(i+1)%atom_1).and. &
				(ligand%bond(i)%atom_2.eq.ligand%bond(i+1)%atom_2)) goto 2
			if(ligand%bond(i)%atom_1 <  ligand%bond(i+1)%atom_1) goto 2
			if((ligand%bond(i)%atom_1.eq.ligand%bond(i+1)%atom_1) .and. &
				(ligand%bond(i)%atom_2<=ligand%bond(i+1)%atom_2)) goto 2

			! If no rule triggered, the bonds are in the wrong order
			mark = 1
			tmp = ligand%bond(i)
			ligand%bond(i) = ligand%bond(i+1)
			ligand%bond(i+1) = tmp
2		end do
		if(mark.eq.0) exit
	end do

	! Reset bond id
	do i = 0,nBonds -1
		ligand%bond(i)%id = i +1
	end do

end subroutine pmfligand_translate

subroutine pmfFixAtomTypes(mol)
	type(tPMFMolecule)	:: mol
	integer							:: i,j,mark,mol2mark
	character(len=2)		:: qtype, tmp

	do i = 0, mol%num_atom -1
		mark = 0
		! set atom MOL2 type
		mol2mark = 0
		do j = 1,size(atomtypetable)
			if(mol%atom(i)%qtype.eq.atomtypetable(j)%other) then
				mol%atom(i)%mol2type = atomtypetable(j)%mol2
				mol2mark = 1
				exit
			end if
		end do
		if(mol2mark.eq.0) then
			write(*,'(a,i5,a,a,a,a,a)') &
			'WARNING: No MOL2 type specified for atom ', &
			mol%atom(i)%id, ' (type ', &
			trim(adjustl(mol%atom(i)%qtype)), &
			') in molecule ', trim(adjustl(mol%name)), &
			'. PMF-type will not be determined correctly.'
			warn = 1
		end if

		! fix qtype to match PMF typing routine
		mol%atom(i)%qtype_long = mol%atom(i)%qtype						! backup qtype
		do j = 1,size(pmf_names)
			if(mol%atom(i)%qtype.eq.pmf_names(j)%qname) then	
				mol%atom(i)%qtype = pmf_names(j)%pmfname
				mark = 1
				exit
			end if
		end do
		if(mark.eq.0) mol%atom(i)%qtype = mol%atom(i)%qtype(1:1)
	end do
end subroutine pmfFixAtomTypes

subroutine pmfValue_Molecule(mol,rule,count,prot,literun)
	type(tPMFMolecule),intent(inout)	:: mol
	type(tPMFRule),pointer						:: rule(:)
	integer														:: count
	integer,optional									:: prot								! use prot=1 to assign ring character to protein atoms
	integer,optional									:: literun						! use literun=1 to skip connection and ring detection

	integer														:: i,j,k,mark,ret,ignored
	type(tPMFAtom)										:: atom
	real															:: t1,t2

	call pmfFixAtomTypes(mol)																! fix atom types
	call pmfDetect_Connections(mol)					

	mol%atom(:)%ring_count = 0															! clear ring data
	mol%atom(:)%ring_id(1) = 0
	mol%atom(:)%ring_id(2) = 0
	mol%atom(:)%ring_id(3) = 0

	call pmfAssign_Rings(mol,input%maximum_ring_size)

	do i = 0,mol%num_atom -1																! fix bond count
		mol%atom(i)%num_bond = 0
		do j = 0,6
			if(mol%atom(i)%bond(j).ne.0) mol%atom(i)%num_bond = mol%atom(i)%num_bond +1
		end do
	end do

	call pmfAssign_Charge(mol)															! assign positive/negative character
	call pmfAssign_Aromatic(mol,req_ring=1)									! assign aromaticity
	call pmfAssign_Hybr(mol)																! assign hybridizing character
	call pmfAssign_Polarity(mol)														! assign polar/nonpolar character
	call pmfAssign_HBchar(mol)															! assign acceptor/donor property
	call pmfAssign_Planar(mol)															! assign planar ring character

	bSecondPass = 0
	do i = 0, mol%num_atom -1
		do j = 1,count																				! check every rule against this atom
			mark = 1
			do k = 1,rule(j)%count
				ret = pmfEval_Rule(mol%atom(i),mol,rule(j)%rule(k))
				if(ret.eq.0) then
					mark = 0
					exit
				end if
			end do
			if(mark.eq.1) then
				mol%atom(i)%ttype = rule(j)%ttype
				exit
			end if
		end do
	end do

	bSecondPass = 1																					! some atoms require a second pass
	do i = 0, mol%num_atom -1
		if(mol%atom(i)%bSecondPass.eq.0) cycle
		do j = 1,count																				! check every rule against this atom
			mark = 1
			do k = 1,rule(j)%count
				ret = pmfEval_Rule(mol%atom(i),mol,rule(j)%rule(k))
				if(ret.eq.0) then

					mark = 0
					exit
				end if
			end do
			if(mark.eq.1) then
				mol%atom(i)%ttype = rule(j)%ttype
				exit
			end if
		end do
	end do

	j = 0
	do i = 0,mol%num_atom -1
		if(mol%atom(i)%ttype.eq.' .') then
			mol%atom(i)%valid = 0
			write(*,'(a,i6,a,a,a,a,a)') &
			'WARNING: Atom ', mol%atom(i)%id, ' (', &
			trim(adjustl(mol%atom(i)%qtype)),') in residue ', &
			trim(adjustl(mol%atom(i)%residue)), &
			' does not match any type rule and will be ignored.'
			j = j +1
		end if
	end do

	ignored = j
	if(j>0) then
		write(*,'(a,i4,a,i4,a,a,a)') 'WARNING: ', j, &
		' atoms (out of ', mol%num_atom, ') in ', &
		adjustl(trim(mol%name)), ' are of unknown type.'
		warn = warn +1
	end if

	! Finally, set ignore flag on atoms of unknown type and, if specifed, atoms in water
	where(mol%atom%ttype.eq.' .') mol%atom%valid = 0
	if(input%ignore_waters.eq.'YES') then
		j = 0
		do i = 0,mol%num_atom -1
			if(mol%atom(i)%valid.eq.0) cycle
			if(mol%atom(i)%residue.eq.'WAT' .OR. mol%atom(i)%residue.eq.'HOH') then
				mol%atom%valid = 0
				j = j +1
			end if
		end do
		ignored = ignored +j
	end if
end subroutine pmfValue_Molecule

subroutine pmfAssign_Planar(mol)
	type(tPMFMolecule)	:: mol
	integer							:: i,j,mark,mark2

	! For each ring, determine if it is aromatic or not. Any atom involved in an aromatic ring is aromatic.

	mol%atom(:)%planar = 0
	do i = 1,mol%num_ring
		mark = 1
		mark2 = 0
		do j = 1,size(mol%ring(i)%path)
			if(mol%atom(mol%ring(i)%path(j)-1)%aromatic.eq.1) mark2 = mark2 +1
			if(mol%atom(mol%ring(i)%path(j)-1)%hybr.ne.'sp2') then
				mark = 0
				exit
			end if
		end do

		if(mark.eq.1) mol%ring(i)%planar = 1
		if(mark2.eq.size(mol%ring(i)%path)) mol%ring(i)%aromatic = 1   

		do j = 1,size(mol%ring(i)%path)
			if(mark.eq.1) mol%atom(mol%ring(i)%path(j)-1)%planar = 1   
		end do
	end do
end subroutine pmfAssign_Planar

integer function pmfEval_Rule(atom,mol,r)
	type(tPMFAtom)								:: atom
	type(tPMFMolecule)						:: mol
	character,intent(in)	:: r   !Removed (len=48) to comile on linux
	character(len=48)							:: rule
	integer												:: invert, mark

	rule = trim(adjustl(r))

	invert = 0
	if(rule(1:3).eq.'not') then															! negate rule?
		invert = 1
		rule = rule(5:48)
	end if

	mark = 0

	if(rule.eq.'aromatic') then
		if(atom%aromatic.eq.1) mark = 1
	elseif(rule.eq.'aliphatic') then
		if(atom%aromatic.eq.0) mark = 1

	elseif(rule.eq.'polar') then
		if(atom%polar.eq.1)	mark = 1
	elseif(rule.eq.'nonpolar') then
		if(atom%polar.eq.0) mark = 1

	elseif(rule.eq.'sp3') then
		if(atom%hybr.eq.'sp3') mark = 1
	elseif(rule.eq.'sp2') then
		if(atom%hybr.eq.'sp2') mark = 1
	elseif(rule.eq.'sp') then
		if(atom%hybr.eq.'sp')  mark = 1

	elseif(rule.eq.'hydrogen bond donor') then
		if(atom%hb.eq.'D') mark = 1
	elseif(rule.eq.'hydrogen bond acceptor') then
		if(atom%hb.eq.'A') mark = 1

	elseif(rule.eq.'in a ring') then
		if(atom%ring.eq.1) mark = 1

	elseif(rule.eq.'in water') then
		if(atom%residue.eq.'WAT' .or. atom%residue.eq.'HOH') mark = 1

	elseif(rule(1:9).eq.'bonded to') then
		if(pmfEval_Bonds(atom,mol,rule(11:48)).eq.1) mark = 1

	elseif(rule.eq.'planar') then
		if(atom%planar.eq.1) mark = 1

	elseif(rule(1:10).eq.'in residue') then
		if(atom%residue.eq.rule(12:14)) mark = 1

	else
		if(index(rule,'-').ne.0) then
			if(atom%charge.eq.-1 .and. atom%qtype.eq.rule(1:1)) mark = 1
		elseif(index(rule,'+').ne.0) then
			if(atom%charge.eq.1  .and. atom%qtype.eq.rule(1:1)) mark = 1
		else
			if(len_trim(rule).eq.1) then
				if(atom%charge.eq.0  .and. atom%qtype.eq.rule(1:1)) mark = 1
			elseif(len_trim(rule).eq.2) then
				if(atom%charge.eq.0  .and. atom%qtype.eq.rule(1:2)) mark = 1
			end if
		end if
	end if

	if(invert.eq.1) mark = abs(1-mark)
	pmfEval_Rule = mark
end function

integer function pmfEval_Bonds(atom,mol,r)
	type(tPMFAtom),intent(inout)	:: atom
	type(tPMFMOlecule),intent(in)	:: mol
	character(len=37),intent(in)	:: r
	character(len=37)							:: rule
	integer												:: i,j,int1,int2,p,filestat,count,mark
	character(len=1),dimension(5) :: tmp
	character(len=1)							:: qtype
	
	rule = r

	pmfEval_Bonds = 0

	if(index(rule,'..')>0) then															! bond def contains interval
		p = index(rule,'..')
		read(rule(p-1:p-1),*,iostat=filestat) int1
		read(rule(p+2:p+2),*,iostat=filestat) int2
		rule = adjustl(trim(rule(p+3:13)))
		count = 0
		do i = 0, atom%num_bond -1
			j = atom%neib(i) -1
			if(mol%atom(j)%qtype.eq.rule) then
				count = count +1
			end if
		end do

		if(count>=int1 .and. count<=int2) then
			pmfEval_Bonds = 1
			return
		end if
	elseif(index(rule,'type')>0) then												! bond def using PMF-type (need second pass)
		if(bSecondPass.eq.0) then
			atom%bSecondPass = 1
			pmfEval_Bonds = 0			
			return
		else
			rule = adjustl(trim(rule(5:13)))
			do i = 0, atom%num_bond -1
				j = atom%neib(i) -1
				if(mol%atom(j)%ttype.eq.rule) then
					pmfEval_Bonds = 1
					return
				end if
			end do
		end if
	elseif(index('123456789',rule(1:1)).ne.0) then					! bond count
		read(rule(1:1),*,iostat=filestat) int1
		rule = adjustl(trim(rule(2:13)))

		count = 0
		do i = 0, atom%num_bond -1
			j = atom%neib(i) -1
			if(mol%atom(j)%qtype.eq.rule) then
				count = count +1
			end if
		end do

		if(count.eq.int1) then
			pmfEval_Bonds = 1
			return
		end if
	elseif(rule(1:16).eq.'atoms other than') then						! other than
		count = 1
		rule = adjustl(trim(rule(17:37)))
		tmp(count)  = rule(1:1)
		
		p = index(rule,'or')
		do while(p.ne.0)
			count = count +1
			rule = adjustl(trim(rule(p+2:37)))
			tmp(count)  = rule(1:1)
			p = index(rule,'or')
		end do

		! make sure at least one bond is to an atom of a type not in the tmp(:)-list
		do i = 0, atom%num_bond -1
			qtype = mol%atom(atom%neib(i)-1)%qtype
			
			mark = 0
			do j = 1, count
				if(qtype.eq.tmp(j)) mark = 1
			end do

			if(mark.eq.0) then
				pmfEval_Bonds = 1
				return
			end if
		end do 
	else
		do i = 0, atom%num_bond -1
			j = atom%neib(i) -1
			if(mol%atom(j)%qtype.eq.rule) then
				pmfEval_Bonds = 1
				return
			end if
		end do
	end if

end function

integer function pmfBonds_Single(atom,mol)								! returns number of single bonds
	type(tPMFAtom),intent(in)			:: atom
	type(tPMFMolecule),intent(in)	:: mol
	integer												:: i

	pmfBonds_Single= 0
	do i = 0, atom%num_bond -1
		if(mol%bond(atom%bond(i)-1)%ttype.eq.'1') pmfBonds_Single= pmfBonds_Single +1
	end do
end function pmfBonds_Single

integer function pmfBonds_Double(atom,mol)								! returns number of double bonds
	type(tPMFAtom),intent(in)			:: atom
	type(tPMFMolecule),intent(in)	:: mol
	integer												:: i

	pmfBonds_Double = 0
	do i = 0, atom%num_bond -1
		if(mol%bond(atom%bond(i)-1)%ttype.eq.'2') pmfBonds_Double = pmfBonds_Double +1
	end do
end function pmfBonds_Double

pure integer function pmfBonds_Triple(atom,mol)						! returns 1 if any bond is triple, otherwise 0
	type(tPMFAtom),intent(in)			:: atom
	type(tPMFMolecule),intent(in)	:: mol
	integer												:: i

	pmfBonds_Triple = 0
	do i = 0, atom%num_neib -1
		if(mol%bond(atom%bond(i)-1)%ttype.eq.'3') pmfBonds_Triple = pmfBonds_Triple +1
	end do
end function pmfBonds_Triple

pure integer function pmfBonds_Aromatic(atom,mol)						! returns 1 if any bond is aromatic, otherwise 0
	type(tPMFAtom),intent(in)			:: atom
	type(tPMFMolecule),intent(in)	:: mol
	integer												:: i

	pmfBonds_Aromatic = 0
	do i = 0, atom%num_neib -1
		if(mol%bond(atom%bond(i)-1)%ttype.eq.'ar') pmfBonds_Aromatic = pmfBonds_Aromatic +1
	end do
end function pmfBonds_Aromatic

integer function pmfBonds_Metal(atom,mol)										! return number of bonds to metals
	type(tPMFAtom),intent(in)			:: atom
	type(tPMFMolecule),intent(in)	:: mol
	integer												:: i

	i = 0
	i = i + pmfEval_Rule(atom,mol,'bonded to Li')
	i = i + pmfEval_Rule(atom,mol,'bonded to Na')
	i = i + pmfEval_Rule(atom,mol,'bonded to K')
	i = i + pmfEval_Rule(atom,mol,'bonded to Rb')
	i = i + pmfEval_Rule(atom,mol,'bonded to Cs')
	i = i + pmfEval_Rule(atom,mol,'bonded to Be')
	i = i + pmfEval_Rule(atom,mol,'bonded to Mg')
	i = i + pmfEval_Rule(atom,mol,'bonded to Ca')
	i = i + pmfEval_Rule(atom,mol,'bonded to Zn')
	i = i + pmfEval_Rule(atom,mol,'bonded to Mn')
	i = i + pmfEval_Rule(atom,mol,'bonded to Fe')
	i = i + pmfEval_Rule(atom,mol,'bonded to V')

	pmfBonds_Metal = i
end function pmfBonds_Metal

subroutine pmfAssign_Hybr(mol)
	type(tPMFMolecule)	:: mol
	integer							:: i,j,k,mark,count,bonds2, bonds3, bonds_ar, all_single
	character(len=3)		:: char
	type(tPMFAtom)			:: atom,atom2
	real								:: angle
	integer							:: amide_flag, atom_1,atom_2
	integer							:: hybr_warn

	! Hybridization is determined from atom MOL2 type only.
	! Translation table for "My Force field" -> MOL2 is read on startup and can be set in input file.
	! Hybridization is assigned to carbons and nitrogens only.
	
	do i = 0,mol%num_atom -1
		char = ''
		select case(mol%atom(i)%mol2type)
			case('C.3')
				char = 'sp3'
			case('C.2')
				char = 'sp2'
			case('C.1')
				char = 'sp'
			case('C.ar')
				char = 'sp2'
			case('C.cat')
				char = 'sp2'

			case('N.3')
				char = 'sp3'
			case('N.2')
				char = 'sp2'
			case('N.1')
				char = 'sp'
			case('N.ar')
				char = 'sp2'
			case('N.am')
				char = 'sp2'
			case('N.pl3')
				char = 'sp2'
			case('N.4')
				char = 'sp3'
		end select

		mol%atom(i)%hybr = char
	end do

	return


!	! Hybridization scheme obtained from:
!	! http://www.sdsc.edu/CCMS/Packages/cambridge/pluto/atom_types.html
!	! Note that hybr. is only assigned for C and N.
!
!	hybr_warn = 0
!	if(input%ff.eq.'gromos') then														! apply gromos rules
!		write(*,'(a)') 'Using the Gromos atom hybridization scheme.'
!		do i = 0, mol%num_atom -1
!			char = ''
!			atom = mol%atom(i)
!			bonds2 = pmfBonds_Double(atom,mol)
!			bonds3 = pmfBonds_Triple(atom,mol)
!			bonds_ar = pmfBonds_Aromatic(atom,mol)
!			all_single = 0
!			if(bonds2.eq.0.and.bonds3.eq.0.and.bonds_ar.eq.0) all_single = 1
!			select case(atom%qtype)
!				case('C')
!					if(atom%num_neib.eq.1) then
!						char = 'sp3'
!					else
!						! Sum all X-N-X angles and calculate mean.
!						angle = 0
!						count = 0
!						do j = 0,atom%num_neib-1
!							do k = j+1,atom%num_neib-1
!								angle = angle + pmfAngle(mol%atom(atom%neib(j)-1)%coor, atom%coor, mol%atom(atom%neib(k)-1)%coor)
!								count = count +1
!							end do
!						end do
!						angle = angle / count
!						
!						if(angle>(109.5-input%hybr_angle_tol) .and. angle<=input%sp2sp3angle) then
!							char = 'sp3'
!						elseif(angle>input%sp2sp3angle .and. angle<(120+input%hybr_angle_tol)) then
!							char = 'sp2'
!						elseif(angle>(180-input%hybr_angle_tol) .and. angle<(180-input%hybr_angle_tol)) then
!							char = 'sp'
!						else
!							char = '-?-'
!							write(*,'(a,i6,a,a,a,f6.2)') 'WARNING: Unable to determine hybridization of atom ', atom%id, ' (type ', trim(adjustl(atom%qtype)), '). Mean of bond angles = ', angle
!							hybr_warn = 1
!						end if
!					end if
!
!
!
!					mark = 0
!					if(all_single.eq.1) then
!						char = 'sp3'
!					elseif(bonds2.eq.1) then
!						char = 'sp2'
!					elseif(bonds2.eq.2) then
!						char = 'sp'
!
!					elseif(atom%num_neib.eq.3 .AND. all_single.eq.1) then
!						char = 'sp3'
!					elseif((atom%num_neib.eq.1 .OR.  atom%num_neib.eq.2) .AND. pmfBonds_Triple(atom,mol).eq.1) then
!						char = 'sp'
!					else
!						char = 'sp2'
!					end if
!
!				case('N')
!					! Nitrogen is a complicated case in Gromos so hybridization is determined from
!					! mean bond angles only. The limit between sp2 and sp3 angles can be set
!					! in input file and defaults to 115 degrees.
!					 
!					if(atom%num_neib.eq.1) then
!						char = 'sp3'
!					else
!						! Sum all X-N-X angles and calculate mean.
!						angle = 0
!						count = 0
!						do j = 0,atom%num_neib-1
!							do k = j+1,atom%num_neib-1
!								angle = angle + pmfAngle(mol%atom(atom%neib(j)-1)%coor, atom%coor, mol%atom(atom%neib(k)-1)%coor)
!								count = count +1
!							end do
!						end do
!						angle = angle / count
!						
!						if(angle>(109.5-input%hybr_angle_tol) .and. angle<=input%sp2sp3angle) then
!							char = 'sp3'
!						elseif(angle>input%sp2sp3angle .and. angle<(120+input%hybr_angle_tol)) then
!							char = 'sp2'
!						elseif(angle>(180-input%hybr_angle_tol) .and. angle<(180-input%hybr_angle_tol)) then
!							char = 'sp'
!						else
!							char = '-?-'
!							write(*,'(a,i6,a,a,a,f6.2)') 'WARNING: Unable to determine hybridization of atom ', atom%id, ' (type ', trim(adjustl(atom%qtype)), '). Mean of bond angles = ', angle
!							hybr_warn = 1
!						end if
!					end if
!			end select
!			mol%atom(i)%hybr = char
!		end do
!	else
!		do i = 0, mol%num_atom -1
!			char = ''
!			atom = mol%atom(i)
!			bonds2 = pmfBonds_Double(atom,mol)
!			bonds3 = pmfBonds_Triple(atom,mol)
!			bonds_ar = pmfBonds_Aromatic(atom,mol)
!			all_single = 0
!			amide_flag = 0
!
!			select case(atom%qtype)
!				case('C')
!					mark = 0
!					if(atom%num_neib.ge.4 .AND. all_single.eq.1) then
!						char = 'sp3'
!					elseif(atom%aromatic.eq.1) then
!						char = 'sp2'
!					elseif((atom%num_bond.eq.1 .OR.  atom%num_bond.eq.2) .AND. pmfBonds_Triple(atom,mol).eq.1) then
!						char = 'sp'
!					else
!						char = 'sp2'
!					end if
!			
!				case('N')
!					count = atom%num_neib - pmfBonds_Metal(atom,mol)						! number of bonds to non metals
!					! determine if this nitrogen is bonded to C=O or C=S
!					if(pmfEval_Rule(atom,mol,'bonded to C').eq.1) then
!						do j = 0, atom%num_neib -1
!							if(mol%atom(atom%neib(j)-1)%qtype.ne.'C') cycle
!							if(pmfBonds_Double(mol%atom(atom%neib(j)-1),mol).ge.1 .and. (pmfEval_Rule(mol%atom(atom%neib(j)-1),mol,'bonded to O') .or. pmfEval_Rule(mol%atom(atom%neib(j)-1),mol,'bonded to S'))) then
!								! this carbon is bonded to an O or an S and has at least one double bond
!								atom2 = mol%atom(atom%neib(j)-1)
!								do k = 0,atom2%num_bond -1
!									if(mol%bond(atom2%bond(k)-1)%ttype.ne.'2') cycle
!									atom_1 = mol%bond(atom2%bond(k)-1)%atom_1
!									atom_2 = mol%bond(atom2%bond(k)-1)%atom_2
!									if(mol%atom(atom_1-1)%qtype.eq.'O' .or. mol%atom(atom_2-1)%qtype.eq.'O') then
!										amide_flag = 1
!										exit
!									elseif(mol%atom(atom_1-1)%qtype.eq.'S' .or. mol%atom(atom_2-1)%qtype.eq.'S') then
!										amide_flag = 1
!										exit
!									end if
!								end do
!								if(amide_flag.eq.1) exit
!							end if
!						end do
!					end if						
!		
!					if(amide_flag.eq.1) write(*,*) 'amide_flag=', amide_flag
!
!					if(count.eq.4 .AND. all_single.eq.1) then
!						char = 'sp3'
!					elseif(atom%aromatic.eq.1) then
!						char = 'sp2'
!					elseif(count.eq.1 .AND. pmfBonds_Triple(atom,mol).eq.1) then
!						char = 'sp'
!					elseif(count.eq.2 .AND. (pmfBonds_Double(atom,mol).eq.2 .OR. (pmfBonds_Single(atom,mol).eq.1 .AND. pmfBonds_Triple(atom,mol).eq.1))) then
!						char = 'sp'
!					elseif(count.eq.3 .AND. amide_flag.eq.1) then
!						char = 'sp2'
!					elseif(count.eq.3) then
!						if(all_single.eq.0) then
!							char = 'sp2'
!						elseif(all_single.eq.1) then
!							if(pmfEval_Rule(atom,mol,'bonded to H').eq.1) then
!								! If one single bond is to an atom that forms a bond of type double, triple, aromatic or 
!								! delocalised .AND. one other single bond is to H then atom_type is N.pl3
!								do j = 0,atom%num_neib-1
!									atom2 = mol%atom(atom%neib(j)-1)
!									if(pmfBonds_Double(atom2,mol).ge.1 .or. pmfBonds_Triple(atom2,mol).ge.1 .or. pmfBonds_Aromatic(atom2,mol).ge.1) then
!										char = 'sp2'
!										exit
!									end if
!								end do
!							end if
!							 
!							if(char.eq.'') then
!								! If one single bond is to an atom that forms a bond of type double, triple, aromatic or 
!								! delocalised .AND. neither of the other single bonds are to H .AND. sum_of_angles around 
!								! N .ge. 350 deg then atom_type is N.pl3									 
!								if(pmfEval_Rule(atom,mol,'not bonded to H').eq.1) then
!									do j = 0,atom%num_neib-1
!										atom2 = mol%atom(atom%neib(j)-1)
!										if(pmfBonds_Double(atom2,mol).ge.1 .or. pmfBonds_Triple(atom2,mol).ge.1 .or. pmfBonds_Aromatic(atom2,mol).ge.1) then
!											char = 'sp2'
!											exit
!										end if
!									end do
!								end if
!							end if
!						end if
!						
!						if(char.eq.'') char = 'sp3'
!					else
!						char = 'sp2'
!					end if
!			end select
!			mol%atom(i)%hybr = char
!		end do
!	end if
!
!	if(hybr_warn.eq.1) warn = warn +1
end subroutine pmfAssign_Hybr

subroutine pmfAssign_Polarity(mol)
	type(tPMFMolecule)	:: mol
	integer							:: i,j

	do i = 0, mol%num_atom -1
		if(mol%atom(i)%qtype.ne.'C') cycle										! skip non-carbon atoms
		mol%atom(i)%polar = 0
		do j = 0,mol%atom(i)%num_neib -1
			if(mol%atom(mol%atom(i)%neib(j)-1)%qtype.ne.'C' .and. &
				 mol%atom(mol%atom(i)%neib(j)-1)%qtype.ne.'H') then
					mol%atom(i)%polar = 1
					exit
			end if
		end do
	end do
end subroutine pmfAssign_Polarity

subroutine pmfAssign_Rings(mol,l)
	type(tPMFMolecule)				:: mol
	integer										:: l																						! ring size limit
	integer										:: limit
	integer										:: i,j
	integer,pointer						:: path(:)																			! list of visited atoms
	integer										:: count																				! number of visited atoms

	limit = l +1																													! to find an n-membered ring, we need limit = n+1
	allocate(path(1:limit))																								! path(:) holds id of atom

	! for each atom: take a walk to neighbouring atoms, but not too far away
	do i = 0, mol%num_atom -1
		path(:) = 0
		call pmfAtom_Walk(i,mol,path,0,limit)
	end do

	! flag atoms in rings
	mol%num_ring = size(mol%ring)
	mol%atom(:)%ring_count = 0
	do i = 1,mol%num_ring
		do j = 1, size(mol%ring(i)%path)
			mol%atom(mol%ring(i)%path(j)-1)%ring_count = mol%atom(mol%ring(i)%path(j)-1)%ring_count +1
			if(mol%atom(mol%ring(i)%path(j)-1)%ring_count<=3) then
				mol%atom(mol%ring(i)%path(j)-1)%ring_id(mol%atom(mol%ring(i)%path(j)-1)%ring_count) = i
			else
				write(*,'(a,i6,a,a,a,a,a)') 'WARNING: Atom ', mol%ring(i)%path(j), ' (', & 
					trim(adjustl(mol%atom(mol%ring(i)%path(j)-1)%qtype)), ' in ', &
					trim(adjustl(mol%atom(mol%ring(i)%path(j)-1)%residue)), &
					') is part of too many rings. Only three rings per atom are considered.'
				warn = warn +1
			end if
			mol%atom(mol%ring(i)%path(j)-1)%ring = 1
		end do
	end do

	! flag bonds between atoms in rings
	where(mol%atom(mol%bond%atom_1-1)%ring.eq.1 .AND. mol%atom(mol%bond%atom_2-1)%ring.eq.1) 
		mol%bond%ring = 1
	elsewhere
		mol%bond%ring = 0
	end where
end subroutine pmfAssign_Rings

recursive subroutine pmfAtom_Walk(start,mol,path,c,limit)
	integer										:: start																				! index of starting atom
	type(tPMFMolecule)				:: mol
	integer,pointer						:: path(:)																			! list of visited atoms
	integer										:: c
	integer										:: limit																				! walk limit
	integer										:: i,j,ret
	integer										:: count,mark

	count = c +1
	path(count) = mol%atom(start)%id																					! path holds id of atom, not index
																																						! atom%neib(:) also holds id of atom
	path(count+1:limit) = 0																										! clear the path ahead

	! if more than 2 steps have been taken and current atom .eq. first atom, we've found a ring
	if(count>2 .AND. path(1).eq.path(count)) call pmfReport_Ring(mol,path,count)

	if(count.eq.limit) return

	! for every neighbour, go for a walk
	do i = 0,mol%atom(start)%num_neib-1
		if(count>1) then
			! this is not the first atom, make sure not to go back where we came from
			if(mol%atom(mol%atom(start)%neib(i)-1)%id.ne.path(count-1)) then				! do not go back where we came from

				! make sure not to go to an atom already in the path, unless it is the first atom
				mark = 1
				do j = 2,count
					if(mol%atom(start)%neib(i).eq.path(j)) then
						mark = 0
						exit
					end if
				end do
					
				if(mark.eq.1) call pmfAtom_Walk(mol%atom(start)%neib(i)-1,mol,path,count,limit)
			end if
		else
			! this is the first atom, no need to check where we came from
			call pmfAtom_Walk(mol%atom(start)%neib(i)-1,mol,path,count,limit)
		end if
	end do
	
end subroutine pmfAtom_Walk

subroutine pmfReport_Ring(mol,path,c)
	type(tPMFMolecule)				:: mol
	integer,pointer						:: path(:)
	integer										:: c,count,s
	integer										:: i,j,k,mark,mark1,mark2
	type(tPMFRing),pointer		:: tmp(:)
	

	count = c -1

	! see if this ring has already been reported
	if(associated(mol%ring)) then
		mark = 0
		do i = 1, size(mol%ring)
			if(count.ne.size(mol%ring(i)%path)) cycle

			! if an element of path(:) is not in mol%ring(i)%path(:) then ring has not been reported
			mark1 = 1
			do j = 1,count
				mark2 = 0
				do k = 1,count
					if(path(j).eq.mol%ring(i)%path(k)) then
						mark2 = 1
						exit
					end if
				end do
				if(mark2.eq.0) then
					mark1 = 0
					exit
				end if
			end do

			! if mark1.eq.1 then found ring matches mol%ring(i)
			if(mark1.eq.1) then
				mark = 1
				exit
			end if
		end do

		if(mark.eq.0) then
			! found ring does not match any of the reported rings, add it
			s = size(mol%ring)
			allocate(tmp(s+1))
			do i = 1,s
				tmp(i)%path	=> mol%ring(i)%path																			! copy reported rings to tmp
			end do
			deallocate(mol%ring)
			mol%ring => tmp																												! update pointer
			
			allocate(mol%ring(s+1)%path(1:count))
			mol%ring(size(mol%ring))%path(1:count) = path(1:count)								! save found ring
		end if
	else
		! first ring
		allocate(mol%ring(1:1))
		allocate(mol%ring(1)%path(1:count))
		mol%ring(1)%path(:) = path(:)
	end if
end subroutine pmfReport_Ring

subroutine pmfInt_PushBack(array, int)
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
end subroutine pmfInt_PushBack		

integer function pmfMolecule_Look_For_A_Ring(molecule, bond_id, atom_path, bond_path, required_size)
! this is a bond-based algorithm to detect rings in a molecule%
! it returns the size of the ring (if zero, the given bond is not in a ring)
! also the corresponding atom path and bond path
	type(tPMFMolecule),intent(inout)	:: molecule
	integer														:: bond_id
	integer,pointer										:: atom_path(:)
	integer,pointer										:: bond_path(:)
	integer														:: required_size,i,ring_size=0,p,next 
	integer,allocatable								::  choice(:,:)

	! initialize the variables first
	allocate(choice(0:molecule%num_bond-1,0:MAX_BOND_NEIB-1),stat=err)
	if(err.ne.0) write(*,'(a)') 'Fatal: Out of memory.'
	choice(:,:) = 0

	! note that atom_path() and bond_path() is not cleaned here
	! they are supposed to be cleaned at the upper level 
	do i = 0,molecule%num_bond -1  
		 if(molecule%bond(i)%valid<=0) goto 1 
		 call pmfMolecule_Reset_Choices(molecule,molecule%bond(i)%id,choice(i,:)) 
1	end do

	! put the given bond at the beginning of the searching chain 
	p=0  
	bond_path(0)=bond_id
	atom_path(0)=molecule%bond(bond_id-1)%atom_1		! bond_id -1

	do
		! search the next possible step 
		do
			next = pmfMolecule_Get_Choices(molecule, choice(bond_path(p)-1,:),atom_path(p),p+1,bond_path) 
			if(next.eq.0) then 	! next step is not available 
				if(p.eq.0) then	! the given bond is not in a ring%
					ring_size=0 
					goto 10
				end if

				! trace back to the previous step 
				call pmfMolecule_Reset_Choices(molecule,bond_path(p),choice(bond_path(p)-1,:)) 
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
		bond_path(p)= next 
		atom_path(p)= pmfMolecule_Two_Bonds_Connection_Check(molecule%bond(bond_path(p-1)-1), molecule%bond(bond_path(p)-1))

		! record the searching history, prevent repeating  
		call pmfMolecule_Clean_Choices(molecule,bond_path(p),choice(bond_path(p-1)-1,:)) 

		! now check if a ring has formed
		if(pmfMolecule_Two_Bonds_Connection_Check(molecule%bond(bond_path(0)-1),molecule%bond(bond_path(p)-1)).eq.atom_path(0)) then
			if(required_size<=0) then ! a ring is found
				ring_size=p+1
				goto 10
      elseif((p+1).eq.required_size) then  ! a ring with required size
				ring_size=p+1
				goto 10	 
			else   ! pre-mature ring, trace back
				call pmfMolecule_Reset_Choices(molecule,bond_path(p), choice(bond_path(p)-1,:)) 
				bond_path(p)=0 
				atom_path(p)=0 
				p = p -1 
				goto 3 
			end if
		else   ! no ring is formed at this step
			if(required_size<=0) then
				goto 3 
			elseif((p+1).eq.required_size) then ! no ring on this path, trace back
				call pmfMolecule_Reset_Choices(molecule,bond_path(p), choice(bond_path(p)-1,:)) 
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
	pmfMolecule_Look_For_A_Ring = ring_size 
end function pmfMolecule_Look_For_A_Ring

subroutine pmfMolecule_Clean_Choices(molecule,tried_choice,choice)
	type(tPMFMolecule)	:: molecule
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

end subroutine pmfMolecule_Clean_Choices

integer function pmfMolecule_Get_Choices(molecule, choice,wrong_end,current_path_length,current_path)
! notice that wrong_end helps to define the direction of the searching path 
	type(tPMFMolecule)	:: molecule
	integer							:: choice(0:MAX_BOND_NEIB -1)
	integer							:: wrong_end,current_path_length
	integer,pointer			:: current_path(:)

	integer							::  i,j,tresult,mark 

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

	pmfMolecule_Get_Choices =  tresult 
end function pmfMolecule_Get_Choices


subroutine pmfMolecule_Reset_Choices(molecule, bond_id,choice)
! find all the neighboring bonds for the given bond 
! it is exclusively used by the ring perception algorithm
	type(tPMFMolecule) :: molecule
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
end subroutine pmfMolecule_Reset_Choices


integer function pmfMolecule_Two_Bonds_Connection_Check(bond1, bond2)
	! this function returns the ID of the joint atom
  type(tPMFBond)		:: bond1, bond2
	integer						:: id 

  id=bond1%atom_1 

  if(id.eq.bond2%atom_1) then
		pmfMolecule_Two_Bonds_Connection_Check = id
		return 
  elseif(id.eq.bond2%atom_2) then
		pmfMolecule_Two_Bonds_Connection_Check = id 
		return
	end if

  id=bond1%atom_2 

  if(id==bond2%atom_1) then
		pmfMolecule_Two_Bonds_Connection_Check = id 
		return
  elseif(id==bond2%atom_2) then
		pmfMolecule_Two_Bonds_Connection_Check = id 
		return
	end if

  pmfMolecule_Two_Bonds_Connection_Check = 0        ! two bonds are not connected
end function pmfMolecule_Two_Bonds_Connection_Check

subroutine pmfAssign_Aromatic(mol,req_ring)
	type(tPMFMolecule)			:: mol
	integer, optional				:: req_ring
	integer									:: i,j

	do i = 0, mol%num_atom -1
		mol%atom(i)%aromatic = 0
		if(present(req_ring).and.req_ring.eq.1 .and. mol%atom(i)%ring.eq.0) cycle									! require atom to be in a ring, if required
		do j = 0,mol%atom(i)%num_bond -1
			if(mol%bond(mol%atom(i)%bond(j)-1)%ttype.eq.'ar') then
				mol%atom(i)%aromatic = 1
				exit
			end if
		end do
	end do
end subroutine pmfAssign_Aromatic

subroutine pmfAssign_Charge(mol)
	type(tPMFMolecule)			:: mol
	integer									:: i,j

	mol%atom(:)%charge = 0
	do i = 0,mol%num_atom -1
		if(mol%atom(i)%qtype.eq.'O'.and.mol%atom(i)%num_bond.eq.1) then						! Oxygens with only 1 single bond are negatively charged
			if(mol%bond(mol%atom(i)%bond(0)-1)%ttype.eq.'1') mol%atom(i)%charge = -1
		elseif(mol%atom(i)%qtype.eq.'N'.and.mol%atom(i)%num_bond.eq.4) then				! Nitrogens with 4 bonds are positively charged
			mol%atom(i)%charge = 1
		end if
	end do
end subroutine pmfAssign_Charge

subroutine pmfAssign_HBchar(mol)
	type(tPMFMolecule)			:: mol
	integer									:: i,j, atom
	character(len=1)				:: mark

	! An atom not bonded to hydrogen is an HB-acceptor. An atom bonded to hydrogen is an HB-donor.
	
	if(input%ff.eq.'gromos') then
		do i = 0, mol%num_atom -1
			mark = 'A'
			if(mol%atom(i)%qtype.eq.'H') then
				mark = ''
			elseif(mol%atom(i)%qtype.eq.'C') then
				! assume any carbon having fewer bonds than its hybridization suggests is bonded to a hydrogen
				if(mol%atom(i)%num_neib<4 .AND. mol%atom(i)%hybr.eq.'sp3') then
					mark = 'D'
				elseif(mol%atom(i)%num_neib<3 .AND. mol%atom(i)%hybr.eq.'sp2') then
					mark = 'D'
				else
					! also check if it actually IS bonded to an H
					do j = 0,mol%atom(i)%num_neib -1
						atom = mol%atom(i)%neib(j) -1
						if(mol%atom(atom)%qtype.eq.'H') then
							mark = 'D'
							exit
						end if
					end do
				end if
			else
				do j = 0,mol%atom(i)%num_neib -1
					atom = mol%atom(i)%neib(j) -1
					if(mol%atom(atom)%qtype.eq.'H') then
						mark = 'D'
						exit
					end if
				end do
			end if
			mol%atom(i)%hb = mark
		end do
	else
		do i = 0, mol%num_atom -1
			mark = 'A'
			if(mol%atom(i)%qtype.eq.'H') then
				mark = ''
			else
				do j = 0,mol%atom(i)%num_neib -1
					atom = mol%atom(i)%neib(j) -1
					if(mol%atom(atom)%qtype.eq.'H') then
						mark = 'D'
						exit
					end if
				end do
			end if
			mol%atom(i)%hb = mark
		end do
	end if
end subroutine pmfAssign_HBchar

subroutine pmfShow_Molecule(mol)
	type(tPMFMolecule)	:: mol
	integer							:: i,j,filestat
	character(len=4)		:: buf
	character(len=16)		:: ring,slask,ar,pol,planar,charge,hb,bondtype,mask
	character(len=64)   :: ringpath,neibs,bonds

	write(*,'(a,a)') 'MOLECULE: ', adjustl(trim(mol%name))
	write(*,'(a)') &
	'  ID topoID MASK TYPE QTYPE   RES   CHARGE  HB  POL HYBR AR  RINGS      PLANAR NEIGHBOURING ATOMS        BONDS TO ATOM'

	do i = 0, mol%num_atom -1
		if(input%show_hydrogens.eq.'NO' .and. mol%atom(i)%ttype(1:1).eq.'H') cycle
		if(input%show_masked.eq.'NO'    .and. mol%atom(i)%mask.eq.0)				 cycle
		
		ring   = ''
		ar		 = ''
		pol    = ''
		planar = ''
		charge = ''
		hb		 = ''
		mask   = ''
		do j = 1, min(mol%atom(i)%ring_count,size(mol%atom(i)%ring_id))
			write(buf,'(i3)') mol%atom(i)%ring_id(j)
			ring((j-1)*4+1:j*4) = buf
		end do		

		neibs  = ''
		bonds  = ''
		do j = 0, mol%atom(i)%num_neib-1
			write(buf,'(i4)') mol%atom(i)%neib(j)
			neibs((j)*5+1:(j+1)*5) = buf
		end do

		do j = 0, mol%atom(i)%num_bond-1
			write(buf,'(i4)') mol%atom(i)%bond(j)
			bonds((j)*5+1:(j+1)*5) = buf
		end do

		if(mol%atom(i)%aromatic.eq.1)	ar		 = 'ar'
		if(mol%atom(i)%planar.eq.1)		planar = 'planar'
		if(mol%atom(i)%polar.eq.1)		pol    = 'pol'
		if(mol%atom(i)%charge.eq.-1)	charge = 'neg'
		if(mol%atom(i)%charge.eq.1)		charge = 'pos'
		if(mol%atom(i)%hb.eq.'D')			hb		 = 'don'
		if(mol%atom(i)%hb.eq.'A')			hb		 = 'acc'
		if(mol%atom(i)%mask.eq.1)			mask   = 'vis.'

!		write(charge,*) mol%atom(i)%charge
	
		write(*,'(t1,i4,t7, i4, t13, a, t18,a2, t23,a, t31,a, t37,a,t45,a,t49,a,t53,a,t58,a,t62,a,t73,a, t80,a,t106,a)') i+1, & 
		  mol%atom(i)%topindex, mask, mol%atom(i)%ttype, mol%atom(i)%qtype_long, mol%atom(i)%residue, &
		  adjustl(trim(charge)), adjustl(trim(hb)), adjustl(trim(pol)), mol%atom(i)%hybr, adjustl(trim(ar)), ring,  &
		  adjustl(trim(planar)), neibs, bonds
								
	end do

	write(*,'(a,i6)') 'Total number of bonds in this molecule = ',mol%num_bond
	if(input%show_bonds.eq.'YES') then
		write(*,'(a)') '  ID TYPE    ATOM-ATOM topoID-topoID RING  NEIGHBOURING BONDS'

		do i = 0,mol%num_bond -1
			if(mol%atom(mol%bond(i)%atom_1-1)%mask.eq.0 .OR. mol%atom(mol%bond(i)%atom_1-1)%mask.eq.0) cycle
			ring = ''
			if(mol%bond(i)%ring.eq.1)			ring		 = 'ring'
			if(mol%bond(i)%ttype.eq.'1')  bondtype = 'single'
			if(mol%bond(i)%ttype.eq.'2')  bondtype = 'double'
			if(mol%bond(i)%ttype.eq.'3')  bondtype = 'triple'
			if(mol%bond(i)%ttype.eq.'ar') bondtype = 'ar'
			if(mol%bond(i)%ttype.eq.'am') bondtype = 'amide'

			neibs = ''
			do j = 0, mol%bond(i)%num_neib-1
				write(buf,'(i4)') mol%bond(i)%neib(j)
				neibs((j)*5+1:(j+1)*5) = buf
			end do

			write(*,'(t1,i4,t6,a,  t14,i4,t19,i4, t25,i4,t31,i4, t38,a, t44,a)') &
						i+1, bondtype, mol%bond(i)%atom_1, mol%bond(i)%atom_2, &
						 mol%atom(mol%bond(i)%atom_1-1)%topindex, &
						 mol%atom(mol%bond(i)%atom_2-1)%topindex, ring, neibs
		end do
	end if

	write(*,'(a,i4)') 'Total number of rings in this molecule = ', mol%num_ring
	if(mol%num_ring>0) then
		write(*,'(a)') '  ID SIZE  AR  PLANAR  ATOMS'
		do i = 1, mol%num_ring
			ringpath = ''
			ar       = ''
			planar	 = ''
			do j = 1, size(mol%ring(i)%path)
				write(buf,'(i4)') mol%ring(i)%path(j)
				ringpath((j-1)*5+1:j*5) = buf
			end do
			if(mol%ring(i)%aromatic.eq.1) ar		 = 'ar'
			if(mol%ring(i)%planar.eq.1)		planar = 'planar'

			write(*,'(t1,i4,t7,i2,t12,a,t16,a,t24,a)') i, size(mol%ring(i)%path),ar, planar,ringpath
		end do
	end if

end subroutine pmfShow_Molecule

subroutine pmfsort_bonds		!identify pairs of atoms involved in q-bonds and/or hydrogen bonding !in start 
	integer						::	atom1,atom2
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
end subroutine pmfsort_bonds

subroutine pmf_calc(iCalc, iFrame)		
	!arguments
	integer, intent(in)	:: iCalc				! calculation index
	integer, intent(in)	:: iFrame				! frame index
	real				:: score								! scoring results
	integer			:: i,mark

	! Update coordinates
	do i = 0,nqat -1
		ligand%atom(i)%coor(0) = xin(3*(iqseq(i+1)-1) +1)
		ligand%atom(i)%coor(1) = xin(3*(iqseq(i+1)-1) +2)
		ligand%atom(i)%coor(2) = xin(3*(iqseq(i+1)-1) +3)
	end do

	do i = 0,protein%num_atom -1
		protein%atom(i)%coor(0) = xin(3*(protein%atom(i)%topindex-1) +1)
		protein%atom(i)%coor(1) = xin(3*(protein%atom(i)%topindex-1) +2)
		protein%atom(i)%coor(2) = xin(3*(protein%atom(i)%topindex-1) +3)
	end do
	
	! NOTE: No need to call pmfValue_Molecule again since it is coordinate-independant

	ligand%atom(:)%abs  = 0
	protein%atom(:)%abs = 0
	score = pmf_Score(ligand,protein,input)	

	call pmflog_frame(iFrame,ligand)

	if(input%show_abs.eq.'YES')		call pmfScoring_Stats(ligand)
	if(input%show_abs.eq.'QATOM') call pmfScoring_Stats(ligand,input%qatom-1)

	if(input%show_total.eq.'YES'.and.input%show_abs.ne.'YES') write(*,'(t30,f8.2)') ligand%score
end subroutine pmf_calc

subroutine pmflog_frame(iFrame,ligand)				! logs scoring results for frame
	integer, intent(in)			:: iFrame						! current frame	
	type(tPMFMolecule)			:: ligand
	type(tPMFScore),pointer	:: new_aPMFScore(:)	! tmp pointer
	
	nPMFScores = nPMFScores +1
	if(nPMFScores>maxPMFScores) then
		maxPMFScores = maxPMFScores + 16										! allocate 16 elements a time
		allocate(new_aPMFScore(maxPMFScores))								! allocate new mem
		new_aPMFScore(1:nPMFScores-1) = pmfscores(1:nPMFScores-1)	! copy old
		
		deallocate(pmfscores)																! return old mem to OS
		pmfscores => new_aPMFScore													! update pointer
	end if

	! finally store values
	pmfscores(nPMFScores)%frame	= iFrame									! iFrame may be redundant
	pmfscores(nPMFScores)%score	= ligand%score
end subroutine pmflog_frame

integer function pmf_add(desc)
	!arguments
	character(*),intent(out)		:: desc
	! locals
	character(len=400)			:: chBuf
	character(len=80)				:: chSmallBuf
	integer									:: i, ats

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
		pmf_add = 0
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
		
	call getlin(fep_file, 'Q-atom (FEP) file: ')
	call getlin(input%ff, 'Force field translation key: ')
	call getlin(chInput,  'Scoring parameters: ')

	pmf_add = Nmasks
	write(desc, 2) masks(Nmasks)%included, nat_pro
2	format('PMF calculation using ',i6,' of ', i6, ' atoms.')
end function pmf_add

subroutine pmf_initialize
	nPMFScores = 0
	maxPMFScores = 16
	allocate(pmfscores(maxPMFScores))

	warn = 0
end subroutine pmf_initialize

subroutine pmf_finalize
	integer	:: i

	if(waterwarn.eq.1) then
		write(*,'(a)') 'WARNING: Water molecules were static in one or more frames.'
		warn = warn +1
	end if

	do i = 0,ligand%num_atom -1
		if(ligand%atom(i)%nonzero_score.eq.0) then
			write(*,'(a,i4,a,a,a,a,a)') 'WARNING: Ligand atom ', ligand%atom(i)%id, ' (', trim(adjustl(ligand%atom(i)%ttype)), ' in residue ', trim(adjustl(ligand%atom(i)%residue)), ') scored zero in all frames.'
			warn = warn +1
		end if
	end do

	if(warn.ne.0) then
		write(*,'(a,i3,a)') 'WARNING: There were ', warn, ' warnings.'
	else
		write(*,'(a)') 'PMF-Score finished succesfully.'
	end if
end subroutine pmf_finalize

subroutine pmf_heading(i)
	integer, intent(in)		:: i
	write(*,'(t33,a5)') 'SCORE'
end subroutine pmf_heading

subroutine pmfscore_mean
	real(8)					:: score
	integer					:: i
	
	score = 0
	do i = 1,nPMFScores
		score = score + pmfscores(i)%score
	end do

	score = score / nPMFScores

	! output mean value
	write(*,'(t1,a, t19,i5, t33,a5)')  'processed frames: ', nPMFScores, 'SCORE'
	write(*,'(a, t30,f8.2)') 'mean score',score

	pmfscores(:)%score  = 0
	pmfscores(:)%frame = 0
	nPMFScores   = 0
end subroutine pmfscore_mean

subroutine pmfDetect_Connections(molecule)
	type(tPMFMolecule),intent(inout)	:: molecule
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

	do i = 0,molecule%num_atom -1
		part_list(i)=0
	end do

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
	
	do i = 0,molecule%num_bond-1 -1 
		do j=i+1,molecule%num_bond  -1
			 if((molecule%bond(i)%valid<=0).or.(molecule%bond(j)%valid<=0)) goto 8 

			 if(pmfMolecule_Two_Bonds_Connection_Check(molecule%bond(i),molecule%bond(j)).ne.0) then
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
end subroutine pmfDetect_Connections

real function pmfAngle_of_Two_Vectors(v1, v2)
	real, dimension(0:2)	:: v1,v2
  real angle
  real l1,l2,tmp1,tmp2

	if((v1(0).eq.v2(0)).and.(v1(1).eq.v2(1)).and.(v1(2).eq.v2(2))) then					! need to treat this case
		pmfAngle_of_Two_Vectors = 0
		return
	end if

	l1=sqrt(v1(0)*v1(0)+v1(1)*v1(1)+v1(2)*v1(2))
	l2=sqrt(v2(0)*v2(0)+v2(1)*v2(1)+v2(2)*v2(2))

	tmp1=v1(0)*v2(0)+v1(1)*v2(1)+v1(2)*v2(2)
	tmp2=l1*l2

	angle=acos(tmp1/tmp2)
	angle=angle/(3.14159265)*180.0

	pmfAngle_of_Two_Vectors = angle  ! return angle in degree, 0-180
end function pmfAngle_of_Two_Vectors

real function pmfAngle(a,b,c)	
	! the a-b-c angle in degree 0-180
	real,dimension(0:2) :: a,b,c
	integer :: i
	real,dimension(0:2)	:: v1,v2

	do i=0,2
		v1(i)=b(i)-a(i)
		v2(i)=b(i)-c(i)
	end do

	pmfAngle=pmfAngle_of_Two_Vectors(v1,v2)
end function pmfAngle

end module
