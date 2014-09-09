! (C) 2014 Uppsala Molekylmekaniska HB, Uppsala, Sweden
! calc_nb.f90
! by Martin Almlöf & Martin Nervall
!       
! nonbonded potential calculation &
! calculation of nb interactions between 'q-atoms' and
! protein residues

module CALC_NB
	use CALC_BASE
	use MASKMANIP
	use PARSE
	use ATOM_MASK
!	use QATOM
	implicit none
!constants
	integer, parameter					::	MAX_LISTS = 20
	integer, parameter					::  MAX_NB_QP = 1
!module variables
	type(MASK_TYPE), private, target	::	masks(2)
	integer, private					::	Nlists = 0
	integer, private                    ::  N_nb_qp = 0
	type NB_LIST_TYPE
		integer							::  number_of_NB=0
		integer, pointer				::	atom1(:), atom2(:)
		real, pointer					::  AA(:), BB(:), qq(:)
	end type NB_LIST_TYPE
	type NBQ_LIST_TYPE
		integer						:: number_of_NB=0
		integer, pointer				:: atom1(:), atom2(:)
		real, pointer					:: AA(:,:), BB(:,:), qq(:,:)
	end type NBQ_LIST_TYPE
	type(NB_LIST_TYPE), private, allocatable::  nb_listan(:)
	type(NB_LIST_TYPE), private, allocatable::  nb_list_res(:)
	type(NBQ_LIST_TYPE), private, allocatable:: nbq_list_res(:)
	type AVERAGES
		real(8)							:: lj, el
	end type AVERAGES
	type(AVERAGES),allocatable			::  aves(:)
	type(AVERAGES),allocatable			::  qp_aves(:)
	integer,allocatable					::	total_frames(:)
	integer								::  total_qp_frames = 0 

	type NB_QP_TYPE
	integer								:: p_first, p_last,q_first, q_last
	end type NB_QP_TYPE
	type(NB_QP_TYPE), private					:: qp_calc
contains

!*******************************************************************

subroutine nb_initialize
	allocate(nb_listan(MAX_LISTS), total_frames(MAX_LISTS), aves(MAX_LISTS))
	aves(:)%lj = 0
	aves(:)%el = 0
	total_frames(:) = 0
end subroutine nb_initialize

!*******************************************************************

subroutine nb_finalize(i)
	integer						::	i


!	Print averages of all trajectories
	write (*,669) i, aves(i)%lj
	write (*,670) i, aves(i)%el
669	format('Average NB_lj energies of NB_calc ',i2 ,' : ',f9.2)
670	format('Average NB_el energies of NB_calc ',i2 ,' : ',f9.2)

end subroutine nb_finalize

!*******************************************************************

integer function nb_add(desc)
	!arguments
	character(*)				::	desc
	integer						::	ats1
	integer						::	ats2
	if(Nlists == MAX_LISTS) then
		write(*,10) MAX_LISTS/2
		nb_add = 0
		return
	end if
10	format('Sorry, the maximum number of nonbond calculations is ',i2)

	write (*, 52)
52	format('N.B. the calculated nonbonded energies are for the trajectory coordinates using topology parameters.')
	write (*, 53)
53	format('i.e. lambda is not taken into account.')


	!make mask for first set of atoms
	write(*, 11)
11	format('Enter mask for first group of atoms')

	call mask_initialize(masks(1))
	ats1 =  maskmanip_make(masks(1))

	if(ats1 == 0 ) then
		call mask_finalize(masks(1))
		nb_add = 0
		return
	end if

	write(*, 13) masks(1)%included
13	format('First mask contains ',i6,' atoms')

	write(*, 14)
14	format('Enter mask for second group of atoms')

	!make mask for second set of atoms
	call mask_initialize(masks(2))
	ats2 =  maskmanip_make(masks(2))

	!discard if no atoms in mask
	if(ats2 == 0) then
		call mask_finalize(masks(1))
		call mask_finalize(masks(2))
		nb_add = 0
		return
	end if

	write(*, 16) masks(2)%included
16	format('Second mask contains ',i6,' atoms')


	Nlists = Nlists + 1
	allocate(nb_listan(Nlists)%atom1(ats1*ats2))
	allocate(nb_listan(Nlists)%atom2(ats1*ats2))
	allocate(nb_listan(Nlists)%AA(ats1*ats2))
	allocate(nb_listan(Nlists)%BB(ats1*ats2))
	allocate(nb_listan(Nlists)%qq(ats1*ats2))

	call nb_make_list(masks(1), masks(2),nb_listan(Nlists),nbq_list_res(Nlists))

	nb_add = Nlists
	write(desc, 20) nb_listan(Nlists)%number_of_NB
20	format('Nonbonded list contains ',i6,' pairs')
 end function nb_add


!*******************************************************************
!* Invoked from Qcalc to calculate nonbonded interactions 
!*  between two masks defined by user.
!*******************************************************************
subroutine nb_calc(i)
integer :: i
real(8) :: NB_Vlj, NB_Vel

	call nb_calc_lists(NB_Vlj, NB_Vel, nb_listan(i),nbq_list_res(i))
 	write(*,99, advance='no') NB_Vlj
 	write(*,100, advance='no') NB_Vel
99  format(f9.2)
100	format(f8.2)

!	Make averages. To be printed after all trajectory files are processed
	aves(i)%lj = (aves(i)%lj*total_frames(i)+NB_Vlj)/(total_frames(i)+1)
	aves(i)%el = (aves(i)%el*total_frames(i)+NB_Vel)/(total_frames(i)+1)

	total_frames(i) = total_frames(i) + 1
end subroutine nb_calc


!*******************************************************************

subroutine nb_calc_lists(NB_Vlj, NB_Vel, nb_list, nbq_list)
	!arguments
	real(8), intent(out)		::  NB_Vlj, NB_Vel
	type(NB_LIST_TYPE), intent(in) :: nb_list
	type(NBQ_LIST_TYPE),intent(in)  :: nbq_list

	!locals
	real(8)						::	r, x1, x2, y1, y2, z1, z2, invr, invr6, invr12
	integer						::  j,k,length,lengthq,tmpstates

	NB_Vlj = 0
	NB_Vel = 0 

	length = nb_list%number_of_NB
	lengthq = nbq_list%number_of_NB
	if(.not.use_fep) then
	do j=1,length
		x1 = xin(3*nb_list%atom1(j)-2)
		x2 = xin(3*nb_list%atom2(j)-2)
		y1 = xin(3*nb_list%atom1(j)-1)
		y2 = xin(3*nb_list%atom2(j)-1)
		z1 = xin(3*nb_list%atom1(j))
		z2 = xin(3*nb_list%atom2(j))
		r = sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
		if (r /= 0) then
			invr = 1/r
			invr6 = invr*invr*invr*invr*invr*invr
			invr12 = invr6*invr6
			
			
			if(ivdw_rule==1) then !geometric comb. rule
				NB_Vlj = NB_Vlj + nb_list%AA(j)*invr12 - nb_list%BB(j)*invr6
            else !arithmetic
				NB_Vlj = NB_Vlj + sqrt(nb_list%BB(j)) * (nb_list%AA(j))**6 * invr6 * ((nb_list%AA(j))**6 * invr6 - 2.0)
			endif			
								
			NB_Vel = NB_Vel + nb_list%qq(j)*invr
		end if
	end do
	else
        do j=1,lengthq
                x1 = xin(3*nbq_list%atom1(j)-2)
                x2 = xin(3*nbq_list%atom2(j)-2)
                y1 = xin(3*nbq_list%atom1(j)-1)
                y2 = xin(3*nbq_list%atom2(j)-1)
                z1 = xin(3*nbq_list%atom1(j))
                z2 = xin(3*nbq_list%atom2(j))
                r = sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
                if (r /= 0) then
                        invr = 1/r
                        invr6 = invr*invr*invr*invr*invr*invr
                        invr12 = invr6*invr6


                        if(ivdw_rule==1) then !geometric comb. rule
				do tmpstates=1,nstates
                                NB_Vlj = NB_Vlj + (lamda(tmpstates) * (nbq_list%AA(tmpstates,j)*invr12 - nbq_list%BB(tmpstates,j)*invr6))
				end do
		        else !arithmetic
				do tmpstates=1,nstates
                                NB_Vlj = NB_Vlj + (lamda(tmpstates) * (sqrt(nbq_list%BB(tmpstates,j)) * (nbq_list%AA(tmpstates,j))**6 * invr6 * ((nbq_list%AA(tmpstates,j))**6 * invr6 - 2.0)))
				end do
                        endif
			do tmpstates=1,nstates
                        NB_Vel = NB_Vel + (lamda(tmpstates) * (nbq_list%qq(tmpstates,j)*invr))
			end do
                end if
        end do
	end if
end subroutine nb_calc_lists

!*******************************************************************
! Make pair lists of interacting atoms in mask1 and mask2
!*******************************************************************

subroutine nb_make_list(mask1, mask2, nb_list, nbq_list)
	integer						::	j, k, l, m, mq, nat1, nat2, LJ_code, nl, qq_pair, atomj, atomk, b
	integer						::  size_groupm,tmpstates
	integer, allocatable	::	group1(:), group2(:)
	integer						::  ljcod(255,255)
	integer,allocatable				::  groupm(:)
	logical						::	NB
	type(MASK_TYPE), intent(in) :: mask1, mask2
	type(NB_LIST_TYPE)			:: nb_list
	type(NBQ_LIST_TYPE)			:: nbq_list
!new real for scaling factor
        real*8                                  :: fullscale=0

	l = 1
	m = 1
	allocate(groupm(nat_pro))
	!make the list of atoms in first mask
	nat1 = mask1%included 
	allocate(group1(nat1))
	do j = 1, nat_pro
		if (mask1%MASK(j)) then
			group1(l) = j
			l = l+1
		end if
	end do

	!make the list of atoms in second mask
	nat2 = mask2%included 
	allocate(group2(nat2))
	do j = 1, nat_pro
		if (mask2%MASK(j)) then
			group2(m) = j
			m = m+1
		end if
	end do

	!make list of mutual atoms
	b=0
	do j = 1, l-1
		do k = 1, m-1
			if (group1(j) == group2(k)) then
				b = b + 1
				groupm(b) = group1(j)
				exit
			end if
		end do
	end do
	size_groupm = b

	ljcod(:,:) = 1
	do m=1,nlj2
		ljcod(lj2(m)%i, lj2(m)%j) = 2
		ljcod(lj2(m)%j, lj2(m)%i) = 2
	end do

	m=1
        mq=1

	do j = 1, nat1
		do k = 1, nat2
                        atomj=0
                        atomk=0
			atomj = group1(j)
			atomk = group2(k)

			NB = .true.
			!check if we're adding a pair which has already been added
			if (atomj .ge. atomk) then
				do l = 1, size_groupm  !cycle through list of mutual atoms
					if (atomk == groupm(l)) then
						NB = .false.
						exit
					end if
				end do
			end if

			!check if the atoms are in a mutual angle
			!if so, they're 1-2 or 1-3 neighbors and no NB terms should be calculated
			!code from calc_geom.f90
			if (NB) then
				do b = 1, nangles
					if( (ang(b)%i == atomj .and. &
						(ang(b)%j == atomk .or. ang(b)%k == atomk)) .or. &
						(ang(b)%j == atomj .and. &
						(ang(b)%i == atomk .or. ang(b)%k == atomk)) .or. & 
						(ang(b)%k == atomj .and. &
						(ang(b)%i == atomk .or. ang(b)%j == atomk))) then
						NB = .false.
						exit
					end if
				end do
			end if


			!some LJ_CODE code taken from md.f90::nbmonitorlist

			!starting guess = use LJ_code matrix for topology atom types


				if (NB) then
					if (use_fep) then
					nbq_list%atom1(mq) = atomj
					nbq_list%atom2(mq) = atomk
					else
					nb_list%atom1(m) = atomj
					nb_list%atom2(m) = atomk
					end if

					LJ_code = ljcod(iac(atomj),iac(atomk)) !iac is from topo.f90 "integer atom code"

	
					!Are atoms of pair in 1-4 position?
!					if(iqatom(atomj) == 0 .and. iqatom(atomk) == 0) then !neither atom is q_atom
					if (abs(atomk-atomj) .le. max_nbr_range ) then
						if (atomj .gt. atomk ) then
							if ( list14(atomj-atomk, atomk) ) LJ_code = 3 !3 means 1-4
																			!list14 contains 1-4 neighbors
																			!within 25 atom#s of eachother
																			!n14long contains 1-4neighbors
																			!which are "further" than 25
																			!atom#s
						else
							if ( list14(atomk-atomj, atomj) ) LJ_code = 3
						end if
					else
						do nl = 1, n14long
							if ((list14long(1, nl) .eq. atomj &
								.and. list14long(2, nl) .eq. atomk) .or. &
								(list14long(1, nl) .eq. atomk &
								.and. list14long(2, nl) .eq. atomj)) then
								LJ_code = 3
							endif
						end do
					endif	!kolla 1-4 interaktioner

!Only troll iqatom list if use_fep is specified and atomj and atomk .le. nprot
                                        if(ivdw_rule==1) then !geometric comb. rule
						!iaclib() from topo.f90, contains atom parameters.
						!argument is "integer atom code"
                                                if (use_fep) then
                                                        if (atomj.le.nat_solute) then
                                                                if (iqatom(atomj).ne.0) then !atom is a Q atom from the FEP file
                                                                        if (atomk.le.nat_solute) then
                                                                                if (iqatom(atomk).ne.0) then !other atom is also a Q atom
                                                                                        do tmpstates=1,nstates
                                                                                        nbq_list%AA(tmpstates,mq)=qavdw(qiac(iqatom(atomj),tmpstates),LJ_code)*qavdw(qiac(iqatom(atomk),tmpstates),LJ_code)
                                                                                        end do
                                                                                else !other atom is not a Q atom, have to use LJ_code = 1 for Q atom instead of possible 2
                                                                                        do tmpstates=1,nstates
                                                                                        if (LJ_code .eq. 2) then
                                                                                                nbq_list%AA(tmpstates,mq)=qavdw(qiac(iqatom(atomj),tmpstates),1)*iaclib(iac(atomk))%avdw(LJ_code)
                                                                                        else 
                                                                                                nbq_list%AA(tmpstates,mq)=qavdw(qiac(iqatom(atomj),tmpstates),LJ_code)*iaclib(iac(atomk))%avdw(LJ_code)
                                                                                        end if
                                                                                        end do
                                                                                end if
                                                                         else !other atom most likely water, not in original topology, so not a qatom
                                                                                do tmpstates=1,nstates
                                                                                if (LJ_code .eq. 2) then
                                                                                        nbq_list%AA(tmpstates,mq)=qavdw(qiac(iqatom(atomj),tmpstates),1)*iaclib(iac(atomk))%avdw(LJ_code)
                                                                                else 
                                                                                        nbq_list%AA(tmpstates,mq)=qavdw(qiac(iqatom(atomj),tmpstates),LJ_code)*iaclib(iac(atomk))%avdw(LJ_code)
                                                                                end if
                                                                                end do
                                                                        end if
                                                                else
                                                                        if  (atomk.le.nat_solute) then
                                                                                if (iqatom(atomk).ne.0) then !other atom is a Q atom
                                                                                        do tmpstates=1,nstates
                                                                                        if (LJ_code .eq. 2) then
                                                                                                nbq_list%AA(tmpstates,mq)=iaclib(iac(atomj))%avdw(LJ_code)*qavdw(qiac(iqatom(atomk),tmpstates),1)
                                                                                        else 
                                                                                                nbq_list%AA(tmpstates,mq)=iaclib(iac(atomj))%avdw(LJ_code)*qavdw(qiac(iqatom(atomk),tmpstates),LJ_code)
                                                                                        end if
                                                                                        end do
                                                                                else !not a qatom
                                                                                        do tmpstates=1,nstates
                                                                                        nbq_list%AA(tmpstates,mq) = iaclib(iac(atomj))%avdw(LJ_code)*iaclib(iac(atomk))%avdw(LJ_code)
                                                                                        end do
                                                                                end if
                                                                        else !other atom most likely water, not in original topology, so not a qatom
                                                                                do tmpstates=1,nstates
                                                                                nbq_list%AA(tmpstates,mq) = iaclib(iac(atomj))%avdw(LJ_code)*iaclib(iac(atomk))%avdw(LJ_code)
                                                                                end do
                                                                        end if
                                                                end if
                                                        else !atom is not in original topology, not a qatom
                                                                if  (atomk.le.nat_solute) then
                                                                        if (iqatom(atomk).ne.0) then !other atom is a Q atom
                                                                                do tmpstates=1,nstates
                                                                                if (LJ_code .eq. 2) then
                                                                                        nbq_list%AA(tmpstates,mq)=iaclib(iac(atomj))%avdw(LJ_code)*qavdw(qiac(iqatom(atomk),tmpstates),1)
                                                                                else 
                                                                                        nbq_list%AA(tmpstates,mq)=iaclib(iac(atomj))%avdw(LJ_code)*qavdw(qiac(iqatom(atomk),tmpstates),LJ_code)
                                                                                end if
                                                                                end do
                                                                        else !not a qatom
                                                                                do tmpstates=1,nstates
                                                                                nbq_list%AA(tmpstates,mq) = iaclib(iac(atomj))%avdw(LJ_code)*iaclib(iac(atomk))%avdw(LJ_code)
                                                                                end do
                                                                        end if
                                                                else !other atom most likely water, not in original topology, so not a qatom
                                                                        do tmpstates=1,nstates
                                                                        nbq_list%AA(tmpstates,mq) = iaclib(iac(atomj))%avdw(LJ_code)*iaclib(iac(atomk))%avdw(LJ_code)
                                                                        end do
                                                                end if
                                                        end if !topology atom or not
                                                else !use_fep
                                                        nb_list%AA(m) = iaclib(iac(atomj))%avdw(LJ_code)*iaclib(iac(atomk))%avdw(LJ_code)
                                                end if
                                        else !arithmetic
                                                if (use_fep) then
                                                        if (atomj.le.nat_solute) then
                                                                if (iqatom(atomj).ne.0) then !atom is a Q atom from the FEP file
                                                                        if (atomk.le.nat_solute) then
                                                                                if (iqatom(atomk).ne.0) then !other atom is also a Q atom
                                                                                        do tmpstates=1,nstates
                                                                                        nbq_list%AA(tmpstates,mq)=qavdw(qiac(iqatom(atomj),tmpstates),LJ_code)+qavdw(qiac(iqatom(atomk),tmpstates),LJ_code)
                                                                                        end do
                                                                                else !other atom is not a Q atom, have to use LJ_code = 1 for Q atom
                                                                                        do tmpstates=1,nstates
                                                                                        if (LJ_code .eq. 2) then
                                                                                                nbq_list%AA(tmpstates,mq)=qavdw(qiac(iqatom(atomj),tmpstates),1)+iaclib(iac(atomk))%avdw(LJ_code)
                                                                                        else
                                                                                                nbq_list%AA(tmpstates,mq)=qavdw(qiac(iqatom(atomj),tmpstates),LJ_code)+iaclib(iac(atomk))%avdw(LJ_code)
                                                                                        end if
                                                                                        end do
                                                                                end if
                                                                        else !other atom most likely water, not in original topology, so not a qatom
                                                                                do tmpstates=1,nstates
                                                                                if (LJ_code .eq. 2) then
                                                                                        nbq_list%AA(tmpstates,mq)=qavdw(qiac(iqatom(atomj),tmpstates),1)+iaclib(iac(atomk))%avdw(LJ_code)
                                                                                else
                                                                                        nbq_list%AA(tmpstates,mq)=qavdw(qiac(iqatom(atomj),tmpstates),LJ_code)+iaclib(iac(atomk))%avdw(LJ_code)
                                                                                end if
                                                                                end do
                                                                        end if
                                                                else 
                                                                        if  (atomk.le.nat_solute) then
                                                                                if (iqatom(atomk).ne.0) then !other atom is a Q atom
                                                                                        do tmpstates=1,nstates
                                                                                        if (LJ_code .eq. 2) then
                                                                                                nbq_list%AA(tmpstates,mq)=iaclib(iac(atomj))%avdw(LJ_code)+qavdw(qiac(iqatom(atomk),tmpstates),1)
                                                                                        else 
                                                                                                nbq_list%AA(tmpstates,mq)=iaclib(iac(atomj))%avdw(LJ_code)+qavdw(qiac(iqatom(atomk),tmpstates),LJ_code)
                                                                                        end if
                                                                                        end do
                                                                                else !not a qatom
                                                                                        do tmpstates=1,nstates
                                                                                                nbq_list%AA(tmpstates,mq) = iaclib(iac(atomj))%avdw(LJ_code)+iaclib(iac(atomk))%avdw(LJ_code)
                                                                                        end do
                                                                                end if
                                                                        else !other atom most likely water, not in original topology, so not a qatom
                                                                                do tmpstates=1,nstates
                                                                                        nbq_list%AA(tmpstates,mq) = iaclib(iac(atomj))%avdw(LJ_code)+iaclib(iac(atomk))%avdw(LJ_code)
                                                                                end do
                                                                        end if
                                                                end if !first atom qatom yes no 
                                                        else !first atom not in orig topology
                                                                if  (atomk.le.nat_solute) then
                                                                        if (iqatom(atomk).ne.0) then !other atom is also a Q atom
                                                                                do tmpstates=1,nstates
                                                                                nbq_list%AA(tmpstates,mq)=qavdw(qiac(iqatom(atomj),tmpstates),LJ_code)+qavdw(qiac(iqatom(atomk),tmpstates),LJ_code)
                                                                                end do
                                                                        else !other atom is not a Q atom, have to use LJ_code = 1 for Q atom
                                                                                do tmpstates=1,nstates
                                                                                        nbq_list%AA(tmpstates,mq) = iaclib(iac(atomj))%avdw(LJ_code)+iaclib(iac(atomk))%avdw(LJ_code)
                                                                                end do
                                                                        end if
                                                                else !other atom most likely water, not in original topology, so not a qatom
                                                                        do tmpstates=1,nstates
                                                                        nbq_list%AA(tmpstates,mq) = iaclib(iac(atomj))%avdw(LJ_code)+iaclib(iac(atomk))%avdw(LJ_code)
                                                                        end do
                                                                end if
                                                        end if !topology or not
                                                else ! use_fep
                                                        nb_list%AA(m) = iaclib(iac(atomj))%avdw(LJ_code)+iaclib(iac(atomk))%avdw(LJ_code)
                                                end if
                                        endif	!geometric / arithmetric
                                        if (use_fep) then
                                                if (atomj.le.nat_solute) then
                                                        if (iqatom(atomj).ne.0) then !atom is a Q atom from the FEP file
                                                                if (atomk.le.nat_solute) then
                                                                        if (iqatom(atomk).ne.0) then !other atom is also a Q atom
                                                                                do tmpstates=1,nstates
                                                                                nbq_list%BB(tmpstates,mq)=qbvdw(qiac(iqatom(atomj),tmpstates),LJ_code)*qbvdw(qiac(iqatom(atomk),tmpstates),LJ_code)
                                                                                end do
                                                                        else !other atom is not a Q atom, have to use LJ_code = 1 for Q atom
                                                                                do tmpstates=1,nstates
                                                                                if (LJ_code .eq. 2) then
                                                                                        nbq_list%BB(tmpstates,mq)=qbvdw(qiac(iqatom(atomj),tmpstates),1)*iaclib(iac(atomk))%bvdw(LJ_code)
                                                                                else
                                                                                        nbq_list%BB(tmpstates,mq)=qbvdw(qiac(iqatom(atomj),tmpstates),LJ_code)*iaclib(iac(atomk))%bvdw(LJ_code)
                                                                                end if
                                                                                end do
                                                                        end if
                                                                else !atom not in original topology
                                                                        do tmpstates=1,nstates
                                                                        if (LJ_code .eq. 2) then
                                                                                nbq_list%BB(tmpstates,mq)=qbvdw(qiac(iqatom(atomj),tmpstates),1)*iaclib(iac(atomk))%bvdw(LJ_code)
                                                                        else
                                                                                nbq_list%BB(tmpstates,mq)=qbvdw(qiac(iqatom(atomj),tmpstates),LJ_code)*iaclib(iac(atomk))%bvdw(LJ_code)
                                                                        end if
                                                                        end do
                                                                end if
                                                        else 
                                                                if (atomk.le.nat_solute) then!first atom not qatom
                                                                        if (iqatom(atomk).ne.0) then !other atom is a Q atom
                                                                                do tmpstates=1,nstates
                                                                                if (LJ_code .eq. 2) then
                                                                                        nbq_list%BB(tmpstates,mq)=iaclib(iac(atomj))%bvdw(LJ_code)*qbvdw(qiac(iqatom(atomk),tmpstates),1)
                                                                                else 
                                                                                        nbq_list%BB(tmpstates,mq)=iaclib(iac(atomj))%bvdw(LJ_code)*qbvdw(qiac(iqatom(atomk),tmpstates),LJ_code)
                                                                                end if
                                                                                end do
                                                                        else !not a qatom
                                                                                do tmpstates=1,nstates
                                                                                        nbq_list%BB(tmpstates,mq) = iaclib(iac(atomj))%bvdw(LJ_code)*iaclib(iac(atomk))%bvdw(LJ_code)
                                                                                end do
                                                                        end if
                                                                else !atom not in orig topology
                                                                        do tmpstates=1,nstates
                                                                        nbq_list%BB(tmpstates,mq) = iaclib(iac(atomj))%bvdw(LJ_code)*iaclib(iac(atomk))%bvdw(LJ_code)
                                                                        end do
                                                                end if
                                                        end if
                                                else !atom not in orig topology
                                                        if (atomk.le.nat_solute) then
                                                                if (iqatom(atomk).ne.0) then !other atom is a Q atom
                                                                        do tmpstates=1,nstates
                                                                        if (LJ_code .eq. 2) then
                                                                                nbq_list%BB(tmpstates,mq)=iaclib(iac(atomj))%bvdw(LJ_code)*qbvdw(qiac(iqatom(atomk),tmpstates),1)
                                                                        else 
                                                                                nbq_list%BB(tmpstates,mq)=iaclib(iac(atomj))%bvdw(LJ_code)*qbvdw(qiac(iqatom(atomk),tmpstates),LJ_code)
                                                                        end if
                                                                        end do
                                                                else !not qatom
                                                                        do tmpstates=1,nstates
                                                                                nbq_list%BB(tmpstates,mq) = iaclib(iac(atomj))%bvdw(LJ_code)*iaclib(iac(atomk))%bvdw(LJ_code)
                                                                        end do
                                                                end if
                                                        else  !atom not in orig topology
                                                                do tmpstates=1,nstates
                                                                nbq_list%BB(tmpstates,mq) = iaclib(iac(atomj))%bvdw(LJ_code)*iaclib(iac(atomk))%bvdw(LJ_code)
                                                                end do
                                                        end if
                                                end if !in topology or not
                                        else !not use_fep
                                                nb_list%BB(m) = iaclib(iac(atomj))%bvdw(LJ_code)*iaclib(iac(atomk))%bvdw(LJ_code)
                                        end if
                                        if (LJ_code .eq. 3) then
                                                fullscale=coulomb_constant*el14_scale
                                        else
                                                fullscale=coulomb_constant
                                        end if
                                        if(use_fep) then
                                                if (atomj.le.nat_solute) then
                                                        if (iqatom(atomj).ne.0) then !atom is a Q atom from the FEP file
                                                                if (atomk.le.nat_solute) then
                                                                        if (iqatom(atomk).ne.0) then !other atom is also a Q atom
                                                                                do tmpstates=1,nstates
                                                                                        nbq_list%qq(tmpstates,mq)=qcrg(iqatom(atomj),tmpstates)*qcrg(iqatom(atomk),tmpstates)*fullscale
                                                                                end do
                                                                        else !other atom is not a Q atom
                                                                                do tmpstates=1,nstates
                                                                                nbq_list%qq(tmpstates,mq)=qcrg(iqatom(atomj),tmpstates)*crg(atomk)*fullscale
                                                                                end do
                                                                        end if
                                                                else !atom not in orig topology
                                                                        do tmpstates=1,nstates
                                                                                nbq_list%qq(tmpstates,mq)=qcrg(iqatom(atomj),tmpstates)*crg(atomk)*fullscale
                                                                        end do
                                                                end if
                                                        else
                                                                if (atomk.le.nat_solute) then
                                                                        if (iqatom(atomk).ne.0) then !other atom is a Q atom
                                                                                do tmpstates=1,nstates
                                                                                nbq_list%qq(tmpstates,mq)=crg(atomj)*qcrg(iqatom(atomk),tmpstates)*fullscale
                                                                                end do
                                                                        else
                                                                                do tmpstates=1,nstates
                                                                                nbq_list%qq(tmpstates,mq) = crg(atomj)*crg(atomk)*fullscale !coulombconstant from topo.f90
                                                                                end do
                                                                        end if
                                                                else !not in topology
                                                                        do tmpstates=1,nstates
                                                                        nbq_list%qq(tmpstates,mq) = crg(atomj)*crg(atomk)*fullscale !coulombconstant from topo.f90
                                                                        end do
                                                                end if
                                                        end if
                                                else !first atom not in topology
                                                        if (atomk.le.nat_solute) then
                                                                if (iqatom(atomk).ne.0) then !other atom is a Q atom
                                                                        do tmpstates=1,nstates
                                                                        nbq_list%qq(tmpstates,mq)=crg(atomj)*qcrg(iqatom(atomk),tmpstates)*fullscale
                                                                        end do
                                                                else
                                                                        do tmpstates=1,nstates
                                                                        nbq_list%qq(tmpstates,mq) = crg(atomj)*crg(atomk)*fullscale !coulombconstant from topo.f90
                                                                        end do
                                                                end if
                                                        else !not in topology
                                                                do tmpstates=1,nstates
                                                                nbq_list%qq(tmpstates,mq) = crg(atomj)*crg(atomk)*fullscale !coulombconstant from topo.f90
                                                                end do
                                                        end if
                                                end if
                                        else !not use_fep
                                                nb_list%qq(m) = crg(atomj)*crg(atomk)*coulomb_constant*el14_scale
                                        end if !use_fep
                                        if (use_fep) then
                                                mq = mq+1  !increase Q index
                                        else
                                                m = m+1 !increase index
                                        end if
				else !NB
				endif

		end do
	end do
	if (use_fep) then
	nbq_list%number_of_NB = mq-1
	else
	nb_list%number_of_NB = m-1
	end if
end subroutine nb_make_list

!*******************************************************************

subroutine nb_heading(i)
	integer						::	i

	write(*,'(a6)', advance='no') 'Vlj'
	write(*,'(a8)', advance='no') 'Vel'
end subroutine nb_heading



!*******************************************************************
!* Add new calculation of nonbonded interactions between
!*  two sets of residues, e.g. protein aa and q-atoms. The
!*  interactions will be calculated for for the first group
!*  one at the time against all residues in the second group.
!*******************************************************************
integer function nb_qp_add(desc)
	!arguments
	character(*)				:: desc
	character(len=200)			:: line
	character*80				:: lambda !name of fep_file if loaded
	integer						:: ires, ats1, ats2, fstat
	integer					:: str_start, str_end, totlen

!	use_fep = .false.

	nb_qp_add = -1
	if (N_nb_qp .ge. MAX_NB_QP) then
		write (*,10) MAX_NB_QP
		return
	end if
10	format('Sorry, the maximum number of residue nonbond calculations is ',i2)
	Nlists = nres
	allocate (qp_aves(nres), nb_list_res(nres))
	allocate (nbq_list_res(nres))
!        allocate(iqatom(nat_solute))
!        iqatom(:)=0
	qp_aves(:)%lj = 0.
	qp_aves(:)%el = 0.
	call mask_initialize(masks(1))
	call mask_initialize(masks(2))
!Added error handling if the user wants to be smart
	write(*,*) ''
	qp_calc%p_first = get_int_arg("Enter first residue of protein:")
	write(*,*) ''
	write(*,*) 'First residue of protein:   ',qp_calc%p_first
	if (qp_calc%p_first .lt. 1) then
		write(*,*) 'First residue not valid as ', qp_calc%p_first
		write(*,*) 'First residue set to 1'
		qp_calc%p_first = 1
	elseif(qp_calc%p_first.gt.nres) then
		write(*,*) 'First residue is larger than number of residues in the topology'
		write(*,*) 'First residue set to ', nres
		qp_calc%p_first = nres
	endif
	qp_calc%p_last = get_int_arg("Enter last residue of protein:")
	write(*,*) ''
	write(*,*) 'Last residue of protein:   ',qp_calc%p_last
	if (qp_calc%p_last .gt. nres ) then
		write(*,*) 'Last residue is larger than number of residues in the topology'
		write(*,*) 'Last residue set to ',nres
		qp_calc%p_last = nres
	elseif(qp_calc%p_last.lt.qp_calc%p_first) then
		write(*,*) 'Last residue is smaller than first residue'
		write(*,*) 'Last residue set to first residue'
		qp_calc%p_last = qp_calc%p_first
	end if
	write(*,*) 'First residue ', qp_calc%p_first
	write(*,*) 'Last residue ', qp_calc%p_last
!changed Paul Bauer to allow FEP file loading
	if(.not.use_fep) then
	write(*,*) 'Enter the name of a FEP file if wanted, or terminate with a .'
	call getlin(fep_file,'')
	if ((fep_file.eq.'.').or.(fep_file.eq.'')) then
		write(*,*) 'No FEP file loaded'
		use_fep = .false.
	else
	        if(.not. prm_open(fep_file)) then
                write(*,'(a,a)') '>>>>> ERROR: Could not open fep file ',trim(fep_file)
                stop
	        end if
		use_fep = .true.
		call get_fep
	end if
!changed Paul Bauer to allow FEP file loading
	if (use_fep) then
!added line for lambda reading, will check with FEP file number of states
!to read the maximum number of states for this kind of FEP calculation
!code from md.f90 file initialisation
	write(*,*)'Give a value for the reaction coordinate (lambda values)'
	call getlin(lambda,'')
!        write(*,*) lambda	
        str_start = 1
        totlen = len_trim(lambda)
        do while(str_start < totlen)
                str_end = string_part(lambda, ' ', str_start)
		states=states+1
                lamda(states)=strtod(lambda(str_start:str_end))
                str_start = str_end + 2
        end do
!	write(*,*)(lamda(str_start),str_start=1,states)
	if (states.ne.nstates) then
		write(*,*)'Number of lambda values not equal to states in FEP file'
		write(*,'(i5,a,i5)') states,' lambda states not equal to FEP file states ', nstates
	end if
	end if
	end if
	write(*, *) 'Enter mask for q-atoms (e.g. <residue  xx yy>.)'
	write(*, *) 'Finalize with end.'
 	ats2 =  maskmanip_make(masks(2))

	if(ats2 == 0 ) then
		call mask_finalize(masks(2))
		return
	end if

	write(*, 13) masks(2)%included
13	format('Q-atom mask contains ',i6,' atoms')


	do ires = qp_calc%p_first, qp_calc%p_last
		call mask_clear(masks(1))
!Done like that to have the resdue numbers in the buffer
!because F90 is to stupid to have strings allocated
!with integers
		write(line, '(a,i6,i6)') 'residue ', ires, ires
		ats1 = mask_add(masks(1), line)
	if (.not.use_fep) then
		allocate(nb_list_res(ires)%atom1(ats1*ats2))
		allocate(nb_list_res(ires)%atom2(ats1*ats2))
		allocate(nb_list_res(ires)%AA(ats1*ats2))
		allocate(nb_list_res(ires)%BB(ats1*ats2))
		allocate(nb_list_res(ires)%qq(ats1*ats2))
                nb_list_res(ires)%atom1(:)=0
                nb_list_res(ires)%atom2(:)=0
                nb_list_res(ires)%AA(:)=0
                nb_list_res(ires)%BB(:)=0
                nb_list_res(ires)%qq(:)=0
        else
		allocate(nbq_list_res(ires)%atom1(ats1*ats2))
		allocate(nbq_list_res(ires)%atom2(ats1*ats2))
		allocate(nbq_list_res(ires)%AA(max_states,ats1*ats2))
		allocate(nbq_list_res(ires)%BB(max_states,ats1*ats2))
		allocate(nbq_list_res(ires)%qq(max_states,ats1*ats2))
                nbq_list_res(ires)%atom1(:)=0
                nbq_list_res(ires)%atom2(:)=0
                nbq_list_res(ires)%AA(:,:)=0
                nbq_list_res(ires)%BB(:,:)=0
                nbq_list_res(ires)%qq(:,:)=0
        end if
		call nb_make_list(masks(1), masks(2), nb_list_res(ires), nbq_list_res(ires))
	end do
	call mask_finalize(masks(1))
	call mask_finalize(masks(2))
	nb_qp_add = 1
end function nb_qp_add


!*******************************************************************
!* Step through one residue at the time and calculate nonbonded
!* interactions with all q-atoms.
!*******************************************************************
subroutine nb_qp_calc()
	!arguments
	integer							::  ires, count, i
	real(8) :: vdw,el


	total_qp_frames = total_qp_frames + 1
	
	do ires = qp_calc%p_first, qp_calc%p_last
		call nb_calc_lists(vdw, el, nb_list_res(ires),nbq_list_res(ires))
		qp_aves(ires)%lj = (qp_aves(ires)%lj*(total_qp_frames-1)+vdw)/(total_qp_frames)
		qp_aves(ires)%el = (qp_aves(ires)%el*(total_qp_frames-1)+el)/(total_qp_frames)
	end do

end subroutine nb_qp_calc

!*******************************************************************
!* Print final averages and clean up.
!*******************************************************************
subroutine nb_qp_finalize()
integer :: ires 
!Print header
	write(*,701) qp_calc%q_first, qp_calc%q_last
	write(*,702) 'Residue','Average LJ','Average el'
701	format('Average nonbonded interaction between residues', i4, ' to ',i4,' and other residues:')
702	format(a15,a15,a15)


!	Print averages of all trajectories
	do ires = qp_calc%p_first, qp_calc%p_last
		write (*,703) ires, qp_aves(ires)%lj, qp_aves(ires)%el
	end do
703	format(i15,f15.2,f15.2)
	deallocate(qp_aves, nb_list_res, nbq_list_res)  
end subroutine nb_qp_finalize

!*******************************************************************
end module CALC_NB
