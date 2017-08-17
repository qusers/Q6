! Q6: A comprehensive simulation package for molecular dynamics simulations and 
! free energy calculations, including empirical valence bond simulations, 
! linear interaction energy calculations, and free energy perturbation.
! 
! Copyright © 2017 Johan Åqvist, John Marelius, Shina Caroline Lynn Kamerlin and Paul Bauer
! 
! This program is free software; you can redistribute it and/or modify it under the 
! terms of the GNU General Public License as published by the Free 
! Software Foundation; either version 2 of the License, or any later version.
! 
! This program is distributed in the hope that it will be useful, 
! but WITHOUT ANY WARRANTY; without even the implied warranty of 
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
! See the GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License along with 
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin 
! Street, Fifth Floor, Boston, MA  02110-1301, USA. Also add information on 
! how to contact you by electronic and paper mail.
! calc_base.f90
! by John Marelius
! shared data for trajectory calculation modules

module CALC_BASE
	use TOPO
	use QATOM
	use PRMFILE
	implicit none
	
	real(kind=prec), allocatable :: xin(:)       ! All functions share this vector which contains a frame (J)
        real(kind=prec), allocatable :: calc_xtop(:) ! Needed because of new data structures
        TYPE(qr_vec),allocatable     :: tmp_coords(:),xin2(:)
   ! real(8), allocatable    ::  xin2d(:,:)       ! Schlitters formula needs all coordinates at the same time
	!for FEP file loading
	real(kind=prec)                         :: lamda(max_states)
	integer                                 :: alloc_stat_gen
	integer                                 :: states
	logical                                 :: use_fep
        real(kind=prec)                         :: totlambda=zero
        real(kind=prec),parameter               :: eps=QREAL_EPS

   ! place FEP file loading routine here so all processes can use it

contains

subroutine get_fep
! local variables
character					::	libtext*80,qaname*2
integer					::	i,j,k,iat
!temp. array for reallocating long-range exclusion list
integer(AI), pointer	::	tempexlong(:,:)
	allocate(iqatom(nat_pro))
        iqatom(:)=0
! --- # states, # q-atoms
if(.not. qatom_load_atoms(fep_file)) then
        call die_general('failure to load Q-atoms from FEP file.')
end if

! set flags
do i=1,nqat
        if(iqseq(i) > 0 .and. iqseq(i) <= nat_solute)  then
                iqatom(iqseq(i)) = i
        else if(iqseq(i) == 0) then
                write(*,10) i
        else
                write(*,20) i, iqseq(i)
                call die_general('invalid q-atom data')
        end if
end do
10	format('>>> WARNING: Q-atom no. ',i2,' is not associated with a topology atom.')
20	format('>>>>> ERROR: Q-atom no. ',i2,' has invalid topology number ',i5)
!allocate memory for qatom charges
allocate(qcrg(nqat,nstates), stat=alloc_stat_gen)
call check_alloc_general(alloc_stat_gen,'Qatom charges')

! --- copy topology charges

do i=1,nqat
        do j=1,nstates
                qcrg(i,j)=crg(iqseq(i))
        end do
end do

!initialize softcore lookup array
allocate (sc_lookup(nqat,natyps+nqat,nstates))
sc_lookup(:,:,:)=0.0

!load rest of fep file
if(.not. qatom_load_fep(fep_file)) then
        call die_general('failure to load FEP file.')
end if

!Adapt LJ parameters to topology
!If arithmetic combination rule take sqrt(epsilon) now
if (qvdw_flag .and. ivdw_rule .eq. 2 ) then
        qbvdw(:,1) = sqrt( qbvdw(:,1) )
        qbvdw(:,3) = sqrt( qbvdw(:,3) )
end if

!remove redefined bonded interactions from topology
if(nqbond > 0 .or. nqangle > 0 .or. nqtor > 0 .or. nqimp > 0 ) then
        write(*,*)
        call centered_heading('Removing redefined interactions from topology','-')
230		format('type',t10,' atom1 atom2 atom3 atom4')
        write(*,230)
231		format(a,t10,4i6)
        !remove bonds that were redefined
        do i=1,nbonds
                do j=1,nqbond
                        if ( (bnd(i)%i==qbnd(j)%i .and. bnd(i)%j==qbnd(j)%j) .or. &
                                (bnd(i)%i==qbnd(j)%j .and. bnd(i)%j==qbnd(j)%i) ) then
                                bnd(i)%cod = 0
                                write (*,231) 'bond',bnd(i)%i,bnd(i)%j
                        end if
                end do
        end do

        !remove angles that were redefined
        do i=1,nangles
                do j=1,nqangle
                        if((ang(i)%i.eq.qang(j)%i .and. ang(i)%j.eq.qang(j)%j .and. &
                                ang(i)%k.eq.qang(j)%k)                          .or. &
                                (ang(i)%i.eq.qang(j)%k .and. ang(i)%j.eq.qang(j)%j .and. &
                                ang(i)%k.eq.qang(j)%i) )                         then

                                ang(i)%cod = 0
                                write (*,231) 'angle',ang(i)%i,ang(i)%j,ang(i)%k
                        end if
                end do
        end do

        !remove torsions that were redefined
        do i=1,ntors
                do j=1,nqtor
			if(( (tor(i)%i.eq.qtor(j)%i .and. tor(i)%j.eq.qtor(j)%j .and. &
				tor(i)%k.eq.qtor(j)%k .and. tor(i)%l.eq.qtor(j)%l) .or. &
				(tor(i)%i.eq.qtor(j)%l .and. tor(i)%j.eq.qtor(j)%k .and. &
				tor(i)%k.eq.qtor(j)%j .and. tor(i)%l.eq.qtor(j)%i )) .and. &
				tor(i)%cod .ne. 0) then
                                tor(i)%cod = 0
                                write (*,231) 'torsion', tor(i)%i,tor(i)%j,tor(i)%k,tor(i)%l
						end if
                end do
        end do


        !remove impropers that were redefined
        select case(ff_type)
        case(FF_CHARMM) !special code for CHARMM
                do i=1,nimps
                        do j=1,nqimp
				if ( ((	(imp(i)%i.eq.qimp(j)%i) .or. &
					(imp(i)%i.eq.qimp(j)%l) .or. &
					(imp(i)%l.eq.qimp(j)%i)) .and. &
				     (	(imp(i)%j.eq.qimp(j)%l) .or. &
					(imp(i)%j.eq.qimp(j)%i) .or. &
					(imp(i)%j.eq.qimp(j)%j) .or. &
					(imp(i)%j.eq.qimp(j)%k)) .and. &
				     (	(imp(i)%k.eq.qimp(j)%i) .or. &
					(imp(i)%k.eq.qimp(j)%j) .or. &
					(imp(i)%k.eq.qimp(j)%k) .or. &
					(imp(i)%k.eq.qimp(j)%l))) .and. &
				     (	imp(i)%cod .ne. 0)) then
                                        imp(i)%cod = 0
                                        write (*,231) &
                                        'improper',imp(i)%i,imp(i)%j,imp(i)%k,imp(i)%l
                                end if
                        end do
                end do

        case default
        do i=1,nimps
                do j=1,nqimp
                        if(( (imp(i)%i.eq.qimp(j)%i) .and. (imp(i)%j.eq.qimp(j)%j) .and. &
                                (imp(i)%k.eq.qimp(j)%k) .and. (imp(i)%l.eq.qimp(j)%l)) .or. &
                                ((imp(i)%i.eq.qimp(j)%k) .and. (imp(i)%j.eq.qimp(j)%j) .and. &
                                (imp(i)%k.eq.qimp(j)%i) .and. (imp(i)%l.eq.qimp(j)%l)) .or. &
                                ((imp(i)%i.eq.qimp(j)%k) .and. (imp(i)%j.eq.qimp(j)%j) .and. &
                                (imp(i)%k.eq.qimp(j)%l) .and. (imp(i)%l.eq.qimp(j)%i)) .or. &
                                ((imp(i)%i.eq.qimp(j)%l) .and. (imp(i)%j.eq.qimp(j)%j) .and. &
                                (imp(i)%k.eq.qimp(j)%i) .and. (imp(i)%l.eq.qimp(j)%k)) .or. &
                                ((imp(i)%i.eq.qimp(j)%l) .and. (imp(i)%j.eq.qimp(j)%j) .and. &
                                (imp(i)%k.eq.qimp(j)%k) .and. (imp(i)%l.eq.qimp(j)%i))) then
                                imp(i)%cod = 0
                                write(*,231)'improper',imp(i)%i,imp(i)%j,imp(i)%k,imp(i)%l
                        end if
                end do
        end do
        end select
end if

!check special exclusions
!modify exclusion lists to inclue special exclusions between Q and non-Q
if(nexspec > 0) then
        allocate(tempexlong(2,nexlong+nexspec))
        tempexlong(:, 1:nexlong) = listexlong(:, 1:nexlong)
        deallocate(listexlong)
        listexlong => tempexlong
end if

do k = 1, nexspec
        i = exspec(k)%i
        j = exspec(k)%j
        if(i < 1 .or. i > nat_pro .or. j < 1 .or. j > nat_pro) then
                write(*, 592) k, i, j
                call die_general('invalid special exclusion data')
        end if
        !if one or more non-Q-atoms modify exclusion lists
        if(iqatom(i)==0 .or. iqatom(j)==0) then
                !With non-Q-atoms involved only accept all or no states
                if(any(exspec(k)%flag(1:nstates))) then
                        if(.not. all(exspec(k)%flag(1:nstates))) then
                                write(*,594) k
                                call die_general('invalid special exclusion data')
                        else !exlcude in all states
                                if(iabs(j-i) <= max_nbr_range) then
                                        if(i < j) then
                                                listex(j-i,i) = .true.
                                        else
                                                listex(i-j,j) = .true.
                                        end if
                                else
                                        nexlong = nexlong + 1
                                        listexlong(1, nexlong) = i
                                        listexlong(2, nexlong) = j
                                end if
                        end if
                end if
        end if
end do
592	format('>>>>> ERROR: Special exclusion pair ',i2,' (',i5,1x,i5,') is invalid')
594	format('>>>>> ERROR: Non-Q-atom special excl. pair ',i2,' must be on in all or no states')
end subroutine get_fep

! need this to keep qcalc working with new data structures, before someone else cleans
! up the whole mess that is qcalc
! returnes a qr_vec object of the correct size for the current operation, to be used in subroutines
! that are no longer retarded
function translate_coords(incoords)
! arguments
real(kind=prec),INTENT(IN) :: incoords(:)
! locals and returns
TYPE(qr_vec) :: translate_coords(SIZE(incoords)/3)
integer      :: length
length = SIZE(incoords)/3
translate_coords(1:length)%x = incoords(1:(length*3)-2:3)
translate_coords(1:length)%y = incoords(2:(length*3)-1:3)
translate_coords(1:length)%z = incoords(3:(length*3)  :3)
end function translate_coords

function btranslate_coords(incoords)
! arguments
TYPE(qr_vec),INTENT(IN) :: incoords(:)
! locals and returns
real(kind=prec) :: btranslate_coords(SIZE(incoords)*3)
integer      :: length
length = SIZE(incoords)
btranslate_coords(1:(length*3)-2:3) = incoords(1:length)%x 
btranslate_coords(2:(length*3)-1:3) = incoords(1:length)%y 
btranslate_coords(3:(length*3)  :3) = incoords(1:length)%z 

end function btranslate_coords

end module CALC_BASE
