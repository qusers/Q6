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
! mask.f90
! by John Marelius
! atom selection mask

module ATOM_MASK
	use TOPO

	implicit none

	! version data
	character(*), private, parameter  ::	MODULE_VERSION = '5.06'
	character(*), private, parameter  ::	MODULE_DATE    = '2014-04-21'

	!private procedures
	private update,update_pretop, set_solute, get_word, finalize_storage

	!declare variables
	logical(1), private, allocatable::	solute(:)

	type, private :: set
		integer				::	first, last
		logical				::	residue_numbering
		logical				::	solute1, solute2
		logical				::	excl1, excl2
		logical				::	restr1, restr2
		logical				::	heavy1, heavy2
		logical             ::  sybyl1, sybyl2
		logical				::  notflag
        character(5)             ::  sybylcode
	end type set

	type, public :: MASK_TYPE
		integer							::	included
		logical(1), pointer	::	mask(:)
	end type MASK_TYPE

contains


subroutine mask_startup
	call topo_startup
end subroutine mask_startup


subroutine mask_shutdown
	call finalize_storage
end subroutine mask_shutdown

subroutine mask_deallocate(m)
!arguments
	type(MASK_TYPE)				::	m
	if(associated(m%mask)) deallocate(m%mask)
end subroutine mask_deallocate

subroutine mask_initialize(m)
!arguments
	type(MASK_TYPE)				::	m

	allocate(m%mask(nat_pro))
	call mask_clear(m)
	if(.not. allocated(solute)) then
		allocate(solute(nat_pro))
		call set_solute
	end if
end subroutine mask_initialize

subroutine mask_finalize(m)
!arguments
	type(MASK_TYPE)				::	m

	if(associated(m%mask)) deallocate(m%mask)
	m%included = 0
end subroutine mask_finalize


subroutine mask_clear(m)
!arguments
	type(MASK_TYPE)				::	m

	!clear mask
	m%mask(:) = .false.
	m%included = 0
end subroutine mask_clear


subroutine mask_all(m)
!arguments
	type(MASK_TYPE)				::	m

	!clear mask
	m%mask(:) = .true.
	m%included = nat_pro
end subroutine mask_all


integer function mask_count(m)
!arguments
	type(MASK_TYPE)				::	m
	mask_count = m%included
end function mask_count

subroutine set_solute
	solute(1:nat_solute) = .true.
	solute(nat_solute+1:nat_pro) = .false.
end subroutine set_solute

!*****************************************************
!Add atoms to mask according to rules in 'line'.
!*****************************************************
integer function mask_add(m, line, pretop)
!arguments
    logical, intent(in),optional::  pretop
	type(MASK_TYPE)				::	m
	character(*)				::	line

	!locals
	character(len=80)			::	word
	integer						::	pos
	type(set)					::	s
	logical						::	empty
	integer						::	readstat

	!clear set description
	s%first = 1
	s%last = nat_pro !default to whole sequence
	s%residue_numbering = .false. !use atom numbering
	s%solute1 = .true.
	s%solute2 = .false.
	s%excl1 = .true.
	s%excl2 = .false.
	s%restr1 = .true.
	s%restr2 = .false.
	s%heavy1 = .true.
	s%heavy2 = .false.
	s%sybyl1 = .true.
	s%sybyl2 = .false.
	s%sybylcode = ""

	s%notflag = .false.
	empty = .true.
	mask_add = 0

	pos = 0
	do while(get_word(line, word, pos))
		call upcase(word)
		select case(word)
		case ('ALL')
			s%notflag = .false.
			empty = .false.
		case ('NOT', '-')
			s%notflag = .true.
		case ('RESIDUE', 'RES', 'RESI')
			s%residue_numbering = .true.
			s%notflag = .false.
		case ('SOLUTE')
			s%solute1 = .not. s%notflag
			s%solute2 = .not. s%notflag
			empty = .false.
			s%notflag = .false.
		case ('EXCLUDED', 'EXCL')
			s%excl1 = .not. s%notflag
			s%excl2 = .not. s%notflag
			empty = .false.
			s%notflag = .false.
		case ('RESTRAINED')
			s%restr1 = .not. s%notflag
			s%restr2 = .not. s%notflag
			empty = .false.
			s%notflag = .false.
		case ('HEAVY')
			s%heavy1 = .not. s%notflag
			s%heavy2 = .not. s%notflag
			empty = .false.
			s%notflag = .false.
		case ('SYBYL')
			s%sybyl1 = .not. s%notflag
			s%sybyl2 = .not. s%notflag
			s%notflag = .false.
			if(.not. get_word(line, word, pos)) then	! get (potential) SYBYL name
				write(*,800)
800				format('>>>>> ERROR: no SYBYL name found')
				exit
			end if
			s%sybylcode = word
			call get_sybylcode(s%sybylcode)
			empty = .false.
		case default
			read(word, *, iostat=readstat) s%first
			if(readstat /= 0) then
				write(*,900) trim(word)
900				format('>>>>> ERROR: unrecognised property ', a)
				empty = .true.
				exit
			elseif(s%residue_numbering .and. (s%first < 1 .or. &
				s%first > nres)) then
				write(*,911) s%first
911				format('>>>>> ERROR: invalid first residue ', i5)
				empty = .true.
				exit
			elseif(s%first < 1 .or. s%first > nat_pro) then
				write(*,910) s%first
910				format('>>>>> ERROR: invalid first atom ', i5)
				empty = .true.
				exit
			elseif(.not. get_word(line, word, pos)) then
				empty = .false. !default to all
				!failed to read last atom, assume last=first
				s%last = s%first
				exit
			else
				empty = .false. !default to all
				read(word, *, iostat=readstat) s%last
				if(readstat /= 0 .or. &
					(s%residue_numbering .and. s%last > nres) .or. &
					s%last < s%first .or. s%last > nat_pro) then
					if(s%residue_numbering) then
						write(*,931) trim(word)
					else
						write(*,930) trim(word)
					end if
930					format('>>>>> ERROR: invalid last atom ', a)
931					format('>>>>> ERROR: invalid last residue ', a)
					empty = .true.
					exit
				else
					exit
				end if
			end if
		end select
	end do
	if (.not. empty .and. present(pretop)) then
		if(pretop) then
		  mask_add=update_pretop(s, m)
		else
		  mask_add=update(s, m)
		end if
	elseif (.not. empty) then
	  mask_add=update(s, m)
	end if
end function mask_add

!*****************************************************

logical function get_word(line, word, pos)
!arguments
	character(*)				::	line
	character(*), intent(out)	::	word
	integer, intent(inout)		::	pos
!locals
	integer, save				::	linelen
	character, save				::	TAB = char(9)
	integer						::	start_field, end_field


	get_word = .false.

	if(pos <= 1) then
		pos = 1
		linelen = len_trim(line)
	end if
	if(linelen == 0) return
	if(pos > linelen) return

	do while(line(pos:pos) == ' ' .or. line(pos:pos) == TAB)
		pos = pos + 1
		if(pos > linelen) then
			return
		end if
	end do

	start_field = pos
	do while(line(pos:pos) /= ' ' .and. line(pos:pos) /= TAB)
		pos = pos + 1
		if(pos > linelen) then
			exit
		end if
	end do
	end_field = pos-1
	pos = pos +1 !update for next call
	word = line(start_field:end_field)
	get_word = .true.

end function get_word

!*****************************************************
! Update mask from rules in s.
! if s%property1 and s%property2 are true, then that
! requested without "not". if both are false then it has
! been requested with "not" flag. else it has not been requested
!*****************************************************
integer function update(s, m)
!arguments
	type(set)					::	s
	type(MASK_TYPE)		::	m
!locals
	integer						::	i
	character(5)			:: sybyl_upcase
	logical(1), allocatable::	sybyl(:)

	update = 0

	if(s%residue_numbering) then
		s%first = res(s%first)%start
		if(s%last == nres) then
			s%last = nat_pro
		else
			s%last = res(s%last+1)%start - 1
		endif
	end if


	allocate(sybyl(nat_pro))
	sybyl(:) = .false.

	if(s%sybyl1 .eqv. s%sybyl2) then  !make sybyl array
		write (*,*) "making sybyl array"
		do i = s%first,s%last
			sybyl_upcase = SYBYL_atom_type(iac(i))
			call upcase(sybyl_upcase)
			if(sybyl_upcase == s%sybylcode) then
				sybyl(i) = .true.
			end if
		end do
	end if



	do i = s%first, s%last
		if( (( solute(i) .eqv. s%solute1) .or. ( solute(i) .eqv. s%solute2)).and.&
		    (( excl(i) .eqv. s%excl1) .or. ( excl(i) .eqv. s%excl2)) .and. &
		    (( (shell(i) .or. excl(i)) .eqv. s%restr1) .or. (( shell(i) .or. excl(i)) .eqv. s%restr2)) .and. &
		    (( heavy(i) .eqv. s%heavy1) .or. ( heavy(i) .eqv. s%heavy2)) .and. &
		    (( sybyl(i) .eqv. s%sybyl1) .or. ( sybyl(i) .eqv. s%sybyl2))) then
			if(.not. m%mask(i)) then
				m%mask(i) = .true.
				update = update + 1
			end if
		end if
	end do
	m%included = m%included + update


	if(allocated(sybyl)) deallocate(sybyl)

end function update

!*****************************************************
! Update mask from rules in s.
! Used before topology is made, thus
! vectors 'excl' or 'shell' cannot be used.
!*****************************************************
integer function update_pretop(s, m)
!arguments
	type(set)					::	s
	type(MASK_TYPE)		::	m
!locals
	integer						::	i

	update_pretop = 0

	if(s%residue_numbering) then
		s%first = res(s%first)%start
		if(s%last == nres) then
			s%last = nat_pro
		else
			s%last = res(s%last+1)%start - 1
		endif
	end if

	do i = s%first, s%last
		if(((solute(i) .eqv. s%solute1) .or. (solute(i) .eqv. s%solute2)).and.&
		   ((heavy(i) .eqv. s%heavy1) .or. (heavy(i) .eqv. s%heavy2))) then
			if(.not. m%mask(i)) then
				m%mask(i) = .true.
				update_pretop = update_pretop + 1
			end if
		end if
	end do
	m%included = m%included + update_pretop

end function update_pretop


!*****************************************************
! Function: get_sybylcode()
! Argument: A (potential) SYBYL atom name
! Returns:  the integer atom code (iac) for that SYBYL name
!           (returns 0 if the argument is not a valid SYBYL name)
!*****************************************************
subroutine get_sybylcode(sybylname)
!arguments
	character(*)			::	sybylname		! (potential) SYBYL name

!locals
	integer					:: i = 1, check_sybylcode			! index
	character(len=5)		:: this_sybyl_caps	! SYBYL_atom_type(i) in CAPS

        check_sybylcode = 0
	do i = 1, natyps
		this_sybyl_caps = SYBYL_atom_type(i)
		call upcase(this_sybyl_caps)
		if(this_sybyl_caps == sybylname) then
			check_sybylcode = i
		end if
	end do

	if(check_sybylcode == 0) then
		write(*,940) sybylname
940				format('>>>>> ERROR: invalid SYBYL code ', a)
	end if

end subroutine get_sybylcode

!*****************************************************

subroutine finalize_storage
	if(allocated(solute)) deallocate(solute)
end subroutine finalize_storage


subroutine mask_get(m, x, xmasked)
	!arguments
	type(MASK_TYPE)			::	m
	TYPE(qr_vec)			::	x(:)
	TYPE(qr_vec)			::	xmasked(:)
	integer				::	i,j

	if(size(x) < nat_pro .or. size(xmasked) < m%included) then
		write(*,900)
900		format('>>>>> ERROR: invalid coordinate array size.')
		return
	end if

	j=0
	do i=1, nat_pro
		if(m%mask(i)) then
			j=j+1
			xmasked(j) = x(i)
		end if
	end do
end subroutine mask_get

!***************************************************
!transfer coords from masked atoms in 'xmasked'
! to topology coord array 'x'
!***************************************************
subroutine mask_put(m, x, xmasked)
	!arguments
	type(MASK_TYPE)			::	m
	TYPE(qr_vec)			::	x(:)
	TYPE(qr_vec)			::	xmasked(:)
	integer				::	i,j

	if(size(x) < nat_pro .or. size(xmasked) < m%included) then
		write(*,900)
900		format('>>>>> ERROR: invalid coordinate array size.')
		return
	end if

	j=0
	do i=1, nat_pro
		if(m%mask(i)) then
			j=j+1
			x(i) = xmasked(j)
		end if
	end do
end subroutine mask_put

!*****************************************************

end module ATOM_MASK

