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
! parse.f90
! by John Marelius
! command parser 

module PARSE
	use MISC
	implicit none

!constants
	character(*), private, parameter	::	MODULE_VERSION = '5.06'
	character(*), private, parameter	::	MODULE_DATE    = '2014-04-21'

	integer, parameter			::	MAX_ARGS = 10

	type substring
		integer					::	istart, iend
	end type substring

	character(200)				::	inbuf
	type(substring)				::	argv(0:MAX_ARGS)
	integer						::	argc = 0, argp = 1
	
	logical,private						:: read_from_file = .false.  !Enables input from file
	integer,private						:: INFILE = 0


contains

subroutine parse_startup()
	argc = 0
	argp = 1
end subroutine parse_startup

!--------------------------------------------------------------------------

logical function parse_open_file(filename)
character(*)	:: filename
integer			:: stat


	parse_open_file = .false.
	read_from_file = .true.

	INFILE = freefile()
	open(unit=INFILE, file = filename, status = 'old', form='formatted', &
		action='read', iostat = stat, access='sequential')
	if (stat /= 0) then
		parse_open_file = .false.
		INFILE = 0
		return
	endif

	parse_open_file = .true.
end function parse_open_file


!--------------------------------------------------------------------------
!Moved from Misc
!--------------------------------------------------------------------------

integer function openit(lun, fil, stat, frm, mode)
!arguments
	INTEGER lun
	CHARACTER( * ) fil, stat, frm
	character(*)				::	mode
!locals
	integer						::	errcode
	
	if(trim(fil).eq.'') write(*,*) 'Invalid filename'

	errcode = -1
!.......................................................................
	do while(errcode /= 0)
		OPEN(unit = lun, file = fil, status = stat, form = frm, &
			iostat = errcode, action=mode)
		if(errcode /= 0) then
			WRITE( * , 10) trim(fil), lun, errcode
			CALL get_string_arg(fil, ' >>>>> Enter a new file name (or nothing to cancel): ')
			if(fil == '' .or. fil == 'nothing') exit
		end if
	end do
	openit = errcode

10	format('>>>>> Failed to open ',a,' (unit ',i2,' error code ',i4,')')
!.......................................................................
end function openit

!--------------------------------------------------------------------------

!split text in inbuf into 
subroutine split
	integer						::	p
	logical						::	ws_flag 
	logical						::	quote_flag 
	character, save				::	TAB = achar(9)
	integer						::	inlen
	integer						::	trimlen(1)

	argc = 0
	argp = 1
	ws_flag = .true.
	quote_flag = .false.
	trimlen = len_trim(inbuf)
	inlen = trimlen(1)
	do p = 1, inlen
		if(inbuf(p:p) == '"') then
			if(quote_flag) then !closing quotation mark
				quote_flag = .false.
				ws_flag = .true.
			else !opening quotation mark
				quote_flag = .true. !set to accept all characters (like ' ') in string
				ws_flag = .true. !set to begin string at next positon
			end if
		elseif((inbuf(p:p) == ' ' .or. inbuf(p:p) == TAB) .and. .not. quote_flag) then
			if(.not. ws_flag) argv(argc)%iend = p-1
			ws_flag = .true.
		elseif(ws_flag) then
			argc = argc + 1
			if(argc > MAX_ARGS) then
				write(*,*) 'WARNING: Command argument overflow'
				exit
			end if
			argv(argc)%istart = p
			argv(argc)%iend = p !step up character by character
			ws_flag = .false.
		else
			argv(argc)%iend = p !step up character by character
			ws_flag = .false.
		end if
	end do
end subroutine split

!--------------------------------------------------------------------------

subroutine getline()
	do
		if (read_from_file) then
			read(INFILE, '(a200)') inbuf
		else  
			read(*, '(a200)') inbuf
		endif
		inbuf = adjustl(inbuf)
		!only exit if not a comment line
		if(scan(inbuf(1:1), '!#*') == 0) exit
	end do
	call split
end subroutine getline

!--------------------------------------------------------------------------

subroutine get_string_arg(arg, prompt)
!arguments
	character(*), intent(out)	::	arg
	character(*), optional, intent(in)	::	prompt
	
	do while(argp > argc)
		if(present(prompt))	write(unit=*, fmt='(a)', advance='no') prompt
		call getline
	end do
	arg = inbuf(argv(argp)%istart :argv(argp)%iend)
	argp = argp + 1
		
end subroutine get_string_arg

!--------------------------------------------------------------------------

subroutine get_line_arg(arg, prompt)
!arguments
	character(*), intent(out)	::	arg
	character(*), optional, intent(in)	::	prompt
	
	do while(argp > argc)
		if(present(prompt))	write(unit=*, fmt='(a)', advance='no') prompt
		call getline
	end do
	!return the rest of the line	
	arg = inbuf(argv(argp)%istart : len_trim(inbuf))
	call parse_reset !reset argument counter - new line next time		
end subroutine get_line_arg

!--------------------------------------------------------------------------

logical function get_string_single_line(arg, prompt)
!arguments
	character(*), intent(out)	::	arg
	character(*), optional, intent(in)	::	prompt
	
	if(argc == 0) then
		if(present(prompt))	write(unit=*, fmt='(a)', advance='no') prompt
		call getline
	end if
	if(argp>argc) then
		arg = ''
		get_string_single_line = .false.
	else
		arg = inbuf(argv(argp)%istart :argv(argp)%iend)
		argp = argp + 1
		get_string_single_line = .true.
	end if		
		
end function get_string_single_line

!--------------------------------------------------------------------------

subroutine parse_reset
	argc = 0
	argp = 1
end subroutine parse_reset

!--------------------------------------------------------------------------

integer function get_int_arg(prompt)
!arguments
	character(*), optional, intent(in)	::	prompt
!locals
	integer						::	value

1	do while(argp > argc)
		if(present(prompt))	write(unit=*, fmt='(a)', advance='no') prompt
		call getline
	end do
	
	read(inbuf(argv(argp)%istart :argv(argp)%iend), fmt=*, err=100) value

	get_int_arg = value

	argp = argp + 1
	return

100	write(*,*) 'Please enter an integer value!'
	call parse_reset
	goto 1
end function get_int_arg

!--------------------------------------------------------------------------

real(kind=prec) function get_real_arg(prompt)
!arguments
	character(*), optional, intent(in)	::	prompt
!locals
	real(kind=prec)						::	value

1	do while(argp > argc)
		if(present(prompt))	write(unit=*, fmt='(a)', advance='no') prompt
		call getline
	end do
	
	read(inbuf(argv(argp)%istart :argv(argp)%iend), fmt=*, err=100) value
	get_real_arg = value

	argp = argp + 1
	return

100	write(*,*) 'Please enter a number!'
	call parse_reset
	goto 1
end function get_real_arg

!--------------------------------------------------------------------------

end module PARSE
