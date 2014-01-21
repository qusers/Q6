!	(C) 2000 Uppsala Molekylmekaniska HB, Uppsala, Sweden

!	prmfile.f90
!	by John Marelius

!	parsing of data files with sections+keywords (input/parameter/library/FEP)

module PRMFILE

!This module reads files with the following format
!	[section_name]
!	key  value
!	key				value  !comment
!	#comment
!	!comment
!   *comment

!	key value
!	[another_section]
!	value value value value
!	value  value value		value
!	value value		value value

!section_name identifies a section in the file which can be opened with 
!prm_open_section. The name is not case-sensitive.
!Blank lines and lines starting with '!', '*' or '#' are ignored. 
!White space anywhere in key-value lines is ignored.
!Values may be followed by comments initated by '!', '*' or '#'
!Keys identify the values by a label. A key name may not contain white-space.


!The parameter file is opened by the function
!prm_open(file) where "file" is the file name. prm_open returns .true. on 
!success or .false. on failure (if the file does not exist or can't be opened)

!Data from the file is retrieved by calling the data retrieval functions 
!which all return .true. on success or .false. on failure
!Entire lines within a section are retrieved by prm_get_line(line) where 
!line is a string of suitable length to hold the expected data

!Key-value pairs can be retrieved by the following functions
!prm_get_value_by_key(key, value)
!	with the input argument "key" identifying the key and the
!	output argument "value" which is an integer, a real or a string
!	(different functions with same name for different data types)

!prm_get_string_string(key, value) 
!	key and value are both strings
!prm_get_string_int(key, value) 
!	key is string and value is integer
!prm_get_string_real(key, value) 
!	key is string and value is real
!prm_get_int_real(key, value) 
!	key is integer and value is real
!prm_get_int_int(key, value) 
!	key and value are both integers
!prm_get_int(value)
!	there is no key, only an integer value


use MISC
implicit none
!constants
	character(*), private, parameter	:: MODULE_VERSION = '5.01'
	character(*), private, parameter	:: MODULE_DATE = '2003-06-02'

!types, NOTE: maximum length of 500 chars per line
type LINE_TYPE
	character(len=900)			::	text
	type(LINE_TYPE), pointer	::	next
end type LINE_TYPE

type SECTION_TYPE
	character(len=80)			::	name
	type(SECTION_TYPE), pointer 	::	next
	type(LINE_TYPE), pointer	::	lp
	integer				::	count, max_enum
end type SECTION_TYPE

!module variables
type(SECTION_TYPE), pointer		::	first_sec, current_sec
type(LINE_TYPE), pointer		::	current_lp
integer					::	section_count

character(len=80)				::	current_title, next_title
integer, private				:: PRM_U = 0
character(200), private			:: PRM_FN
integer, private				:: stat

private							:: find_section
private							:: split
private							:: splitone
private							:: get_strings
private							::	load
contains

subroutine prmfile_startup

end subroutine prmfile_startup

!---------------------------------------------------------------------

logical function prm_get_integer_by_key(key, value, default)
!arguments
	character(*)				::	key
	integer, intent(out)		::	value
	integer, intent(in), optional:: default

!locals
	character(len=80)			::	inkey

	do while (prm_get_string_int(inkey, value))
		if(inkey .eq. key) then
			prm_get_integer_by_key = .true.
			return	
		end if
	end do	
	!go to beginning of section and try again
	call rewind_section
	do while (prm_get_string_int(inkey, value))
		if(inkey .eq. key) then
			prm_get_integer_by_key = .true.
			return	
		end if
	end do	

	if(present(default)) then
		value = default
		prm_get_integer_by_key = .true.
	else
		value = 0
		prm_get_integer_by_key = .false.
	end if
end function prm_get_integer_by_key

!--------------------------------------------------------------

logical function prm_get_real_by_key(key, value, default)
!arguments
	character(*)				::	key
	real, intent(out)			::	value
	real, intent(in), optional	:: default
!locals
	character(len=80)			::	inkey

	do while (prm_get_string_real(inkey, value))
		if(inkey .eq. key) then
			prm_get_real_by_key = .true.
			return	
		end if
	end do	
	!go to beginning of section and try again
	call rewind_section
	do while (prm_get_string_real(inkey, value))
		if(inkey .eq. key) then
			prm_get_real_by_key = .true.
			return	
		end if
	end do	

	if(present(default)) then
		value = default
		prm_get_real_by_key = .true.
	else
		value = 0.
		prm_get_real_by_key = .false.
	end if

end function prm_get_real_by_key

!-----------------------------------------------------------------

logical function prm_get_real8_by_key(key, value, default)
!arguments
	character(*)				::	key
	real(8), intent(out)			::	value
	real(8), intent(in), optional	:: default
!locals
	character(len=80)			::	inkey

	do while (prm_get_string_real8(inkey, value))
		if(inkey .eq. key) then
			prm_get_real8_by_key = .true.
			return	
		end if
	end do	
	!go to beginning of section and try again
	call rewind_section
	do while (prm_get_string_real8(inkey, value))
		if(inkey .eq. key) then
			prm_get_real8_by_key = .true.
			return	
		end if
	end do	

	if(present(default)) then
		value = default
		prm_get_real8_by_key = .true.
	else
		value = 0.
		prm_get_real8_by_key = .false.
	end if

end function prm_get_real8_by_key

!----------------------------------------------------------

logical function prm_get_string_by_key(key, value, default)
!arguments
	character(*)				::	key
	character(*), intent(out)	::	value
	character(*), intent(in), optional	::	default
!locals
	character(len=80)			::	mykey,inkey, invalue

	mykey = trim(key)
	call upcase(mykey)

	do while (prm_get_string_string(inkey, invalue))
		call upcase(inkey)
		if(trim(inkey) .eq. mykey) then
			prm_get_string_by_key = .true.
			value = invalue
			return	
		end if
	end do	
	!go to beginning of section and try again
	call rewind_section
	do while (prm_get_string_string(inkey, invalue))
		call upcase(inkey)
		if(trim(inkey) .eq. mykey) then
			prm_get_string_by_key = .true.
			value = invalue
			return	
		end if
	end do	

	if(present(default)) then
		value = default
		prm_get_string_by_key = .true.
	else
		!don't change value if not found
		prm_get_string_by_key = .false.
	end if

end function prm_get_string_by_key

!----------------------------------------------------------------

logical function prm_get_line_by_key(key, value, default)
!arguments
	character(*)				::	key
	character(*), intent(out)	::	value
	character(*), intent(in), optional	::	default
!locals
	character(len=80)			::	inkey
	character(len=200)			::	invalue

	do while (prm_get_string_line(inkey, invalue))
		if(inkey .eq. key) then
			prm_get_line_by_key = .true.
			value = invalue
			return	
		end if
	end do	
	!go to beginning of section and try again
	call rewind_section
	do while (prm_get_string_line(inkey, invalue))
		if(inkey .eq. key) then
			prm_get_line_by_key = .true.
			value = invalue
			return	
		end if
	end do	

	if(present(default)) then
		value = default
		prm_get_line_by_key = .true.
	else
		!don't change value if not found
		prm_get_line_by_key = .false.
	end if

end function prm_get_line_by_key

!------------------------------------------------------------------

logical function prm_get_logical_by_key(key, value, default)
!arguments
	character(*)				::	key
	logical, intent(out)		::	value
	logical, intent(in), optional::	default

!locals
	character(len=80)			::	inkey
	character(len=80)			::	invalue
	
	do while (prm_get_string_string(inkey, invalue))
		if(inkey .eq. key) then
			goto 100
		end if
	end do	
	!go to beginning of section and try again
	call rewind_section
	do while (prm_get_string_string(inkey, invalue))
		if(inkey .eq. key) then
			goto 100
		end if
	end do	

	if(present(default)) then
		value = default
		prm_get_logical_by_key = .true.
	else
		!don't change value if not found
		prm_get_logical_by_key = .false.
	end if
	return

100	if(invalue == 'on' .or. invalue == 'ON' .or. invalue =='1') then
		value = .true.
		prm_get_logical_by_key = .true.
	elseif(invalue == 'off' .or. invalue == 'OFF' .or. invalue =='0') then
		value = .false.
		prm_get_logical_by_key = .true.
	else
		prm_get_logical_by_key = .false.
	end if
end function prm_get_logical_by_key

!--------------------------------------------------------------------------------

logical function prm_get_line(line)
!arguments
	character(*), intent(out)	::	line
	
	!if there's a line
	if(associated(current_lp)) then
		!retrieve it	
		line = current_lp%text
		!and move one step forward
		current_lp => current_lp%next
		prm_get_line = .true.
	else
		line = ''
		prm_get_line = .false.
	end if

end function prm_get_line

!--------------------------------------------------------------------------

logical function prm_get_field(field, skip)
!arguments
	character(*), intent(out)	::	field
	logical, optional, intent(in)::	skip
!locals
	character(len=500), save	::	line
	integer, save				::	pos, linelen
	type(LINE_TYPE), pointer, save :: my_line
	character, save				::	TAB = char(9)
	integer						::	start_field, end_field
	logical						::	prm_res

	prm_get_field = .false.

	if(present(skip)) then
		if(skip) then
			nullify(my_line)
			pos = 0
			linelen = 0
			return
		end if
	end if
			
	!if other lines were read since last call, get a new one
	if(associated(current_lp) .and. .not. associated(my_line, current_lp)) then
		if(prm_get_line(line)) then
			my_line => current_lp
			pos = 1
			linelen = len_trim(line)
		else
			nullify(my_line)
			pos = 0
			linelen = 0
			return
		end if
	end if

	do while(line(pos:pos) == ' ' .or. line(pos:pos) == TAB)
		pos = pos + 1
		if(pos > linelen) then
			nullify(my_line) !next time new line
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
	field = line(start_field:end_field)
	prm_get_field = .true.
		
end function prm_get_field

!---------------------------------------------------------------------------------
! Finds and return the number of input lines in the section 'section'
!---------------------------------------------------------------------------------
integer function prm_count(section)
!argument
	character(*),intent(in)	::	section

	if(find_section(section)) then
		prm_count = current_sec%count
	else
		prm_count = 0
	end if

end function prm_count

!-------------------------------------------------------------------------------
! Find and return largest number in first column in specified section.
! Optional: Return number of lines in section in 'count_out' 
!-------------------------------------------------------------------------------
integer function prm_max_enum(section, count_out)
!argument
	character(*),intent(in)	::	section
	integer, optional, intent(out):: count_out
!locals
	character(200)				::	line
	integer						::	max_enum, enum, count
	logical						::	dummy

	max_enum = 0
	count = 0
	if(find_section(section)) then
		do while(prm_get_int(enum))
			if(enum > max_enum) max_enum = enum
		end do
		!go back
		call rewind_section
		count = current_sec%count
	end if
	prm_max_enum = max_enum
	if(present(count_out)) count_out = count
end function prm_max_enum

!-------------------------------------------------------------------------------
! Find and return largest number in second column in specified section.
! Optional: Return number of lines in section in 'count_out' 
!-------------------------------------------------------------------------------
integer function prm_max_enum2(section, count_out)
!argument
	character(*),intent(in)	::	section
	integer, optional, intent(out):: count_out
!locals
	character(200)				::	line, tmp_str
	integer						::	max_enum, enum, count
	logical						::	dummy

	max_enum = 0
	count = 0
	if(find_section(section)) then
		do while(prm_get_string_int(tmp_str, enum))
			if(enum > max_enum) max_enum = enum
		end do
		!go back
		call rewind_section
		count = current_sec%count
	end if
	prm_max_enum2 = max_enum
	if(present(count_out)) count_out = count
end function prm_max_enum2

!-----------------------------------------------------------------------------

logical function prm_get_string_string(key, value)
!arguments
	character(*),intent(out)	::	key
	character(*), intent(out)	::	value

	prm_get_string_string = get_strings(key, value)
end function prm_get_string_string

!-----------------------------------------------------------------------------

logical function prm_get_string_line(key, value)
!arguments
	character(*),intent(out)	::	key
	character(*), intent(out)	::	value
!locals
	character(200)				::	line
	logical						::	stat

	key = ''
	value = ''
	prm_get_string_line = .false.

	stat = prm_get_line(line)
	if(.not. stat) return
	if(.not. splitone(line, key, value)) return

	prm_get_string_line = .true.
end function prm_get_string_line

!-----------------------------------------------------------

logical function prm_get_string_int(key, value)
!arguments
	character(*),intent(out)	::	key
	integer, intent(out)		::	value
!locals
	character(200)				::	str_value


	prm_get_string_int = get_strings(key, str_value)
	read(str_value, fmt=*, iostat=stat) value !read integer from value

	if (stat == 0 ) then
		prm_get_string_int = .true.
	else 
	    prm_get_string_int = .false.
	endif

end function prm_get_string_int

!-----------------------------------------------------------------------------

logical function prm_get_string_real(key, value)
!arguments
	character(*),intent(out)	::	key
	real, intent(out)			::	value
!locals
	character(200)				::	str_value


	prm_get_string_real = get_strings(key, str_value)
	read(str_value, fmt=*, iostat=stat) value !read integer from value
	if(stat == 0) then
		prm_get_string_real = .true.
	else 
		prm_get_string_real = .false.
	end if
end function prm_get_string_real

!-----------------------------------------------------------------------------

logical function prm_get_string_real8(key, value)
!arguments
	character(*),intent(out)	::	key
	real(8), intent(out)		::	value
!locals
	character(200)				::	str_value


	prm_get_string_real8 = get_strings(key, str_value)
	read(str_value, fmt=*, iostat=stat) value !read real from str_value

	if (stat == 0) then
		prm_get_string_real8 = .true.
	else
		prm_get_string_real8 = .false.
	endif

end function prm_get_string_real8

!---------------------------------------------------------------------------------
logical function prm_get_int_real(key, value)
!arguments
	integer,intent(out)			::	key
	real, intent(out)			::	value
!locals
	character(200)				::	line


	prm_get_int_real = prm_get_line(line)
	read(line, fmt=*, iostat=stat) key, value !read integer from value
	if(stat == 0) then
		prm_get_int_real = .true.
	else 
		prm_get_int_real = .false.
	end if
end function prm_get_int_real

!---------------------------------------------------------------------------------
logical function prm_get_int_real8(key, value)
!arguments
	integer,intent(out)			::	key
	real(8), intent(out)			::	value
!locals
	character(200)				::	line


	prm_get_int_real8 = prm_get_line(line)
	read(line, fmt=*, iostat=stat) key, value !read integer from value
	if(stat == 0) then
		prm_get_int_real8 = .true.
	else 
		prm_get_int_real8 = .false.
	end if
end function prm_get_int_real8

!----------------------------------------------------------------------------------

logical function prm_get_int_int(key, value)
!arguments
	integer,intent(out)			::	key
	integer, intent(out)		::	value
!locals
	character(200)				::	line


	prm_get_int_int = prm_get_line(line)
	read(line, fmt=*, iostat=stat) key, value !read integer from value
	if(stat == 0) then
		prm_get_int_int = .true.
	else 
		prm_get_int_int = .false.
	end if
end function prm_get_int_int

!--------------------------------------------------------------------------

logical function prm_get_int(value)
!arguments
	integer, intent(out)			::	value
!locals
	character(200)				::	line
	
	prm_get_int = .false.

	if(prm_get_line(line)) then
		read(line, fmt=*, iostat=stat) value !read integer from value
		if(stat == 0) then
			prm_get_int = .true.
		end if
	end if
end function prm_get_int

!----------------------------------------------------------------------------

logical function get_strings(key, value)
!arguments
	character(*),intent(out)	::	key
	character(*), intent(out)	::	value
!locals
	character(200)			::	line
	logical				::	stat

	key = ''
	value = ''
	get_strings = .false.

	stat = prm_get_line(line)
	if(.not. stat) return
	if(.not. split(line, key, value)) return

	get_strings = .true.
end function get_strings

!---------------------------------------------------------------------------------

logical function prm_open_section(section, filename)
!arguments
	character(*),intent(in)	::	section
	character(*),intent(in), optional	::	filename
!locals
	logical						::	used
		
	if(present(filename)) then
		if(filename == PRM_FN) then
			!if the unit has been closed it needs re-opening
			inquire(unit=PRM_U, opened=used)
			if(.not. used) then
				if(.not. prm_open(filename)) then
					prm_open_section = .false.
					return
				end if
			end if
		else
			call prm_close
			if(.not. prm_open(filename)) then
				prm_open_section = .false.
				return
			end if
		end if
	end if
	prm_open_section = find_section(section)	
end function prm_open_section

!---------------------------------------------------------------------------------

subroutine rewind_section
	if(associated(current_sec)) then
		!point to first
		current_lp => current_sec%lp
	else
		call rewind_title
	end if

end subroutine rewind_section

!------------------------------------------------------------------------------

logical function find_section(section)
!arguments
	character(*), intent(in)	::	section
!locals
	character(200)				::	line
	logical						::	rewound
	character(80)				::	ucase_section
	integer						::	sec_len
	integer :: filestat

	type(SECTION_TYPE), pointer::	new_sec

	ucase_section = section
	sec_len = len_trim(section)

	call upcase(ucase_section)

	!loop from current position to end
	new_sec => current_sec
	do while(associated(new_sec))
		if(new_sec%name == ucase_section) goto 100
		new_sec => new_sec%next
	end do
	!loop from start to current position
	new_sec => first_sec
	do while(associated(new_sec) .and. .not. associated(new_sec, current_sec))
		if(new_sec%name == ucase_section) goto 100
		new_sec => new_sec%next
	end do
	find_section = .false.
	return

100	current_sec => new_sec
	current_lp => current_sec%lp
	find_section = .true.

end function find_section

!---------------------------------------------------------------------------------------

logical function prm_get_next_title(title)
!arguments
	character(*), intent(out)	::	title

	title = ''
	prm_get_next_title = .false.

	if(next_title > '') then
		current_title = next_title
		if(load()) then
			title = current_title
			prm_get_next_title = .true.
		end if
	end if

end function prm_get_next_title

!----------------------------------------------------------------------------

subroutine rewind_title
	current_sec => first_sec
	if(associated(current_sec)) then
		current_lp => current_sec%lp
	else
		nullify(current_lp)
	end if
end subroutine rewind_title

!-------------------------------------------------------------------------------

logical function prm_open(filename)
!arguments
	character(*),intent(in)	::	filename

	call prm_close
	PRM_U = freefile()
	open(unit=PRM_U, file = filename, status = 'old', form='formatted', &
		action='read', iostat = stat, access='sequential')
	if(stat /= 0) then
		prm_open = .false.
		PRM_U = 0
	else
		PRM_FN = filename
		if(load()) then
			prm_open = .true.
		else
			prm_open = .false.
			call prm_close
		end if
	end if
end function prm_open

!---------------------------------------------------------------------------------

subroutine clear
!locals
	type(SECTION_TYPE), pointer	::	sp_next
	type(LINE_TYPE), pointer	::	lp_next, lp_current

	!start with first section
	current_sec => first_sec
	!repeat until no more sections
	do while(associated(current_sec))
		sp_next => current_sec%next
		!clear all line entries within the section
		lp_current => current_sec%lp
		do while(associated(lp_current))
			lp_next => lp_current%next
			deallocate(lp_current)
			lp_current => lp_next
		end do
		deallocate(current_sec)
		current_sec => sp_next
	end do
	nullify(first_sec)
end subroutine clear

!--------------------------------------------------------------------------------------

subroutine prm_dump
!locals
	type(SECTION_TYPE), pointer	::	sp_next
	type(LINE_TYPE), pointer	::	lp_next, lp_current

	!start with first section
	call rewind_title
	!repeat until no more sections
	do while(associated(current_sec))
		write(*,100) trim(current_sec%name), current_sec%count

		lp_current => current_sec%lp
		do while(associated(lp_current))
			write(*,200) trim(lp_current%text)
			lp_current => lp_current%next
		end do
		current_sec => current_sec%next
	end do
	call rewind_title

100	format('[',a,'] !',i4,' entries')
200	format(a)
end subroutine prm_dump

!-----------------------------------------------------------------------------------

logical function load()
!locals
	character(len=400)			::	line
	integer						::	stat
	integer						::	end_of_text, comment_start, start, scan_start, comment_pos
	character					::	WS*2, TAB

	TAB = char(9)
	WS(1:1) = ' '
	WS(2:2) = TAB
	
	call clear
	next_title = '' !reset
	section_count = 0
	load = .true.

	do
		line = ''
		read(unit=PRM_U, fmt='(a)', iostat=stat) line
		!check EOF or I/O error
		if(stat > 0) then !read error
			load = .false.
			exit
		elseif(stat < 0) then !EOF
			if(line == '') exit !don't exit if something was on the last line
		end if
		
		start = verify(line, WS)
		if(start == 0) cycle !it was all WS
		!skip blank and comment lines
		if(line(start:start) == '!' .or. &
			line(start:start) == '#' .or. line(start:start) == '*') cycle

		!store next title
		if(line(start:start) == '{') then 
			end_of_text = index(line, '}')
			if(end_of_text > start+1) then
				next_title = line(start+1:end_of_text - 1)
			else
				!invalid title heading
				write(*,'(a)') '>>>>> ERROR: Invalid title heading:'
				write(*,'(a)') line
				load = .false.
			end if
			line = ''
			exit
		elseif(line(start:start) == '[') then !is it a new section?
			end_of_text = index(line, ']')
			if(end_of_text > start+1) then
				!create new section		
				if(section_count == 0) then
					allocate(first_sec)
					current_sec => first_sec
				else
					allocate(current_sec%next)
					current_sec => current_sec%next
				end if
				section_count = section_count +1
				current_sec%name = line(start+1:end_of_text - 1)
				call upcase(current_sec%name)
				current_sec%count = 0
				nullify(current_sec%next) !dissociate pointer
				nullify(current_sec%lp)   !initialise line pointer to null
			else
				!invalid section heading
				write(*,'(a)') '>>>>> ERROR: Invalid section heading:'
				write(*,'(a)') line
				load = .false.
				exit
			end if
		else !if none of the above, then store the line
			!check if in a section, otherwise fail
			if(section_count == 0) then
				load = .false.
				exit
			end if
			if(current_sec%count == 0) then !allocate first
				allocate(current_sec%lp)
				current_lp => current_sec%lp
			else
				allocate(current_lp%next) !allocate subsequent
				current_lp => current_lp%next
			end if
			nullify(current_lp%next) !dissociate pointer
			current_sec%count = current_sec%count + 1
			!find end-of-data in line
			!first avoid finding ! or # within quotes
			!comment_start = scan(line, '!#*')
			scan_start = 1
			comment_start = 0
			do 
				comment_pos = scan(line(scan_start:len_trim(line)), '!#*') + scan_start - 1
				if(comment_pos > scan_start) then
					if(line(comment_pos:comment_pos) == '!' .or. &
						scan(line(comment_pos-1:comment_pos-1), WS) > 0) then
						!it's a real comment
						comment_start = comment_pos
						exit
					else
						!this is # or * after a non-WS character, which in not a comment
						scan_start = comment_pos + 1
					endif
				else
					exit
				endif
			end do
			if(comment_start > 0) then
				!find last non-whitespace
				end_of_text = scan(line(1:comment_start-1), WS, BACK = .TRUE.) - 1
				if(end_of_text <= 0) end_of_text = comment_start - 1
			else
				end_of_text = len_trim(line)+1
			end if

			if(comment_start > 0) then
				!find last non-whitespace
				end_of_text = scan(line(1:comment_start-1), WS, BACK = .TRUE.) - 1
				if(end_of_text <= 0) end_of_text = comment_start - 1
			else
				end_of_text = len_trim(line)+1
			end if
			!store the line without comment/trailing ws
			current_lp%text = line(start:end_of_text)
		end if
	end do
	call rewind_title !reset pointers	
end function load

!------------------------------------------------------------------------------------

subroutine prm_close
	if(PRM_U /= 0) then
		close(PRM_U)
		PRM_U = 0
	end if
	call clear
end subroutine prm_close

!------------------------------------------------------------------------------------

logical function split(line, key, value)
!arguments
	character(*), intent(in)	::	line
	character(*), intent(out)	::	key
	character(*), intent(out)	::	value
	character, save				::	TAB = char(9)
!locals
	integer						::	i
	integer						::	key_start, key_end, value_start, value_end

	split = .false.
	key = ''
	value = ''

!first skip leading ws
	do i = 1, len(line)
		if(line(i:i) /= ' ' .and. line(i:i) /= TAB) exit
	end do
	if(i >= len(line)) return
	key_start = i
!find end of key
	do i = key_start, len(line)
		if(line(i:i) == ' ' .or. line(i:i) == TAB) exit
	end do
	if(i >= len(line)) return
	key_end = i-1
	key = line(key_start:key_end)

!skip ws to start of value
	do i = key_end+2, len(line) !start looking after the first ws character
		if(line(i:i) /= ' ' .and. line(i:i) /= TAB) exit
	end do
	if(i >= len(line)) return
	value_start = i

!find end of value
	do i = value_start, len(line)
		if(line(i:i) == ' ' .or. line(i:i) == TAB .or. line(i:i) == '!' ) exit
	end do
	value_end = i-1
	value = line(value_start:value_end)
	split = .true.

end function split

!------------------------------------------------------------------------------------

logical function splitone(line, key, value)
!arguments
	character(*), intent(in)	::	line
	character(*), intent(out)	::	key
	character(*), intent(out)	::	value
	character, save				::	TAB = char(9)
!locals
	integer						::	i
	integer						::	key_start, key_end, value_start, value_end

	splitone = .false.
	key = ''
	value = ''

!first skip leading ws
	do i = 1, len(line)
		if(line(i:i) /= ' ' .and. line(i:i) /= TAB) exit
	end do
	if(i >= len(line)) return
	key_start = i
!find end of key
	do i = key_start, len(line)
		if(line(i:i) == ' ' .or. line(i:i) == TAB) exit
	end do
	if(i >= len(line)) return
	key_end = i-1
	key = line(key_start:key_end)

!skip ws to start of value
	do i = key_end+2, len(line) !start looking after the first ws character
		if(line(i:i) /= ' ' .and. line(i:i) /= TAB) exit
	end do
	if(i >= len(line)) return
	value_start = i

!find end of data
	do i = value_start, len(line)
		if(line(i:i) == '!' .or. line(i:i) == '#' &
			.or. line(i:i) == '*') exit
	end do
	value_end = i-1
	value = line(value_start:value_end)
	splitone = .true.

end function splitone
	
!------------------------------------------------------------------------------------

end module PRMFILE
