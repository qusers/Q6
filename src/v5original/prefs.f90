! (C) 2000 Uppsala Molekylmekaniska HB, Uppsala, Sweden

! prefs.f90
! by John Marelius

! preference setting handler

module PREFS

 implicit none

!constants
 character(*), private, parameter :: MODULE_VERSION = '5.01'
 character(*), private, parameter :: MODULE_DATE = '2003-06-02'

 integer, private   :: max_prefs, nprefs
 integer, parameter, private :: default_max = 20
 integer, parameter   ::  PREF_LEN=200
 character(len=2), private :: separator

 type, private :: PREF
  character(len=40)  :: name
  logical     :: is_integer, is_real, is_string
  integer, pointer  :: ival
  real, pointer   :: rval
  character(len=PREF_LEN), pointer  :: sval
 end type PREF
 
 type(PREF), private,allocatable :: p(:)

 !private procedures
 private      :: lookup

contains

subroutine pref_initialize(max)
 !args
 integer, optional   :: max

 if(present(max)) then
  max_prefs = max
 else
  max_prefs = default_max
 end if

 allocate(p(max_prefs))
 p(:)%is_integer = .false.
 p(:)%is_real = .false.
 p(:)%is_string = .false.
 nprefs = 0

 separator(1:1) = ' '
 separator(2:2) = char(9)
end subroutine pref_initialize

logical function pref_add(name, ival, rval, sval)
 !args
 character(*), target :: name
 integer, target, optional::  ival
 real, target, optional::  rval
 character(*), target, optional::  sval

 if(nprefs == max_prefs) then
  pref_add = .false.
 else
  nprefs = nprefs + 1
  p(nprefs)%name = name
  if(present(ival)) then
   p(nprefs)%is_integer = .true.
   p(nprefs)%ival => ival
   pref_add = .true.
  elseif(present(rval)) then
   p(nprefs)%is_real = .true.
   p(nprefs)%rval => rval
   pref_add = .true.
  elseif(present(sval)) then
   p(nprefs)%is_string = .true.
   p(nprefs)%sval => sval
   pref_add = .true.
  else
   pref_add =.false.
  end if
 end if
end function pref_add

subroutine pref_list(heading)
 !args
 character(*), optional  :: heading
 !locals
 integer      :: i
 
 if(present(heading)) then
  write(*,'(a)') heading
 end if

10 format(i2,'.',1x,a, t40, i10)
20 format(i2,'.',1x,a, t40, f10.4)
30 format(i2,'.',1x,a, t40, a)

 do i = 1, nprefs
  if(p(i)%is_integer) then
   write(*,10) i,p(i)%name, p(i)%ival  
  elseif(p(i)%is_real) then
   write(*,20) i,p(i)%name, p(i)%rval  
  elseif(p(i)%is_string) then
   write(*,30) i,p(i)%name, p(i)%sval  
  end if
 end do
end subroutine pref_list

logical function pref_set(string, inval)
 !args
 character(*), intent(in) :: string
 character(*), intent(in), optional :: inval
 
 !locals
 integer      :: i
 character(len=40)   :: name
 character(len=PREF_LEN)   :: value
 integer      :: name_end, readstat
 character(len=PREF_LEN)   :: line

 if(present(inval)) then
  name = string
  value = inval
 else
  !split string into name and value 
  line = adjustl(string)
  name_end = scan(line, separator) - 1
  if(name_end == -1) name_end = len_trim(line)

  name = line(1:name_end)
  if(name_end > len_trim(line)-2) then
   value = ''
  else
   value = line(name_end+2:len_trim(line))
  end if
 end if
 pref_set = .false.
 i = lookup(name)
 if(i > 0) then
  if(p(i)%is_integer) then
   read(value, *, iostat=readstat) p(i)%ival
   if(readstat == 0) pref_set = .true.
  elseif(p(i)%is_real) then
   read(value, *, iostat=readstat) p(i)%rval
   if(readstat == 0) pref_set = .true.
  elseif(p(i)%is_string) then
   p(i)%sval = adjustl(trim(value))
   pref_set = .true.
  end if
  if(.not. pref_set) then
   write(*,900) adjustl(trim(value))
900   format('>>>>> ERROR: Invalid value ',a)
  end if
 end if
end function pref_set

integer function lookup(name)
 !args
 character(*), intent(in) :: name
 !locals
 integer       :: i
 integer       :: readstat

 lookup = 0
 do i = 1, nprefs
  if(p(i)%name == name) then
   lookup = i
   exit
  end if
 end do
 if(lookup == 0) then
  !try numeric
  read(name,*, iostat=readstat) i
  if(readstat == 0 .and. i >= 1 .and. i <= nprefs) then
   lookup = i
  else
   write(*,910) trim(name)
910   format('There is no preference named ',a)
  end if
 end if
end function lookup


end module PREFS
