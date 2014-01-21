!	(C) 2000 Uppsala Molekylmekaniska HB, Uppsala, Sweden

!	index.f90
!	by John Marelius

!	string-to-integer index lookup table and functions

module indexer

implicit none

!constants
	character(*), private, parameter	::	MODULE_VERSION = '5.01'
	character(*), private, parameter	::	MODULE_DATE = '2003-06-02'

integer, private ::				count = 0, top = 0
integer, parameter ::	KEYLENGTH = 8
integer, parameter, private	::	RESIZE_INCREMENT = 10
type indexentry
	character(len=KEYLENGTH)::	key
	integer					::	ndx
end type indexentry

type(indexentry), pointer, private :: ndx(:)

integer, private				::	memstat

!declare private procedures
private grow

contains
subroutine index_create(size)
!arguments
	integer, optional			::	size
!locals
	
	call index_clear
	top = RESIZE_INCREMENT !use default starting size
	if(present(size)) then
		if(size > 0) then
			top = size
		end if
	end if

	allocate(ndx(top))
end subroutine index_create

subroutine index_clear

	count = 0
	top = 0
	deallocate(ndx, stat=memstat) !ignore error if not allocated

end subroutine index_clear

subroutine index_shutdown
	call index_clear
end subroutine index_shutdown

subroutine index_startup
end subroutine index_startup

subroutine index_resize(size)
!arguments
	integer						::	size
!locals
	integer						::	mintop
	type(indexentry), pointer :: new_ndx(:)

	if(size <= 0) then
		call index_clear
	else
		!allocate new
		allocate(new_ndx(size))
		!copy
		mintop = min(top, size)
		new_ndx(1:mintop) = ndx(1:mintop)
		!get rid of old
		deallocate(ndx, stat=memstat)
		!point ndx to new_ndx
		ndx => new_ndx
		top = size
	end if

end subroutine index_resize

subroutine grow
	call index_resize(top+RESIZE_INCREMENT)
end subroutine grow

logical function index_add(key, index)
!arguments
	character(*)				::	key
	integer						::	index

!locals
	integer						::	i, j

	do i = 1, count
		if(ndx(i)%key == key) then !cannot redefine
			index_add = .false.
			return
		elseif(lgt(ndx(i)%key, key)) then 
			!grow if necessary
			if(count == top) call grow
			!make a slot by moving subsequent entries
			do j = count, i, -1
				ndx(j+1) = ndx(j)
			end do
			exit
		end if
	end do
	if(count == top) call grow !grow also if inserting last
	ndx(i)%key = key
	ndx(i)%ndx = index
	count = count + 1
	index_add = .true.
end function index_add

logical function index_alias(alias, key)
!arguments
	character(*)				::	alias, key
!locals
	integer						::	index
	index_alias = .false.
	if(index_get(key, index)) then
		if(index_add(alias, index)) then
			index_alias = .true.
		end if
	end if

end function index_alias

logical function index_get(key, i, allow_wildcard)
!arguments
	character(*), intent(in)	::	key
	integer, intent(out)		::	i
	logical, optional			::	allow_wildcard
!locals
	integer						::	hi, lo, mid

!perform standard binary search	
!(Modified from N. C. Shammas/Secrets of the Visual C++ Masters)
	index_get = .false.
	i = 0
	
	if(present(allow_wildcard)) then
		if(allow_wildcard .and. (key == '?' .or. key == '0')) then
			index_get = .true.
			return
		endif
	endif
	lo = 1
	hi = count
	do while(lo <= hi)
		mid = (lo + hi) / 2
		if(ndx(mid)%key == key) then
			i = ndx(mid)%ndx
			index_get = .true.
			return
		elseif(lle(key, ndx(mid)%key)) then 
			hi = mid - 1
		else
			lo = mid + 1
		end if
	end do

end function index_get
	
end module indexer
