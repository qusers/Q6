! (C) 2014 Uppsala Molekylmekaniska HB, Uppsala, Sweden
! sizes.f90
! by John Marelius
! data storage specifications for all Q programs

module	SIZES

        !MN 2002-11-11
        !Set a nice, portable standard for variables
        integer, parameter  :: wp8 = SELECTED_REAL_KIND(15,307) 

	! STORAGE SPECIFICATIONS FOR Q
	! Change according to 
	! 1)	alignment preferences of your machine
	! 2)	maximum number of atoms desired

	! Max number of atoms is limited only by the size of the variables used as indices
	! to the atom arrays.
	! Here we use 16-bit signed integers for indexes to the atom arrays
	! If more than 2**15-1=32767 atoms are required, change max_atom AND AI
	! integer, parameter			::	AI = 2 !integer(AI) = integer(2)
	!	integer, parameter			::	MAX_AI = 2**(8*AI-1)-1
	! Change to this setting to allow up to 2**31-1 = 2147483647 atoms
	 integer, parameter			::	AI = 4 !integer(AI) = integer(4)

	! Size of integer for flag (1 or 0) variables
	! change this for machines do byte-wise memory access slowly
	integer, parameter			::	FLAG = 1
	integer, parameter			::	MAX_FLAG = 2**(8*FLAG-1)-1

	! Size of integer for small range variables 
	! (atom types, qatom/qangle/qtorsion/qimproper types,...)
	! change this for machines do byte-wise memory access slowly
	integer, parameter			::	TINY = 2
	integer, parameter			::	MAX_TINY = 2**(8*TINY-1)-1	

	! Size of integer for medium range variables (angle/torsion types, ...)
	! change this for machines do doublebyte-wise memory access slowly
	integer, parameter			::	SHRT = 2
	integer, parameter			::	MAX_SHRT = 2**(8*SHRT-1)-1

	! NOTE
	! There are still a couple of hard-coded limitations left. 
	! The most important of these are
	! * maxbndtyp	number of entries in bond type library
	! * maxangtyp	number of entries in angle type library
	! * max_atyp	number of entries in atom type library
	! * max_states	number of different states in perturbations
	! * max_qat		number of Q-atoms


    ! NOTE 2  
	! It is important that you change the types in subroutine 'init_nodes'
	! in md.f90 when you change the types above. Otherwise the datastuctures 
	! sent to the slave nodes will be incorrect. This can be dynamicly 
	! handled with later versions of MPI, but is not supported in the 
	! current version.

end module SIZES

