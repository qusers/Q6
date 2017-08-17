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
! math.f90
! by Paul Bauer
! general math functions for vector operations used in the program

module  QMATH

use SIZES

implicit none

! special vector types for reading topologies of different precision
TYPE qr_vecs
real(kind=singleprecision) :: x,y,z
end type qr_vecs

TYPE qr_vecd
real(kind=doubleprecision) :: x,y,z
end type qr_vecd

#ifdef HAVEQUAD
TYPE qr_vecq
real(kind=quadprecision) :: x,y,z
end type qr_vecq
#endif

! default distance type that contains squared and normal distance
! for bond operations
TYPE qr_dist
real(kind=prec) :: r,r2
TYPE(qr_vec)    :: vec
end TYPE qr_dist

! distance type that also contains r6,r12 for vdW
TYPE qr_dist2
real(kind=prec) :: r,r2,r6,r12
TYPE(qr_vec)    :: vec
end TYPE qr_dist2

! and one more that excludes r12, for q nonbonded that need additional calculation of soft core
TYPE qr_dist3
real(kind=prec) :: r,r2,r6
TYPE(qr_vec)    :: vec
end TYPE qr_dist3

! distance and vector for shake calc
! messed up number, so it is now 5
TYPE qr_dist5
real(kind=prec) :: r2
TYPE(qr_vec)    :: vec
end TYPE qr_dist5

real(kind=prec) :: pi, deg2rad, rad2deg

! type to calculate simple statistics
TYPE STAT_TYPE
        real(kind=prec)         :: ave,dev
end TYPE STAT_TYPE

interface operator(+)
        module procedure qvec_add,qarray_add,qvec_realadd
end interface

interface operator(-)
        module procedure qvec_sub,qvec_negate,qarray_sub,qvec_realsub
end interface

interface operator(*)
        module procedure q_realscale,q_arrayscale,q_vecscale,q_realscale2,q_arrayvecscale
end interface

interface operator(/)
        module procedure q_realdiv,q_arraydiv,q_realdiv2,q_vecdiv
end interface

public :: assignment(=)

interface assignment(=)
        module procedure qvec_assign,qarray_assign
end interface

contains

function qvec_add(a,b)
! vector addition, std function used later
! args
TYPE(qr_vec), INTENT (IN) :: a,b 
TYPE(qr_vec) :: qvec_add 
! locals
qvec_add%x = a%x + b%x
qvec_add%y = a%y + b%y
qvec_add%z = a%z + b%z
end function qvec_add

subroutine qvec_assign(a,b)
! vector addition, std function used later
! args
TYPE(qr_vec), INTENT (INOUT) :: a
real(kind=prec),INTENT(IN)   :: b
! locals
a%x = b
a%y = b
a%z = b
end subroutine qvec_assign

subroutine qarray_assign(a,b)
! vector addition, std function used later
! args
TYPE(qr_vec), INTENT (INOUT) :: a(:)
real(kind=prec),INTENT(IN)   :: b
! locals
a(:)%x = b
a(:)%y = b
a(:)%z = b

end subroutine qarray_assign


function qvec_realadd(a,b)
! vector scaling by constant number, operator
! args
TYPE(qr_vec),INTENT (IN)        :: a
real(kind=prec), INTENT(IN)     :: b
TYPE(qr_vec) :: qvec_realadd
qvec_realadd%x = a%x + b
qvec_realadd%y = a%y + b
qvec_realadd%z = a%z + b

end function qvec_realadd

function qarray_add(a,b)
! vector addition for arrays of vectors
TYPE(qr_vec), INTENT(in) :: a(:),b(SIZE(a))
TYPE(qr_vec)             :: qarray_add(SIZE(a))

qarray_add(:)%x = a(:)%x + b(:)%x
qarray_add(:)%y = a(:)%y + b(:)%y
qarray_add(:)%z = a(:)%z + b(:)%z

end function qarray_add

function qvec_sub(a,b)
! vector substraction, std function used later
! args
TYPE(qr_vec), INTENT (IN) :: a,b
TYPE(qr_vec) :: qvec_sub 
! locals
qvec_sub%x = a%x - b%x
qvec_sub%y = a%y - b%y
qvec_sub%z = a%z - b%z
end function qvec_sub

function qvec_realsub(a,b)
! vector scaling by constant number, operator
! args
TYPE(qr_vec),INTENT(IN)         :: a
real(kind=prec),INTENT(IN)      :: b
TYPE(qr_vec)    :: qvec_realsub
! locals
qvec_realsub%x = a%x - b
qvec_realsub%y = a%y - b
qvec_realsub%z = a%z - b

end function qvec_realsub

function qarray_sub(a,b)
! substract vector arrays
! args
TYPE(qr_vec),INTENT(IN) :: a(:), b(SIZE(a))
TYPE(qr_vec)            :: qarray_sub(SIZE(a))

qarray_sub(:)%x = a(:)%x - b(:)%x
qarray_sub(:)%y = a(:)%y - b(:)%y
qarray_sub(:)%z = a(:)%z - b(:)%z

end function qarray_sub
function qvec_negate(a)
! function negates current vector, for operator
! args
TYPE(qr_vec),INTENT (IN) :: a
TYPE(qr_vec) :: qvec_negate
! locals
qvec_negate%x = -a%x
qvec_negate%y = -a%y
qvec_negate%z = -a%z
end function qvec_negate

real(kind=prec) function qvec_square(a)
! returns only square of vector
! used in other functions
! args
TYPE(qr_vec) :: a
qvec_square = a%x**2 + a%y**2 + a%z**2
end function qvec_square

TYPE(qr_dist) function q_dist(a,b)
! returns distance and squared distance
! args
TYPE(qr_vec) :: a,b
! locals
TYPE(qr_vec) :: temp
temp = b - a
q_dist%r2  = one/qvec_square(temp)
q_dist%r   = q_sqrt(q_dist%r2)
q_dist%vec = temp
end function q_dist

TYPE(qr_dist2) function q_dist2(a,b)
! returns distance, squared distance and higher orders for vdW
! args
TYPE(qr_vec) :: a,b
! locals
TYPE(qr_vec) :: temp
temp = b - a
q_dist2%r2  = one/qvec_square(temp)
q_dist2%r   = q_sqrt(q_dist2%r2)
q_dist2%r6  = q_dist2%r2 * q_dist2%r2 * q_dist2%r2
q_dist2%r12 = q_dist2%r6 * q_dist2%r6
q_dist2%vec = temp
end function q_dist2

TYPE(qr_dist3) function q_dist3(a,b)
! returns distance, squared distance and higher orders for vdW
! excluding r12 that is recalculated for q soft cores
! args
TYPE(qr_vec) :: a,b
! locals
TYPE(qr_vec) :: temp
temp = b - a
q_dist3%r2  = one/qvec_square(temp)
q_dist3%r   = q_sqrt(q_dist3%r2)
q_dist3%r6  = q_dist3%r2 * q_dist3%r2 * q_dist3%r2
q_dist3%vec = temp
end function q_dist3

real(kind=prec) function q_dist4(a,b)
! returns only squared distance for list generation
! args
TYPE(qr_vec) :: a,b
! locals
TYPE(qr_vec) :: temp
temp = b - a
q_dist4 = qvec_square(temp)
end function q_dist4

TYPE(qr_dist5) function q_dist5(a,b)
! returns squared distance and vector for shake
! used in shake dot product calculation
! args
TYPE(qr_vec) :: a,b
! locals
q_dist5%vec = b - a
q_dist5%r2  = qvec_square(q_dist5%vec)
end function q_dist5


real(kind=prec) function q_dotprod(a,b)
! returns scalar value of dot product
! implemented the normal fortran dot product
! worked into the new q vector system
! args
TYPE(qr_vec) :: a,b
! locals
! fortran intrinsic return values are not defined, need to calculate by hand
q_dotprod = a%x*b%x + a%y*b%y + a%z*b%z
end function q_dotprod

TYPE(qr_vec) function q_crossprod(a,b)
! returns vector with cross product
! for use with q vector arrays
! args
TYPE(qr_vec) :: a,b
! locals

q_crossprod%x = a%y * b%z - a%z * b%y
q_crossprod%y = a%z * b%x - a%x * b%z
q_crossprod%z = a%x * b%y - a%y * b%x

end function q_crossprod

function q_realscale(a,b)
! function to scale qr_vec types with real numbers
! args
TYPE(qr_vec), INTENT(IN) :: a
real(kind=prec), INTENT(IN) :: b
TYPE(qr_vec) :: q_realscale
q_realscale%x = a%x * b
q_realscale%y = a%y * b
q_realscale%z = a%z * b
end function q_realscale

function q_realscale2(b,a)
! function to scale qr_vec types with real numbers
! args
TYPE(qr_vec), INTENT(IN) :: a
real(kind=prec), INTENT(IN) :: b
TYPE(qr_vec) :: q_realscale2
q_realscale2%x = a%x * b
q_realscale2%y = a%y * b
q_realscale2%z = a%z * b
end function q_realscale2

function q_realdiv(a,b)
! function to divide q vectors, inverse of mult
TYPE(qr_vec), INTENT(IN) :: a
real(kind=prec), INTENT(IN) :: b
TYPE(qr_vec) :: q_realdiv
q_realdiv%x = a%x / b
q_realdiv%y = a%y / b
q_realdiv%z = a%z / b
end function q_realdiv

function q_realdiv2(b,a)
! function to divide q vectors, inverse of mult
TYPE(qr_vec), INTENT(IN) :: a
real(kind=prec), INTENT(IN) :: b
TYPE(qr_vec) :: q_realdiv2
q_realdiv2%x =  b / a%x 
q_realdiv2%y =  b / a%y 
q_realdiv2%z =  b / a%z 
end function q_realdiv2


function q_arrayscale(a,b)
! function to scale qr_vec types with real numbers
! args
TYPE(qr_vec), INTENT(IN) :: a(:)
real(kind=prec), INTENT(IN) :: b
TYPE(qr_vec) :: q_arrayscale(SIZE(a))
q_arrayscale(:)%x = a(:)%x * b
q_arrayscale(:)%y = a(:)%y * b
q_arrayscale(:)%z = a(:)%z * b
end function q_arrayscale

function q_arraydiv(a,b)
! function to divide q vectors, inverse of mult
TYPE(qr_vec), INTENT(IN) :: a(:)
real(kind=prec), INTENT(IN) :: b
TYPE(qr_vec) :: q_arraydiv(SIZE(a))
q_arraydiv(:)%x = a(:)%x / b
q_arraydiv(:)%y = a(:)%y / b
q_arraydiv(:)%z = a(:)%z / b
end function q_arraydiv


function q_vecscale(a,b)
! function to scale one vector with another
! by multiplying elements for each vector unit
! args
TYPE(qr_vec), INTENT(IN) :: a,b
TYPE(qr_vec) ::  q_vecscale
! locals
! none
q_vecscale%x = a%x * b%x
q_vecscale%y = a%y * b%y
q_vecscale%z = a%z * b%z
end function q_vecscale

function q_arrayvecscale(a,b)
! and the same function as above for arrays of qr_vec objects ...
! args
TYPE(qr_vec),INTENT(IN)         :: a(:),b(SIZE(a))
TYPE(qr_vec)                    :: q_arrayvecscale(SIZE(a))
q_arrayvecscale%x = a(:)%x * b(:)%x
q_arrayvecscale%y = a(:)%y * b(:)%y
q_arrayvecscale%z = a(:)%z * b(:)%z

end function q_arrayvecscale

function q_vecdiv(a,b)
! divides components of one vector by those of another
! for operator
TYPE(qr_vec), INTENT(IN) :: a,b
TYPE(qr_vec) :: q_vecdiv
q_vecdiv%x = a%x / b%x
q_vecdiv%y = a%y / b%y
q_vecdiv%z = a%z / b%z
end function q_vecdiv

real(kind=prec) function q_logarithm(a)
! returns results of dlog as real of chosen precision type
! to make compilation independent of variable size
! args
real(kind=prec) :: a
#ifdef QSINGLE
q_logarithm = log(a)
#endif
#ifdef QDOUBLE
q_logarithm = dlog(a)
#endif
#ifdef QQUADRUPLE
q_logarithm = dlog(real(a,kind=doubleprecision))
#endif
end function q_logarithm

real(kind=prec) function q_log2(a)
! returns log to base two of any value
! args
real(kind=prec) :: a
q_log2 = q_logarithm(a) / q_logarithm(2.0_prec)
end function q_log2

real(kind=prec) function q_sqrt(a)
! returns result of dsqrt as real of chosen precision type
! to make compilation independent of variable size
! args
real(kind=prec) :: a
#ifdef QSINGLE
q_sqrt = sqrt(a)
#endif
#ifdef QDOUBLE
q_sqrt = dsqrt(a)
#endif
#ifdef QQUADRUPLE
q_sqrt = dsqrt(real(a,kind=doubleprecision))
#endif
end function q_sqrt

function q_sqrt2(a)
! returns result of dsqrt as real of chosen precision type
! to make compilation independent of variable size
! args
real(kind=prec) :: a(:)
reaL(kind=prec) :: q_sqrt2(size(a))
#ifdef QSINGLE
q_sqrt2(:)  = sqrt(a)
#endif
#ifdef QDOUBLE
q_sqrt2(:)  = dsqrt(a)
#endif
#ifdef QQUADRUPLE
q_sqrt2(:)  = dsqrt(real(a(:),kind=doubleprecision))
#endif
end function q_sqrt2

real(kind=prec) function q_atan(a)
! returns results of datan as real of chosen precision type
! to make compilation independent of variable size
! args
real(kind=prec) :: a
#ifdef QSINGLE
q_atan  = atan(a)
#endif
#ifdef QDOUBLE
q_atan  = datan(a)
#endif
#ifdef QQUADRUPLE
q_atan  = datan(real(a,kind=doubleprecision))
#endif
end function q_atan

real(kind=prec) function q_acos(a)
! returns results of dacos as real of chosen precision type
! to make compilation independent of variable size
! args
real(kind=prec) :: a
#ifdef QSINGLE
q_acos  = acos(a)
#endif
#ifdef QDOUBLE
q_acos  = dacos(a)
#endif
#ifdef QQUADRUPLE
q_acos  = dacos(real(a,kind=doubleprecision))
#endif
end function q_acos

real(kind=prec) function q_cos(a)
! returns results of dcos as real of chosen precision type
! to make compilation independent of variable size
! args
real(kind=prec) :: a
#ifdef QSINGLE
q_cos  = cos(a)
#endif
#ifdef QDOUBLE
q_cos  = dcos(a)
#endif
#ifdef QQUADRUPLE
q_cos  = dcos(real(a,kind=doubleprecision))
#endif
end function q_cos

real(kind=prec) function q_sin(a)
! returns results of dsin as real of chosen precision type
! to make compilation independent of variable size
! args
real(kind=prec) :: a
#ifdef QSINGLE
q_sin  = sin(a)
#endif
#ifdef QDOUBLE
q_sin  = dsin(a)
#endif
#ifdef QQUADRUPLE
q_sin  = dsin(real(a,kind=doubleprecision))
#endif
end function q_sin

real(kind=prec) function q_exp(a)
! returns exponent of value as correct q precision type real
! args
real(kind=prec) :: a
#ifdef QSINGLE
q_exp = exp(a)
#endif
#ifdef QDOUBLE
q_exp = dexp(a)
#endif
#ifdef QQUADRUPLE
q_exp = dexp(real(a,kind=doubleprecision))
#endif
end function q_exp

function q_exp2(a)
! returns exponent of value as correct q precision type real
! args
real(kind=prec) :: a(:)
real(kind=prec) :: q_exp2(size(a))
#ifdef QSINGLE
q_exp2 = exp(a)
#endif
#ifdef QDOUBLE
q_exp2 = dexp(a)
#endif
#ifdef QQUADRUPLE
q_exp2 = dexp(real(a,kind=doubleprecision))
#endif
end function q_exp2


real(kind=prec) function q_epsilon(a)
! from gcc documentation
! returns the smallest number E of the same kind as X such that 1 + E > 1. 
! args
real(kind=prec) :: a
#ifdef QSINGLE
q_epsilon = epsilon(a)
#endif
#ifdef QDOUBLE
q_epsilon = epsilon(a)
#endif
#ifdef QQUADRUPLE
q_epsilon = epsilon(real(a,kind=doubleprecision))
#endif
end function q_epsilon

TYPE(qr_vec) function q_nint(a)
! returns the vector of nearest integers for a given real qr_vec structure
! used when calculating PBC stuff
! args
TYPE(qr_vec) :: a
q_nint%x = nint(a%x)
q_nint%y = nint(a%y)
q_nint%z = nint(a%z)
end function q_nint

function q_ceiling(a)
! returns next highest integer for all components of a qr_vec object
! args
TYPE(qr_vec) :: a
! locals
integer,dimension(3) :: q_ceiling
q_ceiling(1) = ceiling(a%x)
q_ceiling(2) = ceiling(a%y)
q_ceiling(3) = ceiling(a%z)

end function q_ceiling

real(kind=prec) function q_abs(a)
! returns absolute number in requested precision
! always calculated as dabs
! args
real(kind=prec) :: a
#ifdef QSINGLE
q_abs = abs(a)
#endif
#ifdef QDOUBLE
q_abs = dabs(a)
#endif
#ifdef QQUADRUPLE
q_abs = dabs(real(a,kind=doubleprecision))
#endif
end function q_abs

function q_abs2(a)
! returns absolute number in requested precision
! always calculated as dabs
! args
real(kind=prec) :: a(:)
real(kind=prec) :: q_abs2(size(a))
#ifdef QSINGLE
q_abs2(:) = abs(a(:))
#endif
#ifdef QDOUBLE
q_abs2(:) = dabs(a(:))
#endif
#ifdef QQUADRUPLE
q_abs2(:) = dabs(real(a(:),kind=doubleprecision))
#endif
end function q_abs2

subroutine math_initialize
        pi = 4.0_prec*q_atan(one)
        deg2rad = pi/180.0_prec
        rad2deg = 180.0_prec/pi
end subroutine math_initialize

! function to compute standard deviations from sum and sum of squares
TYPE(STAT_TYPE) function get_stddev(val,val2,num)
! arguments
! sum and sum of square
real(kind=prec)         :: val,val2
! number of points
integer                 :: num
! locals
real(kind=prec)         :: tmp

tmp = val/real(num,kind=prec)
get_stddev%ave = tmp
tmp = val2/num - tmp**2
if(tmp.gt.zero) then
        get_stddev%dev = tmp
else
        get_stddev%dev = zero
end if
end function get_stddev


#ifdef DEBUG
!function that returns sum of a given qr_vec array
!useful for printing out total force vectors during debugging
subroutine q_vecsum(a,b,c)
! arguments
TYPE(qr_vec)            :: a(:)
TYPE(qr_vec)            :: vecsum
integer                 :: b
character(*)            :: c
! locals
integer                 :: i

666     format(a,a30,3f20.16)

vecsum = zero
do i = 1 , b
        vecsum = vecsum + a(i)
end do

write(*,666) 'Total force after ',c,vecsum


end subroutine q_vecsum
#endif

end module QMATH

