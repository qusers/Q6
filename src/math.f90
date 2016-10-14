! (C) 2016 Paul Bauer, Beer Ware Licence
! math.f90
! by Paul Bauer
! general math functions for vector operations used in the program

module  QMATH

use SIZES
! this is now the default vector type for all operations
TYPE qr_vec
real(kind=prec) :: x,y,z
end type qr_vec

! special vector types for reading topologies of different precision
TYPE qr_vecs
real(kind=singleprecision) :: x,y,z
end type qr_vecs

TYPE qr_vecd
real(kind=doubleprecision) :: x,y,z
end type qr_vecd

#ifndef PGI
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

contains

TYPE(qr_vec) function qvec_add(a,b)
! vector addition, std function used later
! args
TYPE(qr_vec) :: a,b
! locals
qvec_add%x = a%x + b%x
qvec_add%y = a%y + b%y
qvec_add%z = a%z + b%z
end function qvec_add

TYPE(qr_vec) function qvec_sub(a,b)
! vector substraction, std function used later
! args
TYPE(qr_vec) :: a,b
! locals
qvec_add%x = a%x - b%x
qvec_add%y = a%y - b%y
qvec_add%z = a%z - b%z
end function qvec_sub

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
temp = qvec_sub(b,a)
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
temp = qvec_sub(b,a)
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
temp = qvec_sub(b,a)
q_dist3%r2  = one/qvec_square(temp)
q_dist3%r   = q_sqrt(q_dist2%r2)
q_dist3%r6  = q_dist2%r2 * q_dist2%r2 * q_dist2%r2
q_dist3%vec = temp
end function q_dist3

real(kind=prec) function q_dist4(a,b)
! returns only squared distance for list generation
! args
TYPE(qr_vec) :: a,b
! locals
TYPE(qr_vec) :: temp
temp = qvec_sub(b,a)
q_dist4 = qvec_square(temp)
end function q_dist4

TYPE(qr_dist5) function q_dist5(a,b)
! returns squared distance and vector for shake
! used in shake dot product calculation
! args
TYPE(qr_vec) :: a,b
! locals
q_dist5%vec = qvec_sub(b,a)
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

TYPE(qr_vec) function q_vecscale(a,b)
! function to scale one vector with another
! by multiplying elements for each vector unit
! args
TYPE(qr_vec) :: a,b
! locals
! none
q_vecscale%x = a%x * b%x
q_vecscale%y = a%y * b%y
q_vecscale%z = a%z * b%z
end function q_vecscale

real(kind=prec) function q_logarithm(a)
! returns results of dlog as real of chosen precision type
! to make compilation independent of variable size
! args
real(kind=prec) :: a
! locals
real(kind=doubleprec) :: temp1,temp2

temp1       = a
temp2       = dlog(temp1)
q_loagrithm = temp2

end function q_logarithm

real(kind=prec) function q_sqrt(a)
! returns result of dsqrt as real of chosen precision type
! to make compilation independent of variable size
! args
real(kind=prec) :: a
! locals
real(kind=doubleprec) :: temp1,temp2

temp1  = a
temp2  = dsqrt(temp1)
q_sqrt = temp2

end function q_sqrt

real(kind=prec) function q_atan(a)
! returns results of datan as real of chosen precision type
! to make compilation independent of variable size
! args
real(kind=prec) :: a
! locals
real(kind=doubleprec) :: temp1,temp2

temp1  = a
temp2  = datan(temp1)
q_atan = temp2

end function q_atan

real(kind=prec) function q_acos(a)
! returns results of dacos as real of chosen precision type
! to make compilation independent of variable size
! args
real(kind=prec) :: a
! locals
real(kind=doubleprec) :: temp1,temp2

temp1  = a
temp2  = dacos(temp1)
q_atan = temp2

end function q_acos

real(kind=prec) function q_cos(a)
! returns results of dcos as real of chosen precision type
! to make compilation independent of variable size
! args
real(kind=prec) :: a
! locals
real(kind=doubleprec) :: temp1,temp2

temp1  = a
temp2  = dcos(temp1)
q_atan = temp2

end function q_cos

real(kind=prec) function q_sin(a)
! returns results of dsin as real of chosen precision type
! to make compilation independent of variable size
! args
real(kind=prec) :: a
! locals
real(kind=doubleprec) :: temp1,temp2

temp1  = a
temp2  = dsin(temp1)
q_atan = temp2

end function q_sin


subroutine math_initialize
        pi = 4.0_prec*q_atan(one)
        deg2rad = pi/180.0_prec
        rad2deg = 180.0_prec/pi
end subroutine math_initialize

end module QMATH

