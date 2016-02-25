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

! default distance type that contains squared and normal distance
! for bond operations
TYPE qr_dist
real(kind=prec) :: r,r2
end TYPE qr_dist

! distance type that also contains r6,r12 for vdW
TYPE qr_dist2
real(kind=prec) :: r,r2,r6,r12
end TYPE qr_dist2

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

TYPE(qr_dist) function q_dist(a,b)
! returns distance and squared distance
! args
TYPE(qr_vec) :: a,b
! locals
TYPE(qr_vec) :: temp
temp = qvec_sub(b,a)
q_dist%r2 = temp%x**2 + temp%y**2 + temp%z**2
q_dist%r  = sqrt(q_dist%r2)
end function q_dist

TYPE(qr_dist2) function q_dist2(a,b)
! returns distance, squared distance and higher orders for vdW
! args
TYPE(qr_vec) :: a,b
! locals
TYPE(qr_vec) :: temp
temp = qvec_sub(b,a)
q_dist2%r2  = temp%x**2 + temp%y**2 + temp%z**2
q_dist2%r   = sqrt(q_dist2%r2)
q_dist2%r6  = q_dist2%r2 * q_dist2%r2 * q_dist2%r2
q_dist2%r12 = q_dist2%r6 * q_dist2%r6
end function q_dist2

real(kind=prec) function q_dotprod(a,b)
! returns scalar value of dot product
! implemented the normal fortran dot product
! worked into the new q vector system
! args
TYPE(qr_vec) :: a,b
! locals
real(kind=prec) :: veca(3),vecb(3)
veca(1) = a%x
veca(2) = a%y
veca(3) = a%z
vecb(1) = b%x
vecb(2) = b%y
vecb(3) = b%z
q_dotprod = dot_product(veca,vecb)
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


end module QMATH

