! (C) 2014 Uppsala Molekylmekaniska HB, Uppsala, Sweden
! bonded.f90, based partially on md.f90
! by Johan Åqvist, John Marelius, Anders Kaplan, Isabella Feierberg, Martin Nervall & Martin Almlöf
! new code by Paul Bauer 
! module contains basic calculation of all bonded interaction
! to reduce clutter in md.f90

module BONDED

use QMATH

! variables used in this subroutine

TYPE bond_val
real(kind=prec) :: dist
TYPE(qr_vec)    :: a_vec,b_vec
end TYPE bond_val

TYPE angl_val
real(kind=prec) :: angl
TYPE(qr_vec)    :: a_vec,b_vec,c_vec
end TYPE angl_val

TYPE tors_val
real(kind=prec) :: angl
TYPE(qr_vec)    :: a_vec,b_vec,c_vec,d_vec
end TYPE tors_val

contains

! functions that can be called from anywhere

TYPE(bond_val) function bond(a,b)
! returns bond distance and vector
! args
TYPE(qr_vec) :: a,b
! locals

bond_val%a_vec  = qvec_sub(b,a)
bond_val%b_vec  = -bond_val%a_vec
bond_val%dist = sqrt(bond_val%a_vec%x**2 + bond_val%a_vec%y**2 + bond_val%a_vec%z**2)

end function bond

TYPE(angl_val) function angle(a,b,c)
! returns angle in radians and force vector for constitute atoms
! uses bond_type to get bond vectors and distances directly
! short summary: calculate dot product between the bond vectors
! to obtain the angle between them
! args
TYPE(qr_vec) :: a,b,c
! locals
TYPE(bond_val)  :: tempab,tempbc
real(kind=prec) :: inv_angl,scalar

tempab = bond(a,b)
tempbc = bond(b,c)

scalar = q_dotprod(tempab%a_vec,tempbc%a_vec)
scalar = scalar/(tempab%dist*tempbc%dist)

if ( scalar .gt.  one ) scalar =  one
if ( scalar .lt. -one ) scalar = -one
angle%angl = acos(scalar)

inv_angl = sin(angle%angl)
if ( abs(inv_angl) .lt. 1.e-12_prec ) inv_angl = 1.e-12_prec
inv_angl =  -one / inv_angl

angle%a_vec = inv_ang * ( (tempbc%a_vec/(tempab%dist*tempbc%dist)) - &
               (scalar * tempab%a_vec/tempab%dist**2))
angle%b_vec = inv_ang * ( (tempab%a_vec/(tempab%dist*tempbc%dist)) - &
               (scalar * tempbc%a_vec/tempbc%dist**2))
angle%c_vec = -(angle%a_vec + angle%b_vec)

end function angle

TYPE(tors_val) function torsion(a,b,c,d)
! returns torsion angle and force vectors for constitute atoms
! for a derivation, please check your local copy of a vector 
! math book
! short summary: calculate planes defined by the bond vector
! get orientation by calculating the normal vector defined by those 
! planes vectors and creating the scalar product between the 
! center bond vector and the normal vector
! torsion angle is scalar product of plane vectors divided by 
! absolute value of those vectors
! args
TYPE(qr_vec) :: a,b,c,d
! locals
TYPE(qr_vec)    :: abvec,bcvec,cdvec
TYPE(qr_vec)    :: crossabc,crossbcd,doublecross
TYPE(qr_vec)    :: cavec,dbvec
real(kind=prec) :: scalar,inv_angl,sgn
real(kind=prec) :: abs_abc,abs_bcd,abs2_abc,abs2_bcd
TYPE(qr_vec)    :: vec1,vec2


! get length of the individual bond vectors
abvec = qvec_sub(a,b)
bcvec = qvec_sub(b,c)
cdvec = qvec_sub(c,d)

! get cross product between the vectors
! nedded for later angle calculation
crossabc =  q_crossprod(abvec,bcvec)
crossbcd = -q_crossprod(bcvec,cdvec)

! absolute number of vector equals distance
abs2_abc = (crossabc%x**2 + crossabc%y**2 + crossabc%z**2)
abs2_bcd = (crossbcd%x**2 + crossbcd%y**2 + crossbcd%z**2)

abs_abc  = sqrt(abs2_abc)
abs_bcd  = sqrt(abs2_bcd)

! angle is dotproduct divided by absolute crossproducts
scalar = q_dotprod(crossabc,crossbcd)
scalar = scalar/(abs_abc*abs_bcd)

! checking if angle can be converted
if ( scalar .gt.  one ) scalar =  one
if ( scalar .lt. -one ) scalar = -one

! final angle
torsion%angl = acos(scalar)

! get second crossproduct between two planes
! gives orientation of the angle by calculating
! the normal angle and scalar to bond vector
doublecross  = q_crossprod(crossabc,crossbcd)
sgn          = q_dotprod(bcvec,doublecross)

if ( sgn .lt. zero ) torsion%angl = -torsion%angl

! force vector calculation begins

! first derivative of angle
inv_angl = sin(torsion%angl)
if ( abs(inv_angl) .lt. 1.e-12_prec ) inv_angl = 1.e-12_prec
inv_angl =  -one / inv_angl

! get two vectors that decide force directions on all atoms
! first vector is plane vector 2 divided by distance product
! and then substracted by scalar time vector two divided by its
! own absolute value
! vector two is the same with inverted planes
vec1 = inv_angl * ( (crossbcd/(abs_abc*abs_bcd)) - &
        ( scalar *crossabc/(abs2_abc)))
vec2 = inv_angl * ( (crossabc/(abs_abc*abs_bcd)) - &
        ( scalar *crossbcd/(abs2_bcd)))

! get the vector between indirect bonded atoms
! dbvec needs inverted orientation because it points
! in the other direction
cavec = qvec_sub(bcvec,abvec)
dbvec = qvec_sub(-bcvec,cdvec)

torsion%a_vec = q_crossprod(bcvec,vec1)
torsion%b_vec = qvec_add( &
                 q_crossprod(cavec,vec1) , &
                 q_crossprod(cdvec,vec2))
torsion%c_vec = qvec_sub( &
                 q_crossprod(dbvec,vec2) , &
                 q_crossprod(abvec,vec1))
torsion%d_vec = q_crossprod(bcvec,vec2)

end function torsion

TYPE(tors_val) function improper(a,b,c,d)
! returns improper torsion angle and force vectors for constitute atoms
! for a derivation, please check your local copy of a vector 
! math book
! short summary: calculate planes defined by the bond vector
! args
TYPE(qr_vec) :: a,b,c,d
! locals
TYPE(qr_vec)    :: abvec,bcvec,cdvec
TYPE(qr_vec)    :: crossabc,crossbcd,doublecross
TYPE(qr_vec)    :: cavec,dbvec
real(kind=prec) :: scalar,inv_angl,sgn
real(kind=prec) :: abs_abc,abs_bcd,abs2_abc,abs2_bcd
TYPE(qr_vec)    :: vec1,vec2

! get length of the individual bond vectors
abvec = qvec_sub(a,b)
bcvec = qvec_sub(c,b)
cdvec = qvec_sub(d,c)

! get cross product between the vectors
! nedded for later angle calculation
crossabc =  q_crossprod(abvec,bcvec)
crossbcd = -q_crossprod(bcvec,cdvec)

! absolute number of vector equals distance
abs2_abc = (crossabc%x**2 + crossabc%y**2 + crossabc%z**2)
abs2_bcd = (crossbcd%x**2 + crossbcd%y**2 + crossbcd%z**2)

abs_abc  = sqrt(abs2_abc)
abs_bcd  = sqrt(abs2_bcd)

! angle is dotproduct divided by absolute crossproducts
scalar = q_dotprod(crossabc,crossbcd)
scalar = scalar/(abs_abc*abs_bcd)

! checking if angle can be converted
if ( scalar .gt.  one ) scalar =  one
if ( scalar .lt. -one ) scalar = -one

! final angle
improper%angl = acos(scalar)

! get second crossproduct between two planes
! gives orientation of the angle by calculating
! the normal angle and scalar to bond vector
doublecross  = q_crossprod(crossabc,crossbcd)
sgn          = q_dotprod(bcvec,doublecross)

if ( sgn .lt. zero ) improper%angl = -improper%angl

! force vector calculation begins

! first derivative of angle
inv_angl = sin(improper%angl)
if ( abs(inv_angl) .lt. 1.e-12_prec ) inv_angl = 1.e-12_prec
inv_angl =  -one / inv_angl

! get two vectors that decide force directions on all atoms
! first vector is plane vector 2 divided by distance product
! and then substracted by scalar time vector two divided by its
! own absolute value
! vector two is the same with inverted planes
vec1 = inv_angl * ( (crossbcd/(abs_abc*abs_bcd)) - &
        ( scalar *crossabc/(abs2_abc)))
vec2 = inv_angl * ( (crossabc/(abs_abc*abs_bcd)) - &
        ( scalar *crossbcd/(abs2_bcd)))

! get the vector between indirect bonded atoms
! dbvec needs inverted orientation because it points
! in the other direction
cavec = qvec_sub(bcvec,abvec)
dbvec = qvec_sub(-bcvec,cdvec)

improper%a_vec = q_crossprod(bcvec,vec1)
improper%b_vec = qvec_add( &
                 q_crossprod(cavec,vec1) , &
                 q_crossprod(cdvec,vec2))
improper%c_vec = qvec_sub( &
                 q_crossprod(dbvec,vec2) , &
                 q_crossprod(abvec,vec1))
improper%d_vec = q_crossprod(bcvec,vec2)

end function improper



end module BONDED

