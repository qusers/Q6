! (C) 2014 Uppsala Molekylmekaniska HB, Uppsala, Sweden
! bonded.f90, based partially on md.f90
! by Johan Åqvist, John Marelius, Anders Kaplan, Isabella Feierberg, Martin Nervall & Martin Almlöf
! new code by Paul Bauer 
! module contains basic calculation of all bonded interaction
! to reduce clutter in md.f90

module BONDED

use QMATH
use SIZES

implicit none
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

TYPE(bond_val) function bond_calc(a,b)
! returns bond distance and vector
! args
TYPE(qr_vec) :: a,b
! locals

bond_calc%a_vec  = b - a
bond_calc%b_vec  = -bond_calc%a_vec
bond_calc%dist = q_sqrt(qvec_square(bond_calc%a_vec))

end function bond_calc

TYPE(bond_val) function box_bond_calc(a,b,boxl,invb)
! same as above, but accounts for periodic effects
! only used in calculation of distance restraints
! boxl -> length of periodic box
! invb -> inverse boxlength
TYPE(qr_vec) :: a,b,boxl,invb
! locals

box_bond_calc%a_vec  = b - a
box_bond_calc%a_vec  = box_bond_calc%a_vec - (boxl*q_nint(box_bond_calc%a_vec*invb))
box_bond_calc%b_vec  = -box_bond_calc%a_vec
box_bond_calc%dist = q_sqrt(qvec_square(box_bond_calc%a_vec))

end function box_bond_calc

TYPE(angl_val) function angle_calc(a,b,c)
! returns angle in radians and force vector for constitute atoms
! uses bond_type to get bond vectors and distances directly
! short summary: calculate dot product between the bond vectors
! to obtain the angle between them
! args
TYPE(qr_vec) :: a,b,c
! locals
TYPE(bond_val)  :: tempab,tempbc
real(kind=prec) :: inv_angl,scalar

tempab = bond_calc(b,a)
tempbc = bond_calc(c,b)

scalar = q_dotprod(tempab%a_vec,tempbc%a_vec)
scalar = scalar/(tempab%dist*tempbc%dist)

if ( scalar .gt.  one ) scalar =  one
if ( scalar .lt. -one ) scalar = -one
angle_calc%angl = q_acos(scalar)

inv_angl = q_sin(angle_calc%angl)
if ( abs(inv_angl) .lt. 1.e-12_prec ) inv_angl = 1.e-12_prec
inv_angl =  -one / inv_angl

angle_calc%a_vec = ( (tempbc%a_vec/(tempab%dist*tempbc%dist)) - &
                   ( (tempab%a_vec/tempab%dist**2) * scalar)) * inv_angl
angle_calc%b_vec = ( (tempab%a_vec/(tempab%dist*tempbc%dist)) - &
                   ( (tempbc%a_vec/tempbc%dist**2) * scalar)) * inv_angl
angle_calc%c_vec = -(angle_calc%a_vec + angle_calc%b_vec)

end function angle_calc

TYPE(angl_val) function box_angle_calc(a,b,c,boxl,invb)
! same function as above, but accounts for periodic effects
! only used in calculation of distance restraints
! boxl -> length of periodic box
! invb -> inverse boxlength
TYPE(qr_vec) :: a,b,c,boxl,invb
TYPE(bond_val)  :: tempab,tempbc
real(kind=prec) :: inv_angl,scalar

tempab = box_bond_calc(b,a,boxl,invb)
tempbc = box_bond_calc(c,b,boxl,invb)

scalar = q_dotprod(tempab%a_vec,tempbc%a_vec)
scalar = scalar/(tempab%dist*tempbc%dist)

if ( scalar .gt.  one ) scalar =  one
if ( scalar .lt. -one ) scalar = -one
box_angle_calc%angl = q_acos(scalar)

inv_angl = q_sin(box_angle_calc%angl)
if ( abs(inv_angl) .lt. 1.e-12_prec ) inv_angl = 1.e-12_prec
inv_angl =  -one / inv_angl

box_angle_calc%a_vec = ( (tempbc%a_vec/(tempab%dist*tempbc%dist)) - &
                       ( (tempab%a_vec/tempab%dist**2) * scalar)) * inv_angl
box_angle_calc%b_vec = ( (tempab%a_vec/(tempab%dist*tempbc%dist)) - &
                       ( (tempbc%a_vec/tempbc%dist**2) * scalar)) * inv_angl
box_angle_calc%c_vec = -(box_angle_calc%a_vec + box_angle_calc%b_vec)

end function box_angle_calc

TYPE(tors_val) function torsion_calc(a,b,c,d)
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
abvec = a - b
bcvec = c - b
cdvec = d - c

! get cross product between the vectors
! nedded for later angle calculation
crossabc =  q_crossprod(abvec,bcvec)
crossbcd = -q_crossprod(bcvec,cdvec)

! absolute number of vector equals distance
abs2_abc = qvec_square(crossabc)
abs2_bcd = qvec_square(crossbcd)

abs_abc  = q_sqrt(abs2_abc)
abs_bcd  = q_sqrt(abs2_bcd)

! angle is dotproduct divided by absolute crossproducts
scalar = q_dotprod(crossabc,crossbcd)
scalar = scalar/(abs_abc*abs_bcd)

! checking if angle can be converted
if ( scalar .gt.  one ) scalar =  one
if ( scalar .lt. -one ) scalar = -one

! final angle
torsion_calc%angl = q_acos(scalar)

! get second crossproduct between two planes
! gives orientation of the angle by calculating
! the normal angle and scalar to bond vector
doublecross  = q_crossprod(crossabc,crossbcd)
sgn          = q_dotprod(bcvec,doublecross)

if ( sgn .lt. zero ) torsion_calc%angl = -torsion_calc%angl

! force vector calculation begins

! first derivative of angle
inv_angl = q_sin(torsion_calc%angl)
if ( abs(inv_angl) .lt. 1.e-12_prec ) inv_angl = 1.e-12_prec
inv_angl =  -one / inv_angl

! get two vectors that decide force directions on all atoms
! first vector is plane vector 2 divided by distance product
! and then substracted by scalar time vector two divided by its
! own absolute value
! vector two is the same with inverted planes
vec1 =  ( (crossbcd/(abs_abc*abs_bcd)) - &
        ( (crossabc/(abs2_abc)) * scalar )) * inv_angl
vec2 =  ( (crossabc/(abs_abc*abs_bcd)) - &
        ( (crossbcd/(abs2_bcd)) * scalar )) * inv_angl

! get the vector between indirect bonded atoms
! dbvec needs inverted orientation because it points
! in the other direction
cavec = bcvec  - abvec
dbvec = -bcvec - cdvec

torsion_calc%a_vec = q_crossprod(bcvec,vec1)
torsion_calc%b_vec = q_crossprod(cavec,vec1) + &
                     q_crossprod(cdvec,vec2)
torsion_calc%c_vec = q_crossprod(dbvec,vec2) - &
                     q_crossprod(abvec,vec1)
torsion_calc%d_vec = q_crossprod(bcvec,vec2)

end function torsion_calc

TYPE(tors_val) function improper_calc(a,b,c,d)
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
abvec = a - b
bcvec = c - b
cdvec = d - c

! get cross product between the vectors
! nedded for later angle calculation
crossabc =  q_crossprod(abvec,bcvec)
crossbcd = -q_crossprod(bcvec,cdvec)
! absolute number of vector equals distance
abs2_abc = qvec_square(crossabc)
abs2_bcd = qvec_square(crossbcd)

abs_abc  = q_sqrt(abs2_abc)
abs_bcd  = q_sqrt(abs2_bcd)

! angle is dotproduct divided by absolute crossproducts
scalar = q_dotprod(crossabc,crossbcd)
scalar = scalar/(abs_abc*abs_bcd)

! checking if angle can be converted
if ( scalar .gt.  one ) scalar =  one
if ( scalar .lt. -one ) scalar = -one

! final angle
improper_calc%angl = q_acos(scalar)

! get second crossproduct between two planes
! gives orientation of the angle by calculating
! the normal angle and scalar to bond vector
doublecross  = q_crossprod(crossabc,crossbcd)
sgn          = q_dotprod(bcvec,doublecross)

if ( sgn .lt. zero ) improper_calc%angl = -improper_calc%angl

! force vector calculation begins

! first derivative of angle
inv_angl = q_sin(improper_calc%angl)
if ( abs(inv_angl) .lt. 1.e-12_prec ) inv_angl = 1.e-12_prec
inv_angl =  -one / inv_angl

! get two vectors that decide force directions on all atoms
! first vector is plane vector 2 divided by distance product
! and then substracted by scalar time vector two divided by its
! own absolute value
! vector two is the same with inverted planes
vec1 = ( (crossbcd/(abs_abc*abs_bcd)) - &
        ( (crossabc/(abs2_abc)) * scalar )) * inv_angl
vec2 = ( (crossabc/(abs_abc*abs_bcd)) - &
        ( (crossbcd/(abs2_bcd)) * scalar )) * inv_angl

! get the vector between indirect bonded atoms
! dbvec needs inverted orientation because it points
! in the other direction
cavec = bcvec - abvec
dbvec = -bcvec - cdvec

improper_calc%a_vec = q_crossprod(bcvec,vec1)
improper_calc%b_vec = q_crossprod(cavec,vec1) + &
                      q_crossprod(cdvec,vec2)
improper_calc%c_vec = q_crossprod(dbvec,vec2) - &
                      q_crossprod(abvec,vec1)
improper_calc%d_vec = q_crossprod(bcvec,vec2)

end function improper_calc



end module BONDED

