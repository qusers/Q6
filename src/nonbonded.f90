! (C) 2014 Uppsala Molekylmekaniska HB, Uppsala, Sweden
! nonbonded.f90 
! based on md.f90
! by Johan Åqvist, John Marelius, Anders Kaplan, Isabella Feierberg, Martin Nervall & Martin Almlöf
! subroutines needed for non bonded energy calculations
! by Paul Bauer

module NONBONDED

! used modules
use SIZES
use GLOBALS
use QMATH
use QATOM
implicit none



contains

!-------------------------------------------------------------------------
! nonbond routines are now consolidated to just sphere/box versions
! because combination rules are applied earlier
! all routines here get called by the higher up routines with the proper lists
! to calculate on and return the energies and forces for a given calculation

TYPE(ENERET_TYPE) function nbe(nb)
! protein - protein, protein - water and water - water interactions
! all use same type of nblist -> can use same code
! arguments
TYPE(NB_TYPE)                   :: nb
! locals
integer                         :: i, j
real(kind=prec)                 :: r2,r,r6,r12
TYPE(qr_dist2)                  :: distance

! global variables used:
! x

i = nb%i
j = nb%j

distance = q_dist2(x(i),x(j))
! q_dist2 already gives distances as one/r*
! so no need for further conversion
r2    = distance%r2
r     = distance%r
r6    = distance%r6
r12   = distance%r12

! calculate Vel and dv
nbe%Vel  = nb%elec * r
nbe%V_a  = nb%vdWA *r12 
nbe%V_b  = nb%vdWB *r6
nbe%dv   = r2*( -nbe%Vel -12.0_prec*nbe%V_a +6.0_prec*nbe%V_b )
nbe%vec  = distance%vec
end function nbe

!------------------------------------------------------------------------
TYPE(ENERET_TYPE) function nbe_b(nb,shift)
! same as above, adopted for periodic box shift
! arguments
TYPE(qr_vec)                    :: shift
TYPE(NB_TYPE)                   :: nb
! local variables
integer                         :: i,j
real(kind=prec)                 :: r2,r,r6,r12
TYPE(qr_dist2)                  :: distance

! global variables used:
! x

i  = nb%i
j  = nb%j

! calculate dx, r and r2
distance = q_dist2((x(i)-x(j)),shift)


r2   = distance%r2
r    = distance%r
r6   = distance%r6
r12  = distance%r12

nbe_b%Vel  = nb%elec * r
nbe_b%V_a  = nb%vdWA *r12 
nbe_b%V_b  = nb%vdWB *r6
nbe_b%dv   = r2*( -nbe_b%Vel -12.0_prec*nbe_b%V_a +6.0_prec*nbe_b%V_b )
nbe_b%vec  = distance%vec
end function nbe_b
!----------------------------------------------------------------------

TYPE(ENERET_TYPE) function nbe_qq(nb,lambda)
! returns the nonbonded energies for qq interactions
! arguments
TYPE(NBQ_TYPE)                          :: nb
real(kind=prec)                         :: lambda
! local variables
integer                                 :: iq,i,jq,j
real(kind=prec)                         :: r2,r,r6,r12,r6_hc
TYPE(qr_dist3)                          :: distance

  ! for every pair:
iq   = nb%iq
i    = iqseq(iq)
jq   = nb%jq
j    = iqseq(jq)

! calculate dx and r
distance = q_dist3(x(i),x(j))

r2    = distance%r2
r6_hc = one/distance%r6
r     = distance%r
r6    = r6_hc + nb%score !softcore
r6    = one/r6
r12   = r6*r6

! calculate Vel, V_a, V_b and dv
nbe_qq%Vel  = nb%elec *r
nbe_qq%vec  = distance%vec
if (nb%soft) then
nbe_qq%V_b = zero
nbe_qq%V_a = nb%vdWA*exp(-nb%vdWB/r)
nbe_qq%dv  = r2*( -nbe_qq%Vel -nbe_qq%V_b*nbe_qq%V_a/r )*lambda
else
nbe_qq%V_a  = nb%vdWA *r12 
nbe_qq%V_b  = nb%vdWB *r6
nbe_qq%dv  = r2*( -nbe_qq%Vel -(12.0_prec*nbe_qq%V_a -6.0_prec*nbe_qq%V_b)*r6*r6_hc )*lambda
endif

end function nbe_qq

!-----------------------------------------------------------------------

TYPE(ENERET_TYPE) function nbe_qqb(nb,lambda)
! returns the nonbonded energies for qq interactions with box (only in nb monitor)
! arguments
TYPE(NBQ_TYPE)                          :: nb
real(kind=prec)                         :: lambda
! local variables
integer                                 :: iq,i,jq,j
real(kind=prec)                         :: r2,r,r6,r12,r6_hc
TYPE(qr_dist3)                          :: distance
TYPE(qr_vec)                            :: shift

  ! for every pair:
iq   = nb%iq
i    = iqseq(iq)
jq   = nb%jq
j    = iqseq(jq)

! calculate dx and r
shift = x(i) - x(j)
shift = boxlength*q_nint(shift*inv_boxl)
distance = q_dist3(x(i)-x(j),shift)

r2    = distance%r2
r6_hc = one/distance%r6
r     = distance%r
r6    = r6_hc + nb%score !softcore
r6    = one/r6
r12   = r6*r6

! calculate Vel, V_a, V_b and dv
nbe_qqb%Vel  = nb%elec *r
nbe_qqb%vec  = distance%vec
if (nb%soft) then
nbe_qqb%V_b = zero
nbe_qqb%V_a = nb%vdWA*exp(-nb%vdWB/r)
nbe_qqb%dv  = r2*( -nbe_qqb%Vel -nbe_qqb%V_b*nbe_qqb%V_a/r )*lambda
else
nbe_qqb%V_a  = nb%vdWA *r12 
nbe_qqb%V_b  = nb%vdWB *r6
nbe_qqb%dv  = r2*( -nbe_qqb%Vel -(12.0_prec*nbe_qqb%V_a -6.0_prec*nbe_qqb%V_b)*r6*r6_hc )*lambda
endif

end function nbe_qqb

!-----------------------------------------------------------------------

TYPE(ENERET_TYPE) function nbe_qx(nb,lambda,dist)
! arguments
TYPE(NBQP_TYPE)                 :: nb
TYPE(qr_dist3)                  :: dist
real(kind=prec)                 :: lambda
! local variables
integer                         :: ip,iq,i,j
real(kind=prec)                 :: r2,r,r6,r6_hc,r12


r2    = dist%r2
r6_hc = one/dist%r6
r     = dist%r
r6    = r6_hc + nb%score !softcore
r6    = one/r6
r12   = r6*r6

! calculate qi, Vel, V_a, V_b and dv
nbe_qx%Vel  = nb%elec *r
nbe_qx%V_a  = nb%vdWA *r12
nbe_qx%V_b  = nb%vdWB *r6
nbe_qx%dv   = r2*( -nbe_qx%Vel -(12.0_prec*nbe_qx%V_a -6.0_prec*nbe_qx%V_b)*r6*r6_hc )*lambda

end function nbe_qx

!-----------------------------------------------------------------------

TYPE(ENERET_TYPE) function nbe_qspc(nb,lambda,dist)
! same function as for normal qx return, but does not calculate vdW
! for use with spc water -> might use precomp instruction
! to use same function?
! arguments
TYPE(NBQP_TYPE)                         :: nb
real(kind=prec)                         :: lambda
TYPE(qr_dist)                           :: dist
! locals
real(kind=prec)                         :: r, r2
integer                                 :: iq, i, j

r2   = dist%r2
r    = dist%r

nbe_qspc%Vel = nb%elec *r
nbe_qspc%dv  = -r2*nbe_qspc%Vel*lambda

end function nbe_qspc

!-----------------------------------------------------------------------

TYPE(ENERET_TYPE) function nbe_spc(nb)
! modified nbe_p function to only return electrostatic forces for spc water
! arguments
TYPE(NB_TYPE)                   :: nb
! locals
integer                         :: i, j
real(kind=prec)                 :: r2, r
TYPE(qr_dist)                   :: distance

i    = nb%i
j    = nb%j

distance = q_dist(x(i),x(j))
r2 = distance%r2
r  = distance%r

nbe_spc%vec = distance%vec
nbe_spc%Vel = nb%elec *r
nbe_spc%dv  = r2*( -nbe_spc%Vel)

end function nbe_spc

!-----------------------------------------------------------------------
!******PWadded 2001-10-23
TYPE(ENERET_TYPE) function nbe_spcb(nb,shift)
! same as above, with added correction for periodic box shifts
! arguments
TYPE(NB_TYPE)                   :: nb
TYPE(qr_vec)                    :: shift
! locals
integer                         :: i, j
real(kind=prec)                 :: r2, r
TYPE(qr_dist)                   :: distance

i    = nb%i
j    = nb%j

distance = q_dist(x(i)-x(j),shift)
r2 = distance%r2
r  = distance%r
nbe_spcb%vec = distance%vec
nbe_spcb%Vel = nb%elec *r
nbe_spcb%dv  = r2*( -nbe_spcb%Vel)

end function nbe_spcb

!-----------------------------------------------------------------------

end module NONBONDED




