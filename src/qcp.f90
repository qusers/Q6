! (C) 2014 Uppsala Molekylmekaniska HB, Uppsala, Sweden
! qcp.f90
! Q is written by Johan Åqvist, John Marelius, Anders Kaplan, Isabella Feierberg, Martin Nervall & Martin Almlöf
! Quantum classical path module
! written by Paul Bauer
! based on bisection sampling QCP developed by Dan Thomas Major and Jiali Gao
module QCP

use TOPO
use QATOM
use SIZES
use NRGY
use NONBONDED
use BONDED

! this module does the actual path integral calculations and metropolis sampling of the bead configurations
! code is directly based on the paper, with more bugs than I can count

! global QCP values
! most important: De Brogie Wavelenght, controls step size
real(kind=prec) :: DeB_WL
! other constants, beta, wavelength lambda (wl_lam) and its temp value, boltzmann constant (cboltz), reduced planck constant (hbar)
! conversion of joule <-> calorie (convert)
! path integral factor (FPIFAC???)
real(kind=prec) :: beta, wl_lam, tmp_wl_lam, cboltz, hbar, convert, pi_fac
! array of wavelengths for every atom 
real(kind=prec),allocatable :: wl_lam0(:), wl_lam1(:), wl_lam2(:)
! array of square root of mass ratios for every atom, nned for mass perturbed QCP
! only used if user has set masses with qcp_mass
real(kind=prec),allocatable :: sqmass(:)
contains

subroutine qcp_init
! set up constants after reading in information from inputs
! no matter if input is in MD or Postprocessing
! TODO set constants in MISC
integer :: i

convert    = one / 4.184_prec
cboltz     = 1.38064852E-23_prec
cboltz     = cboltz * convert
beta       = one / (temp * cboltz)
hbar       = 1.054571800E-34_prec
hbar       = hbar * convert
tmp_wl_lam = cboltz * temp / hbar
pi_fac     = 0.5_prec * real(num_Beads, kind=prec) * tmp_wl_lam**2
tmp_wl_lam = hbar / (2.0_prec * tmp_wl_lam * real(num_Beads,kind=prec))

! TODO need allocation check
allocate(wl_lam0(num_QCP),wl_lam1(num_QCP),wl_lam2(num_QCP))

if (use_qcp_mass) then
	allocate(sqmass(num_QCP))
end if

! set wl values for every atom
do i = 1, num_QCP
	wl_lam0(i) = tmp_wl_lam / mass(iqseq(QCP_atoms(i)))
	if(use_qcp_mass) sqmass = sqrt(mass(iqseq(QCP_atoms(i)))/QCP_mass(i))
	wl_lam1(i) = sqrt(wl_lam0(i)/real(qcp_mc_step,kind=prec))
	wl_lam2(i) = 2.0_prec * wl_lam1(i)
end do
! what are TIAC, CHIN, CHAC????
end subroutine qcp_init

subroutine qcp_mc(qcpatom,wavelength1,wavelength2,coords)
! this one does the monte carlo magic to propagate the RP and allow the sampling
! needs its onw random number generator implemented below and shamelessly stolen from the original code
! by D.T. Major and J. Gao
! arguments
integer                                 :: qcpatom
real(kind=prec)                         :: wavelength1,wavelength2
TYPE(qr_vec)                            :: coords(:)
! locals



! QCP does not need to calculate anything related to integrator
! or step propagation, only brute force energy
! pass the array of atoms to MD_NONBONDED energy calculation
! and be happy :)
subroutine qcp_energy



end subroutine qcp_energy

end module QCP
