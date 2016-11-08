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
use POTENE
use GLOBALS
use QMATH
use MPIGLOB
! this module does the actual path integral calculations and metropolis sampling of the bead configurations
! code is directly based on the paper, with more bugs than I can count

implicit none
! global QCP values
real(kind=prec) :: beta, wl_lam, tmp_wl_lam, cboltz, hbar, convert, pi_fac
! array of wavelengths for every atom 
real(kind=prec),allocatable :: wl_lam0(:), wl_lam1(:), wl_lam2(:)
! array of square root of mass ratios for every atom, for mass perturbed QCP (I guess)
! only used if user has set masses with qcp_mass
real(kind=prec),allocatable :: sqmass(:)
! energy storage for beads and particles
! one Q energy array for every FEP lambda?
! added in the end to global EQ as last entry in ene_header%array
TYPE(OQ_ENERGIES),allocatable   :: EQtmp(:),EQsave(:),EQpot_ave(:),qcp_EQbead_old(:,:),qcp_EQbead_new(:,:)
real(kind=prec),allocatable     :: qcp_Ebeta(:),qcp_EQbeta(:,:),total_EQPI(:)
TYPE(ENERGIES)                          :: qcp_E_tot
TYPE(OQ_ENERGIES),allocatable           :: qcp_EQ_tot(:)
real(kind=prec)                         :: total_EPI,Epot_ave,qcp_Ering_ave,Esave,Etmp,qcp_Ebead
real(kind=prec)                         :: qcp_Ering_new,qcp_Ering_old,qcp_Ering
real(kind=prec),allocatable             :: qcp_Ebead_old(:),qcp_Ebead_new(:)
! for gaussian deviation algorithm
integer                         :: gauss_set = 0
real(kind=prec)                 :: gval = zero
! coordinates for ring
TYPE(qr_vec),allocatable        :: qcp_coord(:,:),qcp_coord_old(:,:)

contains

subroutine qcp_init
! set up constants after reading in information from inputs
! no matter if input is in MD or Postprocessing
! arguments
! TODO set constants in MISC
integer :: i

convert    = 4.184_prec / 6.022E23_prec
cboltz     = 1.38064852E-23_prec
hbar       = 1.054571800E-34_prec

!all temp dependent variables set with current instantanious temperature
!during actual QCP calc

! TODO need allocation check
allocate(wl_lam0(qcp_atnum),wl_lam1(qcp_atnum),wl_lam2(qcp_atnum))
! one array each for fep states and qcp steps to save energy
allocate(qcp_EQ(nstates))
! allocate bead coordinate array
allocate(qcp_coord(qcp_atnum,qcp_size),qcp_coord_old(qcp_atnum,qcp_size))

allocate(qcp_Ebeta(maxval(qcp_steps)))
allocate(qcp_EQbeta(maxval(qcp_steps),nstates),EQsave(nstates),total_EQPI(nstates),&
        EQpot_ave(nstates),EQtmp(nstates))

allocate(qcp_EQbead_old(qcp_size,nstates),qcp_EQbead_new(qcp_size,nstates))

allocate(qcp_EQ_tot(nstates))

allocate(x_save(natom),d_save(natom))

allocate(qcp_Ebead_old(qcp_size),qcp_Ebead_new(qcp_size))

! will need this later
qcp_EQ(1:nstates)%lambda = EQ(1:nstates)%lambda

if (use_qcp_mass) then
	allocate(sqmass(qcp_atnum))
end if

! set wl values for every atom
do i = 1, qcp_atnum
call qcp_init_beads(qcp_atom(i),qcp_coord(i,:))
end do

! what are TIAC, CHIN, CHAC????
end subroutine qcp_init

subroutine qcp_init_beads(atom,qcpvec)
! create initial guess for bead positions, wont even try to read stuff in
! again stolen from original code by D.T. Major
! arguments
integer                         :: atom
TYPE(qr_vec)                    :: qcpvec(:)
! locals
! displ -> displacment, dependent on heavy atom or not
real(kind=prec)                 :: angle,percent,displ_l,displ_h,displ,dangl
integer                         :: i, irand,jrand,krand

qcpvec = qcpvec * zero

end subroutine qcp_init_beads

function qcp_bisect(qcpvec,wavelength,level)
! this one takes care of the actual bisection step to generate the bead positions
! shamelessly stolen from original code by D.T. Major and J. Gao
! based on following reference for the crazy
! Levy, P., Compostio Math. 7, 283, 1939
! Ceperly, F.M.; Pollock, E.L. Monte Carlo methods in theoretical physics, p 35. ETS Editrice, Pisa 1992

! arguments
! all the coordinates for the beads of a given PI atom
TYPE(qr_vec)                    :: qcpvec(:)
TYPE(qr_vec)                    :: qcp_bisect(size(qcpvec))
! wavelength, used for step
real(kind=prec)                 :: wavelength
! bisection level, I guess
integer                         :: level

! locals
real(kind=prec),parameter       :: half=0.5_prec
integer                         :: n,i,j,nbsect,bbsect,nextbead,medbead,thisbead
TYPE(qr_vec)                    :: midpoint
real(kind=prec)                 :: lbsect

! for a given particle, go through all beads
! perform number of steps depending on bisection level
! starting at highest and working down
        do i = level, 1, -1
        lbsect   = q_sqrt(real(2**(i-1),kind=prec))
! set number of bisection steps
        nbsect   = 2**(level-i)
! current bisection size
        bbsect   = 2**i
        thisbead = 1
        nextbead = thisbead + bbsect
        medbead  = int(thisbead+nextbead)/2
        do j = 1, nbsect
        if (thisbead.gt.qcp_size) thisbead = thisbead - qcp_size
        if (nextbead.gt.qcp_size) nextbead = nextbead - qcp_size
        if (medbead .gt.qcp_size) medbead  = medbead  - qcp_size
        midpoint = (qcpvec(thisbead) + qcpvec(nextbead)) * half
        qcp_bisect(medbead) = qcp_gauss(wavelength*lbsect,midpoint)
        thisbead = thisbead + bbsect
        nextbead = nextbead + bbsect
        medbead  = medbead  + bbsect
        end do ! nbsect
        end do ! level

end function qcp_bisect

function qcp_center(qcpvec)
! centers all beads of a qcp atom on the middle
! arguments
TYPE(qr_vec)                    :: qcpvec(:)
TYPE(qr_vec)                    :: qcp_center(size(qcpvec))
! locals
integer                         :: i
real(kind=prec)                 :: inv_beads
TYPE(qr_vec)                    :: center

inv_beads = one/(real(qcp_size ,kind=prec))
center = center * zero
do i = 1, qcp_size
center = center + qcpvec(i)
end do
center = center * inv_beads
do i = 1, qcp_size
qcp_center(i) = qcpvec(i) - center
end do

end function qcp_center


TYPE(qr_vec) function qcp_gauss(sigma,coord)
! returns the gaussian distribution weighted coordinate
! reimplemented from code provided by D.T. Major
! from Numerical Recipies, Ch. 7.2, according to Box-Muller algorithm
! separate from general Q gaussian algorith to not be influenced by this
! arguments
TYPE(qr_vec)                    :: coord
real(kind=prec)                 :: sigma
! locals
real(kind=prec)                 :: val1,val2,fac,gdev,rval
real(kind=prec),parameter       :: two=2.0_prec
integer                         :: i

do i = 1, 3

if (gauss_set.eq.0) then

        do 
        val1 = 2.0_prec * qcp_randm() - one
        val2 = 2.0_prec * qcp_randm() - one
        rval = val1**2 + val2**2
        if (rval.gt.zero .and. rval .lt.one) exit
        end do

        fac  = q_sqrt(-two*q_logarithm(rval)/rval)
        gval = val1 * fac
        gdev = val2*fac*sigma
        gauss_set = 1
else
        gdev = gval * sigma
        gauss_set = 0
end if
if (i.eq.1) qcp_gauss%x = coord%x + gdev
if (i.eq.2) qcp_gauss%y = coord%y + gdev
if (i.eq.3) qcp_gauss%z = coord%z + gdev
end do

end function qcp_gauss

subroutine qcp_run(temp,E_init,EQ_init)
! main meaty subroutine that performs the actual run to generate coordinates
! and evaluate energies 
! based on original code by D.T. Major and J.Gao
! after ripping out all the stuff I don't need
! only for bisection now, have to implement rest later
! arguments
real(kind=prec)                         :: temp
TYPE(ENERGIES), INTENT(in)              :: E_init
TYPE(OQ_ENERGIES),INTENT(in)            :: EQ_init(:)
! locals
integer                                 :: ncent,nmove,qcp_loop,istate
integer                                 :: i,ia,j,step
real(kind=prec)                         :: bfactor,deltaH,irand,MCaccept,MCreject,MCrepeat
! in case of MPI this is a mess, because only the pot_ene stuff needs it
! everything else needs to be on the master node to keep it clean

! first set of energies will come from independent calculation, in MD this is ther normal EQ,
! in case of post processing this will be EQ done on the first restart
! first thing needed is saving coordinates and forces to restore them later

x_save = x
d_save = d

! set temperature dependent variables
if (nodeid .eq. 0) then
beta       = convert / (temp * cboltz)
tmp_wl_lam = angstrom * cboltz * temp / hbar
pi_fac     = (0.5_prec * real(qcp_size, kind=prec) * tmp_wl_lam**2 * amu) / convert
tmp_wl_lam = hbar / (2.0_prec * tmp_wl_lam * real(qcp_size,kind=prec) * amu ) 

! set wl values for every atom
do i = 1, qcp_atnum
	wl_lam0(i) = tmp_wl_lam * winv(iqseq(qcp_atom(i)))
	if(use_qcp_mass) sqmass = sqrt(one/winv(iqseq(qcp_atom(i)))/QCP_mass(i))
	wl_lam1(i) = sqrt(wl_lam0(i))
	wl_lam2(i) = 2.0_prec * wl_lam1(i)
end do

end if

do i = 1, qcp_atnum
call qcp_init_beads(qcp_atom(i),qcp_coord(i,:))
end do


qcp_EQ_tot = qcp_EQ_tot * zero
qcp_Ebead_old = zero
qcp_Ebead_new = zero
qcp_Ebead     = zero
Etmp  = zero
Esave = zero
EQtmp = EQtmp * zero
do i = 1, maxval(qcp_steps)
qcp_EQbead_old(i,:) = qcp_EQbead_old(i,:) * zero
qcp_EQbead_new(i,:) = qcp_EQbead_new(i,:) * zero
end do
if (nodeid.eq.0) then
! set number of beads to move, 1 bead for bisection, center all
ncent = qcp_size
nmove = 1
end if

! main loop, first equil, then sampling
! set steps according to qcp_steps
do qcp_loop = 1, 2

if (qcp_loop .eq. 2) then
! sampling needs to set some stuff first

! move all beads, one after the other
! first, each bead for each atom
! only on node 0 for mpi
do i = 1, qcp_size
if (nodeid.eq.0) then
do j = 1, qcp_atnum
ia = qcp_atom(j)
x(ia) = x_save(ia) + qcp_coord(j,i)
end do ! qcp_atnum
end if

! MPI BCast new coordinates to nodes
#ifdef USE_MPI
call MPI_BCast(x,nat3,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
if (ierr .ne. 0) call die('QCP Bcast x')
#endif

! all atoms at first perturbed position -> get energy
! E and EQ and set to zero in pot_ene, don't need to do this manually
call pot_energy(qcp_E,qcp_EQ,.false.)
! rest set to zero
if (nodeid.eq.0) then
Etmp = Etmp + qcp_E%potential
EQtmp(:) = EQtmp(:) + qcp_EQ(:)
qcp_Ebead_old(i) = qcp_E%potential
do  istate = 1, nstates
qcp_EQbead_old(i,istate) = qcp_EQ(istate)
end do
end if
end do ! qcp_size

! what is factor? CHFACN is one for not using strange CHIN method stuff
if (nodeid.eq.0) then
        Esave = Etmp / real(qcp_size,kind=prec)
        Etmp = Esave
end if

if (nodeid.eq.0) then
        EQsave(:) = EQtmp(:) / real(qcp_size,kind=prec)
        EQtmp = EQsave
end if
end if ! qcp_loop .eq.2

if (nodeid .eq.0) then
qcp_Ering = zero
do i = 1, qcp_atnum
qcp_Ering = qcp_Ering + qcp_energy(qcp_atom(i),qcp_coord(i,:),pi_fac)
end do

qcp_Ering_old = qcp_Ering

! set some average values
! why another one of those
total_EPI = zero
qcp_Ering_ave = zero
Epot_ave  = zero
EQpot_ave = EQpot_ave * zero
total_EQPI = zero
! MC stuff
MCrepeat = one
MCaccept = zero
MCreject = zero

end if

do step = 1 , qcp_steps(qcp_loop)
! save bead coords

if (nodeid.eq.0) then

qcp_coord_old = qcp_coord
! call actual bisection step algortihm
! and recenter ring
! for each ring
qcp_Ering = zero
do i = 1, qcp_atnum
qcp_coord(i,:) = qcp_bisect(qcp_coord(i,:),wl_lam1(i),qcp_level)
qcp_coord(i,:) = qcp_center(qcp_coord(i,:))
qcp_Ering = qcp_Ering + qcp_energy(qcp_atom(i),qcp_coord(i,:),pi_fac)
end do
qcp_Ering_new = qcp_Ering

end if

! needed for selection, different for other algorithms
! those would also need Bcasting of energies then
deltaH = zero

bfactor = q_exp(-deltaH*beta)

! get some random number?
irand = qcp_randm()

! this is the monte carlo part
if (bfactor .ge. irand) then
! energy only calculated for sampling, not for equilibration? weird!
        if(qcp_loop .eq. 2) then
                Esave  = Etmp
                EQsave = EQtmp
                if(qcp_veryverbose) then
                write(*,'(a,i4,a,f12.6,a,f12.6)') 'Values at step ',step,' for Esave = ',Esave,' , and Etmp = ',Etmp
                        do istate = 1 ,nstates
                        write(*,'(a,i4,a,i4)') 'Values at step ',step,' state = ',istate
                        write(*,'(a,f12.6,a,f12.6)') 'EQsave = ',EQsave(istate)%total,&
                                ' , EQtmp = ',EQtmp(istate)%total
                        end do
                end if
                do i = 1, qcp_size
                if (nodeid.eq.0) then
                        do j = 1, qcp_atnum
                        ! only atoms with qcp_coords need update, saved in x_save
                        x(qcp_atom(j)) = x_save(qcp_atom(j)) + qcp_coord(j,i)
                        end do ! qcp_atnum
                end if
#ifdef USE_MPI
call MPI_Bcast(x,nat3,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
if (ierr .ne. 0) call die('QCP Bcast x')
#endif

                call pot_energy(qcp_E,qcp_EQ,.false.)
                if (nodeid.eq.0) then
                        qcp_Ebead_new(i) = qcp_E%potential
                        qcp_EQbead_new(i,:) = qcp_EQ(:)
                        Esave = Esave + ((qcp_Ebead_new(i)-qcp_Ebead_old(i))/real(qcp_size,kind=prec))
                        EQsave(:) = EQsave(:) + ((qcp_EQbead_new(i,:) - qcp_EQbead_old(i,:)) &
                                /real(qcp_size,kind=prec))
                        qcp_Ebead_old(i) = qcp_Ebead_new(i)
                        qcp_EQbead_old(i,:) = qcp_EQbead_new(i,:)
                end if
                end do ! qcp_size


                ! sum up some stuff
                if (nodeid .eq.0) then
                        qcp_Ebeta(step) = q_exp(-beta*(Etmp-E_init%potential))
                        if(qcp_veryverbose) then
                                write(*,'(a,i4,a,f12.6)') 'Current beta for step ',step,' = ',qcp_Ebeta(step)
                                write(*,'(a,f12.6)') 'Etmp = ',Etmp
                                write(*,'(a,f12.6)') 'E_init = ',E_init%potential
                        end if
                        do istate = 1, nstates
                        qcp_EQbeta(step,istate) = q_exp(-beta*(EQtmp(istate)%total-EQ_init(istate)%total))
                        if(qcp_veryverbose) then
                                write(*,'(a,i4)') 'Current step = ',step
                                write(*,'(a,i4)') 'Current state = ',istate
                                write(*,'(a,f12.6)') 'Current beta = ',qcp_EQbeta(step,istate)
                                write(*,'(a,f12.6)') 'EQtmp = ',EQtmp(istate)%total
                                write(*,'(a,f12.6)') 'EQ_init = ',EQ_init(istate)%total
                        end if

                        end do
                        total_EPI = total_EPI + MCrepeat * qcp_Ebeta(step)
                        total_EQPI(:) = total_EQPI(:) + MCrepeat * qcp_EQbeta(step,:)
                        Epot_ave  = Epot_ave  + MCrepeat * Etmp
                        EQpot_ave(:) = EQpot_ave(:) + (EQtmp(:) * MCrepeat)
                        Etmp = Esave
                        EQtmp = EQsave
                end if
        end if ! qcp_loop .eq. 2
        if (nodeid .eq.0) then
                MCaccept = MCaccept + MCrepeat
                qcp_Ering_ave = qcp_Ering_ave + MCrepeat * qcp_Ering_new
                qcp_Ering_old = qcp_Ering_new 
                MCrepeat = one
        end if
else ! failed to accept step
        if (nodeid.eq.0) then
                qcp_coord = qcp_coord_old
                MCreject = MCreject + one
                MCrepeat = MCrepeat + one
        end if
end if ! bfactor .ge. irand

! write info about qcp steps if desired -> flag is qcp_veryverbose
if (nodeid.eq.0 .and. qcp_veryverbose .and. qcp_loop .eq. 2) then
        write(*,'(a)') 'Information for individual step of QCP PI'
        write(*,'(a)') 'Current status is free particle sampling'
        write(*,'(a,i4)') 'Current conformation = ',step
        write(*,'(a,f12.6)') 'Current total particle energy = ',Esave
        write(*,'(a,10f12.6)') 'Current Q particle energy, for each state = ',EQsave(:)%total
        write(*,'(a,f12.6)') 'Current ring polymer energy = ',qcp_Ering_new
end if

end do ! step

if (nodeid.eq.0) then
! write out info about QCP success rate and energies
        if (qcp_loop .eq.1) then
                write(*,'(a)') 'QCP summary for free particle equilibration'
        else
                write(*,'(a)') 'QCP summary for free particle sampling'
        end if
        write(*,'(a,i4)') 'Total number of steps is ',step-1
        write(*,'(a,i4)') 'Number of accepted MC steps = ',int(MCaccept)
        write(*,'(a,i4)') 'Number of rejected MC steps = ',int(MCreject)
end if

! make averages for saving to file
if (nodeid.eq.0 )then
        if (MCaccept .eq. zero) then
                MCaccept = real(qcp_steps(qcp_loop),kind=prec)
                qcp_Ering_ave = qcp_Ering_ave/MCaccept
        end if
end if
end do ! qcp_loop

if (nodeid .eq.0) then

total_EPI = total_EPI / MCaccept
Epot_ave  = Epot_ave  / MCaccept
total_EQPI(:) = total_EQPI(:) / MCaccept
EQpot_ave(:) = EQpot_ave(:) / MCaccept

if(qcp_verbose) then
        write(*,'(a,f12.6)') 'Total PI energy = ',total_EPI
        write(*,'(a,f12.6)') 'Final PI potential energy = ',Epot_ave
        do i = 1 , nstates
        write(*,'(a,i4)') 'Total QPI energy for state ',i
        write(*,'(f12.6)') total_EQPI(i)
        end do
        do i = 1 , nstates
        write(*,'(a,i4)') 'istate = ',i
        write(*,'(a,f12.6)') 'Final PI Q potential energy = ',EQpot_ave(i)%total
        end do
end if

! now we have all the energies for the PI sampling, can print info for user and write to global EQ_save 
! qcp_pos is set during loading, need this later for mass perturbation with two arrays
call qatom_savetowrite(EQpot_ave,qcp_pos)

! set everything back to how we found it 
x = x_save
d = d_save

end if
! Bcast when using MPI
#ifdef USE_MPI
call MPI_Bcast(x,nat3,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
if (ierr .ne. 0) call die('QCP Bcast x')
call MPI_Bcast(d,nat3,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
if (ierr .ne. 0) call die('QCP Bcast d')
#endif

! and we are done 

end subroutine qcp_run


real(kind=prec) function qcp_energy(atom,coord,wavelength)
! calculates the ring energies of the particles
! call this for each ring
! arguments
TYPE(qr_vec)                    :: coord(:)
integer                         :: atom
real(kind=prec)                 :: wavelength
! locals
real(kind=prec)                 :: fk
integer                         :: i

qcp_energy = zero
! will need adapting for different QCP types
fk = wavelength/winv(atom)

do i = 1, qcp_size - 1
qcp_energy = qcp_energy + (fk * q_dist4(coord(i),coord(i+1)))
end do

! this only later for closed chain chain polymer
qcp_energy = qcp_energy + (fk * q_dist4(coord(1),coord(qcp_size)))


end function qcp_energy


real(kind=prec) function qcp_randm()
! same random number generator as in charmm code provided by D.T. Major
!
!-----------------------------------------------------------------------
!
!     RANDOM NUMBER GENERATOR
!     SEPTEMBER, 1986: WLJ & J. M. Briggs
!
!     NOTES ON USAGE:
!     THE STATEMENT "COMMON/RANDS/IRN" MUST APPEAR IN EACH
!     PROGRAM,SUBROUTINE... WHERE THE SEED VALUE (IRN) IS
!     BEING WRITTEN OUT OR BEING READ IN. THIS FUNCTION IS
!     USED IN THE SAME MANNER AS THAT USED ON THE GOULD,
!     (I.E. "RANDN = RANU()" ). IRN SHOULD INITIALLY BE
!     AN ODD, I6 INTEGER BUT CAN GET AS BIG AS I7. JMB
!     IMOD-1 UNIQUE RANDOM NUMBERS ARE GENERATED . A SINGLE
!     PRECISION LINEAR CONGRUENTIAL GENERATOR IS USED WITH
!     SHUFFLING OF THE CONSTANT AND OF WHEN THE SHUFFLING IS
!     DONE. CONSEQUENTLY, THE PERIOD IS EXTREMELY LONG - NONE
!     WAS FOUND IN TESTS GENERATING MILLIONS OF RANDOM NUMBERS.
!-----------------------------------------------------------------------
real(kind=prec) :: RNJ,FAC
integer,save :: ICNT=0   ,ICHG=1167     ,ICHG0=1167  ,ICN0=458753759,  &
                IMUL=1173,ICON=458753759,IMOD=1048573, &
                JMUL=1161,JCON=458716759,JMOD=1048573,JRN=124690
integer         :: IRN

  ICNT = ICNT+1
  IF (ICNT.EQ.ICHG) THEN
     !
     !-----------------------------------------------------------------------
     !
     !     CHANGE ICON USING SECONDARY GENERATOR
     !
     !-----------------------------------------------------------------------
     !
     JRN = JRN*JMUL+JCON
     JRN = MOD(JRN,JMOD)
     RNJ = real(JRN,kind=prec)/real(JMOD,kind=prec)
     IF (RNJ.GT.0.5_prec) THEN
        FAC = one + 0.5_prec*RNJ
     ELSE
        FAC = one - 0.5_prec*RNJ
     ENDIF
     FAC = real(ICN0,kind=prec)*FAC
     ICON = INT(FAC)
     !
     !-----------------------------------------------------------------------
     !
     !     CHANGE ICHG USING SECONDARY GENERATOR
     !
     !-----------------------------------------------------------------------
     !
     JRN = JRN*JMUL+JCON
     JRN = MOD(JRN,JMOD)
     RNJ = real(JRN,kind=prec)/real(JMOD,kind=prec)
     IF (RNJ.GT.0.5_prec) THEN
        FAC = one + 0.5_prec*RNJ
     ELSE
        FAC = one - 0.5_prec*RNJ
     ENDIF
     FAC  = real(ICHG0,kind=prec)*FAC
     ICHG = INT(FAC)
     ICNT = 0
  ENDIF
  !
  !-----------------------------------------------------------------------
  !
  !     GENERATE RANDOM NUMBER
  !
  !-----------------------------------------------------------------------
  !
  IRN = IRN*IMUL+ICON
  IRN = MOD(IRN,IMOD)
  qcp_randm = real(IRN,kind=prec)/real(IMOD,kind=prec)

end function qcp_randm

end module QCP
