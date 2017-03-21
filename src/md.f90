! (C) 2014 Uppsala Molekylmekaniska HB, Uppsala, Sweden
! md.f90
! by Johan Åqvist, John Marelius, Anders Kaplan, Isabella Feierberg, Martin Nervall & Martin Almlöf
! molecular dynamics
!TODO: remove default real statment and change real(4) - in accordance with best practice

module MD

! used modules
!use PROFILING
use SIZES
use TRJ
use MPIGLOB
use QATOM
use EXC
use VERSIONS
use QMATH
use GLOBALS
use SIMPREP
use POTENE
use QCP
#if defined (_DF_VERSION_)
use DFPORT
use DFLIB
#endif
!$ use omp_lib
implicit none



!----START OF PUBLIC SUBROUTINES
contains


!----------------------------------------------------------------------


subroutine md_startup
! initialise used modules
call simprep_startup

! read in version info
MD_VERSION = trim(version_pass())
MD_DATE    = trim(date_pass())


end subroutine md_startup


!----------------------------------------------------------------------

subroutine md_shutdown
! call used modules' shutdown subroutines
if(use_qcp) call qcp_shutdown
call simprep_shutdown
call md_deallocate
call topo_deallocate
end subroutine md_shutdown

!-----------------------------------------------------------------------

#ifdef USE_GRID
#ifdef USE_MPI
integer function grid_add(grid1,grid2,length,MPI_Datatype)
! arguments
integer,INTENT (IN)      :: grid1(length)
integer                  :: grid2(length)
integer                         :: length,MPI_Datatype
! locals
integer                 :: j,k

! we have to do something nasty to be able to combine the grid information form
! all the nodes, so here it is
! we search first for a zero element in the final grid, so we not overwrite any
! info that is already there
! as soons as we find the first zero element, we search the slave grid in the
! same row until we find a zero element, and copy all the groups from before
! into the master array
! quick and very dirty way
j = 1
do while (grid2(j) .ne. 0 )
        j = j + 1
end do
k = 1
do while (grid1(k) .ne. 0 )
        grid2(j+k) = grid1(k)
end do

grid_add = -1

end function grid_add
#endif
#endif

!-----------------------------------------------------------------------
logical function initialize()                  
! local variables
character					:: text*200
integer						:: i,j,length,ii
real(kind=prec)						:: stepsize
real(kind=prec)						:: lamda_tmp(max_states)
integer						:: fu, fstat
real(kind=prec)						::	rjunk
integer						::	ijunk

! local parameters
integer						:: num_args
character(200)				:: infilename
logical						::	yes
logical						::	need_restart
character(len=200)				::	instring
logical						::	inlog
integer						::	mask_rows, number
integer,parameter				:: maxmaskrows=30
character(len=200)				:: gc_mask_tmp(maxmaskrows),str,str2
logical                                         :: shake_all,shake_all_solvent,shake_all_solute
logical                                         :: shake_all_hydrogens,shake_all_heavy
character(200)                                  :: qcp_select,qcp_size_name
integer                                         :: timeval(8)
integer                                         :: tgroups
! this subroutine will init:
!  nsteps, stepsize, dt
!  Temp0, tau_T, iseed, Tmaxw
!  use_LRF, NBcycle, Rcpp, Rcww, Rcpw, Rcq
!  shake_solute, shake_solvent, shake_hydrogens, shake_heavy
! fk_pshell
!  fk_wsphere=-1, wpol_restr, wpol_born
!  fkwpol=-1, Dwmz=-1 (values  ized to -1 will
!    be set in water_sphere, once target radius is known)
!  top_file
!  restart, [restart_file]
!  xfin_file
!  itrj_cycle, iene_cycle, iout_cycle, itemp_cycle, [trj_file], [ene_file]
!  fep_file
!  nstates, EQ (allocating memory for EQ)
!  implicit_rstr_from_file, [exrstr_file]
!  nrstr_seq, [rstseq] (allocating memory for rstseq)
!  nrstr_pos, [rstpos] (allocating memory for rstpos)
!  nrstr_dist, [rstdis] (allocating memory for rstdis)
!  nrstr_wall, [rstwal] (allocating memory for rstwal)

! external definition of iargc disabled for gfortran
!integer(4) iargc
!external iargc

! read name of input file from the command line
num_args = command_argument_count()
if (num_args .lt. 1) call die('no input file specified on the command line')
#if defined(CRAY)
call pxfgetarg(num_args, infilename, 200, i)
#elif defined(MPICH)
call getarg(1, infilename)
#else
call getarg(num_args, infilename)
#endif
text = 'Reading input from '//infilename
call centered_heading(trim(text), '-')

initialize = .true. 

if(.not. prm_open_section('md',infilename)) then
call prm_close
! open input file
fu = freefile()
open(unit=fu, file=infilename, action='read', form='formatted', status='old', iostat=fstat)
if (fstat .ne. 0) call die('error opening input file '//infilename)
        initialize = old_initialize(fu)
        close(fu)
return
else

need_restart = .false. !flag for restart file required
if(.not. prm_get_integer_by_key('steps', nsteps)) then
write(*,*) '>>> ERROR: steps not specified (section MD)'
initialize = .false.
end if
if(.not. prm_get_real_by_key('stepsize', stepsize)) then
write(*,*) '>>> ERROR: stepsize not specified (section MD)'
initialize = .false.
end if
write (*,10) nsteps, stepsize
10	format ('Number of MD steps =',i10,'  Stepsize (fs)    =',f10.3)

! convert to internal time units once and for all.
dt=0.020462_prec*stepsize
dt2=0.5_prec*dt
! --- Temperature, Thermostat etc.
if(.not. prm_get_real_by_key('temperature', Temp0)) then
write(*,*) '>>> ERROR: temperature not specified (section MD)'
initialize = .false.
end if

if(.not. prm_get_real_by_key('bath_coupling', tau_T)) then
write(*,*) 'Temperature bath relaxation time tau_T set to default'
tau_T = tau_T_default
end if

write (*,15) Temp0,tau_T
tau_T=0.020462_prec*tau_T
if(Temp0 <= 0) then
write(*,'(a)') &
        '>>> Error: No dynamics at zero temperature!'
initialize = .false.
end if
if(tau_T < dt) then
write(*,'(a)') '>>> Error: tau_t must be >= stepsize.'
initialize = .false.
end if

if(.not. prm_get_string_by_key('thermostat', name_thermostat)) then
write(*,*) 'No thermostat chosen. Berendsen thermostat will be used.'
thermostat = BERENDSEN

else if( name_thermostat == 'langevin' ) then
write(*,*) 'Thermostat chosen: ', name_thermostat
thermostat = LANGEVIN
if(.not. prm_get_real_by_key('langevin_friction', friction)) then
        friction = one/tau_T !***according to GROMACS manual, this is their default value. Need to check - A. Barrozo
!		gkT = 2*friction*Boltz*Temp0/dt !constant to be used to generate the random forces of the thermostat SET LATER	
        write(*,*) 'Langevin thermostat friction constant set to default: 1/tau_T'
end if
if(.not. prm_get_logical_by_key('langevin_random', langevin_gauss_rand )) then
!we use uniform random noise by default for random numbers in the
!langevin thermostat, Paul October 2014
        write(*,*) 'Using uniform random numbers for Langevin thermostat'
else if (langevin_gauss_rand) then
        write(*,*) 'Using Gaussian random numbers for Langevin thermostat'
else if (.not.langevin_gauss_rand) then
         write(*,*) 'Using uniform random numbers for Langevin thermostat'
end if
else if( name_thermostat == 'nose-hoover' ) then
write(*,*) 'Thermostat chosen: Nose-Hoover'
        thermostat = NOSEHOOVER
kbT = Boltz*Temp0

else if( name_thermostat /= 'berendsen' .and. name_thermostat /= 'langevin' .and. name_thermostat /= 'nose_hoover' ) then
write(*,*) '>>> ERROR: this thermostat does not exist in Q.',name_thermostat
initialize = .false.
else
        thermostat = BERENDSEN
        name_thermostat = 'berendsen'
write(*,*) 'Thermostat chosen: ', name_thermostat
end if

if(.not. prm_get_integer_by_key('nhchains', numchain) .and. thermostat == NOSEHOOVER ) then
numchain = 10
write(*,*) 'Nose-Hoover thermostat chain number set to default: 10'
end if

if(.not. prm_get_real_by_key('nose-hoover_mass', nhq) .and. thermostat == NOSEHOOVER ) then
nhq = kbT/(tau_T*tau_T)
write(*,*) 'Nose-Hoover thermostat mass set to default: kbT/tau_T^2' !Based on Martyna, Klein and Tuckerman J.Chem. Phys. 92
end if


if(.not. prm_get_string_by_key('integrator', name_integrator)) then
write(*,*) 'Leap frog integrator used by default.'
integrator = LEAPFROG

else if( name_integrator == 'velocity-verlet' ) then
write(*,*) 'Integrator: ', name_integrator
integrator = VELVERLET
else if( name_integrator /= 'leap-frog' .and. name_integrator /= 'velocity-verlet' ) then
write(*,*) '>>> ERROR: no such integrator exists in Q.',name_integrator
initialize = .false.
else
        name_integrator = 'leap-frog'
        integrator = LEAPFROG
write(*,*) 'Integrator: ', name_integrator
end if
yes = prm_get_logical_by_key('separate_scaling', separate_scaling, .true.)
if(separate_scaling) then
   write(*,'(a)') 'Solute and solvent atoms coupled separately to heat bath.'
   write(*,'(a)') 'May redefine groups in later section'
else 
   write(*,'(a)') 'Solute and solvent atoms coupled together to heat bath.'
   ntgroups = 1
end if

15	format ('Target temperature =',f10.2,'  T-relax time     =',f10.2)

yes = prm_get_integer_by_key('random_seed', iseed, 1) 
if(.not. prm_get_real_by_key('initial_temperature', Tmaxw)) then
iseed = 0 !set iseed = 0 if no initial temp
need_restart = .true.
end if

if (iseed > 0) write (*,16) Tmaxw, iseed
16	format ('Initial velocities will be generated from Maxwell distribution:',&
        /,'Maxwell temperature=',f10.2,' Random number seed=',i10)

! --- shake, LRF
if(.not. prm_get_logical_by_key('shake_solvent', shake_solvent, .true.)) then
write(*,'(a)') '>>> Error: shake_solvent must be on or off.'
initialize = .false.
end if
!write(*,17) 'all solvent bonds', onoff(shake_solvent)
17	format('Shake constraints on ',a,' to ',a,' : ',a3)

if(.not. prm_get_logical_by_key('shake_solute', shake_solute, .false.)) then
write(*,'(a)') '>>> Error: shake_solute must be on or off.'
initialize = .false.
end if 
!write(*,17) 'all solute bonds', onoff(shake_solute)

if(.not. prm_get_logical_by_key('shake_hydrogens', shake_hydrogens, .true.)) then
write(*,'(a)') '>>> Error: shake_hydrogens must be on or off.'
initialize = .false.
end if 

if(.not. prm_get_logical_by_key('shake_heavy', shake_heavy, .false.) ) then
write(*,'(a)') '>>> Error: shake_heavy must be on or off.'
initialize = .false.
end if


if(prm_get_logical_by_key('shake_all',shake_all)) then
if(shake_all) then
       shake_solute    = .true.
       shake_solvent   = .true.
       shake_hydrogens = .true.
       shake_heavy     = .true.
else
       shake_solute    = .false.
       shake_solvent   = .false.
       shake_hydrogens = .false.
       shake_heavy     = .false.
end if
end if

if(prm_get_logical_by_key('shake_all_solvent',shake_all_solvent)) then
if(shake_all_solvent) then
       shake_solvent   = .true.
       shake_hydrogens = .true.
       shake_heavy     = .true.
else
       shake_solvent   = .false.
end if
end if

if(prm_get_logical_by_key('shake_all_solute',shake_all_solute)) then
if(shake_all_solute) then
       shake_solute    = .true.
       shake_hydrogens = .true.
       shake_heavy     = .true.
else
       shake_solute    = .false.
end if
end if

if(prm_get_logical_by_key('shake_all_hydrogens',shake_all_hydrogens)) then
if(shake_all_hydrogens) then
        shake_solute    = .true.
        shake_solvent   = .true.
        shake_hydrogens = .true.
else
        shake_hydrogens = .false.
end if
end if

if(prm_get_logical_by_key('shake_all_heavy',shake_all_heavy)) then
if(shake_all_heavy) then
         shake_solute  = .true.
         shake_solvent = .true.
         shake_heavy   = .true.
else
         shake_heavy   = .false.
end if
end if

if(shake_solvent) then
if(shake_heavy) write(*,17) 'solvent bonds','heavy atoms',onoff(shake_heavy)
if(shake_hydrogens) write(*,17) 'solvent bonds','hydrogens',onoff(shake_hydrogens)
else
write(*,17) 'solvent bonds','any atoms',onoff(shake_solvent)
endif
if(shake_solute) then
if(shake_heavy) write(*,17) 'solute bonds','heavy atoms',onoff(shake_heavy)
if(shake_hydrogens) write(*,17) 'solute bonds','hydrogens',onoff(shake_hydrogens)
else
write(*,17) 'solute bonds','any atoms',onoff(shake_solute)
endif


yes = prm_get_logical_by_key('lrf', use_LRF, .true.)
if(use_LRF) then
       write(*,20) 'LRF Taylor expansion outside cut-off'
else 
       write(*,20) 'standard cut-off'
end if

20	format ('Nonbonded method   = ',a)

yes = prm_get_logical_by_key('force_rms', force_rms, .false.)
if(force_rms) then
write(*,22) 
end if
22	format ('R.M.S. force will be calculated.')

end if

if(.not.prm_open_section('TGROUPS')) then
        if(separate_scaling) then
                ! two groups, solute and solvent
                ntgroups = 2
                ntgroups_kind = DEFAULT_TWO
                allocate(tscale(ntgroups))
        else
                ! one group only for all atoms
                ntgroups = 1
                ntgroups_kind = DEFAULT_ONE
                allocate(tscale(ntgroups))
        end if
        write(*,'(a)') 'Using default temperature coupling groups'
else
        if(.not. separate_scaling) then
                write(*,'(a)') 'Only one temperature coupling group requested, skipping'
                ntgroups = 1
                ntgroups_kind = DEFAULT_ONE
                allocate(tscale(ntgroups))
        else
                ntgroups = prm_count('TGROUPS')
                ! test for bad input, default to solute, solvent groups in this case
                if(ntgroups .le. 1) then
                        write(*,'(a)') 'Weird values for temperature groups, defaulting to two groups'
                        ntgroups = 2
                        allocate(tscale(ntgroups))
                        ntgroups_kind = DEFAULT_TWO
                else
                ! actuallu use the user input for the groups        
                        allocate(tscale(ntgroups))
                        ntgroups_kind = USERSET
                        write(*,'(a,i4)') 'Using this number of temperature control groups: ',ntgroups
                ! check for later to see if all atoms are covered
                        do tgroups = 1 , ntgroups
                                yes = prm_get_line(text)
                                read(text,*, iostat=fstat) tscale(tgroups)%starta,tscale(tgroups)%enda, &
                                        tscale(tgroups)%tname
                                if(fstat /= 0) then
                                        write(*,'(a,i4)') '>>> ERROR: Invalid data for temperature scaling group ',tgroups
                                        initialize = .false.
                                        exit
                                end if
                                ! check if groups overlap and remind the user that this does not work this way
                                ! always this group against previous, might need more general test later
                                ! groups need to be sequential at this point
                                if(tgroups .ge. 2) then
                                        if(tscale(tgroups)%starta .lt. tscale(tgroups-1)%enda) then
                                                write(*,'(a,i4,a,i4)') '>>> ERROR: Overlapping atoms for groups ',tgroups,' and ',tgroups-1
                                                initialize = .false.
                                                exit
                                        end if
                                end if
                                write(*,789) 'Temperature scaling group ',tgroups
                                write(*,789) 'Starting atom = ',tscale(tgroups)%starta
                                write(*,789) 'End atom = ',tscale(tgroups)%enda
                        end do
                end if
        end if
end if
789     format(a,i4)

if(.not. prm_open_section('PBC')) then
box = .false.
write(*,'(a)') 'Boundary: sphere'
else
        box = .true.
        write(*,'(a)') 'Boundary: periodic box'
        if( .not. prm_get_logical_by_key('rigid_box_centre', rigid_box_centre, .false. ) ) then
                write(*,'(a)') '>>> Error: rigid_box_centre must be on or off'
                initialize = .false.
        end if
        write(*,'(a,a3)') 'Rigid box centre ', onoff(rigid_box_centre)
        if( .not. prm_get_logical_by_key('constant_pressure', constant_pressure, .false.) ) then
                write(*,'(a)') '>>> Error: constant_pressure must be on or off'
                initialize = .false.
        end if
        if( constant_pressure ) then
                write(*,'(a)') 'NPT-ensemble'
                volume_try = 0
                volume_acc = 0
                if( .not. prm_get_real_by_key('max_volume_displ', max_vol_displ) ) then
                        initialize = .false.
                        write(*,'(a)') '>>> ERROR: maximum volume displacement not specified (section PBC)'
                else
                        write(*,5) max_vol_displ
                end if
5	format ('Maximum volume displacemet = ', f10.3)
                if( .not. prm_get_integer_by_key('pressure_seed', pressure_seed)) then
                        pressure_seed = 3781
                end if
				write(*, '(a, i4 )' ) 'Pressure seed: ', pressure_seed 
                if( .not. prm_get_real_by_key('pressure', pressure) ) then
                        pressure = one
                end if
                write(*,9) pressure
9	format ('Pressure = ',f10.3,'  bar')
                !convert pressure to strange internal unit
                pressure = pressure * 1.43836e-5_prec
                yes = prm_get_logical_by_key('atom_based_scaling', atom_based_scaling, .false.)
                if (atom_based_scaling) then
                        write (*,'(a)') 'Coordinate scaling on volume changes:  Atom based'
                else
                        write (*,'(a)') 'Coordinate scaling on volume changes:  Molecule based'
                end if 
        else
                write(*,'(a)') 'NVT-ensemble'
                if( prm_get_line_by_key('control_box', instring) ) then
                        read(instring, *) new_boxl
                        control_box = .true.
                        write(*,'(a, 3f10.3)')'Boxsize will be changed to: ', new_boxl
                else
                        control_box = .false.
                end if 
        end if !section constant_pressure 
	yes = prm_get_logical_by_key('put_solvent_back_in_box', put_solvent_back_in_box)
	yes = prm_get_logical_by_key('put_solute_back_in_box', put_solute_back_in_box) 
	if(put_solute_back_in_box .and. put_solvent_back_in_box) then
           write(*,'(a)') 'Solute and solvent molecules will be put back in box.'
	else
		if (put_solute_back_in_box) then
           	write(*,'(a)') 'Only solute molecules will be put back in box.'
		else
			if (put_solvent_back_in_box) then
				write(*,'(a)') 'Only solvent molecules will be put back in box.'
			else
				write(*,'(a)') 'No molecules will be put back in box.'
			end if
		end if
	end if
end if !section PBC

! --- Rcpp, Rcww, Rcpw, Rcq, RcLRF
! change for PBC, q-atoms interact with everything for now, placeholder until we have 
! fixed it to calculate the interaction between the mirror images
if(.not. prm_open_section('cut-offs')) then
write(*,'(a)') 'No cut-offs section, default cut-offs used'
Rcpp = rcpp_default
Rcww = rcww_default
Rcpw = rcpw_default
if (box) then
Rcq = rcq_default_pbc
else
Rcq = rcq_default_sph
end if
RcLRF = rcLRF_default
else
if(.not. prm_get_real_by_key('solute_solute', Rcpp, rcpp_default)) then
        write(*,'(a)') 'solute-solute cut-off set to default'
end if
if(.not. prm_get_real_by_key('solvent_solvent', Rcww, rcww_default)) then
        write(*,'(a)') 'solvent-solvent cut-off set to default'
end if
if(.not. prm_get_real_by_key('solute_solvent', Rcpw, rcpw_default)) then
        write(*,'(a)') 'solute-solvent cut-off set to default'
end if
if (box) then
if(.not. prm_get_real_by_key('q_atom', Rcq, rcq_default_pbc)) then
        write(*,'(a)') 'q-atom cut-off set to default for PBC'
end if
else
if(.not. prm_get_real_by_key('q_atom', Rcq, rcq_default_sph)) then
        write(*,'(a)') 'q-atom cut-off set to default for sphere'
end if
end if
if(use_LRF) then
        if(.not. prm_get_real_by_key('lrf', rcLRF, rcLRF_default)) then
                write(*,'(a)') 'LRF cut-off set to default'
        end if
        if(RcLRF < rcpp .or. RcLRF < rcpw .or. RcLRF < rcww) then
                if (box) then
                      rcLRF = rcLRF_default_pbc
                      write(*,'(a)') 'LRF cut-off set to default for PBC'
                else
                      write(*,'(a)') &
                        '>>> ERROR; LRF cut-off must not be smaller than solute or solvent cut-offs!'
                      initialize = .false.
                end if
        end if
end if
end if

write (*,25) Rcpp,Rcww,Rcpw,Rcq
if(use_LRF) write(*,26) RcLRF
25	format ('Cut-off radii for non-bonded interactions:',/, &
        'Solute-solute:    ',f6.2,/,&
        'Solvent-solvent:  ',f6.2,/,&
        'Solute-solvent:   ',f6.2,/,&
        'Q-atom-non-Q-atom:',f6.2)
26	format ('LRF:              ',f6.2)

30	format ('>>> WARNING: Ignoring obsolete keyword ',a,'.')
! --- simulation sphere

if( .not. box ) then
if(.not. prm_open_section('sphere')) then
                fk_pshell = fk_pshell_default
                print*,'Radius of inner restrained shell set to 85% of exclusion shell radius.'
                rexcl_i = shell_default
                write(*,50) rexcl_i
else
        if(prm_get_line_by_key('centre', instring)) then
                write(*,30) 'centre'
        end if
        ! --- rexcl_o, rexcl_i, fk_pshell
        if(prm_get_real_by_key('radius', rjunk)) then
                write(*,30) 'radius'
        end if
        if(prm_get_real_by_key('shell_radius', rexcl_i)) then  !inner radius of restrained shell
            write(*,50) rexcl_i
            if(rexcl_i < zero) then
              call die('inner radius of restrained shell must be >= 0')
            end if
        else
            print*,'Radius of inner restrained shell set to 85% of exclusion shell radius.'
                                rexcl_i = shell_default
            write(*,50) rexcl_i
        end if
50 format('Radius of inner restrained shell       =    ',f8.3) 
        if(.not. prm_get_real_by_key('shell_force', fk_pshell)) then
                write(*,'(a)') 'Shell force constant set to default'
                fk_pshell = fk_pshell_default
        end if
        if(fk_pshell > 0) then
                write(*,47) fk_pshell
        end if
47		format('Shell restraint force constant         =',f8.2)
        if(.not. prm_get_real_by_key('excluded_force', fk_fix   )) then
                write(*,'(a)') 'Excluded atoms force constant set to default'
                fk_fix = fk_fix_default
        end if
        if(fk_fix > 0) then
                write(*,48) fk_fix
        else
                fk_fix = fk_fix_default
                write(*,'(a)')'Shell restraint force constant can not be less-equal zero, set to default'
                write(*,48) fk_fix
        end if
48              format('Excluded atom restraint force constant         =',f8.2)

        yes = prm_get_logical_by_key('excluded_freeze', freeze, .false.)
        if(freeze) then
                write(*,'(a)') &
                        'Excluded atoms will not move.'
        end if

        yes = prm_get_logical_by_key('exclude_bonded', exclude_bonded, .false.)
        if(exclude_bonded) then
                write(*,'(a)') &
                        'Bonded interactions outside the sphere will be eliminated'
        end if
end if

! --- solvent 
inlog = prm_open_section('solvent')
if(.not. inlog) inlog = prm_open_section('water') !try also the old name
if(.not. inlog) then       !defaults
        fk_wsphere = -one
        Dwmz = -one
        awmz = -one
        wpol_restr = wpol_restr_default
        wpol_born = wpol_restr_default
        fkwpol = -one
else
        if(prm_get_real_by_key('radius', rwat_in)) then
                write(*,'(a,f8.2)') 'Target solvent radius =',rwat_in
        end if
        if(prm_get_line_by_key('centre', instring)) then
                write(*,30) 'centre'
        end if
        if(prm_get_real_by_key('pack', rjunk)) then
                write(*,30) 'pack'
        end if


  if(.not. prm_get_real_by_key('radial_force', fk_wsphere)) then
        write(*,'(a)') 'Solvent radial restraint force constant set to default'
        fk_wsphere = -one ! this will be set in water_sphere, once target radius is known
  end if
  yes=prm_get_logical_by_key('polarisation', wpol_restr, wpol_restr_default)
  !default is on when pol. restr is on, otherwise off
  yes=prm_get_logical_by_key('charge_correction', wpol_born, wpol_restr)
  if(wpol_born .and. .not. wpol_restr) then
        write(*,'(a)') '>>> ERROR: charge_correction on requires polarisation on (section solvent)'
        initialize = .false.
  end if
  if(.not. prm_get_real_by_key('polarisation_force', fkwpol)) then
        write(*,'(a)') 'Solvent polarisation force constant set to default'
        fkwpol = -1 ! this will be set in water_sphere, once target radius is known
  end if
  yes = prm_get_real_by_key('morse_depth', Dwmz, -one)
  yes = prm_get_real_by_key('morse_width', awmz, -one)
  if(prm_get_string_by_key('model', instring)) then
        write(*,30) 'model'
  end if
end if !if (.not. inlog)
end if !if( .not. box )


if(.not. prm_open_section('intervals')) then
write(*,'(a)') 'non-bond list update interval set to default.'
NBcycle = NB_cycle_default
write(*,'(a)') 'energy summary interval set to default.'
iout_cycle = iout_cycle_default
itemp_cycle = iout_cycle_default
iene_cycle = 0 !no energy
itrj_cycle = 0 !no trajectory

ivolume_cycle = ivolume_cycle_default


else
if(.not. prm_get_integer_by_key('non_bond', NBcycle)) then
        write(*,'(a)') 'non-bond list update interval set to default.'
        NBcycle = NB_cycle_default
end if
if(.not. prm_get_integer_by_key('output', iout_cycle)) then
        write(*,'(a)') 'energy summary interval set to default.'
        iout_cycle = iout_cycle_default
end if
if(.not. prm_get_integer_by_key('temperature', itemp_cycle)) then
        write(*,'(a)') 'temperature print-out interval set to default.'
        itemp_cycle = iout_cycle_default
end if
yes = prm_get_integer_by_key('energy', iene_cycle, 0)
yes = prm_get_integer_by_key('trajectory', itrj_cycle, 0)

if( constant_pressure ) then
        if( .not. prm_get_integer_by_key('volume_change', ivolume_cycle) ) then
                write(*,'(a)') 'volume change intervall set to default'
                ivolume_cycle = ivolume_cycle_default
        end if
end if
end if

write(*,84) NBcycle
84	format('Non-bonded pair list update interval   =',i8)
86	format('Energy summary print-out interval      =',i8)
87	format('Temperature print-out interval         =',i8)
88	format('Trajectory write interval              =',i8)
89	format('Energy file write interval             =',i8)
83  format('Volume change interval                 =',i8)

if(iout_cycle > 0) then
write (*,86) iout_cycle
else
write(*,'(a)') 'No energy summaries written.'
iout_cycle = -999999999 ! make sure mod(istep, iout_cycle) never = 0
end if
if(itemp_cycle > 0) then
write (*,87) itemp_cycle
else
write(*,'(a)') 'No temperatures written.'
itemp_cycle = -999999999 ! make sure mod(istep, itemp_cycle) never = 0
end if
if(itrj_cycle > 0) then
write (*,88) itrj_cycle
else
itrj_cycle = -999999999 !no energy
write(*,'(a)') 'No trajectory written.'
end if
if(iene_cycle > 0) then
write (*,89) iene_cycle
else
iene_cycle = -999999999 !no energy
write(*,'(a)') 'No energy file written.'
end if
if( constant_pressure ) then
write(*,83) ivolume_cycle
end if

!read trajectory atom mask
mask_rows = prm_count('trajectory_atoms')
if(itrj_cycle > 0) then
if(mask_rows == 0) then
        write(*,'(a)') 'All atoms will be included in the trajectory.'
        yes = trj_store_mask('all')
else
        do i=1,mask_rows
                yes = prm_get_line(text)
                yes = trj_store_mask(text)
        end do
end if
elseif(mask_rows == 0) then
write(*,'(a)') 'Ignoring section trajectory_atoms.'
end if

if(.not. prm_open_section('files')) then
write(*,'(a)') '>>> ERROR: files section not found.'
initialize = .false.
else
if(.not. prm_get_string_by_key('topology', top_file)) then
        write(*,'(a)') '>>> ERROR: topology not specified (section files)'
        initialize = .false.
end if
write (*,60) trim(top_file)
60		format ('Topology file      = ',a)

if(.not. prm_get_string_by_key('restart', restart_file)) then
        restart = .false.
        if(need_restart) then
                write(*,'(a)') '>>> ERROR: Restart file required when initial temp. not given.'
                initialize = .false.
        end if
else
        restart = .true.
end if

if(restart) then
        write (*,65) trim(restart_file)
else
        write (*,'(a)') 'Initial coordinates taken from topology.'
        if(iseed == 0) then
                write(*,'(a)') &
                        '>>> ERROR: Need a random number seed to generate initial velocities, aborting.'
                initialize = .false.
        end if
end if
65		format ('Initial coord. file= ',a)

if(.not. prm_get_string_by_key('final', xfin_file)) then
        write(*,'(a)') '>>> ERROR: final co-ordinate file not specified (section files, keyword final)'
        initialize = .false.
end if
write (*,80) trim(xfin_file)
80		format ('Final coord. file  = ',a)

if(.not. prm_get_string_by_key('trajectory', trj_file)) then
        if(itrj_cycle > 0) then
                write(*,'(a)') '>>> ERROR: Trajectory file name required to write trajectory!'
                initialize = .false.
        end if
else
        if(itrj_cycle < 0) then
                write(*,*) '>>> Error: Trajectory file given but no output interval'
                initialize = .false.
        end if
        if(itrj_cycle > 0) write (*,90) trim(trj_file)
end if
90		format ('Trajectory file    = ',a)

if(.not. prm_get_string_by_key('energy', ene_file)) then
        if(iene_cycle > 0) then
                write(*,'(a)') '>>> ERROR: Energy file name required to write energies!'
                initialize = .false.
        end if
else

        if(iene_cycle < 0) then

                write(*,'(a)') '>>> ERROR: Energy file given but no energy interval'

                initialize=.false.

        end if
        if(iene_cycle > 0) write (*,94) trim(ene_file)
end if
94		format ('Energy output file = ',a)

if(.not. prm_get_string_by_key('fep', fep_file)) then
        write(*,'(a)') 'No FEP file.'
        !initialize = .false. !This condition IS OK.
        fep_file = ''
else
        write (*,95) trim(fep_file)
95			format ('FEP input file     = ',a,/)
end if
if(.not. prm_get_string_by_key('restraint', exrstr_file)) then
        implicit_rstr_from_file = 0
else
        implicit_rstr_from_file = 1
        write (*,104) trim(exrstr_file)
104			format ('External rstr file = ',a,/)
end if
if(prm_get_string_by_key('water', instring)) then
        write(*,30) 'water'
end if
end if			

! --- states, EQ
nstates = 0
if(prm_open_section('lambdas')) then
do while(prm_get_field(instring))
        nstates = nstates + 1
        read(instring, *, iostat=fstat) lamda_tmp(nstates)
        if(fstat /= 0) then
                write(*,'(a)') '>>> ERROR: Invalid lambda value.'
                initialize = .false.
                exit
        end if
end do
end if
if(nstates == 0 .and. fep_file /= '') then
if(fep_file /= '') then
        write(*,'(a)') 'Defaulting to single FEP state.'
        nstates = 1
        lamda_tmp(1) = one
end if
end if
if(nstates > 0 ) then
if(fep_file == '') then
        write(*,'(a)') '>>> ERROR: FEP file required to use lambdas!'
        initialize = .false.
else
        ! allocate memory for EQ
        allocate(EQ(nstates), stat=alloc_status)
        call check_alloc('Q-atom energy array')

        ! init EQ%lambda
        EQ(1:nstates)%lambda = lamda_tmp(1:nstates)
        write (*,98) (EQ(i)%lambda,i=1,nstates)
98			format ('lambda-values      = ',10f8.5)
end if
else
        ! need dummy EQ allocation for classical MD to work
        allocate(EQ(1),stat=alloc_status)
        call check_alloc('Dummy Q-atom energy array')
end if
!Option to make additional calculation with atom groups excluded from the
!energy calculation to provide 'real' group contribution
!Added Paul Bauer 2014
use_excluded_groups = .false.
ngroups_gc = prm_count('group_contribution')
write (*,'(/,a)') 'List of excluded atom groups'
if ( ngroups_gc .gt. 0 ) then
!allocate new EQ arrays for each group
777 format ('Number of excluded group calculations =',i4)
write (*,777) ngroups_gc
allocate(ST_gc(ngroups_gc), stat=alloc_status)
call check_alloc('gc param store')
778 format ('Groupnumber ',i4,'	contains')
do i=1,ngroups_gc
        write(*,778,advance='no') i
        !read numbers from string
        number = 0
        do while (prm_get_field(instring))
                number = number + 1
                read(instring, *, iostat=fstat) gc_mask_tmp(number)
                if(fstat /= 0) then
                        write(*,'(a)') '>>> ERROR: Invalid mask.'
                        initialize = .false.
                        exit
                end if
        end do

        mask_rows = number - 2
        ST_gc(i)%seltype = trim(gc_mask_tmp(1))
        ST_gc(i)%caltypen = trim(gc_mask_tmp(2))
        if ((trim(ST_gc(i)%seltype) /= 'atom' ).and. (trim(ST_gc(i)%seltype) /= 'residue')) then
779 format (/,'Could not understand name of the selection , ',a)
                write (*,779) trim(ST_gc(i)%seltype)
                initialize = .false.
                exit
        end if
        if ((trim(ST_gc(i)%caltypen) /= 'full').and.(trim(ST_gc(i)%caltypen) /= 'electro') .and. &
                (trim(ST_gc(i)%caltypen) /= 'vdw').and.(trim(ST_gc(i)%caltypen) /= 'all')) then
781 format (/,'Could not understand name of the calculation type , ',a)
                write (*,781) trim(ST_GC(i)%caltypen)
                initialize = .false.
                exit
        end if
        if(iene_cycle > 0) then
        if(mask_rows == 0) then
!You are not as funny as you think ...
                write(*,780)
780 format ('no atoms')
                ST_gc(i)%count=0
        else
                ST_gc(i)%count=mask_rows
                allocate(ST_gc(i)%maskarray(mask_rows),stat=alloc_status)
                call check_alloc('gc maskarray')
                
                do ii=3,mask_rows+2
                        yes = gc_store_mask(ST_gc(i),gc_mask_tmp(ii))
                end do
                use_excluded_groups = .true.
                write (str,'(i4)') i
782 format (i4,' atom groups')
                write(*,782) mask_rows 
        end if
        elseif(mask_rows == 0) then
                write(*,'(a)') 'Ignoring section for excluding atoms.'
        end if
end do

end if

!new section in *inp files to control QCP behaviour
!will only trigger if FEP file is in use -> if more than 0 states
if(nstates > 0 ) then
        if(.not. prm_open_section('QCP')) then
                write(*,'(a)') 'No QCP section found, will not try to use RPMD.'
                use_qcp = .false.
                QCP_N = QCP_OFF
        else
                write(*,'(a)') 'Found QCP section, will use RPMD to describe atoms in Q region.'
                use_qcp = .true.
                if(.not.prm_get_integer_by_key('qcp_seed',qcp_seed)) then
                        write(*,'(a)') 'Using random number for seeding from date_and_time'
                        call date_and_time(values=timeval)
                        qcp_seed = timeval(5)*3600 + timeval(6)*60+timeval(7)+13337
                        qcp_seed = MOD(qcp_seed,10000)
                        if(MOD(qcp_seed,2).eq.0) qcp_seed = qcp_seed + 1
                        write(*,'(a,i4)') 'Using the follwing random number for seeding ',qcp_seed
                else
                        write(*,'(a,i4)') 'Using the follwing random number for seeding ',qcp_seed
                end if
! find out if we use mass perturbation for KIE
                yes = prm_get_logical_by_key('qcp_kie',use_qcp_mass,.false.)
                if(use_qcp_mass) then
                        write(*,'(a)') 'Will perform calculation with mass perturbation for KIE'
                        QCP_N = QCP_ON_KIE
                else
                        write(*,'(a)') 'No mass perturbation'
                        QCP_N = QCP_ON
                end if
! section for methods goes here later
! TODO !
! implement other sampling methods (staging, simple MC, ...)
                yes = prm_get_logical_by_key('qcp_write',qcp_write,.true.)
                if(qcp_write) then
! need file name
                       if(.not.prm_get_string_by_key('qcp_pdb',qcp_pdb_name)) then
                               qcp_write = .false.
                               write(*,'(a)') '>>> WARNING: Need file name for print out of qcp bead coordinates. Disabled write out'
                               write(*,'(a)') 'Keyword: qcp_pdb'
                       else
                               write(*,'(a,a)') 'Writing coordinates to file ',qcp_pdb_name
                       end if
                else
                        write(*,'(a)') 'Not writing out bead coordinates'
                end if
! decide on printout level
               if(.not. prm_get_logical_by_key('qcp_show',qcp_verbose)) then
                       qcp_verbose = .false.
               else
                       if(qcp_verbose) then
                               write(*,'(a)') 'Printing more QCP information'
                       end if
               end if
               if(.not.prm_get_logical_by_key('qcp_debug',qcp_veryverbose)) then
                       qcp_veryverbose = .false.
               else
                       if(qcp_veryverbose) then
                               qcp_verbose = .true.
                               write(*,'(a)') 'Printing all QCP information I can find'
                       end if
               end if
!chose which atoms should be treated as ring polymers
!important later when setting up NB list, RP will be treated different from classical
!this section can be overwritten in the FEP file
		if(.not. prm_get_string_by_key('selection', qcp_select)) then
			write(*,'(a)') 'Will default to Hydrogen atoms only treated as RP!'
			qcp_enum = QCP_HYDROGEN
		else
			call upcase(qcp_select)
			 if ((qcp_select .eq. 'HYDROGEN') .or. &
				(qcp_select .eq. 'HYD') .or. &
				(qcp_select .eq. 'H')) then
				qcp_enum = QCP_HYDROGEN
				write(*,'(a)') 'Only treat Hydrogen atoms as RP!'
			elseif ((qcp_select .eq. 'ALL') .or. &
					(qcp_select .eq. 'QATOM') .or. &
					(qcp_select .eq. 'FEP')) then
				qcp_enum = QCP_ALLATOM
				write(*,'(a)') 'Treat all Q-atoms as RP!'
			elseif (qcp_select .eq. 'INDIVIDUAL') then
				qcp_enum = QCP_FEPATOM
				write(*,'(a)') 'Will use information in FEP file to select QCP atoms'
			else
				write(*,'(a)') ' >>> ERROR: No such QCP atom selection!'
				initialize = .false.
			end if
		end if
!how large should the RP be?
!can again be overwritten in FEP file for each atom itself
		if(.not. prm_get_string_by_key('size',qcp_size_name)) then
			write(*,'(a,i6,a)') 'Will use default sizes for RP, ',QCP_SIZE_DEFAULT, ' ring beads per atom!'
			qcp_size = QCP_SIZE_DEFAULT
		else
			call upcase(qcp_size_name)
			if(qcp_size_name .eq. 'DEFAULT') then
				write(*,'(a,i6,a)') 'Will use default sizes for RP, ',QCP_SIZE_DEFAULT,' ring beads per atom!'
				qcp_size = QCP_SIZE_DEFAULT
			else if (qcp_size_name .eq. 'SMALL') then
				write(*,'(a,i6,a)') 'Will use small RP, ',QCP_SIZE_SMALL,' ring beads per atom!'
				qcp_size = QCP_SIZE_SMALL
			else if (qcp_size_name .eq. 'LARGE') then
				write(*,'(a,i6,a)') 'Will use large RP, ',QCP_SIZE_LARGE,' beads per atom!'
				qcp_size = QCP_SIZE_LARGE
                       else if (qcp_size_name .eq. 'USERDEFINE') then
                                yes = prm_get_integer_by_key('qcp_size_user',qcp_size,QCP_SIZE_DEFAULT)
                                if (qcp_size .lt. 4) then
                                        write(*,'(a,i6)') 'Can not use size of ',qcp_size
                                        write(*,'(a,i6)') 'Resorting to smallest default value of ',QCP_SIZE_VERYSMALL
                                else
                                        
                                        write(*,'(a,i6,a)') 'Will use RP with ',qcp_size,' beads per atom!'
                                end if
			else
				write(*,'(a)') ' >>> ERROR: No such QCP size selection!'
				initialize = .false.
			end if
		end if
! now we know how many beads, set default for qcp_level
! for bisection, this is 2**qcp_level .eq. qcp_size
                qcp_level = nint(q_log2(real(qcp_size,kind=prec)))
                write(*,'(a,i4,a,i4)') 'Setting bisection level to default, log_2(',qcp_size,'), = ',qcp_level
                if (real(qcp_level,kind=prec) - q_log2(real(qcp_size,kind=prec)) .gt. ddiff) then
! whoops, the bead size does not work for bisection
                        write(*,'(a)') 'Bead size is not the result of 2**n. Does not work with bisection strategy'
                        initialize = .false.
                end if
!number of PI steps at each calculation
		if(.not. prm_get_integer_by_key('equilibration',qcp_steps(1))) then
			write(*,'(a,i6,a)') 'Will use default number of PI steps, n = ',QCP_steps_default,' for equilibration'
			qcp_steps(1) = QCP_steps_default
		else
			if(qcp_steps(1) .lt. 1) then
				write(*,'(a)') ' >>> ERROR: Can not use less than 1 PI step per classical step!'
				initialize = .false.
			end if
			write(*,'(a,i4,a)') 'Will use the following number of PI steps, n = ',qcp_steps(1),' for free particle equilibration!'
		end if
                if(.not.prm_get_integer_by_key('sampling',qcp_steps(2))) then
                        write(*,'(a,i6,a)') 'Will use default number of PI steps, n = ',QCP_steps_default,' for sampling'
                        qcp_steps(2) = QCP_steps_default
                else
                        if(qcp_steps(2) .lt.1 ) then
                                write(*,'(a)') ' >>> ERROR: Can not use less than 1 PI step per classical step!'
                                initialize = .false.
                        end if
                        write(*,'(a,i4,a)') 'Will use the following number of PI steps, n = ',qcp_steps(2),' for free particle sampling!'
                end if
	end if
else
	write(*,'(a)') 'No RPMD in classical MD'
        use_qcp = .false.
        QCP_N = QCP_OFF
end if
!	--- restraints:
write (*,'(/,a)') 'Listing of restraining data:'

! --- nrstr_seq, [rstseq]
nrstr_seq = prm_count('sequence_restraints')
109 format (/,'No. of sequence restraints =',i10)
if ( nrstr_seq .gt. 0 ) then
! allocate memory for rstseq
write (*,109) nrstr_seq
allocate(rstseq(nrstr_seq), stat=alloc_status)
call check_alloc('restraint list')
write (*,110)
110		format ('  atom_i  atom_j      fc  H-flag to_centre')
do i=1,nrstr_seq
        ! read rstseq(i)
        yes = prm_get_line(text)
        rstseq(i)%to_centre = 0 
        read(text,*, end=111, err=111) rstseq(i)
111			write(*,112) rstseq(i)
112			format (2i8,f8.2,i8,i10)
end do
end if

! --- nrstr_pos, [rstpos]
nrstr_pos = prm_count('atom_restraints')
115 format (/,'No. of position restratints =',i10)
if ( nrstr_pos .gt. 0 ) then
write (*,115) nrstr_pos
! allocate memory for rstpos
allocate(rstpos(nrstr_pos), stat=alloc_status)
call check_alloc('restraint list')
write (*,120)
120		format ('atom_i      x0      y0      z0     fcx     fcy     fcz   state')
do i=1,nrstr_pos ! read rstpos(i)
        yes = prm_get_line(text)
        read(text,*, iostat=fstat) rstpos(i)%i,rstpos(i)%x, &
                rstpos(i)%fk, rstpos(i)%ipsi
        if(fstat /= 0) then
                write(*,'(a)') '>>> ERROR: Invalid atom restraint data.'
                initialize = .false.
                exit
        end if
        write (*,122) rstpos(i)%i,rstpos(i)%x, &
                rstpos(i)%fk, rstpos(i)%ipsi
end do
122		format (i6,6f8.2,i8)
end if

! --- nrstr_dist, [rstdis]
nrstr_dist = prm_count('distance_restraints')
125	format (/,'No. of distance restraints =',i10)
if ( nrstr_dist .gt. 0 ) then
write (*,125) nrstr_dist
! allocate memory for rstdis
allocate(rstdis(nrstr_dist), stat=alloc_status)
call check_alloc('restraint list')
write (*,130)
130		format ('atom_i atom_j   dist1   dist2   fc        state') 
do i=1,nrstr_dist
        yes=prm_get_line(text)
        ! read rstdis(i)
        if(scan(text, ':') > 0) then !got res:atnr
          !Store in i&j as res:atnr and assign atom nr after topology is read (prep_coord)
          read(text,*, iostat=fstat) rstdis(i)%itext,rstdis(i)%jtext,rstdis(i)%d1,& 
             rstdis(i)%d2, rstdis(i)%fk, rstdis(i)%ipsi
        else !Plain numbers
          read(text,*, iostat=fstat) rstdis(i)%i,rstdis(i)%j,rstdis(i)%d1,&
                rstdis(i)%d2, rstdis(i)%fk, rstdis(i)%ipsi
          rstdis(i)%itext = 'nil'
          rstdis(i)%jtext = 'nil'
        end if
        if(fstat /= 0) then
          write(*,'(a)') '>>> ERROR: Invalid distance restraint data.'
          initialize = .false.
          exit
        end if
        write (*,132) rstdis(i)%i,rstdis(i)%j,rstdis(i)%d1,rstdis(i)%d2,rstdis(i)%fk, &
                rstdis(i)%ipsi
end do
132		format (i6,1x,i6,3f8.2,i8)
end if

! --- nrstr_angl, [rstang]
nrstr_angl = prm_count('angle_restraints')
135     format (/,'No. of angle restraints =',i10)
if ( nrstr_angl .gt. 0 ) then
write (*,135) nrstr_angl
! allocate memory for rstang
allocate(rstang(nrstr_angl), stat=alloc_status)
call check_alloc('restraint list')
write (*,140)
140             format ('atom_i atom_j atom_k   angle   fc        state')
do i=1,nrstr_angl
        yes=prm_get_line(text)
          read(text,*, iostat=fstat) rstang(i)%i,rstang(i)%j,rstang(i)%k,&
                rstang(i)%ang, rstang(i)%fk, rstang(i)%ipsi
        if(fstat /= 0) then
          write(*,'(a)') '>>> ERROR: Invalid angle restraint data.'
          initialize = .false.
          exit
        end if
        write (*,142) rstang(i)%i,rstang(i)%j,rstang(i)%k,rstang(i)%ang,rstang(i)%fk, &
                rstang(i)%ipsi
end do
142             format (i6,1x,i6,1x,i6,2f8.2,i8)
end if


if (.not. box )then
! --- nrstr_wall, [rstwal]
nrstr_wall = prm_count('wall_restraints')
145 format (/,'No. of wall sequence restraints=',i10)
if ( nrstr_wall .gt. 0) then
write (*,145) nrstr_wall
! allocate memory for rstwal
allocate(rstwal(nrstr_wall), stat=alloc_status)
call check_alloc('restraint list')
write (*,150)
150  format ('atom_i atom_j   dist.      fc  aMorse  dMorse  H-flag')
do i=1,nrstr_wall
        ! read rstwal(:)
        yes = prm_get_line(text)
        read(text,*, iostat=fstat) rstwal(i)%i,rstwal(i)%j,rstwal(i)%d,rstwal(i)%fk, &
                rstwal(i)%aMorse, rstwal(i)%dMorse, rstwal(i)%ih
        if(fstat /= 0) then
                write(*,'(a)') '>>> ERROR: Invalid wall restraint data.'
                initialize = .false.
                exit
        end if
        write (*,152) rstwal(i)     
end do
152  format (i6,1x,i6,4f8.2,i8)
end if
end if

call prm_close
end function initialize


!-----------------------------------------------------------------------------------


logical function old_initialize(fu)
!arguments
integer						::	fu
! local variables
integer						::	iuse_indip, shake_flag
character					:: text*200, watmodel*80
integer						:: i,j,length
integer						:: irestart
real(kind=prec)						:: stepsize
real(kind=prec)						:: lamda_tmp(max_states)
integer						:: fstat
integer						::	NBMethod
integer						::	iwpol_restr
real(kind=prec)						::	rjunk

!this is called by initialize to read old-style input file which is 
!alreadu open as unit fu

! this subroutine will init:
!  nsteps, stepsize, dt
!  Temp0, tau_T, iseed, Tmaxw
!  usr_LRF, NBcycle, Rcpp, Rcww, Rcpw, Rcq
!  shake_solvent, shake_solute, shake_hydrogens
!  fk_pshell


!  fk_wsphere=-1, wpol_restr, wpol_born fkwpol=-1, Dwmz=-1, awmz=-1
!    (values initialized to -1 will be set in water_sphere, 
!    once target radius is known)
!  top_file
!  restart, [restart_file]
!  xfin_file
!  itrj_cycle, iene_cycle, iout_cycle, itemp_cycle [trj_file], [ene_file]
!  fep_file
!  nstates, EQ (allocating memory for EQ)
!  implicit_rstr_from_file, [exrstr_file]
!  nrstr_seq, [rstseq] (allocating memory for rstseq)
!  nrstr_pos, [rstpos] (allocating memory for rstpos)
!  nrstr_dist, [rstdis] (allocating memory for rstdis)
!  nrstr_wall, [rstwal] (allocating memory for rstwal)


write(*,1)
1	format('>>> WARNING: Entering unsupported compatibility mode',/ &
   '             to read version 2 input file.',/,&
   '             New features unavailable.')

old_initialize = .true. 

!use default values for new features not in old kind of input.
RcLRF = 999.0_prec
exclude_bonded = .false.
force_rms = .false.
shake_hydrogens = .false.
itemp_cycle = iout_cycle_default
awmz = -1

! --- nsteps, stepsize
!read (fu,*, iostat=stat) nsteps,stepsize
if (.not. prm_get_int_real(nsteps,stepsize)) then
old_initialize = .false.
        call die("Wrong input format.")
end if

write (*,10) nsteps, stepsize
10	format ('Number of MD steps =',i10,'  Stepsize (fs)    =',f10.3)

! convert to internal time units once and for all.
dt=0.020462_prec*stepsize
dt2=0.5_prec*dt

! --- Temp0, tau_T, iseed, Tmaxw
read(fu,'(a80)') text !read line into buffer
!now read buffer (avoid reading more lines from inpunt in search for more values)
read(text,*, err=17, end=17) Temp0,tau_T, iseed,Tmaxw 
17	write (*,15) Temp0,tau_T
if(Temp0 <= 0) then
write(*,'(a)') &
        '>>> Error: No dynamics at zero temperature! Aborting.'
old_initialize = .false.
end if

15	format ('Target temperature =',f10.2,'  T-relax time     =',f10.2)
if (iseed > 0) write (*,16) Tmaxw, iseed
16	format ('Initial velocities will be generated from Maxwell distribution:',&
        /,'Maxwell temperature=',f10.2,' Random number seed=',i10)
tau_T=0.020462_prec*tau_T
if(tau_T < dt) then
write(*,'(a)') '>>> Error: tau_t must be >= stepsize.'
old_initialize = .false.
end if

! --- NBmethod, NBcycle, Rcpp, Rcww, Rcpw, Rcq
read (fu,*) NBmethod,NBcycle,Rcpp,Rcww,Rcpw,Rcq
if(NBMethod == 2) then
use_LRF = .true.
else 
use_LRF = .false.
end if
write (*,20) NBmethod,NBcycle
20	format ('Nonbonded method   =',i10,'  NB update cycle  =',i10,/)
write (*,25) Rcpp,Rcww,Rcpw,Rcq
25	format ('Cutoffs are: Rcpp  =',f6.2,'  Rcww =',f6.2,'  Rcpw =', &
f6.2,'  Rcqp =',f6.2,/)

! --- shake_flag
read (fu,*) shake_flag
shake_solvent = .false.
shake_solute = .false.
if(shake_flag >= 1) then
shake_solvent = .true.
end if
if(shake_flag == 2) then
shake_solute = .true.
end if

write (*,30) shake_flag
30	format ('Shake method       =',i10)

! --- iuse_indip
read (fu,*) 
write (*,35) 
35	format ('Ignoring induced dipole flag.')

! --- protein center: xpcent(:)
read (fu,*) 
write (*,40)
40	format ('Ignoring solute centre.')

! --- rexcl_o, rexcl_i, fk_pshell
read (fu,*) rjunk, rjunk, fk_pshell
write(*,44)
write (*,45) fk_pshell
44	format ('Ignoring exclusion and shell radii.')
45	format ('Restrained shell force const.          =',f8.2)

! --- water center: xwcent(:)
read (fu,*) 
write (*,50) 
50	format ('Ignoring solvent centre.')

! set default values before reading
! done this way because the SGI compiler initialises values to be read to zero

read(fu,'(a80)') text ! read line into buffer
! now read buffer (avoid reading more lines from input in search for more values)
read(text, fmt=*, err=58, end=58) rjunk, rjunk, fk_wsphere, iwpol_restr, fkwpol, Dwmz
goto 59

! set default values:
58	Dwmz = -1
if (fkwpol .eq. 0) then
fkwpol = -1
if (fk_wsphere .eq. 0) then
  fk_wsphere = -1
end if
end if
59	if(iwpol_restr == 0) then
wpol_restr = .false.
wpol_born = .false.
elseif(iwpol_restr == 1) then
wpol_restr = .true.
wpol_born = .true.
elseif(iwpol_restr == 2) then
wpol_restr = .true.
wpol_born = .false.
else
call die('unknown water polarisation restraining mode')
end if
write(*,57)
57	format('Ignoring solvent radius and min. packing distance.')
! --- top_file
read (fu,'(a80)') text
call get_fname (text,length,top_file)
write (*,60) top_file(1:length)
60	format ('Topology file      = ',a)

! --- restart, [restart_file]
read (fu,*) irestart
if ( irestart .eq. 1 ) then
restart = .true.
read (fu,'(a80)') text
call get_fname (text,length,restart_file)
write (*,65) restart_file(1:length)
else
restart = .false.
write (*,'(a)') 'Initial coordinates taken from topology.'
if(iseed == 0) then
        write(*,'(a)') &
                'Error: Need a random number seed to generate initial velocities, aborting.'
        call die('invalid data in input')
end if
end if
65 format ('Initial coord. file= ',a)

! --- xfin_file
read (fu,'(a80)') text
call get_fname (text,length,xfin_file)
write (*,80) xfin_file(1:length)
80 format ('Final coord. file  = ',a,/)

! --- itrj_cycle, iene_cycle, iout_cycle, [trj_file], [ene_file]
read (fu,*) itrj_cycle, iene_cycle, iout_cycle
write (*,85) itrj_cycle, iene_cycle, iout_cycle
85 format ('Trajectory, Energy and Output cycles   =',3i8,/)

if ( itrj_cycle .gt. 0 ) then
read (fu,'(a80)') text
call get_fname (text,length,trj_file)
write (*,90) trj_file(1:length)
else
write (*,'(a)') 'No trajectory written.'
itrj_cycle = -999999999 !make sure mod(istep, itrj_cycle) never = 0
end if
90 format ('Trajectory file    = ',a)

if ( iene_cycle .gt. 0 ) then
read (fu,'(a80)') text
call get_fname (text,length,ene_file)
write (*,94) ene_file(1:length)
else
write (*,'(a)') 'No energy file written'
iene_cycle = -999999999 ! make sure mod(istep, iene_cycle) never = 0
end if
94 format ('Energy output file = ',a)

if(iout_cycle == 0) then
write(*,'(a)') 'No energy summaries written.'
iout_cycle = -999999999 ! make sure mod(istep, iout_cycle) never = 0
end if

! --- fep_file
read (fu,'(a80)') text
call get_fname (text,length,fep_file)
write (*,95) fep_file(1:length)
95 format ('FEP input file     = ',a,/)

! --- nstates, EQ
read (fu,*) nstates, (lamda_tmp(i),i=1,nstates)
if ( nstates .gt. 0 ) then
! allocate memory for EQ
allocate(EQ(nstates), stat=alloc_status)
call check_alloc('Q-atom energy array')

! init EQ%lambda
EQ(1:nstates)%lambda = lamda_tmp(1:nstates)
write (*,98) (EQ(i)%lambda,i=1,nstates)
98  format ('lambda-values      = ',10f8.5)
end if

!	--- restraints:
write (*,'(/,a)') 'Listing of restraining data:'

! --- implicit_rstr_from_file, [exrstr_file]
read (fu,*) implicit_rstr_from_file
write (*,101) implicit_rstr_from_file
101 format ('Read rstr file     =',i10)
if ( implicit_rstr_from_file .eq. 1 ) then
read (fu,'(a80)') text
call get_fname (text,length,exrstr_file)
write (*,104) exrstr_file(1:length)
else
write (*,105)
end if
104 format ('External rstr file = ',a,/)
105 format ('Implicit positional restraints from topology.',/)

! --- nrstr_seq, [rstseq]
read (fu,*) nrstr_seq
write (*,109) nrstr_seq
109 format (/,'No. sequence rstrs =',i10)
if ( nrstr_seq .gt. 0 ) then
! allocate memory for rstseq
allocate(rstseq(nrstr_seq), stat=alloc_status)
call check_alloc('restraint list')
write (*,110)
110 format (1x,'  atom_i  atom_j      fc  H-flag to_centre')
end if
do i=1,nrstr_seq
! read rstseq(i)
read (fu,'(a80)') text
rstseq(i)%to_centre = 0 
read(text,*, end=111, err=111) rstseq(i)
111	write(*,112) rstseq(i)
112 format (2i8,f8.2,i8,i10)
end do

! --- nrstr_pos, [rstpos]
read (fu,*) nrstr_pos
write (*,115) nrstr_pos
115 format (/,'No. position rstrs =',i10)
if ( nrstr_pos .gt. 0 ) then
! allocate memory for rstpos
allocate(rstpos(nrstr_pos), stat=alloc_status)
call check_alloc('restraint list')
write (*,120)
end if
120 format ('atom_i      x0      y0      z0     fcx     fcy     fcz  istate')
do i=1,nrstr_pos
! read rstpos(i)
read (fu,*) rstpos(i)%i,rstpos(i)%x, &
rstpos(i)%fk, rstpos(i)%ipsi
write (*,122) rstpos(i)%i,rstpos(i)%x, &
rstpos(i)%fk, rstpos(i)%ipsi
end do
122 format (i6,6f8.2,i8)

! --- nrstr_dist, [rstdis]
read (fu,*) nrstr_dist
write (*,125) nrstr_dist
125 format ('No. distance rstrs =',i10)
if ( nrstr_dist .gt. 0 ) then
! allocate memory for rstdis
allocate(rstdis(nrstr_dist), stat=alloc_status)
call check_alloc('restraint list')
write (*,130)
end if
130 format ('atom_i atom_j   dist.      fc  istate') 
do i=1,nrstr_dist
! read rstdis(i)
read (fu,*) rstdis(i)%i,rstdis(i)%j,rstdis(i)%d1,rstdis(i)%fk, &
rstdis(i)%ipsi
rstdis(i)%d2 = rstdis(i)%d1 !no flat-bottom
write (*,132) rstdis(i)%i,rstdis(i)%j,rstdis(i)%d1,rstdis(i)%fk, &
rstdis(i)%ipsi
end do
132 format (i6,1x,i6,2f8.2,i8)
! --- nrstr_ang, [rstang]
read (fu,*) nrstr_angl
write (*,135) nrstr_angl
135 format ('No. angle rstrs =',i10)
if ( nrstr_angl .gt. 0 ) then
! allocate memory for rstang
allocate(rstang(nrstr_angl), stat=alloc_status)
call check_alloc('restraint list')
write (*,140)
end if
140 format ('atom_i atom_j atom_k   angle      fc  istate')
do i=1,nrstr_angl
! read rstang(i)
read (fu,*) rstang(i)%i,rstang(i)%j,rstang(i)%k,rstang(i)%ang, &
rstang(i)%fk,rstang(i)%ipsi
write (*,142) rstang(i)%i,rstang(i)%j,rstang(i)%k,rstang(i)%ang, &
rstang(i)%fk,rstang(i)%ipsi
end do
142 format (i6,1x,i6,1x,i6,2f8.2,i8)

! --- nrstr_wall, [rstwal]
read (fu,*) nrstr_wall
write (*,145) nrstr_wall
145 format ('No. wall seq. rstrs=',i10)
if ( nrstr_wall .gt. 0) then


! allocate memory for rstwal
allocate(rstwal(nrstr_wall), stat=alloc_status)
call check_alloc('restraint list')
write (*,150)
end if
150 format ('atom_i atom_j   dist.      fc  H-flag')
do i=1,nrstr_wall
! read rstwal(:)
read (fu,*) rstwal(i)%i,rstwal(i)%j,rstwal(i)%d,rstwal(i)%fk, &
rstwal(i)%ih
write (*,152) rstwal(i)%i,rstwal(i)%j,rstwal(i)%d,rstwal(i)%fk, &
rstwal(i)%ih     
end do
152 format (i6,1x,i6,2f8.2,i8)

read (fu,'(a80)') text
write (*,157) 
157 format ('Ignoring water file.')

! --- determine water model
read (fu,'(a80)') text
write (*,160) 
160 format ('Ignoring water model.')

end function old_initialize
!----------------------------------------------------------

subroutine temperature(tscale,Ekinmax)
! calculate the temperature
!arguments
TYPE(TGROUP_TYPE)                               :: tscale(:)
real(kind=prec)                                 :: Ekinmax

!locals
integer						::	i,tgroups
real(kind=prec)						::	Ekin

tscale(:)%temp  = zero
tscale(:)%tfree = zero
tscale(:)%texcl = zero
Temp  = zero
Tfree = zero
!get kinetic energies for atoms in each temperature group
do tgroups = 1, ntgroups
do i=tscale(tgroups)%starta,tscale(tgroups)%enda
Ekin = 0.5_prec*iaclib(iac(i))%mass*(qvec_square(v(i)))
tscale(tgroups)%temp = tscale(tgroups)%temp + Ekin

!******PWadded if
if( use_PBC .or. ( (.not. use_PBC) .and. (.not. excl(i)) ) ) then
        tscale(tgroups)%tfree = tscale(tgroups)%tfree + Ekin
else
        tscale(tgroups)%texcl = tscale(tgroups)%texcl  + Ekin
end if
if ( Ekin .gt. Ekinmax ) then
        ! hot atom warning
        write (*,180) i,2.0_prec*Ekin/Boltz/3.0_prec
end if
end do
end do
Tfree = sum(tscale(:)%tfree)
Temp  = sum(tscale(:)%temp)

E%kinetic = Temp

Temp  = 2.0_prec*Temp/Boltz/real(Ndegf, kind=prec)
Tfree = 2.0_prec*Tfree/Boltz/real(Ndegfree, kind=prec)

if (detail_temps) then
        do tgroups = 1, ntgroups
                tscale(tgroups)%temp  = 2.0_prec*tscale(tgroups)%temp &
                        /Boltz/real(tscale(tgroups)%Ndegf, kind=prec)
                tscale(tgroups)%tfree = 2.0_prec*tscale(tgroups)%tfree & 
                        /Boltz/real(tscale(tgroups)%Ndegfree, kind=prec)
                if ( tscale(tgroups)%Ndegf .ne. tscale(tgroups)%Ndegfree) &
                        tscale(tgroups)%texcl = 2.0_prec*tscale(tgroups)%texcl &
                        /Boltz/real(tscale(tgroups)%Ndegf - tscale(tgroups)%Ndegfree, kind=prec)
        end do
end if

if (thermostat == BERENDSEN) then
        do tgroups = 1, ntgroups
        if(tscale(tgroups)%tfree.ne.zero) tscale(tgroups)%sfact = &
                Temp0/tscale(tgroups)%tfree - one
        tscale(tgroups)%sfact = q_sqrt( one + dt/tau_T * tscale(tgroups)%sfact )
        end do
end if

180 format ('>>> WARNING: hot atom, i =',i10,' Temp(i)=',f10.2)

end subroutine temperature

!----------------------------------------------------------------------
! Subroutines for velocity update, one for each thermostat
subroutine newvel_ber(step,tempscale,startatom,endatom)!dt_mod,Tscale_solute,i,nat_solute)
real(kind=prec)                 :: step,tempscale
integer                         :: startatom, endatom
!locals
integer                         :: i

do i=startatom,endatom
v(i)  = (v(i) - d(i)*winv(i)*step) * tempscale
xx(i) = x(i)
x(i)  = x(i) + v(i) * dt
if (integrator == VELVERLET) then
        v(i) = (v(i) - d(i)*winv(i)*step) * tempscale
end if
end do

end subroutine newvel_ber

subroutine newvel_lan(step,lfriction,rands,rande)
real(kind=prec)                 :: step,lfriction,randnum
real(kind=prec)                 :: rands(:)
TYPE(qr_vec)                    :: rande(:)
!locals
integer                         :: i,j


if (.not.langevin_gauss_rand) then
!using uniform random numbers, default
do i=1,natom
randnum = randm(iseed)
rande(i)%x = rands(i) * ( randnum - 0.5_prec )
randnum = randm(iseed)
rande(i)%y = rands(i) * ( randnum - 0.5_prec )
randnum = randm(iseed)
rande(i)%z = rands(i) * ( randnum - 0.5_prec )
end do
else
!using gaussian random numbers, not default
do i=1,natom
call q_gauss(zero,rands(i),rande(i),iseed)
end do
end if  

do i=1,natom
v(i)  = (v(i)*lfriction) - (d(i)-rande(i))*winv(i)*step
xx(i) = x(i)
x(i)  = x(i) + v(i) * dt
if (integrator == VELVERLET) then
        v(i) = (v(i)*lfriction) - (d(i)-rande(i))*winv(i)*step
end if
end do


end subroutine newvel_lan

subroutine newvel_nos(step)
real(kind=prec)                 :: step

!locals
integer                         :: i,i3

do i=1,natom
v(i)  = (v(i) - d(i) * winv(i) *step)
xx(i) = x(i)
x(i)  = x(i) + v(i) * dt
if (integrator == VELVERLET) then
        v(i) = (v(i) - d(i) *winv(i) * step)
end if

end do

end subroutine newvel_nos

!-----------------------------------------------------------------------
!******PWchanged 2002-10-01
subroutine md_run

! local variables
integer				:: i,j,k,niter,iii
real(kind=prec)                         :: Tlast
real(kind=prec)                         ::Ekinmax
!new for temperature control
real(kind=prec)				::dv_mod,dv_friction,dv_mod2,randnum
integer				:: n_max = -1,tgroups
real(kind=prec),allocatable	::randva(:)
TYPE(qr_vec),allocatable        ::randfa(:)
!Random variable temperature control array
real(kind=prec)				:: time0, time1, time_per_step, startloop
integer(4)				:: time_completion

!Local variables for file writing to energy files
integer				:: loc_arrays
#if defined(PROFILING)
real(kind=prec)                                         :: start_loop_time1, start_loop_time2

profile(1)%name = 'NB_update'
profile(2)%name = '   nbwwlist_time'
profile(3)%name = '   nbpplist_time'
profile(4)%name = '   nbpwlist_time'
profile(5)%name = '   nbqplist_time'
profile(6)%name = '   nbqwlist_time'
profile(7)%name = 'SHAKE'
profile(8)%name = 'Bonded Terms'
profile(9)%name = 'Restraints'
profile(10)%name = 'Nonbonded Terms'
profile(11)%name = 'Update vel. & coords.'


#endif

#if defined(PROFILING)
#if defined(USE_MPI)
if (nodeid .eq. 0) then
allocate(all_node_times(num_profiling_times*numnodes), stat=alloc_status) !vector for storing all node's node_times, used by mpi_gather at end of md_run
call check_alloc('MPI profiling')
end if
allocate(node_times(num_profiling_times), stat=alloc_status) !each node's profiling times, used at end of md_run by mpi_gather
call check_alloc('MPI profiling')

all_node_times(:) = zero
node_times(:) = zero

#endif
#endif

!Define array length for EQ array in local usage
loc_arrays=ene_header%arrays

!Define number of coord to send/recieve
nat3=natom*3

! calculate maximum temperature
!**MN-> Only master calc. temp for now.
if (nodeid .eq. 0) then
Ekinmax = 1000.0_prec*Ndegf*Boltz*Temp0/2.0_prec/real(natom, kind=prec)

call temperature(tscale,Ekinmax)
!store old Temp
Tlast = Temp
end if

if (nodeid .eq. 0) then
! master node only: print initial temperatures
write (*,*)
write (*,120) 'Initial', Temp, Tfree
        if ( detail_temps ) then
                do tgroups = 1, ntgroups
                        write (*,120) tscale(tgroups)%tname, tscale(tgroups)%temp, tscale(tgroups)%tfree
                end do
        end if
120		format(a7,' temperatures are : Ttot =',f10.2,' Tfree =',f10.2)
write (*,*)

! init timer
time0 = rtime()

        ! Init timer of total loop time
        startloop = rtime()
end if

!set thermostat variables to the values needed in the time steps
!moved here to have definition before we do every step
if (nodeid.eq.0) then
!	allocate temp control arrays
if ( thermostat == LANGEVIN ) then
allocate(randva(natom))
call check_alloc('Random variable temperature control array')
allocate(randfa(natom))
call check_alloc('Random variable temperature control array 2')
end if
! if velocity verlet is used, take half steps
if ( integrator == VELVERLET ) then
dv_mod = dt2
else
dv_mod = dt
end if
if ( thermostat == LANGEVIN ) then
dv_friction = (1.0_prec-friction*dv_mod)
if (.not.langevin_gauss_rand) then
!we use the fast random number generation from unifrom random numbers
!this is the default! -- Paul and Alex October 2014
gkT = 24.0_prec*friction*Boltz*Temp0/dv_mod ! this 24 comes from the fact that we are using uniform random numbers
else
!we use gaussian random numbers
gkT = 2.0_prec*friction*Boltz*Temp0/dv_mod
end if
randva(:)= q_sqrt2 (gkT/winv(:))
end if
end if !nodeid .eq. 0

!***********************************************************************
!	begin MAIN DYNAMICS LOOP (Verlet leap-frog algorithm)
!***********************************************************************

! No loop (only calc. energies) if compiling with the DUM flag
#ifndef DUM
do istep = 0, nsteps-1
#endif
!change volume

if ( mod(istep, NBcycle) .eq. 0 ) then


        ! every NBcycle steps:

                        !Put molecules back in box for nice visualisation, needs to be here to prevent problems with LRF
                        !Update cgp_centers for LRF
                        !only call put_back_in_box if using PBC and either solute or solvent should be put back in box
                        if( use_PBC .and. (put_solute_back_in_box .or. put_solvent_back_in_box) ) then 
                          call put_back_in_box() 
                        end if
!Time estimate removed, can be activated again by passing -DTIME to the compiler
!Paul Bauer October 2014
#ifdef TIME
        if ((nodeid .eq. 0) .and. (istep > 0)) then
                ! print timing info
                call centered_heading('Timing', '-')
                time1 = rtime()
                time_per_step = 1000*(time1-time0)/NBcycle
                time_completion = int(time_per_step*(nsteps-istep)/60000)
                time0 = time1
                write(*,222) time_per_step, time_completion
222    format('Milliseconds per step (wall-clock): ',f5.2,&
                        ' Estimated completion in',i6,' minutes')
        end if
#endif
        ! update lists of nonbonded interaction pairs
        if (nodeid .eq. 0) then
           call centered_heading('Nonbonded pair list generation', '-')
        end if
        call make_pair_lists(Rcq,Rcq**2,RcLRF**2,Rcpp**2,Rcpw**2,Rcww**2)
#if defined(DUMP)
                        write(*,332) 'solute-solute', 'solute-water', 'water-water', 'Q-solute', 'Q-water'
                        write(*,333) nodeid, 'count', nbpp_pair, nbpw_pair, &
                                 &  nbww_pair, nbqp_pair, 3*nqat*nbqw_pair

#if defined(USE_MPI)
                        !reduce totnxx, i.e. collect # pairs found by slave nodes
                        nbxx(1)=nbpp_pair
                        nbxx(2)=nbpw_pair
                        nbxx(3)=nbww_pair
                        nbxx(4)=nbqp_pair
                        nbxx(5)=3*nqat*nbqw_pair

                        call MPI_Reduce(nbxx,nbxx_tot,5,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr) 
                        if (ierr .ne. 0) call die('run/Reduce')

                        if (nodeid .eq. 0) then
                        totnbpp=nbxx_tot(1)
                        totnbpw=nbxx_tot(2)
                        totnbww=nbxx_tot(3)
                        totnbqp=nbxx_tot(4)
                        totnbqw=nbxx_tot(5)
                        write(*,99) 'total', totnbpp,totnbpw,totnbww,totnbqp,totnbqw
                        end if
                        99  format(a10,1x,5(1x,i12))
#endif
                        332	format('node value ',5a13)
                        333	format(i4,1x,a5,1x,5(1x,i12))
#endif



        ! do MC_volume step here, after NBupdate for consistency with LRFs and such.
        ! Implies a limitation in the volume_update interval, i.e. mod(volume_update,nb_update)=0
                        if( use_PBC .and. constant_pressure) then
                           if( mod(istep, ivolume_cycle)==0 .and. istep>0 ) then
                                  call MC_volume
                           end if
                        end if

end if ! every NBcycle steps




! --- start of time step ---



! get potential energy and derivatives from FF
call pot_energy(E,EQ,.true.)
! if we have reached the ene write out, also calculate group contribution
! exclusions and qcp if needed
! exc are only done on master, while qcp is spread to nodes
if ( mod(istep, iene_cycle) == 0 .and. istep > 0) then
        if (use_excluded_groups) call calculate_exc
        if (use_qcp) call qcp_run(Tfree,E,EQ)
end if



if(nodeid .eq. 0) then
        if ( mod(istep,iout_cycle) == 0 .and. monitor_group_pairs > 0) then
           call nonbond_monitor  
        end if

        ! off-diagonals
        if ( noffd .gt. 0 ) call offdiag

#ifndef DUM
        ! update velocities from accelerations,
        ! scale velocities & update positions from velocities
#if defined (PROFILING)
start_loop_time1 = rtime()
#endif

! new selection and subfunction for thermostats
!if Berendsen thermostat was chosen, call two times for diff scaling
if( thermostat == BERENDSEN ) then
        do tgroups = 1, ntgroups
                call newvel_ber(dv_mod,tscale(tgroups)%sfact,tscale(tgroups)%starta,tscale(tgroups)%enda)
        end do
else if (thermostat == LANGEVIN ) then
	call newvel_lan(dv_mod,dv_friction,randva,randfa)
else if (thermostat == NOSEHOOVER ) then
	call nh_prop
	call newvel_nos(dv_mod)
end if
! --- end of thermostat section ---

#if defined (PROFILING)
profile(11)%time = profile(11)%time + rtime() - start_loop_time1
#endif

        ! shake if necessary
        if(shake_constraints > 0) then
                niter=shake(xx, x)
                v(:) = (x(:) - xx(:)) / dt
		if ( thermostat == NOSEHOOVER ) then
			call nh_prop !scaling velocities after applying SHAKE
		end if
        end if

        ! --- end of time step ---
#if defined (PROFILING)
start_loop_time2 = rtime()
#endif

        ! calculate temperature and scaling factor
        call temperature(tscale,Ekinmax)
#if defined (PROFILING)
profile(12)%time = profile(12)%time + rtime() - start_loop_time2
#endif

end if !if(nodeid .eq. 0)


#if defined(USE_MPI)
call MPI_Bcast(x, natom, mpitype_qrvec, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('MD Bcast x')
#endif

        ! print [intermediate] results (master node only)
        if (nodeid .eq. 0) then
                ! trajectory, energy data, output and backup restart file
                if ( mod(istep,itrj_cycle) == 0 .and. istep > 0) then
                ! write_trj: write x to the trajectory file
                        call write_trj
                end if
                ! energies
                if ( mod(istep, iene_cycle) == 0 .and. istep > 0) then
                ! nrgy_put_ene(unit, e2, OFFD): print 'e2'=EQ and OFFD to unit 'unit'=11
                ! save classical full system to EQ save, then write
                        call qatom_savetowrite(EQ,1)
                        call put_ene(11, EQ_save, OFFD,loc_arrays,nstates)
                end if
                ! end-of-line, then call write_out, which will print a report on E and EQ
                if ( mod(istep,iout_cycle) == 0 ) then
                        call write_out
                end if
                ! backup file of coordinates and velocities
                if ( mod(istep,1000) .eq. 0 ) then
                        call write_xfin
                end if
                if ( q_abs(Temp-Tlast)/Temp > TEMP_PRINT_THRESHOLD .or. &
                        (mod(istep, itemp_cycle) == 0 .and. istep > 0)) then
                        ! temperatures
                        Tlast = Temp
                        write(*,201) istep, Temp, Tfree
						if (detail_temps) then
                                                        do tgroups = 1, ntgroups
							        write(*,2020) trim(tscale(tgroups)%tname),tscale(tgroups)%tfree
                                                        end do
!							write(*,2030) Texcl_solute, Texcl_solvent

						end if

                        ! test for NaN
                        if (Temp.ne.Temp) then 
                            call die('a detected NaN.')
                        end if

                end if


        end if ! print results


end do ! time step
201	 format('Temperature at step',i8,':         T_tot=',f10.1,'         T_free=',f10.1)
2020 format('T_free_',a,'=',f10.1)


!***********************************************************************
!	end MAIN DYNAMICS LOOP
!***********************************************************************

! end of Qdum exclusion
#else
end if !(nodeid .eq. 0) from far above
#endif

        ! write final trajectory image when istep = nsteps
#ifndef DUM
if (nodeid .eq. 0) then
    if ( mod(istep,itrj_cycle) == 0) call write_trj
end if
#endif

        ! write output for final step and final coords
call make_pair_lists(Rcq,Rcq**2,RcLRF**2,Rcpp**2,Rcpw**2,Rcww**2)
call pot_energy(E,EQ,.true.)
if (nodeid .eq. 0) then
        write(*,*)
        call write_out
        call write_xfin
end if



if (nodeid .eq. 0) then
!Deallocate temperature control arrays
        if (allocated(randfa)) deallocate(randfa)
        if (allocated(randva)) deallocate(randva)

	time1 = rtime()
	write (*,202) time1 - startloop
	202 format('Total time of main loop:                     ', f15.1,'(s)')
end if
#if defined(PROFILING)
!Print more profiling info

#if defined(USE_MPI)
do i=1,num_profiling_times
	node_times(i) = profile(i)%time
end do
call MPI_GATHER(node_times,num_profiling_times,QMPI_REAL,all_node_times,num_profiling_times,QMPI_REAL,0,MPI_COMM_WORLD,ierr)
if (ierr .ne. 0) call die('md_run/MPI_GATHER profiling times')

if (nodeid .eq. 0) then
write (*,210,advance='no')
do j=0,numnodes-1
	write (*,209,advance='no') j
end do
write(*,*)

do i=1,num_profiling_times
	write (*,207,advance='no') profile(i)%name
	do j=0,numnodes-1
		write (*,208,advance='no') all_node_times(i+j*num_profiling_times)
	end do
	write (*,*) ' (s)'
end do

207	format('Total time of ',A25,T40,': ')
208	format(f10.1,' ')
209	format(I11)
210 format(T30,'node:     ')

write (*,*)

end if

#else

do i=1,num_profiling_times
	write (*,207) profile(i)%name,profile(i)%time
end do
207	format('Total time of ',A25,T40,': ',f15.1,' (s).')
#endif
#endif


end subroutine md_run

!-----------------------------------------------------------------------

subroutine MC_volume()

TYPE(qr_vec)                    :: old_x(nat_pro), old_xx(nat_pro), x_move(nat_pro)
TYPE(qr_vec)                    :: old_boxl, old_inv, cm, old_dMorse_i(max_qat), old_dMorse_j(max_qat)
TYPE(qr_vec)                    :: newvec, movevec

real(kind=prec)									:: old_V, new_V, deltaLength
real(kind=prec)									:: deltaV, deltaE, deltaW
real(kind=prec)									:: box_min
integer									:: starten, slutet, i, j, sw_no !indeces
real(kind=prec)									:: randomno !random number
logical									:: acc
integer									:: longest , niter
real(kind=prec)									:: cubr_vol_ratio
real(kind=prec)									::	old_EMorseD(max_qat)
integer									:: old_nbww_pair, old_nbpw_pair , old_nbpw_cgp_pair 
integer									:: old_nbpp_pair , old_nbpp_cgp_pair , old_nbqp_pair , old_nbqp_cgp_pair
integer									:: old_ww_max, old_pw_max, old_pp_max, old_qp_max

if (nodeid .eq. 0) then
write(*,8) 'Volume change', istep
write(*,*)
write(*,'(a)') '---------- Before move'
8 format('======================== ',A14,' at step ',i6,' ========================')
4 format(16X, 3A10)
6 format(A,T17, 3F10.3)
end if  !(nodeid .eq. 0)

!save the old energies,coordinates and forces
old_x(:) = x(:)
previous_E = E
old_EQ = EQ
old_boxl = boxlength
old_inv = inv_boxl

!qatom stuff
if(nqat.gt.0) then
        old_EMorseD = EMorseD
        old_dMorse_i = dMorse_i
        old_dMorse_j = dMorse_j
end if
!shake stuff
old_xx(:) = xx(:)



if (use_LRF) then
	old_nbww_pair = nbww_pair
	old_nbpw_pair = nbpw_pair
	old_nbpp_pair = nbpp_pair
	old_nbqp_pair = nbqp_pair

	old_nbpw_cgp_pair = nbpw_cgp_pair
	old_nbpp_cgp_pair = nbpp_cgp_pair
	old_nbqp_cgp_pair = nbqp_cgp_pair
!now only reallocate if the old size was too small
	if(calculation_assignment%ww%max.gt.SIZE(old_nbww,1)) then
		if (allocated(old_nbww)) deallocate(old_nbww)
		allocate(old_nbww(calculation_assignment%ww%max))
	endif
        if(calculation_assignment%pw%max.gt.SIZE(old_nbpw,1)) then
                if (allocated(old_nbpw)) deallocate(old_nbpw)
                allocate(old_nbpw(calculation_assignment%pw%max))
        endif
        if(calculation_assignment%pp%max.gt.SIZE(old_nbpp,1)) then
                if (allocated(old_nbpp)) deallocate(old_nbpp)
                allocate(old_nbpp(calculation_assignment%pp%max))
        endif
        if(calculation_assignment%qp%max.gt.SIZE(old_nbqp,1)) then
                if (allocated(old_nbqp)) deallocate(old_nbqp)
                allocate(old_nbqp(calculation_assignment%qp%max,nstates))
        endif
	old_nbww = nbww
	old_nbpw = nbpw
	old_nbpp = nbpp
	old_nbqp = nbqp
	old_nbpw_cgp = nbpw_cgp
	old_nbpp_cgp = nbpp_cgp
	old_nbqp_cgp = nbqp_cgp
	old_ww_max = calculation_assignment%ww%max
	old_pw_max = calculation_assignment%pw%max
	old_pp_max = calculation_assignment%pp%max
	old_qp_max = calculation_assignment%qp%max
	old_lrf(:) = lrf(:)
end if

#if defined(USE_MPI)
!Update modified coordinates  
call MPI_Bcast(x, natom, mpitype_qrvec, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast x')
#endif

call pot_energy(E,EQ,.true.)

if (nodeid .eq. 0 ) then
	old_E = E                !Update to fresh E before changing volume
        if (nqat.gt.0)	old_EQ(1:nstates) = EQ(1:nstates)
	old_V = old_boxl%x*old_boxl%y*old_boxl%z


	!new volume randomized
	randomno = randm(pressure_seed) ! 0<=randomno<=1
	randomno = randomno*2 - 1    !-1 <= randomno <= 1
	deltaV = randomno * max_vol_displ
	new_V = deltaV + old_V
	cubr_vol_ratio = (new_V/old_V)**(one/3.0_prec)
	write(*,4) 'old', 'new', 'delta'
17	format('Volume',8x,f10.3,2x,f10.3,2x,f10.3)
	write(*,17) old_V, new_V, deltaV
	write(*,*)

	!compute new boxlenth and inv_boxl
        boxlength = boxlength * cubr_vol_ratio
	inv_boxl  = one/boxlength
	write(*,10) old_boxl
	write(*,2) boxlength
	write(*,*)
	10 format('Old boxlength', 3f10.3)
	2 format('New boxlength ', 3f10.3)

	!compare cut-offs with new boxsize
	box_min = min( boxlength%x, boxlength%y, boxlength%z )
	!Solute-solute cut-off radii
	if( .not. (box_min .gt. Rcpp*2) ) then
		write(*,*) 'Solute-solute cut-off radii too large', Rcpp
		call die('Solute-solute cut-off radii too large')
	!Solvent-solvent
	else if( .not. (box_min .gt. Rcww*2) ) then
		write(*,*) 'Solvent-solvent cut-off radii too large', Rcww
		call die('Solvent-solvent cut-off radii too large')
	!Solute-solvent
	else if( .not. (box_min .gt. Rcpw*2) ) then
		write(*,*) 'Solute-solvent cut-off radii too large', Rcpw
		call die('Solute-solvent cut-off radii too large')
	!Q-atom
	else if( .not. (box_min .gt. Rcq*2) ) then
		write(*,*) 'Q-atom cut-off radii too large', Rcq
		call die('Q-atom cut-off radii too large')
	!LRF
	else if( .not. (box_min .gt. RcLRF*2) ) then
		write(*,*) 'LRF cut-off radii too large', Rcq
		call die('LRF cut-off radii too large')
	end if




	if (.not. atom_based_scaling) then


		!compute new coordinates after molecules and centre of mass
		do i=1,nmol !looping over all molecules, can also do the last because of nmol +1 ...
			cm = cm * zero
			do j = istart_mol(i),istart_mol(i+1)-1 !loop over all atoms in molecule i
                                cm = cm + x(j) * mass(j)
			end do
!but mol_mass is actually 1/mol_mass
			cm = cm * mol_mass(i) !centre of mass of molecule i

                        movevec = ((cm - boxcentre) * boxlength/old_boxl + boxcentre) - cm

			do j= istart_mol(i),istart_mol(i+1)-1 !move the molecule
                                x(j) = x(j) + movevec
			end do

		end do !over molecules

!we have nmol+1 to be able to iterate over all molecules

	else ! atom_based_scaling = .true.
	!move xx also if shake
		do j=1,natom
                        x_move(j) = (( x(j) - boxcentre) * boxlength/old_boxl + boxcentre) - x(j)
		end do

		x(:) = x(:) + x_move(:)
		
		! shake if necessary
		if(shake_constraints > 0) then
			xx(:) = xx(:) + x_move(:)
			niter=shake(xx, x)
		end if

	end if
end if

#if defined(USE_MPI)
!Update modified coordinates and boxlengths 
call MPI_Bcast(x, natom, mpitype_qrvec, 0, MPI_COMM_WORLD, ierr)
if (ierr .ne. 0) call die('init_nodes/MPI_Bcast x')
call MPI_Bcast(boxlength, 1, mpitype_qrvec, 0, MPI_COMM_WORLD, ierr)
call MPI_Bcast(inv_boxl, 1, mpitype_qrvec, 0, MPI_COMM_WORLD, ierr)
#endif

!Need to update entire LRF... sigh
if (use_LRF) then
        call make_pair_lists(Rcq,Rcq**2,RcLRF**2,Rcpp*2,Rcpw**2,Rcww**2)
end if !use_LRF



!compute the new potential, in parallel if possible
call pot_energy(E,EQ,.true.)

if (nodeid .eq. 0) then
	!Jamfor nya med gamla
	deltaE = E%potential - old_E%potential
	deltaW = deltaE + pressure * deltaV - nmol*Boltz*Temp0*log(new_V/old_V)
	write(*,4) 'old', 'new', 'delta'
18	format('Potential',7x,f14.3,2x,f14.3,2x,f14.3)
	write(*,18) old_E%potential, E%potential, deltaE
	write(*,*)

	!accept or reject
	if( deltaW<=zero ) then
	acc = .true.
	else
	!slumpa tal mellan  0 coh 1
	randomno = randm(pressure_seed)
	if( randomno > exp(- deltaW / Boltz / Temp0) ) then
		acc = .false.
	else
		acc = .true.
	end if
	end if

	volume_try = volume_try + 1
	write(*,'(a)') '---------- After move' 

if( acc ) then
		write(*,'(a)') 'Volume change accepted'
		volume_acc = volume_acc + 1
else
		write(*,'(a)') 'Volume change rejected'

        !put stuff back to what they were before
        x(:) = old_x(:)
        E = previous_E
        EQ(1:nstates) = old_EQ(1:nstates)
        boxlength = old_boxl
        inv_boxl  = old_inv
        if(nqat.gt.0) then
                EQ(1:nstates) = old_EQ(1:nstates)
                EMorseD = old_EMorseD
        	dMorse_i = old_dMorse_i
        	dMorse_j = old_dMorse_j
        end if
	xx(:) = old_xx(:)

	if (use_LRF) then
		nbww_pair = old_nbww_pair
		nbpw_pair = old_nbpw_pair
		nbpp_pair = old_nbpp_pair
		nbqp_pair = old_nbqp_pair

		nbpw_cgp_pair = old_nbpw_cgp_pair
		nbpp_cgp_pair = old_nbpp_cgp_pair
		nbqp_cgp_pair = old_nbqp_cgp_pair
		calculation_assignment%ww%max = old_ww_max
		calculation_assignment%pw%max = old_pw_max
		calculation_assignment%pp%max = old_pp_max
		calculation_assignment%qp%max = old_qp_max
        if(calculation_assignment%ww%max.gt.SIZE(nbww,1)) then
                call reallocate_nonbondlist_ww
        endif
        if(calculation_assignment%pw%max.gt.SIZE(nbpw,1)) then
                call reallocate_nonbondlist_pw
        endif
        if(calculation_assignment%pp%max.gt.SIZE(nbpp,1)) then
                call reallocate_nonbondlist_pp
        endif
        if(calculation_assignment%qp%max.gt.SIZE(nbqp,1)) then
                call reallocate_nonbondlist_qp
        endif
		nbww(1:old_ww_max) = old_nbww(1:old_ww_max)
		nbpw(1:old_pw_max) = old_nbpw(1:old_pw_max)
		nbpp(1:old_pp_max) = old_nbpp(1:old_pp_max)
		nbqp(1:old_qp_max,1:nstates) = old_nbqp(1:old_qp_max,1:nstates)
		nbpw_cgp(:) = old_nbpw_cgp(:)
		nbpp_cgp(:) = old_nbpp_cgp(:)
		nbqp_cgp(:) = old_nbqp_cgp(:)
		lrf(:) = old_lrf(:)
	end if
end if

	write(*,11) boxlength%x*boxlength%y*boxlength%z
	write(*,12) boxlength
	write(*,13) sum(mass(:))/(boxlength%x*boxlength%y*boxlength%z*1E-24_prec*6.02E23_prec)
	11 format('Final volume: ', f10.3)
	12 format('Final boxlength: ', 3f10.3)
	13 format('Final density (g/cm^3): ', f10.3)
	write(*,*)
	write(*,4) 'total', 'accepted', 'ratio'
	write(*,7) 'Attempts', volume_try, volume_acc, real(volume_acc, kind=prec)/volume_try
	write(*,*)
	7 format(A,T17, 2i10, f10.3)
	write(*,'(80a)') '==============================================================================='
end if !(nodeid .eq. 0)
#if defined(USE_MPI)
!Make slave nodes put things back if rejected
call MPI_Bcast(acc, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
if (.not. acc) then
  if (nodeid .ne. 0) then
    x(:) = old_x(:)
    boxlength = old_boxl
    inv_boxl  = old_inv
  end if
end if
#endif
end subroutine MC_volume	

!----------------------------------------------------------------------------------------

end module MD




