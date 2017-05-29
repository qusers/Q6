! (C) 2014 Uppsala Molekylmekaniska HB, Uppsala, Sweden
! minim.f90
! by Paul Bauer, based on
! md.f90
! by Johan Åqvist, John Marelius, Anders Kaplan, Isabella Feierberg, Martin Nervall & Martin Almlöf
! energy minimization

module MINIM

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


subroutine minim_startup
! initialise used modules
call simprep_startup

! read in version info
MD_VERSION = trim(version_pass())
MD_DATE    = trim(date_pass())


end subroutine minim_startup


!----------------------------------------------------------------------

subroutine minim_shutdown
! call used modules' shutdown subroutines
if(use_qcp) call qcp_shutdown
call simprep_shutdown
call md_deallocate
call topo_deallocate
end subroutine minim_shutdown

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
real(kind=prec)                                 :: size_default = one/10000
real(kind=prec)                                 :: Ecut_default = 1.0_prec
logical                                         :: shake_all,shake_all_solvent,shake_all_solute
logical                                         :: shake_all_hydrogens,shake_all_heavy
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

if(.not. prm_open_section('minim',infilename)) then
call prm_close
initialize = .false.
else

need_restart = .false. !flag for restart file required
if(.not. prm_get_integer_by_key('steps', nsteps)) then
write(*,*) '>>> ERROR: steps not specified (section MINIM)'
initialize = .false.
end if
yes = prm_get_real_by_key('initialsize', stepsize,size_default)
write (*,10) nsteps, stepsize
10	format ('Number of Minimization steps =',i10,'  Initial Stepsize    =',f10.3)
dt = stepsize
itdis = one
! update polarization restraints at every step
if(.not. prm_get_real_by_key('Ecut',energy_cutoff)) then
        energy_cutoff = Ecut_default
        write(*,'(a,f8.3)') 'Energy end cutoff set to default, ',energy_cutoff
else if(energy_cutoff .gt.zero) then
        write(*,'(a,f8.3)') 'Energy end cutoff set to ',energy_cutoff
else
        write(*,*) 'No energy end cutoff'
end if

if(.not. prm_get_string_by_key('integrator', name_integrator)) then
        write(*,*) 'Steepest decent used by default.'
        integrator = MIN_STEEP

else if( name_integrator == 'conjugate-gradient' ) then
        write(*,*) 'Integrator: ', name_integrator
        integrator = MIN_CG
else if( name_integrator /= 'steepest-decent' .and. name_integrator /= 'conjugate-gradient' ) then
        write(*,*) '>>> ERROR: no such minimizer exists in Q.',name_integrator
        initialize = .false.
else
        name_integrator = 'steepest-decent'
        integrator = MIN_STEEP
        write(*,*) 'Minimizer: ', name_integrator
end if

if(.not. prm_get_logical_by_key('shake_solvent', shake_solvent, .false.)) then
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

if(.not. prm_get_logical_by_key('shake_hydrogens', shake_hydrogens, .false.)) then
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
constant_pressure = .false.
control_box = .false.

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
        fkwpol = -one ! this will be set in water_sphere, once target radius is known
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
NBcycle = 1 ! update every step during minim
write(*,'(a)') 'energy summary interval set to default.'
iout_cycle = 1 ! write out every step
itemp_cycle = 0 !no temperature in minim
iene_cycle = 0 !no energy
itrj_cycle = 0 !no trajectory

ivolume_cycle = 0 !no volume change in minimization


else
if(.not. prm_get_integer_by_key('non_bond', NBcycle)) then
        write(*,'(a)') 'non-bond list update interval set to default.'
        NBcycle = 1
end if
if(.not. prm_get_integer_by_key('output', iout_cycle)) then
        write(*,'(a)') 'energy summary interval set to default.'
        iout_cycle = 1
end if
        itemp_cycle = 0
yes = prm_get_integer_by_key('energy', iene_cycle, 0)
yes = prm_get_integer_by_key('trajectory', itrj_cycle, 0)
        ivolume_cycle = 0
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
else
        restart = .true.
end if

if(restart) then
        write (*,65) trim(restart_file)
else
        write (*,'(a)') 'Initial coordinates taken from topology.'
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
!new section in *inp files to control QCP behaviour
!will only trigger if FEP file is in use -> if more than 0 states
        use_qcp = .false.
        QCP_N = QCP_OFF
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


! we need some dummy stuff here to handle the temperature
! just using one group and default settings
ntgroups = 1
ntgroups_kind = DEFAULT_ONE
allocate(tscale(ntgroups))

call prm_close
end function initialize


!-----------------------------------------------------------------------------------

integer function maxforce(new_d,oldmax)
! calculate the current maximum force 
! arguments
TYPE(qr_vec)                    :: new_d(:)
real(kind=prec),INTENT(INOUT)   :: oldmax
!locals
real(kind=prec)                 :: maxf
integer                         :: i, atmaxf
real(kind=prec)                 :: f2old,f2new

! check if force on atom is increasing 
maxforce   = 0
maxf       = zero
atmaxf     = -1
do i = 1, nat_pro
        f2new = q_sqrt(qvec_square(new_d(i)))
        if (f2new.gt.maxf) then
                maxf   = f2new
                atmaxf = i
        end if
        if (f2new.gt.oldmax) then
                oldmax = f2new
                maxforce = 1
        end if
end do
! if forces are lower, set oldforce to lowest limit
if(maxforce .eq. 0) oldmax = maxf

if(maxf .lt. energy_cutoff) maxforce = 2

write(*,'(a,f12.3,a,i8)') 'Maximum force of ',maxf,' on atom ',atmaxf

end function maxforce

!----------------------------------------------------------------------
subroutine update_pos(step)
real(kind=prec)                 :: step
!locals
integer                         :: i

do i=1,nat_pro
x(i)  = x(i) - d(i) * step
end do

end subroutine update_pos

!-----------------------------------------------------------------------
subroutine min_run

! local variables
integer				:: i,j,k,niter,retval
!Random variable temperature control array
real(kind=prec)				:: scaling,time0, time1, time_per_step, startloop, force
integer(4)				:: time_completion
!Local variables for file writing to energy files
integer				:: loc_arrays
!Define array length for EQ array in local usage
loc_arrays=ene_header%arrays

!Define number of coord to send/recieve
nat3=natom*3
if (nodeid .eq. 0) then
! master node only: print initial temperatures
! init timer
time0 = rtime()

        ! Init timer of total loop time
        startloop = rtime()
end if

! calculate initial forces etc before main loop

force  = zero
retval = maxforce(d,force)

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

end if ! every NBcycle steps




! --- start of time step ---



! get potential energy and derivatives from FF
call pot_energy(E,EQ,.true.)

if(nodeid .eq. 0) then
        if ( mod(istep,iout_cycle) == 0 .and. monitor_group_pairs > 0) then
           call nonbond_monitor  
        end if

        ! off-diagonals
        if ( noffd .gt. 0 ) call offdiag

#ifndef DUM
        ! update velocities from accelerations,
        ! scale velocities & update positions from velocities
! section for checking max force and position update
retval = maxforce(d,force)
if(retval.eq.1) then
! we have a rise in energy, reset scaling factor for steepest decent
        scaling = dt
        call update_pos(scaling)
else if(retval .eq. 0) then
! increase step size
        scaling = scaling*1.20
        call update_pos(scaling)
end if

        if(constraints .gt.0) niter=constraint(const_method,xx,x)

end if !if(nodeid .eq. 0)

#ifdef USE_MPI
call MPI_Bcast(retval,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if (ierr.ne.0) call die('MINIM Bcast retval')
#endif

if (retval .eq. 2) then
! we have decrease in force and low energy -> stop
if (nodeid.eq.0) write (*,*) 'System has converged'
        exit
end if

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
                        ! test for NaN
                        if (scaling.ne.scaling) then 
                            call die('a detected NaN.')
                    end if

        end if ! print results


end do ! time step

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
        ! set velocities to zero now
        ! TODO
        ! need flag in restart to not allow MD without new velocities
        v(:) = v(:) * zero
        call write_out
        call write_xfin
end if



if (nodeid .eq. 0) then
	time1 = rtime()
	write (*,202) time1 - startloop
	202 format('Total time of main loop:                     ', f15.1,'(s)')
end if

end subroutine min_run

!-----------------------------------------------------------------------

end module MINIM




