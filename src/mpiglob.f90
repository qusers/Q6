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
! mpiglob.f90
! by John Marelius & Anders Kaplan
! global variables for MPI parallell Qdyn

module	MPIGLOB

use NRGY
!$ use omp_lib
#ifdef USE_MPI
include "mpif.h"
#endif

    ! types and variables used for calculation assignment
    type PAIR_ASSIGNMENT_TYPE
            integer					:: start, end
            integer					:: max
    end type PAIR_ASSIGNMENT_TYPE


    type NODE_ASSIGNMENT_TYPE
            type(PAIR_ASSIGNMENT_TYPE) :: pp, pw, ww, qp, qw, shake, &
                                          natom, nat_solute, nat_solvent, &
                                          nmol, ncgp
    end type NODE_ASSIGNMENT_TYPE

 !Gathering
  type MPI_NB_ENERGIES
     sequence
     real(kind=prec)       :: lrf
     type(NB_ENERGIES)    :: pp,pw,ww
  end type MPI_NB_ENERGIES

 !Gathering
  type MPI_NBQ_ENERGIES
     sequence
     type(NB_ENERGIES)    :: qp,qw
  end type MPI_NBQ_ENERGIES

  ! global MPI data
  integer                  :: nodeid, numnodes,ierr
  integer,allocatable      :: mpi_status(:,:)

  !Used for gathering d,E,EQ
  integer			:: reclength=-1
  TYPE(qr_vec),allocatable      :: d_recv(:,:)
  type(MPI_NB_ENERGIES),allocatable  :: E_recv(:),E_send(:)
!additional dimension for new structure of EQ arrays
  type(MPI_NBQ_ENERGIES),allocatable :: EQ_recv(:,:),EQ_send(:)
  integer,allocatable            ::request_recv(:,:)

 !Book keeping of nb-pairs
 integer  :: totnbpp,totnbpw,totnbww,totnbqp,totnbqw
 integer,dimension(5)   ::nbxx,nbxx_tot


 !Balancing nodes
 !Keep track of how many nb-pairs each chargegroup will generate
 integer,allocatable  :: nbpp_per_cgp(:)
 integer,allocatable  :: nbpw_per_cgp(:)
 integer,allocatable  :: nbww_per_cgp(:)
 integer,allocatable  :: nbqp_per_cgp(:)
 integer,allocatable  :: nbqw_per_cgp(:)

! custom Q types to convert MPI_REAL to actual size of Q variables
! will  crash in beautiful ways if not assigned
 integer :: QMPI_REAL = -1
 integer :: QMPI_VECT = -1
 integer :: QMPI_BYTE = -1
 ! stuff for mpi type create
 integer                :: mpitype_batch_lrf,mpi_lrf_add,mpi_lrf_cgp_rep,mpitype_batch_solv_int
 integer                :: mpitype_batch_ppgrid,mpitype_batch_pwgrid,mpitype_batch_wwgrid,mpi_grid_add
 integer                :: mpitype_qrvec,mpitype_batch_cgp,mpitype_batch_iaclib,mpitype_batch_grid
 integer                :: mpitype_pair_assignment,mpitype_node_assignment

#ifdef _OPENMP
 !for omp use
 integer  :: thread_id,threads_num,mp_start,mp_end,mp_counter
 real(kind=prec)  :: mp_real_tmp
 real(kind=prec),allocatable  :: qomp_elec(:),qomp_vdw(:)

 !$omp threadprivate(thread_id, mp_real_tmp, mp_start, mp_end, mp_counter)
#endif

end module MPIGLOB

