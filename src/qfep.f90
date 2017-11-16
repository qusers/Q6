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
! qfep.f90
! by Johan Åqvist, Karin Kolmodin & John Marelius
! Qfep free energy analysis program for FEP, EVB & umbrella sampling

program qfep
use VERSIONS
use NRGY
use PARSE
use QMATH
implicit none
        character(80)   ::      QFEP_VERSION = ''
        character(80)   ::      QFEP_DATE    = ''
        character(10)   ::      QFEP_NAME    = 'qfep'
        character(10)   ::      QFEP_SUFFIX  = ''
	integer		:: mxpts=200,mxbin=10,mxstates=4
	character(80)      ::filnam, line
	integer           ::i,j,ifile,ipt,istate,ibin,nfiles,nstates,ERR, &
	                      nskip,nbins,nmin,idum,noffd,nnoffd,offel,num_offd

	type(OFFDIAG_SAVE),allocatable		:: offd(:)

        real(kind=prec) ::rt,gaprange, &
                        xint,dvg,dGg,scale_Hij, &
                        veff,min
	real(kind=prec),allocatable		:: sumg(:),sumg2(:),avdvg(:),avc1(:),avc2(:),avr(:)
	real(kind=prec),allocatable		:: binsum(:,:)

	integer,allocatable			:: nbinpts(:),ptsum(:)

!Change for new save data type with different array sizes
!	type(OQ_ENERGIES), dimension(mxstates)	:: old_EQ
!	type(OQ_ENERGIES), dimension(mxstates)	:: old_avEQ

	type(Q_ENERGIES),allocatable	:: EQ(:)
	type(Q_ENERGIES),allocatable	:: avEQ(:)
	real(kind=prec),allocatable	:: dvv(:),dGv(:),alfa(:),coeff(:)



	real(kind=prec),allocatable             :: dgf(:,:), dgr(:,:), dgfsum(:,:), dgrsum(:,:), dG(:,:)
	real(kind=prec),allocatable		:: gapmin(:), gapmax(:), sum(:), dv(:), veff1(:), veff2(:)
	real(kind=prec),allocatable		:: A(:,:), mu(:,:), eta(:,:), rxy0(:,:)
	real(kind=prec),allocatable		:: Hij(:,:,:), d(:), e(:)
	type FEP_DATA_TYPE
		integer					::	npts
		real(kind=prec), pointer	:: lambda(:)
		real(kind=prec), pointer	:: v(:,:,:), r(:,:,:) !indices are calulation type, state, point
		real(kind=prec), pointer	:: vg(:,:), gap(:,:), c1(:,:), c2(:,:) !indices are type, point
	end type FEP_DATA_TYPE
	type(FEP_DATA_TYPE), allocatable	::	FEP(:) !index is file
	type(FEP_DATA_TYPE)			::	FEPtmp !temporary storage for one file
							       !now with different types for later

	real(kind=prec), allocatable		:: avdvv(:,:), sumv(:,:), sumv2(:,:)
!Energy file header data type
	type(Q_ENE_HEAD)			:: fileheader
!to keep current residue from header
	integer					:: curres,filestat,jj,canary,myprec
	logical					:: is_old_file = .false.
	character(80)				:: version_precision, file_precision
	integer								::	f,gas=0,error,dummyno !!!!!!!!!masoud
	real(kind=prec)							::	dummy
	character(100)							::	iline !!!!!!! masoud
! reading in version information
        QFEP_VERSION = trim(version_pass())
        QFEP_DATE    = trim(date_pass())


	!header
        call version_check(QFEP_NAME,QFEP_VERSION,QFEP_DATE,QFEP_SUFFIX)
!	write(*,100) QFEP_VERSION,  QFEP_DATE
!	write(*,*)
!100	format('# Qfep',t30,'version ',a,t50,'(modified on ',a,')')

! checking and storing precision
if (prec .eq. singleprecision) then
myprec = 137
version_precision = 'Single'
elseif (prec .eq. doubleprecision) then
myprec = 1337
version_precision = 'Double'
#ifdef HAVEQUAD
elseif (prec .eq. quadprecision) then
myprec = 13337
version_precision = 'Quadruple'
#endif
end if

	!------------------------------------------
	! INPUT OF PARAMETERS
	!------------------------------------------
	call prompt ('--> Number of energy files: ')
	read (*,*) nfiles
	write (*,1) nfiles
1	format('# Number of files                 =',i6)
	call prompt ('--> No. of states, no. of predefined off-diag elements: ')
	read (*,*) nstates, noffd
	write (*,2) nstates, noffd
2	format('# Number of states                 =',i6,/, &
		   '# Number of off-diagonal elements =',i6)
	mxstates = nstates

	!size of secular determinant is allocated

	! Continue to read input
	call prompt ('--> Give kT & no, of pts to skip   & calculation mode: ')
	read(*,'(a)') iline

	!!!!!!!!!!!!!!  masoud reading an extra option for just calculating QQ energies
	!print*, iline
	do i=1,4
		read (iline,*,iostat=error)(dummy,j=1,i)
!		print*, error,"***",i,dummy
		if (error .ne. 0) then
			dummyno=i-1
			exit
		endif
	enddo
	if (dummyno .eq. 2) then
		read (iline,*) rt,nskip
	elseif (dummyno .eq. 3) then
		read (iline,*) rt,nskip,gas
	else
		print*, "not correct argumets for KT, data points to skip or calculation mode"
		stop
	end if
	!!!!!!!!!!!!!!  masoud

!	read (*,*) rt,nskip
	write (*,3) rt,nskip,gas
3	format('# kT                              =',f6.3,/, &
		   '# Number of data points to skip   =',i6,/, '# Only QQ interactions will be considered   = ',i3)

	call prompt ('--> Give number of gap-bins: ')
	read (*,*) nbins
	write (*,4) nbins
4	format('# Number of gap-bins              =',i6)

	mxbin = nbins + 1
! now know the number of bins needed -> allocate dependent variables + 1
! and we have all nstates dependent variables -> allocate nstate arrays
	allocate(sumg(mxbin),sumg2(mxbin),avdvg(mxbin),avc1(mxbin),avc2(mxbin),avr(mxbin), &
		binsum(mxbin,4),nbinpts(mxbin),ptsum(mxbin),STAT=ERR)
	if(ERR /= 0) then
		stop '# >>> ERROR! Qfep5 terminated abnormally: Out of memory when allocating mxbin arrays'
	end if
	
	allocate(EQ(mxstates),avEQ(mxstates),dvv(mxstates),dGv(mxstates),alfa(mxstates),coeff(mxstates),&
		A(mxstates,mxstates), mu(mxstates,mxstates), eta(mxstates,mxstates), rxy0(mxstates,mxstates),&
		avdvv(mxstates,mxbin), sumv(mxstates,mxbin), sumv2(mxstates,mxbin),offd(mxstates), STAT=ERR)
	dvv   = zero
	dGv   = zero
	alfa  = zero
	coeff = zero
	A     = zero
	mu    = zero
	eta   = zero
	rxy0  = zero
	avdvv = zero
	sumv  = zero
	sumv2 = zero
	
	if(ERR /= 0) then
		stop '# >>> ERROR! Qfep5 terminated abnormally: Out of memory when allocating mxstates arrays'
	end if


	call prompt ('--> Give minimum # pts/bin: ')
	read (*,*) nmin
	write (*,5) nmin
5	format('# Minimum number of points per bin=',i6)

	do istate=2,nstates
		write(line,7) istate
		call prompt(line)
		read (*,*) alfa(istate)
		write (*,6) istate,alfa(istate)
	end do
6	format('# Alpha for state ',i2,'              =',f9.2)
7	format('--> Give alpha for state ',i2,':')

	scale_Hij=0.0
	if (noffd /=0) then
		call prompt ('--> Hij scaling:')
		read (*,*) scale_Hij
		write (*,8) scale_Hij
8	format('# Scale factor for Hij            =',f6.2)
! variable to store the actual offdiagonals used later
                num_offd = noffd
	else
		call prompt ('--> Number of off-diagonal elements:')
		read (*,*) nnoffd
		write (*,9) nnoffd
		if(nnoffd >0) write(*,11)
		do offel=1,nnoffd
			call prompt ('--> i, j, A_ij, mu_ij, eta_ij, r_xy0: ')
			read (*,*) i, j, A(i,j), mu(i,j), eta(i,j), rxy0(i,j)
			write(*,12) i, j, A(i,j), mu(i,j), eta(i,j), rxy0(i,j)
		end do
! variable to store the actual offdiagonals used later
                num_offd = nnoffd
9	format('# Number of off-diagonal elements =',i6)
11	format('#   i   j   A(i,j)  mu(i,j) eta(i,j) rxy0(i,j)')
12	format('#',2i4,1x,f9.2,1x,f9.2,1x,f9.2,1x,f9.2)
	end if

	call prompt ('--> linear combination of states defining reaction coord: ')
	read (*,*) (coeff(i),i=1,nstates)
	write(*,13) coeff(1:nstates)
13	format('# Linear combination co-efficients=',8f6.2)


!---------------------------------
! Energy files are opened and read
!---------------------------------
	binsum(:,:)=zero
	f = freefile()

	write(*,*)
	write(*,*)
	write(*,15)
	write(*,16)
15	format('# Part 0: Average energies for all states in all files')
16	format('# file             state   pts   lambda    EQtot   EQbond',&
	'  EQang   EQtor   EQimp    EQel   EQvdW  Eel_qq  EvdW_qq Eel_qp  EvdW_qp Eel_qw EvdW_qw Eqrstr')

	do ifile=1,nfiles
		write(line,14) ifile
14		format('--> Name of file number',i4,':')
		call prompt(line)
		read(*,*) filnam
		write (*,*) ''
		if(openit(f,filnam,'old','unformatted','read') /= 0) then
			stop 'Qfep5 terminated abnormally: Failed to open energy file.'
		end if
                ! reset check for old/new file, if people mix them
                ! bad idea in general, but works with same number of arrays
                is_old_file = .false.
		!get file header for the first time
		read(f, iostat=filestat) canary,fileheader%arrays,fileheader%totresid
		if ((canary .ne. 137).and.(canary.ne.1337).and.(canary.ne.13337)) then
!could be old file type, save the info and continue
!allocate new arrays using the same data layout as the old structures
                        is_old_file = .true.
                        write(*,*) 'Could not read file header! Old energy file?'
			fileheader%arrays = 1
			fileheader%totresid = 1
			fileheader%version=' 5.06'
			if(ifile.eq.1) then
			allocate(fileheader%gcnum(fileheader%arrays),fileheader%types(fileheader%arrays),&
				fileheader%numres(fileheader%totresid),fileheader%resid(fileheader%totresid))
			end if
		else if ( fileheader%arrays .lt. 1) then
			write(*,*) 'Number of different types is < 1. Abort !'
			stop 'Qfep5 failed at loading file header'
		else if (canary .ne.myprec) then
			write(*,*) 'Precision of data in energy file is different from precision of this Q version'
			if (canary.eq.137) file_precision = 'Single'
			if (canary.eq.1337) file_precision = 'Double'
			if (canary.eq.13337) file_precision = 'Quadruple'
			write(*,'(a,a,a,a,a)') 'File precision is ',trim(file_precision),' , while version has ',trim(version_precision),' precision!'
			stop 'Mismatched precisions'
		else
			if(ifile.eq.1) then
			allocate(fileheader%types(fileheader%arrays),fileheader%numres(fileheader%arrays), &
			fileheader%gcnum(fileheader%arrays),fileheader%resid(fileheader%totresid),STAT=ERR)
!			write(*,*) fileheader%arrays,fileheader%totresid
			if (ERR .ne.0) then
			stop 'Qfep5 terminated abnormally: Out of memory when allocating arrays for header strucutre'
			end if
			end if
			rewind(f)
			read(f) canary,fileheader%arrays,fileheader%totresid,fileheader%types(1:fileheader%arrays),&
				fileheader%numres(1:fileheader%arrays),fileheader%resid(1:fileheader%totresid),&
				fileheader%gcnum(1:fileheader%arrays),fileheader%version
		end if
		if(ifile.eq.1) then
		do jj=1,nstates
		allocate(EQ(jj)%qx(fileheader%arrays),EQ(jj)%qq(fileheader%arrays),&
		EQ(jj)%qp(fileheader%arrays), EQ(jj)%qw(fileheader%arrays),&
                EQ(jj)%total(fileheader%arrays),EQ(jj)%q(fileheader%arrays),&
                EQ(jj)%restraint(fileheader%arrays),STAT=ERR)
		if(ERR /= 0) then
			stop 'Qfep5 terminated abnormally: Out of memory when allocating arrays. 11'
		end if
		EQ(jj)%total(:)=zero
		allocate(avEQ(jj)%qx(fileheader%arrays),avEQ(jj)%qq(fileheader%arrays),&
		avEQ(jj)%qp(fileheader%arrays),avEQ(jj)%qw(fileheader%arrays),&
                avEQ(jj)%total(fileheader%arrays),avEQ(jj)%q(fileheader%arrays),&
                avEQ(jj)%restraint(fileheader%arrays),STAT=ERR)

		if(ERR /= 0) then
			stop 'Qfep5 terminated abnormally: Out of memory when allocating arrays. 12'
		end if
		end do
		end if
		if(ifile.eq.1) then
	        !allocate large arrays
	        allocate(FEP(nfiles), FEPtmp%v(fileheader%arrays,nstates,mxpts), &
			FEPtmp%r(fileheader%arrays,nstates,mxpts), FEPtmp%lambda(mxstates), &
	                FEPtmp%vg(fileheader%arrays,mxpts), FEPtmp%gap(fileheader%arrays,mxpts), &
	                FEPtmp%c1(fileheader%arrays,mxpts), FEPtmp%c2(fileheader%arrays,mxpts), &
	                dgf(fileheader%arrays,0:nfiles+1),dgr(fileheader%arrays,0:nfiles+1), &
	                dgfsum(fileheader%arrays,0:nfiles+1),&
			dgrsum(fileheader%arrays,0:nfiles+1),dG(fileheader%arrays,0:nfiles+1), &
			gapmin(fileheader%arrays),gapmax(fileheader%arrays),sum(fileheader%arrays),&
	                veff1(fileheader%arrays),veff2(fileheader%arrays),dv(fileheader%arrays),STAT=ERR)
	        if(ERR /= 0) then
	                stop 'Qfep5 terminated abnormally: Out of memory when allocating arrays. 13'
	        end if
                do i=1,nfiles
			allocate(FEP(i)%lambda(mxstates))
		end do
		
		gapmin(:)=999
		gapmax(:)=-999
		FEPtmp%c1(:,:)= zero
		FEPtmp%c2(:,:)= zero
		FEP(ifile)%lambda(1:mxstates)= zero
		EQ(:)%lambda= zero
		 !size of secular determinant is allocated
	
	        allocate(Hij(nstates,nstates,fileheader%arrays),d(nstates),e(nstates),STAT=ERR)
	        if(ERR /= 0) then
	                write(*,*) 'ERROR: Out of memory when allocation Hij array.'
	                stop 'Qfep5 terminated abnormally: Out of memory.'
	        end if
		end if  !ifile .eq. 1
		write(*,'(a,a)') '# Energy files from Qdyn, version: ',fileheader%version
		!read 1st record to get lambdas
		if(is_old_file) rewind(f)
		idum = get_ene(f, EQ(:), offd, nstates,num_offd,fileheader%arrays)
		if(idum /= 0) then
			!we have a problem
			write(*,'(a,a)') 'ERROR: Unexpected end-of-file in first record of ', &
				trim(filnam)
			if(idum > 0) then
			!the problem is in energies
				write(*,'(a,i1)') 'while reading energies of state ',idum
				write(*,'(a,i1,a)') 'Maybe there are only ',idum-1, ' states?'
				stop 'Qfep5 terminated abnormally: Failed to read energies.'
			else
			!idum < 0 indicates problems with offdiags
				write(*,'(a)') 'while reading off-diagonal elements.'
				write(*,'(a,i1,a)') 'Maybe there are less than ',nnoffd, &
					' off-diagonal elements?'
				stop 'Qfep5 terminated abnormally: Failed to read off-diagonals.'
			end if
		end if

		FEP(ifile)%lambda(1:nstates) = EQ(1:nstates)%lambda

        rewind (f)

		ipt = 0
                do i=1,nstates
                do jj=1,fileheader%arrays !can not use Q energies operator any more
                avEQ(i)%q(jj)%bond     = zero
                avEQ(i)%q(jj)%angle    = zero
                avEQ(i)%q(jj)%torsion  = zero
                avEQ(i)%q(jj)%improper = zero
                avEQ(i)%restraint(jj)  = zero
                avEQ(i)%total(jj)      = zero
                avEQ(i)%qx(jj)%el      = zero
                avEQ(i)%qx(jj)%vdw     = zero
                avEQ(i)%qq(jj)%el      = zero
                avEQ(i)%qq(jj)%vdw     = zero
                avEQ(i)%qp(jj)%el      = zero
                avEQ(i)%qp(jj)%vdw     = zero
                avEQ(i)%qw(jj)%el      = zero
                avEQ(i)%qw(jj)%vdw     = zero
                end do
                end do

!		avEQ(:) = avEQ(:) * 0. !set all fields to zero using multiplication operator
		FEPtmp%gap(:,:) = zero
		if(.not.is_old_file) then
                        read(f) canary,fileheader%arrays,fileheader%totresid,fileheader%types(1:fileheader%arrays),&
                                fileheader%numres(1:fileheader%arrays),fileheader%resid(1:fileheader%totresid),&
                                fileheader%gcnum(1:fileheader%arrays),fileheader%version

		end if
		do while(get_ene(f, EQ(:), offd, nstates, num_offd,fileheader%arrays) == 0) !keep reading till EOF
			ipt = ipt + 1
			if (ipt .eq. (mxpts-1)) then
!ops, we need larger arrays to keep track of the temporary points
				call qfep_reallocate
			end if
	!!!!!!!!!!!!!!  masoud
!if the gas flag is > 0 then total energy will be just q (bonded), qq(nonbonded) and restraint
			if (gas .gt. 0 ) then
			do i=1,nstates
				do j=1,fileheader%arrays
				EQ(i)%total(j)=EQ(i)%q(j)%bond+EQ(i)%q(j)%angle+EQ(i)%q(j)%torsion+ &
					EQ(i)%q(j)%improper+EQ(i)%qq(j)%el+EQ(i)%qq(j)%vdw+EQ(i)%restraint(j)
				end do
!			print*, EQ(i)%total, EQ(i)%q,EQ(i)%qq,EQ(i)%restraint
!			print*, "''''''''''"
!			if (ipt > 10 ) stop
			end do
			end if
	!!!!!!!!!!!!!!   masoud
			if(ipt > nskip) then
				do i=1,nstates
				do jj=1,fileheader%arrays !can not use Q energies operator any more
				avEQ(i)%q(jj)%bond     = avEQ(i)%q(jj)%bond     + EQ(i)%q(jj)%bond
				avEQ(i)%q(jj)%angle    = avEQ(i)%q(jj)%angle    + EQ(i)%q(jj)%angle  
				avEQ(i)%q(jj)%torsion  = avEQ(i)%q(jj)%torsion  + EQ(i)%q(jj)%torsion  
				avEQ(i)%q(jj)%improper = avEQ(i)%q(jj)%improper + EQ(i)%q(jj)%improper  
				avEQ(i)%restraint(jj)  = avEQ(i)%restraint(jj)  + EQ(i)%restraint(jj)  
				avEQ(i)%total(jj)  = avEQ(i)%total(jj)  + EQ(i)%total(jj) 
				avEQ(i)%qx(jj)%el  = avEQ(i)%qx(jj)%el  + EQ(i)%qx(jj)%el
				avEQ(i)%qx(jj)%vdw = avEQ(i)%qx(jj)%vdw + EQ(i)%qx(jj)%vdw
				avEQ(i)%qq(jj)%el  = avEQ(i)%qq(jj)%el  + EQ(i)%qq(jj)%el
				avEQ(i)%qq(jj)%vdw = avEQ(i)%qq(jj)%vdw + EQ(i)%qq(jj)%vdw
				avEQ(i)%qp(jj)%el  = avEQ(i)%qp(jj)%el  + EQ(i)%qp(jj)%el
				avEQ(i)%qp(jj)%vdw = avEQ(i)%qp(jj)%vdw + EQ(i)%qp(jj)%vdw
				avEQ(i)%qw(jj)%el  = avEQ(i)%qw(jj)%el  + EQ(i)%qw(jj)%el
				avEQ(i)%qw(jj)%vdw = avEQ(i)%qw(jj)%vdw + EQ(i)%qw(jj)%vdw
				end do
				end do
			end if


!-------------------------------------------
! Correct H_ii with alfa, and modify H_ij...
!-------------------------------------------

            alfa(1)= zero
			do i=1,nstates
				do j=1,fileheader%arrays
				EQ(i)%total(j)=EQ(i)%total(j)+alfa(i)
				FEPtmp%v(j,i, ipt) = EQ(i)%total(j)
				if (num_offd .ne. 0) then
					FEPtmp%r(j,i,ipt) = offd(i)%rkl
			        end if
				end do
			end do
			do j=1,fileheader%arrays


			do i=1, noffd
				Hij(offd(i)%i, offd(i)%j,j) = offd(i)%Hij
			end do
			end do

			if ( scale_Hij .gt. zero ) then
				Hij(1,2,:) = scale_Hij*Hij(1,2,:)
			else
				do i=1,nstates
					do j=1,nstates
						do jj=1,fileheader%arrays
						if (i==j) then
 							Hij(i,j,jj)=EQ(i)%total(jj)
						else
 							if (A(i,j)==zero) then
 								Hij(i,j,jj) = A(j,i)*q_exp(-mu(j,i)*(offd(1)%rkl-rxy0(j,i)))* &
 								q_exp(-eta(j,i)*(offd(1)%rkl-rxy0(j,i))**2)
 							else
 								Hij(i,j,jj) = A(i,j)*q_exp(-mu(i,j)*(offd(1)%rkl-rxy0(i,j)))* &
 								q_exp(-eta(i,j)*(offd(1)%rkl-rxy0(i,j))**2)
 							end if
						end if
						end do
					end do
 				end do
		            end if

!-----------------------------------------------------------
! Ground state energy is calculated from secular determinant
!-----------------------------------------------------------

			do j=1,fileheader%arrays
				if (nstates==2) then
					FEPtmp%vg(j,ipt)=0.5_prec*(EQ(1)%total(j)+EQ(2)%total(j))-  &
						0.5_prec*q_sqrt( (EQ(1)%total(j)-EQ(2)%total(j))**2 + 4.0_prec*Hij(1,2,j)**2 )
					if(num_offd > 0) then
						FEPtmp%c1(j,ipt)=one/(one+((FEPtmp%vg(j,ipt)-EQ(1)%total(j))/Hij(1,2,j))**2)
						FEPtmp%c2(j,ipt)=one-FEPtmp%c1(j,ipt)
					end if
				else
					call tred2(Hij(:,:,j),nstates,d,e)
					call tqli(d,e,nstates,Hij(:,:,j))
					FEPtmp%vg(j,ipt)=MINVAL(d)
				end if

				do istate=1,nstates
					FEPtmp%gap(j,ipt)=FEPtmp%gap(j,ipt)+FEPtmp%v(j,istate,ipt)*coeff(istate)
				end do


				if(ipt .gt. nskip) then
					if(FEPtmp%gap(j,ipt) .lt. gapmin(j)) gapmin(j)=FEPtmp%gap(j,ipt)
					if(FEPtmp%gap(j,ipt) .gt. gapmax(j)) gapmax(j)=FEPtmp%gap(j,ipt)
				end if
			end do !j->arrays
	         end do  !(ipt)
		 close(f)
		FEP(ifile)%npts = ipt

		!copy FEPtmp to D
		allocate(FEP(ifile)%v(fileheader%arrays,nstates, FEP(ifile)%npts), &
			FEP(ifile)%vg(fileheader%arrays,FEP(ifile)%npts), &
			FEP(ifile)%gap(fileheader%arrays,FEP(ifile)%npts), &
			FEP(ifile)%c1(fileheader%arrays,FEP(ifile)%npts), &
			FEP(ifile)%c2(fileheader%arrays,FEP(ifile)%npts), &
			STAT=ERR)
		if(ERR /= 0) then
			stop 'Qfep5 terminated abnormally: Out of memory when allocating arrays. 2'
		end if
		FEP(ifile)%v(:,:,:) = FEPtmp%v(1:fileheader%arrays,1:nstates, 1:FEP(ifile)%npts)
		FEP(ifile)%vg(:,:) = FEPtmp%vg(1:fileheader%arrays,1:FEP(ifile)%npts)
		FEP(ifile)%gap(:,:) = FEPtmp%gap(1:fileheader%arrays,1:FEP(ifile)%npts)
		FEP(ifile)%c1(:,:) = FEPtmp%c1(1:fileheader%arrays,1:FEP(ifile)%npts)
		FEP(ifile)%c2(:,:) = FEPtmp%c2(1:fileheader%arrays,1:FEP(ifile)%npts)
		if(num_offd > 0) then
			allocate(FEP(ifile)%r(fileheader%arrays,nstates, FEP(ifile)%npts), STAT=ERR)
			if(ERR /= 0) then
				stop 'Qfep5 terminated abnormally: Out of memory when allocating arrays. 1'
			end if
			FEP(ifile)%r(:,:,:) = FEPtmp%r(1:fileheader%arrays,1:nstates, 1:FEP(ifile)%npts)
		end if


		! Average energies for each file is calculated

		if(ipt <= nskip) then !skipped all points in the file
			write(*,900) trim(filnam)
900			format('>>>>> ERROR: number of data sets in ',a,&
					' is less than number of points to skip!')
			stop 'Qfep5 terminated abnormally: Too few data points in file.'
		else
	                do i=1,nstates
	                do jj=1,fileheader%arrays !can not use Q energies operator any more
	                avEQ(i)%q(jj)%bond     = avEQ(i)%q(jj)%bond     * (one/(FEP(ifile)%npts-nskip))
	                avEQ(i)%q(jj)%angle    = avEQ(i)%q(jj)%angle    * (one/(FEP(ifile)%npts-nskip))
	                avEQ(i)%q(jj)%torsion  = avEQ(i)%q(jj)%torsion  * (one/(FEP(ifile)%npts-nskip))
	                avEQ(i)%q(jj)%improper = avEQ(i)%q(jj)%improper * (one/(FEP(ifile)%npts-nskip))
	                avEQ(i)%restraint(jj)  = avEQ(i)%restraint(jj)  * (one/(FEP(ifile)%npts-nskip))
	                avEQ(i)%total(jj)  = avEQ(i)%total(jj)  * (one/(FEP(ifile)%npts-nskip))
	                avEQ(i)%qx(jj)%el  = avEQ(i)%qx(jj)%el  * (one/(FEP(ifile)%npts-nskip))
	                avEQ(i)%qx(jj)%vdw = avEQ(i)%qx(jj)%vdw * (one/(FEP(ifile)%npts-nskip))
	                avEQ(i)%qq(jj)%el  = avEQ(i)%qq(jj)%el  * (one/(FEP(ifile)%npts-nskip))
	                avEQ(i)%qq(jj)%vdw = avEQ(i)%qq(jj)%vdw * (one/(FEP(ifile)%npts-nskip))
	                avEQ(i)%qp(jj)%el  = avEQ(i)%qp(jj)%el  * (one/(FEP(ifile)%npts-nskip))
	                avEQ(i)%qp(jj)%vdw = avEQ(i)%qp(jj)%vdw * (one/(FEP(ifile)%npts-nskip))
	                avEQ(i)%qw(jj)%el  = avEQ(i)%qw(jj)%el  * (one/(FEP(ifile)%npts-nskip))
	                avEQ(i)%qw(jj)%vdw = avEQ(i)%qw(jj)%vdw * (one/(FEP(ifile)%npts-nskip))
	                end do
	                end do
		end if
		do j=1,fileheader%arrays
!Write information for each of the calculations using different exclusion definitions
		do istate=1,nstates
		write(*,17) trim(filnam), istate, FEP(ifile)%npts-nskip, FEP(ifile)%lambda(istate), &
			avEQ(istate)%total(j),avEQ(istate)%q(j)%bond,avEQ(istate)%q(j)%angle,&
			avEQ(istate)%q(j)%torsion, &
			avEQ(istate)%q(j)%improper,avEQ(istate)%qx(j)%el,avEQ(istate)%qx(j)%vdw, &
			avEQ(istate)%qq(j)%el,avEQ(istate)%qq(j)%vdw,avEQ(istate)%qp(j)%el,&
			avEQ(istate)%qp(j)%vdw, &
			avEQ(istate)%qw(j)%el,avEQ(istate)%qw(j)%vdw,avEQ(istate)%restraint(j)

		end do
		end do !fileheader%arrays
17	format(a,1x,t23,1x,i2,1x,i9,1x,f9.6,1x,f8.2,1x,f8.2,1x,f8.2,1x,f8.2,1x,f8.2,1x,f8.2,1x,f8.2,1x,f8.2,1x,f8.2,1x,f8.2,1x,f8.2,1x,f8.2,1x,f8.2,1x,f8.2)
	end do  !ifile

	if(nfiles > 1) then !the following is meaningless for a single file
	do j=1,fileheader%arrays
		dgf(j,:)=    zero
		dgfsum(j,:)= zero
		sum(j)=      zero
		veff1(j)=    zero
		veff2(j)=    zero
		dv(j)=       zero
		do ifile=1,nfiles-1
			!check that number of points > number to skip
			if(nskip >= FEP(ifile)%npts) then
				write(*,999) ifile, nskip, FEP(ifile)%npts
999				format('File',i5,' contains only',i5,' points. Can''t skip',i5)
			end if
			do ipt=nskip+1,FEP(ifile)%npts
				do istate=1,nstates
					veff1(j)=veff1(j)+FEP(ifile)%lambda(istate)* &
					FEP(ifile)%v(j,istate,ipt)
   					veff2(j)=veff2(j)+FEP(ifile+1)%lambda(istate)* &
					FEP(ifile)%v(j,istate,ipt)
   	   			end do
				dv(j)=veff2(j)-veff1(j)
				veff1(j)= zero
				veff2(j)= zero
				sum(j)=sum(j)+q_exp(-dv(j)/rt)
			end do
			sum(j)=sum(j)/real(FEP(ifile)%npts-nskip, kind=prec)
			dgf(j,ifile)=-rt*q_logarithm(sum(j))
			dgfsum(j,ifile+1)=dgfsum(j,ifile)+dgf(j,ifile)
			sum(:)= zero
		end do
		dgrsum(j,:)= zero
		dgr(j,:)=    zero
		sum(j)=      zero
		veff1(j)=    zero
		veff2(j)=    zero
		dv(j)=       zero
		do ifile=nfiles,2,-1
			do ipt=nskip+1,FEP(ifile)%npts
				do istate=1,nstates
      				veff1(j)=veff1(j)+FEP(ifile)%lambda(istate)* &
				FEP(ifile)%v(j,istate,ipt)
     				veff2(j)=veff2(j)+FEP(ifile-1)%lambda(istate)* &
				FEP(ifile)%v(j,istate,ipt)
     			end do
				dv(j)=veff2(j)-veff1(j)
				veff1(j)=zero
				veff2(j)=zero
				sum(j)=sum(j)+q_exp(-dv(j)/rt)
			end do
			sum(j)=sum(j)/real(FEP(ifile)%npts-nskip, kind=prec)
			dgr(j,ifile)=-rt*q_logarithm(sum(j))
			dgrsum(j,ifile-1)=dgrsum(j,ifile)+dgr(j,ifile)
			sum(j)= zero
		end do
21		format('# Part 1: Free energy perturbation summary:')
22		format('# lambda(1)      dGf sum(dGf)      dGr sum(dGr)     <dG>')
        dG(j,1)= zero
        do ifile=2,nfiles
            dG(j,ifile)=dG(j,ifile-1)+0.5_prec*(dgf(j,ifile-1)-dgr(j,ifile))
        end do
	end do !FILEHEADER%ARRAYS
	write(*,*)
	write(*,*)
!	write(*,21)
	curres = 1
	do j=1,fileheader%arrays
66	format(/,'# Calculation for full system')
67	format(/,'# Calculation for system with full exclusion, residues ',2x,10i4)
68	format(/,'# Calculation for system with electrostatic exclusion, residues ',2x,10i4)
69	format(/,'# Calculation for system with vdW exclusion, residues ',2x,10i4)
169     format(/,'# Calculation for QCP, number of atoms ',2x,i4)
170     format(/,'# Calculation for QCP Mass Perturbation, number of atoms ',2x,i4)

	write(*,21)

	select case(fileheader%types(j))
		case(NOGC)
		write(*,66)
		case(FULL)
		write(*,67) (fileheader%resid(jj),jj=curres,curres+fileheader%numres(j)-1)
		case(ELECTRO)
		write(*,68) (fileheader%resid(jj),jj=curres,curres+fileheader%numres(j)-1)
		case(VDW)
		write(*,69) (fileheader%resid(jj),jj=curres,curres+fileheader%numres(j)-1)
                case(QCP_NORM)
                write(*,169) (fileheader%numres(j))
                case(QCP_MASSP)
                write(*,170) (fileheader%numres(j))
		case default
		if (is_old_file) then
		write(*,66)
		else
!Should never reach this, means memory corruption
		write(*,*) 'No such type of the calculation'
		stop 'ABNORMAL TERMINATION OF QFEP'
		end if
	end select
	curres = curres + fileheader%numres(j)
	write(*,22)
        do ifile=1,nfiles
           write (*,23) &
               FEP(ifile)%lambda(1),dgf(j,ifile-1),dgfsum(j,ifile), &
                             dgr(j,ifile+1),dgrsum(j,ifile),dG(j,ifile)
        end do
23		format(2x,f9.6,1x,f9.3,1x,f9.3,1x,f9.3,1x,f9.3,1x,f9.3)

		write (*,*)
		write (*,'(a,f9.2)') '# Min energy-gap is: ',gapmin(j)
		write (*,'(a,f9.2)') '# Max energy-gap is: ',gapmax(j)

	write(*,*)
!	end do !fileheader%array
		!-----------------------------------
		! Reaction free energy is calculated
		!-----------------------------------

        write (*,*)
	    write (*,*)
		write(*,24)
		write(*,25)
24		format('# Part 2: Reaction free energy summary:')
25		format('# Lambda(1)  bin Energy gap      dGa     dGb     dGg    # pts    c1**2    c2**2')
26		format(2x,f9.6,i5,2x,f9.2,2x,f9.2,2x,f9.2,2x,f9.2,2x,i5,1x,f9.3,1x,f9.3)
		gaprange=gapmax(j)-gapmin(j)      !Range of reaction coordinate
		xint=gaprange/real(nbins)   !Divide R.C. into bins
                        binsum(:,:) = zero
                        sumg(:)=      zero
                        sumv(:,:)=    zero
                        ptsum(:)=     zero
                        sumg2(:)=     zero
                        sumv2(:,:)=   zero
                        ibin = 0
                        dGv(:)=       zero
		do ifile=1,nfiles
			avdvv(:,:)= zero
			avdvg(:)=   zero
			sumv=       zero
			sumg=       zero
		    	avc1(:)=    zero
            		avc2(:)=    zero
			avr=        zero
			dvv=        zero
			dvg=        zero
			nbinpts(:)= 0

			do ipt=nskip+1,FEP(ifile)%npts
				ibin=int((FEP(ifile)%gap(j,ipt)-gapmin(j))/xint)+1
				veff= zero
				do istate=1,nstates
					veff=veff+FEP(ifile)%lambda(istate)*FEP(ifile)%v(j,istate,ipt)
				end do  !states
				do istate=1,nstates
				dvv(istate)=FEP(ifile)%v(j,istate,ipt)-veff
				end do !nstates
				dvg=FEP(ifile)%vg(j,ipt)-veff
				avdvv(:,ibin)=avdvv(:,ibin)+dvv(:)
				avdvg(ibin)=avdvg(ibin)+dvg
	         		avc1(ibin)=avc1(ibin)+FEP(ifile)%c1(j,ipt)
                		avc2(ibin)=avc2(ibin)+FEP(ifile)%c2(j,ipt)
				!Only gives first r_xy distance
				if(num_offd > 0) avr(ibin)=avr(ibin)+FEP(ifile)%r(j,1,ipt)
				nbinpts(ibin)=nbinpts(ibin)+1
			end do          !ipt
			do ibin=1,nbins
				if ( nbinpts(ibin) .ne. 0 ) then
			                avc1(ibin)=avc1(ibin)/real(nbinpts(ibin),kind=prec)
			                avc2(ibin)=avc2(ibin)/real(nbinpts(ibin),kind=prec)
					avr(ibin)=avr(ibin)/real(nbinpts(ibin),kind=prec)
					avdvv(:,ibin)=avdvv(:,ibin)/nbinpts(ibin)
					avdvg(ibin)=avdvg(ibin)/nbinpts(ibin)

				end if
			end do !ibin
			do ipt=nskip+1,FEP(ifile)%npts
				ibin=int((FEP(ifile)%gap(j,ipt)-gapmin(j))/xint)+1
				veff= zero
				do istate=1,nstates
					veff=veff+FEP(ifile)%lambda(istate)*FEP(ifile)%v(j,istate,ipt)
				end do  !istate

				do istate=1,nstates
					dvv(istate)=FEP(ifile)%v(j,istate,ipt)-veff-avdvv(istate,ibin)
				end do
				dvg=FEP(ifile)%vg(j,ipt)-veff-avdvg(ibin)
				sumv(:,ibin)=sumv(:,ibin)+q_exp2(-dvv(:)/rt)
				sumg(ibin)=sumg(ibin)+q_exp(-dvg/rt)
			end do   !ipt

			do ibin=1,nbins
	if (nbinpts(ibin).ge.nmin) then
	binsum(ibin,2)=binsum(ibin,2)+avc1(ibin)*nbinpts(ibin)
binsum(ibin,3)=binsum(ibin,3)+avc2(ibin)*nbinpts(ibin)
	binsum(ibin,4)=binsum(ibin,4)+avr(ibin)*nbinpts(ibin)  !Bin-averaged r_xy
	sumv(:,ibin)=sumv(:,ibin)/real(nbinpts(ibin))
sumg(ibin)=sumg(ibin)/real(nbinpts(ibin))

ptsum(ibin)=ptsum(ibin)+nbinpts(ibin)

	do istate=1,nstates
sumv2(istate,ibin)=-rt*q_logarithm(sumv(istate,ibin))+avdvv(istate,ibin)
	end do
sumg2(ibin)=-rt*q_logarithm(sumg(ibin))+avdvg(ibin)
	! These are the diabatic free energy curves
dGv(:)=dG(j,ifile)+sumv2(:,ibin)
	! This is the reaction free energy
dGg=dG(j,ifile)+sumg2(ibin)

binsum(ibin,1)=binsum(ibin,1)+dGg*int(nbinpts(ibin))

	write (*,26) FEP(ifile)%lambda(1),ibin, &
	gapmin(j)+real(ibin,kind=prec)*xint-xint/2., dGv(1),dGv(2),dGg, &
int(nbinpts(ibin)),avc1(ibin),avc2(ibin)
	end if
	end do  !ibin
	end do      !ifile
	write(*,*)
	write(*,*)
	write(*,27)
write(*,28)

	27		format('# Part 3: Bin-averaged summary:')
	28		format('# bin  energy gap  <dGg> <dGg norm> pts  <c1**2> <c2**2> <r_xy>')

	do ibin=1,nbins
	if (ptsum(ibin).ge.nmin) then
	binsum(ibin,1)=binsum(ibin,1)/real(ptsum(ibin),kind=prec) ! Bin-averaged reaction free energy
 	binsum(ibin,2)=binsum(ibin,2)/real(ptsum(ibin),kind=prec) ! Bin-averaged c1**2
 	binsum(ibin,3)=binsum(ibin,3)/real(ptsum(ibin),kind=prec) ! Bin-averaged c2**2
	binsum(ibin,4)=binsum(ibin,4)/real(ptsum(ibin),kind=prec) ! Bin-averaged r_xy
	end if
 	end do
 	min=MINVAL(binsum(:,1))
 	do ibin=1,nbins
		if (ptsum(ibin).ge.nmin) then
29		format(i4,1x,f9.2,1x,f9.2,1x,f9.2,2x,i9,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.2,1x,f8.2,1x,f8.2,1x,f8.2)
 		write(*,29) ibin,gapmin(j)+real(ibin,kind=prec)*xint-xint/2.,binsum(ibin,1),  &
 		binsum(ibin,1)-min,int(ptsum(ibin)),binsum(ibin,2),binsum(ibin,3),binsum(ibin,4)
        end if
	end do !ibin
	end do !fileheader%arrays
	end if !nfiles >1

	!clean up
	do ifile=1,nfiles
		deallocate(FEP(ifile)%v, FEP(ifile)%vg, FEP(ifile)%gap, &
		FEP(ifile)%c1, FEP(ifile)%c2,FEP(ifile)%lambda)
		if(num_offd > 0) deallocate(FEP(ifile)%r)
	end do
	deallocate(FEP)
	deallocate(FEPtmp%v,FEPtmp%vg,FEPtmp%gap,FEPtmp%c1,FEPtmp%c2,FEPtmp%lambda,FEPtmp%r)
	deallocate(Hij,d,e,STAT=ERR)
	if (allocated(fileheader%types))  deallocate(fileheader%types)
	if (allocated(fileheader%numres)) deallocate(fileheader%numres)
	if (allocated(fileheader%gcnum))  deallocate(fileheader%gcnum)
	if (allocated(fileheader%resid))  deallocate(fileheader%resid)
	do istate=1,nstates
		deallocate(EQ(istate)%qx,EQ(istate)%qq,EQ(istate)%qp,EQ(istate)%qw,EQ(istate)%total)
		deallocate(avEQ(istate)%qx,avEQ(istate)%qq,avEQ(istate)%qp,avEQ(istate)%qw,avEQ(istate)%total)
                deallocate(EQ(istate)%q,EQ(istate)%restraint)
                deallocate(avEQ(istate)%q,avEQ(istate)%restraint)
	end do
	deallocate(sumg,sumg2,avdvg,avc1,avc2,avr,binsum,nbinpts,ptsum,EQ,avEQ,dvv,dGv,alfa,coeff,&
		A,mu,eta,rxy0,avdvv,sumv,sumv2,offd)
!.......................................................................

contains

!------------------------------

subroutine qfep_reallocate
!this nasty piece makes sure that the reallocated tmp arrays look right
! locals
type(FEP_DATA_TYPE)			:: FEPlocal
integer					:: oldmxpts

	oldmxpts = mxpts
	mxpts = oldmxpts + int(oldmxpts * 0.25_prec) + 200
	allocate(FEPlocal%v(fileheader%arrays,nstates,oldmxpts), &
		FEPlocal%r(fileheader%arrays,nstates,oldmxpts), FEPlocal%lambda(mxstates), &
		FEPlocal%vg(fileheader%arrays,oldmxpts), FEPlocal%gap(fileheader%arrays,oldmxpts), &
		FEPlocal%c1(fileheader%arrays,oldmxpts), FEPlocal%c2(fileheader%arrays,oldmxpts), &
		STAT=ERR)
                if(ERR /= 0) then
                        stop '# >>> ERROR! Qfep5 terminated abnormally: Out of memory when allocating arrays for reallocation'
                end if
		FEPlocal%c1(:,:)  = zero
		FEPlocal%c2(:,:)  = zero
		FEPlocal%gap(:,:) = zero
	
	FEPlocal%v(1:fileheader%arrays,1:nstates,1:oldmxpts)  = FEPtmp%v(1:fileheader%arrays,1:nstates,1:oldmxpts)
	FEPlocal%r(1:fileheader%arrays,1:nstates,1:oldmxpts)  = FEPtmp%r(1:fileheader%arrays,1:nstates,1:oldmxpts)
	FEPlocal%vg(1:fileheader%arrays,1:oldmxpts)           = FEPtmp%vg(1:fileheader%arrays,1:oldmxpts)
	FEPlocal%gap(1:fileheader%arrays,1:oldmxpts)          = FEPtmp%gap(1:fileheader%arrays,1:oldmxpts)
	FEPlocal%c1(1:fileheader%arrays,1:oldmxpts)           = FEPtmp%c1(1:fileheader%arrays,1:oldmxpts)
	FEPlocal%c2(1:fileheader%arrays,1:oldmxpts)           = FEPtmp%c2(1:fileheader%arrays,1:oldmxpts)
	FEPlocal%lambda(1:nstates) = FEPtmp%lambda(1:nstates)

	deallocate(FEPtmp%lambda)
	deallocate(FEPtmp%v)
	deallocate(FEPtmp%r)
	deallocate(FEPtmp%vg)
	deallocate(FEPtmp%gap)
	deallocate(FEPtmp%c1)
	deallocate(FEPtmp%c2)

	allocate(FEPtmp%v(fileheader%arrays,nstates,mxpts), &
	FEPtmp%r(fileheader%arrays,nstates,mxpts), FEPtmp%lambda(mxstates), &
	FEPtmp%vg(fileheader%arrays,mxpts), FEPtmp%gap(fileheader%arrays,mxpts), &
	FEPtmp%c1(fileheader%arrays,mxpts), FEPtmp%c2(fileheader%arrays,mxpts), &
	STAT=ERR)
	if(ERR /= 0) then
		stop '# >>> ERROR! Qfep5 terminated abnormally: Out of memory when reallocating arrays for energy calculation'
	end if
	FEPtmp%c1(:,:)  = zero
	FEPtmp%c2(:,:)  = zero
	FEPtmp%gap(:,:) = zero
	FEPtmp%v(1:fileheader%arrays,1:nstates,1:oldmxpts)  = FEPlocal%v(1:fileheader%arrays,1:nstates,1:oldmxpts)
	FEPtmp%r(1:fileheader%arrays,1:nstates,1:oldmxpts)  = FEPlocal%r(1:fileheader%arrays,1:nstates,1:oldmxpts)
	FEPtmp%vg(1:fileheader%arrays,1:oldmxpts)           = FEPlocal%vg(1:fileheader%arrays,1:oldmxpts)
	FEPtmp%gap(1:fileheader%arrays,1:oldmxpts)          = FEPlocal%gap(1:fileheader%arrays,1:oldmxpts)
	FEPtmp%c1(1:fileheader%arrays,1:oldmxpts)           = FEPlocal%c1(1:fileheader%arrays,1:oldmxpts)
	FEPtmp%c2(1:fileheader%arrays,1:oldmxpts)           = FEPlocal%c2(1:fileheader%arrays,1:oldmxpts)
	FEPtmp%lambda(1:nstates) = FEPlocal%lambda(1:nstates)

        deallocate(FEPlocal%lambda)
        deallocate(FEPlocal%v)
        deallocate(FEPlocal%r)
        deallocate(FEPlocal%vg)
        deallocate(FEPlocal%gap)
        deallocate(FEPlocal%c1)
        deallocate(FEPlocal%c2)


end subroutine qfep_reallocate

subroutine prompt (outtxt)
	character(*) outtxt
#if defined (__osf__)
	!prompt to STDERR using unit 5=STDIN on OSF/1=DEC UNIX
	integer, parameter			::	f=5
	!write (f,'($,a)') outtxt
#elseif defined (_WIN32)
	!open the ERR file on Win32
	integer, save :: f
	if(f==0) then
		f=17
		open(unit=f,file='ERR')
	end if
	!write (f1,'($,a)') outtxt
#else
	!otherwise prompt to STDOUT
	integer, parameter			::	f=6
	!write (f2,'($,a)') outtxt
#endif
write (f,'(a,$)') outtxt
end subroutine prompt

!------------------------------

SUBROUTINE TRED2(A,N,D,E)
!------------------------------------------------------------
! This subroutine reduces a symmetric matrix to tridiagonal
! form. The tridiagonal matrix can further be diagonalized by
! the subroutine tqli.
! These subroutines were copied by Karin Kolmodin 20 Nov. 1997
! from http://rsc.anu.au/HWS/COURSES/MATHMETH/node70.html
! and rewritten in f90.
!------------------------------------------------------------
	  real(kind=prec),dimension(:)      ::D,E
	  real(kind=prec),dimension(:,:)    ::A
          integer                   ::I,J,K,L,N
	  real(kind=prec)                   ::SCALE,F,G,H,HH

      IF(N.GT.1)THEN
        DO I=N,2,-1
           L=I-1
           H=zero
           SCALE=zero
          IF(L.GT.1)THEN
            DO K=1,L
              SCALE=SCALE+ABS(A(I,K))
            end do
            IF(SCALE.EQ.zero)THEN
              E(I)=A(I,L)
            ELSE
              DO K=1,L
                A(I,K)=A(I,K)/SCALE
                H=H+A(I,K)**2.0_prec
              end do
              F=A(I,L)
              G=-SIGN(q_sqrt(H),F)
              E(I)=SCALE*G
              H=H-F*G
              A(I,L)=F-G
              F=zero
              DO J=1,L
                A(J,I)=A(I,J)/H
                G=zero
                DO K=1,J
                  G=G+A(J,K)*A(I,K)
                end do
                IF(L.GT.J)THEN
                  DO K=J+1,L
                    G=G+A(K,J)*A(I,K)
                  end do
                END IF
                E(J)=G/H
                F=F+E(J)*A(I,J)
              end do
              HH=F/(H+H)
              DO J=1,L
                F=A(I,J)
                G=E(J)-HH*F
                E(J)=G
                DO K=1,J
                  A(J,K)=A(J,K)-F*E(K)-G*A(I,K)
                end do
              end do
            END IF
          ELSE
            E(I)=A(I,L)
          END IF
          D(I)=H
        end do
      END IF
      D(1)=zero
      E(1)=zero
      DO I=1,N
        L=I-1
        IF(D(I).NE.zero)THEN
          DO J=1,L
            G=zero
            DO K=1,L
              G=G+A(I,K)*A(K,J)
            end do
            DO K=1,L
              A(K,J)=A(K,J)-G*A(K,I)
            end do
          end do
        END IF
        D(I)=A(I,I)
        A(I,I)=one
        IF(L.GE.1)THEN
          DO J=1,L
            A(I,J)=zero
            A(J,I)=zero
          end do
        END IF
      end do
END subroutine TRED2

!-----------------------------------

SUBROUTINE TQLI(D,E,N,Z)
!------------------------------------------------------------
! This subroutine diagonalizes a tridiagonal matrix which has
! been prepared by the subroutine tred2.
! These subroutines were copied by Karin Kolmodin 20 Nov. 1997
! from http://rsc.anu.au/HWS/COURSES/MATHMETH/node70.html
! and rewritten in f90
!------------------------------------------------------------
   implicit none
      real(kind=prec),dimension(:)               ::D,E
	  real(kind=prec),dimension(:,:)             ::Z
	  integer                            ::I,N,K,L,M,ITER
	  real(kind=prec)                            ::DD,G,R,S,C,P,F,B

        IF (N.GT.1) THEN
        DO I=2,N
          E(I-1)=E(I)
        end do
        E(N)=zero
        DO L=1,N
          ITER=0
1         DO M=L,N-1
            DD=q_abs(D(M))+q_abs(D(M+1))
            IF (q_abs(E(M))+DD.EQ.DD) GO TO 2
          end do
          M=N
2         IF(M.NE.L)THEN
            IF(ITER.EQ.30)STOP 'too many iterations'
            ITER=ITER+1
            G=(D(L+1)-D(L))/(2.0_prec*E(L))
            R=q_sqrt(G**2.0_prec+one)
            G=D(M)-D(L)+E(L)/(G+SIGN(R,G))
            S=one
            C=one
            P=zero
            DO I=M-1,L,-1
              F=S*E(I)
              B=C*E(I)
              IF(q_abs(F).GE.q_abs(G))THEN
                C=G/F
                R=q_sqrt(C**2.0_prec+one)
                E(I+1)=F*R
                S=one/R
                C=C*S
              ELSE
                S=F/G
                R=q_sqrt(S**2.0_prec+one)
                E(I+1)=G*R
                C=one/R
                S=S*C
              END IF
              G=D(I+1)-P
              R=(D(I)-G)*S+2.0_prec*C*B
              P=S*R
              D(I+1)=G+P
              G=C*R-B
              DO K=1,N
                F=Z(K,I+1)
                Z(K,I+1)=S*Z(K,I)+C*F
                Z(K,I)=C*Z(K,I)-S*F
              end do
            end do
            D(L)=D(L)-P
            E(L)=G
            E(M)=zero
            GO TO 1
          END IF
        end do
      END IF
      RETURN
END subroutine TQLI

!-----------------------------------

end program qfep
