! f95 -fdefault-real-8 -o climate climate.F95
! ncdump -ff -v time out.nc | less

PROGRAM main

USE netcdf

USE state, ONLY: pre1,qbar,the1
USE state, ONLY: ubar,vbar,wbar
USE state, ONLY: util,vtil,wtil
USE state, ONLY: pre2,the1ref
USE state, ONLY: prec

USE state0, ONLY: pre10,qbar0,the10
USE state0, ONLY: ubar0,vbar0,wbar0
USE state0, ONLY: util0,vtil0,wtil0
USE state0, ONLY: pre20,the1ref0
USE state0, ONLY: prec0

USE auxiliaries, ONLY: rhoutil,rhovtil,rhowtil

USE auxiliaries, ONLY: pre1_,qbar_,the1_
USE auxiliaries, ONLY: ubar_,vbar_,wbar_
USE auxiliaries, ONLY: util_,vtil_,wtil_
USE auxiliaries, ONLY: pre2_

USE auxiliaries, ONLY: pre10_,qbar0_,the10_
USE auxiliaries, ONLY: ubar0_,vbar0_,wbar0_
USE auxiliaries, ONLY: util0_,vtil0_,wtil0_
USE auxiliaries, ONLY: pre20_

USE auxiliaries, ONLY: fqbar1,fthe11
USE auxiliaries, ONLY: frhoutil1,frhovtil1,frhowtil1

USE auxiliaries, ONLY: fqbar2,fthe12
USE auxiliaries, ONLY: frhoutil2,frhovtil2,frhowtil2

USE auxiliaries, ONLY: fqbar3,fthe13
USE auxiliaries, ONLY: frhoutil3,frhovtil3,frhowtil3

USE auxiliaries, ONLY: fqbar4,fthe14
USE auxiliaries, ONLY: frhoutil4,frhovtil4,frhowtil4

USE parameters, ONLY: num_ps_steps
USE parameters, ONLY: num_mc_steps
USE parameters, ONLY: print_lag_err
USE parameters, ONLY: print_lag_val
USE parameters, ONLY: nx,ny,nz
USE parameters, ONLY: dt_mc,dt
USE parameters, ONLY: sv_input_start
USE parameters, ONLY: tv_input_start
USE parameters, ONLY: the1ref_input_start
USE parameters, ONLY: the1ref0_input_start
USE parameters, ONLY: flush_time
USE parameters, ONLY: resume
USE parameters, ONLY: delta0,delta

USE psi_solver, ONLY: varpsi
USE phi_solver, ONLY: varphi

USE state, ONLY: check_borders,check_borders_the1ref
USE state0, ONLY: check_borders0,check_borders_the1ref0

USE parameters, ONLY: set_coords

USE evolve, ONLY: evolve_state
USE evolve0, ONLY: evolve_state0

USE psi_solver, ONLY: set_psi,set_psi2,set_psi3
USE phi_solver, ONLY: set_phi,set_phi2

USE utilities, ONLY: rho0,rho1,precip

USE state, ONLY: allocate_state
USE state0, ONLY: allocate_state0
USE psi_solver, ONLY: allocate_varpsi
USE phi_solver, ONLY: allocate_varphi
USE auxiliaries, ONLY: allocate_auxiliaries

IMPLICIT NONE

! netCDF PARAMETERS and VARS:

INTEGER,PARAMETER :: NDIMS=4

INTEGER :: pre1_varid,qbar_varid,the1_varid
INTEGER :: ubar_varid,vbar_varid,wbar_varid
INTEGER :: util_varid,vtil_varid,wtil_varid
INTEGER :: pre2_varid,the1ref_varid
INTEGER :: prec_varid

INTEGER :: pre10_varid,qbar0_varid,the10_varid
INTEGER :: ubar0_varid,vbar0_varid,wbar0_varid
INTEGER :: util0_varid,vtil0_varid,wtil0_varid
INTEGER :: pre20_varid,the1ref0_varid
INTEGER :: prec0_varid

INTEGER :: qbar1_varid,the11_varid
INTEGER :: util1_varid,vtil1_varid,wtil1_varid
INTEGER :: prec1_varid

INTEGER :: time_varid
INTEGER :: ncid,dimids(NDIMS)
INTEGER :: x_dimid,y_dimid,z_dimid,t_dimid
INTEGER :: x1_dimid,y1_dimid,z1_dimid,t1_dimid
INTEGER :: start(NDIMS)=(/1,1,1,1/)

CHARACTER (LEN=*),PARAMETER :: THE1REF_FILE_NAME="/project/amp/pkrause/FINAL/IN/the1ref.nc"
CHARACTER (LEN=*),PARAMETER :: THE1REF0_FILE_NAME="/project/amp/pkrause/FINAL/IN/the1ref0.nc"
CHARACTER (LEN=*),PARAMETER :: INPUT_SV_FILE_NAME="/project/amp/pkrause/FINAL/IN/in.nc"
CHARACTER (LEN=*),PARAMETER :: INPUT_TV_FILE_NAME="/project/amp/pkrause/FINAL/IN/in.nc"
CHARACTER (LEN=*),PARAMETER :: REFERENCE_FILE_NAME="out.nc"
CHARACTER (LEN=*),PARAMETER :: OUTPUT_FILE_NAME="out.nc" ! prior and true values
CHARACTER (LEN=*),PARAMETER :: OUTPUT_FILE_NAME_FLUSHED="out_flushed.nc" ! prior and true values
CHARACTER (LEN=*),PARAMETER :: PERT_FILE_NAME="perturbed_values.nc"

! AUXILIARIES:

INTEGER :: cc=1
INTEGER :: nn,kk,mm
INTEGER :: i,j,k

REAL :: time_ps

! TIME:

REAL :: time ! seconds

! netCDF count:

INTEGER :: count(NDIMS)=(/nx+1,ny+1,nz,1/)

IF (MOD(ny-1,4)/=0) THEN
        WRITE(*,*) 'Call to EXIT: ny-1 must be multiple of 4.'
        CALL exit(0)
END IF

CALL allocate_state()
CALL allocate_state0()

IF (resume==1) THEN
        CALL allocate_varpsi()
        CALL allocate_varphi()
        CALL allocate_auxiliaries()

        CALL set_coords()

        CALL input_resume_values_and_time()
        CALL input_the1ref()
        CALL input_the1ref0()

        ! CREATE ERR FILE:

        OPEN(UNIT=2,FILE='prior_err.dat',ACTION="write",STATUS="old")
        CLOSE(2)

        OPEN(UNIT=2,FILE='prior_err_flushed.dat',ACTION="write",STATUS="old")
        CLOSE(2)
ELSE
        CALL allocate_varpsi()
        CALL allocate_varphi()
        CALL allocate_auxiliaries()

        CALL set_coords()

        CALL input_initial_sample_values()
        CALL input_initial_true_values_and_time()
        CALL input_the1ref()
        CALL input_the1ref0()

        ! CREATE ERR FILE:

        OPEN(UNIT=2,FILE='prior_err.dat',ACTION="write",STATUS="replace")
        CLOSE(2)

        OPEN(UNIT=2,FILE='prior_err_flushed.dat',ACTION="write",STATUS="replace")
        CLOSE(2)

        ! WRITE INITIAL PRIOR AND TRUE VALUES OF ALL VARS:
        ! WRITE INITIAL ERRORS OF TRAILING VARS:

        CALL set_precipitation(pre1,qbar,the1,wbar,wtil,prec)
        CALL set_precipitation(pre10,qbar0,the10,wbar0,wtil0,prec0)
        CALL create_and_write_prior_values_file(0)
        CALL write_prior_rms_errors_of_trailing_vars(0)

        DO k=1,nz
                DO j=0,ny
                        DO i=0,nx
                                pre1_(i,j,k)=pre1(i,j,k)
                                qbar_(i,j,k)=qbar(i,j,k)
                                the1_(i,j,k)=the1(i,j,k)
                                ubar_(i,j,k)=ubar(i,j,k)
                                vbar_(i,j,k)=vbar(i,j,k)
                                wbar_(i,j,k)=wbar(i,j,k)
                                util_(i,j,k)=util(i,j,k)
                                vtil_(i,j,k)=vtil(i,j,k)
                                wtil_(i,j,k)=wtil(i,j,k)
                                pre2_(i,j,k)=pre2(i,j,k)

                                pre10_(i,j,k)=pre10(i,j,k)
                                qbar0_(i,j,k)=qbar0(i,j,k)
                                the10_(i,j,k)=the10(i,j,k)
                                ubar0_(i,j,k)=ubar0(i,j,k)
                                vbar0_(i,j,k)=vbar0(i,j,k)
                                wbar0_(i,j,k)=wbar0(i,j,k)
                                util0_(i,j,k)=util0(i,j,k)
                                vtil0_(i,j,k)=vtil0(i,j,k)
                                wtil0_(i,j,k)=wtil0(i,j,k)
                                pre20_(i,j,k)=pre20(i,j,k)
                        END DO
                END DO
        END DO

        ! FLUSH: BEGIN

        CALL set_psi2(ubar,vbar,wbar)
        CALL set_phi2(pre1,qbar,the1,util,vtil,wtil)

        CALL set_psi(time,pre1,qbar,the1,ubar,vbar,wbar,util,vtil,wtil)

        DO k=1,nz
                DO j=1,ny-1
                        DO i=1,nx-1
                                pre1(i,j,k)=varpsi(i,j,k)
                        END DO
                END DO
        END DO

        DO k=1,nz
                DO j=1,ny-1
                        pre1(0,j,k)=pre1(nx-1,j,k)
                        pre1(nx,j,k)=pre1(1,j,k)
                END DO

                DO i=0,nx
                        pre1(i,0,k)=pre1(i,ny-1,k)
                        pre1(i,ny,k)=pre1(i,1,k)
                END DO
        END DO

        CALL set_phi(delta,time,pre1,qbar,the1,ubar,vbar,wbar,util,vtil,wtil)

        DO k=1,nz
                DO j=1,ny-1
                        DO i=1,nx-1
                                pre2(i,j,k)=varphi(i,j,k)
                        END DO
                END DO
        END DO

        DO k=1,nz
                DO j=1,ny-1
                        pre2(0,j,k)=pre2(nx-1,j,k)
                        pre2(nx,j,k)=pre2(1,j,k)
                END DO

                DO i=0,nx
                        pre2(i,0,k)=pre2(i,ny-1,k)
                        pre2(i,ny,k)=pre2(i,1,k)
                END DO
        END DO

        CALL set_psi2(ubar0,vbar0,wbar0)
        CALL set_phi2(pre10,qbar0,the10,util0,vtil0,wtil0)

        CALL set_psi(time,pre10,qbar0,the10,ubar0,vbar0,wbar0,util0,vtil0,wtil0)

        DO k=1,nz
                DO j=1,ny-1
                        DO i=1,nx-1
                                pre10(i,j,k)=varpsi(i,j,k)
                        END DO
                END DO
        END DO

        DO k=1,nz
                DO j=1,ny-1
                        pre10(0,j,k)=pre10(nx-1,j,k)
                        pre10(nx,j,k)=pre10(1,j,k)
                END DO

                DO i=0,nx
                        pre10(i,0,k)=pre10(i,ny-1,k)
                        pre10(i,ny,k)=pre10(i,1,k)
                END DO
        END DO

        CALL set_phi(delta0,time,pre10,qbar0,the10,ubar0,vbar0,wbar0,util0,vtil0,wtil0)

        DO k=1,nz
                DO j=1,ny-1
                        DO i=1,nx-1
                                pre20(i,j,k)=varphi(i,j,k)
                        END DO
                END DO
        END DO

        DO k=1,nz
                DO j=1,ny-1
                        pre20(0,j,k)=pre20(nx-1,j,k)
                        pre20(nx,j,k)=pre20(1,j,k)
                END DO

                DO i=0,nx
                        pre20(i,0,k)=pre20(i,ny-1,k)
                        pre20(i,ny,k)=pre20(i,1,k)
                END DO
        END DO

        ! FLUSH: END

        CALL set_precipitation(pre1,qbar,the1,wbar,wtil,prec)
        CALL set_precipitation(pre10,qbar0,the10,wbar0,wtil0,prec0)
        CALL create_and_write_prior_values_file(1)
        CALL write_prior_rms_errors_of_trailing_vars(1)

        DO k=1,nz
                DO j=0,ny
                        DO i=0,nx
                                pre1(i,j,k)=pre1_(i,j,k)
                                qbar(i,j,k)=qbar_(i,j,k)
                                the1(i,j,k)=the1_(i,j,k)
                                ubar(i,j,k)=ubar_(i,j,k)
                                vbar(i,j,k)=vbar_(i,j,k)
                                wbar(i,j,k)=wbar_(i,j,k)
                                util(i,j,k)=util_(i,j,k)
                                vtil(i,j,k)=vtil_(i,j,k)
                                wtil(i,j,k)=wtil_(i,j,k)
                                pre2(i,j,k)=pre2_(i,j,k)

                                pre10(i,j,k)=pre10_(i,j,k)
                                qbar0(i,j,k)=qbar0_(i,j,k)
                                the10(i,j,k)=the10_(i,j,k)
                                ubar0(i,j,k)=ubar0_(i,j,k)
                                vbar0(i,j,k)=vbar0_(i,j,k)
                                wbar0(i,j,k)=wbar0_(i,j,k)
                                util0(i,j,k)=util0_(i,j,k)
                                vtil0(i,j,k)=vtil0_(i,j,k)
                                wtil0(i,j,k)=wtil0_(i,j,k)
                                pre20(i,j,k)=pre20_(i,j,k)
                        END DO
                END DO
        END DO

END IF

! START:

DO nn=cc-1,num_ps_steps

        CALL input_reference()

        IF (MOD(time,flush_time)==0) THEN
                ! FLUSH: BEGIN

                CALL set_psi2(ubar,vbar,wbar)
                CALL set_phi2(pre1,qbar,the1,util,vtil,wtil)

                CALL set_psi(time,pre1,qbar,the1,ubar,vbar,wbar,util,vtil,wtil)

                DO k=1,nz
                        DO j=1,ny-1
                                DO i=1,nx-1
                                        pre1(i,j,k)=varpsi(i,j,k)
                                END DO
                        END DO
                END DO

                DO k=1,nz
                        DO j=1,ny-1
                                pre1(0,j,k)=pre1(nx-1,j,k)
                                pre1(nx,j,k)=pre1(1,j,k)
                        END DO

                        DO i=0,nx
                                pre1(i,0,k)=pre1(i,ny-1,k)
                                pre1(i,ny,k)=pre1(i,1,k)
                        END DO
                END DO

                CALL set_phi(delta,time,pre1,qbar,the1,ubar,vbar,wbar,util,vtil,wtil)

                DO k=1,nz
                        DO j=1,ny-1
                                DO i=1,nx-1
                                        pre2(i,j,k)=varphi(i,j,k)
                                END DO
                        END DO
                END DO

                DO k=1,nz
                        DO j=1,ny-1
                                pre2(0,j,k)=pre2(nx-1,j,k)
                                pre2(nx,j,k)=pre2(1,j,k)
                        END DO

                        DO i=0,nx
                                pre2(i,0,k)=pre2(i,ny-1,k)
                                pre2(i,ny,k)=pre2(i,1,k)
                        END DO
                END DO

                CALL set_psi2(ubar0,vbar0,wbar0)
                CALL set_phi2(pre10,qbar0,the10,util0,vtil0,wtil0)

                CALL set_psi(time,pre10,qbar0,the10,ubar0,vbar0,wbar0,util0,vtil0,wtil0)

                DO k=1,nz
                        DO j=1,ny-1
                                DO i=1,nx-1
                                        pre10(i,j,k)=varpsi(i,j,k)
                                END DO
                        END DO
                END DO

                DO k=1,nz
                        DO j=1,ny-1
                                pre10(0,j,k)=pre10(nx-1,j,k)
                                pre10(nx,j,k)=pre10(1,j,k)
                        END DO

                        DO i=0,nx
                                pre10(i,0,k)=pre10(i,ny-1,k)
                                pre10(i,ny,k)=pre10(i,1,k)
                        END DO
                END DO

                CALL set_phi(delta0,time,pre10,qbar0,the10,ubar0,vbar0,wbar0,util0,vtil0,wtil0)

                DO k=1,nz
                        DO j=1,ny-1
                                DO i=1,nx-1
                                        pre20(i,j,k)=varphi(i,j,k)
                                END DO
                        END DO
                END DO

                DO k=1,nz
                        DO j=1,ny-1
                                pre20(0,j,k)=pre20(nx-1,j,k)
                                pre20(nx,j,k)=pre20(1,j,k)
                        END DO

                        DO i=0,nx
                                pre20(i,0,k)=pre20(i,ny-1,k)
                                pre20(i,ny,k)=pre20(i,1,k)
                        END DO
                END DO

                ! FLUSH: END
        END IF

        ! SET POSTERIOR the1, ubar, vbar:

        ! DO k=1,nz
                ! DO j=0,ny
                        ! DO i=0,nx
                                ! the1(i,j,k)=the10(i,j,k)
                                ! ubar(i,j,k)=ubar0(i,j,k)
                                ! vbar(i,j,k)=vbar0(i,j,k)
                        ! END DO
                ! END DO
        ! END DO

        ! EVOLVE STATE: BEGIN

        dt=dt_mc

        DO mm=1,num_mc_steps

                IF (MOD(time,flush_time)==0) THEN ! .AND.MOD(mm,3)/=0) THEN
                        ! FLUSH: BEGIN

                        CALL set_psi2(ubar,vbar,wbar)
                        CALL set_phi2(pre1,qbar,the1,util,vtil,wtil)

                        CALL set_psi(time,pre1,qbar,the1,ubar,vbar,wbar,util,vtil,wtil)

                        DO k=1,nz
                                DO j=1,ny-1
                                        DO i=1,nx-1
                                                pre1(i,j,k)=varpsi(i,j,k)
                                        END DO
                                END DO
                        END DO

                        DO k=1,nz
                                DO j=1,ny-1
                                        pre1(0,j,k)=pre1(nx-1,j,k)
                                        pre1(nx,j,k)=pre1(1,j,k)
                                END DO

                                DO i=0,nx
                                        pre1(i,0,k)=pre1(i,ny-1,k)
                                        pre1(i,ny,k)=pre1(i,1,k)
                                END DO
                        END DO

                        CALL set_phi(delta,time,pre1,qbar,the1,ubar,vbar,wbar,util,vtil,wtil)

                        DO k=1,nz
                                DO j=1,ny-1
                                        DO i=1,nx-1
                                                pre2(i,j,k)=varphi(i,j,k)
                                        END DO
                                END DO
                        END DO

                        DO k=1,nz
                                DO j=1,ny-1
                                        pre2(0,j,k)=pre2(nx-1,j,k)
                                        pre2(nx,j,k)=pre2(1,j,k)
                                END DO

                                DO i=0,nx
                                        pre2(i,0,k)=pre2(i,ny-1,k)
                                        pre2(i,ny,k)=pre2(i,1,k)
                                END DO
                        END DO

                        CALL set_psi2(ubar0,vbar0,wbar0)
                        CALL set_phi2(pre10,qbar0,the10,util0,vtil0,wtil0)

                        CALL set_psi(time,pre10,qbar0,the10,ubar0,vbar0,wbar0,util0,vtil0,wtil0)

                        DO k=1,nz
                                DO j=1,ny-1
                                        DO i=1,nx-1
                                                pre10(i,j,k)=varpsi(i,j,k)
                                        END DO
                                END DO
                        END DO

                        DO k=1,nz
                                DO j=1,ny-1
                                        pre10(0,j,k)=pre10(nx-1,j,k)
                                        pre10(nx,j,k)=pre10(1,j,k)
                                END DO

                                DO i=0,nx
                                        pre10(i,0,k)=pre10(i,ny-1,k)
                                        pre10(i,ny,k)=pre10(i,1,k)
                                END DO
                        END DO

                        CALL set_phi(delta0,time,pre10,qbar0,the10,ubar0,vbar0,wbar0,util0,vtil0,wtil0)

                        DO k=1,nz
                                DO j=1,ny-1
                                        DO i=1,nx-1
                                                pre20(i,j,k)=varphi(i,j,k)
                                        END DO
                                END DO
                        END DO

                        DO k=1,nz
                                DO j=1,ny-1
                                        pre20(0,j,k)=pre20(nx-1,j,k)
                                        pre20(nx,j,k)=pre20(1,j,k)
                                END DO

                                DO i=0,nx
                                        pre20(i,0,k)=pre20(i,ny-1,k)
                                        pre20(i,ny,k)=pre20(i,1,k)
                                END DO
                        END DO

                        ! FLUSH: END
                END IF

                CALL evolve_state(time)
                CALL evolve_state0(time)

                time=time+dt

                ! WRITE PRIOR RMS ERRORS OF TRAILING VARS:

                IF (MOD(mm,print_lag_err)==0) THEN
                        CALL set_precipitation(pre1,qbar,the1,wbar,wtil,prec)
                        CALL set_precipitation(pre10,qbar0,the10,wbar0,wtil0,prec0)
                        CALL write_prior_rms_errors_of_trailing_vars(0)

                        DO k=1,nz
                                DO j=0,ny
                                        DO i=0,nx
                                                pre1_(i,j,k)=pre1(i,j,k)
                                                qbar_(i,j,k)=qbar(i,j,k)
                                                the1_(i,j,k)=the1(i,j,k)
                                                ubar_(i,j,k)=ubar(i,j,k)
                                                vbar_(i,j,k)=vbar(i,j,k)
                                                wbar_(i,j,k)=wbar(i,j,k)
                                                util_(i,j,k)=util(i,j,k)
                                                vtil_(i,j,k)=vtil(i,j,k)
                                                wtil_(i,j,k)=wtil(i,j,k)
                                                pre2_(i,j,k)=pre2(i,j,k)

                                                pre10_(i,j,k)=pre10(i,j,k)
                                                qbar0_(i,j,k)=qbar0(i,j,k)
                                                the10_(i,j,k)=the10(i,j,k)
                                                ubar0_(i,j,k)=ubar0(i,j,k)
                                                vbar0_(i,j,k)=vbar0(i,j,k)
                                                wbar0_(i,j,k)=wbar0(i,j,k)
                                                util0_(i,j,k)=util0(i,j,k)
                                                vtil0_(i,j,k)=vtil0(i,j,k)
                                                wtil0_(i,j,k)=wtil0(i,j,k)
                                                pre20_(i,j,k)=pre20(i,j,k)
                                        END DO
                                END DO
                        END DO

                        ! FLUSH: BEGIN

                        CALL set_psi2(ubar,vbar,wbar)
                        CALL set_phi2(pre1,qbar,the1,util,vtil,wtil)

                        CALL set_psi(time,pre1,qbar,the1,ubar,vbar,wbar,util,vtil,wtil)

                        DO k=1,nz
                                DO j=1,ny-1
                                        DO i=1,nx-1
                                                pre1(i,j,k)=varpsi(i,j,k)
                                        END DO
                                END DO
                        END DO

                        DO k=1,nz
                                DO j=1,ny-1
                                        pre1(0,j,k)=pre1(nx-1,j,k)
                                        pre1(nx,j,k)=pre1(1,j,k)
                                END DO

                                DO i=0,nx
                                        pre1(i,0,k)=pre1(i,ny-1,k)
                                        pre1(i,ny,k)=pre1(i,1,k)
                                END DO
                        END DO

                        CALL set_phi(delta,time,pre1,qbar,the1,ubar,vbar,wbar,util,vtil,wtil)

                        DO k=1,nz
                                DO j=1,ny-1
                                        DO i=1,nx-1
                                                pre2(i,j,k)=varphi(i,j,k)
                                        END DO
                                END DO
                        END DO

                        DO k=1,nz
                                DO j=1,ny-1
                                        pre2(0,j,k)=pre2(nx-1,j,k)
                                        pre2(nx,j,k)=pre2(1,j,k)
                                END DO

                                DO i=0,nx
                                        pre2(i,0,k)=pre2(i,ny-1,k)
                                        pre2(i,ny,k)=pre2(i,1,k)
                                END DO
                        END DO

                        CALL set_psi2(ubar0,vbar0,wbar0)
                        CALL set_phi2(pre10,qbar0,the10,util0,vtil0,wtil0)

                        CALL set_psi(time,pre10,qbar0,the10,ubar0,vbar0,wbar0,util0,vtil0,wtil0)

                        DO k=1,nz
                                DO j=1,ny-1
                                        DO i=1,nx-1
                                                pre10(i,j,k)=varpsi(i,j,k)
                                        END DO
                                END DO
                        END DO

                        DO k=1,nz
                                DO j=1,ny-1
                                        pre10(0,j,k)=pre10(nx-1,j,k)
                                        pre10(nx,j,k)=pre10(1,j,k)
                                END DO

                                DO i=0,nx
                                        pre10(i,0,k)=pre10(i,ny-1,k)
                                        pre10(i,ny,k)=pre10(i,1,k)
                                END DO
                        END DO

                        CALL set_phi(delta0,time,pre10,qbar0,the10,ubar0,vbar0,wbar0,util0,vtil0,wtil0)

                        DO k=1,nz
                                DO j=1,ny-1
                                        DO i=1,nx-1
                                                pre20(i,j,k)=varphi(i,j,k)
                                        END DO
                                END DO
                        END DO

                        DO k=1,nz
                                DO j=1,ny-1
                                        pre20(0,j,k)=pre20(nx-1,j,k)
                                        pre20(nx,j,k)=pre20(1,j,k)
                                END DO

                                DO i=0,nx
                                        pre20(i,0,k)=pre20(i,ny-1,k)
                                        pre20(i,ny,k)=pre20(i,1,k)
                                END DO
                        END DO

                        ! FLUSH: END

                        CALL set_precipitation(pre1,qbar,the1,wbar,wtil,prec)
                        CALL set_precipitation(pre10,qbar0,the10,wbar0,wtil0,prec0)
                        CALL write_prior_rms_errors_of_trailing_vars(1)

                        DO k=1,nz
                                DO j=0,ny
                                        DO i=0,nx
                                                pre1(i,j,k)=pre1_(i,j,k)
                                                qbar(i,j,k)=qbar_(i,j,k)
                                                the1(i,j,k)=the1_(i,j,k)
                                                ubar(i,j,k)=ubar_(i,j,k)
                                                vbar(i,j,k)=vbar_(i,j,k)
                                                wbar(i,j,k)=wbar_(i,j,k)
                                                util(i,j,k)=util_(i,j,k)
                                                vtil(i,j,k)=vtil_(i,j,k)
                                                wtil(i,j,k)=wtil_(i,j,k)
                                                pre2(i,j,k)=pre2_(i,j,k)

                                                pre10(i,j,k)=pre10_(i,j,k)
                                                qbar0(i,j,k)=qbar0_(i,j,k)
                                                the10(i,j,k)=the10_(i,j,k)
                                                ubar0(i,j,k)=ubar0_(i,j,k)
                                                vbar0(i,j,k)=vbar0_(i,j,k)
                                                wbar0(i,j,k)=wbar0_(i,j,k)
                                                util0(i,j,k)=util0_(i,j,k)
                                                vtil0(i,j,k)=vtil0_(i,j,k)
                                                wtil0(i,j,k)=wtil0_(i,j,k)
                                                pre20(i,j,k)=pre20_(i,j,k)
                                        END DO
                                END DO
                        END DO
                END IF

        END DO

        ! EVOLVE STATE: END

        ! WRITE PRIOR AND TRUE VALUES OF ALL VARS:

        IF (MOD(nn,print_lag_val)==0) THEN
                CALL set_precipitation(pre1,qbar,the1,wbar,wtil,prec)
                CALL set_precipitation(pre10,qbar0,the10,wbar0,wtil0,prec0)
                CALL write_prior_and_true_values_of_all_vars_and_time(0)

                DO k=1,nz
                        DO j=0,ny
                                DO i=0,nx
                                        pre1_(i,j,k)=pre1(i,j,k)
                                        qbar_(i,j,k)=qbar(i,j,k)
                                        the1_(i,j,k)=the1(i,j,k)
                                        ubar_(i,j,k)=ubar(i,j,k)
                                        vbar_(i,j,k)=vbar(i,j,k)
                                        wbar_(i,j,k)=wbar(i,j,k)
                                        util_(i,j,k)=util(i,j,k)
                                        vtil_(i,j,k)=vtil(i,j,k)
                                        wtil_(i,j,k)=wtil(i,j,k)
                                        pre2_(i,j,k)=pre2(i,j,k)

                                        pre10_(i,j,k)=pre10(i,j,k)
                                        qbar0_(i,j,k)=qbar0(i,j,k)
                                        the10_(i,j,k)=the10(i,j,k)
                                        ubar0_(i,j,k)=ubar0(i,j,k)
                                        vbar0_(i,j,k)=vbar0(i,j,k)
                                        wbar0_(i,j,k)=wbar0(i,j,k)
                                        util0_(i,j,k)=util0(i,j,k)
                                        vtil0_(i,j,k)=vtil0(i,j,k)
                                        wtil0_(i,j,k)=wtil0(i,j,k)
                                        pre20_(i,j,k)=pre20(i,j,k)
                                END DO
                        END DO
                END DO

                ! FLUSH: BEGIN

                CALL set_psi2(ubar,vbar,wbar)
                CALL set_phi2(pre1,qbar,the1,util,vtil,wtil)

                CALL set_psi(time,pre1,qbar,the1,ubar,vbar,wbar,util,vtil,wtil)

                DO k=1,nz
                        DO j=1,ny-1
                                DO i=1,nx-1
                                        pre1(i,j,k)=varpsi(i,j,k)
                                END DO
                        END DO
                END DO

                DO k=1,nz
                        DO j=1,ny-1
                                pre1(0,j,k)=pre1(nx-1,j,k)
                                pre1(nx,j,k)=pre1(1,j,k)
                        END DO

                        DO i=0,nx
                                pre1(i,0,k)=pre1(i,ny-1,k)
                                pre1(i,ny,k)=pre1(i,1,k)
                        END DO
                END DO

                CALL set_phi(delta,time,pre1,qbar,the1,ubar,vbar,wbar,util,vtil,wtil)

                DO k=1,nz
                        DO j=1,ny-1
                                DO i=1,nx-1
                                        pre2(i,j,k)=varphi(i,j,k)
                                END DO
                        END DO
                END DO

                DO k=1,nz
                        DO j=1,ny-1
                                pre2(0,j,k)=pre2(nx-1,j,k)
                                pre2(nx,j,k)=pre2(1,j,k)
                        END DO

                        DO i=0,nx
                                pre2(i,0,k)=pre2(i,ny-1,k)
                                pre2(i,ny,k)=pre2(i,1,k)
                        END DO
                END DO

                CALL set_psi2(ubar0,vbar0,wbar0)
                CALL set_phi2(pre10,qbar0,the10,util0,vtil0,wtil0)

                CALL set_psi(time,pre10,qbar0,the10,ubar0,vbar0,wbar0,util0,vtil0,wtil0)

                DO k=1,nz
                        DO j=1,ny-1
                                DO i=1,nx-1
                                        pre10(i,j,k)=varpsi(i,j,k)
                                END DO
                        END DO
                END DO

                DO k=1,nz
                        DO j=1,ny-1
                                pre10(0,j,k)=pre10(nx-1,j,k)
                                pre10(nx,j,k)=pre10(1,j,k)
                        END DO

                        DO i=0,nx
                                pre10(i,0,k)=pre10(i,ny-1,k)
                                pre10(i,ny,k)=pre10(i,1,k)
                        END DO
                END DO

                CALL set_phi(delta0,time,pre10,qbar0,the10,ubar0,vbar0,wbar0,util0,vtil0,wtil0)

                DO k=1,nz
                        DO j=1,ny-1
                                DO i=1,nx-1
                                        pre20(i,j,k)=varphi(i,j,k)
                                END DO
                        END DO
                END DO

                DO k=1,nz
                        DO j=1,ny-1
                                pre20(0,j,k)=pre20(nx-1,j,k)
                                pre20(nx,j,k)=pre20(1,j,k)
                        END DO

                        DO i=0,nx
                                pre20(i,0,k)=pre20(i,ny-1,k)
                                pre20(i,ny,k)=pre20(i,1,k)
                        END DO
                END DO

                ! FLUSH: END

                CALL set_precipitation(pre1,qbar,the1,wbar,wtil,prec)
                CALL set_precipitation(pre10,qbar0,the10,wbar0,wtil0,prec0)
                CALL write_prior_and_true_values_of_all_vars_and_time(1)

                DO k=1,nz
                        DO j=0,ny
                                DO i=0,nx
                                        pre1(i,j,k)=pre1_(i,j,k)
                                        qbar(i,j,k)=qbar_(i,j,k)
                                        the1(i,j,k)=the1_(i,j,k)
                                        ubar(i,j,k)=ubar_(i,j,k)
                                        vbar(i,j,k)=vbar_(i,j,k)
                                        wbar(i,j,k)=wbar_(i,j,k)
                                        util(i,j,k)=util_(i,j,k)
                                        vtil(i,j,k)=vtil_(i,j,k)
                                        wtil(i,j,k)=wtil_(i,j,k)
                                        pre2(i,j,k)=pre2_(i,j,k)

                                        pre10(i,j,k)=pre10_(i,j,k)
                                        qbar0(i,j,k)=qbar0_(i,j,k)
                                        the10(i,j,k)=the10_(i,j,k)
                                        ubar0(i,j,k)=ubar0_(i,j,k)
                                        vbar0(i,j,k)=vbar0_(i,j,k)
                                        wbar0(i,j,k)=wbar0_(i,j,k)
                                        util0(i,j,k)=util0_(i,j,k)
                                        vtil0(i,j,k)=vtil0_(i,j,k)
                                        wtil0(i,j,k)=wtil0_(i,j,k)
                                        pre20(i,j,k)=pre20_(i,j,k)
                                END DO
                        END DO
                END DO
        END IF

END DO

CONTAINS

SUBROUTINE check(status)

INTEGER,INTENT(IN) :: status

IF(status /= nf90_noerr) THEN
        PRINT *,trim(nf90_strerror(status))
        STOP "Stopped"
END IF

END SUBROUTINE check

SUBROUTINE input_resume_values_and_time()

CALL check(nf90_open(OUTPUT_FILE_NAME,NF90_NOWRITE,ncid))

CALL check(nf90_inq_varid(ncid,"pre1",pre1_varid))
CALL check(nf90_inq_varid(ncid,"qbar",qbar_varid))
CALL check(nf90_inq_varid(ncid,"the1",the1_varid))
CALL check(nf90_inq_varid(ncid,"ubar",ubar_varid))
CALL check(nf90_inq_varid(ncid,"vbar",vbar_varid))
CALL check(nf90_inq_varid(ncid,"wbar",wbar_varid))
CALL check(nf90_inq_varid(ncid,"util",util_varid))
CALL check(nf90_inq_varid(ncid,"vtil",vtil_varid))
CALL check(nf90_inq_varid(ncid,"wtil",wtil_varid))
CALL check(nf90_inq_varid(ncid,"pre2",pre2_varid))
CALL check(nf90_inq_varid(ncid,"prec",prec_varid))
CALL check(nf90_inq_varid(ncid,"pre10",pre10_varid))
CALL check(nf90_inq_varid(ncid,"qbar0",qbar0_varid))
CALL check(nf90_inq_varid(ncid,"the10",the10_varid))
CALL check(nf90_inq_varid(ncid,"ubar0",ubar0_varid))
CALL check(nf90_inq_varid(ncid,"vbar0",vbar0_varid))
CALL check(nf90_inq_varid(ncid,"wbar0",wbar0_varid))
CALL check(nf90_inq_varid(ncid,"util0",util0_varid))
CALL check(nf90_inq_varid(ncid,"vtil0",vtil0_varid))
CALL check(nf90_inq_varid(ncid,"wtil0",wtil0_varid))
CALL check(nf90_inq_varid(ncid,"pre20",pre20_varid))
CALL check(nf90_inq_varid(ncid,"prec0",prec0_varid))
CALL check(nf90_inq_varid(ncid,"time",time_varid))

CALL check(nf90_inq_dimid(ncid,"t",t_dimid))
CALL check(nf90_inquire_dimension(ncid,t_dimid,len=start(4)))

cc=start(4)+1

CALL check(nf90_get_var(ncid,pre1_varid,pre1,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,qbar_varid,qbar,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,the1_varid,the1,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,ubar_varid,ubar,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,vbar_varid,vbar,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,wbar_varid,wbar,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,util_varid,util,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,vtil_varid,vtil,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,wtil_varid,wtil,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,pre2_varid,pre2,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,prec_varid,prec,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,pre10_varid,pre10,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,qbar0_varid,qbar0,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,the10_varid,the10,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,ubar0_varid,ubar0,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,vbar0_varid,vbar0,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,wbar0_varid,wbar0,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,util0_varid,util0,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,vtil0_varid,vtil0,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,wtil0_varid,wtil0,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,pre20_varid,pre20,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,prec0_varid,prec0,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,time_varid,time,START=start(4:4)))

CALL check(nf90_close(ncid))

END SUBROUTINE input_resume_values_and_time

SUBROUTINE input_reference()

CALL check(nf90_open(OUTPUT_FILE_NAME,NF90_NOWRITE,ncid))

CALL check(nf90_inq_dimid(ncid,"t",t_dimid))
CALL check(nf90_inquire_dimension(ncid,t_dimid,len=start(4)))

CALL check(nf90_close(ncid))

CALL check(nf90_open(REFERENCE_FILE_NAME,NF90_NOWRITE,ncid))

CALL check(nf90_inq_varid(ncid,"pre10",pre10_varid))
CALL check(nf90_inq_varid(ncid,"qbar0",qbar0_varid))
CALL check(nf90_inq_varid(ncid,"the10",the10_varid))
CALL check(nf90_inq_varid(ncid,"ubar0",ubar0_varid))
CALL check(nf90_inq_varid(ncid,"vbar0",vbar0_varid))
CALL check(nf90_inq_varid(ncid,"wbar0",wbar0_varid))
CALL check(nf90_inq_varid(ncid,"util0",util0_varid))
CALL check(nf90_inq_varid(ncid,"vtil0",vtil0_varid))
CALL check(nf90_inq_varid(ncid,"wtil0",wtil0_varid))
CALL check(nf90_inq_varid(ncid,"pre20",pre20_varid))
CALL check(nf90_inq_varid(ncid,"prec0",prec0_varid))

! start(4)=1+(start(4)-1)*2

CALL check(nf90_get_var(ncid,pre10_varid,pre10,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,qbar0_varid,qbar0,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,the10_varid,the10,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,ubar0_varid,ubar0,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,vbar0_varid,vbar0,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,wbar0_varid,wbar0,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,util0_varid,util0,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,vtil0_varid,vtil0,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,wtil0_varid,wtil0,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,pre20_varid,pre20,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,prec0_varid,prec0,START=start,COUNT=count))

CALL check(nf90_close(ncid))

END SUBROUTINE input_reference

SUBROUTINE input_initial_sample_values()

CALL check(nf90_open(INPUT_SV_FILE_NAME,NF90_NOWRITE,ncid))

CALL check(nf90_inq_varid(ncid,"pre1",pre1_varid))
CALL check(nf90_inq_varid(ncid,"qbar",qbar_varid))
CALL check(nf90_inq_varid(ncid,"the1",the1_varid))
CALL check(nf90_inq_varid(ncid,"ubar",ubar_varid))
CALL check(nf90_inq_varid(ncid,"vbar",vbar_varid))
CALL check(nf90_inq_varid(ncid,"wbar",wbar_varid))
CALL check(nf90_inq_varid(ncid,"util",util_varid))
CALL check(nf90_inq_varid(ncid,"vtil",vtil_varid))
CALL check(nf90_inq_varid(ncid,"wtil",wtil_varid))
CALL check(nf90_inq_varid(ncid,"pre2",pre2_varid))
CALL check(nf90_inq_varid(ncid,"prec",prec_varid))

start(4)=sv_input_start
CALL check(nf90_get_var(ncid,pre1_varid,pre1,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,qbar_varid,qbar,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,the1_varid,the1,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,ubar_varid,ubar,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,vbar_varid,vbar,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,wbar_varid,wbar,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,util_varid,util,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,vtil_varid,vtil,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,wtil_varid,wtil,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,pre2_varid,pre2,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,prec_varid,prec,START=start,COUNT=count))

CALL check(nf90_close(ncid))

CALL check_borders()

END SUBROUTINE input_initial_sample_values

SUBROUTINE input_the1ref()

CALL check(nf90_open(THE1REF_FILE_NAME,NF90_NOWRITE,ncid))

CALL check(nf90_inq_varid(ncid,"the10",the1ref_varid))

start(4)=the1ref_input_start
CALL check(nf90_get_var(ncid,the1ref_varid,the1ref,START=start,COUNT=count))

CALL check(nf90_close(ncid))

CALL check_borders_the1ref()

END SUBROUTINE input_the1ref

SUBROUTINE input_the1ref0()

CALL check(nf90_open(THE1REF0_FILE_NAME,NF90_NOWRITE,ncid))

CALL check(nf90_inq_varid(ncid,"the10",the1ref0_varid))

start(4)=the1ref0_input_start
CALL check(nf90_get_var(ncid,the1ref0_varid,the1ref0,START=start,COUNT=count))

CALL check(nf90_close(ncid))

CALL check_borders_the1ref0()

END SUBROUTINE input_the1ref0

SUBROUTINE input_initial_true_values_and_time()

CALL check(nf90_open(INPUT_TV_FILE_NAME,NF90_NOWRITE,ncid))

CALL check(nf90_inq_varid(ncid,"pre10",pre10_varid))
CALL check(nf90_inq_varid(ncid,"qbar0",qbar0_varid))
CALL check(nf90_inq_varid(ncid,"the10",the10_varid))
CALL check(nf90_inq_varid(ncid,"ubar0",ubar0_varid))
CALL check(nf90_inq_varid(ncid,"vbar0",vbar0_varid))
CALL check(nf90_inq_varid(ncid,"wbar0",wbar0_varid))
CALL check(nf90_inq_varid(ncid,"util0",util0_varid))
CALL check(nf90_inq_varid(ncid,"vtil0",vtil0_varid))
CALL check(nf90_inq_varid(ncid,"wtil0",wtil0_varid))
CALL check(nf90_inq_varid(ncid,"pre20",pre20_varid))
CALL check(nf90_inq_varid(ncid,"prec0",prec0_varid))
CALL check(nf90_inq_varid(ncid,"time",time_varid))

start(4)=tv_input_start
CALL check(nf90_get_var(ncid,pre10_varid,pre10,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,qbar0_varid,qbar0,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,the10_varid,the10,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,ubar0_varid,ubar0,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,vbar0_varid,vbar0,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,wbar0_varid,wbar0,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,util0_varid,util0,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,vtil0_varid,vtil0,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,wtil0_varid,wtil0,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,pre20_varid,pre20,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,prec0_varid,prec0,START=start,COUNT=count))
CALL check(nf90_get_var(ncid,time_varid,time,START=start(4:4)))

CALL check(nf90_close(ncid))

CALL check_borders0()

END SUBROUTINE input_initial_true_values_and_time

SUBROUTINE create_and_write_prior_values_file(flushed_state)

INTEGER,INTENT(IN) :: flushed_state

IF (flushed_state==0) THEN
        CALL check(nf90_create(OUTPUT_FILE_NAME,NF90_CLOBBER,ncid))
ELSE
        CALL check(nf90_create(OUTPUT_FILE_NAME_FLUSHED,NF90_CLOBBER,ncid))
END IF

CALL check(nf90_def_dim(ncid,"x",nx+1,x_dimid))
CALL check(nf90_def_dim(ncid,"y",ny+1,y_dimid))
CALL check(nf90_def_dim(ncid,"z",nz,z_dimid))
CALL check(nf90_def_dim(ncid,"t",NF90_UNLIMITED,t_dimid))

dimids=(/ x_dimid,y_dimid,z_dimid,t_dimid /)

CALL check(nf90_def_var(ncid,"pre1",NF90_FLOAT,dimids,pre1_varid))
CALL check(nf90_def_var(ncid,"qbar",NF90_FLOAT,dimids,qbar_varid))
CALL check(nf90_def_var(ncid,"the1",NF90_FLOAT,dimids,the1_varid))
CALL check(nf90_def_var(ncid,"ubar",NF90_FLOAT,dimids,ubar_varid))
CALL check(nf90_def_var(ncid,"vbar",NF90_FLOAT,dimids,vbar_varid))
CALL check(nf90_def_var(ncid,"wbar",NF90_FLOAT,dimids,wbar_varid))
CALL check(nf90_def_var(ncid,"util",NF90_FLOAT,dimids,util_varid))
CALL check(nf90_def_var(ncid,"vtil",NF90_FLOAT,dimids,vtil_varid))
CALL check(nf90_def_var(ncid,"wtil",NF90_FLOAT,dimids,wtil_varid))
CALL check(nf90_def_var(ncid,"pre2",NF90_FLOAT,dimids,pre2_varid))
CALL check(nf90_def_var(ncid,"prec",NF90_FLOAT,dimids,prec_varid))
CALL check(nf90_def_var(ncid,"pre10",NF90_FLOAT,dimids,pre10_varid))
CALL check(nf90_def_var(ncid,"qbar0",NF90_FLOAT,dimids,qbar0_varid))
CALL check(nf90_def_var(ncid,"the10",NF90_FLOAT,dimids,the10_varid))
CALL check(nf90_def_var(ncid,"ubar0",NF90_FLOAT,dimids,ubar0_varid))
CALL check(nf90_def_var(ncid,"vbar0",NF90_FLOAT,dimids,vbar0_varid))
CALL check(nf90_def_var(ncid,"wbar0",NF90_FLOAT,dimids,wbar0_varid))
CALL check(nf90_def_var(ncid,"util0",NF90_FLOAT,dimids,util0_varid))
CALL check(nf90_def_var(ncid,"vtil0",NF90_FLOAT,dimids,vtil0_varid))
CALL check(nf90_def_var(ncid,"wtil0",NF90_FLOAT,dimids,wtil0_varid))
CALL check(nf90_def_var(ncid,"pre20",NF90_FLOAT,dimids,pre20_varid))
CALL check(nf90_def_var(ncid,"prec0",NF90_FLOAT,dimids,prec0_varid))
CALL check(nf90_def_var(ncid,"time",NF90_FLOAT,(/ t_dimid /),time_varid))

CALL check(nf90_enddef(ncid))

start(4)=cc
CALL check(nf90_put_var(ncid,pre1_varid,pre1,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,qbar_varid,qbar,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,the1_varid,the1,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,ubar_varid,ubar,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,vbar_varid,vbar,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,wbar_varid,wbar,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,util_varid,util,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,vtil_varid,vtil,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,wtil_varid,wtil,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,pre2_varid,pre2,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,prec_varid,prec,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,pre10_varid,pre10,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,qbar0_varid,qbar0,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,the10_varid,the10,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,ubar0_varid,ubar0,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,vbar0_varid,vbar0,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,wbar0_varid,wbar0,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,util0_varid,util0,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,vtil0_varid,vtil0,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,wtil0_varid,wtil0,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,pre20_varid,pre20,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,prec0_varid,prec0,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,time_varid,(/ time /),START=start(4:4),COUNT=count(4:4)))
IF (flushed_state==1) THEN
        cc=cc+1
END IF

CALL check(nf90_close(ncid))

END SUBROUTINE create_and_write_prior_values_file

SUBROUTINE write_prior_and_true_values_of_all_vars_and_time(flushed_state)

INTEGER,INTENT(IN) :: flushed_state

IF (flushed_state==0) THEN
        CALL check(nf90_open(OUTPUT_FILE_NAME,NF90_WRITE,ncid))
ELSE
        CALL check(nf90_open(OUTPUT_FILE_NAME_FLUSHED,NF90_WRITE,ncid))
END IF

CALL check(nf90_inq_varid(ncid,"pre1",pre1_varid))
CALL check(nf90_inq_varid(ncid,"qbar",qbar_varid))
CALL check(nf90_inq_varid(ncid,"the1",the1_varid))
CALL check(nf90_inq_varid(ncid,"ubar",ubar_varid))
CALL check(nf90_inq_varid(ncid,"vbar",vbar_varid))
CALL check(nf90_inq_varid(ncid,"wbar",wbar_varid))
CALL check(nf90_inq_varid(ncid,"util",util_varid))
CALL check(nf90_inq_varid(ncid,"vtil",vtil_varid))
CALL check(nf90_inq_varid(ncid,"wtil",wtil_varid))
CALL check(nf90_inq_varid(ncid,"pre2",pre2_varid))
CALL check(nf90_inq_varid(ncid,"prec",prec_varid))
CALL check(nf90_inq_varid(ncid,"pre10",pre10_varid))
CALL check(nf90_inq_varid(ncid,"qbar0",qbar0_varid))
CALL check(nf90_inq_varid(ncid,"the10",the10_varid))
CALL check(nf90_inq_varid(ncid,"ubar0",ubar0_varid))
CALL check(nf90_inq_varid(ncid,"vbar0",vbar0_varid))
CALL check(nf90_inq_varid(ncid,"wbar0",wbar0_varid))
CALL check(nf90_inq_varid(ncid,"util0",util0_varid))
CALL check(nf90_inq_varid(ncid,"vtil0",vtil0_varid))
CALL check(nf90_inq_varid(ncid,"wtil0",wtil0_varid))
CALL check(nf90_inq_varid(ncid,"pre20",pre20_varid))
CALL check(nf90_inq_varid(ncid,"prec0",prec0_varid))
CALL check(nf90_inq_varid(ncid,"time",time_varid))

start(4)=cc
CALL check(nf90_put_var(ncid,pre1_varid,pre1,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,qbar_varid,qbar,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,the1_varid,the1,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,ubar_varid,ubar,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,vbar_varid,vbar,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,wbar_varid,wbar,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,util_varid,util,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,vtil_varid,vtil,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,wtil_varid,wtil,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,pre2_varid,pre2,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,prec_varid,prec,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,pre10_varid,pre10,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,qbar0_varid,qbar0,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,the10_varid,the10,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,ubar0_varid,ubar0,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,vbar0_varid,vbar0,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,wbar0_varid,wbar0,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,util0_varid,util0,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,vtil0_varid,vtil0,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,wtil0_varid,wtil0,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,pre20_varid,pre20,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,prec0_varid,prec0,START=start,COUNT=count))
CALL check(nf90_put_var(ncid,time_varid,(/ time /),START=start(4:4),COUNT=count(4:4)))
IF (flushed_state==1) THEN
        cc=cc+1
END IF

CALL check(nf90_close(ncid))

END SUBROUTINE write_prior_and_true_values_of_all_vars_and_time

SUBROUTINE write_prior_rms_errors_of_trailing_vars(flushed_state)

INTEGER,INTENT(IN) :: flushed_state

INTEGER :: i,j,k

REAL :: err_pre1,err_qbar,err_the1
REAL :: err_ubar,err_vbar,err_wbar
REAL :: err_util,err_vtil,err_wtil
REAL :: err_prec
REAL :: err_wcli

err_pre1=0
err_qbar=0
err_the1=0
err_ubar=0
err_vbar=0
err_wbar=0
err_util=0
err_vtil=0
err_wtil=0
err_prec=0
err_wcli=0
DO k=1,nz
	DO j=1,ny-1
		DO i=1,nx-1
                        err_pre1=err_pre1+(pre1(i,j,k)-pre10(i,j,k))**2
                        err_qbar=err_qbar+(qbar(i,j,k)-qbar0(i,j,k))**2
                        err_the1=err_the1+(the1(i,j,k)-the10(i,j,k))**2
                        err_ubar=err_ubar+(ubar(i,j,k)-ubar0(i,j,k))**2
                        err_vbar=err_vbar+(vbar(i,j,k)-vbar0(i,j,k))**2
                        err_wbar=err_wbar+(wbar(i,j,k)-wbar0(i,j,k))**2
                        err_util=err_util+(util(i,j,k)-util0(i,j,k))**2
                        err_vtil=err_vtil+(vtil(i,j,k)-vtil0(i,j,k))**2
                        err_wtil=err_wtil+(wtil(i,j,k)-wtil0(i,j,k))**2
                        err_prec=err_prec+(prec(i,j,k)-prec0(i,j,k))**2
                        err_wcli=err_wcli+(wbar(i,j,k)+wtil(i,j,k)-wbar0(i,j,k)-wtil0(i,j,k))**2
		END DO
	END DO
END DO
err_pre1=SQRT(err_pre1/((nx-1)*(ny-1)*nz))
err_qbar=SQRT(err_qbar/((nx-1)*(ny-1)*nz))
err_the1=SQRT(err_the1/((nx-1)*(ny-1)*nz))
err_ubar=SQRT(err_ubar/((nx-1)*(ny-1)*nz))
err_vbar=SQRT(err_vbar/((nx-1)*(ny-1)*nz))
err_wbar=SQRT(err_wbar/((nx-1)*(ny-1)*nz))
err_util=SQRT(err_util/((nx-1)*(ny-1)*nz))
err_vtil=SQRT(err_vtil/((nx-1)*(ny-1)*nz))
err_wtil=SQRT(err_wtil/((nx-1)*(ny-1)*nz))
err_prec=SQRT(err_prec/((nx-1)*(ny-1)*nz))
err_wcli=SQRT(err_wcli/((nx-1)*(ny-1)*nz))

IF (flushed_state==0) THEN
        OPEN(UNIT=2,FILE='prior_err.dat',ACTION="write",STATUS="old",POSITION="append")
        WRITE(2,"(11ES15.4)") err_pre1,err_qbar,err_the1,err_ubar,err_vbar,err_wbar,err_util,err_vtil,err_wtil,err_prec,err_wcli
        CLOSE(2)
ELSE
        OPEN(UNIT=2,FILE='prior_err_flushed.dat',ACTION="write",STATUS="old",POSITION="append")
        WRITE(2,"(11ES15.4)") err_pre1,err_qbar,err_the1,err_ubar,err_vbar,err_wbar,err_util,err_vtil,err_wtil,err_prec,err_wcli
        CLOSE(2)
END IF

END SUBROUTINE write_prior_rms_errors_of_trailing_vars

SUBROUTINE set_precipitation(pre1,qbar,the1,wbar,wtil,prec)

REAL,DIMENSION(0:nx,0:ny,nz),INTENT(IN) :: pre1,qbar,the1
REAL,DIMENSION(0:nx,0:ny,nz),INTENT(IN) :: wbar
REAL,DIMENSION(0:nx,0:ny,nz),INTENT(IN) :: wtil

REAL,DIMENSION(0:nx,0:ny,nz),INTENT(OUT) :: prec

INTEGER :: i,j,k

DO k=1,nz
        DO j=1,ny-1
                DO i=1,nx-1
prec(i,j,k)=precip(i,j,k,pre1,wbar,the1,qbar,wtil)
                END DO
        END DO
END DO

DO k=1,nz
        DO j=1,ny-1
                prec(0,j,k)=prec(nx-1,j,k)
                prec(nx,j,k)=prec(1,j,k)
        END DO

        DO i=0,nx
                prec(i,0,k)=prec(i,ny-1,k)
                prec(i,ny,k)=prec(i,1,k)
        END DO
END DO

END SUBROUTINE set_precipitation

END PROGRAM main
