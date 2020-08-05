MODULE evolve0

USE state0, ONLY: pre10,qbar0,the10
USE state0, ONLY: ubar0,vbar0,wbar0
USE state0, ONLY: util0,vtil0,wtil0
USE state0, ONLY: pre20,the1ref0

USE auxiliaries, ONLY: pre1_,qbar_,the1_
USE auxiliaries, ONLY: ubar_,vbar_,wbar_
USE auxiliaries, ONLY: util_,vtil_,wtil_

USE auxiliaries, ONLY: ffqbar,ffthe1
USE auxiliaries, ONLY: ffubar,ffvbar,ffwbar
USE auxiliaries, ONLY: ffrhoutil,ffrhovtil,ffrhowtil

USE auxiliaries, ONLY: ffqbar1,ffthe11
USE auxiliaries, ONLY: ffubar1,ffvbar1,ffwbar1
USE auxiliaries, ONLY: ffrhoutil1,ffrhovtil1,ffrhowtil1

USE auxiliaries, ONLY: ffqbar2,ffthe12
USE auxiliaries, ONLY: ffubar2,ffvbar2,ffwbar2
USE auxiliaries, ONLY: ffrhoutil2,ffrhovtil2,ffrhowtil2

USE auxiliaries, ONLY: ffqbar3,ffthe13
USE auxiliaries, ONLY: ffubar3,ffvbar3,ffwbar3
USE auxiliaries, ONLY: ffrhoutil3,ffrhovtil3,ffrhowtil3

!USE auxiliaries, ONLY: ffubar_,ffvbar_,ffwbar_
!USE auxiliaries, ONLY: ffrhoutil_,ffrhovtil_,ffrhowtil_

!USE auxiliaries, ONLY: ffubar1_,ffvbar1_,ffwbar1_
!USE auxiliaries, ONLY: ffrhoutil1_,ffrhovtil1_,ffrhowtil1_

!USE auxiliaries, ONLY: ffubar2_,ffvbar2_,ffwbar2_
!USE auxiliaries, ONLY: ffrhoutil2_,ffrhovtil2_,ffrhowtil2_

!USE auxiliaries, ONLY: ffubar3_,ffvbar3_,ffwbar3_
!USE auxiliaries, ONLY: ffrhoutil3_,ffrhovtil3_,ffrhowtil3_

USE psi_solver, ONLY: varpsi
USE phi_solver, ONLY: varphi

USE parameters, ONLY: nx,ny,nz
USE parameters, ONLY: dx,dy,dz
USE parameters, ONLY: dt
USE parameters, ONLY: max_ubar,max_vbar,max_wbar
USE parameters, ONLY: max_util,max_vtil,max_wtil
USE parameters, ONLY: delta0

USE utilities, ONLY: rho0,rho1

USE psi_solver, ONLY: set_psi2,set_psi4
USE phi_solver, ONLY: set_phi,set_phi2

USE field, ONLY: set_field1

IMPLICIT NONE
PUBLIC

CONTAINS

SUBROUTINE evolve_state0(time)

REAL,INTENT(IN) :: time

CALL rk_projection_method(time)

END SUBROUTINE evolve_state0

SUBROUTINE rk_projection_method(time)

REAL,INTENT(IN) :: time

INTEGER :: i,j,k

CALL set_field1(delta0,time,pre10,qbar0,the10,ubar0,vbar0,wbar0,util0,vtil0,wtil0,pre20,the1ref0,                                &
                ffqbar,ffthe1,ffubar,ffvbar,ffwbar,ffrhoutil,ffrhovtil,ffrhowtil)

!CALL set_field2(delta0,time,pre10,qbar0,the10,ubar0,vbar0,wbar0,util0,vtil0,wtil0,                                               &
                !ffubar_,ffvbar_,ffwbar_,ffrhoutil_,ffrhovtil_,ffrhowtil_)

DO k=1,nz
        DO j=1,ny-1
                DO i=1,nx-1
			qbar_(i,j,k)=qbar0(i,j,k)+ffqbar(i,j,k)*dt/2
			the1_(i,j,k)=the10(i,j,k)+ffthe1(i,j,k)*dt/2
			ubar_(i,j,k)=ubar0(i,j,k)+ffubar(i,j,k)*dt/2
                        IF (ubar_(i,j,k)>0) THEN
                                ubar_(i,j,k)=MIN(ubar_(i,j,k),max_ubar)
                        ELSE
                                ubar_(i,j,k)=MAX(ubar_(i,j,k),-max_ubar)
                        END IF
			vbar_(i,j,k)=vbar0(i,j,k)+ffvbar(i,j,k)*dt/2
                        IF (vbar_(i,j,k)>0) THEN
                                vbar_(i,j,k)=MIN(vbar_(i,j,k),max_vbar)
                        ELSE
                                vbar_(i,j,k)=MAX(vbar_(i,j,k),-max_vbar)
                        END IF
			wbar_(i,j,k)=wbar0(i,j,k)+ffwbar(i,j,k)*dt/2
                        IF (wbar_(i,j,k)>0) THEN
                                wbar_(i,j,k)=MIN(wbar_(i,j,k),max_wbar)
                        ELSE
                                wbar_(i,j,k)=MAX(wbar_(i,j,k),-max_wbar)
                        END IF
			util_(i,j,k)=(rho0(k)+rho1(i,j,k,pre10,qbar0,the10))*util0(i,j,k)+ffrhoutil(i,j,k)*dt/2
                        IF (util_(i,j,k)>0) THEN
                                util_(i,j,k)=MIN(util_(i,j,k),(rho0(k)+rho1(i,j,k,pre10,qbar0,the10))*max_util)
                        ELSE
                                util_(i,j,k)=MAX(util_(i,j,k),-(rho0(k)+rho1(i,j,k,pre10,qbar0,the10))*max_util)
                        END IF
			vtil_(i,j,k)=(rho0(k)+rho1(i,j,k,pre10,qbar0,the10))*vtil0(i,j,k)+ffrhovtil(i,j,k)*dt/2
                        IF (vtil_(i,j,k)>0) THEN
                                vtil_(i,j,k)=MIN(vtil_(i,j,k),(rho0(k)+rho1(i,j,k,pre10,qbar0,the10))*max_vtil)
                        ELSE
                                vtil_(i,j,k)=MAX(vtil_(i,j,k),-(rho0(k)+rho1(i,j,k,pre10,qbar0,the10))*max_vtil)
                        END IF
			wtil_(i,j,k)=(rho0(k)+rho1(i,j,k,pre10,qbar0,the10))*wtil0(i,j,k)+ffrhowtil(i,j,k)*dt/2
                        IF (wtil_(i,j,k)>0) THEN
                                wtil_(i,j,k)=MIN(wtil_(i,j,k),(rho0(k)+rho1(i,j,k,pre10,qbar0,the10))*max_wtil)
                        ELSE
                                wtil_(i,j,k)=MAX(wtil_(i,j,k),-(rho0(k)+rho1(i,j,k,pre10,qbar0,the10))*max_wtil)
                        END IF
		END DO
	END DO
END DO

DO k=1,nz
        DO j=1,ny-1
                qbar_(0,j,k)=qbar_(nx-1,j,k)
                the1_(0,j,k)=the1_(nx-1,j,k)
                ubar_(0,j,k)=ubar_(nx-1,j,k)
                vbar_(0,j,k)=vbar_(nx-1,j,k)
                wbar_(0,j,k)=wbar_(nx-1,j,k)
                util_(0,j,k)=util_(nx-1,j,k)
                vtil_(0,j,k)=vtil_(nx-1,j,k)
                wtil_(0,j,k)=wtil_(nx-1,j,k)

                qbar_(nx,j,k)=qbar_(1,j,k)
                the1_(nx,j,k)=the1_(1,j,k)
                ubar_(nx,j,k)=ubar_(1,j,k)
                vbar_(nx,j,k)=vbar_(1,j,k)
                wbar_(nx,j,k)=wbar_(1,j,k)
                util_(nx,j,k)=util_(1,j,k)
                vtil_(nx,j,k)=vtil_(1,j,k)
                wtil_(nx,j,k)=wtil_(1,j,k)
        END DO

        DO i=0,nx
                qbar_(i,0,k)=qbar_(i,ny-1,k)
                the1_(i,0,k)=the1_(i,ny-1,k)
                ubar_(i,0,k)=ubar_(i,ny-1,k)
                vbar_(i,0,k)=vbar_(i,ny-1,k)
                wbar_(i,0,k)=wbar_(i,ny-1,k)
                util_(i,0,k)=util_(i,ny-1,k)
                vtil_(i,0,k)=vtil_(i,ny-1,k)
                wtil_(i,0,k)=wtil_(i,ny-1,k)

                qbar_(i,ny,k)=qbar_(i,1,k)
                the1_(i,ny,k)=the1_(i,1,k)
                ubar_(i,ny,k)=ubar_(i,1,k)
                vbar_(i,ny,k)=vbar_(i,1,k)
                wbar_(i,ny,k)=wbar_(i,1,k)
                util_(i,ny,k)=util_(i,1,k)
                vtil_(i,ny,k)=vtil_(i,1,k)
                wtil_(i,ny,k)=wtil_(i,1,k)
        END DO
END DO

CALL set_psi2(ubar_,vbar_,wbar_)

CALL set_psi4(time+dt/2,qbar_,the1_,ubar_,vbar_,wbar_,util_,vtil_,wtil_)

DO k=1,nz
	DO j=1,ny-1
		DO i=1,nx-1
			pre1_(i,j,k)=varpsi(i,j,k)
		END DO
	END DO
END DO

DO k=1,nz
        DO j=1,ny-1
                pre1_(0,j,k)=pre1_(nx-1,j,k)
                pre1_(nx,j,k)=pre1_(1,j,k)
        END DO

        DO i=0,nx
                pre1_(i,0,k)=pre1_(i,ny-1,k)
                pre1_(i,ny,k)=pre1_(i,1,k)
        END DO
END DO

DO k=1,nz
        DO j=1,ny-1
                DO i=1,nx-1
			util_(i,j,k)=util_(i,j,k)/(rho0(k)+rho1(i,j,k,pre1_,qbar_,the1_))
			vtil_(i,j,k)=vtil_(i,j,k)/(rho0(k)+rho1(i,j,k,pre1_,qbar_,the1_))
			wtil_(i,j,k)=wtil_(i,j,k)/(rho0(k)+rho1(i,j,k,pre1_,qbar_,the1_))
		END DO
	END DO
END DO

DO k=1,nz
        DO j=1,ny-1
                util_(0,j,k)=util_(nx-1,j,k)
                vtil_(0,j,k)=vtil_(nx-1,j,k)
                wtil_(0,j,k)=wtil_(nx-1,j,k)

                util_(nx,j,k)=util_(1,j,k)
                vtil_(nx,j,k)=vtil_(1,j,k)
                wtil_(nx,j,k)=wtil_(1,j,k)
        END DO

        DO i=0,nx
                util_(i,0,k)=util_(i,ny-1,k)
                vtil_(i,0,k)=vtil_(i,ny-1,k)
                wtil_(i,0,k)=wtil_(i,ny-1,k)

                util_(i,ny,k)=util_(i,1,k)
                vtil_(i,ny,k)=vtil_(i,1,k)
                wtil_(i,ny,k)=wtil_(i,1,k)
        END DO
END DO

CALL set_phi2(pre1_,qbar_,the1_,util_,vtil_,wtil_)

CALL set_phi(delta0,time+dt/2,pre1_,qbar_,the1_,ubar_,vbar_,wbar_,util_,vtil_,wtil_)

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

CALL set_field1(delta0,time+dt/2,pre1_,qbar_,the1_,ubar_,vbar_,wbar_,util_,vtil_,wtil_,pre20,the1ref0,                           &
                   ffqbar1,ffthe11,ffubar1,ffvbar1,ffwbar1,ffrhoutil1,ffrhovtil1,ffrhowtil1)

!CALL set_field2(delta0,time+dt/2,pre1_,qbar_,the1_,ubar_,vbar_,wbar_,util_,vtil_,wtil_,                                          &
                      !ffubar1_,ffvbar1_,ffwbar1_,ffrhoutil1_,ffrhovtil1_,ffrhowtil1_)

DO k=1,nz
        DO j=1,ny-1
                DO i=1,nx-1
			qbar_(i,j,k)=qbar0(i,j,k)+ffqbar1(i,j,k)*dt/2
			the1_(i,j,k)=the10(i,j,k)+ffthe11(i,j,k)*dt/2
			ubar_(i,j,k)=ubar0(i,j,k)+ffubar1(i,j,k)*dt/2
                        IF (ubar_(i,j,k)>0) THEN
                                ubar_(i,j,k)=MIN(ubar_(i,j,k),max_ubar)
                        ELSE
                                ubar_(i,j,k)=MAX(ubar_(i,j,k),-max_ubar)
                        END IF
			vbar_(i,j,k)=vbar0(i,j,k)+ffvbar1(i,j,k)*dt/2
                        IF (vbar_(i,j,k)>0) THEN
                                vbar_(i,j,k)=MIN(vbar_(i,j,k),max_vbar)
                        ELSE
                                vbar_(i,j,k)=MAX(vbar_(i,j,k),-max_vbar)
                        END IF
			wbar_(i,j,k)=wbar0(i,j,k)+ffwbar1(i,j,k)*dt/2
                        IF (wbar_(i,j,k)>0) THEN
                                wbar_(i,j,k)=MIN(wbar_(i,j,k),max_wbar)
                        ELSE
                                wbar_(i,j,k)=MAX(wbar_(i,j,k),-max_wbar)
                        END IF
			util_(i,j,k)=(rho0(k)+rho1(i,j,k,pre10,qbar0,the10))*util0(i,j,k)+ffrhoutil1(i,j,k)*dt/2
                        IF (util_(i,j,k)>0) THEN
                                util_(i,j,k)=MIN(util_(i,j,k),(rho0(k)+rho1(i,j,k,pre10,qbar0,the10))*max_util)
                        ELSE
                                util_(i,j,k)=MAX(util_(i,j,k),-(rho0(k)+rho1(i,j,k,pre10,qbar0,the10))*max_util)
                        END IF
			vtil_(i,j,k)=(rho0(k)+rho1(i,j,k,pre10,qbar0,the10))*vtil0(i,j,k)+ffrhovtil1(i,j,k)*dt/2
                        IF (vtil_(i,j,k)>0) THEN
                                vtil_(i,j,k)=MIN(vtil_(i,j,k),(rho0(k)+rho1(i,j,k,pre10,qbar0,the10))*max_vtil)
                        ELSE
                                vtil_(i,j,k)=MAX(vtil_(i,j,k),-(rho0(k)+rho1(i,j,k,pre10,qbar0,the10))*max_vtil)
                        END IF
			wtil_(i,j,k)=(rho0(k)+rho1(i,j,k,pre10,qbar0,the10))*wtil0(i,j,k)+ffrhowtil1(i,j,k)*dt/2
                        IF (wtil_(i,j,k)>0) THEN
                                wtil_(i,j,k)=MIN(wtil_(i,j,k),(rho0(k)+rho1(i,j,k,pre10,qbar0,the10))*max_wtil)
                        ELSE
                                wtil_(i,j,k)=MAX(wtil_(i,j,k),-(rho0(k)+rho1(i,j,k,pre10,qbar0,the10))*max_wtil)
                        END IF
		END DO
	END DO
END DO

DO k=1,nz
        DO j=1,ny-1
                qbar_(0,j,k)=qbar_(nx-1,j,k)
                the1_(0,j,k)=the1_(nx-1,j,k)
                ubar_(0,j,k)=ubar_(nx-1,j,k)
                vbar_(0,j,k)=vbar_(nx-1,j,k)
                wbar_(0,j,k)=wbar_(nx-1,j,k)
                util_(0,j,k)=util_(nx-1,j,k)
                vtil_(0,j,k)=vtil_(nx-1,j,k)
                wtil_(0,j,k)=wtil_(nx-1,j,k)

                qbar_(nx,j,k)=qbar_(1,j,k)
                the1_(nx,j,k)=the1_(1,j,k)
                ubar_(nx,j,k)=ubar_(1,j,k)
                vbar_(nx,j,k)=vbar_(1,j,k)
                wbar_(nx,j,k)=wbar_(1,j,k)
                util_(nx,j,k)=util_(1,j,k)
                vtil_(nx,j,k)=vtil_(1,j,k)
                wtil_(nx,j,k)=wtil_(1,j,k)
        END DO

        DO i=0,nx
                qbar_(i,0,k)=qbar_(i,ny-1,k)
                the1_(i,0,k)=the1_(i,ny-1,k)
                ubar_(i,0,k)=ubar_(i,ny-1,k)
                vbar_(i,0,k)=vbar_(i,ny-1,k)
                wbar_(i,0,k)=wbar_(i,ny-1,k)
                util_(i,0,k)=util_(i,ny-1,k)
                vtil_(i,0,k)=vtil_(i,ny-1,k)
                wtil_(i,0,k)=wtil_(i,ny-1,k)

                qbar_(i,ny,k)=qbar_(i,1,k)
                the1_(i,ny,k)=the1_(i,1,k)
                ubar_(i,ny,k)=ubar_(i,1,k)
                vbar_(i,ny,k)=vbar_(i,1,k)
                wbar_(i,ny,k)=wbar_(i,1,k)
                util_(i,ny,k)=util_(i,1,k)
                vtil_(i,ny,k)=vtil_(i,1,k)
                wtil_(i,ny,k)=wtil_(i,1,k)
        END DO
END DO

CALL set_psi2(ubar_,vbar_,wbar_)

CALL set_psi4(time+dt/2,qbar_,the1_,ubar_,vbar_,wbar_,util_,vtil_,wtil_)

DO k=1,nz
	DO j=1,ny-1
		DO i=1,nx-1
			pre1_(i,j,k)=varpsi(i,j,k)
		END DO
	END DO
END DO

DO k=1,nz
        DO j=1,ny-1
                pre1_(0,j,k)=pre1_(nx-1,j,k)
                pre1_(nx,j,k)=pre1_(1,j,k)
        END DO

        DO i=0,nx
                pre1_(i,0,k)=pre1_(i,ny-1,k)
                pre1_(i,ny,k)=pre1_(i,1,k)
        END DO
END DO

DO k=1,nz
        DO j=1,ny-1
                DO i=1,nx-1
			util_(i,j,k)=util_(i,j,k)/(rho0(k)+rho1(i,j,k,pre1_,qbar_,the1_))
			vtil_(i,j,k)=vtil_(i,j,k)/(rho0(k)+rho1(i,j,k,pre1_,qbar_,the1_))
			wtil_(i,j,k)=wtil_(i,j,k)/(rho0(k)+rho1(i,j,k,pre1_,qbar_,the1_))
		END DO
	END DO
END DO

DO k=1,nz
        DO j=1,ny-1
                util_(0,j,k)=util_(nx-1,j,k)
                vtil_(0,j,k)=vtil_(nx-1,j,k)
                wtil_(0,j,k)=wtil_(nx-1,j,k)

                util_(nx,j,k)=util_(1,j,k)
                vtil_(nx,j,k)=vtil_(1,j,k)
                wtil_(nx,j,k)=wtil_(1,j,k)
        END DO

        DO i=0,nx
                util_(i,0,k)=util_(i,ny-1,k)
                vtil_(i,0,k)=vtil_(i,ny-1,k)
                wtil_(i,0,k)=wtil_(i,ny-1,k)

                util_(i,ny,k)=util_(i,1,k)
                vtil_(i,ny,k)=vtil_(i,1,k)
                wtil_(i,ny,k)=wtil_(i,1,k)
        END DO
END DO

CALL set_phi2(pre1_,qbar_,the1_,util_,vtil_,wtil_)

CALL set_phi(delta0,time+dt/2,pre1_,qbar_,the1_,ubar_,vbar_,wbar_,util_,vtil_,wtil_)

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

CALL set_field1(delta0,time+dt/2,pre1_,qbar_,the1_,ubar_,vbar_,wbar_,util_,vtil_,wtil_,pre20,the1ref0,                           &
                   ffqbar2,ffthe12,ffubar2,ffvbar2,ffwbar2,ffrhoutil2,ffrhovtil2,ffrhowtil2)

!CALL set_field2(delta0,time+dt/2,pre1_,qbar_,the1_,ubar_,vbar_,wbar_,util_,vtil_,wtil_,                                          &
                      !ffubar2_,ffvbar2_,ffwbar2_,ffrhoutil2_,ffrhovtil2_,ffrhowtil2_)

DO k=1,nz
        DO j=1,ny-1
                DO i=1,nx-1
			qbar_(i,j,k)=qbar0(i,j,k)+ffqbar2(i,j,k)*dt
			the1_(i,j,k)=the10(i,j,k)+ffthe12(i,j,k)*dt
			ubar_(i,j,k)=ubar0(i,j,k)+ffubar2(i,j,k)*dt
                        IF (ubar_(i,j,k)>0) THEN
                                ubar_(i,j,k)=MIN(ubar_(i,j,k),max_ubar)
                        ELSE
                                ubar_(i,j,k)=MAX(ubar_(i,j,k),-max_ubar)
                        END IF
			vbar_(i,j,k)=vbar0(i,j,k)+ffvbar2(i,j,k)*dt
                        IF (vbar_(i,j,k)>0) THEN
                                vbar_(i,j,k)=MIN(vbar_(i,j,k),max_vbar)
                        ELSE
                                vbar_(i,j,k)=MAX(vbar_(i,j,k),-max_vbar)
                        END IF
			wbar_(i,j,k)=wbar0(i,j,k)+ffwbar2(i,j,k)*dt
                        IF (wbar_(i,j,k)>0) THEN
                                wbar_(i,j,k)=MIN(wbar_(i,j,k),max_wbar)
                        ELSE
                                wbar_(i,j,k)=MAX(wbar_(i,j,k),-max_wbar)
                        END IF
			util_(i,j,k)=(rho0(k)+rho1(i,j,k,pre10,qbar0,the10))*util0(i,j,k)+ffrhoutil2(i,j,k)*dt
                        IF (util_(i,j,k)>0) THEN
                                util_(i,j,k)=MIN(util_(i,j,k),(rho0(k)+rho1(i,j,k,pre10,qbar0,the10))*max_util)
                        ELSE
                                util_(i,j,k)=MAX(util_(i,j,k),-(rho0(k)+rho1(i,j,k,pre10,qbar0,the10))*max_util)
                        END IF
			vtil_(i,j,k)=(rho0(k)+rho1(i,j,k,pre10,qbar0,the10))*vtil0(i,j,k)+ffrhovtil2(i,j,k)*dt
                        IF (vtil_(i,j,k)>0) THEN
                                vtil_(i,j,k)=MIN(vtil_(i,j,k),(rho0(k)+rho1(i,j,k,pre10,qbar0,the10))*max_vtil)
                        ELSE
                                vtil_(i,j,k)=MAX(vtil_(i,j,k),-(rho0(k)+rho1(i,j,k,pre10,qbar0,the10))*max_vtil)
                        END IF
			wtil_(i,j,k)=(rho0(k)+rho1(i,j,k,pre10,qbar0,the10))*wtil0(i,j,k)+ffrhowtil2(i,j,k)*dt
                        IF (wtil_(i,j,k)>0) THEN
                                wtil_(i,j,k)=MIN(wtil_(i,j,k),(rho0(k)+rho1(i,j,k,pre10,qbar0,the10))*max_wtil)
                        ELSE
                                wtil_(i,j,k)=MAX(wtil_(i,j,k),-(rho0(k)+rho1(i,j,k,pre10,qbar0,the10))*max_wtil)
                        END IF
		END DO
	END DO
END DO

DO k=1,nz
        DO j=1,ny-1
                qbar_(0,j,k)=qbar_(nx-1,j,k)
                the1_(0,j,k)=the1_(nx-1,j,k)
                ubar_(0,j,k)=ubar_(nx-1,j,k)
                vbar_(0,j,k)=vbar_(nx-1,j,k)
                wbar_(0,j,k)=wbar_(nx-1,j,k)
                util_(0,j,k)=util_(nx-1,j,k)
                vtil_(0,j,k)=vtil_(nx-1,j,k)
                wtil_(0,j,k)=wtil_(nx-1,j,k)

                qbar_(nx,j,k)=qbar_(1,j,k)
                the1_(nx,j,k)=the1_(1,j,k)
                ubar_(nx,j,k)=ubar_(1,j,k)
                vbar_(nx,j,k)=vbar_(1,j,k)
                wbar_(nx,j,k)=wbar_(1,j,k)
                util_(nx,j,k)=util_(1,j,k)
                vtil_(nx,j,k)=vtil_(1,j,k)
                wtil_(nx,j,k)=wtil_(1,j,k)
        END DO

        DO i=0,nx
                qbar_(i,0,k)=qbar_(i,ny-1,k)
                the1_(i,0,k)=the1_(i,ny-1,k)
                ubar_(i,0,k)=ubar_(i,ny-1,k)
                vbar_(i,0,k)=vbar_(i,ny-1,k)
                wbar_(i,0,k)=wbar_(i,ny-1,k)
                util_(i,0,k)=util_(i,ny-1,k)
                vtil_(i,0,k)=vtil_(i,ny-1,k)
                wtil_(i,0,k)=wtil_(i,ny-1,k)

                qbar_(i,ny,k)=qbar_(i,1,k)
                the1_(i,ny,k)=the1_(i,1,k)
                ubar_(i,ny,k)=ubar_(i,1,k)
                vbar_(i,ny,k)=vbar_(i,1,k)
                wbar_(i,ny,k)=wbar_(i,1,k)
                util_(i,ny,k)=util_(i,1,k)
                vtil_(i,ny,k)=vtil_(i,1,k)
                wtil_(i,ny,k)=wtil_(i,1,k)
        END DO
END DO

CALL set_psi2(ubar_,vbar_,wbar_)

CALL set_psi4(time+dt,qbar_,the1_,ubar_,vbar_,wbar_,util_,vtil_,wtil_)

DO k=1,nz
	DO j=1,ny-1
		DO i=1,nx-1
			pre1_(i,j,k)=varpsi(i,j,k)
		END DO
	END DO
END DO

DO k=1,nz
        DO j=1,ny-1
                pre1_(0,j,k)=pre1_(nx-1,j,k)
                pre1_(nx,j,k)=pre1_(1,j,k)
        END DO

        DO i=0,nx
                pre1_(i,0,k)=pre1_(i,ny-1,k)
                pre1_(i,ny,k)=pre1_(i,1,k)
        END DO
END DO

DO k=1,nz
        DO j=1,ny-1
                DO i=1,nx-1
			util_(i,j,k)=util_(i,j,k)/(rho0(k)+rho1(i,j,k,pre1_,qbar_,the1_))
			vtil_(i,j,k)=vtil_(i,j,k)/(rho0(k)+rho1(i,j,k,pre1_,qbar_,the1_))
			wtil_(i,j,k)=wtil_(i,j,k)/(rho0(k)+rho1(i,j,k,pre1_,qbar_,the1_))
		END DO
	END DO
END DO

DO k=1,nz
        DO j=1,ny-1
                util_(0,j,k)=util_(nx-1,j,k)
                vtil_(0,j,k)=vtil_(nx-1,j,k)
                wtil_(0,j,k)=wtil_(nx-1,j,k)

                util_(nx,j,k)=util_(1,j,k)
                vtil_(nx,j,k)=vtil_(1,j,k)
                wtil_(nx,j,k)=wtil_(1,j,k)
        END DO

        DO i=0,nx
                util_(i,0,k)=util_(i,ny-1,k)
                vtil_(i,0,k)=vtil_(i,ny-1,k)
                wtil_(i,0,k)=wtil_(i,ny-1,k)

                util_(i,ny,k)=util_(i,1,k)
                vtil_(i,ny,k)=vtil_(i,1,k)
                wtil_(i,ny,k)=wtil_(i,1,k)
        END DO
END DO

CALL set_phi2(pre1_,qbar_,the1_,util_,vtil_,wtil_)

CALL set_phi(delta0,time+dt,pre1_,qbar_,the1_,ubar_,vbar_,wbar_,util_,vtil_,wtil_)

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

CALL set_field1(delta0,time+dt,pre1_,qbar_,the1_,ubar_,vbar_,wbar_,util_,vtil_,wtil_,pre20,the1ref0,                             &
                 ffqbar3,ffthe13,ffubar3,ffvbar3,ffwbar3,ffrhoutil3,ffrhovtil3,ffrhowtil3)

!CALL set_field2(delta0,time+dt,pre1_,qbar_,the1_,ubar_,vbar_,wbar_,util_,vtil_,wtil_,                                            &
                    !ffubar3_,ffvbar3_,ffwbar3_,ffrhoutil3_,ffrhovtil3_,ffrhowtil3_)

DO k=1,nz
        DO j=1,ny-1
                DO i=1,nx-1
			qbar_(i,j,k)=qbar0(i,j,k)+dt*(ffqbar(i,j,k)+2*ffqbar1(i,j,k)+2*ffqbar2(i,j,k)+ffqbar3(i,j,k))/6
			the1_(i,j,k)=the10(i,j,k)+dt*(ffthe1(i,j,k)+2*ffthe11(i,j,k)+2*ffthe12(i,j,k)+ffthe13(i,j,k))/6
			ubar_(i,j,k)=ubar0(i,j,k)+dt*(ffubar(i,j,k)+2*ffubar1(i,j,k)+2*ffubar2(i,j,k)+ffubar3(i,j,k))/6
                        IF (ubar_(i,j,k)>0) THEN
                                ubar_(i,j,k)=MIN(ubar_(i,j,k),max_ubar)
                        ELSE
                                ubar_(i,j,k)=MAX(ubar_(i,j,k),-max_ubar)
                        END IF
			vbar_(i,j,k)=vbar0(i,j,k)+dt*(ffvbar(i,j,k)+2*ffvbar1(i,j,k)+2*ffvbar2(i,j,k)+ffvbar3(i,j,k))/6
                        IF (vbar_(i,j,k)>0) THEN
                                vbar_(i,j,k)=MIN(vbar_(i,j,k),max_vbar)
                        ELSE
                                vbar_(i,j,k)=MAX(vbar_(i,j,k),-max_vbar)
                        END IF
			wbar_(i,j,k)=wbar0(i,j,k)+dt*(ffwbar(i,j,k)+2*ffwbar1(i,j,k)+2*ffwbar2(i,j,k)+ffwbar3(i,j,k))/6
                        IF (wbar_(i,j,k)>0) THEN
                                wbar_(i,j,k)=MIN(wbar_(i,j,k),max_wbar)
                        ELSE
                                wbar_(i,j,k)=MAX(wbar_(i,j,k),-max_wbar)
                        END IF
                        util_(i,j,k)=(rho0(k)+rho1(i,j,k,pre10,qbar0,the10))*util0(i,j,k)                                        &
                                +dt*(ffrhoutil(i,j,k)+2*ffrhoutil1(i,j,k)+2*ffrhoutil2(i,j,k)+ffrhoutil3(i,j,k))/6
                        IF (util_(i,j,k)>0) THEN
                                util_(i,j,k)=MIN(util_(i,j,k),(rho0(k)+rho1(i,j,k,pre10,qbar0,the10))*max_util)
                        ELSE
                                util_(i,j,k)=MAX(util_(i,j,k),-(rho0(k)+rho1(i,j,k,pre10,qbar0,the10))*max_util)
                        END IF
                        vtil_(i,j,k)=(rho0(k)+rho1(i,j,k,pre10,qbar0,the10))*vtil0(i,j,k)                                        &
                                +dt*(ffrhovtil(i,j,k)+2*ffrhovtil1(i,j,k)+2*ffrhovtil2(i,j,k)+ffrhovtil3(i,j,k))/6
                        IF (vtil_(i,j,k)>0) THEN
                                vtil_(i,j,k)=MIN(vtil_(i,j,k),(rho0(k)+rho1(i,j,k,pre10,qbar0,the10))*max_vtil)
                        ELSE
                                vtil_(i,j,k)=MAX(vtil_(i,j,k),-(rho0(k)+rho1(i,j,k,pre10,qbar0,the10))*max_vtil)
                        END IF
                        wtil_(i,j,k)=(rho0(k)+rho1(i,j,k,pre10,qbar0,the10))*wtil0(i,j,k)                                        &
                                +dt*(ffrhowtil(i,j,k)+2*ffrhowtil1(i,j,k)+2*ffrhowtil2(i,j,k)+ffrhowtil3(i,j,k))/6
                        IF (wtil_(i,j,k)>0) THEN
                                wtil_(i,j,k)=MIN(wtil_(i,j,k),(rho0(k)+rho1(i,j,k,pre10,qbar0,the10))*max_wtil)
                        ELSE
                                wtil_(i,j,k)=MAX(wtil_(i,j,k),-(rho0(k)+rho1(i,j,k,pre10,qbar0,the10))*max_wtil)
                        END IF
		END DO
	END DO
END DO

DO k=1,nz
        DO j=1,ny-1
                qbar_(0,j,k)=qbar_(nx-1,j,k)
                the1_(0,j,k)=the1_(nx-1,j,k)
                ubar_(0,j,k)=ubar_(nx-1,j,k)
                vbar_(0,j,k)=vbar_(nx-1,j,k)
                wbar_(0,j,k)=wbar_(nx-1,j,k)
                util_(0,j,k)=util_(nx-1,j,k)
                vtil_(0,j,k)=vtil_(nx-1,j,k)
                wtil_(0,j,k)=wtil_(nx-1,j,k)

                qbar_(nx,j,k)=qbar_(1,j,k)
                the1_(nx,j,k)=the1_(1,j,k)
                ubar_(nx,j,k)=ubar_(1,j,k)
                vbar_(nx,j,k)=vbar_(1,j,k)
                wbar_(nx,j,k)=wbar_(1,j,k)
                util_(nx,j,k)=util_(1,j,k)
                vtil_(nx,j,k)=vtil_(1,j,k)
                wtil_(nx,j,k)=wtil_(1,j,k)
        END DO

        DO i=0,nx
                qbar_(i,0,k)=qbar_(i,ny-1,k)
                the1_(i,0,k)=the1_(i,ny-1,k)
                ubar_(i,0,k)=ubar_(i,ny-1,k)
                vbar_(i,0,k)=vbar_(i,ny-1,k)
                wbar_(i,0,k)=wbar_(i,ny-1,k)
                util_(i,0,k)=util_(i,ny-1,k)
                vtil_(i,0,k)=vtil_(i,ny-1,k)
                wtil_(i,0,k)=wtil_(i,ny-1,k)

                qbar_(i,ny,k)=qbar_(i,1,k)
                the1_(i,ny,k)=the1_(i,1,k)
                ubar_(i,ny,k)=ubar_(i,1,k)
                vbar_(i,ny,k)=vbar_(i,1,k)
                wbar_(i,ny,k)=wbar_(i,1,k)
                util_(i,ny,k)=util_(i,1,k)
                vtil_(i,ny,k)=vtil_(i,1,k)
                wtil_(i,ny,k)=wtil_(i,1,k)
        END DO
END DO

CALL set_psi2(ubar_,vbar_,wbar_)

CALL set_psi4(time+dt,qbar_,the1_,ubar_,vbar_,wbar_,util_,vtil_,wtil_)

DO k=1,nz
	DO j=1,ny-1
		DO i=1,nx-1
			pre1_(i,j,k)=varpsi(i,j,k)
		END DO
	END DO
END DO

DO k=1,nz
        DO j=1,ny-1
                pre1_(0,j,k)=pre1_(nx-1,j,k)
                pre1_(nx,j,k)=pre1_(1,j,k)
        END DO

        DO i=0,nx
                pre1_(i,0,k)=pre1_(i,ny-1,k)
                pre1_(i,ny,k)=pre1_(i,1,k)
        END DO
END DO

DO k=1,nz
        DO j=1,ny-1
                DO i=1,nx-1
			util_(i,j,k)=util_(i,j,k)/(rho0(k)+rho1(i,j,k,pre1_,qbar_,the1_))
			vtil_(i,j,k)=vtil_(i,j,k)/(rho0(k)+rho1(i,j,k,pre1_,qbar_,the1_))
			wtil_(i,j,k)=wtil_(i,j,k)/(rho0(k)+rho1(i,j,k,pre1_,qbar_,the1_))
		END DO
	END DO
END DO

DO k=1,nz
        DO j=1,ny-1
                util_(0,j,k)=util_(nx-1,j,k)
                vtil_(0,j,k)=vtil_(nx-1,j,k)
                wtil_(0,j,k)=wtil_(nx-1,j,k)

                util_(nx,j,k)=util_(1,j,k)
                vtil_(nx,j,k)=vtil_(1,j,k)
                wtil_(nx,j,k)=wtil_(1,j,k)
        END DO

        DO i=0,nx
                util_(i,0,k)=util_(i,ny-1,k)
                vtil_(i,0,k)=vtil_(i,ny-1,k)
                wtil_(i,0,k)=wtil_(i,ny-1,k)

                util_(i,ny,k)=util_(i,1,k)
                vtil_(i,ny,k)=vtil_(i,1,k)
                wtil_(i,ny,k)=wtil_(i,1,k)
        END DO
END DO

CALL set_phi2(pre1_,qbar_,the1_,util_,vtil_,wtil_)

CALL set_phi(delta0,time+dt,pre1_,qbar_,the1_,ubar_,vbar_,wbar_,util_,vtil_,wtil_)

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

DO k=1,nz
        DO j=0,ny
                DO i=0,nx
                        pre10(i,j,k)=pre1_(i,j,k)
                        qbar0(i,j,k)=qbar_(i,j,k)
                        the10(i,j,k)=the1_(i,j,k)
                        ubar0(i,j,k)=ubar_(i,j,k)
                        vbar0(i,j,k)=vbar_(i,j,k)
                        wbar0(i,j,k)=wbar_(i,j,k)
                        util0(i,j,k)=util_(i,j,k)
                        vtil0(i,j,k)=vtil_(i,j,k)
                        wtil0(i,j,k)=wtil_(i,j,k)
                END DO
        END DO
END DO

END SUBROUTINE rk_projection_method

END MODULE evolve0
