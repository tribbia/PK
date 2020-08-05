MODULE phi_solver

USE parameters, ONLY: nx,ny,nz
USE parameters, ONLY: dx,dy,dz
USE parameters, ONLY: parg

USE auxiliaries, ONLY: ff1,ff2,ff3

USE utilities, ONLY: omega2,omega3
USE utilities, ONLY: parmu3,parxi3
USE utilities, ONLY: ximudzubar,ximudzvbar
USE utilities, ONLY: xidzutil,xidzvtil

USE utilities, ONLY: rho0,rho1
USE utilities, ONLY: the0

IMPLICIT NONE
PUBLIC

REAL,DIMENSION(nz) :: A,B,C
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: phi
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: varphi
REAL,DIMENSION(2*(30+(nx-1)+(ny-1)+2*nz+MAX(MAX(nx-1,ny-1),nz)+7*(NINT(nx/2)+NINT(ny/2)))) :: W

CONTAINS

SUBROUTINE allocate_varphi

IF (.NOT.ALLOCATED(phi)) THEN
        ALLOCATE(phi(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(varphi)) THEN
        ALLOCATE(varphi(nx+1,ny+1,nz))
END IF

END SUBROUTINE allocate_varphi

SUBROUTINE set_phi(delta,time,pre1,qbar,the1,ubar,vbar,wbar,util,vtil,wtil)

REAL,INTENT(IN) :: delta
REAL,INTENT(IN) :: time

REAL,DIMENSION(0:nx,0:ny,nz),INTENT(IN) :: pre1,qbar,the1
REAL,DIMENSION(0:nx,0:ny,nz),INTENT(IN) :: ubar,vbar,wbar
REAL,DIMENSION(0:nx,0:ny,nz),INTENT(IN) :: util,vtil,wtil

! LOCAL VARABLES:

INTEGER :: IERROR
INTEGER :: i,j,k

REAL :: C1,C2

C1=1.0/dx**2
C2=1.0/dy**2

DO k=2,nz-1
        C(k)=1.0/dz**2
        B(k)=-2.0/dz**2
        A(k)=1.0/dz**2
END DO

C(1)=2.0/dz**2
B(1)=-2.0/dz**2
A(1)=0.0

C(nz)=0.0
B(nz)=-2.0/dz**2
A(nz)=2.0/dz**2

! set right hand side of equation: (varphi on output)

!!! aux core: BEGIN
DO k=2,nz-1
        DO j=1,ny-1
                DO i=1,nx-1
ff1(i,j,k)=                                                                                                                      &
2*(rho0(k)+rho1(i,j,k,pre1,qbar,the1))*(vtil(i,j,k)*omega3(j)-wtil(i,j,k)*omega2(j))                                             &
+parxi3(k)*(ubar(i,j,k+1)-2*ubar(i,j,k)+ubar(i,j,k-1))/(dz*dz)                                                       &
+parxi3(k)*(util(i,j,k+1)-2*util(i,j,k)+util(i,j,k-1))/(dz*dz)

ff2(i,j,k)=                                                                                                                      &
-2*(rho0(k)+rho1(i,j,k,pre1,qbar,the1))*util(i,j,k)*omega3(j)                                                                    &
+parxi3(k)*(vbar(i,j,k+1)-2*vbar(i,j,k)+vbar(i,j,k-1))/(dz*dz)                                                       &
+parxi3(k)*(vtil(i,j,k+1)-2*vtil(i,j,k)+vtil(i,j,k-1))/(dz*dz)

ff3(i,j,k)=                                                                                                                      &
2*(rho0(k)+rho1(i,j,k,pre1,qbar,the1))*util(i,j,k)*omega2(j)                                                                     &
+delta*(rho0(k)+rho1(i,j,k,pre1,qbar,the1))*parg*(rho1(i,j,k,pre1,qbar,the1)/rho0(k))**2                                         &
+parxi3(k)*(wbar(i,j,k+1)-2*wbar(i,j,k)+wbar(i,j,k-1))/(dz*dz)                                                       &
+parxi3(k)*(wtil(i,j,k+1)-2*wtil(i,j,k)+wtil(i,j,k-1))/(dz*dz)
                END DO
        END DO
END DO

k=1
DO j=1,ny-1
        DO i=1,nx-1
ff1(i,j,k)=                                                                                                                      &
2*(rho0(k)+rho1(i,j,k,pre1,qbar,the1))*vtil(i,j,k)*omega3(j)                                                                     &
+((parxi3(k+1)-parmu3(k+1))*(ubar(i,j,k+2)-ubar(i,j,k))/(2*dz)-ximudzubar(i,j,ubar,vbar))/dz                                     &
+(parxi3(k+1)*(util(i,j,k+2)-util(i,j,k))/(2*dz)-xidzutil(i,j,ubar,vbar,util))/dz

ff2(i,j,k)=                                                                                                                      &
-2*(rho0(k)+rho1(i,j,k,pre1,qbar,the1))*util(i,j,k)*omega3(j)                                                                    &
+((parxi3(k+1)-parmu3(k+1))*(vbar(i,j,k+2)-vbar(i,j,k))/(2*dz)-ximudzvbar(i,j,ubar,vbar))/dz                                     &
+(parxi3(k+1)*(vtil(i,j,k+2)-vtil(i,j,k))/(2*dz)-xidzvtil(i,j,ubar,vbar,vtil))/dz

ff3(i,j,k)=                                                                                                                      &
2*(rho0(k)+rho1(i,j,k,pre1,qbar,the1))*util(i,j,k)*omega2(j)                                                                     &
+delta*(rho0(k)+rho1(i,j,k,pre1,qbar,the1))*parg*(rho1(i,j,k,pre1,qbar,the1)/rho0(k))**2                                         &
+parxi3(k)*(wbar(i,j,k+2)-2*wbar(i,j,k+1))/(dz*dz)                                                                   &
+parxi3(k)*(wtil(i,j,k+2)-2*wtil(i,j,k+1))/(dz*dz)
        END DO
END DO

k=nz
DO j=1,ny-1
        DO i=1,nx-1
ff1(i,j,k)=                                                                                                                      &
2*(rho0(k)+rho1(i,j,k,pre1,qbar,the1))*vtil(i,j,k)*omega3(j)                                                                     &
+parxi3(k)*(ubar(i,j,k)-2*ubar(i,j,k-1)+ubar(i,j,k-2))/(dz*dz)                                                       &
+parxi3(k)*(util(i,j,k)-2*util(i,j,k-1)+util(i,j,k-2))/(dz*dz)

ff2(i,j,k)=                                                                                                                      &
-2*(rho0(k)+rho1(i,j,k,pre1,qbar,the1))*util(i,j,k)*omega3(j)                                                                    &
+parxi3(k)*(vbar(i,j,k)-2*vbar(i,j,k-1)+vbar(i,j,k-2))/(dz*dz)                                                       &
+parxi3(k)*(vtil(i,j,k)-2*vtil(i,j,k-1)+vtil(i,j,k-2))/(dz*dz)

ff3(i,j,k)=                                                                                                                      &
2*(rho0(k)+rho1(i,j,k,pre1,qbar,the1))*util(i,j,k)*omega2(j)                                                                     &
+delta*(rho0(k)+rho1(i,j,k,pre1,qbar,the1))*parg*(rho1(i,j,k,pre1,qbar,the1)/rho0(k))**2                                         &
+parxi3(k)*(wbar(i,j,k-2)-2*wbar(i,j,k-1))/(dz*dz)                                                                   &
+parxi3(k)*(wtil(i,j,k-2)-2*wtil(i,j,k-1))/(dz*dz)
        END DO
END DO
!!! aux core: END

!!! aux ghosts: BEGIN
DO k=1,nz
        DO j=1,ny-1
                ff1(0,j,k)=ff1(nx-1,j,k)
                ff2(0,j,k)=ff2(nx-1,j,k)
                ff3(0,j,k)=ff3(nx-1,j,k)

                ff1(nx,j,k)=ff1(1,j,k)
                ff2(nx,j,k)=ff2(1,j,k)
                ff3(nx,j,k)=ff3(1,j,k)
        END DO

        DO i=0,nx
                ff1(i,0,k)=ff1(i,ny-1,k)
                ff2(i,0,k)=ff2(i,ny-1,k)
                ff3(i,0,k)=ff3(i,ny-1,k)

                ff1(i,ny,k)=ff1(i,1,k)
                ff2(i,ny,k)=ff2(i,1,k)
                ff3(i,ny,k)=ff3(i,1,k)
        END DO
END DO
!!! aux ghosts: END

!!! rhs core: BEGIN
DO k=2,nz-1
        DO j=1,ny-1
                DO i=1,nx-1
varphi(i,j,k)=(ff1(i+1,j,k)-ff1(i-1,j,k))/(2*dx)+(ff2(i,j+1,k)-ff2(i,j-1,k))/(2*dy)+(ff3(i,j,k+1)-ff3(i,j,k-1))/(2*dz)
                END DO
        END DO
END DO

k=1
DO j=1,ny-1
        DO i=1,nx-1
varphi(i,j,k)=(ff1(i+1,j,k)-ff1(i-1,j,k))/(2*dx)+(ff2(i,j+1,k)-ff2(i,j-1,k))/(2*dy)+(ff3(i,j,k+1)-ff3(i,j,k))/dz                 &
+(2*(rho0(k)+rho1(i,j,k,pre1,qbar,the1))*util(i,j,k)*omega2(j)                                                                   &
+delta*(rho0(k)+rho1(i,j,k,pre1,qbar,the1))*parg*(rho1(i,j,k,pre1,qbar,the1)/rho0(k))**2                                         &
+parxi3(k)*(wbar(i,j,k+2)-2*wbar(i,j,k+1))/(dz*dz)+parxi3(k)*(wtil(i,j,k+2)-2*wtil(i,j,k+1))/(dz*dz))*2.0/dz
        END DO
END DO

k=nz
DO j=1,ny-1
        DO i=1,nx-1
varphi(i,j,k)=(ff1(i+1,j,k)-ff1(i-1,j,k))/(2*dx)+(ff2(i,j+1,k)-ff2(i,j-1,k))/(2*dy)+(ff3(i,j,k)-ff3(i,j,k-1))/dz                 &
-(2*(rho0(k)+rho1(i,j,k,pre1,qbar,the1))*util(i,j,k)*omega2(j)                                                                   &
+delta*(rho0(k)+rho1(i,j,k,pre1,qbar,the1))*parg*(rho1(i,j,k,pre1,qbar,the1)/rho0(k))**2                                         &
+parxi3(k)*(wbar(i,j,k-2)-2*wbar(i,j,k-1))/(dz*dz)+parxi3(k)*(wtil(i,j,k-2)-2*wtil(i,j,k-1))/(dz*dz))*2.0/dz
        END DO
END DO
!!! rhs core: END

! set phi:

CALL pois3d(0,nx-1,C1,0,ny-1,C2,1,nz,A,B,C,nx+1,ny+1,varphi,IERROR,W)

END SUBROUTINE set_phi

SUBROUTINE set_phi2(pre1,qbar,the1,util,vtil,wtil)

REAL,DIMENSION(0:nx,0:ny,nz),INTENT(IN) :: pre1,qbar,the1
REAL,DIMENSION(0:nx,0:ny,nz),INTENT(INOUT) :: util,vtil,wtil

! LOCAL VARABLES:

INTEGER :: IERROR
INTEGER :: i,j,k

REAL :: C1,C2

C1=1.0/dx**2
C2=1.0/dy**2

DO k=2,nz-1
        C(k)=1.0/dz**2
        B(k)=-2.0/dz**2
        A(k)=1.0/dz**2
END DO

C(1)=2.0/dz**2
B(1)=-2.0/dz**2
A(1)=0.0

C(nz)=0.0
B(nz)=-2.0/dz**2
A(nz)=2.0/dz**2

! set right hand side of equation: (varphi on output)

!!! rhs core: BEGIN
DO k=2,nz-1
        DO j=1,ny-1
                DO i=1,nx-1
varphi(i,j,k)=((rho0(k)+rho1(i+1,j,k,pre1,qbar,the1))*util(i+1,j,k)-(rho0(k)+rho1(i-1,j,k,pre1,qbar,the1))*util(i-1,j,k))/(2*dx) &
             +((rho0(k)+rho1(i,j+1,k,pre1,qbar,the1))*vtil(i,j+1,k)-(rho0(k)+rho1(i,j-1,k,pre1,qbar,the1))*vtil(i,j-1,k))/(2*dy) &
         +((rho0(k+1)+rho1(i,j,k+1,pre1,qbar,the1))*wtil(i,j,k+1)-(rho0(k-1)+rho1(i,j,k-1,pre1,qbar,the1))*wtil(i,j,k-1))/(2*dz)
                END DO
        END DO
END DO

k=1
DO j=1,ny-1
        DO i=1,nx-1
varphi(i,j,k)=((rho0(k)+rho1(i+1,j,k,pre1,qbar,the1))*util(i+1,j,k)-(rho0(k)+rho1(i-1,j,k,pre1,qbar,the1))*util(i-1,j,k))/(2*dx) &
             +((rho0(k)+rho1(i,j+1,k,pre1,qbar,the1))*vtil(i,j+1,k)-(rho0(k)+rho1(i,j-1,k,pre1,qbar,the1))*vtil(i,j-1,k))/(2*dy) &
             +((rho0(k+1)+rho1(i,j,k+1,pre1,qbar,the1))*wtil(i,j,k+1)-(rho0(k)+rho1(i,j,k,pre1,qbar,the1))*wtil(i,j,k))/dz
        END DO
END DO

k=nz
DO j=1,ny-1
        DO i=1,nx-1
varphi(i,j,k)=((rho0(k)+rho1(i+1,j,k,pre1,qbar,the1))*util(i+1,j,k)-(rho0(k)+rho1(i-1,j,k,pre1,qbar,the1))*util(i-1,j,k))/(2*dx) &
             +((rho0(k)+rho1(i,j+1,k,pre1,qbar,the1))*vtil(i,j+1,k)-(rho0(k)+rho1(i,j-1,k,pre1,qbar,the1))*vtil(i,j-1,k))/(2*dy) &
             +((rho0(k)+rho1(i,j,k,pre1,qbar,the1))*wtil(i,j,k)-(rho0(k-1)+rho1(i,j,k-1,pre1,qbar,the1))*wtil(i,j,k-1))/dz
        END DO
END DO
!!! rhs core: END

! set phi:

CALL pois3d(0,nx-1,C1,0,ny-1,C2,1,nz,A,B,C,nx+1,ny+1,varphi,IERROR,W)

DO k=1,nz
	DO j=1,ny-1
		DO i=1,nx-1
			phi(i,j,k)=varphi(i,j,k)
		END DO
	END DO
END DO

DO k=1,nz
        DO j=1,ny-1
                phi(0,j,k)=phi(nx-1,j,k)
                phi(nx,j,k)=phi(1,j,k)
        END DO

        DO i=0,nx
                phi(i,0,k)=phi(i,ny-1,k)
                phi(i,ny,k)=phi(i,1,k)
        END DO
END DO

DO k=2,nz-1
	DO j=1,ny-1
		DO i=1,nx-1
util(i,j,k)=util(i,j,k)-((phi(i+1,j,k)-phi(i-1,j,k))/(2*dx))/(rho0(k)+rho1(i,j,k,pre1,qbar,the1))
vtil(i,j,k)=vtil(i,j,k)-((phi(i,j+1,k)-phi(i,j-1,k))/(2*dy))/(rho0(k)+rho1(i,j,k,pre1,qbar,the1))
wtil(i,j,k)=wtil(i,j,k)-((phi(i,j,k+1)-phi(i,j,k-1))/(2*dz))/(rho0(k)+rho1(i,j,k,pre1,qbar,the1))
		END DO
	END DO
END DO

k=1
DO j=1,ny-1
        DO i=1,nx-1
util(i,j,k)=util(i,j,k)-((phi(i+1,j,k)-phi(i-1,j,k))/(2*dx))/(rho0(k)+rho1(i,j,k,pre1,qbar,the1))
vtil(i,j,k)=vtil(i,j,k)-((phi(i,j+1,k)-phi(i,j-1,k))/(2*dy))/(rho0(k)+rho1(i,j,k,pre1,qbar,the1))
        END DO
END DO

k=nz
DO j=1,ny-1
        DO i=1,nx-1
util(i,j,k)=util(i,j,k)-((phi(i+1,j,k)-phi(i-1,j,k))/(2*dx))/(rho0(k)+rho1(i,j,k,pre1,qbar,the1))
vtil(i,j,k)=vtil(i,j,k)-((phi(i,j+1,k)-phi(i,j-1,k))/(2*dy))/(rho0(k)+rho1(i,j,k,pre1,qbar,the1))
        END DO
END DO

DO k=1,nz
        DO j=1,ny-1
                util(0,j,k)=util(nx-1,j,k)
                vtil(0,j,k)=vtil(nx-1,j,k)
                wtil(0,j,k)=wtil(nx-1,j,k)

                util(nx,j,k)=util(1,j,k)
                vtil(nx,j,k)=vtil(1,j,k)
                wtil(nx,j,k)=wtil(1,j,k)
        END DO

        DO i=0,nx
                util(i,0,k)=util(i,ny-1,k)
                vtil(i,0,k)=vtil(i,ny-1,k)
                wtil(i,0,k)=wtil(i,ny-1,k)

                util(i,ny,k)=util(i,1,k)
                vtil(i,ny,k)=vtil(i,1,k)
                wtil(i,ny,k)=wtil(i,1,k)
        END DO
END DO

END SUBROUTINE set_phi2

END MODULE phi_solver
