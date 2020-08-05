MODULE psi_solver

USE parameters, ONLY: nx,ny,nz
USE parameters, ONLY: dx,dy,dz
USE parameters, ONLY: parg

USE auxiliaries, ONLY: ff1,ff2,ff3

USE utilities, ONLY: omega2,omega3
USE utilities, ONLY: parmu3
USE utilities, ONLY: mudzubar,mudzvbar

USE utilities, ONLY: rho0
USE utilities, ONLY: the0

IMPLICIT NONE
PUBLIC

REAL,DIMENSION(nz) :: A,B,C
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: psi
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: varpsi
REAL,DIMENSION(2*(30+(nx-1)+(ny-1)+2*nz+MAX(MAX(nx-1,ny-1),nz)+7*(NINT(nx/2)+NINT(ny/2)))) :: W

CONTAINS

SUBROUTINE allocate_varpsi

IF (.NOT.ALLOCATED(psi)) THEN
        ALLOCATE(psi(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(varpsi)) THEN
        ALLOCATE(varpsi(nx+1,ny+1,nz))
END IF

END SUBROUTINE allocate_varpsi

SUBROUTINE set_psi(time,pre1,qbar,the1,ubar,vbar,wbar,util,vtil,wtil)

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

! set right hand side of equation: (varpsi on output)

!!! aux core: BEGIN
DO k=2,nz-1
        DO j=1,ny-1
                DO i=1,nx-1
ff1(i,j,k)=                                                                                                                      &
-rho0(k)*(ubar(i,j,k)+util(i,j,k))*(ubar(i+1,j,k)-ubar(i-1,j,k))/(2*dx)                                                          &
-rho0(k)*(vbar(i,j,k)+vtil(i,j,k))*(ubar(i,j+1,k)-ubar(i,j-1,k))/(2*dy)                                                          &
-rho0(k)*(wbar(i,j,k)+wtil(i,j,k))*(ubar(i,j,k+1)-ubar(i,j,k-1))/(2*dz)                                                          &
+2*rho0(k)*(vbar(i,j,k)*omega3(j)-wbar(i,j,k)*omega2(j))                                                                         &
+parmu3(k)*(ubar(i,j,k+1)-2*ubar(i,j,k)+ubar(i,j,k-1))/(dz*dz)

ff2(i,j,k)=                                                                                                                      &
-rho0(k)*(ubar(i,j,k)+util(i,j,k))*(vbar(i+1,j,k)-vbar(i-1,j,k))/(2*dx)                                                          &
-rho0(k)*(vbar(i,j,k)+vtil(i,j,k))*(vbar(i,j+1,k)-vbar(i,j-1,k))/(2*dy)                                                          &
-rho0(k)*(wbar(i,j,k)+wtil(i,j,k))*(vbar(i,j,k+1)-vbar(i,j,k-1))/(2*dz)                                                          &
-2*rho0(k)*ubar(i,j,k)*omega3(j)                                                                                                 &
+parmu3(k)*(vbar(i,j,k+1)-2*vbar(i,j,k)+vbar(i,j,k-1))/(dz*dz)

ff3(i,j,k)=                                                                                                                      &
-rho0(k)*(ubar(i,j,k)+util(i,j,k))*(wbar(i+1,j,k)-wbar(i-1,j,k))/(2*dx)                                                          &
-rho0(k)*(vbar(i,j,k)+vtil(i,j,k))*(wbar(i,j+1,k)-wbar(i,j-1,k))/(2*dy)                                                          &
-rho0(k)*(wbar(i,j,k)+wtil(i,j,k))*(wbar(i,j,k+1)-wbar(i,j,k-1))/(2*dz)                                                          &
+2*rho0(k)*ubar(i,j,k)*omega2(j)                                                                                                 &
+rho0(k)*parg*the1(i,j,k)/the0(k)                                                                                                &
+parmu3(k)*(wbar(i,j,k+1)-2*wbar(i,j,k)+wbar(i,j,k-1))/(dz*dz)
                END DO
        END DO
END DO

k=1
DO j=1,ny-1
        DO i=1,nx-1
ff1(i,j,k)=                                                                                                                      &
-rho0(k)*(ubar(i,j,k)+util(i,j,k))*(ubar(i+1,j,k)-ubar(i-1,j,k))/(2*dx)                                                          &
-rho0(k)*(vbar(i,j,k)+vtil(i,j,k))*(ubar(i,j+1,k)-ubar(i,j-1,k))/(2*dy)                                                          &
+2*rho0(k)*vbar(i,j,k)*omega3(j)                                                                                                 &
+(parmu3(k+1)*(ubar(i,j,k+2)-ubar(i,j,k))/(2*dz)-mudzubar(i,j,ubar,vbar))/dz

ff2(i,j,k)=                                                                                                                      &
-rho0(k)*(ubar(i,j,k)+util(i,j,k))*(vbar(i+1,j,k)-vbar(i-1,j,k))/(2*dx)                                                          &
-rho0(k)*(vbar(i,j,k)+vtil(i,j,k))*(vbar(i,j+1,k)-vbar(i,j-1,k))/(2*dy)                                                          &
-2*rho0(k)*ubar(i,j,k)*omega3(j)                                                                                                 &
+(parmu3(k+1)*(vbar(i,j,k+2)-vbar(i,j,k))/(2*dz)-mudzvbar(i,j,ubar,vbar))/dz

ff3(i,j,k)=                                                                                                                      &
2*rho0(k)*ubar(i,j,k)*omega2(j)                                                                                                  &
+rho0(k)*parg*the1(i,j,k)/the0(k)                                                                                                &
+parmu3(k)*(wbar(i,j,k+2)-2*wbar(i,j,k+1))/(dz*dz)
        END DO
END DO

k=nz
DO j=1,ny-1
        DO i=1,nx-1
ff1(i,j,k)=                                                                                                                      &
-rho0(k)*(ubar(i,j,k)+util(i,j,k))*(ubar(i+1,j,k)-ubar(i-1,j,k))/(2*dx)                                                          &
-rho0(k)*(vbar(i,j,k)+vtil(i,j,k))*(ubar(i,j+1,k)-ubar(i,j-1,k))/(2*dy)                                                          &
+2*rho0(k)*vbar(i,j,k)*omega3(j)                                                                                                 &
+parmu3(k)*(ubar(i,j,k)-2*ubar(i,j,k-1)+ubar(i,j,k-2))/(dz*dz)

ff2(i,j,k)=                                                                                                                      &
-rho0(k)*(ubar(i,j,k)+util(i,j,k))*(vbar(i+1,j,k)-vbar(i-1,j,k))/(2*dx)                                                          &
-rho0(k)*(vbar(i,j,k)+vtil(i,j,k))*(vbar(i,j+1,k)-vbar(i,j-1,k))/(2*dy)                                                          &
-2*rho0(k)*ubar(i,j,k)*omega3(j)                                                                                                 &
+parmu3(k)*(vbar(i,j,k)-2*vbar(i,j,k-1)+vbar(i,j,k-2))/(dz*dz)

ff3(i,j,k)=                                                                                                                      &
2*rho0(k)*ubar(i,j,k)*omega2(j)                                                                                                  &
+parmu3(k)*(wbar(i,j,k-2)-2*wbar(i,j,k-1))/(dz*dz)
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
varpsi(i,j,k)=(ff1(i+1,j,k)-ff1(i-1,j,k))/(2*dx)+(ff2(i,j+1,k)-ff2(i,j-1,k))/(2*dy)+(ff3(i,j,k+1)-ff3(i,j,k-1))/(2*dz)
                END DO
        END DO
END DO

k=1
DO j=1,ny-1
        DO i=1,nx-1
varpsi(i,j,k)=(ff1(i+1,j,k)-ff1(i-1,j,k))/(2*dx)+(ff2(i,j+1,k)-ff2(i,j-1,k))/(2*dy)+(ff3(i,j,k+1)-ff3(i,j,k))/dz                 &
+(2*rho0(k)*ubar(i,j,k)*omega2(j)+rho0(k)*parg*the1(i,j,k)/the0(k)+parmu3(k)*(wbar(i,j,k+2)-2*wbar(i,j,k+1))/(dz*dz))*2.0/dz
        END DO
END DO

k=nz
DO j=1,ny-1
        DO i=1,nx-1
varpsi(i,j,k)=(ff1(i+1,j,k)-ff1(i-1,j,k))/(2*dx)+(ff2(i,j+1,k)-ff2(i,j-1,k))/(2*dy)+(ff3(i,j,k)-ff3(i,j,k-1))/dz                 &
-(2*rho0(k)*ubar(i,j,k)*omega2(j)+parmu3(k)*(wbar(i,j,k-2)-2*wbar(i,j,k-1))/(dz*dz))*2.0/dz
        END DO
END DO
!!! rhs core: END

! set psi:

CALL pois3d(0,nx-1,C1,0,ny-1,C2,1,nz,A,B,C,nx+1,ny+1,varpsi,IERROR,W)

END SUBROUTINE set_psi

SUBROUTINE set_psi2(ubar,vbar,wbar)

REAL,DIMENSION(0:nx,0:ny,nz),INTENT(INOUT) :: ubar,vbar,wbar

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

! set right hand side of equation: (varpsi on output)

!!! rhs core: BEGIN
DO k=2,nz-1
        DO j=1,ny-1
                DO i=1,nx-1
varpsi(i,j,k)=(rho0(k)*ubar(i+1,j,k)-rho0(k)*ubar(i-1,j,k))/(2*dx)+(rho0(k)*vbar(i,j+1,k)-rho0(k)*vbar(i,j-1,k))/(2*dy)          &
+(rho0(k+1)*wbar(i,j,k+1)-rho0(k-1)*wbar(i,j,k-1))/(2*dz)
                END DO
        END DO
END DO

k=1
DO j=1,ny-1
        DO i=1,nx-1
varpsi(i,j,k)=(rho0(k)*ubar(i+1,j,k)-rho0(k)*ubar(i-1,j,k))/(2*dx)+(rho0(k)*vbar(i,j+1,k)-rho0(k)*vbar(i,j-1,k))/(2*dy)          &
+(rho0(k+1)*wbar(i,j,k+1)-rho0(k)*wbar(i,j,k))/dz
        END DO
END DO

k=nz
DO j=1,ny-1
        DO i=1,nx-1
varpsi(i,j,k)=(rho0(k)*ubar(i+1,j,k)-rho0(k)*ubar(i-1,j,k))/(2*dx)+(rho0(k)*vbar(i,j+1,k)-rho0(k)*vbar(i,j-1,k))/(2*dy)          &
+(rho0(k)*wbar(i,j,k)-rho0(k-1)*wbar(i,j,k-1))/dz
        END DO
END DO
!!! rhs core: END

! set psi:

CALL pois3d(0,nx-1,C1,0,ny-1,C2,1,nz,A,B,C,nx+1,ny+1,varpsi,IERROR,W)

DO k=1,nz
	DO j=1,ny-1
		DO i=1,nx-1
			psi(i,j,k)=varpsi(i,j,k)
		END DO
	END DO
END DO

DO k=1,nz
        DO j=1,ny-1
                psi(0,j,k)=psi(nx-1,j,k)
                psi(nx,j,k)=psi(1,j,k)
        END DO

        DO i=0,nx
                psi(i,0,k)=psi(i,ny-1,k)
                psi(i,ny,k)=psi(i,1,k)
        END DO
END DO

DO k=2,nz-1
	DO j=1,ny-1
		DO i=1,nx-1
ubar(i,j,k)=ubar(i,j,k)-((psi(i+1,j,k)-psi(i-1,j,k))/(2*dx))/rho0(k)
vbar(i,j,k)=vbar(i,j,k)-((psi(i,j+1,k)-psi(i,j-1,k))/(2*dy))/rho0(k)
wbar(i,j,k)=wbar(i,j,k)-((psi(i,j,k+1)-psi(i,j,k-1))/(2*dz))/rho0(k)
		END DO
	END DO
END DO

k=1
DO j=1,ny-1
        DO i=1,nx-1
ubar(i,j,k)=ubar(i,j,k)-((psi(i+1,j,k)-psi(i-1,j,k))/(2*dx))/rho0(k)
vbar(i,j,k)=vbar(i,j,k)-((psi(i,j+1,k)-psi(i,j-1,k))/(2*dy))/rho0(k)
        END DO
END DO

k=nz
DO j=1,ny-1
        DO i=1,nx-1
ubar(i,j,k)=ubar(i,j,k)-((psi(i+1,j,k)-psi(i-1,j,k))/(2*dx))/rho0(k)
vbar(i,j,k)=vbar(i,j,k)-((psi(i,j+1,k)-psi(i,j-1,k))/(2*dy))/rho0(k)
        END DO
END DO

DO k=1,nz
        DO j=1,ny-1
                ubar(0,j,k)=ubar(nx-1,j,k)
                vbar(0,j,k)=vbar(nx-1,j,k)
                wbar(0,j,k)=wbar(nx-1,j,k)

                ubar(nx,j,k)=ubar(1,j,k)
                vbar(nx,j,k)=vbar(1,j,k)
                wbar(nx,j,k)=wbar(1,j,k)
        END DO

        DO i=0,nx
                ubar(i,0,k)=ubar(i,ny-1,k)
                vbar(i,0,k)=vbar(i,ny-1,k)
                wbar(i,0,k)=wbar(i,ny-1,k)

                ubar(i,ny,k)=ubar(i,1,k)
                vbar(i,ny,k)=vbar(i,1,k)
                wbar(i,ny,k)=wbar(i,1,k)
        END DO
END DO

END SUBROUTINE set_psi2

SUBROUTINE set_psi3(ubar,vbar,wbar)

REAL,DIMENSION(0:nx,0:ny,nz),INTENT(INOUT) :: ubar,vbar,wbar

! LOCAL VARABLES:

INTEGER :: IERROR
INTEGER :: i,j,k

REAL :: C1,C2

C1=0
C2=0

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

! set right hand side of equation: (varpsi on output)

!!! rhs core: BEGIN
DO k=2,nz-1
        DO j=1,ny-1
                DO i=1,nx-1
varpsi(i,j,k)=(ubar(i+1,j,k)-ubar(i-1,j,k))/(2*dx)+(vbar(i,j+1,k)-vbar(i,j-1,k))/(2*dy)+(wbar(i,j,k+1)-wbar(i,j,k-1))/(2*dz)
                END DO
        END DO
END DO

k=1
DO j=1,ny-1
        DO i=1,nx-1
varpsi(i,j,k)=(ubar(i+1,j,k)-ubar(i-1,j,k))/(2*dx)+(vbar(i,j+1,k)-vbar(i,j-1,k))/(2*dy)+(wbar(i,j,k+1)-wbar(i,j,k))/dz
        END DO
END DO

k=nz
DO j=1,ny-1
        DO i=1,nx-1
varpsi(i,j,k)=(ubar(i+1,j,k)-ubar(i-1,j,k))/(2*dx)+(vbar(i,j+1,k)-vbar(i,j-1,k))/(2*dy)+(wbar(i,j,k)-wbar(i,j,k-1))/dz
        END DO
END DO
!!! rhs core: END

! set psi:

CALL pois3d(0,nx-1,C1,0,ny-1,C2,1,nz,A,B,C,nx+1,ny+1,varpsi,IERROR,W)

DO k=1,nz
	DO j=1,ny-1
		DO i=1,nx-1
			psi(i,j,k)=varpsi(i,j,k)
		END DO
	END DO
END DO

DO k=1,nz
        DO j=1,ny-1
                psi(0,j,k)=psi(nx-1,j,k)
                psi(nx,j,k)=psi(1,j,k)
        END DO

        DO i=0,nx
                psi(i,0,k)=psi(i,ny-1,k)
                psi(i,ny,k)=psi(i,1,k)
        END DO
END DO

DO k=2,nz-1
	DO j=1,ny-1
		DO i=1,nx-1
wbar(i,j,k)=wbar(i,j,k)-(psi(i,j,k+1)-psi(i,j,k-1))/(2*dz)
		END DO
	END DO
END DO

DO k=1,nz
        DO j=1,ny-1
                wbar(0,j,k)=wbar(nx-1,j,k)

                wbar(nx,j,k)=wbar(1,j,k)
        END DO

        DO i=0,nx
                wbar(i,0,k)=wbar(i,ny-1,k)

                wbar(i,ny,k)=wbar(i,1,k)
        END DO
END DO

END SUBROUTINE set_psi3

SUBROUTINE set_psi4(time,qbar,the1,ubar,vbar,wbar,rhoutil,rhovtil,rhowtil)

REAL,INTENT(IN) :: time

REAL,DIMENSION(0:nx,0:ny,nz),INTENT(IN) :: qbar,the1
REAL,DIMENSION(0:nx,0:ny,nz),INTENT(IN) :: ubar,vbar,wbar
REAL,DIMENSION(0:nx,0:ny,nz),INTENT(IN) :: rhoutil,rhovtil,rhowtil

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

! set right hand side of equation: (varpsi on output)

!!! aux core: BEGIN
DO k=2,nz-1
        DO j=1,ny-1
                DO i=1,nx-1
ff1(i,j,k)=                                                                                                                      &
-(rho0(k)*ubar(i,j,k)+rhoutil(i,j,k))*(ubar(i+1,j,k)-ubar(i-1,j,k))/(2*dx)                                                       &
-(rho0(k)*vbar(i,j,k)+rhovtil(i,j,k))*(ubar(i,j+1,k)-ubar(i,j-1,k))/(2*dy)                                                       &
-(rho0(k)*wbar(i,j,k)+rhowtil(i,j,k))*(ubar(i,j,k+1)-ubar(i,j,k-1))/(2*dz)                                                       &
+2*rho0(k)*(vbar(i,j,k)*omega3(j)-wbar(i,j,k)*omega2(j))                                                                         &
+parmu3(k)*(ubar(i,j,k+1)-2*ubar(i,j,k)+ubar(i,j,k-1))/(dz*dz)

ff2(i,j,k)=                                                                                                                      &
-(rho0(k)*ubar(i,j,k)+rhoutil(i,j,k))*(vbar(i+1,j,k)-vbar(i-1,j,k))/(2*dx)                                                       &
-(rho0(k)*vbar(i,j,k)+rhovtil(i,j,k))*(vbar(i,j+1,k)-vbar(i,j-1,k))/(2*dy)                                                       &
-(rho0(k)*wbar(i,j,k)+rhowtil(i,j,k))*(vbar(i,j,k+1)-vbar(i,j,k-1))/(2*dz)                                                       &
-2*rho0(k)*ubar(i,j,k)*omega3(j)                                                                                                 &
+parmu3(k)*(vbar(i,j,k+1)-2*vbar(i,j,k)+vbar(i,j,k-1))/(dz*dz)

ff3(i,j,k)=                                                                                                                      &
-(rho0(k)*ubar(i,j,k)+rhoutil(i,j,k))*(wbar(i+1,j,k)-wbar(i-1,j,k))/(2*dx)                                                       &
-(rho0(k)*vbar(i,j,k)+rhovtil(i,j,k))*(wbar(i,j+1,k)-wbar(i,j-1,k))/(2*dy)                                                       &
-(rho0(k)*wbar(i,j,k)+rhowtil(i,j,k))*(wbar(i,j,k+1)-wbar(i,j,k-1))/(2*dz)                                                       &
+2*rho0(k)*ubar(i,j,k)*omega2(j)                                                                                                 &
+rho0(k)*parg*the1(i,j,k)/the0(k)                                                                                                &
+parmu3(k)*(wbar(i,j,k+1)-2*wbar(i,j,k)+wbar(i,j,k-1))/(dz*dz)
                END DO
        END DO
END DO

k=1
DO j=1,ny-1
        DO i=1,nx-1
ff1(i,j,k)=                                                                                                                      &
-(rho0(k)*ubar(i,j,k)+rhoutil(i,j,k))*(ubar(i+1,j,k)-ubar(i-1,j,k))/(2*dx)                                                       &
-(rho0(k)*vbar(i,j,k)+rhovtil(i,j,k))*(ubar(i,j+1,k)-ubar(i,j-1,k))/(2*dy)                                                       &
+2*rho0(k)*vbar(i,j,k)*omega3(j)                                                                                                 &
+(parmu3(k+1)*(ubar(i,j,k+2)-ubar(i,j,k))/(2*dz)-mudzubar(i,j,ubar,vbar))/dz

ff2(i,j,k)=                                                                                                                      &
-(rho0(k)*ubar(i,j,k)+rhoutil(i,j,k))*(vbar(i+1,j,k)-vbar(i-1,j,k))/(2*dx)                                                       &
-(rho0(k)*vbar(i,j,k)+rhovtil(i,j,k))*(vbar(i,j+1,k)-vbar(i,j-1,k))/(2*dy)                                                       &
-2*rho0(k)*ubar(i,j,k)*omega3(j)                                                                                                 &
+(parmu3(k+1)*(vbar(i,j,k+2)-vbar(i,j,k))/(2*dz)-mudzvbar(i,j,ubar,vbar))/dz

ff3(i,j,k)=                                                                                                                      &
2*rho0(k)*ubar(i,j,k)*omega2(j)                                                                                                  &
+rho0(k)*parg*the1(i,j,k)/the0(k)                                                                                                &
+parmu3(k)*(wbar(i,j,k+2)-2*wbar(i,j,k+1))/(dz*dz)
        END DO
END DO

k=nz
DO j=1,ny-1
        DO i=1,nx-1
ff1(i,j,k)=                                                                                                                      &
-(rho0(k)*ubar(i,j,k)+rhoutil(i,j,k))*(ubar(i+1,j,k)-ubar(i-1,j,k))/(2*dx)                                                       &
-(rho0(k)*vbar(i,j,k)+rhovtil(i,j,k))*(ubar(i,j+1,k)-ubar(i,j-1,k))/(2*dy)                                                       &
+2*rho0(k)*vbar(i,j,k)*omega3(j)                                                                                                 &
+parmu3(k)*(ubar(i,j,k)-2*ubar(i,j,k-1)+ubar(i,j,k-2))/(dz*dz)

ff2(i,j,k)=                                                                                                                      &
-(rho0(k)*ubar(i,j,k)+rhoutil(i,j,k))*(vbar(i+1,j,k)-vbar(i-1,j,k))/(2*dx)                                                       &
-(rho0(k)*vbar(i,j,k)+rhovtil(i,j,k))*(vbar(i,j+1,k)-vbar(i,j-1,k))/(2*dy)                                                       &
-2*rho0(k)*ubar(i,j,k)*omega3(j)                                                                                                 &
+parmu3(k)*(vbar(i,j,k)-2*vbar(i,j,k-1)+vbar(i,j,k-2))/(dz*dz)

ff3(i,j,k)=                                                                                                                      &
2*rho0(k)*ubar(i,j,k)*omega2(j)                                                                                                  &
+parmu3(k)*(wbar(i,j,k-2)-2*wbar(i,j,k-1))/(dz*dz)
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
varpsi(i,j,k)=(ff1(i+1,j,k)-ff1(i-1,j,k))/(2*dx)+(ff2(i,j+1,k)-ff2(i,j-1,k))/(2*dy)+(ff3(i,j,k+1)-ff3(i,j,k-1))/(2*dz)
                END DO
        END DO
END DO

k=1
DO j=1,ny-1
        DO i=1,nx-1
varpsi(i,j,k)=(ff1(i+1,j,k)-ff1(i-1,j,k))/(2*dx)+(ff2(i,j+1,k)-ff2(i,j-1,k))/(2*dy)+(ff3(i,j,k+1)-ff3(i,j,k))/dz                 &
+(2*rho0(k)*ubar(i,j,k)*omega2(j)+rho0(k)*parg*the1(i,j,k)/the0(k)+parmu3(k)*(wbar(i,j,k+2)-2*wbar(i,j,k+1))/(dz*dz))*2.0/dz
        END DO
END DO

k=nz
DO j=1,ny-1
        DO i=1,nx-1
varpsi(i,j,k)=(ff1(i+1,j,k)-ff1(i-1,j,k))/(2*dx)+(ff2(i,j+1,k)-ff2(i,j-1,k))/(2*dy)+(ff3(i,j,k)-ff3(i,j,k-1))/dz                 &
-(2*rho0(k)*ubar(i,j,k)*omega2(j)+parmu3(k)*(wbar(i,j,k-2)-2*wbar(i,j,k-1))/(dz*dz))*2.0/dz
        END DO
END DO
!!! rhs core: END

! set psi:

CALL pois3d(0,nx-1,C1,0,ny-1,C2,1,nz,A,B,C,nx+1,ny+1,varpsi,IERROR,W)

END SUBROUTINE set_psi4

END MODULE psi_solver
