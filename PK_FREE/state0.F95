MODULE state0

USE parameters, ONLY: nx,ny,nz

IMPLICIT NONE
PUBLIC

! STATE COMPONENTS:

REAL,ALLOCATABLE,DIMENSION(:,:,:) :: pre10,qbar0,the10
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: ubar0,vbar0,wbar0
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: util0,vtil0,wtil0
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: pre20,the1ref0
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: prec0

CONTAINS

SUBROUTINE allocate_state0()

IF (.NOT.ALLOCATED(pre10)) THEN
        ALLOCATE(pre10(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(qbar0)) THEN
        ALLOCATE(qbar0(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(the10)) THEN
        ALLOCATE(the10(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ubar0)) THEN
        ALLOCATE(ubar0(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(vbar0)) THEN
        ALLOCATE(vbar0(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(wbar0)) THEN
        ALLOCATE(wbar0(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(util0)) THEN
        ALLOCATE(util0(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(vtil0)) THEN
        ALLOCATE(vtil0(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(wtil0)) THEN
        ALLOCATE(wtil0(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(pre20)) THEN
        ALLOCATE(pre20(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(the1ref0)) THEN
        ALLOCATE(the1ref0(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(prec0)) THEN
        ALLOCATE(prec0(0:nx,0:ny,nz))
END IF

END SUBROUTINE allocate_state0

SUBROUTINE check_borders0()

INTEGER :: i,j,k

REAL :: aux_pre1,aux_qbar,aux_the1
REAL :: aux_ubar,aux_vbar,aux_wbar
REAL :: aux_util,aux_vtil,aux_wtil
REAL :: aux_pre2

aux_wbar=0
aux_wtil=0
DO j=0,ny
        DO i=0,nx
                aux_wbar=aux_wbar+ABS(wbar0(i,j,1))
                aux_wtil=aux_wtil+ABS(wtil0(i,j,1))
        END DO
END DO

IF (aux_wbar>0) THEN
        WRITE(*,*) 'Call to EXIT: wbar0 must be 0 at z=0.'
        CALL exit(0)
END IF

IF (aux_wtil>0) THEN
        WRITE(*,*) 'Call to EXIT: wtil0 must be 0 at z=0.'
        CALL exit(0)
END IF

aux_the1=0
aux_wbar=0
aux_wtil=0
DO j=0,ny
        DO i=0,nx
                aux_the1=aux_the1+ABS(the10(i,j,nz))
                aux_wbar=aux_wbar+ABS(wbar0(i,j,nz))
                aux_wtil=aux_wtil+ABS(wtil0(i,j,nz))
        END DO
END DO

IF (aux_the1>0) THEN
        WRITE(*,*) 'Call to EXIT: the10 must be 0 at z=Lz.'
        CALL exit(0)
END IF

IF (aux_wbar>0) THEN
        WRITE(*,*) 'Call to EXIT: wbar0 must be 0 at z=Lz.'
        CALL exit(0)
END IF

IF (aux_wtil>0) THEN
        WRITE(*,*) 'Call to EXIT: wtil0 must be 0 at z=Lz.'
        CALL exit(0)
END IF

aux_pre1=0
aux_qbar=0
aux_the1=0
aux_ubar=0
aux_vbar=0
aux_wbar=0
aux_util=0
aux_vtil=0
aux_wtil=0
aux_pre2=0
DO k=1,nz
        DO j=0,ny
                aux_pre1=aux_pre1+ABS(pre10(0,j,k)-pre10(nx-1,j,k))
                aux_pre1=aux_pre1+ABS(pre10(1,j,k)-pre10(nx,j,k))
                aux_qbar=aux_qbar+ABS(qbar0(0,j,k)-qbar0(nx-1,j,k))
                aux_qbar=aux_qbar+ABS(qbar0(1,j,k)-qbar0(nx,j,k))
                aux_the1=aux_the1+ABS(the10(0,j,k)-the10(nx-1,j,k))
                aux_the1=aux_the1+ABS(the10(1,j,k)-the10(nx,j,k))
                aux_ubar=aux_ubar+ABS(ubar0(0,j,k)-ubar0(nx-1,j,k))
                aux_ubar=aux_ubar+ABS(ubar0(1,j,k)-ubar0(nx,j,k))
                aux_vbar=aux_vbar+ABS(vbar0(0,j,k)-vbar0(nx-1,j,k))
                aux_vbar=aux_vbar+ABS(vbar0(1,j,k)-vbar0(nx,j,k))
                aux_wbar=aux_wbar+ABS(wbar0(0,j,k)-wbar0(nx-1,j,k))
                aux_wbar=aux_wbar+ABS(wbar0(1,j,k)-wbar0(nx,j,k))
                aux_util=aux_util+ABS(util0(0,j,k)-util0(nx-1,j,k))
                aux_util=aux_util+ABS(util0(1,j,k)-util0(nx,j,k))
                aux_vtil=aux_vtil+ABS(vtil0(0,j,k)-vtil0(nx-1,j,k))
                aux_vtil=aux_vtil+ABS(vtil0(1,j,k)-vtil0(nx,j,k))
                aux_wtil=aux_wtil+ABS(wtil0(0,j,k)-wtil0(nx-1,j,k))
                aux_wtil=aux_wtil+ABS(wtil0(1,j,k)-wtil0(nx,j,k))
                aux_pre2=aux_pre2+ABS(pre20(0,j,k)-pre20(nx-1,j,k))
                aux_pre2=aux_pre2+ABS(pre20(1,j,k)-pre20(nx,j,k))
        END DO
END DO

IF (aux_pre1>0) THEN
        WRITE(*,*) 'Call to EXIT: pre10 must be periodic in x.'
        CALL exit(0)
END IF

IF (aux_qbar>0) THEN
        WRITE(*,*) 'Call to EXIT: qbar0 must be periodic in x.'
        CALL exit(0)
END IF

IF (aux_the1>0) THEN
        WRITE(*,*) 'Call to EXIT: the10 must be periodic in x.'
        CALL exit(0)
END IF

IF (aux_ubar>0) THEN
        WRITE(*,*) 'Call to EXIT: ubar0 must be periodic in x.'
        CALL exit(0)
END IF

IF (aux_vbar>0) THEN
        WRITE(*,*) 'Call to EXIT: vbar0 must be periodic in x.'
        CALL exit(0)
END IF

IF (aux_wbar>0) THEN
        WRITE(*,*) 'Call to EXIT: wbar0 must be periodic in x.'
        CALL exit(0)
END IF

IF (aux_util>0) THEN
        WRITE(*,*) 'Call to EXIT: util0 must be periodic in x.'
        CALL exit(0)
END IF

IF (aux_vtil>0) THEN
        WRITE(*,*) 'Call to EXIT: vtil0 must be periodic in x.'
        CALL exit(0)
END IF

IF (aux_wtil>0) THEN
        WRITE(*,*) 'Call to EXIT: wtil0 must be periodic in x.'
        CALL exit(0)
END IF

IF (aux_pre2>0) THEN
        WRITE(*,*) 'Call to EXIT: pre20 must be periodic in x.'
        CALL exit(0)
END IF

aux_pre1=0
aux_qbar=0
aux_the1=0
aux_ubar=0
aux_vbar=0
aux_wbar=0
aux_util=0
aux_vtil=0
aux_wtil=0
aux_pre2=0
DO k=1,nz
        DO i=0,nx
                aux_pre1=aux_pre1+ABS(pre10(i,0,k)-pre10(i,ny-1,k))
                aux_pre1=aux_pre1+ABS(pre10(i,1,k)-pre10(i,ny,k))
                aux_qbar=aux_qbar+ABS(qbar0(i,0,k)-qbar0(i,ny-1,k))
                aux_qbar=aux_qbar+ABS(qbar0(i,1,k)-qbar0(i,ny,k))
                aux_the1=aux_the1+ABS(the10(i,0,k)-the10(i,ny-1,k))
                aux_the1=aux_the1+ABS(the10(i,1,k)-the10(i,ny,k))
                aux_ubar=aux_ubar+ABS(ubar0(i,0,k)-ubar0(i,ny-1,k))
                aux_ubar=aux_ubar+ABS(ubar0(i,1,k)-ubar0(i,ny,k))
                aux_vbar=aux_vbar+ABS(vbar0(i,0,k)-vbar0(i,ny-1,k))
                aux_vbar=aux_vbar+ABS(vbar0(i,1,k)-vbar0(i,ny,k))
                aux_wbar=aux_wbar+ABS(wbar0(i,0,k)-wbar0(i,ny-1,k))
                aux_wbar=aux_wbar+ABS(wbar0(i,1,k)-wbar0(i,ny,k))
                aux_util=aux_util+ABS(util0(i,0,k)-util0(i,ny-1,k))
                aux_util=aux_util+ABS(util0(i,1,k)-util0(i,ny,k))
                aux_vtil=aux_vtil+ABS(vtil0(i,0,k)-vtil0(i,ny-1,k))
                aux_vtil=aux_vtil+ABS(vtil0(i,1,k)-vtil0(i,ny,k))
                aux_wtil=aux_wtil+ABS(wtil0(i,0,k)-wtil0(i,ny-1,k))
                aux_wtil=aux_wtil+ABS(wtil0(i,1,k)-wtil0(i,ny,k))
                aux_pre2=aux_pre2+ABS(pre20(i,0,k)-pre20(i,ny-1,k))
                aux_pre2=aux_pre2+ABS(pre20(i,1,k)-pre20(i,ny,k))
        END DO
END DO

IF (aux_pre1>0) THEN
        WRITE(*,*) 'Call to EXIT: pre10 must be periodic in y.'
        CALL exit(0)
END IF

IF (aux_qbar>0) THEN
        WRITE(*,*) 'Call to EXIT: qbar0 must be periodic in y.'
        CALL exit(0)
END IF

IF (aux_the1>0) THEN
        WRITE(*,*) 'Call to EXIT: the10 must be periodic in y.'
        CALL exit(0)
END IF

IF (aux_ubar>0) THEN
        WRITE(*,*) 'Call to EXIT: ubar0 must be periodic in y.'
        CALL exit(0)
END IF

IF (aux_vbar>0) THEN
        WRITE(*,*) 'Call to EXIT: vbar0 must be periodic in y.'
        CALL exit(0)
END IF

IF (aux_wbar>0) THEN
        WRITE(*,*) 'Call to EXIT: wbar0 must be periodic in y.'
        CALL exit(0)
END IF

IF (aux_util>0) THEN
        WRITE(*,*) 'Call to EXIT: util0 must be periodic in y.'
        CALL exit(0)
END IF

IF (aux_vtil>0) THEN
        WRITE(*,*) 'Call to EXIT: vtil0 must be periodic in y.'
        CALL exit(0)
END IF

IF (aux_wtil>0) THEN
        WRITE(*,*) 'Call to EXIT: wtil0 must be periodic in y.'
        CALL exit(0)
END IF

IF (aux_pre2>0) THEN
        WRITE(*,*) 'Call to EXIT: pre20 must be periodic in y.'
        CALL exit(0)
END IF

END SUBROUTINE check_borders0

SUBROUTINE check_borders_the1ref0()

INTEGER :: i,j,k

REAL :: aux_the1ref

aux_the1ref=0
DO j=0,ny
        DO i=0,nx
                aux_the1ref=aux_the1ref+ABS(the1ref0(i,j,nz))
        END DO
END DO

IF (aux_the1ref>0) THEN
        WRITE(*,*) 'Call to EXIT: the1ref0 must be 0 at z=Lz.'
        CALL exit(0)
END IF

aux_the1ref=0
DO k=1,nz
        DO j=0,ny
                aux_the1ref=aux_the1ref+ABS(the1ref0(0,j,k)-the1ref0(nx-1,j,k))
                aux_the1ref=aux_the1ref+ABS(the1ref0(1,j,k)-the1ref0(nx,j,k))
        END DO
END DO

IF (aux_the1ref>0) THEN
        WRITE(*,*) 'Call to EXIT: the1ref0 must be periodic in x.'
        CALL exit(0)
END IF

aux_the1ref=0
DO k=1,nz
        DO i=0,nx
                aux_the1ref=aux_the1ref+ABS(the1ref0(i,0,k)-the1ref0(i,ny-1,k))
                aux_the1ref=aux_the1ref+ABS(the1ref0(i,1,k)-the1ref0(i,ny,k))
        END DO
END DO

IF (aux_the1ref>0) THEN
        WRITE(*,*) 'Call to EXIT: the1ref0 must be periodic in y.'
        CALL exit(0)
END IF

END SUBROUTINE check_borders_the1ref0

END MODULE state0
