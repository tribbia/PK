MODULE state

USE parameters, ONLY: nx,ny,nz

IMPLICIT NONE
PUBLIC

! STATE COMPONENTS:

REAL,ALLOCATABLE,DIMENSION(:,:,:) :: pre1,qbar,the1
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: ubar,vbar,wbar
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: util,vtil,wtil
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: pre2,the1ref
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: prec

CONTAINS

SUBROUTINE allocate_state()

IF (.NOT.ALLOCATED(pre1)) THEN
        ALLOCATE(pre1(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(qbar)) THEN
        ALLOCATE(qbar(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(the1)) THEN
        ALLOCATE(the1(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ubar)) THEN
        ALLOCATE(ubar(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(vbar)) THEN
        ALLOCATE(vbar(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(wbar)) THEN
        ALLOCATE(wbar(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(util)) THEN
        ALLOCATE(util(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(vtil)) THEN
        ALLOCATE(vtil(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(wtil)) THEN
        ALLOCATE(wtil(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(pre2)) THEN
        ALLOCATE(pre2(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(the1ref)) THEN
        ALLOCATE(the1ref(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(prec)) THEN
        ALLOCATE(prec(0:nx,0:ny,nz))
END IF

END SUBROUTINE allocate_state

SUBROUTINE check_borders()

INTEGER :: i,j,k

REAL :: aux_pre1,aux_qbar,aux_the1
REAL :: aux_ubar,aux_vbar,aux_wbar
REAL :: aux_util,aux_vtil,aux_wtil
REAL :: aux_pre2

aux_wbar=0
aux_wtil=0
DO j=0,ny
        DO i=0,nx
                aux_wbar=aux_wbar+ABS(wbar(i,j,1))
                aux_wtil=aux_wtil+ABS(wtil(i,j,1))
        END DO
END DO

IF (aux_wbar>0) THEN
        WRITE(*,*) 'Call to EXIT: wbar must be 0 at z=0.'
        CALL exit(0)
END IF

IF (aux_wtil>0) THEN
        WRITE(*,*) 'Call to EXIT: wtil must be 0 at z=0.'
        CALL exit(0)
END IF

aux_the1=0
aux_wbar=0
aux_wtil=0
DO j=0,ny
        DO i=0,nx
                aux_the1=aux_the1+ABS(the1(i,j,nz))
                aux_wbar=aux_wbar+ABS(wbar(i,j,nz))
                aux_wtil=aux_wtil+ABS(wtil(i,j,nz))
        END DO
END DO

IF (aux_the1>0) THEN
        WRITE(*,*) 'Call to EXIT: the1 must be 0 at z=Lz.'
        CALL exit(0)
END IF

IF (aux_wbar>0) THEN
        WRITE(*,*) 'Call to EXIT: wbar must be 0 at z=Lz.'
        CALL exit(0)
END IF

IF (aux_wtil>0) THEN
        WRITE(*,*) 'Call to EXIT: wtil must be 0 at z=Lz.'
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
                aux_pre1=aux_pre1+ABS(pre1(0,j,k)-pre1(nx-1,j,k))
                aux_pre1=aux_pre1+ABS(pre1(1,j,k)-pre1(nx,j,k))
                aux_qbar=aux_qbar+ABS(qbar(0,j,k)-qbar(nx-1,j,k))
                aux_qbar=aux_qbar+ABS(qbar(1,j,k)-qbar(nx,j,k))
                aux_the1=aux_the1+ABS(the1(0,j,k)-the1(nx-1,j,k))
                aux_the1=aux_the1+ABS(the1(1,j,k)-the1(nx,j,k))
                aux_ubar=aux_ubar+ABS(ubar(0,j,k)-ubar(nx-1,j,k))
                aux_ubar=aux_ubar+ABS(ubar(1,j,k)-ubar(nx,j,k))
                aux_vbar=aux_vbar+ABS(vbar(0,j,k)-vbar(nx-1,j,k))
                aux_vbar=aux_vbar+ABS(vbar(1,j,k)-vbar(nx,j,k))
                aux_wbar=aux_wbar+ABS(wbar(0,j,k)-wbar(nx-1,j,k))
                aux_wbar=aux_wbar+ABS(wbar(1,j,k)-wbar(nx,j,k))
                aux_util=aux_util+ABS(util(0,j,k)-util(nx-1,j,k))
                aux_util=aux_util+ABS(util(1,j,k)-util(nx,j,k))
                aux_vtil=aux_vtil+ABS(vtil(0,j,k)-vtil(nx-1,j,k))
                aux_vtil=aux_vtil+ABS(vtil(1,j,k)-vtil(nx,j,k))
                aux_wtil=aux_wtil+ABS(wtil(0,j,k)-wtil(nx-1,j,k))
                aux_wtil=aux_wtil+ABS(wtil(1,j,k)-wtil(nx,j,k))
                aux_pre2=aux_pre2+ABS(pre2(0,j,k)-pre2(nx-1,j,k))
                aux_pre2=aux_pre2+ABS(pre2(1,j,k)-pre2(nx,j,k))
        END DO
END DO

IF (aux_pre1>0) THEN
        WRITE(*,*) 'Call to EXIT: pre1 must be periodic in x.'
        CALL exit(0)
END IF

IF (aux_qbar>0) THEN
        WRITE(*,*) 'Call to EXIT: qbar must be periodic in x.'
        CALL exit(0)
END IF

IF (aux_the1>0) THEN
        WRITE(*,*) 'Call to EXIT: the1 must be periodic in x.'
        CALL exit(0)
END IF

IF (aux_ubar>0) THEN
        WRITE(*,*) 'Call to EXIT: ubar must be periodic in x.'
        CALL exit(0)
END IF

IF (aux_vbar>0) THEN
        WRITE(*,*) 'Call to EXIT: vbar must be periodic in x.'
        CALL exit(0)
END IF

IF (aux_wbar>0) THEN
        WRITE(*,*) 'Call to EXIT: wbar must be periodic in x.'
        CALL exit(0)
END IF

IF (aux_util>0) THEN
        WRITE(*,*) 'Call to EXIT: util must be periodic in x.'
        CALL exit(0)
END IF

IF (aux_vtil>0) THEN
        WRITE(*,*) 'Call to EXIT: vtil must be periodic in x.'
        CALL exit(0)
END IF

IF (aux_wtil>0) THEN
        WRITE(*,*) 'Call to EXIT: wtil must be periodic in x.'
        CALL exit(0)
END IF

IF (aux_pre2>0) THEN
        WRITE(*,*) 'Call to EXIT: pre2 must be periodic in x.'
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
                aux_pre1=aux_pre1+ABS(pre1(i,0,k)-pre1(i,ny-1,k))
                aux_pre1=aux_pre1+ABS(pre1(i,1,k)-pre1(i,ny,k))
                aux_qbar=aux_qbar+ABS(qbar(i,0,k)-qbar(i,ny-1,k))
                aux_qbar=aux_qbar+ABS(qbar(i,1,k)-qbar(i,ny,k))
                aux_the1=aux_the1+ABS(the1(i,0,k)-the1(i,ny-1,k))
                aux_the1=aux_the1+ABS(the1(i,1,k)-the1(i,ny,k))
                aux_ubar=aux_ubar+ABS(ubar(i,0,k)-ubar(i,ny-1,k))
                aux_ubar=aux_ubar+ABS(ubar(i,1,k)-ubar(i,ny,k))
                aux_vbar=aux_vbar+ABS(vbar(i,0,k)-vbar(i,ny-1,k))
                aux_vbar=aux_vbar+ABS(vbar(i,1,k)-vbar(i,ny,k))
                aux_wbar=aux_wbar+ABS(wbar(i,0,k)-wbar(i,ny-1,k))
                aux_wbar=aux_wbar+ABS(wbar(i,1,k)-wbar(i,ny,k))
                aux_util=aux_util+ABS(util(i,0,k)-util(i,ny-1,k))
                aux_util=aux_util+ABS(util(i,1,k)-util(i,ny,k))
                aux_vtil=aux_vtil+ABS(vtil(i,0,k)-vtil(i,ny-1,k))
                aux_vtil=aux_vtil+ABS(vtil(i,1,k)-vtil(i,ny,k))
                aux_wtil=aux_wtil+ABS(wtil(i,0,k)-wtil(i,ny-1,k))
                aux_wtil=aux_wtil+ABS(wtil(i,1,k)-wtil(i,ny,k))
                aux_pre2=aux_pre2+ABS(pre2(i,0,k)-pre2(i,ny-1,k))
                aux_pre2=aux_pre2+ABS(pre2(i,1,k)-pre2(i,ny,k))
        END DO
END DO

IF (aux_pre1>0) THEN
        WRITE(*,*) 'Call to EXIT: pre1 must be periodic in y.'
        CALL exit(0)
END IF

IF (aux_qbar>0) THEN
        WRITE(*,*) 'Call to EXIT: qbar must be periodic in y.'
        CALL exit(0)
END IF

IF (aux_the1>0) THEN
        WRITE(*,*) 'Call to EXIT: the1 must be periodic in y.'
        CALL exit(0)
END IF

IF (aux_ubar>0) THEN
        WRITE(*,*) 'Call to EXIT: ubar must be periodic in y.'
        CALL exit(0)
END IF

IF (aux_vbar>0) THEN
        WRITE(*,*) 'Call to EXIT: vbar must be periodic in y.'
        CALL exit(0)
END IF

IF (aux_wbar>0) THEN
        WRITE(*,*) 'Call to EXIT: wbar must be periodic in y.'
        CALL exit(0)
END IF

IF (aux_util>0) THEN
        WRITE(*,*) 'Call to EXIT: util must be periodic in y.'
        CALL exit(0)
END IF

IF (aux_vtil>0) THEN
        WRITE(*,*) 'Call to EXIT: vtil must be periodic in y.'
        CALL exit(0)
END IF

IF (aux_wtil>0) THEN
        WRITE(*,*) 'Call to EXIT: wtil must be periodic in y.'
        CALL exit(0)
END IF

IF (aux_pre2>0) THEN
        WRITE(*,*) 'Call to EXIT: pre2 must be periodic in y.'
        CALL exit(0)
END IF

END SUBROUTINE check_borders

SUBROUTINE check_borders_the1ref()

INTEGER :: i,j,k

REAL :: aux_the1ref

aux_the1ref=0
DO j=0,ny
        DO i=0,nx
                aux_the1ref=aux_the1ref+ABS(the1ref(i,j,nz))
        END DO
END DO

IF (aux_the1ref>0) THEN
        WRITE(*,*) 'Call to EXIT: the1ref must be 0 at z=Lz.'
        CALL exit(0)
END IF

aux_the1ref=0
DO k=1,nz
        DO j=0,ny
                aux_the1ref=aux_the1ref+ABS(the1ref(0,j,k)-the1ref(nx-1,j,k))
                aux_the1ref=aux_the1ref+ABS(the1ref(1,j,k)-the1ref(nx,j,k))
        END DO
END DO

IF (aux_the1ref>0) THEN
        WRITE(*,*) 'Call to EXIT: the1ref must be periodic in x.'
        CALL exit(0)
END IF

aux_the1ref=0
DO k=1,nz
        DO i=0,nx
                aux_the1ref=aux_the1ref+ABS(the1ref(i,0,k)-the1ref(i,ny-1,k))
                aux_the1ref=aux_the1ref+ABS(the1ref(i,1,k)-the1ref(i,ny,k))
        END DO
END DO

IF (aux_the1ref>0) THEN
        WRITE(*,*) 'Call to EXIT: the1ref must be periodic in y.'
        CALL exit(0)
END IF

END SUBROUTINE check_borders_the1ref

END MODULE state
