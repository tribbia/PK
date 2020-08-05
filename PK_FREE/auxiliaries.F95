MODULE auxiliaries

USE parameters, ONLY: nx,ny,nz

IMPLICIT NONE
PUBLIC

REAL,ALLOCATABLE,DIMENSION(:,:,:) :: pre1_,qbar_,the1_
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: ubar_,vbar_,wbar_
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: util_,vtil_,wtil_
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: pre2_

REAL,ALLOCATABLE,DIMENSION(:,:,:) :: pre10_,qbar0_,the10_
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: ubar0_,vbar0_,wbar0_
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: util0_,vtil0_,wtil0_
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: pre20_

REAL,ALLOCATABLE,DIMENSION(:,:,:) :: rhoutil,rhovtil,rhowtil

REAL,ALLOCATABLE,DIMENSION(:,:,:) :: dqbar,dthe1
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: drhoutil,drhovtil,drhowtil

REAL,ALLOCATABLE,DIMENSION(:,:,:) :: fqbar1,fthe11
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: frhoutil1,frhovtil1,frhowtil1

REAL,ALLOCATABLE,DIMENSION(:,:,:) :: fqbar2,fthe12
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: frhoutil2,frhovtil2,frhowtil2

REAL,ALLOCATABLE,DIMENSION(:,:,:) :: fqbar3,fthe13
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: frhoutil3,frhovtil3,frhowtil3

REAL,ALLOCATABLE,DIMENSION(:,:,:) :: fqbar4,fthe14
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: frhoutil4,frhovtil4,frhowtil4

REAL,ALLOCATABLE,DIMENSION(:,:,:) :: ff1,ff2,ff3

REAL,ALLOCATABLE,DIMENSION(:,:,:) :: ffqbar,ffthe1
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: ffubar,ffvbar,ffwbar
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: ffrhoutil,ffrhovtil,ffrhowtil

REAL,ALLOCATABLE,DIMENSION(:,:,:) :: ffqbar1,ffthe11
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: ffubar1,ffvbar1,ffwbar1
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: ffrhoutil1,ffrhovtil1,ffrhowtil1

REAL,ALLOCATABLE,DIMENSION(:,:,:) :: ffqbar2,ffthe12
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: ffubar2,ffvbar2,ffwbar2
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: ffrhoutil2,ffrhovtil2,ffrhowtil2

REAL,ALLOCATABLE,DIMENSION(:,:,:) :: ffqbar3,ffthe13
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: ffubar3,ffvbar3,ffwbar3
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: ffrhoutil3,ffrhovtil3,ffrhowtil3

REAL,ALLOCATABLE,DIMENSION(:,:,:) :: ffubar_,ffvbar_,ffwbar_
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: ffrhoutil_,ffrhovtil_,ffrhowtil_

REAL,ALLOCATABLE,DIMENSION(:,:,:) :: ffubar1_,ffvbar1_,ffwbar1_
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: ffrhoutil1_,ffrhovtil1_,ffrhowtil1_

REAL,ALLOCATABLE,DIMENSION(:,:,:) :: ffubar2_,ffvbar2_,ffwbar2_
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: ffrhoutil2_,ffrhovtil2_,ffrhowtil2_

REAL,ALLOCATABLE,DIMENSION(:,:,:) :: ffubar3_,ffvbar3_,ffwbar3_
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: ffrhoutil3_,ffrhovtil3_,ffrhowtil3_

CONTAINS

SUBROUTINE allocate_auxiliaries()

IF (.NOT.ALLOCATED(pre1_)) THEN
        ALLOCATE(pre1_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(qbar_)) THEN
        ALLOCATE(qbar_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(the1_)) THEN
        ALLOCATE(the1_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ubar_)) THEN
        ALLOCATE(ubar_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(vbar_)) THEN
        ALLOCATE(vbar_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(wbar_)) THEN
        ALLOCATE(wbar_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(util_)) THEN
        ALLOCATE(util_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(vtil_)) THEN
        ALLOCATE(vtil_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(wtil_)) THEN
        ALLOCATE(wtil_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(pre2_)) THEN
        ALLOCATE(pre2_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(pre10_)) THEN
        ALLOCATE(pre10_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(qbar0_)) THEN
        ALLOCATE(qbar0_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(the10_)) THEN
        ALLOCATE(the10_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ubar0_)) THEN
        ALLOCATE(ubar0_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(vbar0_)) THEN
        ALLOCATE(vbar0_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(wbar0_)) THEN
        ALLOCATE(wbar0_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(util0_)) THEN
        ALLOCATE(util0_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(vtil0_)) THEN
        ALLOCATE(vtil0_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(wtil0_)) THEN
        ALLOCATE(wtil0_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(pre20_)) THEN
        ALLOCATE(pre20_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(rhoutil)) THEN
        ALLOCATE(rhoutil(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(rhovtil)) THEN
        ALLOCATE(rhovtil(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(rhowtil)) THEN
        ALLOCATE(rhowtil(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(dqbar)) THEN
        ALLOCATE(dqbar(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(dthe1)) THEN
        ALLOCATE(dthe1(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(drhoutil)) THEN
        ALLOCATE(drhoutil(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(drhovtil)) THEN
        ALLOCATE(drhovtil(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(drhowtil)) THEN
        ALLOCATE(drhowtil(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(fqbar1)) THEN
        ALLOCATE(fqbar1(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(fthe11)) THEN
        ALLOCATE(fthe11(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(frhoutil1)) THEN
        ALLOCATE(frhoutil1(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(frhovtil1)) THEN
        ALLOCATE(frhovtil1(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(frhowtil1)) THEN
        ALLOCATE(frhowtil1(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(fqbar2)) THEN
        ALLOCATE(fqbar2(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(fthe12)) THEN
        ALLOCATE(fthe12(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(frhoutil2)) THEN
        ALLOCATE(frhoutil2(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(frhovtil2)) THEN
        ALLOCATE(frhovtil2(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(frhowtil2)) THEN
        ALLOCATE(frhowtil2(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(fqbar3)) THEN
        ALLOCATE(fqbar3(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(fthe13)) THEN
        ALLOCATE(fthe13(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(frhoutil3)) THEN
        ALLOCATE(frhoutil3(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(frhovtil3)) THEN
        ALLOCATE(frhovtil3(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(frhowtil3)) THEN
        ALLOCATE(frhowtil3(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(fqbar4)) THEN
        ALLOCATE(fqbar4(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(fthe14)) THEN
        ALLOCATE(fthe14(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(frhoutil4)) THEN
        ALLOCATE(frhoutil4(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(frhovtil4)) THEN
        ALLOCATE(frhovtil4(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(frhowtil4)) THEN
        ALLOCATE(frhowtil4(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ff1)) THEN
        ALLOCATE(ff1(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ff2)) THEN
        ALLOCATE(ff2(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ff3)) THEN
        ALLOCATE(ff3(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffqbar)) THEN
        ALLOCATE(ffqbar(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffthe1)) THEN
        ALLOCATE(ffthe1(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffubar)) THEN
        ALLOCATE(ffubar(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffvbar)) THEN
        ALLOCATE(ffvbar(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffwbar)) THEN
        ALLOCATE(ffwbar(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffrhoutil)) THEN
        ALLOCATE(ffrhoutil(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffrhovtil)) THEN
        ALLOCATE(ffrhovtil(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffrhowtil)) THEN
        ALLOCATE(ffrhowtil(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffqbar1)) THEN
        ALLOCATE(ffqbar1(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffthe11)) THEN
        ALLOCATE(ffthe11(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffubar1)) THEN
        ALLOCATE(ffubar1(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffvbar1)) THEN
        ALLOCATE(ffvbar1(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffwbar1)) THEN
        ALLOCATE(ffwbar1(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffrhoutil1)) THEN
        ALLOCATE(ffrhoutil1(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffrhovtil1)) THEN
        ALLOCATE(ffrhovtil1(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffrhowtil1)) THEN
        ALLOCATE(ffrhowtil1(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffqbar2)) THEN
        ALLOCATE(ffqbar2(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffthe12)) THEN
        ALLOCATE(ffthe12(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffubar2)) THEN
        ALLOCATE(ffubar2(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffvbar2)) THEN
        ALLOCATE(ffvbar2(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffwbar2)) THEN
        ALLOCATE(ffwbar2(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffrhoutil2)) THEN
        ALLOCATE(ffrhoutil2(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffrhovtil2)) THEN
        ALLOCATE(ffrhovtil2(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffrhowtil2)) THEN
        ALLOCATE(ffrhowtil2(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffqbar3)) THEN
        ALLOCATE(ffqbar3(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffthe13)) THEN
        ALLOCATE(ffthe13(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffubar3)) THEN
        ALLOCATE(ffubar3(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffvbar3)) THEN
        ALLOCATE(ffvbar3(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffwbar3)) THEN
        ALLOCATE(ffwbar3(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffrhoutil3)) THEN
        ALLOCATE(ffrhoutil3(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffrhovtil3)) THEN
        ALLOCATE(ffrhovtil3(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffrhowtil3)) THEN
        ALLOCATE(ffrhowtil3(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffubar_)) THEN
        ALLOCATE(ffubar_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffvbar_)) THEN
        ALLOCATE(ffvbar_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffwbar_)) THEN
        ALLOCATE(ffwbar_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffrhoutil_)) THEN
        ALLOCATE(ffrhoutil_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffrhovtil_)) THEN
        ALLOCATE(ffrhovtil_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffrhowtil_)) THEN
        ALLOCATE(ffrhowtil_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffubar1_)) THEN
        ALLOCATE(ffubar1_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffvbar1_)) THEN
        ALLOCATE(ffvbar1_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffwbar1_)) THEN
        ALLOCATE(ffwbar1_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffrhoutil1_)) THEN
        ALLOCATE(ffrhoutil1_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffrhovtil1_)) THEN
        ALLOCATE(ffrhovtil1_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffrhowtil1_)) THEN
        ALLOCATE(ffrhowtil1_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffubar2_)) THEN
        ALLOCATE(ffubar2_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffvbar2_)) THEN
        ALLOCATE(ffvbar2_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffwbar2_)) THEN
        ALLOCATE(ffwbar2_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffrhoutil2_)) THEN
        ALLOCATE(ffrhoutil2_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffrhovtil2_)) THEN
        ALLOCATE(ffrhovtil2_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffrhowtil2_)) THEN
        ALLOCATE(ffrhowtil2_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffubar3_)) THEN
        ALLOCATE(ffubar3_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffvbar3_)) THEN
        ALLOCATE(ffvbar3_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffwbar3_)) THEN
        ALLOCATE(ffwbar3_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffrhoutil3_)) THEN
        ALLOCATE(ffrhoutil3_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffrhovtil3_)) THEN
        ALLOCATE(ffrhovtil3_(0:nx,0:ny,nz))
END IF

IF (.NOT.ALLOCATED(ffrhowtil3_)) THEN
        ALLOCATE(ffrhowtil3_(0:nx,0:ny,nz))
END IF

END SUBROUTINE allocate_auxiliaries

END MODULE auxiliaries
