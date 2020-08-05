MODULE utilities

USE parameters, ONLY: nx,ny,nz
USE parameters, ONLY: dx,dy,dz
USE parameters, ONLY: parpi
USE parameters, ONLY: Lx,Ly,Lz
USE parameters, ONLY: y,z
USE parameters, ONLY: year,omega
USE parameters, ONLY: parg,parc,dryR
USE parameters, ONLY: kappa,parcp
USE parameters, ONLY: rho00,pre00
USE parameters, ONLY: parT0,parH
USE parameters, ONLY: parLv,parRv,refEs,refT
USE parameters, ONLY: Cd_uvbar,Cd_uvtil
USE parameters, ONLY: C_the1,C_qbar
USE parameters, ONLY: nu12_0,nu3_0
USE parameters, ONLY: mu12_0,mu3_0
USE parameters, ONLY: xi12_0,xi3_0
USE parameters, ONLY: up12_0,up3_0
USE parameters, ONLY: va12_0,va3_0
USE parameters, ONLY: z00
USE parameters, ONLY: sigma

IMPLICIT NONE
PUBLIC

! INTERFACE rho1
        ! MODULE PROCEDURE rho1_global_pre1
        ! MODULE PROCEDURE rho1_local_pre1
! END INTERFACE

CONTAINS

REAL FUNCTION omega2(j)

        INTEGER,INTENT(IN) :: j

        omega2=omega*ABS(COS(2*parpi*y(j)/Ly))

END FUNCTION omega2

REAL FUNCTION omega3(j)

        INTEGER,INTENT(IN) :: j

        omega3=omega*SIN(2*parpi*y(j)/Ly)

END FUNCTION omega3

REAL FUNCTION rho0(k)

        INTEGER,INTENT(IN) :: k

        rho0=rho00*EXP(-z(k)/parH) ! Kg/m^3

END FUNCTION rho0

REAL FUNCTION dzrho0(k)

        INTEGER,INTENT(IN) :: k

        dzrho0=-rho0(k)/parH

END FUNCTION dzrho0

REAL FUNCTION pre0(k)

        INTEGER,INTENT(IN) :: k

        pre0=pre00*EXP(-z(k)/parH) ! Pa

END FUNCTION pre0

REAL FUNCTION the0(k)

        INTEGER,INTENT(IN) :: k

        the0=parT0*(parc/pre0(k))**kappa ! K

END FUNCTION the0

REAL FUNCTION dzthe0(k)

        INTEGER,INTENT(IN) :: k

        dzthe0=kappa*the0(k)/parH

END FUNCTION dzthe0

REAL FUNCTION buoyancy_vint(i,j,k,the1)

        INTEGER,INTENT(IN) :: i,j,k

        REAL,DIMENSION(0:nx,0:ny,nz),INTENT(IN) :: the1

        INTEGER :: kk

        buoyancy_vint=0.0
        DO kk=1,k
                buoyancy_vint=buoyancy_vint+dz*parg*the1(i,j,kk)/the0(kk)
        END DO

END FUNCTION buoyancy_vint

REAL FUNCTION rho1(i,j,k,pre1,qbar,the1)

        INTEGER,INTENT(IN) :: i,j,k

        REAL,DIMENSION(0:nx,0:ny,nz),INTENT(IN) :: pre1,qbar,the1

        REAL :: Tv

        Tv=(((pre0(k)+pre1(i,j,k))/parc)**kappa)*(the0(k)+the1(i,j,k))*(1.0+0.6*qbar(i,j,k))

        rho1=pre1(i,j,k)/(dryR*Tv)+rho0(k)*(parT0/Tv-1.0)

END FUNCTION rho1

REAL FUNCTION precip(i,j,k,pre1,wbar,the1,qbar,wtil)

        INTEGER,INTENT(IN) :: i,j,k

        REAL,DIMENSION(0:nx,0:ny,nz),INTENT(IN) :: pre1,qbar,the1
        REAL,DIMENSION(0:nx,0:ny,nz),INTENT(IN) :: wbar,wtil

        REAL :: pressure,theta,Tbar,Es,qs,aux

        precip=0.0

        IF (z(k)>=z00.AND.wtil(i,j,k)<=0.AND.(wbar(i,j,k)+wtil(i,j,k))>0) THEN
                pressure=pre0(k)+pre1(i,j,k)
                theta=the0(k)+the1(i,j,k)
                Tbar=theta*((pressure/parc)**kappa)
                Es=refEs*EXP((1.0/refT-1.0/Tbar)*parLv/parRv)
                qs=0.392*Es/pressure
                IF (qbar(i,j,k)>=qs) THEN
                aux=((qs*parg)/(dryR*(1.0+0.6*qbar(i,j,k))))*(parLv*dryR-parcp*parRv*Tbar)/(parcp*parRv*(Tbar**2)+qs*(parLv**2))
                        precip=aux*(wbar(i,j,k)+wtil(i,j,k))
                END IF
        END IF

END FUNCTION precip

REAL FUNCTION rad(i,j,k,pre1,the1,wbar,wtil,precipitation)

        INTEGER,INTENT(IN) :: i,j,k

        REAL,DIMENSION(0:nx,0:ny,nz),INTENT(IN) :: pre1,the1
        REAL,DIMENSION(0:nx,0:ny,nz),INTENT(IN) :: wbar,wtil

        REAL,INTENT(IN) :: precipitation

        REAL,PARAMETER :: r0=1.5/(24*3600)
        REAL,PARAMETER :: dzrad=-0.7/(16e3*24*3600)

        rad=(4*sigma/parcp)*(parT0**3)*the1(i,j,k)*(((pre0(k)+pre1(i,j,k))/parc)**kappa)

        IF (precipitation>0) THEN
        ! IF (z(k)<4e3.AND.wbar(i,j,k)+wtil(i,j,k)>1e-5) THEN
                rad=rad+((parc/(pre0(k)+pre1(i,j,k)))**kappa)*(r0+dzrad*z(k))
        END IF

END FUNCTION rad

REAL FUNCTION parmu12(k)

        INTEGER,INTENT(IN) :: k

        parmu12=mu12_0

END FUNCTION parmu12

REAL FUNCTION parmu3(k)

        INTEGER,INTENT(IN) :: k

        parmu3=mu3_0

END FUNCTION parmu3

REAL FUNCTION parnu12(k)

        INTEGER,INTENT(IN) :: k

        parnu12=nu12_0/rho0(k) ! m^2/s

END FUNCTION parnu12

REAL FUNCTION parnu3(k)

        INTEGER,INTENT(IN) :: k

        parnu3=nu3_0/rho0(k) ! m^2/s

END FUNCTION parnu3

REAL FUNCTION parxi12(k)

        INTEGER,INTENT(IN) :: k

        parxi12=xi12_0

END FUNCTION parxi12

REAL FUNCTION parxi3(k)

        INTEGER,INTENT(IN) :: k

        parxi3=xi3_0

        ! IF ((k-1)*dz<=2e3) THEN
                ! parxi3=xi3_0
        ! ELSE
                ! parxi3=0.0
        ! END IF

END FUNCTION parxi3

REAL FUNCTION parup12(k)

        INTEGER,INTENT(IN) :: k

        parup12=up12_0

END FUNCTION parup12

REAL FUNCTION parup3(k)

        INTEGER,INTENT(IN) :: k

        parup3=up3_0

END FUNCTION parup3

REAL FUNCTION parva12(k)

        INTEGER,INTENT(IN) :: k

        parva12=va12_0

END FUNCTION parva12

REAL FUNCTION parva3(k)

        INTEGER,INTENT(IN) :: k

        parva3=va3_0

END FUNCTION parva3

REAL FUNCTION mudzubar(i,j,ubar,vbar)

        INTEGER,INTENT(IN) :: i,j

        REAL,DIMENSION(0:nx,0:ny,nz),INTENT(IN) :: ubar,vbar

        mudzubar=Cd_uvbar*SQRT(ubar(i,j,1)**2+vbar(i,j,1)**2)*ubar(i,j,1)

END FUNCTION mudzubar

REAL FUNCTION mudzvbar(i,j,ubar,vbar)

        INTEGER,INTENT(IN) :: i,j

        REAL,DIMENSION(0:nx,0:ny,nz),INTENT(IN) :: ubar,vbar

        mudzvbar=Cd_uvbar*SQRT(ubar(i,j,1)**2+vbar(i,j,1)**2)*vbar(i,j,1)

END FUNCTION mudzvbar

REAL FUNCTION nudzubar(i,j,ubar,vbar)

        INTEGER,INTENT(IN) :: i,j

        REAL,DIMENSION(0:nx,0:ny,nz),INTENT(IN) :: ubar,vbar

        nudzubar=mudzubar(i,j,ubar,vbar)/rho0(1)

END FUNCTION nudzubar

REAL FUNCTION nudzvbar(i,j,ubar,vbar)

        INTEGER,INTENT(IN) :: i,j

        REAL,DIMENSION(0:nx,0:ny,nz),INTENT(IN) :: ubar,vbar

        nudzvbar=mudzvbar(i,j,ubar,vbar)/rho0(1)

END FUNCTION nudzvbar

REAL FUNCTION ximudzubar(i,j,ubar,vbar)

        INTEGER,INTENT(IN) :: i,j

        REAL,DIMENSION(0:nx,0:ny,nz),INTENT(IN) :: ubar,vbar

        ximudzubar=mudzubar(i,j,ubar,vbar)*(parxi3(1)-parmu3(1))/parxi3(1)

END FUNCTION ximudzubar

REAL FUNCTION ximudzvbar(i,j,ubar,vbar)

        INTEGER,INTENT(IN) :: i,j

        REAL,DIMENSION(0:nx,0:ny,nz),INTENT(IN) :: ubar,vbar

        ximudzvbar=mudzvbar(i,j,ubar,vbar)*(parxi3(1)-parmu3(1))/parxi3(1)

END FUNCTION ximudzvbar

REAL FUNCTION xidzutil(i,j,ubar,vbar,util)

        INTEGER,INTENT(IN) :: i,j

        REAL,DIMENSION(0:nx,0:ny,nz),INTENT(IN) :: ubar,vbar,util

        xidzutil=Cd_uvtil*SQRT(ubar(i,j,1)**2+vbar(i,j,1)**2)*util(i,j,1)

END FUNCTION xidzutil

REAL FUNCTION xidzvtil(i,j,ubar,vbar,vtil)

        INTEGER,INTENT(IN) :: i,j

        REAL,DIMENSION(0:nx,0:ny,nz),INTENT(IN) :: ubar,vbar,vtil

        xidzvtil=Cd_uvtil*SQRT(ubar(i,j,1)**2+vbar(i,j,1)**2)*vtil(i,j,1)

END FUNCTION xidzvtil

REAL FUNCTION updzthe1(i,j,ubar,vbar,pre1,the1,time)

        INTEGER,INTENT(IN) :: i,j

        REAL,DIMENSION(0:nx,0:ny,nz),INTENT(IN) :: ubar,vbar
        REAL,DIMENSION(0:nx,0:ny,nz),INTENT(IN) :: pre1,the1
        REAL,INTENT(IN) :: time

        updzthe1=C_the1*SQRT(ubar(i,j,1)**2+vbar(i,j,1)**2)*(the1(i,j,1)-ground_the1(i,j,time,pre1))

END FUNCTION updzthe1

REAL FUNCTION vadzqbar(i,j,ubar,vbar,qbar,time)

        INTEGER,INTENT(IN) :: i,j

        REAL,DIMENSION(0:nx,0:ny,nz),INTENT(IN) :: ubar,vbar,qbar
        REAL,INTENT(IN) :: time

        vadzqbar=C_qbar*SQRT(ubar(i,j,1)**2+vbar(i,j,1)**2)*(qbar(i,j,1)-ground_qbar(i,j,time))

END FUNCTION vadzqbar

REAL FUNCTION ground_the1(i,j,time,pre1)

        INTEGER,INTENT(IN) :: i,j

        REAL,INTENT(IN) :: time

        REAL,DIMENSION(0:nx,0:ny,nz),INTENT(IN) :: pre1

        INTEGER :: k=1
        REAL,PARAMETER :: amp1=10 ! (1-EXP(-(.00005*time)**2))*81*2/3.1416

        ground_the1=((20+amp1*(1+COS((2*parpi/(nx-1))*(i+((nx-1)/year)*(time+year/2)))))*EXP(-0.5*((j-0.5*ny)/(ny/4.0))**2)+267) &
                    *((parc/(pre0(k)+pre1(i,j,k)))**kappa)-the0(k)

END FUNCTION ground_the1

REAL FUNCTION ground_qbar(i,j,time)

        INTEGER,INTENT(IN) :: i,j

        REAL,INTENT(IN) :: time

        REAL,PARAMETER :: amp2=ny/8.0

        ground_qbar=10e-3+15e-3*EXP(-0.5*((j-(0.5*ny+amp2*SIN(2*parpi*0.0/year)))/(ny/3.0))**2)                                  &
                    *(0.5+0.25*(1+COS((16*parpi/(nx-1))*(i-1))))

END FUNCTION ground_qbar

! REAL FUNCTION top_the1(i,j,time,pre1)

        ! INTEGER,INTENT(IN) :: i,j

        ! REAL,INTENT(IN) :: time

        ! REAL,DIMENSION(0:nx,0:ny,nz),INTENT(IN) :: pre1

        ! INTEGER :: jj,k=nz
        ! REAL,PARAMETER :: amp1=51
        ! REAL,PARAMETER :: amp2=ny/16.0

        ! SELECT CASE (j)

        ! CASE (1:(ny-1)/4)

        ! jj=(ny-1)/2-j+2
        ! top_the1=(250+amp1*EXP(-0.5*(z(k)/(10e3/2.2))**2)*EXP(-0.5*((jj-(0.5*ny+amp2*SIN(2*parpi*time/year)))/(ny/11.25))**2)    &
                    ! *EXP(-0.5*((i-0.5*nx)/(nx/1.2))**2))*((parc/(pre0(k)+pre1(i,jj,k)))**kappa)-the0(k)

        ! CASE ((ny-1)/4+1:3*(ny-1)/4+1)

        ! top_the1=(250+amp1*EXP(-0.5*(z(k)/(10e3/2.2))**2)*EXP(-0.5*((j-(0.5*ny+amp2*SIN(2*parpi*time/year)))/(ny/11.25))**2)     &
                    ! *EXP(-0.5*((i-0.5*nx)/(nx/1.2))**2))*((parc/(pre0(k)+pre1(i,j,k)))**kappa)-the0(k)

        ! CASE (3*(ny-1)/4+2:)

        ! jj=3*(ny-1)/2-j+2
        ! top_the1=(250+amp1*EXP(-0.5*(z(k)/(10e3/2.2))**2)*EXP(-0.5*((jj-(0.5*ny+amp2*SIN(2*parpi*time/year)))/(ny/11.25))**2)    &
                    ! *EXP(-0.5*((i-0.5*nx)/(nx/1.2))**2))*((parc/(pre0(k)+pre1(i,jj,k)))**kappa)-the0(k)

        ! END SELECT

! END FUNCTION top_the1

END MODULE utilities
