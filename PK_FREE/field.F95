MODULE field

USE parameters, ONLY: nx,ny,nz
USE parameters, ONLY: dx,dy,dz
USE parameters, ONLY: parLv,parcp
USE parameters, ONLY: parc,parg,parK
USE parameters, ONLY: kappa

USE utilities, ONLY: rho1,precip,rad
USE utilities, ONLY: omega2,omega3
USE utilities, ONLY: pre0
USE utilities, ONLY: rho0
USE utilities, ONLY: the0,dzthe0
USE utilities, ONLY: parmu12,parmu3
USE utilities, ONLY: parnu12,parnu3
USE utilities, ONLY: parup12,parup3
USE utilities, ONLY: parva12,parva3
USE utilities, ONLY: parxi12,parxi3

USE utilities, ONLY: nudzubar,nudzvbar
USE utilities, ONLY: ximudzubar,ximudzvbar
USE utilities, ONLY: xidzutil,xidzvtil
USE utilities, ONLY: updzthe1,vadzqbar

IMPLICIT NONE
PUBLIC

CONTAINS

SUBROUTINE set_field1(delta,time,pre1,qbar,the1,ubar,vbar,wbar,util,vtil,wtil,pre2,the1ref,                                      &
                     fqbar,fthe1,fubar,fvbar,fwbar,frhoutil,frhovtil,frhowtil)

REAL,INTENT(IN) :: delta
REAL,INTENT(IN) :: time

REAL,DIMENSION(0:nx,0:ny,nz),INTENT(IN) :: pre1,qbar,the1
REAL,DIMENSION(0:nx,0:ny,nz),INTENT(IN) :: ubar,vbar,wbar
REAL,DIMENSION(0:nx,0:ny,nz),INTENT(IN) :: util,vtil,wtil
REAL,DIMENSION(0:nx,0:ny,nz),INTENT(IN) :: pre2,the1ref

REAL,DIMENSION(0:nx,0:ny,nz),INTENT(OUT) :: fqbar,fthe1
REAL,DIMENSION(0:nx,0:ny,nz),INTENT(OUT) :: fubar,fvbar,fwbar
REAL,DIMENSION(0:nx,0:ny,nz),INTENT(OUT) :: frhoutil,frhovtil,frhowtil

! LOCAL VARIABLES:

INTEGER :: i,j,k

REAL :: precipitation,heat

!!! core: BEGIN
DO k=2,nz-1
        DO j=1,ny-1
                DO i=1,nx-1
precipitation=precip(i,j,k,pre1,wbar,the1,qbar,wtil)
heat=((parc/(pre0(k)+pre1(i,j,k)))**kappa)*precipitation*parLv/parcp

fqbar(i,j,k)=                                                                                                                    &
-(ubar(i,j,k)+util(i,j,k))*(qbar(i+1,j,k)-qbar(i-1,j,k))/(2*dx)                                                                  &
-(vbar(i,j,k)+vtil(i,j,k))*(qbar(i,j+1,k)-qbar(i,j-1,k))/(2*dy)                                                                  &
-(wbar(i,j,k)+wtil(i,j,k))*(qbar(i,j,k+1)-qbar(i,j,k-1))/(2*dz)                                                                  &
-precipitation                                                                                                                   &
+parva12(k)*(qbar(i+1,j,k)-2*qbar(i,j,k)+qbar(i-1,j,k))/(dx*dx)                                                                  &
+parva12(k)*(qbar(i,j+1,k)-2*qbar(i,j,k)+qbar(i,j-1,k))/(dy*dy)                                                                  &
+parva3(k)*(qbar(i,j,k+1)-2*qbar(i,j,k)+qbar(i,j,k-1))/(dz*dz)

fthe1(i,j,k)=                                                                                                                    &
-(ubar(i,j,k)+util(i,j,k))*(the1(i+1,j,k)-the1(i-1,j,k))/(2*dx)                                                                  &
-(vbar(i,j,k)+vtil(i,j,k))*(the1(i,j+1,k)-the1(i,j-1,k))/(2*dy)                                                                  &
-(wbar(i,j,k)+wtil(i,j,k))*(the1(i,j,k+1)-the1(i,j,k-1))/(2*dz)                                                                  &
-(wbar(i,j,k)+wtil(i,j,k))*dzthe0(k)                                                                                             &
+heat-rad(i,j,k,pre1,the1,wbar,wtil,precipitation)                                                                               &
+parup12(k)*(the1(i+1,j,k)-2*the1(i,j,k)+the1(i-1,j,k))/(dx*dx)                                                                  &
+parup12(k)*(the1(i,j+1,k)-2*the1(i,j,k)+the1(i,j-1,k))/(dy*dy)                                                                  &
+parup3(k)*(the1(i,j,k+1)-2*the1(i,j,k)+the1(i,j,k-1))/(dz*dz)                                                                   &
-parK*(the1(i,j,k)-the1ref(i,j,k))

fubar(i,j,k)=                                                                                                                    &
-(ubar(i,j,k)+util(i,j,k))*(ubar(i+1,j,k)-ubar(i-1,j,k))/(2*dx)                                                                  &
-(vbar(i,j,k)+vtil(i,j,k))*(ubar(i,j+1,k)-ubar(i,j-1,k))/(2*dy)                                                                  &
-(wbar(i,j,k)+wtil(i,j,k))*(ubar(i,j,k+1)-ubar(i,j,k-1))/(2*dz)                                                                  &
+2*(vbar(i,j,k)*omega3(j)-wbar(i,j,k)*omega2(j))                                                                                 &
-((pre1(i+1,j,k)-pre1(i-1,j,k))/(2*dx))/rho0(k)                                                                                  &
+parnu12(k)*(ubar(i+1,j,k)-2*ubar(i,j,k)+ubar(i-1,j,k))/(dx*dx)                                                                  &
+parnu12(k)*(ubar(i,j+1,k)-2*ubar(i,j,k)+ubar(i,j-1,k))/(dy*dy)                                                                  &
+parnu3(k)*(ubar(i,j,k+1)-2*ubar(i,j,k)+ubar(i,j,k-1))/(dz*dz)

fvbar(i,j,k)=                                                                                                                    &
-(ubar(i,j,k)+util(i,j,k))*(vbar(i+1,j,k)-vbar(i-1,j,k))/(2*dx)                                                                  &
-(vbar(i,j,k)+vtil(i,j,k))*(vbar(i,j+1,k)-vbar(i,j-1,k))/(2*dy)                                                                  &
-(wbar(i,j,k)+wtil(i,j,k))*(vbar(i,j,k+1)-vbar(i,j,k-1))/(2*dz)                                                                  &
-2*ubar(i,j,k)*omega3(j)                                                                                                         &
-((pre1(i,j+1,k)-pre1(i,j-1,k))/(2*dy))/rho0(k)                                                                                  &
+parnu12(k)*(vbar(i+1,j,k)-2*vbar(i,j,k)+vbar(i-1,j,k))/(dx*dx)                                                                  &
+parnu12(k)*(vbar(i,j+1,k)-2*vbar(i,j,k)+vbar(i,j-1,k))/(dy*dy)                                                                  &
+parnu3(k)*(vbar(i,j,k+1)-2*vbar(i,j,k)+vbar(i,j,k-1))/(dz*dz)

fwbar(i,j,k)=                                                                                                                    &
-(ubar(i,j,k)+util(i,j,k))*(wbar(i+1,j,k)-wbar(i-1,j,k))/(2*dx)                                                                  &
-(vbar(i,j,k)+vtil(i,j,k))*(wbar(i,j+1,k)-wbar(i,j-1,k))/(2*dy)                                                                  &
-(wbar(i,j,k)+wtil(i,j,k))*(wbar(i,j,k+1)-wbar(i,j,k-1))/(2*dz)                                                                  &
+2*ubar(i,j,k)*omega2(j)                                                                                                         &
+parg*the1(i,j,k)/the0(k)                                                                                                        &
-((pre1(i,j,k+1)-pre1(i,j,k-1))/(2*dz))/rho0(k)                                                                                  &
+parnu12(k)*(wbar(i+1,j,k)-2*wbar(i,j,k)+wbar(i-1,j,k))/(dx*dx)                                                                  &
+parnu12(k)*(wbar(i,j+1,k)-2*wbar(i,j,k)+wbar(i,j-1,k))/(dy*dy)                                                                  &
+parnu3(k)*(wbar(i,j,k+1)-2*wbar(i,j,k)+wbar(i,j,k-1))/(dz*dz)

frhoutil(i,j,k)=                                                                                                                 &
-((vbar(i,j+1,k)*(rho0(k)+rho1(i,j+1,k,pre1,qbar,the1))*util(i,j+1,k)-vtil(i,j+1,k)*(rho0(k)+rho1(i,j+1,k,pre1,qbar,the1))       &
*ubar(i,j+1,k))                                                                                                                  &
 -(vbar(i,j-1,k)*(rho0(k)+rho1(i,j-1,k,pre1,qbar,the1))*util(i,j-1,k)-vtil(i,j-1,k)*(rho0(k)+rho1(i,j-1,k,pre1,qbar,the1))       &
*ubar(i,j-1,k)))/(2*dy)                                                                                                          &
-((wbar(i,j,k+1)*(rho0(k+1)+rho1(i,j,k+1,pre1,qbar,the1))*util(i,j,k+1)-wtil(i,j,k+1)*(rho0(k+1)+rho1(i,j,k+1,pre1,qbar,the1))   &
*ubar(i,j,k+1))                                                                                                                  &
 -(wbar(i,j,k-1)*(rho0(k-1)+rho1(i,j,k-1,pre1,qbar,the1))*util(i,j,k-1)-wtil(i,j,k-1)*(rho0(k-1)+rho1(i,j,k-1,pre1,qbar,the1))   &
*ubar(i,j,k-1)))/(2*dz)                                                                                                          &
+2*(rho0(k)+rho1(i,j,k,pre1,qbar,the1))*(vtil(i,j,k)*omega3(j)-wtil(i,j,k)*omega2(j))                                            &
+parxi12(k)*(ubar(i+1,j,k)-2*ubar(i,j,k)+ubar(i-1,j,k))/(dx*dx)                                                     &
+parxi12(k)*(ubar(i,j+1,k)-2*ubar(i,j,k)+ubar(i,j-1,k))/(dy*dy)                                                     &
+parxi3(k)*(ubar(i,j,k+1)-2*ubar(i,j,k)+ubar(i,j,k-1))/(dz*dz)                                                       &
+parxi12(k)*(util(i+1,j,k)-2*util(i,j,k)+util(i-1,j,k))/(dx*dx)                                                                  &
+parxi12(k)*(util(i,j+1,k)-2*util(i,j,k)+util(i,j-1,k))/(dy*dy)                                                                  &
+parxi3(k)*(util(i,j,k+1)-2*util(i,j,k)+util(i,j,k-1))/(dz*dz)                                                                   &
-(pre2(i+1,j,k)-pre2(i-1,j,k))/(2*dx)

frhovtil(i,j,k)=                                                                                                                 &
-((ubar(i+1,j,k)*(rho0(k)+rho1(i+1,j,k,pre1,qbar,the1))*vtil(i+1,j,k)-util(i+1,j,k)*(rho0(k)+rho1(i+1,j,k,pre1,qbar,the1))       &
*vbar(i+1,j,k))                                                                                                                  &
 -(ubar(i-1,j,k)*(rho0(k)+rho1(i-1,j,k,pre1,qbar,the1))*vtil(i-1,j,k)-util(i-1,j,k)*(rho0(k)+rho1(i-1,j,k,pre1,qbar,the1))       &
*vbar(i-1,j,k)))/(2*dx)                                                                                                          &
-((wbar(i,j,k+1)*(rho0(k+1)+rho1(i,j,k+1,pre1,qbar,the1))*vtil(i,j,k+1)-wtil(i,j,k+1)*(rho0(k+1)+rho1(i,j,k+1,pre1,qbar,the1))   &
*vbar(i,j,k+1))                                                                                                                  &
 -(wbar(i,j,k-1)*(rho0(k-1)+rho1(i,j,k-1,pre1,qbar,the1))*vtil(i,j,k-1)-wtil(i,j,k-1)*(rho0(k-1)+rho1(i,j,k-1,pre1,qbar,the1))   &
*vbar(i,j,k-1)))/(2*dz)                                                                                                          &
-2*(rho0(k)+rho1(i,j,k,pre1,qbar,the1))*util(i,j,k)*omega3(j)                                                                    &
+parxi12(k)*(vbar(i+1,j,k)-2*vbar(i,j,k)+vbar(i-1,j,k))/(dx*dx)                                                     &
+parxi12(k)*(vbar(i,j+1,k)-2*vbar(i,j,k)+vbar(i,j-1,k))/(dy*dy)                                                     &
+parxi3(k)*(vbar(i,j,k+1)-2*vbar(i,j,k)+vbar(i,j,k-1))/(dz*dz)                                                       &
+parxi12(k)*(vtil(i+1,j,k)-2*vtil(i,j,k)+vtil(i-1,j,k))/(dx*dx)                                                                  &
+parxi12(k)*(vtil(i,j+1,k)-2*vtil(i,j,k)+vtil(i,j-1,k))/(dy*dy)                                                                  &
+parxi3(k)*(vtil(i,j,k+1)-2*vtil(i,j,k)+vtil(i,j,k-1))/(dz*dz)                                                                   &
-(pre2(i,j+1,k)-pre2(i,j-1,k))/(2*dy)

frhowtil(i,j,k)=                                                                                                                 &
-((ubar(i+1,j,k)*(rho0(k)+rho1(i+1,j,k,pre1,qbar,the1))*wtil(i+1,j,k)-util(i+1,j,k)*(rho0(k)+rho1(i+1,j,k,pre1,qbar,the1))       &
*wbar(i+1,j,k))                                                                                                                  &
 -(ubar(i-1,j,k)*(rho0(k)+rho1(i-1,j,k,pre1,qbar,the1))*wtil(i-1,j,k)-util(i-1,j,k)*(rho0(k)+rho1(i-1,j,k,pre1,qbar,the1))       &
*wbar(i-1,j,k)))/(2*dx)                                                                                                          &
-((vbar(i,j+1,k)*(rho0(k)+rho1(i,j+1,k,pre1,qbar,the1))*wtil(i,j+1,k)-vtil(i,j+1,k)*(rho0(k)+rho1(i,j+1,k,pre1,qbar,the1))       &
*wbar(i,j+1,k))                                                                                                                  &
 -(vbar(i,j-1,k)*(rho0(k)+rho1(i,j-1,k,pre1,qbar,the1))*wtil(i,j-1,k)-vtil(i,j-1,k)*(rho0(k)+rho1(i,j-1,k,pre1,qbar,the1))       &
*wbar(i,j-1,k)))/(2*dy)                                                                                                          &
+2*(rho0(k)+rho1(i,j,k,pre1,qbar,the1))*util(i,j,k)*omega2(j)                                                                    &
+delta*(rho0(k)+rho1(i,j,k,pre1,qbar,the1))*parg*(rho1(i,j,k,pre1,qbar,the1)/rho0(k))**2                                         &
+parxi12(k)*(wbar(i+1,j,k)-2*wbar(i,j,k)+wbar(i-1,j,k))/(dx*dx)                                                     &
+parxi12(k)*(wbar(i,j+1,k)-2*wbar(i,j,k)+wbar(i,j-1,k))/(dy*dy)                                                     &
+parxi3(k)*(wbar(i,j,k+1)-2*wbar(i,j,k)+wbar(i,j,k-1))/(dz*dz)                                                       &
+parxi12(k)*(wtil(i+1,j,k)-2*wtil(i,j,k)+wtil(i-1,j,k))/(dx*dx)                                                                  &
+parxi12(k)*(wtil(i,j+1,k)-2*wtil(i,j,k)+wtil(i,j-1,k))/(dy*dy)                                                                  &
+parxi3(k)*(wtil(i,j,k+1)-2*wtil(i,j,k)+wtil(i,j,k-1))/(dz*dz)                                                                   &
-(pre2(i,j,k+1)-pre2(i,j,k-1))/(2*dz)
                END DO
        END DO
END DO

k=1
DO j=1,ny-1
        DO i=1,nx-1
precipitation=precip(i,j,k,pre1,wbar,the1,qbar,wtil)
heat=((parc/(pre0(k)+pre1(i,j,k)))**kappa)*precipitation*parLv/parcp

fqbar(i,j,k)=                                                                                                                    &
-(ubar(i,j,k)+util(i,j,k))*(qbar(i+1,j,k)-qbar(i-1,j,k))/(2*dx)                                                                  &
-(vbar(i,j,k)+vtil(i,j,k))*(qbar(i,j+1,k)-qbar(i,j-1,k))/(2*dy)                                                                  &
-precipitation                                                                                                                   &
+parva12(k)*(qbar(i+1,j,k)-2*qbar(i,j,k)+qbar(i-1,j,k))/(dx*dx)                                                                  &
+parva12(k)*(qbar(i,j+1,k)-2*qbar(i,j,k)+qbar(i,j-1,k))/(dy*dy)                                                                  &
+(parva3(k+1)*(qbar(i,j,k+2)-qbar(i,j,k))/(2*dz)-vadzqbar(i,j,ubar,vbar,qbar,time))/dz

fthe1(i,j,k)=                                                                                                                    &
-(ubar(i,j,k)+util(i,j,k))*(the1(i+1,j,k)-the1(i-1,j,k))/(2*dx)                                                                  &
-(vbar(i,j,k)+vtil(i,j,k))*(the1(i,j+1,k)-the1(i,j-1,k))/(2*dy)                                                                  &
+heat-rad(i,j,k,pre1,the1,wbar,wtil,precipitation)                                                                               &
+parup12(k)*(the1(i+1,j,k)-2*the1(i,j,k)+the1(i-1,j,k))/(dx*dx)                                                                  &
+parup12(k)*(the1(i,j+1,k)-2*the1(i,j,k)+the1(i,j-1,k))/(dy*dy)                                                                  &
+(parup3(k+1)*(the1(i,j,k+2)-the1(i,j,k))/(2*dz)-updzthe1(i,j,ubar,vbar,pre1,the1,time))/dz                                      &
-parK*(the1(i,j,k)-the1ref(i,j,k))

fubar(i,j,k)=                                                                                                                    &
-(ubar(i,j,k)+util(i,j,k))*(ubar(i+1,j,k)-ubar(i-1,j,k))/(2*dx)                                                                  &
-(vbar(i,j,k)+vtil(i,j,k))*(ubar(i,j+1,k)-ubar(i,j-1,k))/(2*dy)                                                                  &
+2*vbar(i,j,k)*omega3(j)                                                                                                         &
-((pre1(i+1,j,k)-pre1(i-1,j,k))/(2*dx))/rho0(k)                                                                                  &
+parnu12(k)*(ubar(i+1,j,k)-2*ubar(i,j,k)+ubar(i-1,j,k))/(dx*dx)                                                                  &
+parnu12(k)*(ubar(i,j+1,k)-2*ubar(i,j,k)+ubar(i,j-1,k))/(dy*dy)                                                                  &
+(parnu3(k+1)*(ubar(i,j,k+2)-ubar(i,j,k))/(2*dz)-nudzubar(i,j,ubar,vbar))/dz

fvbar(i,j,k)=                                                                                                                    &
-(ubar(i,j,k)+util(i,j,k))*(vbar(i+1,j,k)-vbar(i-1,j,k))/(2*dx)                                                                  &
-(vbar(i,j,k)+vtil(i,j,k))*(vbar(i,j+1,k)-vbar(i,j-1,k))/(2*dy)                                                                  &
-2*ubar(i,j,k)*omega3(j)                                                                                                         &
-((pre1(i,j+1,k)-pre1(i,j-1,k))/(2*dy))/rho0(k)                                                                                  &
+parnu12(k)*(vbar(i+1,j,k)-2*vbar(i,j,k)+vbar(i-1,j,k))/(dx*dx)                                                                  &
+parnu12(k)*(vbar(i,j+1,k)-2*vbar(i,j,k)+vbar(i,j-1,k))/(dy*dy)                                                                  &
+(parnu3(k+1)*(vbar(i,j,k+2)-vbar(i,j,k))/(2*dz)-nudzvbar(i,j,ubar,vbar))/dz

fwbar(i,j,k)=0.0

frhoutil(i,j,k)=                                                                                                                 &
-((vbar(i,j+1,k)*(rho0(k)+rho1(i,j+1,k,pre1,qbar,the1))*util(i,j+1,k)-vtil(i,j+1,k)*(rho0(k)+rho1(i,j+1,k,pre1,qbar,the1))       &
*ubar(i,j+1,k))                                                                                                                  &
 -(vbar(i,j-1,k)*(rho0(k)+rho1(i,j-1,k,pre1,qbar,the1))*util(i,j-1,k)-vtil(i,j-1,k)*(rho0(k)+rho1(i,j-1,k,pre1,qbar,the1))       &
*ubar(i,j-1,k)))/(2*dy)                                                                                                          &
-(wbar(i,j,k+1)*(rho0(k+1)+rho1(i,j,k+1,pre1,qbar,the1))*util(i,j,k+1)-wtil(i,j,k+1)*(rho0(k+1)+rho1(i,j,k+1,pre1,qbar,the1))    &
*ubar(i,j,k+1))/dz                                                                                                               &
+2*(rho0(k)+rho1(i,j,k,pre1,qbar,the1))*vtil(i,j,k)*omega3(j)                                                                    &
+parxi12(k)*(ubar(i+1,j,k)-2*ubar(i,j,k)+ubar(i-1,j,k))/(dx*dx)                                                     &
+parxi12(k)*(ubar(i,j+1,k)-2*ubar(i,j,k)+ubar(i,j-1,k))/(dy*dy)                                                     &
+((parxi3(k+1)-parmu3(k+1))*(ubar(i,j,k+2)-ubar(i,j,k))/(2*dz)-ximudzubar(i,j,ubar,vbar))/dz                                     &
+parxi12(k)*(util(i+1,j,k)-2*util(i,j,k)+util(i-1,j,k))/(dx*dx)                                                                  &
+parxi12(k)*(util(i,j+1,k)-2*util(i,j,k)+util(i,j-1,k))/(dy*dy)                                                                  &
+(parxi3(k+1)*(util(i,j,k+2)-util(i,j,k))/(2*dz)-xidzutil(i,j,ubar,vbar,util))/dz                                                &
-(pre2(i+1,j,k)-pre2(i-1,j,k))/(2*dx)

frhovtil(i,j,k)=                                                                                                                 &
-((ubar(i+1,j,k)*(rho0(k)+rho1(i+1,j,k,pre1,qbar,the1))*vtil(i+1,j,k)-util(i+1,j,k)*(rho0(k)+rho1(i+1,j,k,pre1,qbar,the1))       &
*vbar(i+1,j,k))                                                                                                                  &
 -(ubar(i-1,j,k)*(rho0(k)+rho1(i-1,j,k,pre1,qbar,the1))*vtil(i-1,j,k)-util(i-1,j,k)*(rho0(k)+rho1(i-1,j,k,pre1,qbar,the1))       &
*vbar(i-1,j,k)))/(2*dx)                                                                                                          &
-(wbar(i,j,k+1)*(rho0(k+1)+rho1(i,j,k+1,pre1,qbar,the1))*vtil(i,j,k+1)-wtil(i,j,k+1)*(rho0(k+1)+rho1(i,j,k+1,pre1,qbar,the1))    &
*vbar(i,j,k+1))/dz                                                                                                               &
-2*(rho0(k)+rho1(i,j,k,pre1,qbar,the1))*util(i,j,k)*omega3(j)                                                                    &
+parxi12(k)*(vbar(i+1,j,k)-2*vbar(i,j,k)+vbar(i-1,j,k))/(dx*dx)                                                     &
+parxi12(k)*(vbar(i,j+1,k)-2*vbar(i,j,k)+vbar(i,j-1,k))/(dy*dy)                                                     &
+((parxi3(k+1)-parmu3(k+1))*(vbar(i,j,k+2)-vbar(i,j,k))/(2*dz)-ximudzvbar(i,j,ubar,vbar))/dz                                     &
+parxi12(k)*(vtil(i+1,j,k)-2*vtil(i,j,k)+vtil(i-1,j,k))/(dx*dx)                                                                  &
+parxi12(k)*(vtil(i,j+1,k)-2*vtil(i,j,k)+vtil(i,j-1,k))/(dy*dy)                                                                  &
+(parxi3(k+1)*(vtil(i,j,k+2)-vtil(i,j,k))/(2*dz)-xidzvtil(i,j,ubar,vbar,vtil))/dz                                                &
-(pre2(i,j+1,k)-pre2(i,j-1,k))/(2*dy)

frhowtil(i,j,k)=0.0
        END DO
END DO

k=nz
DO j=1,ny-1
        DO i=1,nx-1
precipitation=precip(i,j,k,pre1,wbar,the1,qbar,wtil)
heat=((parc/(pre0(k)+pre1(i,j,k)))**kappa)*precipitation*parLv/parcp

fqbar(i,j,k)=                                                                                                                    &
-(ubar(i,j,k)+util(i,j,k))*(qbar(i+1,j,k)-qbar(i-1,j,k))/(2*dx)                                                                  &
-(vbar(i,j,k)+vtil(i,j,k))*(qbar(i,j+1,k)-qbar(i,j-1,k))/(2*dy)                                                                  &
-precipitation                                                                                                                   &
+parva12(k)*(qbar(i+1,j,k)-2*qbar(i,j,k)+qbar(i-1,j,k))/(dx*dx)                                                                  &
+parva12(k)*(qbar(i,j+1,k)-2*qbar(i,j,k)+qbar(i,j-1,k))/(dy*dy)                                                                  &
+parva3(k)*(qbar(i,j,k)-2*qbar(i,j,k-1)+qbar(i,j,k-2))/(dz*dz)

fthe1(i,j,k)=0.0

fubar(i,j,k)=                                                                                                                    &
-(ubar(i,j,k)+util(i,j,k))*(ubar(i+1,j,k)-ubar(i-1,j,k))/(2*dx)                                                                  &
-(vbar(i,j,k)+vtil(i,j,k))*(ubar(i,j+1,k)-ubar(i,j-1,k))/(2*dy)                                                                  &
+2*vbar(i,j,k)*omega3(j)                                                                                                         &
-((pre1(i+1,j,k)-pre1(i-1,j,k))/(2*dx))/rho0(k)                                                                                  &
+parnu12(k)*(ubar(i+1,j,k)-2*ubar(i,j,k)+ubar(i-1,j,k))/(dx*dx)                                                                  &
+parnu12(k)*(ubar(i,j+1,k)-2*ubar(i,j,k)+ubar(i,j-1,k))/(dy*dy)                                                                  &
+parnu3(k)*(ubar(i,j,k)-2*ubar(i,j,k-1)+ubar(i,j,k-2))/(dz*dz)

fvbar(i,j,k)=                                                                                                                    &
-(ubar(i,j,k)+util(i,j,k))*(vbar(i+1,j,k)-vbar(i-1,j,k))/(2*dx)                                                                  &
-(vbar(i,j,k)+vtil(i,j,k))*(vbar(i,j+1,k)-vbar(i,j-1,k))/(2*dy)                                                                  &
-2*ubar(i,j,k)*omega3(j)                                                                                                         &
-((pre1(i,j+1,k)-pre1(i,j-1,k))/(2*dy))/rho0(k)                                                                                  &
+parnu12(k)*(vbar(i+1,j,k)-2*vbar(i,j,k)+vbar(i-1,j,k))/(dx*dx)                                                                  &
+parnu12(k)*(vbar(i,j+1,k)-2*vbar(i,j,k)+vbar(i,j-1,k))/(dy*dy)                                                                  &
+parnu3(k)*(vbar(i,j,k)-2*vbar(i,j,k-1)+vbar(i,j,k-2))/(dz*dz)

fwbar(i,j,k)=0.0

frhoutil(i,j,k)=                                                                                                                 &
-((vbar(i,j+1,k)*(rho0(k)+rho1(i,j+1,k,pre1,qbar,the1))*util(i,j+1,k)-vtil(i,j+1,k)*(rho0(k)+rho1(i,j+1,k,pre1,qbar,the1))       &
*ubar(i,j+1,k))                                                                                                                  &
 -(vbar(i,j-1,k)*(rho0(k)+rho1(i,j-1,k,pre1,qbar,the1))*util(i,j-1,k)-vtil(i,j-1,k)*(rho0(k)+rho1(i,j-1,k,pre1,qbar,the1))       &
*ubar(i,j-1,k)))/(2*dy)                                                                                                          &
+(wbar(i,j,k-1)*(rho0(k-1)+rho1(i,j,k-1,pre1,qbar,the1))*util(i,j,k-1)-wtil(i,j,k-1)*(rho0(k-1)+rho1(i,j,k-1,pre1,qbar,the1))    &
*ubar(i,j,k-1))/dz                                                                                                               &
+2*(rho0(k)+rho1(i,j,k,pre1,qbar,the1))*vtil(i,j,k)*omega3(j)                                                                    &
+parxi12(k)*(ubar(i+1,j,k)-2*ubar(i,j,k)+ubar(i-1,j,k))/(dx*dx)                                                     &
+parxi12(k)*(ubar(i,j+1,k)-2*ubar(i,j,k)+ubar(i,j-1,k))/(dy*dy)                                                     &
+parxi3(k)*(ubar(i,j,k)-2*ubar(i,j,k-1)+ubar(i,j,k-2))/(dz*dz)                                                       &
+parxi12(k)*(util(i+1,j,k)-2*util(i,j,k)+util(i-1,j,k))/(dx*dx)                                                                  &
+parxi12(k)*(util(i,j+1,k)-2*util(i,j,k)+util(i,j-1,k))/(dy*dy)                                                                  &
+parxi3(k)*(util(i,j,k)-2*util(i,j,k-1)+util(i,j,k-2))/(dz*dz)                                                                   &
-(pre2(i+1,j,k)-pre2(i-1,j,k))/(2*dx)

frhovtil(i,j,k)=                                                                                                                 &
-((ubar(i+1,j,k)*(rho0(k)+rho1(i+1,j,k,pre1,qbar,the1))*vtil(i+1,j,k)-util(i+1,j,k)*(rho0(k)+rho1(i+1,j,k,pre1,qbar,the1))       &
*vbar(i+1,j,k))                                                                                                                  &
 -(ubar(i-1,j,k)*(rho0(k)+rho1(i-1,j,k,pre1,qbar,the1))*vtil(i-1,j,k)-util(i-1,j,k)*(rho0(k)+rho1(i-1,j,k,pre1,qbar,the1))       &
*vbar(i-1,j,k)))/(2*dx)                                                                                                          &
+(wbar(i,j,k-1)*(rho0(k-1)+rho1(i,j,k-1,pre1,qbar,the1))*vtil(i,j,k-1)-wtil(i,j,k-1)*(rho0(k-1)+rho1(i,j,k-1,pre1,qbar,the1))    &
*vbar(i,j,k-1))/dz                                                                                                               &
-2*(rho0(k)+rho1(i,j,k,pre1,qbar,the1))*util(i,j,k)*omega3(j)                                                                    &
+parxi12(k)*(vbar(i+1,j,k)-2*vbar(i,j,k)+vbar(i-1,j,k))/(dx*dx)                                                     &
+parxi12(k)*(vbar(i,j+1,k)-2*vbar(i,j,k)+vbar(i,j-1,k))/(dy*dy)                                                     &
+parxi3(k)*(vbar(i,j,k)-2*vbar(i,j,k-1)+vbar(i,j,k-2))/(dz*dz)                                                       &
+parxi12(k)*(vtil(i+1,j,k)-2*vtil(i,j,k)+vtil(i-1,j,k))/(dx*dx)                                                                  &
+parxi12(k)*(vtil(i,j+1,k)-2*vtil(i,j,k)+vtil(i,j-1,k))/(dy*dy)                                                                  &
+parxi3(k)*(vtil(i,j,k)-2*vtil(i,j,k-1)+vtil(i,j,k-2))/(dz*dz)                                                                   &
-(pre2(i,j+1,k)-pre2(i,j-1,k))/(2*dy)

frhowtil(i,j,k)=0.0
        END DO
END DO
!!! core: END

END SUBROUTINE set_field1

SUBROUTINE set_field2(delta,time,pre1,qbar,the1,ubar,vbar,wbar,util,vtil,wtil,fubar,fvbar,fwbar,frhoutil,frhovtil,frhowtil)

REAL,INTENT(IN) :: delta
REAL,INTENT(IN) :: time

REAL,DIMENSION(0:nx,0:ny,nz),INTENT(IN) :: pre1,qbar,the1
REAL,DIMENSION(0:nx,0:ny,nz),INTENT(IN) :: ubar,vbar,wbar
REAL,DIMENSION(0:nx,0:ny,nz),INTENT(IN) :: util,vtil,wtil

REAL,DIMENSION(0:nx,0:ny,nz),INTENT(OUT) :: fubar,fvbar,fwbar
REAL,DIMENSION(0:nx,0:ny,nz),INTENT(OUT) :: frhoutil,frhovtil,frhowtil

! LOCAL VARIABLES:

INTEGER :: i,j,k

!!! core: BEGIN
DO k=2,nz-1
        DO j=1,ny-1
                DO i=1,nx-1
fubar(i,j,k)=                                                                                                                    &
-(ubar(i,j,k)+util(i,j,k))*(ubar(i+1,j,k)-ubar(i-1,j,k))/(2*dx)                                                                  &
-(vbar(i,j,k)+vtil(i,j,k))*(ubar(i,j+1,k)-ubar(i,j-1,k))/(2*dy)                                                                  &
-(wbar(i,j,k)+wtil(i,j,k))*(ubar(i,j,k+1)-ubar(i,j,k-1))/(2*dz)                                                                  &
+2*(vbar(i,j,k)*omega3(j)-wbar(i,j,k)*omega2(j))                                                                                 &
+parnu12(k)*(ubar(i+1,j,k)-2*ubar(i,j,k)+ubar(i-1,j,k))/(dx*dx)                                                                  &
+parnu12(k)*(ubar(i,j+1,k)-2*ubar(i,j,k)+ubar(i,j-1,k))/(dy*dy)                                                                  &
+parnu3(k)*(ubar(i,j,k+1)-2*ubar(i,j,k)+ubar(i,j,k-1))/(dz*dz)

fvbar(i,j,k)=                                                                                                                    &
-(ubar(i,j,k)+util(i,j,k))*(vbar(i+1,j,k)-vbar(i-1,j,k))/(2*dx)                                                                  &
-(vbar(i,j,k)+vtil(i,j,k))*(vbar(i,j+1,k)-vbar(i,j-1,k))/(2*dy)                                                                  &
-(wbar(i,j,k)+wtil(i,j,k))*(vbar(i,j,k+1)-vbar(i,j,k-1))/(2*dz)                                                                  &
-2*ubar(i,j,k)*omega3(j)                                                                                                         &
+parnu12(k)*(vbar(i+1,j,k)-2*vbar(i,j,k)+vbar(i-1,j,k))/(dx*dx)                                                                  &
+parnu12(k)*(vbar(i,j+1,k)-2*vbar(i,j,k)+vbar(i,j-1,k))/(dy*dy)                                                                  &
+parnu3(k)*(vbar(i,j,k+1)-2*vbar(i,j,k)+vbar(i,j,k-1))/(dz*dz)

fwbar(i,j,k)=                                                                                                                    &
-(ubar(i,j,k)+util(i,j,k))*(wbar(i+1,j,k)-wbar(i-1,j,k))/(2*dx)                                                                  &
-(vbar(i,j,k)+vtil(i,j,k))*(wbar(i,j+1,k)-wbar(i,j-1,k))/(2*dy)                                                                  &
-(wbar(i,j,k)+wtil(i,j,k))*(wbar(i,j,k+1)-wbar(i,j,k-1))/(2*dz)                                                                  &
+2*ubar(i,j,k)*omega2(j)                                                                                                         &
+parg*the1(i,j,k)/the0(k)                                                                                                        &
+parnu12(k)*(wbar(i+1,j,k)-2*wbar(i,j,k)+wbar(i-1,j,k))/(dx*dx)                                                                  &
+parnu12(k)*(wbar(i,j+1,k)-2*wbar(i,j,k)+wbar(i,j-1,k))/(dy*dy)                                                                  &
+parnu3(k)*(wbar(i,j,k+1)-2*wbar(i,j,k)+wbar(i,j,k-1))/(dz*dz)

frhoutil(i,j,k)=                                                                                                                 &
-((vbar(i,j+1,k)*(rho0(k)+rho1(i,j+1,k,pre1,qbar,the1))*util(i,j+1,k)-vtil(i,j+1,k)*(rho0(k)+rho1(i,j+1,k,pre1,qbar,the1))       &
*ubar(i,j+1,k))                                                                                                                  &
 -(vbar(i,j-1,k)*(rho0(k)+rho1(i,j-1,k,pre1,qbar,the1))*util(i,j-1,k)-vtil(i,j-1,k)*(rho0(k)+rho1(i,j-1,k,pre1,qbar,the1))       &
*ubar(i,j-1,k)))/(2*dy)                                                                                                          &
-((wbar(i,j,k+1)*(rho0(k+1)+rho1(i,j,k+1,pre1,qbar,the1))*util(i,j,k+1)-wtil(i,j,k+1)*(rho0(k+1)+rho1(i,j,k+1,pre1,qbar,the1))   &
*ubar(i,j,k+1))                                                                                                                  &
 -(wbar(i,j,k-1)*(rho0(k-1)+rho1(i,j,k-1,pre1,qbar,the1))*util(i,j,k-1)-wtil(i,j,k-1)*(rho0(k-1)+rho1(i,j,k-1,pre1,qbar,the1))   &
*ubar(i,j,k-1)))/(2*dz)                                                                                                          &
+2*(rho0(k)+rho1(i,j,k,pre1,qbar,the1))*(vtil(i,j,k)*omega3(j)-wtil(i,j,k)*omega2(j))                                            &
+parxi12(k)*(ubar(i+1,j,k)-2*ubar(i,j,k)+ubar(i-1,j,k))/(dx*dx)                                                     &
+parxi12(k)*(ubar(i,j+1,k)-2*ubar(i,j,k)+ubar(i,j-1,k))/(dy*dy)                                                     &
+parxi3(k)*(ubar(i,j,k+1)-2*ubar(i,j,k)+ubar(i,j,k-1))/(dz*dz)                                                       &
+parxi12(k)*(util(i+1,j,k)-2*util(i,j,k)+util(i-1,j,k))/(dx*dx)                                                                  &
+parxi12(k)*(util(i,j+1,k)-2*util(i,j,k)+util(i,j-1,k))/(dy*dy)                                                                  &
+parxi3(k)*(util(i,j,k+1)-2*util(i,j,k)+util(i,j,k-1))/(dz*dz)

frhovtil(i,j,k)=                                                                                                                 &
-((ubar(i+1,j,k)*(rho0(k)+rho1(i+1,j,k,pre1,qbar,the1))*vtil(i+1,j,k)-util(i+1,j,k)*(rho0(k)+rho1(i+1,j,k,pre1,qbar,the1))       &
*vbar(i+1,j,k))                                                                                                                  &
 -(ubar(i-1,j,k)*(rho0(k)+rho1(i-1,j,k,pre1,qbar,the1))*vtil(i-1,j,k)-util(i-1,j,k)*(rho0(k)+rho1(i-1,j,k,pre1,qbar,the1))       &
*vbar(i-1,j,k)))/(2*dx)                                                                                                          &
-((wbar(i,j,k+1)*(rho0(k+1)+rho1(i,j,k+1,pre1,qbar,the1))*vtil(i,j,k+1)-wtil(i,j,k+1)*(rho0(k+1)+rho1(i,j,k+1,pre1,qbar,the1))   &
*vbar(i,j,k+1))                                                                                                                  &
 -(wbar(i,j,k-1)*(rho0(k-1)+rho1(i,j,k-1,pre1,qbar,the1))*vtil(i,j,k-1)-wtil(i,j,k-1)*(rho0(k-1)+rho1(i,j,k-1,pre1,qbar,the1))   &
*vbar(i,j,k-1)))/(2*dz)                                                                                                          &
-2*(rho0(k)+rho1(i,j,k,pre1,qbar,the1))*util(i,j,k)*omega3(j)                                                                    &
+parxi12(k)*(vbar(i+1,j,k)-2*vbar(i,j,k)+vbar(i-1,j,k))/(dx*dx)                                                     &
+parxi12(k)*(vbar(i,j+1,k)-2*vbar(i,j,k)+vbar(i,j-1,k))/(dy*dy)                                                     &
+parxi3(k)*(vbar(i,j,k+1)-2*vbar(i,j,k)+vbar(i,j,k-1))/(dz*dz)                                                       &
+parxi12(k)*(vtil(i+1,j,k)-2*vtil(i,j,k)+vtil(i-1,j,k))/(dx*dx)                                                                  &
+parxi12(k)*(vtil(i,j+1,k)-2*vtil(i,j,k)+vtil(i,j-1,k))/(dy*dy)                                                                  &
+parxi3(k)*(vtil(i,j,k+1)-2*vtil(i,j,k)+vtil(i,j,k-1))/(dz*dz)

frhowtil(i,j,k)=                                                                                                                 &
-((ubar(i+1,j,k)*(rho0(k)+rho1(i+1,j,k,pre1,qbar,the1))*wtil(i+1,j,k)-util(i+1,j,k)*(rho0(k)+rho1(i+1,j,k,pre1,qbar,the1))       &
*wbar(i+1,j,k))                                                                                                                  &
 -(ubar(i-1,j,k)*(rho0(k)+rho1(i-1,j,k,pre1,qbar,the1))*wtil(i-1,j,k)-util(i-1,j,k)*(rho0(k)+rho1(i-1,j,k,pre1,qbar,the1))       &
*wbar(i-1,j,k)))/(2*dx)                                                                                                          &
-((vbar(i,j+1,k)*(rho0(k)+rho1(i,j+1,k,pre1,qbar,the1))*wtil(i,j+1,k)-vtil(i,j+1,k)*(rho0(k)+rho1(i,j+1,k,pre1,qbar,the1))       &
*wbar(i,j+1,k))                                                                                                                  &
 -(vbar(i,j-1,k)*(rho0(k)+rho1(i,j-1,k,pre1,qbar,the1))*wtil(i,j-1,k)-vtil(i,j-1,k)*(rho0(k)+rho1(i,j-1,k,pre1,qbar,the1))       &
*wbar(i,j-1,k)))/(2*dy)                                                                                                          &
+2*(rho0(k)+rho1(i,j,k,pre1,qbar,the1))*util(i,j,k)*omega2(j)                                                                    &
+delta*(rho0(k)+rho1(i,j,k,pre1,qbar,the1))*parg*(rho1(i,j,k,pre1,qbar,the1)/rho0(k))**2                                         &
+parxi12(k)*(wbar(i+1,j,k)-2*wbar(i,j,k)+wbar(i-1,j,k))/(dx*dx)                                                     &
+parxi12(k)*(wbar(i,j+1,k)-2*wbar(i,j,k)+wbar(i,j-1,k))/(dy*dy)                                                     &
+parxi3(k)*(wbar(i,j,k+1)-2*wbar(i,j,k)+wbar(i,j,k-1))/(dz*dz)                                                       &
+parxi12(k)*(wtil(i+1,j,k)-2*wtil(i,j,k)+wtil(i-1,j,k))/(dx*dx)                                                                  &
+parxi12(k)*(wtil(i,j+1,k)-2*wtil(i,j,k)+wtil(i,j-1,k))/(dy*dy)                                                                  &
+parxi3(k)*(wtil(i,j,k+1)-2*wtil(i,j,k)+wtil(i,j,k-1))/(dz*dz)
                END DO
        END DO
END DO

k=1
DO j=1,ny-1
        DO i=1,nx-1
fubar(i,j,k)=                                                                                                                    &
-(ubar(i,j,k)+util(i,j,k))*(ubar(i+1,j,k)-ubar(i-1,j,k))/(2*dx)                                                                  &
-(vbar(i,j,k)+vtil(i,j,k))*(ubar(i,j+1,k)-ubar(i,j-1,k))/(2*dy)                                                                  &
+2*vbar(i,j,k)*omega3(j)                                                                                                         &
+parnu12(k)*(ubar(i+1,j,k)-2*ubar(i,j,k)+ubar(i-1,j,k))/(dx*dx)                                                                  &
+parnu12(k)*(ubar(i,j+1,k)-2*ubar(i,j,k)+ubar(i,j-1,k))/(dy*dy)                                                                  &
+(parnu3(k+1)*(ubar(i,j,k+2)-ubar(i,j,k))/(2*dz)-nudzubar(i,j,ubar,vbar))/dz

fvbar(i,j,k)=                                                                                                                    &
-(ubar(i,j,k)+util(i,j,k))*(vbar(i+1,j,k)-vbar(i-1,j,k))/(2*dx)                                                                  &
-(vbar(i,j,k)+vtil(i,j,k))*(vbar(i,j+1,k)-vbar(i,j-1,k))/(2*dy)                                                                  &
-2*ubar(i,j,k)*omega3(j)                                                                                                         &
+parnu12(k)*(vbar(i+1,j,k)-2*vbar(i,j,k)+vbar(i-1,j,k))/(dx*dx)                                                                  &
+parnu12(k)*(vbar(i,j+1,k)-2*vbar(i,j,k)+vbar(i,j-1,k))/(dy*dy)                                                                  &
+(parnu3(k+1)*(vbar(i,j,k+2)-vbar(i,j,k))/(2*dz)-nudzvbar(i,j,ubar,vbar))/dz

fwbar(i,j,k)=0.0

frhoutil(i,j,k)=                                                                                                                 &
-((vbar(i,j+1,k)*(rho0(k)+rho1(i,j+1,k,pre1,qbar,the1))*util(i,j+1,k)-vtil(i,j+1,k)*(rho0(k)+rho1(i,j+1,k,pre1,qbar,the1))       &
*ubar(i,j+1,k))                                                                                                                  &
 -(vbar(i,j-1,k)*(rho0(k)+rho1(i,j-1,k,pre1,qbar,the1))*util(i,j-1,k)-vtil(i,j-1,k)*(rho0(k)+rho1(i,j-1,k,pre1,qbar,the1))       &
*ubar(i,j-1,k)))/(2*dy)                                                                                                          &
-(wbar(i,j,k+1)*(rho0(k+1)+rho1(i,j,k+1,pre1,qbar,the1))*util(i,j,k+1)-wtil(i,j,k+1)*(rho0(k+1)+rho1(i,j,k+1,pre1,qbar,the1))    &
*ubar(i,j,k+1))/dz                                                                                                               &
+2*(rho0(k)+rho1(i,j,k,pre1,qbar,the1))*vtil(i,j,k)*omega3(j)                                                                    &
+parxi12(k)*(ubar(i+1,j,k)-2*ubar(i,j,k)+ubar(i-1,j,k))/(dx*dx)                                                     &
+parxi12(k)*(ubar(i,j+1,k)-2*ubar(i,j,k)+ubar(i,j-1,k))/(dy*dy)                                                     &
+((parxi3(k+1)-parmu3(k+1))*(ubar(i,j,k+2)-ubar(i,j,k))/(2*dz)-ximudzubar(i,j,ubar,vbar))/dz                                     &
+parxi12(k)*(util(i+1,j,k)-2*util(i,j,k)+util(i-1,j,k))/(dx*dx)                                                                  &
+parxi12(k)*(util(i,j+1,k)-2*util(i,j,k)+util(i,j-1,k))/(dy*dy)                                                                  &
+(parxi3(k+1)*(util(i,j,k+2)-util(i,j,k))/(2*dz)-xidzutil(i,j,ubar,vbar,util))/dz

frhovtil(i,j,k)=                                                                                                                 &
-((ubar(i+1,j,k)*(rho0(k)+rho1(i+1,j,k,pre1,qbar,the1))*vtil(i+1,j,k)-util(i+1,j,k)*(rho0(k)+rho1(i+1,j,k,pre1,qbar,the1))       &
*vbar(i+1,j,k))                                                                                                                  &
 -(ubar(i-1,j,k)*(rho0(k)+rho1(i-1,j,k,pre1,qbar,the1))*vtil(i-1,j,k)-util(i-1,j,k)*(rho0(k)+rho1(i-1,j,k,pre1,qbar,the1))       &
*vbar(i-1,j,k)))/(2*dx)                                                                                                          &
-(wbar(i,j,k+1)*(rho0(k+1)+rho1(i,j,k+1,pre1,qbar,the1))*vtil(i,j,k+1)-wtil(i,j,k+1)*(rho0(k+1)+rho1(i,j,k+1,pre1,qbar,the1))    &
*vbar(i,j,k+1))/dz                                                                                                               &
-2*(rho0(k)+rho1(i,j,k,pre1,qbar,the1))*util(i,j,k)*omega3(j)                                                                    &
+parxi12(k)*(vbar(i+1,j,k)-2*vbar(i,j,k)+vbar(i-1,j,k))/(dx*dx)                                                     &
+parxi12(k)*(vbar(i,j+1,k)-2*vbar(i,j,k)+vbar(i,j-1,k))/(dy*dy)                                                     &
+((parxi3(k+1)-parmu3(k+1))*(vbar(i,j,k+2)-vbar(i,j,k))/(2*dz)-ximudzvbar(i,j,ubar,vbar))/dz                                     &
+parxi12(k)*(vtil(i+1,j,k)-2*vtil(i,j,k)+vtil(i-1,j,k))/(dx*dx)                                                                  &
+parxi12(k)*(vtil(i,j+1,k)-2*vtil(i,j,k)+vtil(i,j-1,k))/(dy*dy)                                                                  &
+(parxi3(k+1)*(vtil(i,j,k+2)-vtil(i,j,k))/(2*dz)-xidzvtil(i,j,ubar,vbar,vtil))/dz

frhowtil(i,j,k)=0.0
        END DO
END DO

k=nz
DO j=1,ny-1
        DO i=1,nx-1
fubar(i,j,k)=                                                                                                                    &
-(ubar(i,j,k)+util(i,j,k))*(ubar(i+1,j,k)-ubar(i-1,j,k))/(2*dx)                                                                  &
-(vbar(i,j,k)+vtil(i,j,k))*(ubar(i,j+1,k)-ubar(i,j-1,k))/(2*dy)                                                                  &
+2*vbar(i,j,k)*omega3(j)                                                                                                         &
+parnu12(k)*(ubar(i+1,j,k)-2*ubar(i,j,k)+ubar(i-1,j,k))/(dx*dx)                                                                  &
+parnu12(k)*(ubar(i,j+1,k)-2*ubar(i,j,k)+ubar(i,j-1,k))/(dy*dy)                                                                  &
+parnu3(k)*(ubar(i,j,k)-2*ubar(i,j,k-1)+ubar(i,j,k-2))/(dz*dz)

fvbar(i,j,k)=                                                                                                                    &
-(ubar(i,j,k)+util(i,j,k))*(vbar(i+1,j,k)-vbar(i-1,j,k))/(2*dx)                                                                  &
-(vbar(i,j,k)+vtil(i,j,k))*(vbar(i,j+1,k)-vbar(i,j-1,k))/(2*dy)                                                                  &
-2*ubar(i,j,k)*omega3(j)                                                                                                         &
+parnu12(k)*(vbar(i+1,j,k)-2*vbar(i,j,k)+vbar(i-1,j,k))/(dx*dx)                                                                  &
+parnu12(k)*(vbar(i,j+1,k)-2*vbar(i,j,k)+vbar(i,j-1,k))/(dy*dy)                                                                  &
+parnu3(k)*(vbar(i,j,k)-2*vbar(i,j,k-1)+vbar(i,j,k-2))/(dz*dz)

fwbar(i,j,k)=0.0

frhoutil(i,j,k)=                                                                                                                 &
-((vbar(i,j+1,k)*(rho0(k)+rho1(i,j+1,k,pre1,qbar,the1))*util(i,j+1,k)-vtil(i,j+1,k)*(rho0(k)+rho1(i,j+1,k,pre1,qbar,the1))       &
*ubar(i,j+1,k))                                                                                                                  &
 -(vbar(i,j-1,k)*(rho0(k)+rho1(i,j-1,k,pre1,qbar,the1))*util(i,j-1,k)-vtil(i,j-1,k)*(rho0(k)+rho1(i,j-1,k,pre1,qbar,the1))       &
*ubar(i,j-1,k)))/(2*dy)                                                                                                          &
+(wbar(i,j,k-1)*(rho0(k-1)+rho1(i,j,k-1,pre1,qbar,the1))*util(i,j,k-1)-wtil(i,j,k-1)*(rho0(k-1)+rho1(i,j,k-1,pre1,qbar,the1))    &
*ubar(i,j,k-1))/dz                                                                                                               &
+2*(rho0(k)+rho1(i,j,k,pre1,qbar,the1))*vtil(i,j,k)*omega3(j)                                                                    &
+parxi12(k)*(ubar(i+1,j,k)-2*ubar(i,j,k)+ubar(i-1,j,k))/(dx*dx)                                                     &
+parxi12(k)*(ubar(i,j+1,k)-2*ubar(i,j,k)+ubar(i,j-1,k))/(dy*dy)                                                     &
+parxi3(k)*(ubar(i,j,k)-2*ubar(i,j,k-1)+ubar(i,j,k-2))/(dz*dz)                                                       &
+parxi12(k)*(util(i+1,j,k)-2*util(i,j,k)+util(i-1,j,k))/(dx*dx)                                                                  &
+parxi12(k)*(util(i,j+1,k)-2*util(i,j,k)+util(i,j-1,k))/(dy*dy)                                                                  &
+parxi3(k)*(util(i,j,k)-2*util(i,j,k-1)+util(i,j,k-2))/(dz*dz)

frhovtil(i,j,k)=                                                                                                                 &
-((ubar(i+1,j,k)*(rho0(k)+rho1(i+1,j,k,pre1,qbar,the1))*vtil(i+1,j,k)-util(i+1,j,k)*(rho0(k)+rho1(i+1,j,k,pre1,qbar,the1))       &
*vbar(i+1,j,k))                                                                                                                  &
 -(ubar(i-1,j,k)*(rho0(k)+rho1(i-1,j,k,pre1,qbar,the1))*vtil(i-1,j,k)-util(i-1,j,k)*(rho0(k)+rho1(i-1,j,k,pre1,qbar,the1))       &
*vbar(i-1,j,k)))/(2*dx)                                                                                                          &
+(wbar(i,j,k-1)*(rho0(k-1)+rho1(i,j,k-1,pre1,qbar,the1))*vtil(i,j,k-1)-wtil(i,j,k-1)*(rho0(k-1)+rho1(i,j,k-1,pre1,qbar,the1))    &
*vbar(i,j,k-1))/dz                                                                                                               &
-2*(rho0(k)+rho1(i,j,k,pre1,qbar,the1))*util(i,j,k)*omega3(j)                                                                    &
+parxi12(k)*(vbar(i+1,j,k)-2*vbar(i,j,k)+vbar(i-1,j,k))/(dx*dx)                                                     &
+parxi12(k)*(vbar(i,j+1,k)-2*vbar(i,j,k)+vbar(i,j-1,k))/(dy*dy)                                                     &
+parxi3(k)*(vbar(i,j,k)-2*vbar(i,j,k-1)+vbar(i,j,k-2))/(dz*dz)                                                       &
+parxi12(k)*(vtil(i+1,j,k)-2*vtil(i,j,k)+vtil(i-1,j,k))/(dx*dx)                                                                  &
+parxi12(k)*(vtil(i,j+1,k)-2*vtil(i,j,k)+vtil(i,j-1,k))/(dy*dy)                                                                  &
+parxi3(k)*(vtil(i,j,k)-2*vtil(i,j,k-1)+vtil(i,j,k-2))/(dz*dz)

frhowtil(i,j,k)=0.0
        END DO
END DO
!!! core: END

END SUBROUTINE set_field2

SUBROUTINE set_field3(time,pre1,qbar,the1,ubar,vbar,wbar,util,vtil,wtil,the1ref,fqbar,fthe1)

REAL,INTENT(IN) :: time

REAL,DIMENSION(0:nx,0:ny,nz),INTENT(IN) :: pre1,qbar,the1
REAL,DIMENSION(0:nx,0:ny,nz),INTENT(IN) :: ubar,vbar,wbar
REAL,DIMENSION(0:nx,0:ny,nz),INTENT(IN) :: util,vtil,wtil
REAL,DIMENSION(0:nx,0:ny,nz),INTENT(IN) :: the1ref

REAL,DIMENSION(0:nx,0:ny,nz),INTENT(OUT) :: fqbar,fthe1

! LOCAL VARIABLES:

INTEGER :: i,j,k

REAL :: precipitation,heat

!!! core: BEGIN
DO k=2,nz-1
        DO j=1,ny-1
                DO i=1,nx-1
precipitation=precip(i,j,k,pre1,wbar,the1,qbar,wtil)
heat=((parc/(pre0(k)+pre1(i,j,k)))**kappa)*precipitation*parLv/parcp

fqbar(i,j,k)=                                                                                                                    &
-(ubar(i,j,k)+util(i,j,k))*(qbar(i+1,j,k)-qbar(i-1,j,k))/(2*dx)                                                                  &
-(vbar(i,j,k)+vtil(i,j,k))*(qbar(i,j+1,k)-qbar(i,j-1,k))/(2*dy)                                                                  &
-(wbar(i,j,k)+wtil(i,j,k))*(qbar(i,j,k+1)-qbar(i,j,k-1))/(2*dz)                                                                  &
-precipitation                                                                                                                   &
+parva12(k)*(qbar(i+1,j,k)-2*qbar(i,j,k)+qbar(i-1,j,k))/(dx*dx)                                                                  &
+parva12(k)*(qbar(i,j+1,k)-2*qbar(i,j,k)+qbar(i,j-1,k))/(dy*dy)                                                                  &
+parva3(k)*(qbar(i,j,k+1)-2*qbar(i,j,k)+qbar(i,j,k-1))/(dz*dz)

fthe1(i,j,k)=                                                                                                                    &
-(ubar(i,j,k)+util(i,j,k))*(the1(i+1,j,k)-the1(i-1,j,k))/(2*dx)                                                                  &
-(vbar(i,j,k)+vtil(i,j,k))*(the1(i,j+1,k)-the1(i,j-1,k))/(2*dy)                                                                  &
-(wbar(i,j,k)+wtil(i,j,k))*(the1(i,j,k+1)-the1(i,j,k-1))/(2*dz)                                                                  &
-(wbar(i,j,k)+wtil(i,j,k))*dzthe0(k)                                                                                             &
+heat-rad(i,j,k,pre1,the1,wbar,wtil,precipitation)                                                                               &
+parup12(k)*(the1(i+1,j,k)-2*the1(i,j,k)+the1(i-1,j,k))/(dx*dx)                                                                  &
+parup12(k)*(the1(i,j+1,k)-2*the1(i,j,k)+the1(i,j-1,k))/(dy*dy)                                                                  &
+parup3(k)*(the1(i,j,k+1)-2*the1(i,j,k)+the1(i,j,k-1))/(dz*dz)                                                                   &
-parK*(the1(i,j,k)-the1ref(i,j,k))
                END DO
        END DO
END DO

k=1
DO j=1,ny-1
        DO i=1,nx-1
precipitation=precip(i,j,k,pre1,wbar,the1,qbar,wtil)
heat=((parc/(pre0(k)+pre1(i,j,k)))**kappa)*precipitation*parLv/parcp

fqbar(i,j,k)=                                                                                                                    &
-(ubar(i,j,k)+util(i,j,k))*(qbar(i+1,j,k)-qbar(i-1,j,k))/(2*dx)                                                                  &
-(vbar(i,j,k)+vtil(i,j,k))*(qbar(i,j+1,k)-qbar(i,j-1,k))/(2*dy)                                                                  &
-precipitation                                                                                                                   &
+parva12(k)*(qbar(i+1,j,k)-2*qbar(i,j,k)+qbar(i-1,j,k))/(dx*dx)                                                                  &
+parva12(k)*(qbar(i,j+1,k)-2*qbar(i,j,k)+qbar(i,j-1,k))/(dy*dy)                                                                  &
+(parva3(k+1)*(qbar(i,j,k+2)-qbar(i,j,k))/(2*dz)-vadzqbar(i,j,ubar,vbar,qbar,time))/dz

fthe1(i,j,k)=                                                                                                                    &
-(ubar(i,j,k)+util(i,j,k))*(the1(i+1,j,k)-the1(i-1,j,k))/(2*dx)                                                                  &
-(vbar(i,j,k)+vtil(i,j,k))*(the1(i,j+1,k)-the1(i,j-1,k))/(2*dy)                                                                  &
+heat-rad(i,j,k,pre1,the1,wbar,wtil,precipitation)                                                                               &
+parup12(k)*(the1(i+1,j,k)-2*the1(i,j,k)+the1(i-1,j,k))/(dx*dx)                                                                  &
+parup12(k)*(the1(i,j+1,k)-2*the1(i,j,k)+the1(i,j-1,k))/(dy*dy)                                                                  &
+(parup3(k+1)*(the1(i,j,k+2)-the1(i,j,k))/(2*dz)-updzthe1(i,j,ubar,vbar,pre1,the1,time))/dz                                      &
-parK*(the1(i,j,k)-the1ref(i,j,k))
        END DO
END DO

k=nz
DO j=1,ny-1
        DO i=1,nx-1
precipitation=precip(i,j,k,pre1,wbar,the1,qbar,wtil)
heat=((parc/(pre0(k)+pre1(i,j,k)))**kappa)*precipitation*parLv/parcp

fqbar(i,j,k)=                                                                                                                    &
-(ubar(i,j,k)+util(i,j,k))*(qbar(i+1,j,k)-qbar(i-1,j,k))/(2*dx)                                                                  &
-(vbar(i,j,k)+vtil(i,j,k))*(qbar(i,j+1,k)-qbar(i,j-1,k))/(2*dy)                                                                  &
-precipitation                                                                                                                   &
+parva12(k)*(qbar(i+1,j,k)-2*qbar(i,j,k)+qbar(i-1,j,k))/(dx*dx)                                                                  &
+parva12(k)*(qbar(i,j+1,k)-2*qbar(i,j,k)+qbar(i,j-1,k))/(dy*dy)                                                                  &
+parva3(k)*(qbar(i,j,k)-2*qbar(i,j,k-1)+qbar(i,j,k-2))/(dz*dz)

fthe1(i,j,k)=0.0
        END DO
END DO
!!! core: END

END SUBROUTINE set_field3

END MODULE field
