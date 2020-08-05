MODULE parameters

IMPLICIT NONE
PUBLIC

! MATH CONSTANTS:

REAL,PARAMETER :: parpi=3.1416

! MODEL PARAMETERS:

REAL,PARAMETER :: z00=0.8e3 ! meters
REAL,PARAMETER :: year=24*3600.0 ! seconds
REAL,PARAMETER :: Lx=40000e3,Ly=40000e3,Lz=18e3 ! meters
REAL,PARAMETER :: omega=7.3e-5 ! rad/s
REAL,PARAMETER :: radius=6350e3 ! m
REAL,PARAMETER :: beta=2.0*omega/radius
REAL,PARAMETER :: parg=10.0 ! m/s^2
REAL,PARAMETER :: parc=1e5 ! Pa
REAL,PARAMETER :: parcv=718.0 ! J/(kg*K)
REAL,PARAMETER :: parcp=1005.0 ! J/(kg*K)
REAL,PARAMETER :: parK=0.0 !1.0/(0.05*year) ! 1/s
REAL,PARAMETER :: dryR=287 ! J/(kg*K)
REAL,PARAMETER :: kappa=0.286
REAL,PARAMETER :: rho00=1.2
REAL,PARAMETER :: pre00=1e5
REAL,PARAMETER :: parT0=pre00/(dryR*rho00)
REAL,PARAMETER :: parH=dryR*parT0/parg
REAL,PARAMETER :: parLv=2.3e6 ! J/kg
REAL,PARAMETER :: parRv=462.0 ! J/(kg*K)
REAL,PARAMETER :: refEs=2337.0 ! Pa=N/m^2
REAL,PARAMETER :: refT=293.0 ! K
REAL,PARAMETER :: Cd_uvbar=0.03
REAL,PARAMETER :: Cd_uvtil=0.03
REAL,PARAMETER :: C_the1=0.03
REAL,PARAMETER :: C_qbar=0.03
REAL,PARAMETER :: nu12_0=1e1
REAL,PARAMETER :: nu3_0=1e1
REAL,PARAMETER :: mu12_0=1e1 ! Kg/(m.s)
REAL,PARAMETER :: mu3_0=1e1 ! Kg/(m.s)
REAL,PARAMETER :: xi12_0=1e2 ! MIN(dx**2,dy**2)/(5*24*3600)
REAL,PARAMETER :: xi3_0=1e2 ! (dz**2)/(5*24*3600)
REAL,PARAMETER :: up12_0=1e1 ! m^2/s
REAL,PARAMETER :: up3_0=1e1 ! m^2/s
REAL,PARAMETER :: va12_0=1e1 ! m^2/s
REAL,PARAMETER :: va3_0=1e1 ! m^2/s
REAL,PARAMETER :: flush_time=1e50 !10*60.0 ! 1/12sc
REAL,PARAMETER :: sigma=5.7e-8 ! W/(m^2*K^4)
REAL,PARAMETER :: max_wbar=1.0 ! m/s
REAL,PARAMETER :: max_wtil=1.0 ! m/s

! NUMERICAL PARAMETERS:

INTEGER,PARAMETER :: nx=33,ny=33,nz=131
INTEGER,PARAMETER :: num_ps_steps=2*20 ! final time 40sc
INTEGER,PARAMETER :: print_lag_err=1 !12 ! 1/60sc
INTEGER,PARAMETER :: print_lag_val=1 ! 0.25sc
INTEGER,PARAMETER :: resume=0

REAL,PARAMETER :: dx=Lx/(nx-1),dy=Ly/(ny-1),dz=Lz/(nz-1) ! meters
REAL,PARAMETER :: dt_mc=20.0 ! MIN(dx,dy)/SQRT(parg*Lz/2) ! seconds
INTEGER,PARAMETER :: num_mc_steps=(year/dt_mc)/20 ! 0.25sc
REAL,PARAMETER :: delta0=0.0
REAL,PARAMETER :: delta=0.0

INTEGER, PARAMETER :: sv_input_start=1
INTEGER, PARAMETER :: tv_input_start=1
INTEGER, PARAMETER :: the1ref_input_start=1
INTEGER, PARAMETER :: the1ref0_input_start=1

REAL,DIMENSION(0:nx) :: x
REAL,DIMENSION(0:ny) :: y
REAL,DIMENSION(nz) :: z

REAL :: dt

CONTAINS

SUBROUTINE set_coords()

INTEGER :: i,j,k

DO i=0,nx
        x(i)=-0.5*Lx+(i-1)*dx
END DO

DO j=0,ny
        y(j)=-0.5*Ly+(j-1)*dy
END DO

DO k=1,nz
        z(k)=(k-1)*dz
END DO

END SUBROUTINE set_coords

END MODULE parameters
