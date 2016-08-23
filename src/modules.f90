module parameters

  implicit none

  
  character(120) :: filename
  integer :: tmpint
  character(80) :: modelname
  
  integer, parameter :: times = 1 ! this can make the dx,dz,dt finer  

  ! switch OPT / CONV, use Optimally accurate operators (Takeuchi, Geller, 1998) or [1 -2 1] scheme
  logical :: optimise

  ! switch video
  logical :: videoornot 

  ! writing strains
  logical :: writingStrain

  ! Constants
  double precision, parameter :: pi=3.1415926535897932d0 
  double precision, parameter :: ZERO = 0.d0


  !###############################################################################################
  ! Line of RECEIVERS
  integer :: iReceiverStart,iReceiverInterval,nReceiver
  integer :: izReceiverStart
  integer, allocatable :: nrx(:),nrz(:) ! Receiver positions

  ! Line of SOURCES
  integer :: iSourceStart,iSourceInterval,nSource
  integer :: izSourceStart
  integer, allocatable, dimension(:) :: iisx, iisz

  integer :: iSource, iReceiver
  integer :: maxnx,maxnz,maxnt

  ! Ricker wavelets source
  double precision f0,t0    ! Dominant frequency and excitation time
  !double precision tp,ts

  !###############################################################################################
    
  ! parameters for the gridding
  double precision dt,dx,dz     ! time step, step over OX and OZ
  
  ! parameters for the wavefield
  integer nt,nx,nz,it,ist,isx,isz,ix,iz,recl_size
  ! Attention ! nx and nz are modified with absorbing boundaries

  ! Wavefield for time steps N, N-1 and N-2
  double precision, allocatable, dimension(:,:) :: ux,uz,ux1,ux2,uz1,uz2
  double precision, allocatable, dimension(:,:) :: e1,e2,e3,e4,e5,e6,e7,e8
  double precision, allocatable, dimension(:,:) :: e13,e14,e15,e16,e17,e18,e19,e20
  double precision, allocatable, dimension(:,:) :: f1,f2,f3,f4,f5,f6,f7,f8
  double precision, allocatable, dimension(:,:) :: f13,f14,f15,f16,f17,f18,f19,f20
  double precision, allocatable, dimension(:,:) :: work

  ! for discontinuities
  double precision, allocatable, dimension(:,:) :: ee12,ee34,ee56,ee65,ee78,ee87
  double precision, allocatable, dimension(:,:) :: ff12,ff34,ff56,ff65,ff78,ff87

  ! parameter for the structure
  character(80) :: vpfile, vsfile, rhofile          ! modelname
  double precision, allocatable, dimension(:,:) :: rho,lam,mu,fx,fz,vs,vp

  ! Courant number
  double precision :: cp                            ! maxvalue of vp
  double precision :: Courant_number

  
  ! parameter for the receiver
  integer :: ir,j
  real, allocatable, dimension(:,:) :: synx,synz    !seismograms
  real, allocatable, dimension(:) :: time
  character(200) :: outfile
 
  
  ! parameter for the waveform
  double precision t

  !parameter for video  
  real, allocatable, dimension(:,:) :: video
  real, allocatable, dimension(:,:) :: snapux,snapuz
  integer :: IT_DISPLAY                             ! show pic every IT_DISPLAY steps
 
  integer(2) head(1:120)
  character(80) :: routine


  !###############################################################################################
  ! PML Boundary
  logical, parameter :: USE_PML_XMIN = .true.
  logical, parameter :: USE_PML_XMAX = .true.
  logical, parameter :: USE_PML_YMIN = .true.
  logical, parameter :: USE_PML_YMAX = .true.
  
  integer, parameter :: NPOINTS_PML = 100*times       ! thickness of the PML layer in grid points
  double precision, parameter :: NPOWER = 2.d0        ! power to compute d0 profile
  double precision, parameter :: K_MAX_PML = 1.d0     ! from Gedney page 8.11
  double precision :: ALPHA_MAX_PML

  ! Cerjan Boundary
  double precision, parameter :: CerjanRate = 0.0015          ! Decay rate
  double precision, allocatable, dimension(:,:) :: weightBC   ! decaying mask of size of the wavefield
  integer :: lmargin(1:2),rmargin(1:2)                        ! Thickness


  !###############################################################################################
  ! ENERGY

  double precision epsilon_xx,epsilon_yy,epsilon_xy     ! for evolution of total energy in the medium
  double precision, allocatable, dimension(:) :: total_energy_kinetic,total_energy_potential


  !###############################################################################################
  ! MARKERS

  ! for water
  integer, allocatable, dimension(:,:) :: liquidmarkers   ! are non-zero where is a liquid

  ! for discontinuities
  integer, allocatable, dimension(:,:) :: markers
  integer :: nDiscon                           ! number of discontinuities
  integer :: lengthDiscon                      ! with x,z coordinates
  double precision, allocatable :: dscr(:,:,:) ! discontinuity coordinates
  double precision :: tmpvaluex,tmpvaluez

  ! for free surface
  integer,allocatable, dimension(:,:) :: zerodisplacement
  integer :: lengthFreeSurface ! with x,z coordinates
  double precision, allocatable :: free(:,:)
 
  
  !###############################################################################################
  ! FWI
  ! for waveform inversion
  real(kind(0e0)), allocatable, dimension(:,:):: singleStrainDiagonal,singleStrainShear,tmpsingleStrain

  
  character(340) :: commandline
  

  

end module parameters


module paramFrechet
  
  implicit none
  integer :: i1Source, i2Source ! 2 sources for cross correlations
  integer :: isx1,isx2,isz1,isz2,it1,it2
  real, allocatable, dimension(:,:) :: singleStrainForward,singleStrainBack
  real, allocatable, dimension(:,:) :: singleKernelP,singleKernelS
  double precision, allocatable, dimension (:,:) :: strainForward,strainBack
  double precision, allocatable, dimension (:,:) :: kernelP,kernelS
  
end module paramFrechet
