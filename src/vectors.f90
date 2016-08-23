subroutine vectorAllocate
  use parameters
  implicit none
   
  ! Array allocation

  allocate(nrx(1:nReceiver)) ! receivers
  allocate(nrz(1:nReceiver))
  allocate(iisx(1:nSource)) ! Sources
  allocate(iisz(1:nSource))

  allocate(ux(maxnx+1,maxnz+1),uz(maxnx+1,maxnz+1))
  allocate(ux1(maxnx+1,maxnz+1),ux2(maxnx+1,maxnz+1))
  allocate(uz1(maxnx+1,maxnz+1),uz2(maxnx+1,maxnz+1))
  allocate(e1(maxnx+1,maxnz+1), e2(maxnx+1,maxnz+1))
  allocate(e3(maxnx+1,maxnz+1), e4(maxnx+1,maxnz+1))
  allocate(e5(maxnx+1,maxnz+1), e6(maxnx+1,maxnz+1))
  allocate(e7(maxnx+1,maxnz+1), e8(maxnx+1,maxnz+1))
  allocate(e13(maxnx+1,maxnz+1),e14(maxnx+1,maxnz+1))
  allocate(e15(maxnx+1,maxnz+1),e16(maxnx+1,maxnz+1))
  allocate(e17(maxnx+1,maxnz+1),e18(maxnx+1,maxnz+1))
  allocate(e19(maxnx+1,maxnz+1),e20(maxnx+1,maxnz+1))
  allocate(f1(maxnx+1,maxnz+1), f2(maxnx+1,maxnz+1))
  allocate(f3(maxnx+1,maxnz+1), f4(maxnx+1,maxnz+1))
  allocate(f5(maxnx+1,maxnz+1), f6(maxnx+1,maxnz+1))
  allocate(f7(maxnx+1,maxnz+1), f8(maxnx+1,maxnz+1))
  allocate(f13(maxnx+1,maxnz+1),f14(maxnx+1,maxnz+1))
  allocate(f15(maxnx+1,maxnz+1),f16(maxnx+1,maxnz+1))
  allocate(f17(maxnx+1,maxnz+1),f18(maxnx+1,maxnz+1))
  allocate(f19(maxnx+1,maxnz+1),f20(maxnx+1,maxnz+1))
  allocate(work(maxnx+1,32)) ! NF, is it nz or nx ??
  
  allocate(ee12(maxnx+1,maxnz+1),ee34(maxnx+1,maxnz+1),ee56(maxnx+1,maxnz+1))
  allocate(ee65(maxnx+1,maxnz+1),ee78(maxnx+1,maxnz+1),ee87(maxnx+1,maxnz+1))
  allocate(ff12(maxnx+1,maxnz+1),ff34(maxnx+1,maxnz+1),ff56(maxnx+1,maxnz+1))
  allocate(ff65(maxnx+1,maxnz+1),ff78(maxnx+1,maxnz+1),ff87(maxnx+1,maxnz+1))

  allocate(rho(maxnx+1,maxnz+1))
  allocate(lam(maxnx+1,maxnz+1),mu(maxnx+1,maxnz+1))
  allocate(fx(maxnx+1,maxnz+1),fz(maxnx+1,maxnz+1))
  allocate(vs(maxnx+1,maxnz+1),vp(maxnx+1,maxnz+1))
  
  
  allocate(synx(0:maxnt,1:nReceiver),synz(0:maxnt,1:nReceiver),time(0:maxnt)) ! synthetics

  allocate(video(maxnx+1,maxnz+1))

  allocate(snapux(maxnx+1,maxnz+1),snapuz(maxnx+1,maxnz+1))
  
  allocate(weightBC(maxnx+1,maxnz+1))

  allocate(liquidmarkers(maxnx+1,maxnz+1))
  allocate(markers(maxnx+1,maxnz+1))
  allocate(zerodisplacement(maxnx+1,maxnz+1))
  
  allocate(total_energy_kinetic(maxnt+1),total_energy_potential(maxnt+1))
    
  allocate(singleStrainDiagonal(maxnx+1,maxnz+1))
  allocate(singleStrainShear(maxnx+1,maxnz+1))
  allocate(tmpsingleStrain(1:nx+1,1:nz+1))
  
end subroutine vectorAllocate


subroutine vectorAllocateFrechet
  use paramFrechet
  use parameters
  implicit none
  
  ! For Frechet calculation, we don't need to think about absorbing boundaries
  maxnx=nx
  maxnz=nz

  ! Array allocation

  allocate(nrx(1:nReceiver)) ! receivers
  allocate(nrz(1:nReceiver))
  allocate(iisx(1:nSource)) ! Sources
  allocate(iisz(1:nSource))
  
  ! nx and nz here are not considering boundary conditions
  allocate(singleStrainForward(1:nx+1,1:nz+1))
  allocate(singleStrainBack(1:nx+1,1:nz+1))
  allocate(strainForward(1:nx+1,1:nz+1))
  allocate(strainBack(1:nx+1,1:nz+1))
  
  allocate(kernelP(1:nx+1,1:nz+1))
  allocate(kernelS(1:nx+1,1:nz+1))
  allocate(singleKernelP(1:nx+1,1:nz+1))
  allocate(singleKernelS(1:nx+1,1:nz+1))
  

  allocate(synx(0:maxnt,1:nReceiver),synz(0:maxnt,1:nReceiver),time(0:maxnt)) ! synthetics

  allocate(video(maxnx+1,maxnz+1))

  allocate(snapux(maxnx+1,maxnz+1),snapuz(maxnx+1,maxnz+1))
  
end subroutine vectorAllocateFrechet




subroutine calStrainDiagonal(nx,nz,ux,uz,lmargin,rmargin,singleStrainDiagonal)
  implicit none
  integer :: maxnz,nx,nz,ix,iz
  double precision :: ux(1:nx+1,1:nz+1),uz(1:nx+1,1:nz+1)
  integer :: lmargin(1:2), rmargin(1:2)
  real(kind(0e0)) :: singleStrainDiagonal(1:nx+1,1:nz+1)
  double precision :: straintmp
  double precision, parameter :: onetwelfth = 0.0833333333333333d0
  
  integer :: nxstart,nxend,nzstart,nzend
  ! we calculate only for the box of interest (without absorbing boundaries)
  ! and we suppose that we have lmargin != 0 and rmargin != 0
  ! i.e., we use the four-point first derivative operators with one extra point to the left

  nxstart=lmargin(1)+1
  nxend=nx+1-rmargin(1)

  nzstart=lmargin(2)+1
  nzend=nz+1-rmargin(2)

  
  do ix=nxstart,nxend
     do iz=nzstart,nzend
        straintmp=0.d0
        straintmp=(5.d0*ux(ix+1,iz)+3.d0*ux(ix,iz)-9.d0*ux(ix-1,iz)+ux(ix-2,iz))*onetwelfth
        straintmp=straintmp+(5.d0*uz(ix,iz+1)+3.d0*uz(ix,iz)-9.d0*uz(ix,iz-1)+uz(ix,iz-2))*onetwelfth
        singleStrainDiagonal(ix,iz)=straintmp
     enddo
  enddo
 
     
end subroutine calStrainDiagonal




subroutine calStrainShear(nx,nz,ux,uz,lmargin,rmargin,singleStrain)
  implicit none
  integer :: maxnz,nx,nz,ix,iz
  double precision :: ux(1:nx+1,1:nz+1),uz(1:nx+1,1:nz+1)
  integer :: lmargin(1:2), rmargin(1:2)
  real(kind(0e0)) :: singleStrain(1:nx+1,1:nz+1)
  double precision :: straintmp
  double precision, parameter :: onetwelfth = 0.0833333333333333d0
  
  integer :: nxstart,nxend,nzstart,nzend
  ! we calculate only for the box of interest (without absorbing boundaries)
  ! and we suppose that we have lmargin != 0 and rmargin != 0
  ! i.e., we use the four-point first derivative operators with one extra point to the left

  nxstart=lmargin(1)+1
  nxend=nx+1-rmargin(1)

  nzstart=lmargin(2)+1
  nzend=nz+1-rmargin(2)

  
  do ix=nxstart,nxend
     do iz=nzstart,nzend
        straintmp=0.d0
        straintmp=(5.d0*uz(ix+1,iz)+3.d0*uz(ix,iz)-9.d0*uz(ix-1,iz)+uz(ix-2,iz))*onetwelfth
        straintmp=straintmp+(5.d0*ux(ix,iz+1)+3.d0*ux(ix,iz)-9.d0*ux(ix,iz-1)+ux(ix,iz-2))*onetwelfth
        singleStrain(ix,iz)=straintmp
     enddo
  enddo
 
     
end subroutine calStrainShear
