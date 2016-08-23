! Reads parameters from file and creates output folders

subroutine paramMultiReader
  use parameters
  implicit none

  ! Reading Inf File
110 format(a80)
  read(5,110) modelname
  read(5,110) vpfile
  read(5,110) vsfile
  read(5,110) rhofile
  read(5,'(L1)') optimise
  read(5,'(L1)') videoornot
  read(5,'(L1)') writingStrain
  read(5,*) IT_DISPLAY
  read(5,*) iSourceStart,iSourceInterval,nSource
  read(5,*) izSourceStart
  read(5,*) iReceiverStart,iReceiverInterval,nReceiver
  read(5,*) izReceiverStart
  read(5,*) nt,nx,nz
  read(5,*) dt,dx,dz
  read(5,*) f0, t0


  nt=nt*times
  nx=nx*times
  nz=nz*times
  dt=dt/dble(times)
  dx=dx/dble(times)
  dz=dz/dble(times)
  
  maxnt = nt
  maxnx = nx+(lmargin(1)+rmargin(1))
  maxnz = nz+(lmargin(2)+rmargin(2))


  
 
  call system('mkdir ./inffile')
   
  commandline="mkdir synthetics"
  call system(commandline)
  commandline="mkdir snapshots"
  call system(commandline)
  commandline="mkdir videos"
  call system(commandline)
  commandline="mkdir synthetics/"//trim(modelname)
  call system(commandline)
  commandline="mkdir videos/"//trim(modelname)
  call system(commandline)
  commandline="mkdir strains"
  call system(commandline)
  commandline="mkdir strains/"//trim(modelname)
  call system(commandline)
end subroutine paramMultiReader



subroutine paramFrechetReader
  use parameters
  use paramFrechet
  implicit none

  ! Reading Inf File
110 format(a80)
  read(5,110) modelname
  read(5,110) vpfile
  read(5,110) vsfile
  read(5,110) rhofile
  read(5,'(L1)') optimise
  read(5,'(L1)') videoornot
  read(5,'(L1)') writingStrain
  read(5,*) IT_DISPLAY
  read(5,*) iSourceStart,iSourceInterval,nSource
  read(5,*) izSourceStart
  read(5,*) iReceiverStart,iReceiverInterval,nReceiver
  read(5,*) izReceiverStart
  read(5,*) nt,nx,nz
  read(5,*) dt,dx,dz
  read(5,*) f0, t0
  ! Frechet needs some more infos here
  read(5,*) i1Source, i2Source
 


  nt=nt*times
  nx=nx*times
  nz=nz*times
  dt=dt/dble(times)
  dx=dx/dble(times)
  dz=dz/dble(times)
  
  maxnt = nt
  maxnx = nx+(lmargin(1)+rmargin(1))
  maxnz = nz+(lmargin(2)+rmargin(2))


  
 
  call system('mkdir ./inffile')
   
  commandline="mkdir synthetics"
  call system(commandline)
  commandline="mkdir snapshots"
  call system(commandline)
  commandline="mkdir kernelPsnapshots" ! we need kernelsnapshots folders too
  call system(commandline)
  commandline="mkdir kernelSsnapshots"
  call system(commandline)
  commandline="mkdir videos"
  call system(commandline)
  commandline="mkdir synthetics/"//trim(modelname)
  call system(commandline)
  commandline="mkdir videos/"//trim(modelname)
  call system(commandline)
  commandline="mkdir strains"
  call system(commandline)
  commandline="mkdir strains/"//trim(modelname)
  call system(commandline)
end subroutine paramFrechetReader



! Creates file for each source with source and corresponding receivers information
! saves in ../bin/inffile/

subroutine ReceiverSourcePositions

  use parameters
  implicit none

  
  do iReceiver = 1, nReceiver
     nrx(iReceiver)=(iReceiverStart-1)*times+1+iReceiverInterval*times*(iReceiver-1)
     nrz(iReceiver)=(izReceiverStart-1)*times+1
  enddo
  
  do iSource = 1, nSource
     iisx(iSource)=(iSourceStart-1)*times+1+iSourceInterval*times*(iSource-1)
     iisz(iSource)=(izSourceStart-1)*times+1
     write(filename, '(I5,".",I5,".inf")') iisx(iSource),iisz(iSource)
     do j=1, 12
        if(filename(j:j).eq.' ') filename(j:j)='0'
     enddo
    filename= './inffile/'//trim(modelname)//'.'//filename
    
       

     open(1, file=filename, form='formatted')
     write(1,'(a)') 'c grids'
     write(1,*) nt, nx, nz
     write(1,*) dt, dx, dz
     write(1,'(a)') 'c modelname'
     write(1,'(a)') trim(modelname)
     write(1,'(a)') trim(vpfile)
     write(1,'(a)') trim(vsfile)
     write(1,'(a)') trim(rhofile)

     write(1,'(a)') 'c source position (in grids) '
     write(1,*) iisx(iSource),iisz(iSource)
     write(1,'(a)') 'c source time function (Ricker wavelet)'
     write(1,*) f0,t0     
     write(1,'(a)') 'c receivers information'
     write(1,*) nReceiver
     do iReceiver = 1, nReceiver
        write(1,*) nrx(iReceiver), nrz(iReceiver)
     enddo
     

     !write(1,*) 'c'
     write(1,'(a)') 'end'     
  enddo

end subroutine ReceiverSourcePositions




! READ FROM file2d into array of size maxnx x maxnz

subroutine calstruct( maxnx,maxnz,file2d, nx,nz,rho )
  implicit none
  integer maxnx,maxnz,nx,nz
  double precision rho(1:maxnx+1,1:maxnz+1)
  real(kind(1.e0)),allocatable :: rrho(:,:)
  integer i,j,k,nox(6),noz(6)
  double precision x,z,xmax,zmax,trho,coef1,coef2
  integer recl_size
  character*80 file2d
  recl_size=kind(1.e0)*(nx+1)*(nz+1)
  
  allocate(rrho(1:nx+1,1:nz+1))
  open (1,file=file2d,form='unformatted',access='direct',recl=recl_size)
  read(1,rec=1) rrho(1:nx+1,1:nz+1)
  close(1)
  rho(1:nx+1,1:nz+1)=1.d-3*rrho(1:nx+1,1:nz+1)
  
  deallocate(rrho)
  return
end subroutine calstruct


! Create mu and lam arrays of size nx+1 x nz+1
! INPUT: rho, vp, vs arrays 
! OUTPUT: mu, lam, liquidmarkers arrays

subroutine calstruct2(maxnx,maxnz,nx,nz,rho,vp,vs,lam,mu,liquidmarkers)
  implicit none
  
  integer i,j,maxnz,maxnx,nx,nz
  double precision rho(maxnx+1,maxnz+1),vp(maxnx+1,maxnz+1),vs(maxnx+1,maxnz+1)
  double precision lam(maxnx+1,maxnz+1),mu(maxnx+1,maxnz+1)
  integer liquidmarkers(maxnx+1,maxnz+1)

  liquidmarkers = 0
  lam = 0.d0
  mu = 0.d0
  
  do i=1,nx+1
     do j=1,nz+1
        if(vs(i,j).eq.0.d0) then
           liquidmarkers(i,j)=1
           !NF should take out this now
           !vs(i,j)=vp(i,j)/1.7d0
        endif
           
        mu(i,j)=rho(i,j)*vs(i,j)*vs(i,j)
        lam(i,j)=rho(i,j)*vp(i,j)*vp(i,j)-2*mu(i,j)

        
     enddo
  enddo
  return
end subroutine calstruct2
  

! Extend model with absorbing layers
  !      L____NX____R
  ! U   |  ________  |
  ! NZ  | |        | |
  !     | |________| |
  ! D   |____________|
! INPUT: nx,nz,lmargin,rmargin. Thikness of abs boundaries 
! OUTPUT: rho,lam,mu,markers,liquidmarkers,zerodisplacement. Extended by these boundaries model

subroutine calstructBC(maxnx,maxnz,nx,nz,rho,lam,mu,markers,liquidmarkers,zerodisplacement,lmargin,rmargin)
  implicit none
  integer :: i,j,maxnx,maxnz,nx,nz,nnx,nnz
  double precision ::lam(maxnx+1,maxnz+1),mu(maxnx+1,maxnz+1),rho(maxnx+1,maxnz+1)  
  integer :: rmargin(1:2), lmargin(1:2)
  integer :: markers(maxnx+1,maxnz+1) ! discontinuities
  integer :: zerodisplacement(maxnx+1,maxnz+1) ! above free surface
  integer :: liquidmarkers(maxnx+1,maxnz+1)
  ! real(kind(0d0)), dimension(maxnz+1,maxnz+1) ::mmu,rrho,llam
  double precision, allocatable :: mmu(:,:), rrho(:,:), llam(:,:)
  integer, allocatable :: mmarkers(:,:),lliquidmarkers(:,:)
  integer, allocatable :: zzerodisplacement(:,:)
  
  allocate(rrho(1:maxnx+1,1:maxnz+1))
  allocate(mmu(1:maxnx+1,1:maxnz+1))
  allocate(llam(1:maxnx+1,1:maxnz+1))
  allocate(mmarkers(1:maxnx+1,1:maxnz+1))
  allocate(lliquidmarkers(1:maxnx+1,1:maxnz+1))
  allocate(zzerodisplacement(1:maxnx+1,1:maxnz+1))

  mmu = 0.d0
  rrho = 0.d0
  mmarkers = 0
  lliquidmarkers = 0
  llam = 0.d0
  zzerodisplacement = 0


  ! Main body of rectangle
  !    ____________
  !   |  ________  |
  !   | |    x   | |
  !   | |________| |
  !   |____________|

  llam(1+lmargin(1):nx+1+lmargin(1),1+lmargin(2):nz+1+lmargin(2)) = lam(1:nx+1,1:nz+1)
  mmu(1+lmargin(1):nx+1+lmargin(1),1+lmargin(2):nz+1+lmargin(2)) = mu(1:nx+1,1:nz+1)
  rrho(1+lmargin(1):nx+1+lmargin(1),1+lmargin(2):nz+1+lmargin(2)) = rho(1:nx+1,1:nz+1)

  mmarkers(1+lmargin(1):nx+1+lmargin(1),1+lmargin(2):nz+1+lmargin(2)) = markers(1:nx+1,1:nz+1)
  lliquidmarkers(1+lmargin(1):nx+1+lmargin(1),1+lmargin(2):nz+1+lmargin(2)) =  &
       liquidmarkers(1:nx+1,1:nz+1)
  zzerodisplacement(1+lmargin(1):nx+1+lmargin(1),1+lmargin(2):nz+1+lmargin(2)) =  &
       zerodisplacement(1:nx+1,1:nz+1)


  ! 4 corners
  !    ____________
  !   |x ________ x|
  !   | |        | |
  !   | |________| |
  !   |x__________x|

  llam(1:lmargin(1),1:lmargin(2))=lam(1,1)
  mmu(1:lmargin(1),1:lmargin(2))=mu(1,1)
  rrho(1:lmargin(1),1:lmargin(2))=rho(1,1)
  zzerodisplacement(1:lmargin(1),1:lmargin(2))=zerodisplacement(1,1)


  llam(1:lmargin(1),1+nz+1+lmargin(2):rmargin(2)+nz+1+lmargin(2))=lam(1,nz+1)
  mmu(1:lmargin(1),1+nz+1+lmargin(2):rmargin(2)+nz+1+lmargin(2))=mu(1,nz+1)
  rrho(1:lmargin(1),1+nz+1+lmargin(2):rmargin(2)+nz+1+lmargin(2))=rho(1,nz+1)
  zzerodisplacement(1:lmargin(1),1+nz+1+lmargin(2):rmargin(2)+nz+1+lmargin(2)) &
       = zerodisplacement(1,nz+1)
  !print *, llam(1,nz+lmargin(2)+5),mu(1,nz+1),rho(1,nz+1)

  llam(1+nx+1+lmargin(1):rmargin(1)+nx+1+lmargin(1),1:lmargin(2))=lam(nx+1,1)
  mmu(1+nx+1+lmargin(1):rmargin(1)+nx+1+lmargin(1),1:lmargin(2))=mu(nx+1,1)
  rrho(1+nx+1+lmargin(1):rmargin(1)+nx+1+lmargin(1),1:lmargin(2))=rho(nx+1,1)
  zzerodisplacement(1+nx+1+lmargin(1):rmargin(1)+nx+1+lmargin(1),1:lmargin(2)) &
       = zerodisplacement(nx+1,1)

  llam(1+nx+1+lmargin(1):rmargin(1)+nx+1+lmargin(1),1+nz+1+lmargin(2):rmargin(2)+nz+1+lmargin(2)) &
       = lam(nx+1,nz+1)
  mmu(1+nx+1+lmargin(1):rmargin(1)+nx+1+lmargin(1),1+nz+1+lmargin(2):rmargin(2)+nz+1+lmargin(2)) &
       = mu(nx+1,nz+1)
  rrho(1+nx+1+lmargin(1):rmargin(1)+nx+1+lmargin(1),1+nz+1+lmargin(2):rmargin(2)+nz+1+lmargin(2)) &
       = rho(nx+1,nz+1)
  zzerodisplacement(1+nx+1+lmargin(1):rmargin(1)+nx+1+lmargin(1),1+nz+1+lmargin(2):rmargin(2)+nz+1+lmargin(2)) &
       = zerodisplacement(nx+1,nz+1)
  

  ! 4 rectangles

  do i = 1,lmargin(1)
     llam(i,1+lmargin(2):nz+1+lmargin(2)) = lam(1,1:nz+1)
     mmu(i,1+lmargin(2):nz+1+lmargin(2)) = mu(1,1:nz+1)
     rrho(i,1+lmargin(2):nz+1+lmargin(2)) = rho(1,1:nz+1)
     zzerodisplacement(i,1+lmargin(2):nz+1+lmargin(2))  &
          = zerodisplacement(1,1:nz+1)
  enddo
  
  do i = 1+nx+1+lmargin(1),rmargin(1)+nx+1+lmargin(1)
     llam(i,1+lmargin(2):nz+1+lmargin(2)) = lam(nx+1,1:nz+1)
     mmu(i,1+lmargin(2):nz+1+lmargin(2)) = mu(nx+1,1:nz+1)
     rrho(i,1+lmargin(2):nz+1+lmargin(2)) = rho(nx+1,1:nz+1)
     zzerodisplacement(i,1+lmargin(2):nz+1+lmargin(2)) &
          = zerodisplacement(nx+1,1:nz+1)
  enddo
  
  do i = 1,lmargin(2)
     llam(1+lmargin(1):nx+1+lmargin(2),i)=lam(1:nx+1,1)
     mmu(1+lmargin(1):nx+1+lmargin(2),i)=mu(1:nx+1,1)
     rrho(1+lmargin(1):nx+1+lmargin(2),i)=rho(1:nx+1,1)
     zzerodisplacement(1+lmargin(1):nx+1+lmargin(2),i) &
          =zerodisplacement(1:nx+1,1)
  enddo

  do i = 1+nz+1+lmargin(2),rmargin(2)+nz+1+lmargin(2)
     llam(1+lmargin(1):nx+1+lmargin(2),i) = lam(1:nx+1,nz+1)
     mmu(1+lmargin(1):nx+1+lmargin(2),i) = mu(1:nx+1,nz+1)
     rrho(1+lmargin(1):nx+1+lmargin(2),i) = rho(1:nx+1,nz+1)
     zzerodisplacement(1+lmargin(1):nx+1+lmargin(2),i) &
          = zerodisplacement(1:nx+1,nz+1)

  enddo

  nnx=rmargin(1)+nx+lmargin(1)
  nnz=rmargin(2)+nz+lmargin(2)

  nx=nnx
  nz=nnz
  lam=0.d0
  rho=0.d0
  mu=0.d0
  liquidmarkers = 0
  zerodisplacement = 0
  markers = 0

  lam(1:nx+1,1:nz+1) = llam(1:nx+1,1:nz+1)
  rho(1:nx+1,1:nz+1) = rrho(1:nx+1,1:nz+1)
  mu(1:nx+1,1:nz+1) = mmu(1:nx+1,1:nz+1)
  markers(1:nx+1,1:nz+1)=mmarkers(1:nx+1,1:nz+1)
  zerodisplacement(1:nx+1,1:nz+1)=zzerodisplacement(1:nx+1,1:nz+1)
  liquidmarkers(1:nx+1,1:nz+1)=lliquidmarkers(1:nx+1,1:nz+1)
  !print *, nx,nz
  !write(12,*) rho(:,:)
  !write(13,*) lam(:,:)
  !write(14,*) mmu(1:nx+1,1:nz+1)
  !stop
  
  deallocate(llam)
  deallocate(rrho)
  deallocate(mmu)
  deallocate(mmarkers)
  deallocate(lliquidmarkers)
  deallocate(zzerodisplacement)
end subroutine calstructBC


! Fill up array with zeros
! INPUT: nx,nz,ux
! OUTPUT: ux(1:nx+1, 1:nx+1) = 0

subroutine datainit( nx,nz,ux )
  
  integer nx,nz
  double precision ux(nx+1,*)
  integer i,j
  
  do j=1,nz+1
     do i=1,nx+1
        ux(i,j) = 0.d0
     enddo
  enddo
  
  return
end subroutine datainit



! All-in-one implementation of Cerjan boundaries
! Returns attenuated wavefield 

subroutine  compNRBC2(ux,ux1,ux2,uz,uz1,uz2, rrate, lmargin, rmargin,nnx,nnz)

  ! Cerjan boundary conditions (2D)
  implicit none
  integer :: nnx, nnz

  real*8, intent(inout) :: ux2(nnx,nnz),uz2(nnx,nnz)
  real*8, intent(inout) :: ux1(nnx,nnz),uz1(nnx,nnz)
  real*8, intent(inout) :: ux(nnx,nnz),uz(nnx,nnz)
  real*8, intent(in) :: rrate
  integer, dimension(3), intent(in) :: lmargin, rmargin
  integer ix, iy, iz
  integer i, j, k
  real*8 r
  
  do iz = 1, nnz
     !do iy = 1, nny
     do ix = 1, nnx
        
        i = 0
        j = 0
        k = 0
           
        if (ix < lmargin(1) + 1) i = lmargin(1) + 1 - ix
        !   if (iy < lmargin(2) + 1) j = lmargin(2) + 1 - iy
        if (iz < lmargin(2) + 1) k = lmargin(2) + 1 - iz
        if (nnx - rmargin(1) < ix) i = ix - nnx + rmargin(1)
        !if (nny - rmargin(2) < iy) j = iy - nny + rmargin(2)
        if (nnz - rmargin(2) < iz) k = iz - nnz + rmargin(2)
           
        if (i == 0 .and. j == 0 .and. k == 0) cycle
        
        r = rrate * rrate * dble( i * i + j * j + k * k )
        r = exp( - r )
        

        ux2(:,:) = ux2(:,:) * r
        ux1(:,:) = ux1(:,:) * r
        ux(:,:) = ux(:,:) * r

        uz2(:,:) = uz2(:,:) * r
        uz1(:,:) = uz1(:,:) * r
        uz(:,:) = uz(:,:) * r
        
     enddo
     !enddo
  enddo
  
end subroutine compNRBC2


! Cerjan boundary attenuation coefficients
! Returns value of decaying rate at each grid point

subroutine  compNRBCpre(r,rrate, lmargin, rmargin,nnx,nnz)

  implicit none
  integer :: nnx, nnz
  ! Cerjan boundary conditions (2D)
  double precision :: r(nnx,nnz)
  real*8, intent(in) :: rrate
  integer, dimension(3), intent(in) :: lmargin, rmargin
  integer ix, iy, iz
  integer i, j, k
  double precision :: rr
  
  
  do iz = 1, nnz
     !do iy = 1, nny
     do ix = 1, nnx
        
        i = 0
        j = 0
        k = 0
           
        if (ix < lmargin(1) + 1) i = lmargin(1) + 1 - ix
        !   if (iy < lmargin(2) + 1) j = lmargin(2) + 1 - iy
        if (iz < lmargin(2) + 1) k = lmargin(2) + 1 - iz
   
        if (nnx - rmargin(1) < ix) i = ix - nnx + rmargin(1)
        !if (nny - rmargin(2) < iy) j = iy - nny + rmargin(2)
        if (nnz - rmargin(2) < iz) k = iz - nnz + rmargin(2)
           
        if (i == 0 .and. j == 0 .and. k == 0) cycle
        
        rr = rrate * rrate * dble( i * i + j * j + k * k )
        r(ix,iz) = exp( - rr )
        
        !if(r(ix,iz).ne.1.d0) then
        !   print *, ix,iz,r(ix,iz)
        !endif
        
        !print *, ix,iy,r
        !ux2(:,:) = ux2(:,:) * r
        !ux1(:,:) = ux1(:,:) * r
        !ux(:,:) = ux(:,:) * r

        !uz2(:,:) = uz2(:,:) * r
        !uz1(:,:) = uz1(:,:) * r
        !uz(:,:) = uz(:,:) * r
        
     enddo

     
     
     !enddo
  enddo
  
end subroutine compNRBCpre




