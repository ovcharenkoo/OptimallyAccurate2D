subroutine MizutaniIsoFree(coef,rho0,rho1,lam0,lam1,mu0,mu1,ik,jk,dx,dz,eta,normal)

  ! isotropic modified operators for discontinuities

  implicit none 
  integer :: i,info 
  double precision :: coef(1:6,1:2),eta(0:1,1:2),normal(1:2)
  double precision :: rho0,rho1,lam0,lam1,mu0,mu1,dx,dz
  integer :: ik,jk
  double precision :: nvx,nvy,nvx2,nvy2,nvxnvy
  double precision :: C1(1:4), C0(1:4)
  double precision :: eta0x,eta1x,eta0y,eta1y

  ! intermediate values
  double precision :: xe1ux,ye1uy,xe2ux_dx,xe2ux_dy,ye2uy_dx,ye2uy_dy
  double precision :: xe3ux_dxx,xe3ux_dyy,xe3ux_dxy,ye3uy_dxx,ye3uy_dyy,ye3uy_dxy 
  double precision :: xe4ux_dxa,xe4ux_dya,xe4uy_dxa,xe4uy_dya 
  double precision :: ye4ux_dxa,ye4ux_dya,ye4uy_dxa,ye4uy_dya
  double precision :: xe4ux_dxb,xe4ux_dyb,xe4uy_dxb,xe4uy_dyb 
  double precision :: ye4ux_dxb,ye4ux_dyb,ye4uy_dxb,ye4uy_dyb
  double precision :: xe5ux_dxxa,xe5ux_dxya,xe5uy_dxxa,xe5uy_dxya,xe5ux_dyya,xe5uy_dyya
  double precision :: ye5ux_dxxa,ye5ux_dxya,ye5uy_dxxa,ye5uy_dxya,ye5ux_dyya,ye5uy_dyya
  double precision :: xe5ux_dxxb,xe5ux_dxyb,xe5uy_dxxb,xe5uy_dxyb,xe5ux_dyyb,xe5uy_dyyb
  double precision :: ye5ux_dxxb,ye5ux_dxyb,ye5uy_dxxb,ye5uy_dxyb,ye5ux_dyyb,ye5uy_dyyb 
  double precision :: xe6ux_dxxa,xe6uy_dxya,xe6ux_dyya,ye6ux_dxya,ye6uy_dxxa,ye6uy_dyya
  double precision :: xe6ux_dxxb,xe6uy_dxyb,xe6ux_dyyb,ye6ux_dxyb,ye6uy_dxxb,ye6uy_dyyb 
 
  ! Matrices
  double precision, dimension(1:12,1:12) :: A0,B0,A1,B1,M0,Mm0,M1,M2,Bm0,Mm1
  double precision :: r, s
  

  ! A1 is always zero for the free surface

  nvx=normal(1)
  nvy=normal(2)
  nvx2=nvx*nvx
  nvy2=nvy*nvy
  nvxnvy=nvx*nvy


  eta0x = dble(ik)*eta(0,1)
  eta1x =-dble(ik)*eta(1,1)

  eta0y = dble(jk)*eta(0,2)
  eta1y =-dble(jk)*eta(1,2)  

  C0 = 0.d0
  C1 = 0.d0

  A0 = 0.d0
  A1 = 0.d0
  B0 = 0.d0
  B1 = 0.d0
  
  !print *, eta0x,eta1x,eta0y,eta1y
  !print *, nvx,nvy

  
  ! isotropic, PSV

  C0(1) = lam0 + 2.d0*mu0
  C0(2) = lam0
  C0(3) = C0(1)
  C0(4) = mu0

 
  C1(1) = lam1 + 2.d0*mu1
  C1(2) = lam1
  C1(3) = C1(1)
  C1(4) = mu1
  
  

  ! Hereafter we follow the Oleg Ovcharenko's formulation
  
  
  ! BC eq.1
  xe1ux = 1.d0
  ye1uy = 1.d0
  
  ! BC eq.2
  xe2ux_dx = nvy
  xe2ux_dy = -nvx
  ye2uy_dx = nvy
  ye2uy_dy = -nvx
  
  ! BC eq.3
  xe3ux_dxx = nvy2
  xe3ux_dyy = nvx2
  xe3ux_dxy = -2.d0*nvxnvy
  ye3uy_dxx = nvy2
  ye3uy_dyy = nvx2
  ye3uy_dxy = -2.d0*nvxnvy

  ! BC eq.4

  ! Above
  xe4ux_dxa = nvx*C1(1)
  xe4ux_dya = nvy*C1(4)
  xe4uy_dxa = nvy*C1(4)
  xe4uy_dya = nvx*C1(2)
  
  ye4ux_dxa = nvy*C1(2)
  ye4ux_dya = nvx*C1(4)
  ye4uy_dxa = nvx*C1(4)
  ye4uy_dya = nvy*C1(3)

  ! Below
  xe4ux_dxb = nvx*C0(1)
  xe4ux_dyb = nvy*C0(4)
  xe4uy_dxb = nvy*C0(4)
  xe4uy_dyb = nvx*C0(2)
  
  ye4ux_dxb = nvy*C0(2)
  ye4ux_dyb = nvx*C0(4)
  ye4uy_dxb = nvx*C0(4)
  ye4uy_dyb = nvy*C0(3)

  ! BC eq.5
  
  ! Above
  xe5ux_dxxa = nvxnvy*C1(1)
  xe5ux_dxya = nvy2*C1(4)-nvx2*C1(1)
  xe5uy_dxxa = nvy2*C1(4)
  xe5uy_dxya = nvxnvy*C1(2)-nvxnvy*C1(4)
  xe5ux_dyya = - nvxnvy*C1(4)
  xe5uy_dyya = - nvx2*C1(2)

  ye5ux_dxxa = nvy2*C1(2)
  ye5ux_dxya = nvxnvy*C1(4)-nvxnvy*C1(2)
  ye5uy_dxxa = nvxnvy*C1(4)
  ye5uy_dxya = nvy2*C1(3)-nvx2*C1(4)
  ye5ux_dyya =  - nvx2*C1(4)
  ye5uy_dyya =  - nvxnvy*C1(3)
  
  ! Below
  xe5ux_dxxb = nvxnvy*C0(1)
  xe5ux_dxyb = nvy2*C0(4)-nvx2*C0(1)
  xe5uy_dxxb = nvy2*C0(4)
  xe5uy_dxyb = nvxnvy*C0(2)-nvxnvy*C0(4)
  xe5ux_dyyb = - nvxnvy*C0(4)
  xe5uy_dyyb = - nvx2*C0(2)
  
  ye5ux_dxxb = nvy2*C0(2)
  ye5ux_dxyb = nvxnvy*C0(4)-nvxnvy*C0(2)
  ye5uy_dxxb = nvxnvy*C0(4)
  ye5uy_dxyb = nvy2*C0(3)-nvx2*C0(4)
  ye5ux_dyyb =  - nvx2*C0(4)
  ye5uy_dyyb =  - nvxnvy*C0(3)

  ! BC eq.6
  ! Above C1
  xe6ux_dxxa = C1(1)/rho1
  xe6uy_dxya = (C1(2)+C1(4))/rho1
  xe6ux_dyya = C1(4)/rho1
  
  ye6ux_dxya = (C1(4)+C1(2))/rho1
  ye6uy_dxxa = C1(4)/rho1
  ye6uy_dyya = C1(3)/rho1

  ! Below C0
  xe6ux_dxxb = C0(1)/rho0
  xe6uy_dxyb = (C0(2)+C0(4))/rho0
  xe6ux_dyyb = C0(4)/rho0
  
  ye6ux_dxyb = (C0(4)+C0(2))/rho0
  ye6uy_dxxb = C0(4)/rho0
  ye6uy_dyyb = C0(3)/rho0
    
    
  ! Full matrices for left and right points
  
  ! B0

  B0 = 0.d0
  
  B0( 1, 1) = xe1ux

  B0( 2, 2) = xe2ux_dx
  B0( 2, 3) = xe2ux_dy
 
  B0( 3, 4) = xe3ux_dxx
  B0( 3, 5) = xe3ux_dyy
  B0( 3, 6) = xe3ux_dxy
  
  B0( 4, 2) = xe4ux_dxb
  B0( 4, 3) = xe4ux_dyb 
  B0( 4, 8) = xe4uy_dxb
  B0( 4, 9) = xe4uy_dyb

  B0( 5, 4) = xe5ux_dxxb
  B0( 5, 5) = xe5ux_dyyb
  B0( 5, 6) = xe5ux_dxyb
  B0( 5,10) = xe5uy_dxxb
  B0( 5,11) = xe5uy_dyyb
  B0( 5,12) = xe5uy_dxyb
  
  B0( 6, 4) = xe6ux_dxxb
  B0( 6, 5) = xe6ux_dyyb
  B0( 6,12) = xe6uy_dxyb
  
  B0( 7, 7) = ye1uy

  B0( 8, 8) = ye2uy_dx
  B0( 8, 9) = ye2uy_dy

  B0( 9,10) = ye3uy_dxx
  B0( 9,11) = ye3uy_dyy
  B0( 9,12) = ye3uy_dxy

  B0(10, 2) = ye4ux_dxb
  B0(10, 3) = ye4ux_dyb
  B0(10, 8) = ye4uy_dxb
  B0(10, 9) = ye4uy_dyb

  B0(11, 4) = ye5ux_dxxb 
  B0(11, 5) = ye5ux_dyyb
  B0(11, 6) = ye5ux_dxyb
  B0(11,10) = ye5uy_dxxb
  B0(11,11) = ye5uy_dyyb
  B0(11,12) = ye5uy_dxyb
  
  B0(12, 6) = ye6ux_dxyb
  B0(12,10) = ye6uy_dxxb
  B0(12,11) = ye6uy_dyyb
  

  ! A0

  r = eta0x
  s = eta0y

  A0 = 0.d0
  
  A0( 1, 1) = 1.d0
  A0( 1, 2) = r*dx
  A0( 1, 3) = s*dz
  A0( 1, 4) = (r*dx)*(r*dx)*5.d-1
  A0( 1, 5) = (s*dz)*(s*dz)*5.d-1
  A0( 1, 6) = r*s*dx*dz

  A0( 2, 2) = 1.d0
  A0( 2, 4) = r*dx
  A0( 2, 6) = s*dz

  A0( 3, 3) = 1.d0
  A0( 3, 5) = s*dz
  A0( 3, 6) = r*dx
  
  A0( 4, 4) = 1.d0
  A0( 5, 5) = 1.d0
  A0( 6, 6) = 1.d0
  
  A0(7:12,7:12) = A0(1:6,1:6)

         

  ! B1

  B1 = 0.d0
  
  B1( 1, 1) = xe1ux

  B1( 2, 2) = xe2ux_dx
  B1( 2, 3) = xe2ux_dy
 
  B1( 3, 4) = xe3ux_dxx
  B1( 3, 5) = xe3ux_dyy
  B1( 3, 6) = xe3ux_dxy
  
  B1( 4, 2) = xe4ux_dxa
  B1( 4, 3) = xe4ux_dya 
  B1( 4, 8) = xe4uy_dxa
  B1( 4, 9) = xe4uy_dya

  B1( 5, 4) = xe5ux_dxxa
  B1( 5, 5) = xe5ux_dyya
  B1( 5, 6) = xe5ux_dxya
  B1( 5,10) = xe5uy_dxxa
  B1( 5,11) = xe5uy_dyya
  B1( 5,12) = xe5uy_dxya
  
  B1( 6, 4) = xe6ux_dxxa
  B1( 6, 5) = xe6ux_dyya
  B1( 6,12) = xe6uy_dxya
  
  B1( 7, 7) = ye1uy

  B1( 8, 8) = ye2uy_dx
  B1( 8, 9) = ye2uy_dy

  B1( 9,10) = ye3uy_dxx
  B1( 9,11) = ye3uy_dyy
  B1( 9,12) = ye3uy_dxy

  B1(10, 2) = ye4ux_dxa
  B1(10, 3) = ye4ux_dya
  B1(10, 8) = ye4uy_dxa
  B1(10, 9) = ye4uy_dya

  B1(11, 4) = ye5ux_dxxa 
  B1(11, 5) = ye5ux_dyya
  B1(11, 6) = ye5ux_dxya
  B1(11,10) = ye5uy_dxxa
  B1(11,11) = ye5uy_dyya
  B1(11,12) = ye5uy_dxya
  
  B1(12, 6) = ye6ux_dxya
  B1(12,10) = ye6uy_dxxa
  B1(12,11) = ye6uy_dyya


      


  ! A1


  r = eta1x;
  s = eta1y;
  


  A1 = 0.d0
  
  A1( 1, 1) = 1.d0
  A1( 1, 2) = r*dx
  A1( 1, 3) = s*dz
  A1( 1, 4) = (r*dx)*(r*dx)*5.d-1
  A1( 1, 5) = (s*dz)*(s*dz)*5.d-1
  A1( 1, 6) = r*s*dx*dz
  
  A1( 2, 2) = 1.d0
  A1( 2, 4) = r*dx
  A1( 2, 6) = s*dz
  
  A1( 3, 3) = 1.d0
  A1( 3, 5) = s*dz
  A1( 3, 6) = r*dx
  
  A1( 4, 4) = 1.d0
  A1( 5, 5) = 1.d0
  A1( 6, 6) = 1.d0
  
  A1(7:12,7:12) = A1(1:6,1:6)
  

  

  ! temporary matrices
  
  

  !print *, 'A1',A1
  
 
  


  
  M0=matmul(B0,A0)
  M1=matmul(B1,A1)
  Mm1=0.d0
  call inverseLU(12,M1,Mm1)
  M2=matmul(Mm1,M0)

  !Mm0=0.d0
  !call svdinverse(12,12,M0,Mm0,5*12,info)
  !call inverseLU(12,M0,Mm0)

  !call inverse(12,M0,12,Mm0)
  !print *, "hello"
  
  !M1=matmul(B1,A1)

  !print *, 'Mm0', Mm0(2,:)
  !print *, 'identity' 
  !A0=matmul(Mm0,A0)
  !do i = 1, 12
  !   print *, A0(i,i)
  !enddo
     
  
  !print *, 'M1', M1(1,:)
  ! coefs 

  coef(1:6,1:2) = 0.d0
  
  !M2=matmul(Mm0,M1)

  !print *, 'M2',M2(7,:)
 
  !print *, M2(1,1:6)
  !stop

  coef(1:6,1)=M2(1,1:6)
  coef(1:6,2)=M2(7,7:12)

  !print *, "das ist gut"

  !print *, "ik, jk", ik, jk

  !if(info.ne.0) print *, "bad configuration=", normal,eta
end subroutine MizutaniIsoFree


subroutine MizutaniIso(coef,rho0,rho1,lam0,lam1,mu0,mu1,ik,jk,dx,dz,eta,normal)

  ! isotropic modified operators for discontinuities

  implicit none 
  integer :: i,info 
  double precision :: coef(1:6,1:2),eta(0:1,1:2),normal(1:2)
  double precision :: rho0,rho1,lam0,lam1,mu0,mu1,dx,dz
  integer :: ik,jk
  double precision :: nvx,nvy,nvx2,nvy2,nvxnvy
  double precision :: C1(1:4), C0(1:4)
  double precision :: eta0x,eta1x,eta0y,eta1y

  ! intermediate values
  double precision :: xe1ux,ye1uy,xe2ux_dx,xe2ux_dy,ye2uy_dx,ye2uy_dy
  double precision :: xe3ux_dxx,xe3ux_dyy,xe3ux_dxy,ye3uy_dxx,ye3uy_dyy,ye3uy_dxy 
  double precision :: xe4ux_dxa,xe4ux_dya,xe4uy_dxa,xe4uy_dya 
  double precision :: ye4ux_dxa,ye4ux_dya,ye4uy_dxa,ye4uy_dya
  double precision :: xe4ux_dxb,xe4ux_dyb,xe4uy_dxb,xe4uy_dyb 
  double precision :: ye4ux_dxb,ye4ux_dyb,ye4uy_dxb,ye4uy_dyb
  double precision :: xe5ux_dxxa,xe5ux_dxya,xe5uy_dxxa,xe5uy_dxya,xe5ux_dyya,xe5uy_dyya
  double precision :: ye5ux_dxxa,ye5ux_dxya,ye5uy_dxxa,ye5uy_dxya,ye5ux_dyya,ye5uy_dyya
  double precision :: xe5ux_dxxb,xe5ux_dxyb,xe5uy_dxxb,xe5uy_dxyb,xe5ux_dyyb,xe5uy_dyyb
  double precision :: ye5ux_dxxb,ye5ux_dxyb,ye5uy_dxxb,ye5uy_dxyb,ye5ux_dyyb,ye5uy_dyyb 
  double precision :: xe6ux_dxxa,xe6uy_dxya,xe6ux_dyya,ye6ux_dxya,ye6uy_dxxa,ye6uy_dyya
  double precision :: xe6ux_dxxb,xe6uy_dxyb,xe6ux_dyyb,ye6ux_dxyb,ye6uy_dxxb,ye6uy_dyyb 
 
  ! Matrices
  double precision, dimension(1:12,1:12) :: A0,B0,A1,B1,M0,Mm0,M1,M2,Bm0,Mm1
  double precision :: r, s
  
  nvx=normal(1)
  nvy=normal(2)
  nvx2=nvx*nvx
  nvy2=nvy*nvy
  nvxnvy=nvx*nvy


  eta0x = dble(ik)*eta(0,1)
  eta1x =-dble(ik)*eta(1,1)

  eta0y = dble(jk)*eta(0,2)
  eta1y =-dble(jk)*eta(1,2)  

  C0 = 0.d0
  C1 = 0.d0

  A0 = 0.d0
  A1 = 0.d0
  B0 = 0.d0
  B1 = 0.d0
  
  !print *, eta0x,eta1x,eta0y,eta1y
  !print *, nvx,nvy

  
  ! isotropic, PSV

  C0(1) = lam0 + 2.d0*mu0
  C0(2) = lam0
  C0(3) = C0(1)
  C0(4) = mu0

 
  C1(1) = lam1 + 2.d0*mu1
  C1(2) = lam1
  C1(3) = C1(1)
  C1(4) = mu1
  
  

  ! Hereafter we follow the Oleg Ovcharenko's formulation
  
  
  ! BC eq.1
  xe1ux = 1.d0
  ye1uy = 1.d0
  
  ! BC eq.2
  xe2ux_dx = nvy
  xe2ux_dy = -nvx
  ye2uy_dx = nvy
  ye2uy_dy = -nvx
  
  ! BC eq.3
  xe3ux_dxx = nvy2
  xe3ux_dyy = nvx2
  xe3ux_dxy = -2.d0*nvxnvy
  ye3uy_dxx = nvy2
  ye3uy_dyy = nvx2
  ye3uy_dxy = -2.d0*nvxnvy

  ! BC eq.4

  ! Above
  xe4ux_dxa = nvx*C1(1)
  xe4ux_dya = nvy*C1(4)
  xe4uy_dxa = nvy*C1(4)
  xe4uy_dya = nvx*C1(2)
  
  ye4ux_dxa = nvy*C1(2)
  ye4ux_dya = nvx*C1(4)
  ye4uy_dxa = nvx*C1(4)
  ye4uy_dya = nvy*C1(3)

  ! Below
  xe4ux_dxb = nvx*C0(1)
  xe4ux_dyb = nvy*C0(4)
  xe4uy_dxb = nvy*C0(4)
  xe4uy_dyb = nvx*C0(2)
  
  ye4ux_dxb = nvy*C0(2)
  ye4ux_dyb = nvx*C0(4)
  ye4uy_dxb = nvx*C0(4)
  ye4uy_dyb = nvy*C0(3)

  ! BC eq.5
  
  ! Above
  xe5ux_dxxa = nvxnvy*C1(1)
  xe5ux_dxya = nvy2*C1(4)-nvx2*C1(1)
  xe5uy_dxxa = nvy2*C1(4)
  xe5uy_dxya = nvxnvy*C1(2)-nvxnvy*C1(4)
  xe5ux_dyya = - nvxnvy*C1(4)
  xe5uy_dyya = - nvx2*C1(2)

  ye5ux_dxxa = nvy2*C1(2)
  ye5ux_dxya = nvxnvy*C1(4)-nvxnvy*C1(2)
  ye5uy_dxxa = nvxnvy*C1(4)
  ye5uy_dxya = nvy2*C1(3)-nvx2*C1(4)
  ye5ux_dyya =  - nvx2*C1(4)
  ye5uy_dyya =  - nvxnvy*C1(3)
  
  ! Below
  xe5ux_dxxb = nvxnvy*C0(1)
  xe5ux_dxyb = nvy2*C0(4)-nvx2*C0(1)
  xe5uy_dxxb = nvy2*C0(4)
  xe5uy_dxyb = nvxnvy*C0(2)-nvxnvy*C0(4)
  xe5ux_dyyb = - nvxnvy*C0(4)
  xe5uy_dyyb = - nvx2*C0(2)
  
  ye5ux_dxxb = nvy2*C0(2)
  ye5ux_dxyb = nvxnvy*C0(4)-nvxnvy*C0(2)
  ye5uy_dxxb = nvxnvy*C0(4)
  ye5uy_dxyb = nvy2*C0(3)-nvx2*C0(4)
  ye5ux_dyyb =  - nvx2*C0(4)
  ye5uy_dyyb =  - nvxnvy*C0(3)

  ! BC eq.6
  ! Above C1
  xe6ux_dxxa = C1(1)/rho1
  xe6uy_dxya = (C1(2)+C1(4))/rho1
  xe6ux_dyya = C1(4)/rho1
  
  ye6ux_dxya = (C1(4)+C1(2))/rho1
  ye6uy_dxxa = C1(4)/rho1
  ye6uy_dyya = C1(3)/rho1

  ! Below C0
  xe6ux_dxxb = C0(1)/rho0
  xe6uy_dxyb = (C0(2)+C0(4))/rho0
  xe6ux_dyyb = C0(4)/rho0
  
  ye6ux_dxyb = (C0(4)+C0(2))/rho0
  ye6uy_dxxb = C0(4)/rho0
  ye6uy_dyyb = C0(3)/rho0
    
    
  ! Full matrices for left and right points
  
  ! B0

  B0 = 0.d0
  
  B0( 1, 1) = xe1ux

  B0( 2, 2) = xe2ux_dx
  B0( 2, 3) = xe2ux_dy
 
  B0( 3, 4) = xe3ux_dxx
  B0( 3, 5) = xe3ux_dyy
  B0( 3, 6) = xe3ux_dxy
  
  B0( 4, 2) = xe4ux_dxb
  B0( 4, 3) = xe4ux_dyb 
  B0( 4, 8) = xe4uy_dxb
  B0( 4, 9) = xe4uy_dyb

  B0( 5, 4) = xe5ux_dxxb
  B0( 5, 5) = xe5ux_dyyb
  B0( 5, 6) = xe5ux_dxyb
  B0( 5,10) = xe5uy_dxxb
  B0( 5,11) = xe5uy_dyyb
  B0( 5,12) = xe5uy_dxyb
  
  B0( 6, 4) = xe6ux_dxxb
  B0( 6, 5) = xe6ux_dyyb
  B0( 6,12) = xe6uy_dxyb
  
  B0( 7, 7) = ye1uy

  B0( 8, 8) = ye2uy_dx
  B0( 8, 9) = ye2uy_dy

  B0( 9,10) = ye3uy_dxx
  B0( 9,11) = ye3uy_dyy
  B0( 9,12) = ye3uy_dxy

  B0(10, 2) = ye4ux_dxb
  B0(10, 3) = ye4ux_dyb
  B0(10, 8) = ye4uy_dxb
  B0(10, 9) = ye4uy_dyb

  B0(11, 4) = ye5ux_dxxb 
  B0(11, 5) = ye5ux_dyyb
  B0(11, 6) = ye5ux_dxyb
  B0(11,10) = ye5uy_dxxb
  B0(11,11) = ye5uy_dyyb
  B0(11,12) = ye5uy_dxyb
  
  B0(12, 6) = ye6ux_dxyb
  B0(12,10) = ye6uy_dxxb
  B0(12,11) = ye6uy_dyyb
  

  ! A0

  r = eta0x
  s = eta0y

  A0 = 0.d0
  
  A0( 1, 1) = 1.d0
  A0( 1, 2) = r*dx
  A0( 1, 3) = s*dz
  A0( 1, 4) = (r*dx)*(r*dx)*5.d-1
  A0( 1, 5) = (s*dz)*(s*dz)*5.d-1
  A0( 1, 6) = r*s*dx*dz

  A0( 2, 2) = 1.d0
  A0( 2, 4) = r*dx
  A0( 2, 6) = s*dz

  A0( 3, 3) = 1.d0
  A0( 3, 5) = s*dz
  A0( 3, 6) = r*dx
  
  A0( 4, 4) = 1.d0
  A0( 5, 5) = 1.d0
  A0( 6, 6) = 1.d0
  
  A0(7:12,7:12) = A0(1:6,1:6)

         

  ! B1

  B1 = 0.d0
  
  B1( 1, 1) = xe1ux

  B1( 2, 2) = xe2ux_dx
  B1( 2, 3) = xe2ux_dy
 
  B1( 3, 4) = xe3ux_dxx
  B1( 3, 5) = xe3ux_dyy
  B1( 3, 6) = xe3ux_dxy
  
  B1( 4, 2) = xe4ux_dxa
  B1( 4, 3) = xe4ux_dya 
  B1( 4, 8) = xe4uy_dxa
  B1( 4, 9) = xe4uy_dya

  B1( 5, 4) = xe5ux_dxxa
  B1( 5, 5) = xe5ux_dyya
  B1( 5, 6) = xe5ux_dxya
  B1( 5,10) = xe5uy_dxxa
  B1( 5,11) = xe5uy_dyya
  B1( 5,12) = xe5uy_dxya
  
  B1( 6, 4) = xe6ux_dxxa
  B1( 6, 5) = xe6ux_dyya
  B1( 6,12) = xe6uy_dxya
  
  B1( 7, 7) = ye1uy

  B1( 8, 8) = ye2uy_dx
  B1( 8, 9) = ye2uy_dy

  B1( 9,10) = ye3uy_dxx
  B1( 9,11) = ye3uy_dyy
  B1( 9,12) = ye3uy_dxy

  B1(10, 2) = ye4ux_dxa
  B1(10, 3) = ye4ux_dya
  B1(10, 8) = ye4uy_dxa
  B1(10, 9) = ye4uy_dya

  B1(11, 4) = ye5ux_dxxa 
  B1(11, 5) = ye5ux_dyya
  B1(11, 6) = ye5ux_dxya
  B1(11,10) = ye5uy_dxxa
  B1(11,11) = ye5uy_dyya
  B1(11,12) = ye5uy_dxya
  
  B1(12, 6) = ye6ux_dxya
  B1(12,10) = ye6uy_dxxa
  B1(12,11) = ye6uy_dyya


      


  ! A1


  r = eta1x;
  s = eta1y;
  


  A1 = 0.d0
  
  A1( 1, 1) = 1.d0
  A1( 1, 2) = r*dx
  A1( 1, 3) = s*dz
  A1( 1, 4) = (r*dx)*(r*dx)*5.d-1
  A1( 1, 5) = (s*dz)*(s*dz)*5.d-1
  A1( 1, 6) = r*s*dx*dz
  
  A1( 2, 2) = 1.d0
  A1( 2, 4) = r*dx
  A1( 2, 6) = s*dz
  
  A1( 3, 3) = 1.d0
  A1( 3, 5) = s*dz
  A1( 3, 6) = r*dx
  
  A1( 4, 4) = 1.d0
  A1( 5, 5) = 1.d0
  A1( 6, 6) = 1.d0
  
  A1(7:12,7:12) = A1(1:6,1:6)
  

  

  ! temporary matrices
  
  

  !print *, 'A1',A1
  
 
  


  
  M0=matmul(B0,A0)
  M1=matmul(B1,A1)
  Mm1=0.d0
  call inverseLU(12,M1,Mm1)
  M2=matmul(Mm1,M0)

  !Mm0=0.d0
  !call svdinverse(12,12,M0,Mm0,5*12,info)
  !call inverseLU(12,M0,Mm0)

  !call inverse(12,M0,12,Mm0)
  !print *, "hello"
  
  !M1=matmul(B1,A1)

  !print *, 'Mm0', Mm0(2,:)
  !print *, 'identity' 
  !A0=matmul(Mm0,A0)
  !do i = 1, 12
  !   print *, A0(i,i)
  !enddo
     
  
  !print *, 'M1', M1(1,:)
  ! coefs 

  coef(1:6,1:2) = 0.d0
  
  !M2=matmul(Mm0,M1)

  !print *, 'M2',M2(7,:)
 
  !print *, M2(1,1:6)
  !stop

  coef(1:6,1)=M2(1,1:6)
  coef(1:6,2)=M2(7,7:12)

  !print *, "das ist gut"

  !print *, "ik, jk", ik, jk

  !if(info.ne.0) print *, "bad configuration=", normal,eta
end subroutine MizutaniIso


