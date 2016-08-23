program generator
 implicit none

 integer, parameter     :: prec = kind (1.0) 
 integer                :: i,j
 integer                :: level,tmp

 integer, parameter :: times = 1

 integer, parameter                             :: nbre_layers = 3
 integer, dimension(nbre_layers)                :: layer_thikness
 real(kind=prec), dimension(nbre_layers)        :: vp
 real(kind=prec), dimension(nbre_layers)        :: vs
 real(kind=prec), allocatable                   :: fullvp (:,:)
 real(kind=prec), allocatable                   :: fullvs (:,:)
 real(kind=prec), allocatable                   :: fullrho(:,:),tmpM(:,:)
 real(kind=prec)                                :: dval,pval,xover1,xover2,grad,const
 real(kind=prec), parameter                     :: vpwater = 1.5e3
 integer                                        :: NZ, NX
 integer ::ix,iz,tmpint
 double precision :: boundary,dx,tmpfloat
 integer                                        :: NX_TOTAL, NZ_TOTAL
 integer                                        :: recl_size
!************************************************************************	

! PROGRAM GENERATES VELOCITY MODEL (VP & VS)

!************************************************************************
! Model initialization	
 NX_TOTAL = (601-1)*times+1
 NZ_TOTAL = (161-1)*times+1

!*** Thikness, Vp, Vs *******************************************		

 layer_thikness(1)=71;vp(1)=1500;vs(1)=1000;
 layer_thikness(2)=91;vp(2)=3000;vs(2)=1700;
 layer_thikness(3)=162;vp(3)=2000;vs(3)=1170;
! layer_thikness(4)=180;vp(4)=2700;vs(4)=1700;
! layer_thikness(5)=200;vp(5)=3000;vs(5)=1900;

! layer_thikness(1)=40;vp(1)=2200;vs(1)=1400;
! layer_thikness(2)=90;vp(2)=2200;vs(2)=1400;
! layer_thikness(3)=140;vp(3)=2200;vs(3)=1400;
! layer_thikness(4)=180;vp(4)=2200;vs(4)=1400;
! layer_thikness(5)=200;vp(5)=2200;vs(5)=1400;
 
 

 
 



 allocate (fullvp (NX_TOTAL, NZ_TOTAL) )
 allocate (fullvs (NX_TOTAL, NZ_TOTAL) )
 allocate (fullrho(NX_TOTAL, NZ_TOTAL) )
 !allocate(tmpM(NZ_TOTAL,NX_TOTAL))
 if(0.eq.1) then
 open(1,file='./models/model_cp_400x200.txt',form='formatted')
 open(2,file='./models/model_cs_400x200.txt',form='formatted')
 open(3,file='./models/model_rho_400x200.txt',form='formatted')
 read(1,*) tmpM
 tmpM=transpose(tmpM)
 fullvp(1:NX_TOTAL,1:NZ_TOTAL)=tmpM(1:NX_TOTAL,NZ_TOTAL:1:-1)
 read(2,*) tmpM
 tmpM=transpose(tmpM)
 fullvs(1:NX_TOTAL,1:NZ_TOTAL)=tmpM(1:NX_TOTAL,NZ_TOTAL:1:-1)
 read(3,*) tmpM
 tmpM=transpose(tmpM)
 fullrho(1:NX_TOTAL,1:NZ_TOTAL)=tmpM(1:NX_TOTAL,NZ_TOTAL:1:-1)

 !write(14,*)fullvp
 do ix=1,NX_TOTAL
    do iz=1,NZ_TOTAL
       if(fullvs(ix,iz)<400.d0) then
          fullvs(ix,iz) = fullvp(ix,iz)/1.7d0
       endif
    enddo
 enddo
 endif


 if(0.eq.1) then

 dx = 2.d-2

 do ix = 1,NX_TOTAL
    boundary = 199.d0*dx-dx*dble(ix-1)*199.d0/399.d0
    tmpint =nint(boundary/dx)
    fullvp(ix,1:tmpint) = 2200.e0
    fullvs(ix,1:tmpint) = 1400.e0
    fullvp(ix,tmpint+1:NZ_TOTAL) = 2500.e0
    fullvs(ix,tmpint+1:NZ_TOTAL) = 2100.e0
 enddo

 endif


 if(1.eq.0) then

 dx = 2.d-2

 do ix = 1,NX_TOTAL
    boundary = 199.d0*dx*dble(times)-dx*dble(ix-1)*199.d0/399.d0
    tmpint =nint(boundary/dx)
    !fullvp(ix,1:tmpint) = 2290.e0
    !fullvs(ix,1:tmpint) = 3000.e0
  ! 
    !fullvp(ix,tmpint+1:NZ_TOTAL) = 2200.e0
    !fullvs(ix,tmpint+1:NZ_TOTAL) = 1400.e0
    ! giza inversed (above)
    
    fullvp(ix,1:tmpint)=2200.e0
    fullvs(ix,1:tmpint)=1400.e0
    fullrho(ix,1:tmpint)=2119.e0
    
    fullvp(ix,tmpint+1:NZ_TOTAL)=3000.e0
    fullvs(ix,tmpint+1:NZ_TOTAL)=1900.e0
    fullrho(ix,tmpint+1:NZ_TOTAL)=2290.e0

! below
!rho2290
!vp 3000
!vs 1900
!for above:
!rho 2119
!Vp 2200
!Vs 1400



   
    
 enddo
 
 write(14,*)fullvp(1:400,1:200)

 endif




 if(0.eq.1) then
    dx=2.d-2
    fullvp(:,:) = 2200.e0
    fullvs(:,:) = 1400.e0
 
    do iz=50,150
       do ix = 100,200
          tmpfloat=(dble(ix)-150.d0)**2+(dble(iz)-100.d0)**2
          if(tmpfloat.le.2500.d0) then
             fullvp(ix,iz)=3000.e0
             fullvs(ix,iz)=1900.e0

          endif
       enddo
    enddo
     write(14,*)fullvp(1:400,1:200)

 endif


 recl_size = prec * NX_TOTAL * NZ_TOTAL

!****************************************************************	   
 
 open (1,file='../models/2d_guy.vp',form='unformatted',access='direct',recl=recl_size)
 open (2,file='../models/2d_guy.vs',form='unformatted',access='direct',recl=recl_size) 
 open (3,file='../models/2d_guy.rho',form='unformatted',access='direct',recl=recl_size)

 if(0.eq.0) then
    tmp=1
    do j=1, NZ_TOTAL
       level=tmp
       if (j .gt. layer_thikness(level)) level=level+1
       tmp=level
       do i=1, NX_TOTAL
          fullvp(i,j) = vp(level)
          fullvs(i,j) = vs(level)
       enddo
    enddo
 endif
 write(1,rec=1) fullvp(:,:)
 write(2,rec=1) fullvs(:,:)

 if(0.eq.0) then
 do j=1,NZ_TOTAL
     do i=1,NX_TOTAL
        pval = fullvp(i,j)
        if (pval.le.vpwater) then
           dval = 1000.E0
        elseif (pval > vpwater .and. pval < 2000.E0) then
           dval = 2351.E0-(7497.E0)*(pval/1000.E0)**(-4.656E0)
        elseif (pval >= 2000.E0 .and. pval <= 2150.E0) then
           xover1 = 2351.E0-(7497.E0)*(2000.E0/1000.E0)**(-4.656E0)
           xover2 = 1740.E0*(2150.E0/1000.E0)**(0.25E0)
           grad = 150.E0/(xover2-xover1)
           const = 2000.E0-(xover1*grad)
           dval = (pval-const)/grad
        elseif (pval > 2150) then
           dval = 1740.E0*(pval/1000.E0)**(0.25E0)
        endif
        fullrho(i,j) = dval
     enddo
  enddo

  endif


  write(3,rec=1) fullrho(:,:)

  

 close (1,status='keep')
 close (2,status='keep')
 close (3,status='keep')
end program generator
