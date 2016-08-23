subroutine disconConfig
  use parameters
  implicit none

  
  ! Discontinuity configuration


  ! diagonal discontinuity
  nDiscon = 0


  if(0.eq.0) then
  
  nDiscon = 1
  lengthDiscon = 40*nx+1
  
  if(nDiscon.ne.0) then
     allocate(dscr(1:2,1:lengthDiscon,1:nDiscon))
     do ix =1,lengthDiscon
        dscr(1,ix,1) = dble(ix-1)*dx/40.d0
        dscr(2,ix,1) = 199.d0*dble(times)*dz-dscr(1,ix,1)*199.d0/399.d0
     enddo
  endif
  markers(1:maxnz,1:maxnz) = 0
  markers(1:nx+1,1:nz+1) = 1 ! for the moment NF will search for all the points (of course it is not good)

  markers(1:maxnz,1:maxnz) = 0
  do ix = 1,nx+1
     tmpint=nint(199.d0-dble(ix-1)*199.d0/399.d0)
     if(tmpint-3.ge.1) then
        markers(ix,tmpint-3:tmpint+3) = 1
     else if(tmpint-2.ge.1) then
        markers(ix,tmpint-2:tmpint+3) = 1
     else
        markers(ix,1:tmpint+3) = 1
     endif
  enddo
 
  endif

  ! Circle discontinuity
  
  if(1.eq.0) then
     
     nDiscon = 2
     lengthDiscon = 40*100+1
     
     if(nDiscon.ne.0) then
        allocate(dscr(1:2,1:lengthDiscon,1:nDiscon))
        do ix =1,lengthDiscon
           dscr(1,ix,1) = 99.d0*dx+dble(ix-1)*dx/40.d0
           dscr(1,ix,2) = dscr(1,ix,1)
           dscr(2,ix,1) = 99.d0+sqrt(2500.d0*dx**2-(dble(ix-1)*dx/40.d0-dx*149.d0)**2)
           !dscr(2,ix,1) = 99.d0
           dscr(2,ix,2) = 99.d0-sqrt(2500.d0*dx**2-(dble(ix-1)*dx/40.d0-dx*149.d0)**2)
        enddo
     endif
     markers(1:nx+1,1:nz+1)=1
  endif


  ! Oleg discontinuties

  if(0.eq.1) then
     nDiscon=1
     lengthDiscon = 160000
     allocate(dscr(1:2,1:lengthDiscon,1:nDiscon))
     open(1, file='grid_x_Oleg.txt')
     do ix = 1,400
        read(1,*)  dscr(1,400*(ix-1)+1:400*ix,1)
     enddo
!123  format(400(F.2,1x))
     close(1)
     open(1, file='grid_y_Oleg.txt')
     do ix = 1,400
        read(1,*)  dscr(2,400*(ix-1)+1:400*ix,1)
     enddo
     close(1)
     


  endif
  

end subroutine disconConfig


subroutine freeConfig
  
  use parameters
  implicit none

  

  ! Free surface configuration

  lengthFreeSurface = 0
  if(lengthFreeSurface.ne.0) then
     allocate(free(1:2,1:lengthFreeSurface))
     do ix=1,lengthFreeSurface
        free(1,ix) = dble(ix-1)*dx/40.d0
        free(2,ix) = 5.d-1+3.d-1*sin(free(1,ix)/(dble(nx)*dx)*pi*4)
     enddo
     zerodisplacement(1:maxnz,1:maxnz)=0
     do ix=1,nx+1
        tmpint=nint((5.d-1+3.d-1*sin((dble(ix-1)*dx)/(dble(nx)*dx)*pi*4))/dx)
        zerodisplacement(ix,1:tmpint) = 1
        vp(ix,1:tmpint)=0.d0
        vs(ix,1:tmpint)=0.d0
     enddo

     
  endif

end subroutine freeConfig
