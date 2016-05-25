program opt1d
  
  ! Computation of the synthetic seismograms in the time domain
  !using the normal operators.
  !
  !					originally from	1995.10  N.Takeuchi
  !                                                     2016. 5  N.Fuji
  
  implicit none
  integer maxnz
  real*8 pi
  integer, parameter :: maxnz = 10000 
  double precision, parameter :: pi=3.1415926535897932d0 
  
  !  parameters for the gridding
  integer nt,nz,it,ist,isz
  real*8 dt,dz
  ! parameter for the structure
  double precision rrho(4),kkappa(4)
  double precision rho(maxnz+1),kappa(maxnz+1)
  ! parameter for the wavefield
  double precision t
  double precision u(maxnz+1),u1(maxnz+1),u2(maxnz+1),f(maxnz+1)
  double precision e1(maxnz+1),e2(maxnz+1),e3(maxnz+1)
  double precision f1(maxnz+1),f2(maxnz+1),f3(maxnz+1)
  ! parameter for the source
  double precision tp,ts
  ! parameter for the receiver
  integer nr
  
  ! reading the parameter files
  call pinput( maxnz,nt,nz,dt,dz,rho,kappa,tp,ts,nr )
  !Initializing the data
  call datainit( maxnz,u )
  call datainit( maxnz,u1 )
  call datainit( maxnz,u2 )
  call datainit( maxnz,f )
  call datainit( maxnz,rho )
  call datainit( maxnz,kappa )
  ! computing the intermediate parameters
  call calstruct( rrho,dz,nz,rho )
  call calstruct( kkappa,dz,nz,kappa )
  call cales( nz,rho,kappa,dt,dz,e1,e2,e3,f1,f2,f3 )
  ist = dnint( 2 * tp / dt )
  isz = nz / 2 + 1
  
  t = 0.d0
  write(6,*) real(t),real(u(nr))
  do it=0,nt
     call calf( it,t,ist,isz,dt,dz,rho(isz),tp,ts,f )
     ! evaluating the next step
     call calstep( nz,e1,e2,e3,f1,f2,f3,u,u1,u2,f )
     ! increment of t
     t = t + dt
     if ( mod(it,10).eq.9 ) write(6,*) real(t),real(u(nr))
     write(6,*) real(t),real(u(nr))
  enddo
  
end program opt1d



subroutine pinput (maxnz,nt,nz,dt,dz,rho,kappa,tp,ts,nr )

  implicit none
  double precision f0,t0
  integer maxnz,nt,nx,nz,nrx,nrz
  double precision dt,dx,dz !rho(*),lam(*),mu(*),tp,ts
  double precision rho(maxnz+1),kappa(maxnz+1),vp(maxnz+1),tmpvector(maxnz+1)
  character*80 tmpfile,dummy
  character*80 vpfile,vsfile,rhofile
  tmpfile='tmpfileforwork'
  
  
  ! temporary file open
  open( unit=11, file=tmpfile, status='unknown' )
  !writing to the temporary file
100 continue
  read(5,110) dummy
110 format(a80)
  if ( dummy(1:1).eq.'c' ) goto 100
  if ( dummy(1:3).eq.'end' ) goto 120
  write(11,110) dummy
  goto 100
120 continue
  ! temporary file close
  close(11)
  
  ! temporary file open
  open( unit=11, file=tmpfile, status='unknown' )
  ! reading the parameter
  read(11,*) nt,nx,nz
  if ( nx.gt.maxnz ) pause 'nx is too large (pinput).'
  if ( nz.gt.maxnz ) pause 'nz is too large (pinput).'
  read(11,*) dt,dx,dz
111 format(a80)
  read(11,111) vpfile
  read(11,111) vsfile
  read(11,111) rhofile
  read(11,*) f0,t0
  read(11,*) nrx,nrz
  ! temporary file close
  close(11)
  !


  call calstruct2D1D(maxnz,rhofile,dx,dz,nx,nz,rho)
  call calstruct2D1D(maxnz,vpfile,dx,dz,nx,nz,vp)
  call calstruct2D1D(maxnz,vsfile,dx,dz,nx,nz,tmpvector)
  call calstruct2(maxnz,1,nz,rho,vp,vs,lam,mu)
  return
end subroutine pinput


subroutine calstruct2(maxnz,nx,nz,rho,vp,vs,lam,mu)
  implicit none
  
  integer i,j,maxnz,nx,nz
  double precision rho(maxnz+1,maxnz+1),vp(maxnz+1,maxnz+1),vs(maxnz+1,maxnz+1)
  double precision lam(maxnz+1,maxnz+1),mu(maxnz+1,maxnz+1)

  do i=1,nx+1
     do j=1,nz+1
        mu(i,j)=rho(i,j)*vs(i,j)*vs(i,j)
        lam(i,j)=rho(i,j)*vp(i,j)*vp(i,j)-2*mu(i,j)
     enddo
  enddo
end subroutine calstruct2

subroutine calstruct2D1D( maxnz,file2d,dx,dz,nx,nz,rho )
  implicit none
  integer maxnz,nx,nz
  double precision dx,dz,rho(maxnz+1)
  real(kind(1.0)) rrho(maxnz+1,maxnz+1)
  integer i,j,k,nox(6),noz(6)
  double precision x,z,xmax,zmax,trho,coef1,coef2
  integer recl_size
  character*80 file2d
  recl_size=kind(1.0)*(nx+1)*(nz+1)

  open (1,file=file2d,form='unformatted',access='direct',recl=recl_size)
  read(1,rec=1) rrho(1:nx+1,1:nz+1)
  close(1)
  rho(1:nz+1)=1.d-3*rrho(1,1:nz+1)
  
  return
end subroutine calstruct2D1D


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine pinput_old( maxnz,nt,nz,dt,dz,rho,kappa,tp,ts,nr )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	integer maxnz,nt,nz,nr
	double precision dt,dz,rho(4),kappa(4),tp,ts
	character*80 tmpfile,dummy
c
	data tmpfile / '/tmp/work' /
c
c temporary file open
	open( unit=11, file=tmpfile, status='unknown' )
c writing to the temporary file
  100	continue
	  read(5,110) dummy
  110	  format(a80)
	  if ( dummy(1:1).eq.'c' ) goto 100
	  if ( dummy(1:3).eq.'end' ) goto 120
	  write(11,*) dummy
	  goto 100
  120	continue
c temporary file close
	close(11)
c 
c temporary file open
	open( unit=11, file=tmpfile, status='unknown' )
c reading the parameter
	read(11,*) nt,nz
	if ( nz.gt.maxnz ) pause 'nz is too large (pinput).'
	read(11,*) dt,dz
	read(11,*) rho(1),rho(2),rho(3),rho(4)
	read(11,*) kappa(1),kappa(2),kappa(3),kappa(4)
	read(11,*) tp,ts
	read(11,*) nr
c temporary file close
	close(11)
c
	return
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine datainit( maxnz,u )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	integer maxnz
	double precision u(*)
	integer i
c
	do 100 i=1,maxnz+1
	  u(i) = 0.d0
  100	continue
c
	return
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine calstruct( rrho,dz,nz,rho )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	integer nz
	double precision rrho(4),dz,rho(*)
	integer i,j
	double precision r,rmax,coef,trho
c
	rmax = dz * nz
	do 100 i=1,nz+1
	  r = dble(i-1) * dz
          trho = 0.d0
	  do 110 j=1,4
	    if ( j.eq.1 ) then
	      coef = 1.d0
	    else
	      coef = coef * ( r / rmax )
	    endif
	    trho = trho + rrho(j) * coef
  110	  continue
          rho(i) = trho
  100	continue
c
	return
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine cales( nz,rho,kappa,dt,dz,e1,e2,e3,f1,f2,f3 )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	integer nz
	double precision rho(*),kappa(*),dt,dz,e1(*),e2(*),e3(*),f1(*),f2(*),f3(*)
	integer iz
	double precision dt2,dz2
c
	dt2 = dt * dt
	dz2 = dz * dz
c
	e1(1) = 0.d0
	e2(1) = 2.d0 - ( kappa(1) + kappa(2) )
     &	               / rho(1) * dt2 / dz2
	e3(1) = ( kappa(1) + kappa(2) ) / rho(1) * dt2 / dz2
	f1(1) = 0.d0
	f2(1) = e2(1) / 12.d0
	f3(1) = ( e3(1) - 2.d0 ) / 12.d0
	do 100 iz=2,nz
	  e1(iz) = ( kappa(iz-1) + kappa(iz) )
     &	           / ( 2.d0 * rho(iz) ) * dt2 / dz2
	  e2(iz) = 2.d0
     &	           - ( kappa(iz-1) + 2.d0 * kappa(iz) + kappa(iz+1) )
     &	             / ( 2.d0 * rho(iz) ) * dt2 / dz2
	  e3(iz) = ( kappa(iz) + kappa(iz+1) )
     &	           / ( 2.d0 * rho(iz) ) * dt2 / dz2
	  f1(iz) = ( e1(iz) - 1.d0 ) / 12.d0
	  f2(iz) = e2(iz) / 12.d0
	  f3(iz) = ( e3(iz) - 1.d0 ) / 12.d0
  100	continue
	e1(nz+1) = ( kappa(nz) + kappa(nz+1) )
     &	           / rho(nz+1) * dt2 / dz2
	e2(nz+1) = 2.d0 - ( kappa(nz) + kappa(nz+1) )
     &	                  / rho(nz+1) * dt2 / dz2
	e3(nz+1) = 0.d0
	f1(nz+1) = ( e1(nz+1) - 2.d0 ) / 12.d0
	f2(nz+1) = e2(nz+1) / 12.d0
	f3(nz+1) = 0.d0
c
	return
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine calf( it,t,ist,isz,dt,dz,rho,tp,ts,f )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	double precision pi
	parameter ( pi=3.1415926535897932d0 )
	integer it,ist,isz
	double precision t,dt,dz,rho,tp,ts,f(*)
	double precision b
c
	if ( it.le.ist ) then
	  b = pi * ( t - ts ) / tp
	  f(isz) = dsqrt(pi) / 2.d0 * (b*b-0.5d0) * dexp(-b*b) / dz
	  if ( (it.eq.0).or.(it.eq.ist) ) f(isz) = f(isz) / 2.d0
	  f(isz) = f(isz) * dt * dt / rho
	else
	  f(isz) = 0.d0
	endif
c
	return
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine calstep( nz,e1,e2,e3,f1,f2,f3,u,u1,u2,f )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	integer nz
	double precision e1(*),e2(*),e3(*),f1(*),f2(*),f3(*)
	double precision u(*),u1(*),u2(*),f(*)
	integer iz
	double precision tmp1,tmp2,tmp3
c
c evaluating the u using the unmodified operators
	u(1) = - u2(1)
     &	       + e2(1) * u1(1)
     &	       + e3(1) * u1(2)
     &	       + f(1)
	do 100 iz=2,nz
	  u(iz) = - u2(iz)
     &	          + e1(iz) * u1(iz-1)
     &	          + e2(iz) * u1(iz)
     &	          + e3(iz) * u1(iz+1)
     &	          + f(iz)
  100	continue
	u(nz+1) = - u2(nz+1)
     &	          + e1(nz+1) * u1(nz)
     &	          + e2(nz+1) * u1(nz+1)
     &	          + f(nz+1)
c computing du using 1-D Born approximation
	tmp2 = ( u(1) - u1(1) - u1(1) + u2(1) )
	do 110 iz=1,nz
	  tmp3 = ( u(iz+1) - u1(iz+1) - u1(iz+1) + u2(iz+1) )
	  u(iz) = u(iz)
     &	          + f1(iz) * tmp1
     &	          + f2(iz) * tmp2
     &	          + f3(iz) * tmp3
	  tmp1 = tmp2
	  tmp2 = tmp3
  110	continue
	u(nz+1) = u(nz+1)
     &	          + f1(nz+1) * tmp1
     &	          + f2(nz+1) * tmp2
c swapping u1 & u2
	do 120 iz=1,nz+1
	  u2(iz) = u1(iz)
	  u1(iz) = u(iz)
  120	continue
c
	return
	end