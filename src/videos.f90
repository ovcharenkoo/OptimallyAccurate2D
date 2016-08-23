subroutine create_color_image(image_data_2D,NX,NY,it,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
     NPOINTS_PML,USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX,field_number)
	
  implicit none
  
  !       non linear display to enhance small amplitudes for graphics
  double precision, parameter :: POWER_DISPLAY = 1.d0
  
  !       amplitude threshold above which we draw the color point
  double precision, parameter :: cutvect = 0.01d0
  
  !       use black or white background for points that are below the threshold
  logical, parameter :: WHITE_BACKGROUND = .true.
  
  !       size of cross and square in pixels drawn to represent the source and the receivers
  integer, parameter :: width_cross=5,thickness_cross=1
  integer, parameter :: size_square=3
  
  integer NX,NY,it,field_number,ISOURCE,JSOURCE,NPOINTS_PML,nrec
  logical USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX
  
  double precision, dimension(NX,NY) :: image_data_2D
  
  integer, dimension(nrec) :: ix_rec,iy_rec
  logical, parameter :: kurama =.true.
  integer :: ix,iy,irec
  
  character(len=100) :: file_name,system_command1,system_command2,system_command3
	
  integer :: R, G, B
  double precision :: R1, G1, B1
  
  double precision :: normalized_value,max_amplitude
  
  !       open image file and create system command to convert image to more convenient format
  !       use the "convert" command from ImageMagick http://www.imagemagick.org
  if(field_number == 1) then
     write(file_name,"('image',i6.6,'_Ux.pnm')") it
     write(system_command1, "('convert image',i6.6,'_Ux.pnm snapshots/imageUx',i6.6,'.png')") it,it
     write(system_command2, "('rm image',i6.6,'_Ux.pnm')") it
  else if(field_number == 2) then
     write(file_name,"('image',i6.6,'_Uz.pnm')") it
     write(system_command1,"('convert image',i6.6,'_Uz.pnm snapshots/imageUz',i6.6,'.png')") it,it
     write(system_command2,"('rm image',i6.6,'_Uz.pnm')") it
  endif
  
  open(unit=27, file=file_name, status='unknown')
  
  write(27,"('P3')")	! write image in PNM P3 format
  
  write(27,*) NX,NY	! write image size
  write(27,*) '255'	! maximum value of each pixel color
	
  !       compute maximum amplitude
  max_amplitude = maxval(abs(image_data_2D))
  
  !       image starts in upper-left corner in PNM format
  do iy=1,NY,1
     do ix=1,NX
        
        !       define data as vector component normalized to [-1:1] and rounded to nearest integer
        !       keeping in mind that amplitude can be negative
        normalized_value = image_data_2D(ix,iy) / max_amplitude *1.5d0
        
        !       suppress values that are outside [-1:+1] to avoid small edge effects
        if(normalized_value < -1.d0) normalized_value = -1.d0
        if(normalized_value > 1.d0) normalized_value = 1.d0
        
        !       draw an orange cross to represent the source
        if((ix >= ISOURCE - width_cross .and. ix <= ISOURCE + width_cross .and. &
             iy >= JSOURCE - thickness_cross .and. iy <= JSOURCE + thickness_cross) .or. &
             (ix >= ISOURCE - thickness_cross .and. ix <= ISOURCE + thickness_cross .and. &
             iy >= JSOURCE - width_cross .and. iy <= JSOURCE + width_cross)) then
           R = 255
           G = 157
           B = 0
           
           !       display two-pixel-thick black frame around the image
        else if(ix <= 2 .or. ix >= NX-1 .or. iy <= 2 .or. iy >= NY-1) then
           R = 0
           G = 0
           B = 0
           
           !       display edges of the PML layers
        else if((USE_PML_XMIN .and. ix == NPOINTS_PML) .or. &
             (USE_PML_XMAX .and. ix == NX - NPOINTS_PML) .or. &
             (USE_PML_YMIN .and. iy == NPOINTS_PML) .or. &
             (USE_PML_YMAX .and. iy == NY - NPOINTS_PML)) then
           R = 255
           G = 150
           B = 0
           
           !       suppress all the values that are below the threshold
        else if(abs(image_data_2D(ix,iy)) <= max_amplitude * cutvect) then
           
           !       use a black or white background for points that are below the threshold
           if(WHITE_BACKGROUND) then
              R = 255
              G = 255
              B = 255
           else
              R = 0
              G = 0
              B = 0
           endif
           
           !       represent regular image points using red if value is positive, blue if negative
        elseif(kurama) then
           normalized_value=normalized_value**POWER_DISPLAY
           call plotcolor(normalized_value,R1,G1,B1)
           R = nint(R1)
           G = nint(G1)
           B = nint(B1)
        else
           if(normalized_value >= 0.d0) then
              R = 255
              G = nint(255.d0-255.d0*normalized_value**POWER_DISPLAY)
              B = nint(255.d0-255.d0*normalized_value**POWER_DISPLAY)
           else
              R = nint(255.d0-255.d0*abs(normalized_value)**POWER_DISPLAY)
              G = nint(255.d0-255.d0*abs(normalized_value)**POWER_DISPLAY)
              B = 255
           endif
        endif
        
        !       draw a green square to represent the receivers
        do irec = 1,nrec
           if((ix >= ix_rec(irec) - size_square .and. ix <= ix_rec(irec) + size_square .and. &
                iy >= iy_rec(irec) - size_square .and. iy <= iy_rec(irec) + size_square) .or. &
                (ix >= ix_rec(irec) - size_square .and. ix <= ix_rec(irec) + size_square .and. &
                iy >= iy_rec(irec) - size_square .and. iy <= iy_rec(irec) + size_square)) then
              !       use dark green color
              R = 30
              G = 180
              B = 60
           endif
        enddo
        
        !       write color pixel
        write(27,"(i3,' ',i3,' ',i3)") R,G,B
        
     enddo
  enddo
  
  !       close file
close(27)

!       call the system to convert image to GIF (can be commented out if "call system" is missing in your compiler)
call system(system_command1)
call system(system_command2)
end subroutine create_color_image


subroutine create_color_kernel(image_data_2D,NX,NY,it,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, field_number, max_amplitude,folder)
	
  implicit none
  
  character(len=100) :: folder

  !       non linear display to enhance small amplitudes for graphics
  double precision, parameter :: POWER_DISPLAY = 1.d0
  
  !       amplitude threshold above which we draw the color point
  double precision, parameter :: cutvect = 0.01d0
  
  !       use black or white background for points that are below the threshold
  logical, parameter :: WHITE_BACKGROUND = .true.
  
  !       size of cross and square in pixels drawn to represent the source and the receivers
  integer, parameter :: width_cross=5,thickness_cross=1
  integer, parameter :: size_square=3
  
  integer NX,NY,it,field_number,ISOURCE,JSOURCE,NPOINTS_PML,nrec
  logical USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX
  
  double precision, dimension(NX,NY) :: image_data_2D
  
  integer, dimension(nrec) :: ix_rec,iy_rec
  logical, parameter :: kurama =.true.
  integer :: ix,iy,irec
  
  character(len=100) :: file_name,system_command1,system_command2,system_command3
	
  integer :: R, G, B
  double precision :: R1, G1, B1
  
  double precision :: normalized_value,max_amplitude
  

  USE_PML_XMIN = .false.
  USE_PML_XMAX = .false. 
  USE_PML_YMIN = .false.
  USE_PML_YMAX = .false.
  NPOINTS_PML = 0

  !       open image file and create system command to convert image to more convenient format
  !       use the "convert" command from ImageMagick http://www.imagemagick.org
  if(field_number == 1) then
     write(file_name,"('image',i6.6,'_Kx.pnm')") it
     write(system_command1, "('convert image',i6.6,'_Kx.pnm ',a,'/imageKx',i6.6,'.png')") it,trim(folder),it
     write(system_command2, "('rm image',i6.6,'_Kx.pnm')") it
  else if(field_number == 2) then
     write(file_name,"('image',i6.6,'_Kz.pnm')") it
     write(system_command1,"('convert image',i6.6,'_Kz.pnm ',a,'/imageKz',i6.6,'.png')") it,trim(folder),it
     write(system_command2,"('rm image',i6.6,'_Kz.pnm')") it
  endif
  
  open(unit=27, file=file_name, status='unknown')
  
  write(27,"('P3')")	! write image in PNM P3 format
  
  write(27,*) NX,NY	! write image size
  write(27,*) '255'	! maximum value of each pixel color
	
  !       compute maximum amplitude
  !max_amplitude = maxval(abs(image_data_2D))
  
  !       image starts in upper-left corner in PNM format
  do iy=1,NY,1
     do ix=1,NX
        
        !       define data as vector component normalized to [-1:1] and rounded to nearest integer
        !       keeping in mind that amplitude can be negative
        normalized_value = image_data_2D(ix,iy) / max_amplitude 
        
        !       suppress values that are outside [-1:+1] to avoid small edge effects
        if(normalized_value < -1.d0) normalized_value = -1.d0
        if(normalized_value > 1.d0) normalized_value = 1.d0
        
        !       draw an orange cross to represent the source
        if((ix >= ISOURCE - width_cross .and. ix <= ISOURCE + width_cross .and. &
             iy >= JSOURCE - thickness_cross .and. iy <= JSOURCE + thickness_cross) .or. &
             (ix >= ISOURCE - thickness_cross .and. ix <= ISOURCE + thickness_cross .and. &
             iy >= JSOURCE - width_cross .and. iy <= JSOURCE + width_cross)) then
           R = 255
           G = 157
           B = 0
           
           !       display two-pixel-thick black frame around the image
        else if(ix <= 2 .or. ix >= NX-1 .or. iy <= 2 .or. iy >= NY-1) then
           R = 0
           G = 0
           B = 0
           
           !       display edges of the PML layers
        else if((USE_PML_XMIN .and. ix == NPOINTS_PML) .or. &
             (USE_PML_XMAX .and. ix == NX - NPOINTS_PML) .or. &
             (USE_PML_YMIN .and. iy == NPOINTS_PML) .or. &
             (USE_PML_YMAX .and. iy == NY - NPOINTS_PML)) then
           R = 255
           G = 150
           B = 0
           
           !       suppress all the values that are below the threshold
        else if(abs(image_data_2D(ix,iy)) <= max_amplitude * cutvect) then
           
           !       use a black or white background for points that are below the threshold
           if(WHITE_BACKGROUND) then
              R = 255
              G = 255
              B = 255
           else
              R = 0
              G = 0
              B = 0
           endif
           
           !       represent regular image points using red if value is positive, blue if negative
        elseif(kurama) then
           normalized_value=normalized_value**POWER_DISPLAY
           call plotcolor(normalized_value,R1,G1,B1)
           R = nint(R1)
           G = nint(G1)
           B = nint(B1)
        else
           if(normalized_value >= 0.d0) then
              R = 255
              G = nint(255.d0-255.d0*normalized_value**POWER_DISPLAY)
              B = nint(255.d0-255.d0*normalized_value**POWER_DISPLAY)
           else
              R = nint(255.d0-255.d0*abs(normalized_value)**POWER_DISPLAY)
              G = nint(255.d0-255.d0*abs(normalized_value)**POWER_DISPLAY)
              B = 255
           endif
        endif
        
        !       draw a green square to represent the receivers
        do irec = 1,nrec
           if((ix >= ix_rec(irec) - size_square .and. ix <= ix_rec(irec) + size_square .and. &
                iy >= iy_rec(irec) - size_square .and. iy <= iy_rec(irec) + size_square) .or. &
                (ix >= ix_rec(irec) - size_square .and. ix <= ix_rec(irec) + size_square .and. &
                iy >= iy_rec(irec) - size_square .and. iy <= iy_rec(irec) + size_square)) then
              !       use dark green color
              R = 30
              G = 180
              B = 60
           endif
        enddo
        
        !       write color pixel
        write(27,"(i3,' ',i3,' ',i3)") R,G,B
        
     enddo
  enddo
  
  !       close file
close(27)

!       call the system to convert image to GIF (can be commented out if "call system" is missing in your compiler)
call system(system_command1)
call system(system_command2)
end subroutine create_color_kernel









subroutine plotcolor(v,r1,g1,b1)
  !======================================================================
  ! Interpolate r1 g1 b1
  !               2016. 5. Kurama OKUBO IPGP
  !======================================================================
  implicit none
  real(8) v,v1,nl
  double precision r1, g1, b1
  real(8),allocatable,dimension(:) :: x,r,g,b
  
  integer n,i
  
  !read colormap data
  open (17, file='../colormap/colormap.dat', status='old')
  read (17, *) nl
  
  n = int(nl)
  allocate( x(n) )
  allocate( r(n) )
  allocate( g(n) )
  allocate( b(n) )
  
  do i = 1, int(nl)
     read (17, *) x(i), r(i), g(i), b(i)
  end do
  close(17)
  !rescale normalized value from [-1 1] to [0 1]
  !Depending on colormap.dat
  
  !when colormap domain is [0 1]
  !v1 = (1.0d0+v)/2.0d0
  
  !when colormap domain is [-1 1]
  v1 = v
  
  call colormap(v1,n,x,r,g,b,r1,g1,b1)
  
end subroutine plotcolor

subroutine colormap(v,n,x,r,g,b,r1,g1,b1)
  !======================================================================
  !v: normalized value [-1 1]
  !n: number of the data in colormap
  !x: ampritude in dataset
  !r: dataset of R
  !g: dataset of G
  !b: dataset of B
  !r1: interpolated R
  !g1 interpolated G
  !b1: interpolated B
  !
  !                  2016.5. Kurama OKUBO IPGP
  !----------------------------------------------------------------------
  
  implicit none
  integer i,n
  real(8) v,x(n),r(n),g(n),b(n)
  double precision ispline, r1, g1, b1
  real(8) b_spline(n),c_spline(n),d_spline(n)

  !for R
  call spline (x, r, b_spline, c_spline, d_spline,n)
  r1 = ispline(v, x, r, b_spline, c_spline, d_spline, n)

  !for G
  call spline (x, g, b_spline, c_spline, d_spline,n)
  g1 = ispline(v, x, g, b_spline, c_spline, d_spline, n)

  !for B
  call spline (x, b, b_spline, c_spline, d_spline,n)
  b1 = ispline(v, x, b, b_spline, c_spline, d_spline, n)


end subroutine colormap


subroutine spline (x, y, b, c, d, n)
  !======================================================================
  !  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
  !  for cubic spline interpolation
  !  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
  !  for  x(i) <= x <= x(i+1)
  !  Alex G: January 2010
  !----------------------------------------------------------------------
  !  input..
  !  x = the arrays of data abscissas (in strictly increasing order)
  !  y = the arrays of data ordinates
  !  n = size of the arrays xi() and yi() (n>=2)
  !  output..
  !  b, c, d  = arrays of spline coefficients
  !  comments ...
  !  spline.f90 program is based on fortran version of program spline.f
  !  the accompanying function fspline can be used for interpolation
  !======================================================================
  implicit none
  integer n
  double precision x(n), y(n), b(n), c(n), d(n)
  integer i, j, gap
  double precision h
  
  gap = n-1
  ! check input
  if ( n < 2 ) return
  if ( n < 3 ) then
     b(1) = (y(2)-y(1))/(x(2)-x(1))   ! linear interpolation
     c(1) = 0.
     d(1) = 0.
     b(2) = b(1)
     c(2) = 0.
     d(2) = 0.
     return
  end if
  !
  ! step 1: preparation
  !
  d(1) = x(2) - x(1)
  c(2) = (y(2) - y(1))/d(1)
  do i = 2, gap
     d(i) = x(i+1) - x(i)
     b(i) = 2.0*(d(i-1) + d(i))
     c(i+1) = (y(i+1) - y(i))/d(i)
     c(i) = c(i+1) - c(i)
  end do
  !
  ! step 2: end conditions
  !
  b(1) = -d(1)
  b(n) = -d(n-1)
  c(1) = 0.0
  c(n) = 0.0
  if(n /= 3) then
     c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
     c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
     c(1) = c(1)*d(1)**2/(x(4)-x(1))
     c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
  end if
  !
  ! step 3: forward elimination
  !
  do i = 2, n
     h = d(i-1)/b(i-1)
     b(i) = b(i) - h*d(i-1)
     c(i) = c(i) - h*c(i-1)
  end do
  !
  ! step 4: back substitution
  !
  c(n) = c(n)/b(n)
  do j = 1, gap
     i = n-j
     c(i) = (c(i) - d(i)*c(i+1))/b(i)
  end do
  !
  ! step 5: compute spline coefficients
  !
  b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
  do i = 1, gap
     b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
     d(i) = (c(i+1) - c(i))/d(i)
     c(i) = 3.*c(i)
  end do
  c(n) = 3.0*c(n)
  d(n) = d(n-1)
end subroutine spline

function ispline(u, x, y, b, c, d, n)
  !======================================================================
  ! function ispline evaluates the cubic spline interpolation at point z
  ! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
  ! where  x(i) <= u <= x(i+1)
  !----------------------------------------------------------------------
  ! input..
  ! u       = the abscissa at which the spline is to be evaluated
  ! x, y    = the arrays of given data points
  ! b, c, d = arrays of spline coefficients computed by spline
  ! n       = the number of data points
  ! output:
  ! ispline = interpolated value at point u
  !=======================================================================
  implicit none
  double precision ispline
  integer n
  double precision  u, x(n), y(n), b(n), c(n), d(n)
  integer i, j, k
  double precision dx
  
  ! if u is ouside the x() interval take a boundary value (left or right)
  if(u <= x(1)) then
     ispline = y(1)
     return
  end if
  if(u >= x(n)) then
     ispline = y(n)
     return
  end if
  
  !*
  !  binary search for for i, such that x(i) <= u <= x(i+1)
  !*
  i = 1
  j = n+1
  do while (j > i+1)
     k = (i+j)/2
     if(u < x(k)) then
        j=k
     else
        i=k
     end if
  end do
  !*
  !  evaluate spline interpolation
  !*
  dx = u - x(i)
  ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
end function ispline
