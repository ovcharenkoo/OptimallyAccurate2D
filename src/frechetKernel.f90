program frechetKernel

  

  ! Computation of Frechet derivatives for Vp and Vs
  ! Strain wavefields should be calculated beforehand with OPT2D.f90
  !

  !
  !
  !                                            2016.6. N. Fuji (IPGP)
  !
  !


  use parameters
  use paramFrechet
  implicit none
  character(100) :: tmpfolder
  
  call paramFrechetReader
  call vectorAllocateFrechet
  call ReceiverSourcePositions
  
  isx1 = iisx(i1Source)
  isz1 = iisz(i1Source)
  isx2 = iisx(i2Source)
  isz2 = iisz(i2Source)

  do it = 0, nt
     time(it)=dt*dble(it)
  enddo

  recl_size=kind(1.e0)*(nx+1)*(nz+1)
  
  do it = 0,nt,IT_DISPLAY

     kernelP = 0.d0
     kernelS = 0.d0
     do it1 = 0,nt,IT_DISPLAY
     
        it2 = -it1+it
        
        

        if(it2.gt.nt) cycle
        if(it2.lt.0) cycle

        ! for Vp sensitivity kernel

        if(optimise) then
           write(outfile,'("strainD",I5,".",I5,".",I5,".OPT_dat") ') it1,isx1,isz1
        else
           write(outfile,'("strainD",I5,".",I5,".",I5,".CON_dat") ') it1,isx1,isz1
        endif
        do j=1,24
           if(outfile(j:j).eq.' ') outfile(j:j)='0'
        enddo
        
        outfile = './strains/'//trim(modelname)//'/'//outfile
        open(1,file=outfile,form='unformatted',access='direct',recl=recl_size)
        read(1,rec=1) singleStrainForward(1:nx+1,1:nz+1)
        close(1)
        StrainForward(:,:) = singleStrainForward(:,:)
        
        
        if(optimise) then
           write(outfile,'("strainD",I5,".",I5,".",I5,".OPT_dat") ') it2,isx2,isz2
        else
           write(outfile,'("strainD",I5,".",I5,".",I5,".CON_dat") ') it2,isx2,isz2
        endif
        do j=1,24
           if(outfile(j:j).eq.' ') outfile(j:j)='0'
        enddo

       
        
        outfile = './strains/'//trim(modelname)//'/'//outfile
        open(1,file=outfile,form='unformatted',access='direct',recl=recl_size)
        read(1,rec=1) singleStrainBack(1:nx+1,1:nz+1)
        close(1)
        StrainBack(:,:) = singleStrainBack(:,:)
        

        !open(1,file='tmp.dat',form='unformatted',access='direct',recl=recl_size)
        !write(1,rec=1) singleStrainBack
        !close(1)
        !stop

        !kernelP= kernelP+IT_DISPLAY*dble(dt)*(StrainForward*StrainBack)
        
        do iz = 1,nz+1
           do ix = 1,nx+1
              kernelP(ix,iz)=kernelP(ix,iz)+IT_DISPLAY*dble(dt)*StrainForward(ix,iz)*StrainBack(ix,iz)
           enddo
        enddo


        ! for Vs sensitivity kernel

        if(optimise) then
           write(outfile,'("strainS",I5,".",I5,".",I5,".OPT_dat") ') it1,isx1,isz1
        else
           write(outfile,'("strainS",I5,".",I5,".",I5,".CON_dat") ') it1,isx1,isz1
        endif
        do j=1,24
           if(outfile(j:j).eq.' ') outfile(j:j)='0'
        enddo
        
        outfile = './strains/'//trim(modelname)//'/'//outfile
        open(1,file=outfile,form='unformatted',access='direct',recl=recl_size)
        read(1,rec=1) singleStrainForward(1:nx+1,1:nz+1)
        close(1)
        StrainForward(:,:) = singleStrainForward(:,:)
        
        
        if(optimise) then
           write(outfile,'("strainS",I5,".",I5,".",I5,".OPT_dat") ') it2,isx2,isz2
        else
           write(outfile,'("strainS",I5,".",I5,".",I5,".CON_dat") ') it2,isx2,isz2
        endif
        do j=1,24
           if(outfile(j:j).eq.' ') outfile(j:j)='0'
        enddo

       
        
        outfile = './strains/'//trim(modelname)//'/'//outfile
        open(1,file=outfile,form='unformatted',access='direct',recl=recl_size)
        read(1,rec=1) singleStrainBack(1:nx+1,1:nz+1)
        close(1)
        StrainBack(:,:) = singleStrainBack(:,:)
        

        !open(1,file='tmp.dat',form='unformatted',access='direct',recl=recl_size)
        !write(1,rec=1) singleStrainBack
        !close(1)
        !stop

        !kernelP= kernelP+IT_DISPLAY*dble(dt)*(StrainForward*StrainBack)
        
        do iz = 1,nz+1
           do ix = 1,nx+1
              kernelS(ix,iz)=kernelS(ix,iz)+IT_DISPLAY*dble(dt)*StrainForward(ix,iz)*StrainBack(ix,iz)
           enddo
        enddo





        ! NF for debugging
        !if(videoornot) then
        !   call create_color_kernel(StrainBack,nx+1,nz+1,it2,isx1,isz1,iisx(2:2),iisz(2:2),1,2,1.d-3)
           
        
        !endif
     enddo   
     
     
     
     if(videoornot) then
        tmpfolder="kernelPsnapshots"
        call create_color_kernel(kernelP,nx+1,nz+1,it,isx1,isz1,iisx(i2Source:i2Source),iisz(i2Source:i2Source),1,2,5.d-9,tmpfolder)
        tmpfolder="kernelSsnapshots"
        call create_color_kernel(kernelS,nx+1,nz+1,it,isx1,isz1,iisx(i2Source:i2Source),iisz(i2Source:i2Source),1,2,5.d-9,tmpfolder)
     endif
     
  enddo
  
  if(videoornot) then
     
     if(optimise) then
        write(outfile,'("frechetP",".",I3,".",I3,".",I3,".",I3,".OPT.mp4") ') isx1,isz1,isx2,isz2
     else
        write(outfile,'("frechetP",".",I3,".",I3,".",I3,".",I3,".CON.mp4") ') isx1,isz1,isx2,isz2
     endif
     do j=1,24
        if(outfile(j:j).eq.' ') outfile(j:j)='0'
     enddo
     
     outfile = './videos/'//trim(modelname)//'/'//trim(outfile)
     
     


      
     commandline="ffmpeg -framerate 5 -pattern_type glob -i 'kernelPsnapshots/*.png' "
     commandline=trim(commandline)//" -c:v libx264 -pix_fmt yuv420p -vf 'scale=trunc(iw/"
     commandline=trim(commandline)//"2)*2:trunc(ih/2"//")*2' "
     commandline=trim(commandline)//" "//trim(outfile)
     print *, commandline
     call system(commandline)
     



     
     if(optimise) then
        write(outfile,'("frechetS",".",I3,".",I3,".",I3,".",I3,".OPT.mp4") ') isx1,isz1,isx2,isz2
     else
        write(outfile,'("frechetS",".",I3,".",I3,".",I3,".",I3,".CON.mp4") ') isx1,isz1,isx2,isz2
     endif
     do j=1,24
        if(outfile(j:j).eq.' ') outfile(j:j)='0'
     enddo
     
     outfile = './videos/'//trim(modelname)//'/'//trim(outfile)
     
     


      
     commandline="ffmpeg -framerate 5 -pattern_type glob -i 'kernelSsnapshots/*.png' "
     commandline=trim(commandline)//" -c:v libx264 -pix_fmt yuv420p -vf 'scale=trunc(iw/"
     commandline=trim(commandline)//"2)*2:trunc(ih/2"//")*2' "
     commandline=trim(commandline)//" "//trim(outfile)
     print *, commandline
     call system(commandline)
     


     
  endif
  
     



end program frechetKernel

