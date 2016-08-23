module snapshot_mod

  ! Copyright 2016 Thomas Belahi
  ! Thomas Belahi, Institut de Physique du Globe de Paris

  
  implicit none
  
  integer :: snafla,snasta,snainc,snanum
  
contains
  
  function snap_sample(n,nshots) result(do_snap)
    integer, intent(in) :: n   ! current time sample
    integer, intent(in) :: nshots
    logical :: do_snap
    
    do_snap = .false.
    
    !if (nshots /= 1 .or. snafla == 0 ) return
    if (n < snasta .or. n > snasta+(snanum-1)*snainc) return
    
    do_snap = (modulo(n-snasta,snainc) == 0)
  end function snap_sample
  
  
  subroutine write_snapshots_header(x_start, x_end, z_start, z_end,sn)
    ! subroutine to create a header file that contains information
    ! on the dimension of the binary direct access files snapshots_s01.txx,... saved by
    ! the write_full_fields subroutine
    
    ! the header file is saved as a namelist
    
    use allglobals
    use fd_mod
    
    ! declarations
    integer, intent(in)                :: x_start, x_end, z_start, z_end, sn
    real(kind=kind(1.d0))              :: time_step
    integer                            :: ntimesteps, starting_snapshot, snap_increment, snap_number
    integer                            :: ios, recl_size
    character*256                      :: nameoutput
    character*3                        :: snum
    real(kind=kind(0.e0)), allocatable :: tmp_single(:,:)
    
    
    
    ! declare the namelist
    namelist /HEADER/ x_start, x_end, z_start, z_end, recl_size, ntimesteps,&
         time_step, starting_snapshot, snap_increment, snap_number
    
    !write sn to snum to convert shot number to a string
    write(snum,'(i3.3)')  sn
    !output path for the file
    nameoutput = outpath(1:outlen)//'/snapshots.header'//snum
    
    ! define the size of a snapshot in a binary file
    allocate(tmp_single(x_start:x_end,z_start:z_end))
    tmp_single=0
    inquire(iolength=recl_size) tmp_single(:,:)
    deallocate(tmp_single)
    
    !convert some useful parameters
    ntimesteps        = k
    time_step         = dt
    starting_snapshot = snasta
    snap_increment    = snainc
    snap_number       = snanum
    
    ! write the namelist /HEADER/ to the file fullfield.header
    open(unit=10, file=nameoutput, status='replace', form='formatted', access='sequential', action='write', iostat=ios)
    if (ios/=0) then ! problè´≤me è´¢ l'ouverture
       stop "erreur è´¢ l'ouverture du fichier snapshot.header"
    else
       write(unit=10, nml=HEADER)
       close(unit=10)
    endif
    
  end subroutine write_snapshots_header
  
  
  subroutine write_snapshot(isforward,n, x1, x2, z1, z2, sn)
    use fd_mod
    use allglobals
    logical, intent(in) :: isforward
    integer, intent(in) :: n
    integer, intent(in) :: x1,x2,z1,z2, sn
    
    character :: nameo*256,namep*243
    character:: nnum*6
    character :: snum*3
    integer :: plen, matSizeSingle, nrec
    real(kind=kind(1.d0)) :: t
    real(kind=kind(0.e0)), allocatable :: tmp_single(:,:)
    
    !write sn to snum to convert shot number to a string
    write(snum, '(i3.3)') sn
    
    t = n * dt
    
    if (isforward) then
       namep = outpath(1:outlen)//'/snapshots_s'//snum
       write(nnum,'(i5.5)') n
    else
       namep = outpath(1:outlen)//'/snapshots_back'
       write(nnum,'(i5.5)') k - n
    endif
    
    plen = len_trim(namep)
    
    !define the rec position for the snapshots
    nrec = (n-snasta)/snainc + 1
    
    !allocate the moche tabelau bien mutable qu'on va faire varier comme des sales
    allocate(tmp_single(x1:x2,z1:z2))
    
    !determine the size of a matrix to be passed to recl for writing the binary file
    ! nicely made platform independent and compiler indepent (i.e. Fortran is shit !!)
    inquire(iolength=matSizeSingle) tmp_single(:,:)
    
    nameo=namep(1:plen)//'.txx'
    open(1, file=nameo, status='unknown', form='unformatted', access='direct',  recl=matSizeSingle)
    !convert matrix to single
    tmp_single = txx(x1:x2,z1:z2)
    write(1,rec=nrec) tmp_single
    close(1)
    
    nameo=namep(1:plen)//'.tzz'
    open(1, file=nameo, status='unknown', form='unformatted', access='direct',  recl=matSizeSingle)
    !convert matrix to single
    tmp_single = tzz(x1:x2,z1:z2)
    write(1,rec=nrec) tmp_single
    close(1)
    
    nameo=namep(1:plen)//'.txz'
    open(1, file=nameo, status='unknown', form='unformatted', access='direct',  recl=matSizeSingle)
    !convert matrix to single
    tmp_single = txz(x1:x2,z1:z2)
    write(1,rec=nrec) tmp_single
    close(1)
    
    nameo=namep(1:plen)//'.u'
    open(1, file=nameo, status='unknown', form='unformatted', access='direct',  recl=matSizeSingle)
    !convert matrix to single
    tmp_single = u(x1:x2,z1:z2)
    write(1,rec=nrec) tmp_single
    close(1)
    
    nameo=namep(1:plen)//'.v'
    open(1, file=nameo, status='unknown', form='unformatted', access='direct',  recl=matSizeSingle)
    !convert matrix to single
    tmp_single = v(x1:x2,z1:z2)
    write(1,rec=nrec) tmp_single
    close(1)
    
    nameo=namep(1:plen)//'.rxxu'
    open(1, file=nameo, status='unknown', form='unformatted', access='direct',  recl=matSizeSingle)
    !convert matrix to single
    tmp_single = rxxu(x1:x2,z1:z2)
    write(1,rec=nrec) tmp_single
    close(1)
    
    nameo=namep(1:plen)//'.rzzu'
    open(1, file=nameo, status='unknown', form='unformatted', access='direct',  recl=matSizeSingle)
    !convert matrix to single
    tmp_single = rzzu(x1:x2,z1:z2)
    write(1,rec=nrec) tmp_single
    close(1)
    
    nameo=namep(1:plen)//'.rxzu'
    open(1, file=nameo, status='unknown', form='unformatted', access='direct',  recl=matSizeSingle)
    !convert matrix to single
    tmp_single = rxzu(x1:x2,z1:z2)
    write(1,rec=nrec) tmp_single
    close(1)
    
    deallocate(tmp_single)
  end subroutine write_snapshot
  
  subroutine read_snapshot(isforward,n, x1, x2, z1, z2, sn)
    !il suffit de reprendre write_snapshot et de remplacer write par read et
    !d'inverse l'è´±tape de "casting to single" et de write en read puis
    !"casting to double"
    use fd_mod
    use allglobals
    logical, intent(in) :: isforward
    integer, intent(in) :: n
    integer, intent(in) :: x1,x2,z1,z2, sn
    
    character :: nameo*256,namep*243
    character:: nnum*6
    character :: snum*3
    integer :: plen, matSizeSingle, nrec
    real(kind=kind(1.d0)) :: t
    real(kind=kind(0.e0)), allocatable :: tmp_single(:,:)
    
    !write sn to snum to convert shot number to a string
    write(snum, '(i3.3)') sn
    
    t = n * dt
    
    if (isforward) then
       namep = outpath(1:outlen)//'/snapshots_s'//snum
       write(nnum,'(i5.5)') n
    else
       namep = outpath(1:outlen)//'/snapshots_s'//snum
    endif
    
    plen = len_trim(namep)
    
    !define the rec position for the snapshots
    nrec = (n-snasta)/snainc + 1
    
    !allocate the moche tabelau bien mutable qu'on va faire varier comme des sales
    allocate(tmp_single(x1:x2,z1:z2))
    
    !determine the size of a matrix to be passed to recl for writing the binary file
    ! nicely made platform independent and compiler indepent (i.e. Fortran is shit !!)
    inquire(iolength=matSizeSingle) tmp_single(:,:)
    
    nameo=namep(1:plen)//'.txx'
    open(1, file=nameo, status='unknown', form='unformatted', access='direct',  recl=matSizeSingle)
    read(1,rec=nrec) tmp_single
    !convert matrix to double
    rectxx(x1:x2,z1:z2) = tmp_single(x1:x2,z1:z2)
    close(1)
    
    nameo=namep(1:plen)//'.tzz'
    open(1, file=nameo, status='unknown', form='unformatted', access='direct',  recl=matSizeSingle)
    read(1,rec=nrec) tmp_single
    !convert matrix to double
    rectzz(x1:x2,z1:z2) = tmp_single(x1:x2,z1:z2)
    close(1)
    
    nameo=namep(1:plen)//'.txz'
    open(1, file=nameo, status='unknown', form='unformatted', access='direct',  recl=matSizeSingle)
    read(1,rec=nrec) tmp_single
    !convert matrix to double
    rectxz = tmp_single(x1:x2,z1:z2)
    close(1)
    
    nameo=namep(1:plen)//'.u'
    open(1, file=nameo, status='unknown', form='unformatted', access='direct',  recl=matSizeSingle)
    read(1,rec=nrec) tmp_single
    !convert matrix to double
    recu(x1:x2,z1:z2) = tmp_single(x1:x2,z1:z2)
    close(1)
    
    nameo=namep(1:plen)//'.v'
    open(1, file=nameo, status='unknown', form='unformatted', access='direct',  recl=matSizeSingle)
    read(1,rec=nrec) tmp_single
    !convert matrix to double
    recv(x1:x2,z1:z2) = tmp_single(x1:x2,z1:z2)
    close(1)
    
    nameo=namep(1:plen)//'.rxxu'
    open(1, file=nameo, status='unknown', form='unformatted', access='direct',  recl=matSizeSingle)
    read(1,rec=nrec) tmp_single
    !convert matrix to double
    recrxx(x1:x2,z1:z2) = tmp_single(x1:x2,z1:z2)
    close(1)
    
    nameo=namep(1:plen)//'.rzzu'
    open(1, file=nameo, status='unknown', form='unformatted', access='direct',  recl=matSizeSingle)
    read(1,rec=nrec) tmp_single
    !convert matrix to double
    recrzz(x1:x2,z1:z2) = tmp_single(x1:x2,z1:z2)
    close(1)
    
    nameo=namep(1:plen)//'.rxzu'
    open(1, file=nameo, status='unknown', form='unformatted', access='direct',  recl=matSizeSingle)
    read(1,rec=nrec) tmp_single
    !convert matrix to double
    recrxz(x1:x2,z1:z2) = tmp_single(x1:x2,z1:z2)
    close(1)
    
    deallocate(tmp_single)
  end subroutine read_snapshot
  
end module snapshot_mod
