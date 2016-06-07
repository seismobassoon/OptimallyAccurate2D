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
     
     do it1 = 0,nt,IT_DISPLAY
     
        !it2 = nt-it1+it
        
        it2=it1

        if(it2.gt.nt) cycle
        if(it2.lt.0) cycle

        if(optimise) then
           write(outfile,'("strain",I5,".",I5,".",I5,".OPT_dat") ') it1,isx1,isz1
        else
           write(outfile,'("strain",I5,".",I5,".",I5,".CON_dat") ') it1,isx1,isz1
        endif
        do j=1,24
           if(outfile(j:j).eq.' ') outfile(j:j)='0'
        enddo
        
        outfile = './strains/'//trim(modelname)//'/'//outfile
        open(1,file=outfile,form='unformatted',access='direct',recl=recl_size)
        read(1,rec=1) singleStrainForward(1:nx+1,1:nz+1)
        
        StrainForward(:,:) = singleStrainForward(:,:)
        
        
        if(optimise) then
           write(outfile,'("strain",I5,".",I5,".",I5,".OPT_dat") ') it2,isx2,isz2
        else
           write(outfile,'("strain",I5,".",I5,".",I5,".CON_dat") ') it2,isx2,isz2
        endif
        do j=1,24
           if(outfile(j:j).eq.' ') outfile(j:j)='0'
        enddo
        
        outfile = './strains/'//trim(modelname)//'/'//outfile
        open(1,file=outfile,form='unformatted',access='direct',recl=recl_size)
        read(1,rec=1) singleStrainBack
        
        StrainBack(:,:) = singleStrainBack(:,:)

        !kernelP= kernelP+IT_DISPLAY*dble(dt)*(StrainForward*StrainBack)
        ! NF for debugging
        if(videoornot) then
           call create_color_kernel(StrainBack,nx,nz,it2,isx1,isz1,iisx(2:2),iisz(2:2),1,2,1.d-3)
           
        endif
     enddo   
     stop
     
     
     if(videoornot) then
        call create_color_kernel(kernelP,nx,nz,it,isx1,isz1,iisx(2:2),iisz(2:2),1,2,1.d-4)
        
     endif
     
  enddo
  
  if(videoornot) then
     
     if(optimise) then
        write(outfile,'("frechet",".",I5,".",I5,".",I5,".",I5,".OPT.mp4") ') isx1,isz1,isx2,isz2
     else
        write(outfile,'("frechet",".",I5,".",I5,".",I5,".",I5,".CON.mp4") ') isx1,isz1,isx2,isz2
     endif
     do j=1,24
        if(outfile(j:j).eq.' ') outfile(j:j)='0'
     enddo
     
     outfile = './videos/'//trim(modelname)//'/'//outfile
     
     
     commandline="ffmpeg -framerate 5 -pattern_type glob -i 'kernelsnapshots/*.png' -c:v libx264 -pix_fmt yuv420p "//outfile
     
     call system(commandline)
     
  endif
  
     



end program frechetKernel

