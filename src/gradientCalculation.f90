subroutine gradientCalculation

  

  !
  !    gradient calculation for Vp and Vs, modified from frechetKernel.f90
  !
  !      Here in this subroutine, kernelP and kernelS mean gradient direction
  !
  !
  !                                            2017.1. N. Fuji (IPGP)
  !                   
  !


  use parameters
  use paramFWI
  implicit none
  character(100) :: tmpfolder
  
  !call paramFrechetReader
  !call vectorAllocateFrechet
  !call ReceiverSourcePositions
  


 
  
  kernelP = 0.d0
  kernelS = 0.d0
  
  
  do iSource = 1,nSource
     i1Source = iSource
     i2Source = iSource
     

     
     isx1 = iisx(i1Source)
     isz1 = iisz(i1Source)
     isx2 = iisx(i2Source)
     isz2 = iisz(i2Source)
     

     do it = 0, nt
        time(it)=dt*dble(it)
     enddo
 

     recl_size=kind(1.e0)*(nx+1)*(nz+1)
     !print *, recl_size
     
     do it = 0,nt,IT_DISPLAY
        
        
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
              write(outfile,'("strainD",I5,".",I5,".",I5,".OPT_del") ') it2,isx2,isz2
           else
              write(outfile,'("strainD",I5,".",I5,".",I5,".CON_del") ') it2,isx2,isz2
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
                 kernelP(ix,iz)=kernelP(ix,iz) &
                      +2.d0*rho(ix,iz)*vp(ix,iz)*IT_DISPLAY*dble(dt)*StrainForward(ix,iz)*StrainBack(ix,iz)
                      !  +IT_DISPLAY*dble(dt)*StrainForward(ix,iz)*StrainBack(ix,iz)
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
              write(outfile,'("strainS",I5,".",I5,".",I5,".OPT_del") ') it2,isx2,isz2
           else
              write(outfile,'("strainS",I5,".",I5,".",I5,".CON_del") ') it2,isx2,isz2
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
                 kernelS(ix,iz)=kernelS(ix,iz) &
                      +2.d0*rho(ix,iz)*vs(ix,iz)*IT_DISPLAY*dble(dt)*StrainForward(ix,iz)*StrainBack(ix,iz)
              enddo
           enddo
           


           
        enddo
        
      
     enddo
  
     

     
  enddo
  
  !tmpfolder="gradientlPsnapshots"
  !call create_color_kernel(kernelP,nx+1,nz+1,iterationIndex,isx1,isz1,iisx(i2Source:i2Source), &
  !     iisz(i2Source:i2Source),1,2,5.d-9,tmpfolder)
  !tmpfolder="gradientSsnapshots"
  !call create_color_kernel(kernelS,nx+1,nz+1,iterationIndex,isx1,isz1,iisx(i2Source:i2Source), &
  !     iisz(i2Source:i2Source),1,2,5.d-9,tmpfolder)
  
end subroutine gradientCalculation
