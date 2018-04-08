subroutine backpropagation
  use parameters
  use paramFWI
  implicit none
  double precision :: normaliseX(1:nReceiver), normaliseZ(1:nReceiver)
  
  ! record size
  recl_size=(nx+1-rmargin(1)-lmargin(1))*(nz+1-rmargin(2)-lmargin(2))*kind(0e0)
  recl_size_syn=(maxnt+1)*(nReceiver+1)*kind(0e0)

  ! Smoothed version of CONV/OPT operators

  call cales( nx,nz,rho,lam,mu,dt,dx,dz, &
       e1, e2, e3, e4, e5, e6, e7, e8, &
       e13,e14,e15,e16,e17,e18,e19,e20, &
       f1, f2, f3, f4, f5, f6, f7, f8, &
       f13,f14,f15,f16,f17,f18,f19,f20 )

  ! discontinuities
  
  ee12 = 0.d0
  ee34 = 0.d0
  ee56 = 0.d0
  ee65 = 0.d0
  ee78 = 0.d0
  ee87 = 0.d0

  ff12 = 0.d0
  ff34 = 0.d0
  ff56 = 0.d0
  ff65 = 0.d0
  ff78 = 0.d0
  ff87 = 0.d0
  
  if(nDiscon.ne.0) then
 
     ! changing dscr by putting lmargin(1) and (2)
     tmpvaluex=dble(lmargin(1))*dx
     tmpvaluez=dble(lmargin(2))*dz
     do ix=1,nDiscon
        do iz=1,lengthDiscon
           dscr(1,iz,ix)=dscr(1,iz,ix)+tmpvaluex
           dscr(2,iz,ix)=dscr(2,iz,ix)+tmpvaluez
        enddo
     enddo
          
     call cales_discon( nx,nz,rho,lam,mu,dt,dx,dz,e1, e2, e3, e4, e5, e6, e7, e8,&
     e13,e14,e15,e16,e17,e18,e19,e20, &
     f1, f2, f3, f4, f5, f6, f7, f8, &
     f13,f14,f15,f16,f17,f18,f19,f20, & 
     ! hereafter are new variables for cales_discon
     ee12,ee34,ee56,ee65,ee78,ee87, &
     ff12,ff34,ff56,ff65,ff78,ff87, &
     markers,nDiscon,lengthDiscon,dscr)

  endif
  
  if(lengthFreeSurface.ne.0) then
          
     call cales_free( maxnx,nx,nz,rho,lam,mu,dt,dx,dz,e1, e2, e3, e4, e5, e6, e7, e8,&
     e13,e14,e15,e16,e17,e18,e19,e20, &
     f1, f2, f3, f4, f5, f6, f7, f8, &
     f13,f14,f15,f16,f17,f18,f19,f20, & 
     ! hereafter are new variables for cales_discon
     ee12,ee34,ee56,ee65,ee78,ee87, &
     ff12,ff34,ff56,ff65,ff78,ff87, &
     zerodisplacement,lengthFreeSurface,free)

     
  endif



  if(lengthFreeSurface.ne.0) then
     ! changing free by putting lmargin(1) and (2)
     tmpvaluex=dble(lmargin(1))*dx
     tmpvaluez=dble(lmargin(2))*dz
     
     do ix=1,lengthFreeSurface
           free(1,ix)=free(1,ix)+tmpvaluex
           free(2,ix)=free(2,ix)+tmpvaluez           
     enddo
  endif

  ! for Cerjan absorbing boundary


  weightBC=1.d0
     
  call compNRBCpre(weightBC(1:nx+1,1:nz+1),CerjanRate,lmargin,rmargin,nx+1,nz+1)
  


  do iSource = 1, nSource

    
     isx=iisx(iSource)+lmargin(1)
     isz=iisz(iSource)+lmargin(2)
     
     ist=nt/4
    



     ! Here in this routine, I do not calculate for the sources themselves I back propagate delta d
     ! 

     ! Reading OBS data


     if(trim(extentionOBSx).ne."9999") then 
        write(outfile,'(I5,".",I5,".") ')  &
             isx-lmargin(1),isz-lmargin(2)
        
        do j=1,12
           if(outfile(j:j).eq.' ') outfile(j:j)='0'
        enddo
        
      
 
        outfile = trim(outfile)//trim(extentionOBSx)       
        outfile = trim(obsdir)//'/'//trim(outfile)
       
       
        open(1,file=outfile,form='unformatted',access='direct',recl=recl_size_syn)
        read(1,rec=1) obsx(0:maxnt,1:nReceiver)
        close(1)
   
     endif

     if(trim(extentionOBSz).ne."9999") then
     
        
        write(outfile,'(I5,".",I5,".") ') &
             isx-lmargin(1),isz-lmargin(2)
        
        
        do j=1,12
           if(outfile(j:j).eq.' ') outfile(j:j)='0'
        enddo
        
        outfile = trim(outfile)//trim(extentionOBSz)
        outfile = trim(obsdir)//'/'//trim(outfile)
        
        open(1,file=outfile,form='unformatted',access='direct',recl=recl_size_syn)
        read(1,rec=1) obsz(0:maxnt,1:nReceiver)
        close(1)
        
     endif

     ! Reading OBS data done


     ! Reading SYN data

  
     if(optimise) then
        write(outfile,'(I5,".",I5,".OPT_UX") ')  &
             isx-lmargin(1),isz-lmargin(2)
     else
        write(outfile,'(I5,".",I5,".CON_UX") ') &
             isx-lmargin(1),isz-lmargin(2)
     endif
     
     do j=1,12
        if(outfile(j:j).eq.' ') outfile(j:j)='0'
     enddo
     
     outfile = './synthetics/'//trim(modelname)//'/'//outfile
     
     open(1,file=outfile,form='unformatted',access='direct',recl=recl_size_syn)
     read(1,rec=1) synx(0:maxnt,1:nReceiver)
     close(1)

     
     
     if(optimise) then
        write(outfile,'(I5,".",I5,".OPT_UZ") ') &
             isx-lmargin(1),isz-lmargin(2)
     else
        write(outfile,'(I5,".",I5,".CON_UZ") ') &
             isx-lmargin(1),isz-lmargin(2)
     endif
     
     do j=1,12
        if(outfile(j:j).eq.' ') outfile(j:j)='0'
     enddo
     
     outfile = './synthetics/'//trim(modelname)//'/'//outfile
   
     
     open(1,file=outfile,form='unformatted',access='direct',recl=recl_size_syn)
     read(1,rec=1) synz(0:maxnt,1:nReceiver)
     close(1)
   
     ! Reading SYN data done
     

     ! taking waveform difference
     delx(:,:)=obsx(:,:)-synx(:,:)
     delz(:,:)=obsz(:,:)-synz(:,:)

     ! reverse delx,delz

     delx(0:maxnt,:)=delx(maxnt:0:-1,:)
     delz(0:maxnt,:)=delz(maxnt:0:-1,:)

     do ir=1,nReceiver
        normaliseX(ir)=maxval(abs(delx(0:maxnt/2,ir)))
        normaliseZ(ir)=maxval(abs(delz(0:maxnt/2,ir)))
        
        !print *, normaliseX,normaliseZ
     enddo

     !stop
    ! do it=1,nt
    ! write(13,*) synz(it,5)
    ! enddo
    ! do it=1,nt
    !   write(14,*) obsz(it,5)
    !enddo

    !do it=1,nt
    !   write(15,*) delz(it,5)
    !enddo
    ! stop


     ! for video (without boundary)

     
     ALPHA_MAX_PML = 2.d0*PI*(f0/2.d0) ! from Festa and Vilotte
     
     ! Initializing the data

     ux=0.d0
     uz=0.d0
     ux1=0.d0
     uz1=0.d0
     ux2=0.d0
     uz2=0.d0
     work=0.d0
     

     
 
     ! R. Courant et K. O. Friedrichs et H. Lewy (1928)
     cp=maxval(vp)
     Courant_number = cp * dt * sqrt(1.d0/dx**2 + 1.d0/dz**2)
     print *, 'backpropagating'
     print *, 'Courant number is', Courant_number
     

     
     t=0.d0
     time(0)=t


     

     do it=0,nt
      
        
           

        fx=0.d0
        fz=0.d0



        ! adjoint sources
        do ir=1,nReceiver
        
           fx(nrx(ir),nrz(ir)) = delx(it,ir)!/normaliseZ(ir)
           fz(nrx(ir),nrz(ir)) = delz(it,ir)!/normaliseZ(ir)         

           !print *, fz(nrx(ir),nrz(ir)), nrx(ir),nrz(ir)
        enddo
       
        !do ir=1,nReceiver
        !   print *, fz(nrx(ir),nrz(ir)),nrx(ir),nrz(ir)
        !enddo


        !else
           call calstep_discon( nx,nz, &
                e1, e2, e3, e4, e5, e6, e7, e8, &
                e13,e14,e15,e16,e17,e18,e19,e20, &
                f1, f2, f3, f4, f5, f6, f7, f8, &
                f13,f14,f15,f16,f17,f18,f19,f20, &
                ux,uz,ux1,ux2,uz1,uz2,isx,isz,fx,fz, &
                work(1,1), work(1,5), work(1,9),work(1,13), &
                work(1,17),work(1,18),work(1,20),work(1,21), &
                work(1,23),work(1,24),work(1,28),work(1,29), optimise, & 
                ! Hereafter are new variables for cales_discon
                ee12,ee34,ee56,ee65,ee78,ee87, &
                ff12,ff34,ff56,ff65,ff78,ff87)
           

        !endif


           if(lengthFreeSurface.ne.0) then
              do iz=1,nz+1
                 do ix=1,nx+1
                    if(zerodisplacement(ix,iz).eq.1) then
                       ux(ix,iz) = 0.d0
                       uz(ix,iz) = 0.d0
                    endif
                 enddo
              enddo
           endif


           ! increment of t
        t = t + dt
        time(it)=t
        !write(14,*) real(t),real(ux(nrx,nrz)),real(uz(nrx,nrz))
       
        
     
        
        ! applying Cerjan boundary
        
        do iz=1,nz+1
           do ix=1,nx+1
              uz(ix,iz)=uz(ix,iz)*weightBC(ix,iz)
              ux1(ix,iz)=ux1(ix,iz)*weightBC(ix,iz)
              uz1(ix,iz)=uz1(ix,iz)*weightBC(ix,iz)
              ux2(ix,iz)=ux2(ix,iz)*weightBC(ix,iz)
              uz2(ix,iz)=uz2(ix,iz)*weightBC(ix,iz)
           enddo
        enddo
        
        
        
        ! calculating strains
        
        
        if(writingStrain.and.(mod(it,IT_DISPLAY).eq.0)) then
           singleStrainDiagonal=0.e0
           tmpsingleStrain=0.e0
           call calStrainDiagonal(nx,nz,ux,uz,lmargin,rmargin,singleStrainDiagonal)
           call calStrainShear(nx,nz,ux,uz,lmargin,rmargin,singleStrainShear)


           if(optimise) then
              write(outfile,'("strainD",I5,".",I5,".",I5,".OPT_del") ') it,isx-lmargin(1),isz-lmargin(2)
           else
              write(outfile,'("strainD",I5,".",I5,".",I5,".CON_del") ') it,isx-lmargin(1),isz-lmargin(2)
           endif
           do j=1,24
              if(outfile(j:j).eq.' ') outfile(j:j)='0'
           enddo
           
           outfile = './strains/'//trim(modelname)//'/'//outfile
           open(1,file=outfile,form='unformatted',access='direct',recl=recl_size)

           tmpsingleStrain(1:nx+1-lmargin(1)-rmargin(1),1:nz+1-lmargin(2)-rmargin(2)) = &
                singleStrainDiagonal(lmargin(1)+1:nx+1-rmargin(1),lmargin(2)+1:nz+1-rmargin(2))
           write(1,rec=1)  tmpsingleStrain
           close(1,status='keep')


           
           if(optimise) then
              write(outfile,'("strainS",I5,".",I5,".",I5,".OPT_del") ') it,isx-lmargin(1),isz-lmargin(2)
           else
              write(outfile,'("strainS",I5,".",I5,".",I5,".CON_del") ') it,isx-lmargin(1),isz-lmargin(2)
           endif
           do j=1,24
              if(outfile(j:j).eq.' ') outfile(j:j)='0'
           enddo
           
           outfile = './strains/'//trim(modelname)//'/'//outfile
           open(1,file=outfile,form='unformatted',access='direct',recl=recl_size)

           tmpsingleStrain(1:nx+1-lmargin(1)-rmargin(1),1:nz+1-lmargin(2)-rmargin(2)) = &
                singleStrainShear(lmargin(1)+1:nx+1-rmargin(1),lmargin(2)+1:nz+1-rmargin(2))
           write(1,rec=1)  tmpsingleStrain
           close(1,status='keep')




           
        endif


        
        !write(*,*) it, ' of ', nt
        if(mod(it,IT_DISPLAY) == 0)then
          
           if(videoornot) then
              call create_color_image(uz(1:nx+1,1:nz+1),nx+1,nz+1,it,isx,isz, &
                   nrx(1:nReceiver),nrz(1:nReceiver),nReceiver, &
                   NPOINTS_PML,USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX,2)
              
              
           endif
           
         
           
        endif

        
      
     enddo

     !write(18,*) singleStrainDiagonal(:,:)


     if(videoornot) then
        
        if(optimise) then
           write(outfile,'("video",".",I5,".",I5,".OPT.mp4") ') isx-lmargin(1),isz-lmargin(2)
        else
           write(outfile,'("video",".",I5,".",I5,".CON.mp4") ') isx-lmargin(1),isz-lmargin(2)
        endif
        do j=1,24
           if(outfile(j:j).eq.' ') outfile(j:j)='0'
        enddo
        
        outfile = './videos/'//trim(modelname)//'/'//trim(outfile)
        
        
        commandline="ffmpeg -framerate 5 -pattern_type glob -i 'snapshots/*.png' -c:v libx264 -pix_fmt yuv420p -vf 'scale=trunc(iw/"//"2)*2:trunc(ih/2"//")*2' "//trim(outfile)
       
        call system(commandline)
     
     endif
     
     

     
  enddo
end subroutine backpropagation
