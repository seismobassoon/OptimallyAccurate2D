program multipleSourcesFWI2D

  ! Computation of the synthetic seismograms in the time domain
  ! using the optimally accurate operators.
  ! 2D PSV heterogeneous medium
  ! CPML or Cerjan boundary conditions
  !
  !					originally from 1997.6  N. Takeuchi
  !                                                     2016.5. N. Fuji
  !                         discon operators (matlab) : 2016.5. O. Ovcharenko
  !                                         colorbars : 2016.5. K. Okubo
  !  
  !                                          cleaning : 2016.6. N. Fuji   
  !                                waveform inversion : 2017.1. N. Fuji
  !                                          Hessian  : 2017.5. N. Fuji

  use parameters
  use paramFWI
  implicit none


  ! Cerjan boundary
  lmargin(1)=NPOINTS_PML
  rmargin(1)=NPOINTS_PML
  lmargin(2)=NPOINTS_PML
  rmargin(2)=NPOINTS_PML
  
  call paramFWIReader

  call vectorAllocate
  
  call vectorAllocateFWI

  call disconConfig ! discontinuity configuration
  
  call ReceiverSourcePositions 

  !!!! for each source we calculate synthetics
  
  ! reading intermediate parameters (vp,vs,rho)
  
  call calstruct( maxnx,maxnz,rhofile,nx,nz,rho )
  call calstruct( maxnx,maxnz,vpfile, nx,nz,vp )
  call calstruct( maxnx,maxnz,vsfile, nx,nz,vs )
    
  call freeConfig

  ! calculate lamda and mu
  call calstruct2(maxnx,maxnz,nx,nz,rho,vp,vs,lam,mu,liquidmarkers)

  ! structuring absorbing boundary

  call calstructBC(maxnx, maxnz,nx,nz,rho,lam,mu,markers,liquidmarkers,zerodisplacement,lmargin,rmargin)

  
  do ir= 1, nReceiver
     nrx(ir)=nrx(ir)+lmargin(1)
     nrz(ir)=nrz(ir)+lmargin(2)
  enddo
  



  ! first forward modelling
    
  iterationIndex=0

  ! forced to write strains

  writingStrain = .true.


  call forwardmodelling


  

 

  ! NF assumes that sources and receivers are placed at the same points

  
  boxnx=nx-rmargin(1)-lmargin(1)
  boxnz=nz-rmargin(2)-lmargin(2)
  !call FourierAllocate

  do while (iterationIndex<numberIteration) 

     iterationIndex=iterationIndex+1
     
     print *, "iterationIndex = ", iterationIndex
     
     call backpropagation

     ! FFT and deconvolution with Ricker wavelet
     ! It allocates also Frechet derivatives
  
     ! NF should comment out

     !call FourierAll

     ! kernelP/S are A^T \delta d
     


     kernelP=0.d0
     kernelS=0.d0
     
     nx=boxnx
     nz=boxnz

     !call approximatedHessian
     
     ! Here we have already ata and atd (i.e. we can do anything we want!)
     ! However, note that ata here is AtA conjugate transpose!!
   

     
     ! NF should use CG inversion scheme from old libraries
     

     !call invbyCG
     
     
     call gradientCalculation 

     recl_size=kind(1.e0)*(boxnx+1)*(boxnz+1)
     
     
     singleStrainForward(:,:)=kernelP(:,:)

     write(outfile,'("./iteratedModels/",I3.3,".vpgrad")'),iterationIndex
     open(1,file=outfile,form='unformatted',access='direct',recl=recl_size)
     write(1,rec=1) singleStrainForward
     close(1)

     singleStrainForward(:,:)=kernelS(:,:)
 
     write(outfile,'("./iteratedModels/",I3.3,".vsgrad")'),iterationIndex
     open(1,file=outfile,form='unformatted',access='direct',recl=recl_size)
     write(1,rec=1) singleStrainForward
     close(1)

     


     if(1.eq.1) then ! we just change kernelP from iz=10 for the singularity problem

        do iz=1,9
           vp(1:nx+1,iz) = vp(1:nx+1,iz) + steplengthVp * kernelP(1:nx+1,iz) * dexp(dble(-11+iz)*0.53)
           vs(1:nx+1,iz) = vp(1:nx+1,iz) + steplengthVs * kernelS(1:nx+1,iz) * dexp(dble(-11+iz)*0.53)
        enddo

        vp(1:nx+1,10:nz+1) = vp(1:nx+1,10:nz+1) + steplengthVp * kernelP(1:nx+1,10:nz+1)
        vs(1:nx+1,10:nz+1) = vs(1:nx+1,10:nz+1) + steplengthVs * kernelS(1:nx+1,10:nz+1)
        call calstruct2(maxnx,maxnz,nx,nz,rho,vp,vs,lam,mu,liquidmarkers)
        call calstructBC(maxnx,maxnz,nx,nz,rho,lam,mu,markers,liquidmarkers,zerodisplacement,lmargin,rmargin)
        print *, "small perturbation"
        

        synx2(:,:)=synx(:,:)
        synz2(:,:)=synz(:,:)


        synx=0.e0
        synz=0.e0

        call forwardmodelling
        
        numeratorG = 0.d0
        denominatorG = 0.d0
        
        
        
        
        synx(:,:) = synx(:,:)-synx2(:,:)
        synz(:,:) = synz(:,:)-synz2(:,:)
        

        delx(0:maxnt,:)=delx(maxnt:0:-1,:)
        delz(0:maxnt,:)=delz(maxnt:0:-1,:)


        ! here, syn is no more syn !!!

        numeratorG = sum(synx(:,:)*delx(:,:))+sum(synz(:,:)*delz(:,:))
        denominatorG = sum(synx(:,:)*synx(:,:))+sum(synz(:,:)*synz(:,:))

        print *, "num, dem", -numeratorG, denominatorG

        alphaVp = -numeratorG/denominatorG*steplengthVp
        alphaVs = -numeratorG/denominatorG*steplengthVs
     
        print *, "alphaVp/Vs = ",  alphaVp,alphaVs
     

     endif

     
     ! new model construction


     
     vp(1:boxnx+1,1:boxnz+1) = vp(1:boxnx+1,1:boxnz+1) + (alphaVp-steplengthVp)*kernelP(1:boxnx+1,1:boxnz+1)
     vs(1:boxnx+1,1:boxnz+1) = vs(1:boxnx+1,1:boxnz+1) + (alphaVs-steplengthVs)*kernelS(1:boxnx+1,1:boxnz+1)



     singleStrainForward(:,:)=vp(1:boxnx+1,1:boxnz+1)*1.e3
     
     write(outfile,'("./iteratedModels/",I3.3,".vpmodel")'),iterationIndex
     open(1,file=outfile,form='unformatted',access='direct',recl=recl_size)
     write(1,rec=1) singleStrainForward
     close(1)

 
     singleStrainForward(:,:)=vs(1:boxnx+1,1:boxnz+1)*1.e3

     write(outfile,'("./iteratedModels/",I3.3,".vsmodel")'),iterationIndex
     open(1,file=outfile,form='unformatted',access='direct',recl=recl_size)
     write(1,rec=1) singleStrainForward
     close(1)


     ! rho model construction from Vp model

     vs = 0.d0

     call vp2rho(boxnx+1,boxnz+1,vp,vs)

     
     singleStrainForward(:,:)=vs(1:boxnx+1,1:boxnz+1)*1.e3

     write(outfile,'("./iteratedModels/",I3.3,".rhomodel")'),iterationIndex
     open(1,file=outfile,form='unformatted',access='direct',recl=recl_size)
     write(1,rec=1) singleStrainForward
     close(1)
     
     
     write(vpfile,'("./iteratedModels/",I3.3,".vpmodel")'),iterationIndex
     write(vsfile,'("./iteratedModels/",I3.3,".vsmodel")'),iterationIndex
     write(rhofile,'("./iteratedModels/",I3.3,".rhomodel")'),iterationIndex


     maxnx=nx+lmargin(1)+rmargin(1)
     maxnz=nz+lmargin(2)+rmargin(2)

     print *, "this is line 232 maxnx,maxnz,nx,nz=",maxnx,maxnz,nx,nz


     
     call calstruct( maxnx,maxnz,rhofile,nx,nz,rho )
     call calstruct( maxnx,maxnz,vpfile, nx,nz,vp )
     call calstruct( maxnx,maxnz,vsfile, nx,nz,vs )
     
     call freeConfig


     call calstruct2(maxnx,maxnz,nx,nz,rho,vp,vs,lam,mu,liquidmarkers)
     call calstructBC(maxnx,maxnz,nx,nz,rho,lam,mu,markers,liquidmarkers,zerodisplacement,lmargin,rmargin)
     call forwardmodelling
     
     iterationIndex=iterationIndex+1

  enddo

  call FourierDeallocate

     
     
  

end program multipleSourcesFWI2D




