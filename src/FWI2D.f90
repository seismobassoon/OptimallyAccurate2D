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
  
  
  !write(12,*) rho
  !write(13,*) vp
  !write(14,*) vs

  ! structuring absorbing boundary

  call calstructBC(maxnx, maxnz,nx,nz,rho,lam,mu,markers,liquidmarkers,zerodisplacement,lmargin,rmargin)

  
  do ir= 1, nReceiver
     nrx(ir)=nrx(ir)+lmargin(1)
     nrz(ir)=nrz(ir)+lmargin(2)
  enddo
  



  ! first forward modelling
    
  iterationIndex=0
  
  print *, "first forward"

  call forwardmodelling
  

  do while (iterationIndex<numberIteration) 
     iterationIndex=iterationIndex+1
     call backpropagation
  
     !print *, nx, maxnx, maxnz,nz

     nx=nx-rmargin(1)-lmargin(1)
     nz=nz-rmargin(1)-lmargin(1)
     
     kernelP=0.d0
     kernelS=0.d0

     call gradientCalculation
     

     recl_size=kind(1.e0)*(nx+1)*(nz+1)
     
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

     print *, "finished gradient"


     vp(1:nx+1,1:nz+1) = vp(1:nx+1,1:nz+1) + steplengthVp * kernelP(1:nx+1,1:nz+1)
     vs(1:nx+1,1:nz+1) = vs(1:nx+1,1:nz+1) + steplengthVs * kernelS(1:nx+1,1:nz+1)
     call calstruct2(maxnx,maxnz,nx,nz,rho,vp,vs,lam,mu,liquidmarkers)
     call calstructBC(maxnx, maxnz,nx,nz,rho,lam,mu,markers,liquidmarkers,zerodisplacement,lmargin,rmargin)
     print *, "small perturbation"

     synx=0.e0
     synz=0.e0

     call forwardmodelling

     numeratorG = 0.d0
     denominatorG = 0.d0


     

     synx(:,:) = obsx(:,:)-synx(:,:)
     synz(:,:) = obsz(:,:)-synz(:,:)

     ! here, syn is no more syn !!!

     numeratorG = sum(synx(:,:)*delx(:,:))+sum(synz(:,:)*delz(:,:))
     denominatorG = sum(synx(:,:)*synx(:,:))+sum(synz(:,:)*synz(:,:))

     print *, "num, dem", numeratorG, denominatorG

     alphaVp = numeratorG/denominatorG*steplengthVp
     alphaVs = numeratorG/denominatorG*steplengthVs
     
     print *, "alphaVp/Vs = ",  alphaVp,alphaVs

     ! new model construction


     nx=nx-rmargin(1)-lmargin(1)
     nz=nz-rmargin(1)-lmargin(1)
     
     vp(1:nx+1,1:nz+1) = vp(1:nx+1,1:nz+1) + (alphaVp-steplengthVp) * kernelP(1:nx+1,1:nz+1)
     vs(1:nx+1,1:nz+1) = vs(1:nx+1,1:nz+1) + (alphaVs-steplengthVs) * kernelS(1:nx+1,1:nz+1)


     recl_size=kind(1.e0)*(nx+1)*(nz+1)

     singleStrainForward(:,:)=vp(1:nx+1,1:nz+1)
     
     write(outfile,'("./iteratedModels/",I3.3,".vpmodel")'),iterationIndex
     open(1,file=outfile,form='unformatted',access='direct',recl=recl_size)
     write(1,rec=1) singleStrainForward
     close(1)

 
     singleStrainForward(:,:)=vs(1:nx+1,1:nz+1)

     write(outfile,'("./iteratedModels/",I3.3,".vsmodel")'),iterationIndex
     open(1,file=outfile,form='unformatted',access='direct',recl=recl_size)
     write(1,rec=1) singleStrainForward
     close(1)


     call calstruct2(maxnx,maxnz,nx,nz,rho,vp,vs,lam,mu,liquidmarkers)
     call calstructBC(maxnx, maxnz,nx,nz,rho,lam,mu,markers,liquidmarkers,zerodisplacement,lmargin,rmargin)
     call forwardmodelling
     
     

  enddo
     

     
     
  

end program multipleSourcesFWI2D




