
subroutine cales_discon( maxnz,nx,nz,rho,lam,mu,dt,dx,dz,e1, e2, e3, e4, e5, e6, e7, e8,&
     e13,e14,e15,e16,e17,e18,e19,e20, &
     f1, f2, f3, f4, f5, f6, f7, f8, &
     f13,f14,f15,f16,f17,f18,f19,f20, & ! hereafter are new variables for cales_discon
     markers,nDiscon,lengthDiscon,dscr)

  ! in the vicinity of boundaries :
  ! we use Operators for discontinuities (Mizutani 2002; Ovcharenko 2015)
  
  integer maxnz,nx,nz
  double precision rho(maxnz+1,*),lam(maxnz+1,*),mu(maxnz+1,*)
  double precision dt,dx,dz
  double precision  e1(maxnz+1,*), e2(maxnz+1,*), e3(maxnz+1,*)
  double precision  e4(maxnz+1,*), e5(maxnz+1,*), e6(maxnz+1,*)
  double precision  e7(maxnz+1,*), e8(maxnz+1,*)
  double precision e13(maxnz+1,*),e14(maxnz+1,*),e15(maxnz+1,*)
  double precision e16(maxnz+1,*),e17(maxnz+1,*),e18(maxnz+1,*)
  double precision e19(maxnz+1,*),e20(maxnz+1,*)
  double precision  f1(maxnz+1,*), f2(maxnz+1,*), f3(maxnz+1,*)
  double precision  f4(maxnz+1,*), f5(maxnz+1,*), f6(maxnz+1,*)
  double precision  f7(maxnz+1,*), f8(maxnz+1,*)
  double precision f13(maxnz+1,*),f14(maxnz+1,*),f15(maxnz+1,*)
  double precision f16(maxnz+1,*),f17(maxnz+1,*),f18(maxnz+1,*)
  double precision f19(maxnz+1,*),f20(maxnz+1,*)
  integer ix,iz
  double precision dt2,dx2,dz2,dxdz
  
  ! interface markers

  integer :: ik(1:9),jk(1:9),ctr
  integer :: markers(maxnz+1,maxnz+1)
  double precision :: pt0x,pt0z,pt1x,pt1z

  integer :: nDiscon,lengthDiscon ! number of discontinuities
  double precision :: dscr(1:2,1:lengthDiscon,1:nDiscon)

  double precision :: dDiagonal,dDiagonal2 ! sqrt(dx^2+dz^2)
  double precision :: eps ! zero tolerance

  double precision :: xi,zi,distan2 ! intersecion coordinates
  integer :: iLengthDiscon,iDiscon,iInterSection(1:2),err
  double precision :: eta(0:1,1:2)
  double precision :: normal(1:2) ! normal vector
  double precision :: coeftmp(1:4,1:9) ! temporal coef for two points

  integer :: nointersections

  ! temporary small matrices
  
  double precision, dimension (3,3) :: tmpM3,pre_dx2,pre_dy2
  double precision, dimension (4,6) :: tmpM46
  double precision, dimension (6,6) :: tmpM6,tmppM6
  double precision, dimension (6,4) :: tmpM64,pre_dxdy
  
  


  ! Verify that all the coordinates are already with lmargins


  dt2 = dt * dt
  dx2 = dx * dx
  dz2 = dz * dz
  dxdz = dx * dz
  
  dDiagonal2= dx2+dz2
  dDiagonal = sqrt(dDiagonal2)
  
  eps=1.d-8

  ik = 0
  jk = 0

  ctr = 0
 
  do ix = 1,-1,-1
     do iz = 1,-1,-1
        ctr = ctr + 1
        ik(ctr) = ix
        jk(ctr) = iz
     enddo
  enddo
  
  ctr = 0 


  ! modified operators for discontinuities
  
  do iz=2,nz
     do ix=2,nx
        if(markers(ix,iz)>0) then

           pt0x=dble(ix-1)*dx
           pt0z=dble(iz-1)*dz

      
           coeftmp(1:6,1:2,1:9) = 0.d0


           nointersections = 1


           ! ctr = 1 right-top ix+1,iz+1
           ctr = 1
           distan2 = dDiagonal2
           
           
           pt1x = pt0x + dx
           pt1z = pt0z + dz


           call findNearestPoint(pt0x,pt0z,pt1x,pt1z,distan2,xi,zi,lengthDiscon,nDiscon,iInterSection,err,dscr)
           if(err.eq.0) then ! if there's no intersection and we take ordinary operators 
              call xiziEta(xi,zi,pt0x,pt0z,dx,dz,eta)
              call NormalFinder(normal,lengthDiscon,nDiscon,iInterSection,dscr)
              nointersections = nointersections * 0
              else
                 normal=0.d0
                 eta(1,1) = dble(abs(ik(ctr)))
                 eta(1,2) = dble(abs(jk(ctr)))
                 eta(0,1) = 1.d0-eta(1,1)
                 eta(0,2) = 1.d0-eta(1,2)                 
           endif
           
           call MizutaniIso(coeftmp(1:6,1:2,ctr),rho(ix,iz),rho(ix+ik(ctr),iz+jk(ctr)), &
                lam(ix,iz),lam(ix+ik(ctr),iz+jk(ctr)), &
                mu(ix,iz),mu(ix+ik(ctr),iz+jk(ctr)),ik(ctr),jk(ctr),dx,dz,eta,normal)       
           
           ! ctr = 2 right-centre ix+1,iz
           ctr = 2
           distan2 = dx2

           pt1x = pt0x + dx
           pt1z = pt0z

           call findNearestPoint(pt0x,pt0z,pt1x,pt1z,distan2,xi,zi,lengthDiscon,nDiscon,iInterSection,err,dscr)
           if(err.eq.0) then ! if there's no intersection and we take ordinary operators 
              call xiziEta(xi,zi,pt0x,pt0z,dx,dz,eta)
              call NormalFinder(normal,lengthDiscon,nDiscon,iInterSection,dscr)
              nointersections = nointersections * 0
              else
                 normal=0.d0
                 eta(1,1) = dble(abs(ik(ctr)))
                 eta(1,2) = dble(abs(jk(ctr)))
                 eta(0,1) = 1.d0-eta(1,1)
                 eta(0,2) = 1.d0-eta(1,2)                 
           endif
           
           call MizutaniIso(coeftmp(1:6,1:2,ctr),rho(ix,iz),rho(ix+ik(ctr),iz+jk(ctr)), &
                lam(ix,iz),lam(ix+ik(ctr),iz+jk(ctr)), &
                mu(ix,iz),mu(ix+ik(ctr),iz+jk(ctr)),ik(ctr),jk(ctr),dx,dz,eta,normal) 
           


           ! ctr = 3 right-bottom ix+1,iz-1
           ctr = 3
           distan2 = dDiagonal2
           
           pt1x = pt0x + dx
           pt1z = pt0z - dz
           
           call findNearestPoint(pt0x,pt0z,pt1x,pt1z,distan2,xi,zi,lengthDiscon,nDiscon,iInterSection,err,dscr)
           if(err.eq.0) then ! if there's no intersection and we take ordinary operators 
              call xiziEta(xi,zi,pt0x,pt0z,dx,dz,eta)
              call NormalFinder(normal,lengthDiscon,nDiscon,iInterSection,dscr)
              nointersections = nointersections * 0
              else
                 normal=0.d0
                 eta(1,1) = dble(abs(ik(ctr)))
                 eta(1,2) = dble(abs(jk(ctr)))
                 eta(0,1) = 1.d0-eta(1,1)
                 eta(0,2) = 1.d0-eta(1,2)                 
           endif
           
           call MizutaniIso(coeftmp(1:6,1:2,ctr),rho(ix,iz),rho(ix+ik(ctr),iz+jk(ctr)), &
                lam(ix,iz),lam(ix+ik(ctr),iz+jk(ctr)), &
                mu(ix,iz),mu(ix+ik(ctr),iz+jk(ctr)),ik(ctr),jk(ctr),dx,dz,eta,normal) 
           



           



           
           ! ctr = 4 centre-top ix,iz+1
           ctr = 4
           distan2 = dz2

           pt1x = pt0x 
           pt1z = pt0z + dz

             
           call findNearestPoint(pt0x,pt0z,pt1x,pt1z,distan2,xi,zi,lengthDiscon,nDiscon,iInterSection,err,dscr)
           if(err.eq.0) then ! if there's no intersection and we take ordinary operators 
              call xiziEta(xi,zi,pt0x,pt0z,dx,dz,eta)
              call NormalFinder(normal,lengthDiscon,nDiscon,iInterSection,dscr)
              nointersections = nointersections * 0
              else
                 normal=0.d0
                 eta(1,1) = dble(abs(ik(ctr)))
                 eta(1,2) = dble(abs(jk(ctr)))
                 eta(0,1) = 1.d0-eta(1,1)
                 eta(0,2) = 1.d0-eta(1,2)                 
           endif
           
           call MizutaniIso(coeftmp(1:6,1:2,ctr),rho(ix,iz),rho(ix+ik(ctr),iz+jk(ctr)), &
                lam(ix,iz),lam(ix+ik(ctr),iz+jk(ctr)), &
                mu(ix,iz),mu(ix+ik(ctr),iz+jk(ctr)),ik(ctr),jk(ctr),dx,dz,eta,normal) 


           
           
           ! ctr = 5 centre ix,iz
           ctr = 5
           distan2 = 0.d0

           pt1x = pt0x 
           pt1z = pt0z 
           
           
           ! ctr = 6 centre-bottom ix,iz-1
           ctr = 6
           distan2 = dz2
           
           pt1x = pt0x 
           pt1z = pt0z - dz

           call findNearestPoint(pt0x,pt0z,pt1x,pt1z,distan2,xi,zi,lengthDiscon,nDiscon,iInterSection,err,dscr)
           if(err.eq.0) then ! if there's no intersection and we take ordinary operators 
              call xiziEta(xi,zi,pt0x,pt0z,dx,dz,eta)
              call NormalFinder(normal,lengthDiscon,nDiscon,iInterSection,dscr)
              nointersections = nointersections * 0
              else
                 normal=0.d0
                 eta(1,1) = dble(abs(ik(ctr)))
                 eta(1,2) = dble(abs(jk(ctr)))
                 eta(0,1) = 1.d0-eta(1,1)
                 eta(0,2) = 1.d0-eta(1,2)                 
           endif
           
           call MizutaniIso(coeftmp(1:6,1:2,ctr),rho(ix,iz),rho(ix+ik(ctr),iz+jk(ctr)), &
                lam(ix,iz),lam(ix+ik(ctr),iz+jk(ctr)), &
                mu(ix,iz),mu(ix+ik(ctr),iz+jk(ctr)),ik(ctr),jk(ctr),dx,dz,eta,normal) 







           
           
           ! ctr = 7 left-top ix-1,iz+1
           ctr = 7
           distan2 = dDiagonal2
            
           pt1x = pt0x - dx
           pt1z = pt0z + dz

           call findNearestPoint(pt0x,pt0z,pt1x,pt1z,distan2,xi,zi,lengthDiscon,nDiscon,iInterSection,err,dscr)
           if(err.eq.0) then ! if there's no intersection and we take ordinary operators 
              call xiziEta(xi,zi,pt0x,pt0z,dx,dz,eta)
              call NormalFinder(normal,lengthDiscon,nDiscon,iInterSection,dscr)
              nointersections = nointersections * 0
              else
                 normal=0.d0
                 eta(1,1) = dble(abs(ik(ctr)))
                 eta(1,2) = dble(abs(jk(ctr)))
                 eta(0,1) = 1.d0-eta(1,1)
                 eta(0,2) = 1.d0-eta(1,2)                 
           endif
           
           call MizutaniIso(coeftmp(1:6,1:2,ctr),rho(ix,iz),rho(ix+ik(ctr),iz+jk(ctr)), &
                lam(ix,iz),lam(ix+ik(ctr),iz+jk(ctr)), &
                mu(ix,iz),mu(ix+ik(ctr),iz+jk(ctr)),ik(ctr),jk(ctr),dx,dz,eta,normal) 


            
           ! ctr = 8 left-centre ix-1,iz
           ctr = 8
           distan2 = dx2
           
           pt1x = pt0x - dx
           pt1z = pt0z 

           call findNearestPoint(pt0x,pt0z,pt1x,pt1z,distan2,xi,zi,lengthDiscon,nDiscon,iInterSection,err,dscr)
           if(err.eq.0) then ! if there's no intersection and we take ordinary operators 
              call xiziEta(xi,zi,pt0x,pt0z,dx,dz,eta)
              call NormalFinder(normal,lengthDiscon,nDiscon,iInterSection,dscr)
              nointersections = nointersections * 0
              else
                 normal=0.d0
                 eta(1,1) = dble(abs(ik(ctr)))
                 eta(1,2) = dble(abs(jk(ctr)))
                 eta(0,1) = 1.d0-eta(1,1)
                 eta(0,2) = 1.d0-eta(1,2)                 
           endif
           
           call MizutaniIso(coeftmp(1:6,1:2,ctr),rho(ix,iz),rho(ix+ik(ctr),iz+jk(ctr)), &
                lam(ix,iz),lam(ix+ik(ctr),iz+jk(ctr)), &
                mu(ix,iz),mu(ix+ik(ctr),iz+jk(ctr)),ik(ctr),jk(ctr),dx,dz,eta,normal) 

           
           
            
           ! ctr = 9 left-bottom ix-1,iz-1
           ctr = 9
           distan2 = dDiagonal2

           pt1x = pt0x - dx
           pt1z = pt0z - dz

           call findNearestPoint(pt0x,pt0z,pt1x,pt1z,distan2,xi,zi,lengthDiscon,nDiscon,iInterSection,err,dscr)
           if(err.eq.0) then ! if there's no intersection and we take ordinary operators 
              call xiziEta(xi,zi,pt0x,pt0z,dx,dz,eta)
              call NormalFinder(normal,lengthDiscon,nDiscon,iInterSection,dscr)
              nointersections = nointersections * 0
              else
                 normal=0.d0
                 eta(1,1) = dble(abs(ik(ctr)))
                 eta(1,2) = dble(abs(jk(ctr)))
                 eta(0,1) = 1.d0-eta(1,1)
                 eta(0,2) = 1.d0-eta(1,2)                 
           endif
           
           call MizutaniIso(coeftmp(1:6,1:2,ctr),rho(ix,iz),rho(ix+ik(ctr),iz+jk(ctr)), &
                lam(ix,iz),lam(ix+ik(ctr),iz+jk(ctr)), &
                mu(ix,iz),mu(ix+ik(ctr),iz+jk(ctr)),ik(ctr),jk(ctr),dx,dz,eta,normal) 

           


           ! Now we have all the elements for coeftmp(1:6,1:2,1:9) 
           
           if(nointersections.eq.1) cycle ! for smoothed points
           
           ! for derivatives of ux


           tmpM3 = 0.d0
           
           tmpM3(1,1) = coeftmp(1,1,2)
           tmpM3(1,2) = coeftmp(2,1,2)
           tmpM3(1,3) = coeftmp(4,1,2)

           tmpM3(2,1) = coeftmp(1,1,5)
           tmpM3(2,2) = coeftmp(2,1,5)
           tmpM3(2,3) = coeftmp(4,1,5)

           tmpM3(3,1) = coeftmp(1,1,8)
           tmpM3(3,2) = coeftmp(2,1,8)
           tmpM3(3,3) = coeftmp(4,1,8)

           call inverse(3,tmpM3,3,pre_dx2)


           tmpM3 = 0.d0

           tmpM3(1,1) = coeftmp(1,1,4)
           tmpM3(1,2) = coeftmp(3,1,4)
           tmpM3(1,3) = coeftmp(5,1,4)

           tmpM3(2,1) = coeftmp(1,1,5)
           tmpM3(2,2) = coeftmp(3,1,5)
           tmpM3(2,3) = coeftmp(5,1,5)

           tmpM3(3,1) = coeftmp(1,1,6)
           tmpM3(3,2) = coeftmp(3,1,6)
           tmpM3(3,3) = coeftmp(5,1,6)

           call inverse(3,tmpM3,3,pre_dy2)

           tmpM46 = 0.d0
           
           tmpM46(1,1:6) = coeftmp(1:6,1,1)
           tmpM46(2,1:6) = coeftmp(1:6,1,3)
           tmpM46(3,1:6) = coeftmp(1:6,1,7)
           tmpM46(4,1:6) = coeftmp(1:6,1,9)

           tmpM6 = 0.d0
           tmpM64 = transpose(tmpM64)
           tmpM6 = matmul(tmpM64,tmpM46)
          
           tmppM6 =0.d0
           
           call inverse(6,tmpM6,6,tmppM6)
           pre_dxdy=matmul(tmppM6,tmpM64)

           
           

           ! modified operators for interfaces
           !                   (Oleg Ovcharenko)

           e1(ix,iz) = dt2 / rho(ix,iz) &
                * (  lam(ix,iz)  &
                + 2.d0 *  mu(ix,iz)  ) &
                * pre_dx2(3,1) &
                / ( dx2 )

           e2(ix,iz) = dt2 / rho(ix,iz) &
                * (  lam(ix,iz)  &
                + 2.d0 *  mu(ix,iz)  ) &
                * pre_dx2(3,3) &
                / ( dx2 )
           
           ee12(ix,iz) = dt2 / rho(ix,iz) &
                * (  lam(ix,iz)  &
                + 2.d0 *  mu(ix,iz)  ) &
                * (pre_dx2(3,1)+pre_dx2(3,2)+pre_dx2(3,3)) &
                / ( dx2 )
           
           
           e3(ix,iz) = dt2 / rho(ix,iz) &
                * mu(ix,iz)  &
                * pre_dy2(3,1) &
                / ( dz2 )


           e4(ix,iz) = dt2 / rho(ix,iz) &
                * mu(ix,iz)  &
                * pre_dy2(3,3) &
                / ( dz2 )
           
           ee34(ix,iz) = dt2 / rho(ix,iz) &
                * mu(ix,iz)  &
                * (pre_dy2(3,1)+pre_dy2(3,2)+pre_dy2(3,3)) &
                / ( dz2 )

           

           f5(ix,iz) = dt2 / rho(ix,iz) * mu(ix,iz) &
                * (-pre_dxdy(6,3)) &
                / ( 4.d0 * dxdz )
           f6(ix,iz) = dt2 / rho(ix,iz) * mu(ix,iz) &
                * pre_dxdy(6,1) &
                / ( 4.d0 * dxdz )
           ff56(ix,iz) =  dt2 / rho(ix,iz) * mu(ix,iz) &
                * (pre_dxdy(6,3)+pre_dxdy(6,4)) &
                / ( 4.d0 * dxdz )
           ff65(ix,iz) =  dt2 / rho(ix,iz) * mu(ix,iz) &
                * (pre_dxdy(6,1)+pre_dxdy(6,2)) &
                / ( 4.d0 * dxdz )
                
                
           f7(ix,iz) = dt2 / rho(ix,iz) * lam(ix,iz) &
                * (-pre_dxdy(6,2)) &
                / ( 4.d0 * dxdz )                     
           f8(ix,iz) = dt2 / rho(ix,iz) * lam(ix,iz) &
                * pre_dxdy(6,1) &
                / ( 4.d0 * dxdz )
           f78(ix,iz) =  dt2 / rho(ix,iz) * lam(ix,iz) &
                * (pre_dxdy(6,2)+pre_dxdy(6,4)) &
                / ( 4.d0 * dxdz )
           f87(ix,iz) =  dt2 / rho(ix,iz) * lam(ix,iz) &
                * (pre_dxdy(6,1)+pre_dxdy(6.3)) &
                / ( 4.d0 * dxdz )
                


           
           ! for derivatives of uy


           tmpM3 = 0.d0
           
           tmpM3(1,1) = coeftmp(1,2,2)
           tmpM3(1,2) = coeftmp(2,2,2)
           tmpM3(1,3) = coeftmp(4,2,2)

           tmpM3(2,1) = coeftmp(1,2,5)
           tmpM3(2,2) = coeftmp(2,2,5)
           tmpM3(2,3) = coeftmp(4,2,5)

           tmpM3(3,1) = coeftmp(1,2,8)
           tmpM3(3,2) = coeftmp(2,2,8)
           tmpM3(3,3) = coeftmp(4,2,8)

           call inverse(3,tmpM3,3,pre_dx2)


           tmpM3 = 0.d0

           tmpM3(1,1) = coeftmp(1,2,4)
           tmpM3(1,2) = coeftmp(3,2,4)
           tmpM3(1,3) = coeftmp(5,2,4)

           tmpM3(2,1) = coeftmp(1,2,5)
           tmpM3(2,2) = coeftmp(3,2,5)
           tmpM3(2,3) = coeftmp(5,2,5)

           tmpM3(3,1) = coeftmp(1,2,6)
           tmpM3(3,2) = coeftmp(3,2,6)
           tmpM3(3,3) = coeftmp(5,2,6)

           call inverse(3,tmpM3,3,pre_dy2)

           tmpM46 = 0.d0
           
           tmpM46(1,1:6) = coeftmp(1:6,2,1)
           tmpM46(2,1:6) = coeftmp(1:6,2,3)
           tmpM46(3,1:6) = coeftmp(1:6,2,7)
           tmpM46(4,1:6) = coeftmp(1:6,2,9)

           tmpM6 = 0.d0
           tmpM64 = transpose(tmpM64)
           tmpM6 = matmul(tmpM64,tmpM46)
          
           tmppM6 =0.d0
           
           call inverse(6,tmpM6,6,tmppM6)
           pre_dxdy=matmul(tmppM6,tmpM64)

           f1(ix,iz) = dt2 / rho(ix,iz) &
                *  mu(ix,iz)  &
                * pre_dx2(3,1) &
                / ( dx2 )

           
           f2(ix,iz) = dt2 / rho(ix,iz) &
                *  mu(ix,iz)  &
                * pre_dx2(3,3) &
                / ( dx2 )

           ff12(ix,iz) = dt2 / rho(ix,iz) &
                *  mu(ix,iz)  &
                * (pre_dx2(3,1)+pre_dx2(3,2)+pre_dx2(3,3)) &
                / ( dx2 )


           f3(ix,iz) = dt2 / rho(ix,iz) &
                * ( lam(ix,iz)  &
                + 2.d0 *  mu(ix,iz)  ) &
                * pre_dy2(3,1) &
                / ( 2.d0 * dz2 )


           f4(ix,iz) = dt2 / rho(ix,iz) &
                * ( lam(ix,iz) &
                + 2.d0 *  mu(ix,iz) ) &
                * pre_dy2(3,3) &
                / ( 2.d0 * dz2 )

           ff34(ix,iz) = dt2 / rho(ix,iz) &
                * ( lam(ix,iz) &
                + 2.d0 *  mu(ix,iz) ) &
                * (pre_dy2(3,1)+pre_dy2(3,2)+pre_dy2(3,3)) &
                / ( 2.d0 * dz2 )
      

           e5(ix,iz) = dt2 / rho(ix,iz) * lam(ix,iz) &
                * (-pre_dxdy(6,3)) &
                / ( 4.d0 * dxdz )

           e6(ix,iz) = dt2 / rho(ix,iz) * lam(ix,iz) &
                * pre_dxdy(6,1) &
                / ( 4.d0 * dxdz )

           ee56(ix,iz) = dt2 / rho(ix,iz) * lam(ix,iz) &
                * (pre_dxdy(6,3)+pre_dxdy(6,4)) &
                / ( 4.d0 * dxdz )

           ee65(ix,iz) = dt2 / rho(ix,iz) * lam(ix,iz) &
                * (pre_dxdy(6,1)+pre_dxdy(6,2)) &
                / ( 4.d0 * dxdz )

           e7(ix,iz) = dt2 / rho(ix,iz) * mu(ix,iz) &
                * (-pre_dxdy(6,2)) &
                / ( 4.d0 * dxdz )

           e8(ix,iz) = dt2 / rho(ix,iz) * mu(ix,iz) &
                * pre_dxdy(6,1) &
                / ( 4.d0 * dxdz )

           ee78(ix,iz) = dt2 / rho(ix,iz) * mu(ix,iz) &
                * (pre_dxdy(6,2)+pre_dxdy(6,4)) &
                / ( 4.d0 * dxdz )
           
           ee87(ix,iz) = dt2 / rho(ix,iz) * mu(ix,iz) &
                * (pre_dxdy(6,1)+pre_dxdy(6,3)) &
                / ( 4.d0 * dxdz )



        endif
     enddo
  enddo











  return
end subroutine cales_discon






subroutine findNearestPoint(x1,z1,x2,z2,distan2,xi,zi,eta,lengthDiscon,nDiscon,iInterSection,err,dscr)
  implicit none
  double precision :: x1,z1,x2,z2,distan2
  double precision :: xi, zi,eta(0:1,1:2)
  
  integer :: nDiscon,lengthDiscon ! number of discontinuities
  double precision :: dscr(1:2,1:lengthDiscon,1:nDiscon)
  double precision :: distance2(1:lengthDiscon,1:nDiscon)
  integer :: iLengthDiscon,iDiscon,iInterSection(1:2),err
  double precision :: x,z
  
  err=0
  xi=0.d0
  zi=0.d0
  distance2=0.d0

  do iDiscon = 1,nDiscon
     do iLengthDiscon = 1,lengthDiscon
        
        x=dscr(1,iLengthDiscon,iDiscon)
        z=dscr(2,iLengthDiscon,iDiscon)
        
        distance2(iLengthDiscon,iDiscon)=(x-x1)*(x-x1)+(z-z1)*(z-z1)+(x-x2)*(x-x2)+(z-z2)*(z-z2)
        
     enddo
  enddo


  iInterSection(2)=minloc(distance2)
  
  if(distance2(iInterSection(1),iInterSection(2))>distan2) then
     err = 1
  else
     xi = dscr(1,iInterSection(1),iInterSection(2))
     zi = dscr(2,iInterSection(1),iInterSection(2))
     
  endif
end subroutine findNearestPoint
     
        
        
  

  
subroutine xiziEta(xi,zi,x,z,dx,dz,eta)
  implicit none
  double precision :: xi,zi,x,z,dx,dz,eta(0:1,1:2)
  
  eta=0.d0
 

  eta(0,1) = abs(xi-x)/dx
  eta(1,1) = 1.d0-eta(0,1)

  eta(0,2) = abs(zi-z)/dz
  eta(1,2) - 1.d0-eta(0,2)
end subroutine xiziEta




subroutine NormalFinder(normal,lengthDiscon,nDiscon,iInterSection,dscr)
  implicit none
  double precision :: normal(1:2)  
  integer :: nDiscon,lengthDiscon ! number of discontinuities
  double precision :: dscr(1:2,1:lengthDiscon,1:nDiscon)
  integer :: iInterSection(1:2)
  double precision :: xi, zi, xj, zj, xk, zk
  double precision :: s12, s23
  double precision :: dxds, dzds

  normal=0.d0

  if(iInterSection(1).eq.1) then
     
     xi = dscr(1,2,iInterSection(2))
     zi = dscr(2,2,iInterSection(2))
     xj = dscr(1,1,iInterSection(2))
     zj = dscr(2,1,iInterSection(2))
     xk = dscr(1,3,iInterSection(2))
     zk = dscr(2,3,iInterSection(2))
     
  elseif(iInterSection(1).eq.lengthDiscon) then
     xi = dscr(1,iInterSection(1)-1,iInterSection(2))
     zi = dscr(2,iInterSection(1)-1,iInterSection(2))
     xj = dscr(1,iInterSection(1)-2,iInterSection(2))
     zj = dscr(2,iInterSection(1)-2,iInterSection(2))
     xk = dscr(1,iInterSection(1),iInterSection(2))
     zk = dscr(2,iInterSection(1),iInterSection(2))
     
  else
     xi = dscr(1,iInterSection(1),iInterSection(2))
     zi = dscr(2,iInterSection(1),iInterSection(2))
     xj = dscr(1,iInterSection(1)-1,iInterSection(2))
     zj = dscr(2,iInterSection(1)-1,iInterSection(2))
     xk = dscr(1,iInterSection(1)+1,iInterSection(2))
     zk = dscr(2,iInterSection(1)+1,iInterSection(2))
  endif

  s12 = sqrt((xi-xj)*(xi-xj)+(zi-zj)*(zi-zj))
  s23 = sqrt((xi-xk)*(xi-xk)+(xi-xk)*(xi-xk))
  
  dxds = (s23*s23*(xi-xj)+s12*s12*(xk-xi))/(s12*s23*(s12+s23))
  dzds = (s23*s23*(zi-zj)+s12*s12*(zk-zi))/(s12*s23*(s12+s23))

  normal(1)=-dzds
  normal(2)=dxds

end subroutine NormalFinder
  
  
  
  

  