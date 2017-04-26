


subroutine cales( nx,nz,rho,lam,mu,dt,dx,dz,e1, e2, e3, e4, e5, e6, e7, e8,&
     e13,e14,e15,e16,e17,e18,e19,e20, &
     f1, f2, f3, f4, f5, f6, f7, f8, &
     f13,f14,f15,f16,f17,f18,f19,f20 )
  implicit none
  integer nx,nz
  double precision rho(nx+1,nz+1),lam(nx+1,nz+1),mu(nx+1,nz+1)
  double precision dt,dx,dz
  double precision  e1(nx+1,nz+1), e2(nx+1,nz+1), e3(nx+1,nz+1)
  double precision  e4(nx+1,nz+1), e5(nx+1,nz+1), e6(nx+1,nz+1)
  double precision  e7(nx+1,nz+1), e8(nx+1,nz+1)
  double precision e13(nx+1,nz+1),e14(nx+1,nz+1),e15(nx+1,nz+1)
  double precision e16(nx+1,nz+1),e17(nx+1,nz+1),e18(nx+1,nz+1)
  double precision e19(nx+1,nz+1),e20(nx+1,nz+1)
  double precision  f1(nx+1,nz+1), f2(nx+1,nz+1), f3(nx+1,nz+1)
  double precision  f4(nx+1,nz+1), f5(nx+1,nz+1), f6(nx+1,nz+1)
  double precision  f7(nx+1,nz+1), f8(nx+1,nz+1)
  double precision f13(nx+1,nz+1),f14(nx+1,nz+1),f15(nx+1,nz+1)
  double precision f16(nx+1,nz+1),f17(nx+1,nz+1),f18(nx+1,nz+1)
  double precision f19(nx+1,nz+1),f20(nx+1,nz+1)
  integer ix,iz
  double precision dt2,dx2,dz2,dxdz
  ! for smoothed part of the model :
  !  we use Zahradnik operators and optimally accurate operators

  print *,nx,nz 
  
  dt2 = dt * dt
  dx2 = dx * dx
  dz2 = dz * dz
  dxdz = dx * dz
  
  do iz=2,nz
     do ix=2,nx
        e1(ix,iz) = dt2 / rho(ix,iz) &
             * ( ( lam(ix-1,iz) + lam(ix,iz) ) &
             + 2.d0 * ( mu(ix-1,iz) + mu(ix,iz) ) ) &
             / ( 2.d0 * dx2 )
        e2(ix,iz) = dt2 / rho(ix,iz) &
             * ( ( lam(ix,iz) + lam(ix+1,iz) ) &
             + 2.d0 * ( mu(ix,iz) + mu(ix+1,iz) ) ) &
             / ( 2.d0 * dx2 )
        e3(ix,iz) = dt2 / rho(ix,iz) &
             * ( mu(ix,iz-1) + mu(ix,iz) ) &
             / ( 2.d0 * dz2 )
        e4(ix,iz) = dt2 / rho(ix,iz) &
             * ( mu(ix,iz) + mu(ix,iz+1) ) &
             / ( 2.d0 * dz2 )
        e5(ix,iz) = dt2 / rho(ix,iz) * lam(ix-1,iz) &
             / ( 4.d0 * dxdz )
        e6(ix,iz) = dt2 / rho(ix,iz) * lam(ix+1,iz) &
             / ( 4.d0 * dxdz )
        e7(ix,iz) = dt2 / rho(ix,iz) * mu(ix,iz-1) &
             / ( 4.d0 * dxdz )
        e8(ix,iz) = dt2 / rho(ix,iz) * mu(ix,iz+1) &
             / ( 4.d0 * dxdz )
        e13(ix,iz) = dt2 / rho(ix,iz) * lam(ix-1,iz) &
             * ( -5.d0 ) / ( 1728.d0 * dxdz )
        e14(ix,iz) = dt2 / rho(ix,iz) * lam(ix,iz) &
             * ( -3.d0 ) / ( 1728.d0 * dxdz )
        e15(ix,iz) = dt2 / rho(ix,iz) * lam(ix+1,iz) &
           * ( +9.d0 ) / ( 1728.d0 * dxdz )
        if ( ix+2.le.nx+1 ) then
           e16(ix,iz) = dt2 / rho(ix,iz) * lam(ix+2,iz) &
                * ( -1.d0 ) / ( 1728.d0 * dxdz )
        else
           e16(ix,iz) = 0.d0
        endif
        e17(ix,iz) = dt2 / rho(ix,iz) * mu(ix,iz-1) &
             * ( -5.d0 ) / ( 1728.d0 * dxdz )
        e18(ix,iz) = dt2 / rho(ix,iz) * mu(ix,iz) &
             * ( -3.d0 ) / ( 1728.d0 * dxdz )
        e19(ix,iz) = dt2 / rho(ix,iz) * mu(ix,iz+1) &
             * ( +9.d0 ) / ( 1728.d0 * dxdz )
        if ( iz+2.le.nz+1) then
           e20(ix,iz) = dt2 / rho(ix,iz) * mu(ix,iz+2) &
                * ( -1.d0 ) / ( 1728.d0 * dxdz )
        else
           e20(ix,iz) = 0.d0
        endif
        f1(ix,iz) = dt2 / rho(ix,iz) &
             * ( mu(ix-1,iz) + mu(ix,iz) ) &
             / ( 2.d0 * dx2 )
        f2(ix,iz) = dt2 / rho(ix,iz) &
             * ( mu(ix,iz) + mu(ix+1,iz) ) &
             / ( 2.d0 * dx2 )
        f3(ix,iz) = dt2 / rho(ix,iz) &
             * ( ( lam(ix,iz-1) + lam(ix,iz) ) &
             + 2.d0 * ( mu(ix,iz-1) + mu(ix,iz) ) ) &
             / ( 2.d0 * dz2 )
        f4(ix,iz) = dt2 / rho(ix,iz) &
             * ( ( lam(ix,iz) + lam(ix,iz+1) ) &
           + 2.d0 * ( mu(ix,iz) + mu(ix,iz+1) ) ) &
           / ( 2.d0 * dz2 )
        f5(ix,iz) = dt2 / rho(ix,iz) * mu(ix-1,iz) &
             / ( 4.d0 * dxdz )
        f6(ix,iz) = dt2 / rho(ix,iz) * mu(ix+1,iz) &
             / ( 4.d0 * dxdz )
        f7(ix,iz) = dt2 / rho(ix,iz) * lam(ix,iz-1) &
             / ( 4.d0 * dxdz )
        f8(ix,iz) = dt2 / rho(ix,iz) * lam(ix,iz+1) &
             / ( 4.d0 * dxdz )
        if ( ix-2.ge.1 ) then
           f13(ix,iz) = dt2 / rho(ix,iz) * mu(ix-2,iz) &
                * (  1.d0 ) / ( 1728.d0 * dxdz )
        else
           f13(ix,iz) = 0.d0
        endif
        f14(ix,iz) = dt2 / rho(ix,iz) * mu(ix-1,iz) &
             * ( -9.d0 ) / ( 1728.d0 * dxdz )
        f15(ix,iz) = dt2 / rho(ix,iz) * mu(ix,iz) &
             * (  3.d0 ) / ( 1728.d0 * dxdz )
        f16(ix,iz) = dt2 / rho(ix,iz) * mu(ix+1,iz) &
             * (  5.d0 ) / ( 1728.d0 * dxdz )
        if ( iz-2.ge.1 ) then
           f17(ix,iz) = dt2 / rho(ix,iz) * lam(ix,iz-2) &
                * (  1.d0 ) / ( 1728.d0 * dxdz )
        else
           f17(ix,iz) = 0.d0
        endif
        f18(ix,iz) = dt2 / rho(ix,iz) * lam(ix,iz-1) &
             * ( -9.d0 ) / ( 1728.d0 * dxdz )
        f19(ix,iz) = dt2 / rho(ix,iz) * lam(ix,iz) &
             * (  3.d0 ) / ( 1728.d0 * dxdz )
        f20(ix,iz) = dt2 / rho(ix,iz) * lam(ix,iz+1) &
             * (  5.d0 ) / ( 1728.d0 * dxdz )

        
     enddo
  enddo


  return
end subroutine cales


  

subroutine calstep( nx,nz, &
     e1, e2, e3, e4, e5, e6, e7, e8, &
     e13,e14,e15,e16,e17,e18,e19,e20, &
     f1, f2, f3, f4, f5, f6, f7, f8, &
     f13,f14,f15,f16,f17,f18,f19,f20, &
     ux,uz,ux1,ux2,uz1,uz2,isx,isz,fx,fz, &
     work1,work2,work3,work4, &
     work5,work6,work7,work8, &
     work9,work10,work11,work12,optimise, & ! Hereafter are new variables for cales_discon
     ee12,ee34,ee56,ee65,ee78,ee87, &
     ff12,ff34,ff56,ff65,ff78,ff87)

  integer nx,nz,isx,isz
  double precision ux(nx+1,nz+1),ux1(nx+1,nz+1),ux2(nx+1,nz+1)
  double precision uz(nx+1,nz+1),uz1(nx+1,nz+1),uz2(nx+1,nz+1)
  double precision  e1(nx+1,nz+1), e2(nx+1,nz+1), e3(nx+1,nz+1)
  double precision  e4(nx+1,nz+1), e5(nx+1,nz+1), e6(nx+1,nz+1)
  double precision  e7(nx+1,nz+1), e8(nx+1,nz+1)
  double precision e13(nx+1,nz+1),e14(nx+1,nz+1),e15(nx+1,nz+1)
  double precision e16(nx+1,nz+1),e17(nx+1,nz+1),e18(nx+1,nz+1)
  double precision e19(nx+1,nz+1),e20(nx+1,nz+1)
  double precision  f1(nx+1,nz+1), f2(nx+1,nz+1), f3(nx+1,nz+1)
  double precision  f4(nx+1,nz+1), f5(nx+1,nz+1), f6(nx+1,nz+1)
  double precision  f7(nx+1,nz+1), f8(nx+1,nz+1)
  double precision f13(nx+1,nz+1),f14(nx+1,nz+1),f15(nx+1,nz+1)
  double precision f16(nx+1,nz+1),f17(nx+1,nz+1),f18(nx+1,nz+1)
  double precision f19(nx+1,nz+1),f20(nx+1,nz+1)
  double precision fx(nx+1,nz+1),fz(nx+1,nz+1)
  double precision work1(nx+1,-2:1),work2(nx+1,-1:2)
  double precision work3(nx+1,-2:1),work4(nx+1,-1:2)
  double precision work5(nx+1),work6(nx+1,0:1)
  double precision work7(nx+1),work8(nx+1,0:1)
  double precision work9(nx+1),work10(nx+1,-2:1)
  double precision work11(nx+1),work12(nx+1,-1:2)
  integer ix,iz,iz1,iz2,ix11,ix12,ix21,ix22
  logical optimise


  
  double precision, dimension(nx+1,nz+1) :: ee12,ee34,ee56,ee65,ee78,ee87
  double precision, dimension(nx+1,nz+1) :: ff12,ff34,ff56,ff65,ff78,ff87

  ! predicting the wavefield
 do iz=2,nz
    do ix=2,nx
       ux(ix,iz) = 2.d0 * ux1(ix,iz) - ux2(ix,iz) &
            + e1(ix,iz) * ( ux1(ix-1,iz) - ux1(ix,iz) ) &
            + e2(ix,iz) * ( ux1(ix+1,iz) - ux1(ix,iz) ) &
            + e3(ix,iz) * ( ux1(ix,iz-1) - ux1(ix,iz) ) &
            + e4(ix,iz) * ( ux1(ix,iz+1) - ux1(ix,iz) ) &
            - e5(ix,iz) * ( uz1(ix-1,iz+1) - uz1(ix-1,iz-1) ) &
            + e6(ix,iz) * ( uz1(ix+1,iz+1) - uz1(ix+1,iz-1) ) &
            - e7(ix,iz) * ( uz1(ix+1,iz-1) - uz1(ix-1,iz-1) ) &
            + e8(ix,iz) * ( uz1(ix+1,iz+1) - uz1(ix-1,iz+1) ) & ! Hereafter discon
            + ee12(ix,iz) * ux1(ix,iz) &
            + ee34(ix,iz) * ux1(ix,iz) &
            + ee56(ix,iz) * uz1(ix-1,iz-1) &
            + ee65(ix,iz) * uz1(ix+1,iz-1) &
            + ee78(ix,iz) * uz1(ix-1,iz-1) &
            + ee87(ix,iz) * uz1(ix+1,iz-1)
       



       uz(ix,iz) = 2.d0 * uz1(ix,iz) - uz2(ix,iz) &
            + f1(ix,iz) * ( uz1(ix-1,iz) - uz1(ix,iz) ) &
            + f2(ix,iz) * ( uz1(ix+1,iz) - uz1(ix,iz) ) &
            + f3(ix,iz) * ( uz1(ix,iz-1) - uz1(ix,iz) ) &
            + f4(ix,iz) * ( uz1(ix,iz+1) - uz1(ix,iz) ) &
            - f5(ix,iz) * ( ux1(ix-1,iz+1) - ux1(ix-1,iz-1) ) &
            + f6(ix,iz) * ( ux1(ix+1,iz+1) - ux1(ix+1,iz-1) ) &
            - f7(ix,iz) * ( ux1(ix+1,iz-1) - ux1(ix-1,iz-1) ) &
            + f8(ix,iz) * ( ux1(ix+1,iz+1) - ux1(ix-1,iz+1) ) & ! Hereafter discon
            + ff12(ix,iz) * uz1(ix,iz) &
            + ff34(ix,iz) * uz1(ix,iz) &
            + ff56(ix,iz) * ux1(ix-1,iz-1) &
            + ff65(ix,iz) * ux1(ix+1,iz-1) &
            + ff78(ix,iz) * ux1(ix-1,iz-1) &
            + ff87(ix,iz) * ux1(ix+1,iz-1)

    enddo
 enddo
 ux(isx,isz) = ux(isx,isz) + fx(isx,isz)
 uz(isx,isz) = uz(isx,isz) + fz(isx,isz)

 
 if(optimise) then
    ! correcting the wavefield
    !
    do ix=2,nx
       iz1 = 2
       iz2 = 3
       work1(ix,-2) = 0.d0
       work1(ix,-1) = 0.d0
       work1(ix,0) = 0.d0
       work1(ix,1) = ux(ix,iz1) - 2.d0 * ux1(ix,iz1) + ux2(ix,iz1)
       work2(ix,-1) = 0.d0
       work2(ix,0) = 0.d0
       work2(ix,1) = uz(ix,iz1) - 2.d0 * uz1(ix,iz1) + uz2(ix,iz1)
       work2(ix,2) = uz(ix,iz2) - 2.d0 * uz1(ix,iz2) + uz2(ix,iz2)
       work3(ix,-2) = 0.d0
       work3(ix,-1) = 0.d0
       work3(ix,0) = 0.d0
       work3(ix,1) = work1(ix,1) + 12.d0 * ux1(ix,iz1)
       work4(ix,-1) = 0.d0
       work4(ix,0) = 0.d0
       work4(ix,1) = work2(ix,1) + 12.d0 * uz1(ix,iz1)
       work4(ix,2) = work2(ix,2) + 12.d0 * uz1(ix,iz2)
    enddo
    
    do ix=1,nx+1
       ix11 = max0( ix-1,1 )
       ix12 = min0( ix+1,nx+1 )
       ix21 = max0( ix-2,1 )
       ix22 = min0( ix+2,nx+1 )
       work6(ix,0) = 0.d0
       work6(ix,1) = &
            (           ( -work3(ix11,1) ) &
            + 10.d0 * ( -work3(ix,  1) ) & 
            +         ( -work3(ix12,1) ) &
            )
       work8(ix,0) = 0.d0
       work8(ix,1) = &
            (           ( -work4(ix11,1) ) &
            + 10.d0 * ( -work4(  ix,1) ) &
            +         ( -work4(ix12,1) ) &
            )
       work10(ix,-2) = 0.d0
       work10(ix,-1) = 0.d0
       work10(ix,0) = 0.d0
       work10(ix,1)   = (          work3(ix21,1) - 9.d0 * work3(ix11,1) &
            + 3.d0 * work3(  ix,1) + 5.d0 * work3(ix12,1) )
       work12(ix,-1) = 0.d0
       work12(ix,0) = 0.d0
       work12(ix,1) = ( - 5.d0 * work4(ix11,1) - 3.d0 * work4(  ix,1) &
            + 9.d0 * work4(ix12,1) -        work4(ix22,1) )
       work12(ix,2) = ( - 5.d0 * work4(ix11,2) - 3.d0 * work4(  ix,2) &
            + 9.d0 * work4(ix12,2) -        work4(ix22,2) )
       
    enddo
    
    do iz=2,nz
       iz1 = iz + 1
       iz2 = min0( iz+2, nz+1 )
       do ix=2,nx
          work1(ix,-2) = work1(ix,-1)
          work1(ix,-1) = work1(ix,0)
          work1(ix,0) = work1(ix,1)
          work1(ix,1) = ux(ix,iz1) - 2.d0 * ux1(ix,iz1) + ux2(ix,iz1)
          work2(ix,-1) = work2(ix,0)
          work2(ix,0) = work2(ix,1)
          work2(ix,1) = work2(ix,2)
          work2(ix,2) = uz(ix,iz2) - 2.d0 * uz1(ix,iz2) + uz2(ix,iz2)
          work3(ix,-2) = work3(ix,-1)
          work3(ix,-1) = work3(ix,0)
          work3(ix,0) = work3(ix,1)
          work3(ix,1) = work1(ix,1) + 12.d0 * ux1(ix,iz1)
          work4(ix,-1) = work4(ix,0)
          work4(ix,0) = work4(ix,1)
          work4(ix,1) = work4(ix,2)
          work4(ix,2) = work2(ix,2) + 12.d0 * uz1(ix,iz2)
       enddo
       do ix=1,nx+1
          ix11 = max0( ix-1,1 )
          ix12 = min0( ix+1,nx+1 )
          ix21 = max0( ix-2,1 )
          ix22 = min0( ix+2,nx+1 )
          work5(ix) =   (           ( work3(ix11,-1)-work3(ix,-1) ) &
               + 10.d0 * ( work3(ix11, 0)-work3(ix, 0) ) &
               +         ( work3(ix11, 1)-work3(ix, 1) ) )
          work6(ix,0) = work6(ix,1)
          work6(ix,1) =  (  ( work3(ix11,0)-work3(ix11,1) ) &
               + 10.d0 * ( work3(  ix,0)-work3(ix,  1) ) &
               +         ( work3(ix12,0)-work3(ix12,1) ) )
          
          work7(ix) =( ( work4(ix11,-1)-work4(ix,-1) ) &
               + 10.d0 * ( work4(ix11, 0)-work4(ix, 0) ) &
               +         ( work4(ix11, 1)-work4(ix, 1) ))
          
          work8(ix,0) = work8(ix,1)
          work8(ix,1) =  (           ( work4(ix11,0)-work4(ix11,1) ) &
               + 10.d0 * ( work4(  ix,0)-work4(  ix,1) ) &
               +         ( work4(ix12,0)-work4(ix12,1) ))
          
          work9(ix) = (          work3(ix,-2) - 9.d0 * work3(ix,-1) &
               + 3.d0 * work3(ix,0)  + 5.d0 * work3(ix,1))
          
          work10(ix,-2) = work10(ix,-1)
          work10(ix,-1) = work10(ix,0)
          work10(ix,0) = work10(ix,1)
          work10(ix,1) = ( work3(ix21,1) - 9.d0 * work3(ix11,1) &
               + 3.d0 * work3(  ix,1) + 5.d0 * work3(ix12,1) )
          
          work11(ix) = ( - 5.d0 * work4(ix,-1)  - 3.d0 * work4(ix,0) &
               + 9.d0 * work4(ix, 1)  -        work4(ix,2) )
          
          work12(ix,-1) = work12(ix,0)
          work12(ix,0) = work12(ix,1)
          work12(ix,1) = work12(ix,2)
          work12(ix,2) = ( - 5.d0 * work4(ix11,2) - 3.d0 * work4(  ix,2) &
               + 9.d0 * work4(ix12,2) -        work4(ix22,2))
          
       enddo
       
       do ix=2,nx
          ix21 = max0( ix-2,1 )
          ix22 = min0( ix+2,nx+1 )
          ux(ix,iz) = ux(ix,iz) &
               + ( &
               - (           (   work1(ix-1,-1) + work1(ix-1,1) &
               + work1(ix+1,-1) + work1(ix+1,1) ) &
               + 10.d0 * (   work1(ix-1, 0) + work1(  ix,-1) &
               + work1(  ix, 1) + work1(ix+1, 0) ) &
               + 100.d0 * work1(ix,0) ) &
               + e1(ix,iz) * work5(ix) &
               - e2(ix,iz) * work5(ix+1) &
               + e3(ix,iz) * work6(ix,0) &
               - e4(ix,iz) * work6(ix,1) &
               ) / 144.d0 &
               + e13(ix,iz) * work11(ix-1) &
               + e14(ix,iz) * work11(ix) &
               + e15(ix,iz) * work11(ix+1) &
               + e16(ix,iz) * work11(ix22) &
               + e17(ix,iz) * work12(ix,-1) &
               + e18(ix,iz) * work12(ix,0) &
               + e19(ix,iz) * work12(ix,1) &
               + e20(ix,iz) * work12(ix,2)
          uz(ix,iz) = uz(ix,iz) &
               + ( &
               - (           (   work2(ix-1,-1) + work2(ix-1,1) &
               + work2(ix+1,-1) + work2(ix+1,1) ) &
               + 10.d0 * (   work2(ix-1, 0) + work2(  ix,-1) &
               + work2(  ix, 1) + work2(ix+1, 0) ) &
               + 100.d0 * work2(ix,0) ) &
               + f1(ix,iz) * work7(ix) &
               - f2(ix,iz) * work7(ix+1) &
               + f3(ix,iz) * work8(ix,0) &
               - f4(ix,iz) * work8(ix,1) &
               ) / 144.d0 &
               + f13(ix,iz) * work9(ix21) &
               + f14(ix,iz) * work9(ix-1) &
               + f15(ix,iz) * work9(ix) &
               + f16(ix,iz) * work9(ix+1) &
               + f17(ix,iz) * work10(ix,-2) &
               + f18(ix,iz) * work10(ix,-1) &
               + f19(ix,iz) * work10(ix,0) &
               + f20(ix,iz) * work10(ix,1)
       enddo
    enddo
     ux(isx,isz) = ux(isx,isz) + fx(isx,isz)
     uz(isx,isz) = uz(isx,isz) + fz(isx,isz)     
 endif




 ! swapping u1 & u2 
 do iz=2,nz
    do ix=2,nx
       ux2(ix,iz) = ux1(ix,iz)
       ux1(ix,iz) =  ux(ix,iz)
       uz2(ix,iz) = uz1(ix,iz)
       uz1(ix,iz) =  uz(ix,iz)
    enddo
 enddo


 
 return
end subroutine calstep
	

