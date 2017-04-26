program generator

  ! Shihao Yuan and Nobuaki Fuji 2016

  
 implicit none

 integer, parameter     :: prec = kind (1.e0) 
 integer                :: i,j
 integer                :: level,tmp

 integer, parameter :: times = 10

 integer, parameter                             :: nbre_layers = 5
 integer, dimension(nbre_layers)                :: layer_thikness
 real(kind=prec), dimension(nbre_layers)        :: vp
 real(kind=prec), dimension(nbre_layers)        :: vs
 real(kind=prec), dimension(nbre_layers)        :: rho
 real(kind=prec), allocatable                   :: fullvp (:,:)
 real(kind=prec), allocatable                   :: fullvs (:,:)
 real(kind=prec), allocatable                   :: fullrho(:,:),tmpM(:,:)
 real(kind=prec)                                :: dval,pval,xover1,xover2,grad,const
 real(kind=prec), parameter                     :: vpwater = 1.5e3
 integer                                        :: NZ, NX
 integer ::ix,iz,tmpint
 double precision :: boundary,dx,tmpfloat
 integer                                        :: NX_TOTAL, NZ_TOTAL
 integer                                        :: recl_size
!************************************************************************	

! PROGRAM GENERATES VELOCITY MODEL (VP & VS)

!************************************************************************
! Model initialization	
 NX_TOTAL = (1589-1)*times+1
 NZ_TOTAL = (1601-1)*times+1

 print *, NX_TOTAL, NZ_TOTAL

 !NX_TOTAL = (601-1)*times+1
 !NZ_TOTAL = (161-1)*times+1

!*** Thikness, Vp, Vs *******************************************		

 layer_thikness(1)=501;vp(1)=10000.d0;vs(1)=5000.d0;rho(1)=4000.d0;
 layer_thikness(2)=801;vp(2)=8000.d0;vs(2)=3000.d0;rho(2)=3000.d0;
 layer_thikness(3)=1101;vp(3)=10000.d0;vs(3)=5000.d0;rho(3)=4000.d0;
 layer_thikness(4)=1401;vp(4)=8000.d0;vs(4)=3000.d0;rho(4)=3000.d0;
 layer_thikness(5)=1601;vp(5)=10000.d0;vs(5)=5000.d0;rho(5)=4000.d0;

 allocate (fullvp (NX_TOTAL, NZ_TOTAL) )
 allocate (fullvs (NX_TOTAL, NZ_TOTAL) )
 allocate (fullrho(NX_TOTAL, NZ_TOTAL) )

 recl_size = prec * NX_TOTAL * NZ_TOTAL

print *, recl_size, prec

!****************************************************************	   
 
 open (1,file='../models/2d_anais.vp',form='unformatted',access='direct',recl=recl_size)
 open (2,file='../models/2d_anais.vs',form='unformatted',access='direct',recl=recl_size) 
 open (3,file='../models/2d_anais.rho',form='unformatted',access='direct',recl=recl_size)




 if(0.eq.0) then
    tmp=1
    do j=1, NZ_TOTAL
       level=tmp
       if (j .gt. layer_thikness(level)) level=level+1
       tmp=level
       do i=1, NX_TOTAL
          fullvp(i,j) = vp(level)
          fullvs(i,j) = vs(level)
          fullrho(i,j) = rho(level)
       enddo
    enddo
 endif
 write(1,rec=1) fullvp(:,:)
 write(2,rec=1) fullvs(:,:)
 write(3,rec=1) fullrho(:,:)


 close (1,status='keep')
 close (2,status='keep')
 close (3,status='keep')
end program generator
