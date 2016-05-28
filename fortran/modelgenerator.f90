program generator
 implicit none

 integer, parameter     :: prec = kind (1.0) 
 integer                :: i,j
 integer                :: level,tmp


 integer, parameter                             :: nbre_layers = 5
 integer, dimension(nbre_layers)                :: layer_thikness
 real(kind=prec), dimension(nbre_layers)        :: vp
 real(kind=prec), dimension(nbre_layers)        :: vs
 real(kind=prec), allocatable                   :: fullvp (:,:)
 real(kind=prec), allocatable                   :: fullvs (:,:)
 real(kind=prec), allocatable                   :: fullrho(:,:)
 real(kind=prec)                                :: dval,pval,xover1,xover2,grad,const
 real(kind=prec), parameter                     :: vpwater = 1.5e3
 integer                                        :: NZ, NX
 integer                                        :: NX_TOTAL, NZ_TOTAL
 integer                                        :: recl_size
!************************************************************************	

! PROGRAM GENERATES VELOCITY MODEL (VP & VS)

!************************************************************************
! Model initialization	
 NX_TOTAL = 400
 NZ_TOTAL = 200

!*** Thikness, Vp, Vs *******************************************		

 layer_thikness(1)=40;vp(1)=2200;vs(1)=1400;
 layer_thikness(2)=90;vp(2)=2300;vs(2)=1450;
 layer_thikness(3)=140;vp(3)=2500;vs(3)=1550;
 layer_thikness(4)=180;vp(4)=2700;vs(4)=1700;
 layer_thikness(5)=200;vp(5)=3000;vs(5)=1900;

! layer_thikness(1)=40;vp(1)=2200;vs(1)=1400;
! layer_thikness(2)=90;vp(2)=2200;vs(2)=1400;
! layer_thikness(3)=140;vp(3)=2200;vs(3)=1400;
! layer_thikness(4)=180;vp(4)=2200;vs(4)=1400;
! layer_thikness(5)=200;vp(5)=2200;vs(5)=1400;
 



 allocate (fullvp (NX_TOTAL, NZ_TOTAL) )
 allocate (fullvs (NX_TOTAL, NZ_TOTAL) )
 allocate (fullrho(NX_TOTAL, NZ_TOTAL) )
 recl_size = prec * NX_TOTAL * NZ_TOTAL

!****************************************************************	   
 
 open (1,file='./2d_homo.vp',form='unformatted',access='direct',recl=recl_size)
 open (2,file='./2d_homo.vs',form='unformatted',access='direct',recl=recl_size) 
 open (3,file='./2d_homo.rho',form='unformatted',access='direct',recl=recl_size)

 tmp=1
 do j=1, NZ_TOTAL
   level=tmp
   if (j .gt. layer_thikness(level)) level=level+1
   tmp=level
   do i=1, NX_TOTAL
     fullvp(i,j) = vp(level)
     fullvs(i,j) = vs(level)
   enddo
 enddo

 write(1,rec=1) fullvp(:,:)
 write(2,rec=1) fullvs(:,:)


 do j=1,NZ_TOTAL
     do i=1,NX_TOTAL
        pval = fullvp(i,j)
        if (pval.le.vpwater) then
           dval = 1000.E0
        elseif (pval > vpwater .and. pval < 2000.E0) then
           dval = 2351.E0-(7497.E0)*(pval/1000.E0)**(-4.656E0)
        elseif (pval >= 2000.E0 .and. pval <= 2150.E0) then
           xover1 = 2351.E0-(7497.E0)*(2000.E0/1000.E0)**(-4.656E0)
           xover2 = 1740.E0*(2150.E0/1000.E0)**(0.25E0)
           grad = 150.E0/(xover2-xover1)
           const = 2000.E0-(xover1*grad)
           dval = (pval-const)/grad
        elseif (pval > 2150) then
           dval = 1740.E0*(pval/1000.E0)**(0.25E0)
        endif
        fullrho(i,j) = dval
     enddo
  enddo
  write(3,rec=1) fullrho(:,:)

 close (1,status='keep')
 close (2,status='keep')
 close (3,status='keep')
end program generator
