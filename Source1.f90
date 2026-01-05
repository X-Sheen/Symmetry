! Particle cluster & RP instibility
! Color_gradient model
! Wetting boundary condition by Takashi Akai 2018  
! 1/8 symmetric about X,Y,Z
! From Xu Zhiyuan
! 21-JAN-2021
  module const
implicit none
integer, parameter:: xl=152
integer, parameter:: yl=152
integer, parameter:: zl=152
double precision, dimension(0:18):: wi
integer, dimension(0:18):: cix, ciy, ciz, opc,symx,symy,symz
double precision:: sigma, Rad, Rad_c
double precision:: nu_r, nu_b, rho0
double precision, parameter:: EVAP_LIM = 1.0d-8
double precision, parameter:: beta = 0.7d0   !!!!!!界面厚度
real*8, dimension(0:8):: weight8(0:8)=(/0.0d0, 4.0d0/45.0d0, 1.0d0/21.0d0, 2.0d0/105.0d0, 5.0d0/504.0d0, 1.0d0/315.0d0, &
                        &1.0d0/630.0d0, 0.0d0, 1.0d0/5040.0d0/)
real*8, parameter:: PI=4.0d0*datan(1.0d0)
real*8, parameter:: cita=90.0d0/180.0d0*PI       !!!Contact angle
real*8, parameter:: alpha=2.0d0/180.0d0*pi       !!!Related to the total number of Lagrange points
real*8 :: gravity
integer, parameter:: par=12
integer:: La_m,La_n,La_num
end module const

program main
use const
implicit none
integer t,i, j, x, y, z
double precision, dimension(0:18,0:18):: M, M_inv
integer, allocatable, dimension(:,:,:):: walls
integer, allocatable, dimension(:,:,:):: walls_old,walls_all,pin,pout
integer, allocatable, dimension(:,:):: next_x, next_y, next_z
double precision, allocatable, dimension(:,:,:):: rho, ux, uy, uz
double precision, allocatable, dimension(:,:,:):: rho_N, rho_R, rho_B
double precision, allocatable, dimension(:,:,:):: n_x, n_y, n_z
double precision, allocatable, dimension(:,:,:):: gradx, grady, gradz
double precision, allocatable, dimension(:,:,:,:):: fi,ri,bi,ffi
logical, allocatable, dimension(:,:,:):: flag
double precision,allocatable,dimension(:,:) :: Lax,Lay,Laz,La_rhoN,La_nx,La_ny,La_nz,ds
logical,allocatable,dimension(:,:) :: La_flag

real*8 :: center(par,3),center_old(par,3),center0(par,3)
real*8 :: Fpx(par),Fpy(par),Fpz(par),Tpx(par),Tpy(par),Tpz(par),Fg
real*8 :: Frex(par),Frey(par),Frez(par)
real*8 :: tension(par,3),tension_t(par,3)
real*8 :: upx(par),upy(par),upz(par),rpx(par),rpy(par),rpz(par)
real*8 :: rho_p,mass,inertia_moment
real*8 :: sum_ri,cita_star,error,dis_c
real*8 :: tau_R,tau_B
real*8 :: Oh,k_star,lamda_d,lamda
real*8 :: tr,tmax,tsmall
real*8 :: Wcap_old
real*8 :: Shb(par)
!set parameters
Oh=0.1d0       !1.0d0!0.1d0  

lamda=10.0d0   ! viscosity ratio ratio

lamda_d=(real(yl)-2)*2.0d0
Rad_c=30.0d0

! Fluid
sigma=0.05d0
rho0=0.05d0!0.05d0!0.1d0

nu_R=Oh/sqrt(rho0/(Rad_c*sigma))
nu_B=nu_R/lamda
tau_R=3.0d0*nu_R+0.5d0
tau_B=3.0d0*nu_B+0.5d0
!rho0=Oh**2*sigma*Rad_c/nu_R**2
write(*,*)'nu_R=',nu_R
! Particle
Rad=7.5d0!10.0d0
gravity=0.00098

rho_p=1.0d0*rho0!100.0d0*rho0  !!!!!!!!!!!!!!! density of particle

mass=4.0d0/3.0d0*rho_p*PI*Rad**3
inertia_moment=2.0d0/5.0d0*mass*rad**2   
!Fg=-(rho_p)*mass*gravity
Fg=-pi*2.0d0*Rad*sigma

! Lagrange points
  La_m=2.0d0*pi/alpha  
  La_n=pi/alpha
  La_num=La_m*(La_n+1)

! Output
tr=sqrt(Rad_c**3*rho0/sigma)
tmax=tr*200.0d0  !tr*300.0d0               !!!!!!
tsmall=tr*30.0d0 !tr*100.0d0               !!!!!!
write(*,*)'tmax=',tmax
write(*,*)'tsmall=',tsmall
write(*,*)'tr=',tr
! 验算参数
k_star=2.0d0*PI*Rad_c/lamda_d
Oh=nu_R*rho0/sqrt(rho0*sigma*Rad_c)
write(*,*)'k_star=',k_star
write(*,*)'lamda_d=',lamda_d
write(*,*)'Oh=',Oh

data M /1,-30,12,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
		1,-11,-4,1,-4,0,0,0,0,2,-4,0,0,0,0,0,0,0,0, &
		1,-11,-4,-1,4,0,0,0,0,2,-4,0,0,0,0,0,0,0,0, &
		1,-11,-4,0,0,1,-4,0,0,-1,2,1,-2,0,0,0,0,0,0, &
		1,-11,-4,0,0,-1,4,0,0,-1,2,1,-2,0,0,0,0,0,0, &
		1,-11,-4,0,0,0,0,1,-4,-1,2,-1,2,0,0,0,0,0,0, &
		1,-11,-4,0,0,0,0,-1,4,-1,2,-1,2,0,0,0,0,0,0, &
		1,8,1,1,1,1,1,0,0,1,1,1,1,1,0,0,1,-1,0, &
		1,8,1,-1,-1,1,1,0,0,1,1,1,1,-1,0,0,-1,-1,0, &
		1,8,1,1,1,-1,-1,0,0,1,1,1,1,-1,0,0,1,1,0, &
		1,8,1,-1,-1,-1,-1,0,0,1,1,1,1,1,0,0,-1,1,0, & 
		1,8,1,1,1,0,0,1,1,1,1,-1,-1,0,0,1,-1,0,1, &
		1,8,1,-1,-1,0,0,1,1,1,1,-1,-1,0,0,-1,1,0,1, &
		1,8,1,1,1,0,0,-1,-1,1,1,-1,-1,0,0,-1,-1,0,-1, & 
		1,8,1,-1,-1,0,0,-1,-1,1,1,-1,-1,0,0,1,1,0,-1, &
		1,8,1,0,0,1,1,1,1,-2,-2,0,0,0,1,0,0,1,-1, & 
		1,8,1,0,0,-1,-1,1,1,-2,-2,0,0,0,-1,0,0,-1,-1, &
		1,8,1,0,0,1,1,-1,-1,-2,-2,0,0,0,-1,0,0,1,1, & 
		1,8,1,0,0,-1,-1,-1,-1,-2,-2,0,0,0,1,0,0,-1,1/

!inverse of M
do i = 0, 18
	do j = 0, 18
		M_inv(i,j) = M(j,i) / DOT_PRODUCT(M(j,:),M(j,:))
	end do !j
end do	!i	

cix = (/0,1,-1,0,0,0,0,1,-1,1,-1,1,-1,1,-1,0,0,0,0/)
ciy = (/0,0,0,1,-1,0,0,1,1,-1,-1,0,0,0,0,1,-1,1,-1/)
ciz = (/0,0,0,0,0,1,-1,0,0,0,0,1,1,-1,-1,1,1,-1,-1/)
opc = (/0,2,1,4,3,6,5,10,9,8,7,14,13,12,11,18,17,16,15/)
symx= (/0,2,1,3,4,5,6,8,7,10,9,12,11,14,13,15,16,17,18/)! 关于YOZ平面对称
symy= (/0,1,2,4,3,5,6,9,10,7,8,11,12,13,14,16,15,18,17/)! 关于XOZ平面对称 
symz= (/0,1,2,3,4,6,5,7,8,9,10,13,14,11,12,17,18,15,16/)! 关于XOY平面对称 

wi(0) = 1.0d0/3.0d0
wi(1:6) = 1.0d0/18.0d0
wi(7:18) = 1.0d0/36.0d0

! allocate memory
allocate(next_x(1:18,1:xl))
allocate(next_y(1:18,1:yl))
allocate(next_z(1:18,1:zl))
allocate(walls(zl,yl,xl))
allocate(walls_all(zl,yl,xl))
allocate(walls_old(zl,yl,xl))
allocate(pin(zl,yl,xl))
allocate(pout(zl,yl,xl))
allocate(rho(zl,yl,xl))
allocate(rho_r(zl,yl,xl))
allocate(rho_b(zl,yl,xl))
allocate(rho_N(zl,yl,xl))
allocate(ux(zl,yl,xl))
allocate(uy(zl,yl,xl))
allocate(uz(zl,yl,xl))
allocate(flag(zl,yl,xl))
allocate(gradx(zl,yl,xl))
allocate(grady(zl,yl,xl))
allocate(gradz(zl,yl,xl))
allocate(n_x(zl,yl,xl))
allocate(n_y(zl,yl,xl))
allocate(n_z(zl,yl,xl))
allocate(fi(0:18,zl,yl,xl))
allocate(ffi(0:18,zl,yl,xl))
allocate(ri(0:18,zl,yl,xl))
allocate(bi(0:18,zl,yl,xl))
allocate(Lax(par,La_num))
allocate(Lay(par,La_num))
allocate(Laz(par,La_num))
allocate(La_rhoN(par,La_num))
allocate(La_nx(par,La_num))
allocate(La_ny(par,La_num))
allocate(La_nz(par,La_num))
allocate(La_flag(par,La_num))
allocate(ds(par,La_num))

walls=0
walls_all=0
pin=0
pout=0
upx=0.0d0
upy=0.0d0
upz=0.0d0
rpx=0.0d0
rpy=0.0d0
rpz=0.0d0
Fpx=0.0d0
Fpy=0.0d0
Fpz=0.0d0
Tpx=0.0d0
Tpy=0.0d0
Tpz=0.0d0
Frex=0.0d0
Frey=0.0d0
Frez=0.0d0
tension=0.0d0
tension_t=0.0d0
Wcap_old=0.0d0
open(unit=400,file='parameter.dat',status='unknown')
write(400,*)'VARIABLES=k_star,Oh,cita,Rad_c,lamda_d,sigma,nu_R,rho0,Rad,par,tr'
write(400,"(11(1xF12.6))")k_star,Oh,cita,Rad_c,lamda_d,sigma,nu_R,rho0,Rad,par*1.0,tr
close(400)

call tables(next_x, next_y, next_z)
call initia_center(center,center0,lamda_d)   
!center=center*200.0d0
!center0=center
call set_walls(walls,walls_old,walls_all,pin,pout,center,next_x,next_y,next_z)
call initialization(lamda_d,rho,rho_N,rho_R,rho_B,ux,uy,uz,fi,ffi,ri,bi)
call derivatives(walls_all,pin,pout,center,next_x,next_y,next_z,rho_N,gradx,grady,gradz,n_x,n_y,n_z,flag)
do t = 0, int(tmax)
	call calmacros(fi,ri,bi,rho,ux,uy,uz,rho_R,rho_B,rho_N,walls_all)

    if(t<=int(tsmall).and.mod(t,20)==0) then
        write(*,*)t,t/tr,sum(Shb)
         call output_result2(center,t,tr,upx,upy,upz,rpx,rpy,rpz,Fpx,Fpy,Fpz,Tpx,Tpy,Tpz,tension,tension_t,Frex,Frey,Frez)      
         call output3(rho_N,t,tr,walls_all,pout,lamda_d)
         call output4(rho_N,t,tr,walls_all,pout,center0,ux,uy,uz,rho,upx,upy,upz,rpx,rpy,rpz,tension,tension_t,mass,inertia_moment,Wcap_old,Shb)
    end if

    if(t<=int(tsmall).and.mod(t,100)==0) then
      call output_results(t,walls_all,rho_N,ux,uy,uz,rho)
      
    elseif(mod(t,5000)==0) then
      call output_results(t,walls_all,rho_N,ux,uy,uz,rho)
      call output_result2(center,t,tr,upx,upy,upz,rpx,rpy,rpz,Fpx,Fpy,Fpz,Tpx,Tpy,Tpz,tension,tension_t,Frex,Frey,Frez)
      call output3(rho_N,t,tr,walls_all,pout,lamda_d)
      call output4(rho_N,t,tr,walls_all,pout,center0,ux,uy,uz,rho,upx,upy,upz,rpx,rpy,rpz,tension,tension_t,mass,inertia_moment,Wcap_old,Shb)
      write(*,*)t,t/tr,sum(Shb)
    end if
    
!if(1==2)then   
      call Lagrange_points(walls,center,rho_N,gradx,grady,gradz,n_x,n_y,n_z,flag,Lax,Lay,Laz,La_rhoN,La_nx,La_ny,La_nz,La_flag,ds,Shb)
      call capillary_force2(center,Lax,Lay,Laz,La_nx,La_ny,La_nz,La_rhoN,La_flag,ds,tension,tension_t)
      call momentum(walls,pout,center,Fpx,Fpy,Fpz,Tpx,Tpy,Tpz,upx,upy,upz,rpx,rpy,rpz,ffi,fi)
      call repulsive(center,Frex,Frey,Frez,Fg)
      call movement(center,center_old,center0,Fpx,Fpy,Fpz,Frex,Frey,Frez,Tpx,Tpy,Tpz,upx,upy,upz,rpx,rpy,rpz,Fg,mass,inertia_moment,tension,tension_t)
      call set_walls(walls,walls_old,walls_all,pin,pout,center,next_x,next_y,next_z)
      call refilling(next_x,next_y,next_z,pout,walls_old,center_old,upx,upy,upz,rpx,rpy,rpz,rho,rho_N,rho_R,rho_B,ux,uy,uz,ri,bi,fi) 
!end if
    call derivatives(walls_all,pin,pout,center,next_x,next_y,next_z,rho_N,gradx,grady,gradz,n_x,n_y,n_z,flag)
	call collision(walls_all,pout,next_x,next_y,next_z,flag,n_x,n_y,n_z,gradx,grady,gradz,rho,rho_N,ux,uy,uz,fi,ffi,M,M_inv)
	call Recolor(flag,fi,ri,bi,rho_R,rho_B,rho,n_x,n_y,n_z,walls_all)
	call boundary(next_x,next_y,next_z,ri,bi)
    call streaming(ri, next_x, next_y, next_z)
	call streaming(bi, next_x, next_y, next_z)
    call curve(walls,pout,ri,bi,upx,upy,upz,rpx,rpy,rpz,center,rho_R,rho_B)  
end do !t
    end program main
      
subroutine angle_test(t,rho_N,center,dis_c,error,cita_star)
use const
implicit none
real*8, intent(in):: rho_N(zl,yl,xl)
real*8 :: center(par,3),dis_c,cita_star,error
real*8 :: dis,z0,a
integer:: t,x,y,z,zn

x=(1+xl)/2
y=(1+yl)/2
do z=int(center(1,3)+Rad+3.0),zl
  zn=z+1
    if(rho_N(z,y,x)*rho_N(zn,y,z)<0.0d0) then   !!! 跨界面
      z0=dble(z)-rho_N(z,y,x)/(rho_N(zn,y,x)-rho_N(z,y,x))
      cycle
    end if
end do

dis=z0-center(1,3)-Rad
error=(dis-dis_c)/dis_c
a=asin(Rad/dis*sin(cita))
cita_star=2.0d0*a/PI*180.0d0


if(isNAN(cita_star))then
  cita_star=0.0d0
end if

!write(*,*)'a=',a,dis,sin(cita),Rad/dis*sin(cita)
open(unit=100, file='cita.dat',status='unknown')
if(t==0)then
  write(100,*)'TITLE="Contact_Angles"'
  write(100,*)'VARIABLES=t,cita_star,error'
end if

write(100,"(1(1xI4),2(1xF12.8))") t,cita_star,error

return
end subroutine angle_test

subroutine tables(next_x, next_y, next_z)
use const, only: xl, yl, zl, cix, ciy, ciz
implicit none
integer, dimension(18,xl), intent(OUT):: next_x
integer, dimension(18,yl), intent(OUT):: next_y
integer, dimension(18,zl), intent(OUT):: next_z
integer i, j
do i = 1, xl
        do j = 1, 18
                next_x(j,i) = i + cix(j)
                if (next_x(j,i) > xl) next_x(j,i) = 1
                if (next_x(j,i) < 1)  next_x(j,i) = xl
        end do !j
end do !i
do i = 1, yl
        do j = 1, 18
                next_y(j,i) = i + ciy(j)
                if (next_y(j,i) > yl) next_y(j,i) = 1
                if (next_y(j,i) < 1)  next_y(j,i) = yl
        end do !j
end do !i
do i = 1, zl
        do j = 1, 18
                next_z(j,i) = i + ciz(j)
                if (next_z(j,i) > zl) next_z(j,i) = 1
                if (next_z(j,i) < 1)  next_z(j,i) = zl
        end do !j
end do !i
return
end subroutine tables

subroutine set_walls(walls,walls_old,walls_all,pin,pout,center,next_x,next_y,next_z)
use const
implicit none
integer:: next_x(18,xl), next_y(18,yl), next_z(18,zl)
integer:: walls(zl,yl,xl),walls_old(zl,yl,xl),walls_all(zl,yl,xl),pin(zl,yl,xl),pout(zl,yl,xl)
integer:: x,y,z,xf,xxf,yf,yyf,zf,i,pp
real*8 :: distance,center(par,3)
walls_old=walls_all
walls_all=0
walls=0
pin=0
pout=0

!$OMP PARALLEL DO PRIVATE(distance), SCHEDULE(GUIDED)
do pp=1,par
  do x=1,xl
    do y=1,yl
      do z=1,zl
        distance= (dble(x)-center(pp,1))**2+(dble(y)-center(pp,2))**2+(dble(z)-center(pp,3))**2
        if(distance<=Rad**2)then
          walls(z,y,x)=pp ! Solid
          walls_all(z,y,x)=10
        end if
      end do
    end do
  end do
end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(xf,yf,zf), SCHEDULE(GUIDED)
do pp=1,par
  do x=2,xl-1
    do y=2,yl-1
      do z=2,zl-1
        if(walls_all(z,y,x)==0)then
          do i=1,18
            xf=next_x(i,x)           
            yf=next_y(i,y)       
            zf=next_z(i,z)
            if(walls(zf,yf,xf)==pp)then ! Forward to solid point, particle pp
              pout(z,y,x)=pp         ! Live fluid point
              pin(zf,yf,xf)=pp       ! Live solid point
            end if
          end do
        end if
      end do
    end do
  end do
end do
!$OMP END PARALLEL DO          
walls(:,:,xl)=walls(:,:,xl-1)
walls_all(:,:,xl)=walls_all(:,:,xl-1)
pin(:,:,xl)=pin(:,:,xl-1)
pout(:,:,xl)=pout(:,:,xl-1)

walls(:,yl,:)=walls(:,yl-1,:)
walls_all(:,yl,:)=walls_all(:,yl-1,:)
pin(:,yl,:)=pin(:,yl-1,:)
pout(:,yl,:)=pout(:,yl-1,:)

walls(:,1,:)=walls(:,2,:)
walls_all(:,1,:)=walls_all(:,2,:)
pin(:,1,:)=pin(:,2,:)
pout(:,1,:)=pout(:,2,:)

walls(zl,:,:)=walls(zl-1,:,:)
walls_all(zl,:,:)=walls_all(zl-1,:,:)
pin(zl,:,:)=pin(zl-1,:,:)
pout(zl,:,:)=pout(zl-1,:,:)
return
end subroutine set_walls

subroutine equilf(meq, fk, rho, ux, uy, uz, fx, fy, fz)
implicit none
double precision:: usqr
double precision, intent(IN):: rho, ux, uy, uz
double precision, intent(IN):: fx, fy, fz
double precision:: uf
double precision, intent(OUT):: meq(0:18), fk(0:18)

usqr = ux*ux+uy*uy+uz*uz
meq(0) = 1.0d0
meq(1) = 19.0d0*usqr-11.0d0
meq(2) = 3.0d0-5.5d0*usqr
meq(3) = ux
meq(4) = -ux/1.5d0
meq(5) = uy
meq(6) = -uy/1.5d0
meq(7) = uz
meq(8) = -uz/1.5d0
meq(9) = 3.0d0*ux*ux-usqr
meq(10) = -0.5d0*(3.0d0*ux*ux-usqr)
meq(11) = uy*uy-uz*uz
meq(12) = -0.5d0*(uy*uy-uz*uz)
meq(13) = ux*uy
meq(14) = uy*uz
meq(15) = ux*uz
meq(16) = 0.0d0
meq(17) = 0.0d0
meq(18) = 0.0d0

meq = rho*meq

uf = ux*fx+uy*fy+uz*fz
fk(0) = 0.0d0
fk(1) = 38.0d0*uf
fk(2) = -11.0d0*uf
fk(3) = fx
fk(4) = -fx/1.5d0
fk(5) = fy
fk(6) = -fy/1.5d0
fk(7) = fz
fk(8) = -fz/1.5d0
fk(9) = 6.0d0*ux*fx-2.0d0*uf
fk(10) = -(3.0d0*ux*fx-uf)
fk(11) = 2.0d0*(uy*fy-uz*fz)
fk(12) = -(uy*fy-uz*fz)
fk(13) = ux*fy+uy*fx
fk(14) = uy*fz+uz*fy
fk(15) = ux*fz+uz*fx
fk(16) = 0.0d0
fk(17) = 0.0d0
fk(18) = 0.0d0

return
end subroutine equilf

subroutine derivatives(walls_all,pin,pout,center,next_x,next_y,next_z,rho_N,gradx,grady,gradz,n_x,n_y,n_z,flag)
use const
implicit none
integer x, y, z, i
integer xn, yn, zn
double precision:: val, fnorm
integer :: walls_all(zl,yl,xl),pin(zl,yl,xl),pout(zl,yl,xl)
double precision :: rho_N(zl,yl,xl),gradx(zl,yl,xl),grady(zl,yl,xl),gradz(zl,yl,xl)
double precision :: n_x(zl,yl,xl), n_y(zl,yl,xl), n_z(zl,yl,xl)
integer :: next_x(18,xl), next_y(18,yl), next_z(18,zl)
logical, dimension(zl,yl,xl) :: flag
real*8 :: center(par,3),nw(3), n_temp(3)
integer :: pp
real*8:: wi_sum, coeff(3), z_check(2),theta
real*8:: temp, sol(2), a(2), b(2), vec1(3), vec2(3), db1, db2
! 补充固体内部 rho_N 值
!$OMP PARALLEL DO PRIVATE(wi_sum, xn, yn, zn), SCHEDULE(GUIDED)
do x=2,xl-1
    do y=2,yl-1
        do z=2,zl-1
            if(pin(z,y,x)>0) then ! Solid
                wi_sum=0.0d0
                rho_N(z,y,x)=0.0d0
                do i=1,18
                    xn=next_x(i,x)
                    yn=next_y(i,y)
                    zn=next_z(i,z)
                    if(walls_all(zn,yn,xn)<=0) then ! Fluid
                        wi_sum=wi_sum+wi(i)
                        rho_N(z,y,x)=rho_N(z,y,x)+wi(i)*rho_N(zn,yn,xn)
                    end if
                end do
                rho_N(z,y,x)=rho_N(z,y,x)/wi_sum
                if(wi_sum==0.0)then
                  write(*,*)'warnning!!! devision wi_sum=0',x,y,z
                end if
            end if
        end do
    end do
end do
!$OMP END PARALLEL DO
call symmetry_X_phi(rho_N)
call symmetry_Y_phi(rho_N)
call symmetry_Z_phi(rho_N)

!计算导数及界面法矢量
!$OMP PARALLEL DO PRIVATE(xn, yn, zn, val, fnorm), SCHEDULE(GUIDED)
do x = 2, xl-1
	do y = 2, yl-1
		do z = 2, zl-1
            gradx(z,y,x)=0.0d0
            grady(z,y,x)=0.0d0
            gradz(z,y,x)=0.0d0
            n_x(z,y,x)=0.0d0
            n_y(z,y,x)=0.0d0
            n_z(z,y,x)=0.0d0
            flag(z,y,x)=.false.   
			if (walls_all(z,y,x) <= 0) then     ! for fluid nodes
				do i = 1, 18
                    xn = next_x(i,x)
					yn = next_y(i,y)
					zn = next_z(i,z)                    
					val = wi(i) * rho_N(zn,yn,xn)
					gradx(z,y,x) = gradx(z,y,x) + val * cix(i)
					grady(z,y,x) = grady(z,y,x) + val * ciy(i)
					gradz(z,y,x) = gradz(z,y,x) + val * ciz(i)
				end do !i
				gradx(z,y,x) = 3.0d0 * gradx(z,y,x)
				grady(z,y,x) = 3.0d0 * grady(z,y,x)
				gradz(z,y,x) = 3.0d0 * gradz(z,y,x)
				fnorm = dsqrt(gradx(z,y,x)**2 + grady(z,y,x)**2 + gradz(z,y,x)**2)
				if (fnorm > EVAP_LIM) then
					n_x(z,y,x) = gradx(z,y,x) / fnorm
					n_y(z,y,x) = grady(z,y,x) / fnorm
					n_z(z,y,x) = gradz(z,y,x) / fnorm
					flag(z,y,x) = .true.
				endif				
			endif		
		end do !z
	end do !y
end do !x
!$OMP END PARALLEL DO

! wetting boundary condition
!$OMP PARALLEL DO PRIVATE(pp,nw,n_temp,theta,a,b,vec1,vec2,db1,db2,fnorm), SCHEDULE(GUIDED)
do x=2,xl-1
    do y=2,yl-1
        do z=2,zl-1
          if(pout(z,y,x)>=1.and.flag(z,y,x))then
            pp=pout(z,y,x)
            nw(1)=dble(x)-center(pp,1)
            nw(2)=dble(y)-center(pp,2)
            nw(3)=dble(z)-center(pp,3)
            call unitization(nw)
            n_temp=(/n_x(z,y,x), n_y(z,y,x), n_z(z,y,x)/)
            if (  maxval ( abs( nw(:)-n_temp(:) ) ) <1.0d-5 ) then
                theta=1.0d-5                                               !!!! ???
            else
                theta=dacos(dot_product(nw, n_temp))
            endif
            
            a(1)=dcos(cita)-dsin(cita)*dcos(theta)/dsin(theta)
            a(2)=dcos(cita)+dsin(cita)*dcos(theta)/dsin(theta)
            b(1)=dsin(cita)/dsin(theta)
            b(2)=-dsin(cita)/dsin(theta)
           
            vec1=a(1)*nw+b(1)*(/n_x(z,y,x), n_y(z,y,x), n_z(z,y,x)/)
            vec2=a(2)*nw+b(2)*(/n_x(z,y,x), n_y(z,y,x), n_z(z,y,x)/)           
                       
            db1=(vec1(1)-n_x(z,y,x))**2+(vec1(2)-n_y(z,y,x))**2+(vec1(3)-n_z(z,y,x))**2
            db2=(vec2(1)-n_x(z,y,x))**2+(vec2(2)-n_y(z,y,x))**2+(vec2(3)-n_z(z,y,x))**2
                        
            if(db1<db2) then
                n_x(z,y,x)=vec1(1)
                n_y(z,y,x)=vec1(2)
                n_z(z,y,x)=vec1(3)
                
                fnorm = dsqrt(gradx(z,y,x)**2 + grady(z,y,x)**2 + gradz(z,y,x)**2)
                gradx(z,y,x)=fnorm*n_x(z,y,x)
                grady(z,y,x)=fnorm*n_y(z,y,x)
                gradz(z,y,x)=fnorm*n_z(z,y,x)
            else
                n_x(z,y,x)=vec2(1)
                n_y(z,y,x)=vec2(2)
                n_z(z,y,x)=vec2(3)
                
                fnorm = dsqrt(gradx(z,y,x)**2 + grady(z,y,x)**2 + gradz(z,y,x)**2)
                gradx(z,y,x)=fnorm*n_x(z,y,x)
                grady(z,y,x)=fnorm*n_y(z,y,x)
                gradz(z,y,x)=fnorm*n_z(z,y,x)
            end if
        end if
        end do
    end do
end do
!$OMP END PARALLEL DO

call symmetry_X_u(n_x,n_y,n_z)
call symmetry_Y_u(n_x,n_y,n_z)
call symmetry_Z_u(n_x,n_y,n_z)
call symmetry_X_u(gradx,grady,gradz)
call symmetry_Y_u(gradx,grady,gradz)
call symmetry_Z_u(gradx,grady,gradz)

return
end subroutine derivatives

subroutine collision(walls_all,pout,next_x,next_y,next_z,flag,n_x,n_y,n_z,gradx,grady,gradz,rho,rho_N,ux,uy,uz,fi,ffi,M,M_inv)
use const
implicit none
integer x, y, z, i
integer xn, yn, zn
integer :: walls_all(zl,yl,xl),pout(zl,yl,xl)
integer :: next_x(18,xl), next_y(18,yl), next_z(18,zl)
logical :: flag(zl,yl,xl)
double precision :: n_x(zl,yl,xl), n_y(zl,yl,xl), n_z(zl,yl,xl)
double precision :: gradx(zl,yl,xl), grady(zl,yl,xl), gradz(zl,yl,xl)
double precision :: rho(zl,yl,xl), rho_N(zl,yl,xl)
double precision :: ux(zl,yl,xl), uy(zl,yl,xl), uz(zl,yl,xl)
double precision :: fi(0:18,1:zl,1:yl,1:xl),ffi(0:18,1:zl,1:yl,1:xl)
double precision :: M(0:18,0:18), M_inv(0:18,0:18)
double precision:: mi(0:18), meq(0:18), fk(0:18), si(0:18), force(3)
double precision:: dnxdx, dnxdy, dnxdz, dnydx, dnydy, dnydz, dnzdx, dnzdy, dnzdz
double precision:: tmp, kurv, val, nu, su, sq

!$OMP PARALLEL DO PRIVATE(dnxdx, dnxdy, dnxdz, dnydx, dnydy, dnydz, dnzdx, dnzdy, dnzdz,&
!$OMP xn, yn, zn, tmp, kurv, val, nu, su, sq, mi, meq, fk, si, force), SCHEDULE(GUIDED)
do x = 2, xl-1
	do y = 2, yl-1
		do z = 2, zl-1
            if(walls_all(z,y,x)>0) cycle ! Fluid point
			force = 0.0d0
			if(flag(z,y,x).and.pout(z,y,x)<=0.and.z>=3.and.x>=3)then ! Interface, but not live point  !!! bottom wall 
            !if(flag(z,y,x))then ! Interface, but not live point
				dnxdx = 0.0d0
				dnxdy = 0.0d0
				dnxdz = 0.0d0
				dnydx = 0.0d0
				dnydy = 0.0d0
				dnydz = 0.0d0
				dnzdx = 0.0d0
				dnzdy = 0.0d0
				dnzdz = 0.0d0
				do i = 1, 18
					xn = next_x(i,x)
					yn = next_y(i,y)
					zn = next_z(i,z)               
					tmp = wi(i)*n_x(zn,yn,xn)
					dnxdx = dnxdx + tmp*cix(i)
					dnxdy = dnxdy + tmp*ciy(i)
					dnxdz = dnxdz + tmp*ciz(i)
					tmp = wi(i)*n_y(zn,yn,xn)
					dnydx = dnydx + tmp*cix(i)
					dnydy = dnydy + tmp*ciy(i)
					dnydz = dnydz + tmp*ciz(i)
					tmp = wi(i)*n_z(zn,yn,xn)
					dnzdx = dnzdx + tmp*cix(i)
					dnzdy = dnzdy + tmp*ciy(i)
					dnzdz = dnzdz + tmp*ciz(i)
				end do !i
				dnxdx = 3.0d0*dnxdx
				dnxdy = 3.0d0*dnxdy
				dnxdz = 3.0d0*dnxdz
				dnydx = 3.0d0*dnydx
				dnydy = 3.0d0*dnydy
				dnydz = 3.0d0*dnydz
				dnzdx = 3.0d0*dnzdx
				dnzdy = 3.0d0*dnzdy
				dnzdz = 3.0d0*dnzdz
				!interface curvature
				kurv = n_x(z,y,x)*n_y(z,y,x)*(dnydx+dnxdy)+n_x(z,y,x)*n_z(z,y,x)*(dnzdx+dnxdz)+n_y(z,y,x)*n_z(z,y,x)*(dnydz+dnzdy)&
				 -(n_y(z,y,x)**2+n_z(z,y,x)**2)*dnxdx-(n_x(z,y,x)**2+n_z(z,y,x)**2)*dnydy-(n_x(z,y,x)**2+n_y(z,y,x)**2)*dnzdz
				val = 0.5d0*sigma*kurv
				force(1) =val*gradx(z,y,x)
				force(2) =val*grady(z,y,x)
				force(3) =val*gradz(z,y,x)
			endif
			! update velocity
			ux(z,y,x) = ux(z,y,x)+0.5d0*force(1)/rho(z,y,x)
			uy(z,y,x) = uy(z,y,x)+0.5d0*force(2)/rho(z,y,x)
			uz(z,y,x) = uz(z,y,x)+0.5d0*force(3)/rho(z,y,x)
			! relaxation time
			if (rho_N(z,y,x) <= -1.0d0) then
				nu = nu_b
			elseif (rho_N(z,y,x) >= 1.0d0) then
				nu = nu_r
			else
				nu = 2.0d0/((1.0d0+rho_N(z,y,x))/nu_r+(1.0d0-rho_N(z,y,x))/nu_b)
			endif
			! TRT collision model
			su = 1.0d0/(3.0d0*nu+0.5d0)
			sq = 8.0d0*(2.0d0-su)/(8.0d0-su)
			si = (/1.0d0,su,su,1.0d0,sq,1.0d0,sq,1.0d0,sq,su,su,su,su,su,su,su,sq,sq,sq/) 
			mi = MATMUL(M,fi(:,z,y,x))
			call equilf(meq, fk, rho(z,y,x), ux(z,y,x), uy(z,y,x), uz(z,y,x), force(1), force(2), force(3))
			mi = mi-si*(mi-meq)+(1.0d0-0.5d0*si)*fk
			fi(:,z,y,x) = MATMUL(M_inv,mi)
            ffi(:,z,y,x)=fi(:,z,y,x)
		end do !z
	end do !y
end do !x
!$OMP END PARALLEL DO
return
end subroutine collision

subroutine Recolor(flag,fi,ri,bi,rho_R,rho_B,rho,n_x,n_y,n_z,walls_all)
use const
implicit none
integer x, y, z, i 
logical, dimension(zl,yl,xl) :: flag
double precision, dimension(0:18,1:zl,1:yl,1:xl):: fi, ri, bi
double precision, dimension(zl,yl,xl) :: rho_R, rho_B
double precision, dimension(zl,yl,xl) :: n_x, n_y, n_z
double precision:: rho(zl,yl,xl)
integer :: walls_all(zl,yl,xl)

!$OMP PARALLEL DO SCHEDULE(GUIDED)
do x = 2, xl-1
	do y = 2, yl-1
		do z = 2, zl-1
            if(walls_all(z,y,x)>=1) cycle
			ri(:,z,y,x) = rho_R(z,y,x) / rho(z,y,x) * fi(:,z,y,x)
			if(flag(z,y,x)) then
				do i = 1, 18
					ri(i,z,y,x)=ri(i,z,y,x)+rho_r(z,y,x)*rho_b(z,y,x)/rho(z,y,x)*beta*wi(i)&
					 *(n_x(z,y,x)*cix(i)+n_y(z,y,x)*ciy(i)+n_z(z,y,x)*ciz(i))
				end do !i
			endif
			if (rho_R(z,y,x) <= EVAP_LIM) then
				ri(:,z,y,x) = 0.0d0
			endif
			if (rho_B(z,y,x) <= EVAP_LIM) then
				ri(:,z,y,x) = fi(:,z,y,x)
			endif
			bi(:,z,y,x) = fi(:,z,y,x) - ri(:,z,y,x) 
		end do !z
	end do !y
end do !x
!$OMP END PARALLEL DO 
return
end subroutine Recolor


subroutine streaming(fi, next_x, next_y, next_z)
use const, only: xl, yl, zl, opc
implicit none
integer, intent(in):: next_x(18,xl), next_y(18,yl), next_z(18,zl)
double precision, dimension(0:18,zl,yl,xl), intent(out):: fi
double precision, dimension(0:18,zl,yl,xl):: ft
integer:: x, y, z, i

!$OMP PARALLEL DO SCHEDULE(GUIDED)
do x=1,xl
    do y=1,yl
        do z=1,zl
            do i=1,18
                ft(i,z,y,x)=fi(i, next_z(opc(i),z), next_y(opc(i),y), next_x(opc(i),x))
            end do
        end do
    end do
end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO SCHEDULE(GUIDED)
do x=1,xl
    do y=1,yl
        do z=1,zl
            do i=1,18
                fi(i,z,y,x)=ft(i,z,y,x)
            end do
        end do
    end do
end do
!$OMP END PARALLEL DO

return
end subroutine streaming

subroutine calmacros(fi,ri,bi,rho,ux,uy,uz,rho_R,rho_B,rho_N,walls_all)
use const
implicit none
double precision, dimension(0:18,1:zl,1:yl,1:xl):: fi
double precision, dimension(0:18,1:zl,1:yl,1:xl) :: ri, bi
double precision, dimension(zl,yl,xl) :: rho, ux, uy, uz
double precision, dimension(zl,yl,xl) :: rho_N
double precision, dimension(zl,yl,xl) :: rho_R, rho_B
integer :: walls_all(zl,yl,xl)
integer x, y, z
!$OMP PARALLEL DO SCHEDULE(GUIDED)
do x = 2, xl-1
  do y = 2, yl-1
    do z = 2, zl-1
      if(walls_all(z,y,x)>=1)then
          rho_N(z,y,x)=0.0d0
          ux(z,y,x)=0.0d0
          uy(z,y,x)=0.0d0
          uz(z,y,x)=0.0d0
          rho(z,y,x)=rho0
      else
          ! Only fluid point
        fi(:,z,y,x) = ri(:,z,y,x) + bi(:,z,y,x)
        rho_R(z,y,x) = sum(ri(:,z,y,x))
        rho_B(z,y,x) = sum(bi(:,z,y,x))
        rho(z,y,x) = rho_R(z,y,x)+rho_B(z,y,x)
        ux(z,y,x) = sum(fi(:,z,y,x)*cix)/rho(z,y,x)
        uy(z,y,x) = sum(fi(:,z,y,x)*ciy)/rho(z,y,x)
        uz(z,y,x) = sum(fi(:,z,y,x)*ciz)/rho(z,y,x)
        rho_N(z,y,x) = (rho_R(z,y,x)-rho_B(z,y,x))/rho(z,y,x)
      end if
    end do !z
  end do !y
end do !x
!$OMP END PARALLEL DO
!!!  Walls
rho_N(1,:,:)=rho_N(2,:,:)  !!!   z=1 is wall
rho_N(:,:,1)=rho_N(:,:,2)  !!!   x=1 is wall
!!! symmetry
call symmetry_X_phi(rho_N)
call symmetry_X_phi(rho_R)
call symmetry_X_phi(rho_B)
call symmetry_X_phi(rho)
call symmetry_X_u(ux,uy,uz)
call symmetry_X_fi(fi)

call symmetry_Y_phi(rho_N)
call symmetry_Y_phi(rho_R)
call symmetry_Y_phi(rho_B)
call symmetry_Y_phi(rho)
call symmetry_Y_u(ux,uy,uz)
call symmetry_Y_fi(fi)

call symmetry_Z_phi(rho_N)
call symmetry_Z_phi(rho_R)
call symmetry_Z_phi(rho_B)
call symmetry_Z_phi(rho)
call symmetry_Z_u(ux,uy,uz)
call symmetry_Z_fi(fi)

return
end subroutine calmacros

subroutine equilf_SRT(feq, rho, ux, uy, uz)
use const, only: cix, ciy, ciz, wi
implicit none
double precision:: usqr, uc
double precision, intent(IN):: rho, ux, uy, uz
double precision, intent(OUT):: feq(0:18)
integer i

usqr = ux*ux + uy*uy + uz*uz
do i = 0, 18
	uc = ux*cix(i) + uy*ciy(i) + uz*ciz(i)
	feq(i) = rho*wi(i)*(1.0d0+3.0d0*uc+4.5d0*uc*uc-1.5d0*usqr)
end do !i

return
end subroutine equilf_SRT

subroutine initialization(lamda_d,rho,rho_N,rho_R,rho_B,ux,uy,uz,fi,ffi,ri,bi)
use const
implicit none
integer x, y, z
double precision, dimension(zl,yl,xl) :: rho_r, rho_b, rho_N
double precision, dimension(zl,yl,xl) :: ux, uy, uz, rho
double precision :: fi(0:18,1:zl,1:yl,1:xl),ffi(0:18,1:zl,1:yl,1:xl)
double precision :: ri(0:18,1:zl,1:yl,1:xl), bi(0:18,1:zl,1:yl,1:xl)
double precision:: x0,y0,z0, sqr,dis,dis_x,lamda_d,star,plus
ux = 0.0d0
uy = 0.0d0
uz = 0.0d0
rho = rho0
rho_r = 0.0d0
rho_b = 0.0d0
x0=real(xl)-0.5d0
y0=real(yl)-0.5d0
z0=real(zl)-0.5d0
do x=1,xl
  do y=1,yl
    do z=1,zl
      dis=(dble(x)-x0)**2+(dble(z)-z0)**2
      dis_x=0.1d0*Rad_c*cos(2.0d0*PI*(dble(y)-1.5d0)/lamda_d)
      star=Rad_c+dis_x-sqrt(dis)
      rho_N(z,y,x)=tanh(star*6.0d0*beta*0.134d0) ! D3Q19: k=0.134, D2Q9: k=0.1504
      plus=sigma/(Rad_c+dis_x)*3.0d0
      
      if(rho_N(z,y,x)>0.0d0)then
        rho(z,y,x)=rho0+plus
      else
        rho(z,y,x)=rho0
      endif
      
      rho_R(z,y,x)=rho(z,y,x)*(rho_N(z,y,x)+1.0d0)/2.0d0
      rho_B(z,y,x)=rho(z,y,x)-rho_R(z,y,x)
      
      !rho_R(z,y,x)=rho0*(rho_N(z,y,x)+1.0d0)/2.0d0
      !rho_B(z,y,x)=rho0-rho_R(z,y,x)
      call equilf_SRT(ri(:,z,y,x), rho_R(z,y,x), ux(z,y,x), uy(z,y,x), uz(z,y,x))
	  call equilf_SRT(bi(:,z,y,x), rho_B(z,y,x), ux(z,y,x), uy(z,y,x), uz(z,y,x))
    end do
  end do
end do
fi=ri+bi
ffi=fi
return 
end subroutine initialization

subroutine output_results(t,walls_all,rho_N,ux,uy,uz,rho)
use const
implicit none
integer :: walls_all(zl,yl,xl)
double precision, dimension(zl,yl,xl) ::  rho_N,rho
double precision, dimension(zl,yl,xl) :: ux, uy, uz
integer :: t
integer x, y, z
character(len=50):: filename

write(filename,*) t

open(unit=100, file='OUT_'//trim(adjustl(filename))//'.dat',status='unknown')
write(100,*)'TITLE="Contact_Angles"'
write(100,*)'VARIABLES=X,Y,Z,WALL,RHO_N,U,V,W,RHO'
write(100,*)'ZONE F=POINT I=', xl, 'J=', yl, 'K=', zl

do z = 1, zl
    do y = 1, yl
        do x = 1, xl
            write(100,"(4(1xI6),5(1xF12.8))") x, y, z, walls_all(z,y,x),rho_N(z,y,x), ux(z,y,x), uy(z,y,x), uz(z,y,x),rho(z,y,x)
        end do !x
    end do !y
end do !z
close(100)
return
  end subroutine output_results

  subroutine initia_center(center,center0,lamda_d)
  use const
  implicit none
  real*8::center(par,3),center0(par,3),a,b,x0,y0,z0,dis,dis_x,lamda_d,star,tab
  integer::pp,num_col,num,col(3),row(3)

  col=(/4,8,12/)
  row=(/4,4,4/)
  
    tab=dble(yl-2)/(dble(par)/3.0-0.5)  !!交错
    !tab=dble(yl-2)/(dble(par)/3.0-1.0)  !!全对齐
    
    !col=(/6,12/)
    !row=(/6,6/)
    !tab=dble(yl-2)/(dble(par)/2.0-1.0)  !! 象限点Only 全对齐
  
  x0=real(xl)-0.5d0
  y0=real(yl)-0.5d0
  z0=real(zl)-0.5d0
  
  num_col=1
  num=0
  do pp=1,par
      num=num+1
      if(pp>col(num_col))then
          num_col=num_col+1
      end if
      !write(*,*)'center,num=',num_col
      a=0.0d0-dble(num_col-1)*pi/4.0d0   ! 8方位
      !a=0.0d0-dble(num_col-1)*pi/2.0d0    ! 象限点4方位
      if(mod(num_col,2)==0)then  !偶数行
        b=1.5d0+(real(num)-1)*tab
      else
        b=1.5d0+(real(num)-0.5d0)*tab
      end if
      !b=1.5d0+(real(num)-1.0d0)*tab   !!! 全对齐

      if(num==row(num_col))then
          num=0
      end if
  
      dis_x=0.1d0*Rad_c*cos(2.0d0*PI*(dble(b)-1.5)/lamda_d)
      star=Rad_c+dis_x !扰动后的液柱半径
      star=sqrt(rad**2+star**2+2.0d0*rad*star*cos(cita))
      center(pp,1)=x0-star*cos(a)
      center(pp,2)=b
      center(pp,3)=z0+star*sin(a)
  end do
  center0=center
  return
    end subroutine initia_center
    
 
  subroutine unitization(a)
  implicit none
  real*8:: a(3),star
  star=a(1)**2+a(2)**2+a(3)**2
  star=sqrt(star)
  a(1)=a(1)/star
  a(2)=a(2)/star
  a(3)=a(3)/star
  return
  end subroutine unitization
  
  subroutine curve(walls,pout,ri,bi,upx,upy,upz,rpx,rpy,rpz,center,rho_R,rho_B)   
  use const
  implicit none
  real*8 :: ri(0:18,zl,yl,xl),bi(0:18,zl,yl,xl)
  real*8 :: upx(par),upy(par),upz(par),rpx(par),rpy(par),rpz(par),center(par,3),rho_R(zl,yl,xl),rho_B(zl,yl,xl)
  integer:: walls(zl,yl,xl),pout(zl,yl,xl)
  integer:: x,y,z,i,xb,yb,zb,pp
  real*8 :: xr,yr,zr,star,uw(3)
  !$OMP PARALLEL DO PRIVATE(pp,xb,yb,zb,xr,yr,zr,uw,star), SCHEDULE(GUIDED)
  do x=2,xl-1
    do y=2,yl-1
      do z=2,zl-1
        if(pout(z,y,x)>=1)then ! live fluid point
          pp=pout(z,y,x)
          do i=1,18
            xb=x-cix(i)
            yb=y-ciy(i)
            zb=z-ciz(i)
            if(walls(zb,yb,xb)==pp)then
              xr=dble(x)-center(pp,1)
              yr=dble(y)-center(pp,2)
              zr=dble(z)-center(pp,3)
              uw(1)=upx(pp)+zr*rpy(pp)-yr*rpz(pp)
              uw(2)=upy(pp)+xr*rpz(pp)-zr*rpx(pp)
              uw(3)=upz(pp)+yr*rpx(pp)-xr*rpy(pp)
              !write(*,*)'uw(2)',uw(2)
              star=uw(1)*dble(cix(i))+uw(2)*dble(ciy(i))+uw(3)*dble(ciz(i))
              ri(i,z,y,x)=ri(opc(i),zb,yb,xb)+6.0d0*wi(i)*rho_R(z,y,x)*star
              bi(i,z,y,x)=bi(opc(i),zb,yb,xb)+6.0d0*wi(i)*rho_B(z,y,x)*star
            end if
          end do
        end if
      end do
    end do
  end do
  !$OMP END PARALLEL DO
  
  call symmetry_X_fi(ri)
  call symmetry_X_fi(bi)
  call symmetry_Y_fi(ri)
  call symmetry_Y_fi(bi)
  call symmetry_Z_fi(ri)
  call symmetry_Z_fi(bi)
  return
  end subroutine curve
  
    subroutine curve2(walls,pout,ri,bi,upx,upy,upz,rpx,rpy,rpz,center,rho_R,rho_B)   
  use const
  implicit none
  real*8 :: ri(0:18,zl,yl,xl),bi(0:18,zl,yl,xl)
  real*8 :: upx(par),upy(par),upz(par),rpx(par),rpy(par),rpz(par),center(par,3),rho_R(zl,yl,xl),rho_B(zl,yl,xl)
  integer:: walls(zl,yl,xl),pout(zl,yl,xl)
  integer:: x,y,z,i,xb,yb,zb,xf,yf,zf,pp
  real*8 :: x0,y0,z0,xw,yw,zw,xr,yr,zr,star,star1,star2,qq,dis,uw(3)
  !$OMP PARALLEL DO PRIVATE(pp,xb,yb,zb,xf,yf,zf,x0,y0,z0,star,star1,star2,xw,yw,zw,dis,qq,xr,yr,zr,uw), SCHEDULE(GUIDED)
  do x=2,xl-1
    do y=2,yl-1
      do z=2,zl-1
        if(pout(z,y,x)>=1)then ! live fluid point
          pp=pout(z,y,x)
          do i=1,18
            xb=x-cix(i)
            yb=y-ciy(i)
            zb=z-ciz(i)
            if(walls(zb,yb,xb)==pp)then !后退点是固体点
              xf=x+cix(i)
              yf=y+ciy(i)
              zf=z+ciz(i)
              ! 求qq
              x0=dble(x)-center(pp,1)
              y0=dble(y)-center(pp,2)
              z0=dble(z)-center(pp,3)   ! 坐标平移
              star=cix(i)*x0+ciy(i)*y0+ciz(i)*z0
              star1=cix(i)**2*(Rad**2-y0**2-z0**2)+ciy(i)**2*(Rad**2-x0**2-z0**2) &
                    +ciz(i)**2*(Rad**2-x0**2-y0**2)
              star1=star1+2.0d0*(cix(i)*ciy(i)*x0*y0+ciy(i)*ciz(i)*y0*z0+ciz(i)*cix(i)*z0*x0)
              star1=sqrt(star1)
              star2=cix(i)**2+ciy(i)**2+ciz(i)**2
              xw=x0+cix(i)*(star1-star)/star2     ! +\-star1
              yw=y0+ciy(i)*(star1-star)/star2
              zw=z0+ciz(i)*(star1-star)/star2
              dis=(xw-x0)**2+(yw-y0)**2+(zw-z0)**2
              qq=sqrt(dis)/sqrt(star2)    ! 后退qq
              if(qq>1.0.or.qq<0.0) write(*,*)'warnning curve2 qq', qq
              xr=dble(x)-cix(i)*qq-center(pp,1)
              yr=dble(y)-ciy(i)*qq-center(pp,2)
              zr=dble(z)-ciz(i)*qq-center(pp,3)
              uw(1)=upx(pp)+zr*rpy(pp)-yr*rpz(pp)
              uw(2)=upy(pp)+xr*rpz(pp)-zr*rpx(pp)
              uw(3)=upz(pp)+yr*rpx(pp)-xr*rpy(pp)
              star=uw(1)*dble(cix(i))+uw(2)*dble(ciy(i))+uw(3)*dble(ciz(i))
              ri(i,z,y,x)=1.0d0/(1.0d0+qq)*(qq*(ri(i,zf,yf,xf)&
                          +ri(opc(i),zb,yb,xb))+(1.0d0-qq)*ri(opc(i),z,y,x)&
                          +6.0d0*wi(i)*rho_R(z,y,x)*star)   !迁移后 
              bi(i,z,y,x)=1.0d0/(1.0d0+qq)*(qq*(bi(i,zf,yf,xf)&
                          +bi(opc(i),zb,yb,xb))+(1.0d0-qq)*bi(opc(i),z,y,x)&
                          +6.0d0*wi(i)*rho_B(z,y,x)*star)   !迁移后 
            end if
          end do
        end if
      end do
    end do
  end do
  !$OMP END PARALLEL DO
  call symmetry_X_fi(ri)
  call symmetry_X_fi(bi)
  call symmetry_Y_fi(ri)
  call symmetry_Y_fi(bi)
  call symmetry_Z_fi(ri)
  call symmetry_Z_fi(bi)
  return
  end subroutine curve2
  
      subroutine capillary_force(pout,rho_N,gradx,grady,gradz,n_x,n_y,n_z,flag,center,tension,tension_t) 
    use const
    implicit none
    real*8 :: rho_N(zl,yl,xl),gradx(zl,yl,xl),grady(zl,yl,xl),gradz(zl,yl,xl),n_x(zl,yl,xl),n_y(zl,yl,xl),n_z(zl,yl,xl),center(par,3)
    real*8 :: tension(par,3),tension_t(par,3)
    real*8 :: xr(3),temp(3),temp_t(3),star,drhoN_star,nw(3),nt(3),interface_extend(3)
    integer :: pout(zl,yl,xl),x,y,z,pp
    logical :: flag(zl,yl,xl)
    tension=0.0d0
    tension_t=0.0d0
!$OMP PARALLEL DO PRIVATE(pp,nw,xr,interface_extend,temp,temp_t,nt,drhoN_star,star) REDUCTION(+: tension,tension_t), SCHEDULE(GUIDED)
    do x=2,xl-1
      do y=2,yl-1
        do z=2,zl-1
          if(pout(z,y,x)>=1)then ! live fluid point
            pp=pout(z,y,x)
            if(rho_N(z,y,x)>-0.95d0.and.rho_N(z,y,x)<0.95d0)then ! interface point
                nw(1)=dble(x)-center(pp,1)
                nw(2)=dble(y)-center(pp,2)
                nw(3)=dble(z)-center(pp,3)
                xr=nw
                call unitization(nw)       ! unitization 
                interface_extend(1)=-n_x(z,y,x)
                interface_extend(2)=-n_y(z,y,x)
                interface_extend(3)=-n_z(z,y,x)  ! figure to outside
                call VecCross(interface_extend,nw,temp)
                call VecCross(temp,interface_extend,nt)
                call unitization(nt)     ! unitization        
                drhoN_star=0.0d0   
                if(pout(z-1,y,x)==pp)then
                  drhoN_star=drhoN_star+0.5d0*abs(gradz(z,y,x))
                end if
                if(pout(z+1,y,x)==pp)then
                  drhoN_star=drhoN_star+0.5d0*abs(gradz(z,y,x))
                end if    
                if(pout(z,y-1,x)==pp)then
                  drhoN_star=drhoN_star+0.5d0*abs(grady(z,y,x))
                end if
                if(pout(z,y+1,x)==pp)then
                  drhoN_star=drhoN_star+0.5d0*abs(grady(z,y,x))
                end if
                if(pout(z,y,x-1)==pp)then
                  drhoN_star=drhoN_star+0.5d0*abs(gradx(z,y,x))
                end if
                if(pout(z,y,x+1)==pp)then
                  drhoN_star=drhoN_star+0.5d0*abs(gradx(z,y,x))
                end if
                star=0.75d0*(1.0d0-rho_N(z,y,x)**2)*drhoN_star
                if(cita==PI/2.0)then
                  temp=star*nw*sigma
                else
                  temp=star*nt*sigma
                end if
                tension(pp,:)=tension(pp,:)+temp
                call VecCross(xr,temp,temp_t)
                tension_t(pp,:)=tension_t(pp,:)+temp_t
            end if
          end if
        end do
      end do
    end do
      !$OMP END PARALLEL DO
    return
  end subroutine capillary_force
  
    subroutine VecCross(aa,bb,cc)
    implicit none
    real*8 aa(3),bb(3),cc(3)
    cc(1)=   aa(2)*bb(3)-aa(3)*bb(2)
    cc(2)=-( aa(1)*bb(3)-aa(3)*bb(1) )
    cc(3)=   aa(1)*bb(2)-aa(2)*bb(1)
    return
  end subroutine VecCross  
  
  
  subroutine momentum(walls,pout,center,Fpx,Fpy,Fpz,Tpx,Tpy,Tpz,upx,upy,upz,rpx,rpy,rpz,ffi,fi)
use const
implicit none
integer:: pout(zl,yl,xl),walls(zl,yl,xl)
real*8 :: ffi(0:18,zl,yl,xl),fi(0:18,zl,yl,xl)
real*8 :: center(par,3),Fpx(par),Fpy(par),Fpz(par),Tpx(par),Tpy(par),Tpz(par),upx(par),upy(par),upz(par),rpx(par),rpy(par),rpz(par)
real*8 :: fmomx,fmomy,fmomz,torquex,torquey,torquez,xr,yr,zr,uw(3)
real*8 :: x0,y0,z0,star,star1,star2,xw,yw,zw,dis,qq
integer:: x,y,z,i,xb,yb,zb,pp
Fpx=0.0d0
Fpy=0.0d0
Fpz=0.0d0
Tpx=0.0d0
Tpy=0.0d0
Tpz=0.0d0
!$OMP PARALLEL DO PRIVATE(pp,xb,yb,zb,fmomx,fmomy,fmomz,torquex,torquey,torquez,xr,yr,zr,uw, &
!$OMP x0,y0,z0,star,star1,star2,xw,yw,zw,dis,qq) REDUCTION(+:Fpx,Fpy,Fpz,Tpx,Tpy,Tpz), SCHEDULE(GUIDED)
do x=2,xl-1
  do y=2,yl-1
    do z=2,zl-1
      if(pout(z,y,x)>=1)then ! live fluid point
        pp=pout(z,y,x)
        do i=1,18
          xb=x-cix(i)
          yb=y-ciy(i)
          zb=z-ciz(i)
          fmomx=0.0d0
          fmomy=0.0d0
          fmomz=0.0d0
          torquex=0.0d0
          torquey=0.0d0
          torquez=0.0d0
          if(walls(zb,yb,xb)==pp)then  ! back to a solid point
            if(1==2)then
             ! 求qq
              x0=dble(x)-center(pp,1)
              y0=dble(y)-center(pp,2)
              z0=dble(z)-center(pp,3)   ! 坐标平移
              star=cix(i)*x0+ciy(i)*y0+ciz(i)*z0
              star1=cix(i)**2*(Rad**2-y0**2-z0**2)+ciy(i)**2*(Rad**2-x0**2-z0**2) &
                    +ciz(i)**2*(Rad**2-x0**2-y0**2)
              star1=star1+2.0d0*(cix(i)*ciy(i)*x0*y0+ciy(i)*ciz(i)*y0*z0+ciz(i)*cix(i)*z0*x0)
              star1=sqrt(star1)
              star2=cix(i)**2+ciy(i)**2+ciz(i)**2
              xw=x0+cix(i)*(star1-star)/star2     ! +\-star1
              yw=y0+ciy(i)*(star1-star)/star2
              zw=z0+ciz(i)*(star1-star)/star2
              dis=(xw-x0)**2+(yw-y0)**2+(zw-z0)**2
              qq=sqrt(dis)/sqrt(star2)    ! 后退qq
              if(qq>1.0.or.qq<0.0) write(*,*)'warnning curve2 qq', qq
              xr=dble(x)-cix(i)*qq-center(pp,1)
              yr=dble(y)-ciy(i)*qq-center(pp,2)
              zr=dble(z)-ciz(i)*qq-center(pp,3)
            end if
            
              xr=dble(x)-center(pp,1)
              yr=dble(y)-center(pp,2)
              zr=dble(z)-center(pp,3)
              
              uw(1)=upx(pp)+zr*rpy(pp)-yr*rpz(pp)
              uw(2)=upy(pp)+xr*rpz(pp)-zr*rpx(pp)
              uw(3)=upz(pp)+yr*rpx(pp)-xr*rpy(pp)
              
            fmomx=((dble(cix(i))-uw(1))*fi(i,z,y,x)-(dble(cix(opc(i)))-uw(1))*ffi(opc(i),z,y,x))
            fmomy=((dble(ciy(i))-uw(2))*fi(i,z,y,x)-(dble(ciy(opc(i)))-uw(2))*ffi(opc(i),z,y,x))
            fmomz=((dble(ciz(i))-uw(3))*fi(i,z,y,x)-(dble(ciz(opc(i)))-uw(3))*ffi(opc(i),z,y,x))
            torquex=yr*fmomz-zr*fmomy
            torquey=zr*fmomx-xr*fmomz
            torquez=xr*fmomy-yr*fmomx
            Fpx(pp)=Fpx(pp)-fmomx
            Fpy(pp)=Fpy(pp)-fmomy
            Fpz(pp)=Fpz(pp)-fmomz
            Tpx(pp)=Tpx(pp)-torquex
            Tpy(pp)=Tpy(pp)-torquey
            Tpz(pp)=Tpz(pp)-torquez
          end if
        end do
      end if
    end do
  end do
end do
!$OMP END PARALLEL DO
return
  end subroutine momentum
  
 subroutine movement(center,center_old,center0,Fpx,Fpy,Fpz,Frex,Frey,Frez,Tpx,Tpy,Tpz,upx,upy,upz,rpx,rpy,rpz,Fg,mass,inertia_moment,tension,tension_t)
use const
implicit none
real*8 :: Fpx(par),Fpy(par),Fpz(par),Tpx(par),Tpy(par),Tpz(par),upx(par),upy(par),upz(par),rpx(par),rpy(par),rpz(par),Fg
real*8 :: Frex(par),Frey(par),Frez(par)
real*8 :: tension(par,3),tension_t(par,3)
real*8 :: mass,inertia_moment,center(par,3),center_old(par,3),center0(par,3),temp
integer:: pp
!对称
do pp=1,par
  if(abs(center0(pp,1)-(dble(xl)-0.5d0))<0.1)then   ! 初始在YOZ面上
    !write(*,*)'move',pp
    Fpx(pp)=0.0d0
    Fpy(pp)=2.0d0*Fpy(pp)
    Fpz(pp)=2.0d0*Fpz(pp)
    Tpx(pp)=2.0d0*Tpx(pp)
    Tpy(pp)=0.0d0
    Tpz(pp)=0.0d0
    Frex(pp)=0.0d0
    tension(pp,1)=0.0d0
    tension(pp,2)=tension(pp,2)*2.0d0
    tension(pp,3)=tension(pp,3)*2.0d0
    tension_t(pp,1)=2.0d0*tension_t(pp,1)
    tension_t(pp,2)=0.0d0
    tension_t(pp,3)=0.0d0
  end if
  if(abs(center0(pp,2)-(dble(yl)-0.5d0))<0.1.or.abs(center0(pp,2)-1.5d0)<0.1)then   ! 初始在中间XOZ面上
    Fpx(pp)=2.0d0*Fpx(pp)
    Fpy(pp)=0.0d0
    Fpz(pp)=2.0d0*Fpz(pp)
    Tpx(pp)=0.0d0
    Tpy(pp)=2.0d0*Tpy(pp)
    Tpz(pp)=0.0d0
    Frey(pp)=0.0d0
    tension(pp,1)=2.0d0*tension(pp,1)
    tension(pp,2)=0.0d0
    tension(pp,3)=2.0d0*tension(pp,3)
    tension_t(pp,1)=0.0d0
    tension_t(pp,2)=0.0d0!2.0d0*tension_t(pp,2)
    tension_t(pp,3)=0.0d0
  end if
  if(abs(center0(pp,3)-(dble(zl)-0.5d0))<0.1)then   ! 初始在中间XOY面上
    Fpx(pp)=2.0d0*Fpx(pp)
    Fpy(pp)=2.0d0*Fpy(pp)
    Fpz(pp)=0.0d0
    Tpx(pp)=0.0d0
    Tpy(pp)=0.0d0
    Tpz(pp)=2.0d0*Tpz(pp)
    Frez(pp)=0.0d0
    tension(pp,1)=tension(pp,1)*2.0d0
    tension(pp,2)=tension(pp,2)*2.0d0
    tension(pp,3)=0.0d0
    tension_t(pp,1)=0.0d0
    tension_t(pp,2)=0.0d0
    tension_t(pp,3)=2.0d0*tension_t(pp,3)
  end if
end do

do pp=1,par
 if(pp==par/3+1)then
    temp=(tension(pp,1)+tension(pp,3))/2.0d0
    tension(pp,1)=temp
    tension(pp,3)=temp
    tension_t(pp,2)=0.0d0
    temp=(Fpx(pp)+Fpz(pp))/2.0d0
    Fpx(pp)=temp
    Fpz(pp)=temp
    Tpy(pp)=0.0d0
    temp=(Frex(pp)+Frez(pp))/2.0d0
    Frex(pp)=temp
    Frez(pp)=temp    
  end if
end do

upx(:)=upx(:)+(Fpx(:)+Frex(:)+tension(:,1))/mass
upy(:)=upy(:)+(Fpy(:)+Frey(:)+tension(:,2))/mass
upz(:)=upz(:)+(Fpz(:)+Frez(:)+tension(:,3))/mass
rpx(:)=rpx(:)+(Tpx(:)+tension_t(:,1))/inertia_moment
rpy(:)=rpy(:)+(Tpy(:)+tension_t(:,2))/inertia_moment
rpz(:)=rpz(:)+(Tpz(:)+tension_t(:,3))/inertia_moment
center_old=center
center(:,1)=center(:,1)+upx(:)+0.5d0*(Fpx(:)+Frex(:)+tension(:,1))/mass
center(:,2)=center(:,2)+upy(:)+0.5d0*(Fpy(:)+Frey(:)+tension(:,2))/mass
center(:,3)=center(:,3)+upz(:)+0.5d0*(Fpz(:)+Frez(:)+tension(:,3))/mass

return
  end subroutine movement
  
 subroutine refilling(next_x,next_y,next_z,pout,walls_old,center_old,upx,upy,upz,rpx,rpy,rpz,rho,rho_N,rho_R,rho_B,ux,uy,uz,ri,bi,fi) 
 use const
 implicit none
 integer :: pout(zl,yl,xl),walls_old(zl,yl,xl),next_x(18,xl),next_y(18,yl),next_z(18,zl)
 real*8 :: center_old(par,3),upx(par),upy(par),upz(par),rpx(par),rpy(par),rpz(par)
 real*8 :: rho(zl,yl,xl),rho_N(zl,yl,xl),rho_R(zl,yl,xl),rho_B(zl,yl,xl),ux(zl,yl,xl),uy(zl,yl,xl),uz(zl,yl,xl)
 real*8 :: ri(0:18,zl,yl,xl),bi(0:18,zl,yl,xl),fi(0:18,zl,yl,xl)
 real*8 :: xr,yr,zr,sum_r,sum_b,num
 integer:: x,y,z,i,pp,xn,yn,zn
 !$OMP PARALLEL DO PRIVATE(pp,xr,yr,zr,xn,yn,zn,sum_r,sum_b,num), SCHEDULE(GUIDED)
 do x=2,xl-1
   do y=2,yl-1
     do z=2,zl-1
       if(pout(z,y,x)>=1.and.walls_old(z,y,x)>=1)then
         pp=pout(z,y,x)
         xr=dble(x)-center_old(pp,1)
         yr=dble(y)-center_old(pp,2)
         zr=dble(z)-center_old(pp,3)
         ux(z,y,x)=upx(pp)+zr*rpy(pp)-yr*rpz(pp)
         uy(z,y,x)=upy(pp)+xr*rpz(pp)-zr*rpx(pp)
         uz(z,y,x)=upz(pp)+yr*rpx(pp)-xr*rpy(pp)
         sum_r=0.0d0
         sum_b=0.0d0
         num=0.0d0
         do i=1,18
           zn=next_z(i,z)
           yn=next_y(i,y)
           xn=next_x(i,x)
           if(walls_old(zn,yn,xn)<=0)then
             sum_r=sum_r+rho_r(next_z(i,z),next_y(i,y),next_x(i,x))
             sum_b=sum_b+rho_b(next_z(i,z),next_y(i,y),next_x(i,x))
             num=num+1.0d0
           end if
         end do
         rho_R(z,y,x)=sum_r/num
         rho_B(z,y,x)=sum_b/num
         rho(z,y,x)=rho_R(z,y,x)+rho_B(z,y,x)
         rho_N(z,y,x)=(rho_R(z,y,x)-rho_B(z,y,x))/rho(z,y,x)         
         call equilf_SRT(ri(:,z,y,x), rho_r(z,y,x), ux(z,y,x), uy(z,y,x), uz(z,y,x))
         call equilf_SRT(bi(:,z,y,x), rho_b(z,y,x), ux(z,y,x), uy(z,y,x), uz(z,y,x))
         fi(:,z,y,x)=ri(:,z,y,x)+bi(:,z,y,x)
       end if
     end do
   end do
 end do
 !$OMP END PARALLEL DO
call symmetry_X_phi(rho_N)
call symmetry_X_phi(rho_R)
call symmetry_X_phi(rho_B)
call symmetry_X_phi(rho)
call symmetry_X_u(ux,uy,uz)
call symmetry_X_fi(fi)

call symmetry_Y_phi(rho_N)
call symmetry_Y_phi(rho_R)
call symmetry_Y_phi(rho_B)
call symmetry_Y_phi(rho)
call symmetry_Y_u(ux,uy,uz)
call symmetry_Y_fi(fi)

call symmetry_Z_phi(rho_N)
call symmetry_Z_phi(rho_R)
call symmetry_Z_phi(rho_B)
call symmetry_Z_phi(rho)
call symmetry_Z_u(ux,uy,uz)
call symmetry_Z_fi(fi)
 return
  end subroutine refilling
  
    subroutine output_result2(center,t,tr,upx,upy,upz,rpx,rpy,rpz,Fpx,Fpy,Fpz,Tpx,Tpy,Tpz,tension,tension_t,Frex,Frey,Frez)
  use const
  implicit none
  real*8 :: upx(par),upy(par),upz(par),rpx(par),rpy(par),rpz(par)
  real*8 :: Fpx(par),Fpy(par),Fpz(par),Tpx(par),Tpy(par),Tpz(par)
  real*8 :: tension(par,3),tension_t(par,3)
  real*8 :: Frex(par),Frey(par),Frez(par)
  real*8 :: center(par,3),tr
  character(len=50)::filename1,filename2
  integer:: pp,t
  do pp=1,par
    write(filename1,*)t
    write(filename2,*)pp
    open(unit=pp,file='P_'//trim(adjustl(filename2))//'.dat',status='unknown')
    if(t==0)then
      write(pp,*)'VARIABLES=t,t*,cx,cy,cz,upx,upy,upz,rpx,rpy,rpz,Fpx,Fpy,Fpz,Tpx,Tpy,Tpz,Fcx,Fcy,Fcz,Tcx,Tcy,Tcz,Frex,Frey,Frez'
    end if
    write(pp,"(1(1xI6),25(1xF12.6))")t,t/tr,center(pp,1),center(pp,2),center(pp,3),upx(pp),upy(pp),upz(pp),rpx(pp),rpy(pp),rpz(pp),Fpx(pp),Fpy(pp),Fpz(pp),Tpx(pp),Tpy(pp),Tpz(pp),&
                                        tension(pp,1),tension(pp,2),tension(pp,3),tension_t(pp,1),tension_t(pp,2),tension_t(pp,3),Frex(pp),Frey(pp),Frez(pp)
  end do
  return
  end subroutine output_result2

subroutine symmetry_X_fi(fi)
use const
implicit none
real*8 :: fi(0:18,zl,yl,xl)
integer:: x, y, z, i
!symmetry about Y0Z
do i=0,18
  fi(i,:,:,xl)=fi(symx(i),:,:,xl-1)
end do
return
end subroutine symmetry_X_fi 
  
subroutine symmetry_X_phi(phi)
use const
implicit none
real*8 :: phi(zl,yl,xl)
!对称边界 about Y0Z
phi(:,:,xl)=phi(:,:,xl-1)
return
  end subroutine symmetry_X_phi
  
subroutine symmetry_X_u(ux,uy,uz)
use const
implicit none
real*8 :: ux(zl,yl,xl),uy(zl,yl,xl),uz(zl,yl,xl)
!对称边界 about Y0Z
ux(:,:,xl)=-ux(:,:,xl-1)
uy(:,:,xl)= uy(:,:,xl-1)
uz(:,:,xl)= uz(:,:,xl-1)
return
  end subroutine  symmetry_X_u
  
subroutine symmetry_Y_fi(fi)
use const
implicit none
real*8 :: fi(0:18,zl,yl,xl)
integer:: x, y, z, i
!symmetry about Y0Z
  do i=0,18
    fi(i,:,yl,:)=fi(symy(i),:,yl-1,:)
    fi(i,:,1,:) =fi(symy(i),:,2,:)
  end do
return
  end subroutine symmetry_Y_fi
  
  subroutine symmetry_Y_phi(phi)
use const
implicit none
real*8 :: phi(zl,yl,xl)
!对称边界 about X0Z
phi(:,yl,:)=phi(:,yl-1,:)
phi(:,1,:) =phi(:,2,:)
return
  end subroutine symmetry_Y_phi
  
  subroutine symmetry_Y_u(ux,uy,uz)
use const
implicit none
real*8 :: ux(zl,yl,xl),uy(zl,yl,xl),uz(zl,yl,xl)
!对称边界 about Y0Z
ux(:,yl,:)= ux(:,yl-1,:)
uy(:,yl,:)=-uy(:,yl-1,:)
uz(:,yl,:)= uz(:,yl-1,:)
ux(:,1,:)= ux(:,2,:)
uy(:,1,:)=-uy(:,2,:)
uz(:,1,:)= uz(:,2,:)
return
  end subroutine  symmetry_Y_u
  
subroutine symmetry_Z_fi(fi)
use const
implicit none
real*8 :: fi(0:18,zl,yl,xl)
integer:: x, y, z, i
!symmetry about XOY
  do i=0,18
    fi(i,zl,:,:)=fi(symz(i),zl-1,:,:)
  end do
return
  end subroutine symmetry_Z_fi
  
subroutine symmetry_Z_phi(phi)
use const
implicit none
real*8 :: phi(zl,yl,xl)
!对称边界 about X0Z
phi(zl,:,:)=phi(zl-1,:,:)
return
end subroutine symmetry_Z_phi
  
  subroutine symmetry_Z_u(ux,uy,uz)
use const
implicit none
real*8 :: ux(zl,yl,xl),uy(zl,yl,xl),uz(zl,yl,xl)
!对称边界 about Y0Z
ux(zl,:,:)= ux(zl-1,:,:)
uy(zl,:,:)= uy(zl-1,:,:)
uz(zl,:,:)=-uz(zl-1,:,:)
return
  end subroutine  symmetry_Z_u  
  
subroutine boundary(next_x,next_y,next_z,ri,bi)
use const
implicit none
integer:: next_x(18,xl),next_y(18,yl),next_z(18,zl)
real*8 :: ri(0:18,zl,yl,xl),bi(0:18,zl,yl,xl)
integer x,y,z,i,xn,yn,zn

call symmetry_X_fi(ri)
call symmetry_X_fi(bi)
call symmetry_Y_fi(ri)
call symmetry_Y_fi(bi)
call symmetry_Z_fi(ri)
call symmetry_Z_fi(bi)
! Bottom wall
z=2
!$OMP PARALLEL DO PRIVATE(xn,yn,zn), SCHEDULE(GUIDED)
do x=1,xl
  do y=1, yl
    do i = 1, 18
      xn=next_x(i,x)
      yn=next_y(i,y)
      zn=next_z(i,z)
      if(zn==1)then
        ri(opc(i),zn,yn,xn)=ri(i,z,y,x)
        bi(opc(i),zn,yn,xn)=bi(i,z,y,x)
      end if
    end do
  end do
end do
!$OMP END PARALLEL DO
! Back wall
x=2
!$OMP PARALLEL DO PRIVATE(xn,yn,zn), SCHEDULE(GUIDED)
do y=1,yl
  do z=1,zl
    do i = 1, 18
      xn=next_x(i,x)
      yn=next_y(i,y)
      zn=next_z(i,z)
      if(xn==1)then
        ri(opc(i),zn,yn,xn)=ri(i,z,y,x)
        bi(opc(i),zn,yn,xn)=bi(i,z,y,x)
      end if
    end do
  end do
end do
!$OMP END PARALLEL DO
return
end subroutine boundary

    subroutine repulsive(center,Frex,Frey,Frez,Fg)
  use const
  implicit none
  real*8 :: center(par,3),Frex(par),Frey(par),Frez(par),Fg
  real*8 :: dis_left,dis_right,dis_top,dis_bottom,dis_front,dis_back,dis(par,par),range,order
  integer :: i,j,pp
  range=3.0d0
  order=0.01d0
  Frex=0.0d0
  Frey=0.0d0
  Frez=0.0d0
  !$OMP PARALLEL DO PRIVATE(dis_top,dis_bottom,dis_left,dis_right,dis_front,dis_back) REDUCTION(+:Frez,Frey,Frex), SCHEDULE(GUIDED)
  do pp=1,par
    dis_top   =dble(zl)-0.5d0-center(pp,3)-rad
    dis_bottom=center(pp,3)-1.5d0-rad   
    dis_left  =center(pp,2)-1.5d0-rad
    dis_right =dble(yl)-1.5d0-center(pp,2)-rad
    dis_front =center(pp,1)-1.5d0-rad
    dis_back  =dble(xl)-1.5d0-center(pp,1)-rad
    
    if(dis_top<2.0)then
        !write(*,*)'warnning dis_top<2.0',pp,center(pp,3)
    end if
    if(dis_bottom<2.0)then
        write(*,*)'warnning dis_bottom<2.0',pp,center(pp,3)
    end if
    if(dis_left<2.0)then
        !write(*,*)'warnning dis_left<2.0',pp,center(pp,2)
    end if
    if(dis_right<2.0)then
        !write(*,*)'warnning dis_right<2.0',pp,center(pp,2)
    end if
    if(dis_front<2.0)then
        write(*,*)'warnning dis_front<2.0',pp,center(pp,1)
    end if
    if(dis_back<2.0)then
        !write(*,*)'warnning dis_back<2.0',pp,center(pp,1)
    end if
    
    if(dis_top<range)then
      Frez(pp)=Frez(pp)+2.0d0*(Fg/order)*(dis_top-range)**2/range**2   !"-"
    end if
    if(dis_bottom<range)then
      Frez(pp)=Frez(pp)-2.0d0*(Fg/order)*(dis_bottom-range)**2/range**2 !"+"
    end if
    if(dis_left<range)then
      Frey(pp)=Frey(pp)-2.0d0*(Fg/order)*(dis_left-range)**2/range**2 !"+"
    end if
    if(dis_right<range)then
      Frey(pp)=Frey(pp)+2.0d0*(Fg/order)*(dis_right-range)**2/range**2 !"-"
    end if
    if(dis_front<range)then
      Frex(pp)=Frex(pp)-2.0d0*(Fg/order)*(dis_front-range)**2/range**2 !"+"
    end if
    if(dis_back<range)then
      Frex(pp)=Frex(pp)+2.0d0*(Fg/order)*(dis_back-range)**2/range**2 ! "-"
    end if
  end do
    !$OMP END PARALLEL DO
    
  !$OMP PARALLEL DO REDUCTION(+:Frez,Frey,Frex),SCHEDULE(GUIDED)
  do i=1,par
      do j=1,par
          if(i==j)cycle
          dis(i,j)=dsqrt((center(i,1)-center(j,1))**2+(center(i,2)-center(j,2))**2+(center(i,3)-center(j,3))**2)-rad*2.0d0
          if(dis(i,j)<2.0d0)then
              write(*,*)'warnning dis_ij<2.0',i,j,dis(i,j)
          end if
          if(dis(i,j)<range)then
              Frex(i)=Frex(i)-(Fg/order)*(dis(i,j)-range)**2/range**2*(center(i,1)-center(j,1))/dis(i,j)
              Frey(i)=Frey(i)-(Fg/order)*(dis(i,j)-range)**2/range**2*(center(i,2)-center(j,2))/dis(i,j)
              Frez(i)=Frez(i)-(Fg/order)*(dis(i,j)-range)**2/range**2*(center(i,3)-center(j,3))/dis(i,j)
          end if
      end do
  end do
  !$OMP END PARALLEL DO
  return 
    end subroutine repulsive

    subroutine output3(rho_N,t,tr,walls_all,pout,lamda_d)
    use const
    implicit none
    real*8 :: rho_N(zl,yl,xl)
    integer :: i,x,y,z,zn,yn,xn,t,walls_all(zl,yl,xl),y0,pout(zl,yl,xl)
    real*8 :: top(yl),bottom(yl),h(yl),hmin,tr,area,energy
    real*8 :: area_R,area_B,lamda_d

    x=xl-1
    hmin=real(zl)
    bottom=1.0
    top=real(zl)-0.5
!$OMP PARALLEL DO PRIVATE(zn),SCHEDULE(GUIDED)
    do y=2,yl-1
        do z=2,zl-1
            if(pout(z,y,x)==0.and.walls_all(z,y,x)==0)then
                zn=z+1
                if(pout(zn,y,x)==0.and.walls_all(zn,y,x)==0)then
                    if(rho_N(z,y,x)*rho_N(zn,y,x)<0.0d0)then  ! 跨界面
                        !write(*,*)'rho_N',rho_N(z,y,x),rho_N(zn,y,x),walls_all(z,y,x),walls_all(z,y,x)
                        bottom(y)=real(zn)-rho_N(zn,y,x)/(rho_N(zn,y,x)-rho_N(z,y,x))
                    end if
                end if
            end if
        end do
        !write(*,*)'bottom,top',bottom(y),top(y)
    end do
  !$OMP END PARALLEL DO
      
    do y=2,yl-1
        h(y)=top(y)-bottom(y)
        if(h(y)<hmin)then
            hmin=h(y)
            y0=y
        end if
    end do
    
    open(unit=300,file='hmin.dat',status='unknown')
    if(t==0)then
        write(300,*)'VARIABLES=t,t*,h,y0'
    end if
    write(300,"(1(1xI8),3(1xF12.6))")t,t/tr,hmin/Rad_c,y0/lamda_d
    
  return
  end subroutine output3
  
      subroutine output4(rho_N,t,tr,walls_all,pout,center0,ux,uy,uz,rho,upx,upy,upz,rpx,rpy,rpz,tension,tension_t,mass,inertia_moment,Wcap_old,Shb)
    use const
    implicit none
    real*8 :: rho_N(zl,yl,xl),ux(zl,yl,xl),uy(zl,yl,xl),uz(zl,yl,xl),rho(zl,yl,xl)
    integer :: i,x,y,z,zn,yn,xn,t,walls_all(zl,yl,xl),pout(zl,yl,xl),pp
    real*8 :: upx(par),upy(par),upz(par),rpx(par),rpy(par),rpz(par),tension(par,3),tension_t(par,3),mass,inertia_moment
    real*8 :: center0(par,3),tr,Er,Epv,Esur,Emov,Wcap(par),Wcap_old,Emov_p(par),Wcap_sum,Emov_psum,area,Shb(par)
    
    Er=pi*Rad_c**2*sigma
    Epv=0.0d0
    Esur=0.0d0
    Emov=0.0d0
    area=0.0d0
      !$OMP PARALLEL DO PRIVATE(xn,yn,zn) REDUCTION(+:Epv,Emov,area),SCHEDULE(GUIDED)
    do x=2,xl-1
      do y=2,yl-1
        do z=2,zl-1
          if(walls_all(z,y,z)<=0)then
            if(rho_N(z,y,x)>0.0d0)then
              Epv=Epv+(rho(z,y,x)-rho0)/3.0d0                                   ! 液柱内压力能
              Emov=Emov+0.5d0*rho(z,y,z)*(ux(z,y,x)**2+uy(z,y,x)**2+uz(z,y,z)**2)
            end if
            do i=1,18
              xn=x+cix(i)
              yn=y+ciy(i)
              zn=z+ciz(i)
              if(rho_N(z,y,x)*rho_N(zn,yn,xn)<0.0.and.walls_all(zn,yn,xn)<=0)then  ! 跨界面
                area=area+1.0
                exit
              end if
            end do
          end if
        end do
      end do
    end do
  !$OMP END PARALLEL DO
Epv=Epv/Er*8.0d0
Emov=Emov/Er*8.0d0
Esur=area/2.0d0*sigma/Er*8.0d0

Emov_p=0.0d0
Wcap=0.0d0
Wcap_sum=0.0d0
Emov_psum=0.0d0
!$OMP PARALLEL DO  REDUCTION(+:Emov_p,Wcap),SCHEDULE(GUIDED)
do pp=1,par
  Emov_p(pp)=0.5*mass*(upx(pp)**2+upy(pp)**2+upz(pp)**2) &  
            +0.5*inertia_moment*(rpx(pp)**2+rpy(pp)**2+rpz(pp)**2)
  Wcap(pp)=tension(pp,1)*upx(pp)+tension(pp,2)*upy(pp)+tension(pp,3)*upz(pp) &
          +tension_t(pp,1)*rpx(pp)+tension_t(pp,2)*rpy(pp)+tension_t(pp,3)*rpz(pp)
  if(abs(center0(pp,2)-(dble(yl)-0.5d0))<0.1.or.abs(center0(pp,2)-1.5d0)<0.1)then   ! 初始在中间XOZ面上
    Emov_p(pp)=Emov_p(pp)*0.5d0
    Wcap(pp)=Wcap(pp)*0.5d0
    Shb(pp)=Shb(pp)*0.5d0
  end if
  if(abs(center0(pp,1)-(dble(xl)-0.5d0))<0.1)then
    Emov_p(pp)=Emov_p(pp)*0.5d0
    Wcap(pp)=Wcap(pp)*0.5d0
    Shb(pp)=Shb(pp)*0.5d0
  end if
  if(abs(center0(pp,3)-(dble(zl)-0.5d0))<0.1)then
    Emov_p(pp)=Emov_p(pp)*0.5d0
    Wcap(pp)=Wcap(pp)*0.5d0
    Shb(pp)=Shb(pp)*0.5d0
  end if
end do
!$OMP END PARALLEL DO
Emov_psum=sum( Emov_p)/Er*8.0d0
Wcap_sum=Wcap_old+sum(Wcap)/Er*8.0d0
Wcap_old=Wcap_sum
Shb=Shb/(pi*Rad_c*Rad_c)*8.0d0

    open(unit=400,file='Energy.dat',status='unknown')
    if(t==0)then
        write(400,*)'VARIABLES=t,t*,Esur,Shb,Efree,Emov,Wcap,Emov_p,Esum'
    end if
    write(400,"(1(1xI8),8(1xF12.6))")t,t/tr,Esur,sum(Shb),Esur-Epv,Emov,Wcap_sum,Emov_psum,(Esur-Epv+Emov+Emov_psum)
  return
    end subroutine output4
  
  subroutine Lagrange_points(walls,center,rho_N,gradx,grady,gradz,n_x,n_y,n_z,flag,Lax,Lay,Laz,La_rhoN,La_nx,La_ny,La_nz,La_flag,ds,Shb)
  use const
  implicit none
  integer :: walls(zl,yl,xl)
  double precision,dimension(zl,yl,xl) :: rho_N,gradx,grady,gradz,n_x,n_y,n_z
  logical,dimension(zl,yl,xl):: flag
  double precision,dimension(par,La_num):: Lax,Lay,Laz,La_rhoN,La_nx,La_ny,La_nz,ds
  logical,dimension(par,La_num):: La_flag
  real*8 :: center(par,3),aa,bb,dis,weight,weight_sum,star,nw(3)
  real*8 :: n_temp(3),theta,a(2),b(2),vec1(3),vec2(3),db1,db2,min,Shb(par)
  integer:: x,y,z,pp,num,x0,y0,z0,xf,yf,zf,ex(8),ey(8),ez(8),i,xx,yy,zz,na,nb
  ex=(/0,1,1,0,0,1,1,0/)
  ey=(/0,0,1,1,0,0,1,1/)
  ez=(/0,0,0,0,1,1,1,1/)
  La_flag=.false.
  La_rhoN=0.0d0
  La_nx=0.0d0
  La_ny=0.0d0
  La_nz=0.0d0
  Shb=0.0d0
  !$OMP PARALLEL DO PRIVATE(num,aa,bb,x0,y0,z0,weight_sum,star,xf,yf,zf,dis,weight,nw,min,xx,yy,zz,n_temp,&
  !$OMP                     theta,a,b,vec1,vec2,db1,db2),SCHEDULE(GUIDED)
  do pp=1,par
    num=1
    do na=1,La_m        !aa=0,2.0*pi-alpha,alpha ! La_m
      do nb=1,La_n+1    !bb=0,pi,alpha         ! La_n+1    
        aa=dble(na-1)*alpha
        bb=dble(nb-1)*alpha
        Lax(pp,num)=center(pp,1)+Rad*sin(bb)*cos(aa)
        Lay(pp,num)=center(pp,2)+Rad*sin(bb)*sin(aa)
        Laz(pp,num)=center(pp,3)+Rad*cos(bb)
        ds(pp,num)=Rad**2*sin(bb)*alpha**2
        !!! rhoN of Lagrange points
        if(Lax(pp,num)<xl-0.5d0.and.Lay(pp,num)>=1.5d0.and.Lay(pp,num)<yl-0.5d0.and.Laz(pp,num)<zl-0.5d0)then !!! within the domain
          x0=int(Lax(pp,num))
          y0=int(Lay(pp,num))
          z0=int(Laz(pp,num))
          weight_sum=0.0d0
          star=0.0d0
          La_rhoN(pp,num)=0.0d0
          do i=1,8        ! cubic
            xf=x0+ex(i)
            yf=y0+ey(i)
            zf=z0+ez(i)
            if(walls(zf,yf,xf)<=0)then ! Fluid
              star=rho_N(zf,yf,xf)+gradx(zf,yf,xf)*(Lax(pp,num)-dble(xf))+grady(zf,yf,xf)*(Lay(pp,num)-dble(yf)) &
                                  +gradz(zf,yf,xf)*(Laz(pp,num)-dble(zf))
              dis=(dble(xf)-Lax(pp,num))**2+(dble(yf)-Lay(pp,num))**2+(dble(zf)-Laz(pp,num))**2
              weight=1.0d0/(dis)
              La_rhoN(pp,num)=La_rhoN(pp,num)+weight*star
              weight_sum=weight_sum+weight
            end if
          end do
          La_rhoN(pp,num)=La_rhoN(pp,num)/weight_sum
          if(La_rhoN(pp,num)<-1.0)then
            !write(*,*)'warnning La_rhoN',La_rhoN(pp,num)
            La_rhoN(pp,num)=-1.0d0
          end if
          if(La_rhoN(pp,num)>1.0)then
            !write(*,*)'warnning La_rhoN',La_rhoN(pp,num)
            La_rhoN(pp,num)=1.0d0
          end if
          if(abs(La_rhoN(pp,num))<0.99d0) La_flag(pp,num)=.true.
        end if

        if(La_rhoN(pp,num)>0.0d0)then
          Shb(pp)=Shb(pp)+ds(pp,num)      ! Aera immersed in the liquid column
        end if
        
        !!! nx,ny,nz of Lagrange points,related to the contact angle
        if(La_flag(pp,num))then
          nw(1)=Lax(pp,num)-center(pp,1)
          nw(2)=Lay(pp,num)-center(pp,2)
          nw(3)=Laz(pp,num)-center(pp,3)
          call unitization(nw)

          min=100.0d0
          xx=x0
          yy=y0
          zz=z0
          do i=1,8        ! cubic
            xf=x0+ex(i)
            yf=y0+ey(i)
            zf=z0+ez(i)
            if(xf<=xl.and.yf>=1.and.yf<=yl.and.zf<=zl)then   ! within the domain
              if(walls(zf,yf,xf)<=0)then                 ! fluid point
                !if(abs(rho_N(zf,yf,xf))<1.0d0)then     ! within the interface
                if(abs(n_x(zf,yf,xf)**2+n_y(zf,yf,xf)**2+n_z(zf,yf,xf)**2)>0.0d0)then! within the interface
                  dis=(dble(xf)-Lax(pp,num))**2+(dble(yf)-Lay(pp,num))**2+(dble(zf)-Laz(pp,num))**2
                  if(sqrt(dis)<min)then
                    min=sqrt(dis)
                    xx=xf
                    yy=yf
                    zz=zf
                  end if
                end if
              end if
            end if
          end do
          La_nx(pp,num)=n_x(zz,yy,xx)
          La_ny(pp,num)=n_y(zz,yy,xx)
          La_nz(pp,num)=n_z(zz,yy,xx)
          
          if(min==100.0d0)then
              write(*,*)'warnning no point',pp,x0,y0,z0,La_rhoN(pp,num)
              write(*,*) rho_N(z0,y0,x0),walls(z0,y0,x0)
              write(*,*) rho_N(z0,y0,x0+1),walls(z0,y0,x0+1)
              write(*,*) rho_N(z0,y0+1,x0+1),walls(z0,y0+1,x0+1)
              write(*,*) rho_N(z0,y0+1,x0+1),walls(z0,y0+1,x0+1)
              write(*,*) rho_N(z0+1,y0,x0),walls(z0+1,y0,x0)
              write(*,*) rho_N(z0+1,y0,x0+1),walls(z0+1,y0,x0+1)
              write(*,*) rho_N(z0+1,y0+1,x0),walls(z0+1,y0+1,x0)
              write(*,*) rho_N(z0+1,y0+1,x0+1),walls(z0+1,y0+1,x0+1)
              !La_nx(pp,num)=0.5d0*nw(2)*nw(3)   !给一个与壁面法矢量垂直的向量
              !La_ny(pp,num)=-nw(1)*nw(3)
              !La_nz(pp,num)=0.5d0*nw(1)*nw(2)
          endif
           
          n_temp=(/La_nx(pp,num), La_ny(pp,num), La_nz(pp,num)/)
          call unitization(n_temp)
          La_nx(pp,num)=n_temp(1)
          La_ny(pp,num)=n_temp(2)
          La_nz(pp,num)=n_temp(3)

          if (  maxval ( abs( nw(:)-n_temp(:) ) ) <1.0d-5 ) then
            theta=1.0d-5                                               !!!! ???
          else
            theta=dacos(dot_product(nw, n_temp))
          endif

          theta=dacos(dot_product(nw, n_temp))

          a(1)=dcos(cita)-dsin(cita)*dcos(theta)/dsin(theta)
          a(2)=dcos(cita)+dsin(cita)*dcos(theta)/dsin(theta)
          b(1)=dsin(cita)/dsin(theta)
          b(2)=-dsin(cita)/dsin(theta)

          vec1=a(1)*nw+b(1)*(/n_temp(1), n_temp(2), n_temp(3)/)
          vec2=a(2)*nw+b(2)*(/n_temp(1), n_temp(2), n_temp(3)/)           
          db1=(vec1(1)-n_temp(1))**2+(vec1(2)-n_temp(2))**2+(vec1(3)-n_temp(3))**2
          db2=(vec2(1)-n_temp(1))**2+(vec2(2)-n_temp(2))**2+(vec2(3)-n_temp(3))**2
          if(db1<db2) then
            La_nx(pp,num)=vec1(1)
            La_ny(pp,num)=vec1(2)
            La_nz(pp,num)=vec1(3)
          else
            La_nx(pp,num)=vec2(1)
            La_ny(pp,num)=vec2(2)
            La_nz(pp,num)=vec2(3)
          end if
          if(abs(La_nx(pp,num))>1.0)then
            write(*,*)'warnning,La_nx',La_nx(pp,num),n_temp(1)*nw(1)+n_temp(2)*nw(2)+n_temp(3)*nw(3)
          end if
          if(abs(La_ny(pp,num))>1.0)then
            write(*,*)'warnning,La_ny',La_ny(pp,num),n_temp(1)*nw(1)+n_temp(2)*nw(2)+n_temp(3)*nw(3)
          end if
          if(abs(La_nz(pp,num))>1.0)then
            write(*,*)'warnning,La_nz',La_nz(pp,num),n_temp(1)*nw(1)+n_temp(2)*nw(2)+n_temp(3)*nw(3)
          end if
        end if

        num=num+1
      end do
    end do
    !write(*,*)'check,La_num',La_num,num-1    
  end do
  !$OMP END PARALLEL DO
  
  !if(1==2)then
    open(unit=500,file='Laglange.dat',status='unknown')
    write(500,*)'VARIABLES=x,y,z,rhoN,nx,ny,nz'
    pp=6
    do num=1,La_num
      write(500,"(7(1xF12.6))")Lax(pp,num),Lay(pp,num),Laz(pp,num),La_rhoN(pp,num),La_nx(pp,num),La_ny(pp,num),La_nz(pp,num)
    end do
    close(500)
  !end if
  
  return
  end subroutine Lagrange_points
  
  subroutine capillary_force2(center,Lax,Lay,Laz,La_nx,La_ny,La_nz,La_rhoN,La_flag,ds,tension,tension_t)
  use const
  implicit none  
  real*8 :: center(par,3),tension(par,3),tension_t(par,3)
  double precision,dimension(par,La_num):: Lax,Lay,Laz,La_rhoN,La_nx,La_ny,La_nz,ds
  logical,dimension(par,La_num):: La_flag
  real*8 :: nw(3),xr(3),interface_extend(3),temp(3),temp_t(3),nt(3),projector,star,k
  integer:: pp,num
  
  k=0.134d0*1.43d0   !!! for D3Q19 color gradient model
  tension=0.0d0
  tension_t=0.0d0
!$OMP PARALLEL DO PRIVATE(nw,xr,interface_extend,temp,projector,nt,star,temp_t) REDUCTION(+:tension,tension_t),SCHEDULE(GUIDED)
do pp=1,par
  do num=1,La_num
    if(La_flag(pp,num))then
      nw(1)=Lax(pp,num)-center(pp,1)
      nw(2)=Lay(pp,num)-center(pp,2)
      nw(3)=Laz(pp,num)-center(pp,3)
      xr=nw
      call unitization(nw)
      interface_extend(1)=-La_nx(pp,num)
      interface_extend(2)=-La_ny(pp,num)
      interface_extend(3)=-La_nz(pp,num)  ! figure to bule
      call VecCross(interface_extend,nw,temp)
      projector=temp(1)**2+temp(2)**2+temp(3)**2
      projector=sqrt(projector)
      
      call VecCross(temp,interface_extend,nt)
      call unitization(nt)     ! direction of the capillary force
      star=4.5d0*k*beta*(1.0d0-La_rhoN(pp,num)**2)**2*ds(pp,num)*projector    !!! capillary force model by Haoran Liu et al. in 2020
      if(isNAN(star))then
        write(*,*)'star_cap2',La_rhoN(pp,num),projector
      end if
      
      if(cita==PI/2.0)then
        temp=star*nw*sigma
      else
        temp=star*nt*sigma
      end if
      tension(pp,:)=tension(pp,:)+temp
      call VecCross(xr,temp,temp_t)
      tension_t(pp,:)=tension_t(pp,:)+temp_t
    end if
  end do
end do
  !$OMP END PARALLEL DO
  return
  end subroutine capillary_force2