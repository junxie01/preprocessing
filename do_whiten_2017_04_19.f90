subroutine do_whiten(ns,dom,sachead,fre,sig,seisout)
use sacio
implicit none
integer FFTW_FORWARD,FFTW_BACKWARD
parameter (FFTW_FORWARD=-1,FFTW_BACKWARD=1)
integer FFTW_ESTIMATE,FFTW_MEASURE
parameter (FFTW_ESTIMATE=0,FFTW_MEASURE=1)
type(sac_head) :: sachead
integer,parameter:: nmax=4000000
integer nerr,i,nk,npts,npow,ns,nn
real*4 fre(4),dt,amax
real*8 dom,f,plan
real*4 sig(nmax,3)
real*4 ss1(nmax),ss2(nmax),ss3(nmax)
complex seisout(nmax,3)
complex s1(nmax),s2(nmax),s3(nmax)
complex sa(nmax),sb(nmax),sc(nmax)
complex saf(nmax),sbf(nmax),scf(nmax)
!
sa=cmplx(0,0)
sb=cmplx(0,0)
sc=cmplx(0,0)
s1=cmplx(0,0)
s2=cmplx(0,0)
s3=cmplx(0,0)
dt=sachead%delta
npts=sachead%npts
seisout=0;ns=2;npow=1
do while(ns.lt.npts)
   ns=ns*2
   npow=npow+1
enddo
nk=ns/2+1
dom=dble(1.0/dt/ns)
!write(*,*)'ns=',ns,'dt=',dt
sa(1:npts)=cmplx(sig(1:npts,1),0)
sb(1:npts)=cmplx(sig(1:npts,2),0)
sc(1:npts)=cmplx(sig(1:npts,3),0)
!write(*,*)"Do fft"
!call clogc(npow,sa,1,dt)
!call clogc(npow,sb,1,dt)
!call clogc(npow,sc,1,dt)
call dfftw_plan_dft_1d(plan,ns,sa,saf,FFTW_FORWARD, FFTW_ESTIMATE)
call dfftw_execute(plan)
call dfftw_destroy_plan(plan)

call dfftw_plan_dft_1d(plan,ns,sb,sbf,FFTW_FORWARD, FFTW_ESTIMATE)
call dfftw_execute(plan)
call dfftw_destroy_plan(plan)

call dfftw_plan_dft_1d(plan,ns,sc,scf,FFTW_FORWARD, FFTW_ESTIMATE)
call dfftw_execute(plan)
call dfftw_destroy_plan(plan)
!sa(nk+1:ns)=cmplx(0,0)
!sb(nk+1:ns)=cmplx(0,0)
!sc(nk+1:ns)=cmplx(0,0)
! smooth the spectrum
!write(*,*)"Smooth the spectrum"
nn=int(0.005d0/dom)
call smoothf(dom,nk,nn,fre,saf,ss1)
call smoothf(dom,nk,nn,fre,sbf,ss2)
call smoothf(dom,nk,nn,fre,scf,ss3)
!write(*,*)"Smooth the spectrum done!",dom,nk,npow
do i=1,nk 
   amax=(ss1(i)+ss2(i)+ss3(i))/3.0
   amax=max(ss1(i),ss2(i),ss3(i))
   amax=ss1(i)
   f = dble((i-1))*dom
   if( f .ge. dble(fre(1)) .and. f .le. dble(fre(4)) ) then
      s1(i)=sa(i)/amax
      s2(i)=sb(i)/amax
      s3(i)=sc(i)/amax
   endif
enddo
! taper the spectrum
!write(*,*)"Taper the spectrum!"
!,dom,nk,npow
call  taperf(s1,fre,dom,nk,npow)
call  taperf(s2,fre,dom,nk,npow)
call  taperf(s3,fre,dom,nk,npow)
!s1(nk+1:ns)=cmplx(0,0)
!s2(nk+1:ns)=cmplx(0,0)
!s3(nk+1:ns)=cmplx(0,0)
!write(*,*)"Do ifft"
!call  clogc(npow,s1,-1,dt)
!call  clogc(npow,s2,-1,dt)
!call  clogc(npow,s3,-1,dt)
seisout(1:ns,1) = s1(1:ns)
seisout(1:ns,2) = s2(1:ns)
seisout(1:ns,3) = s3(1:ns)
end subroutine
