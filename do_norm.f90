subroutine do_norm(sig,fre,sachead)
use sacio
implicit none
type(sac_head):: sachead
integer,parameter :: nmax=4000000
integer i,nerr,npts,halfwidth
real dum,fre(4)
real temp,dt
real x1(nmax),x2(nmax),x3(nmax)
real sig(nmax,3)
!
dt=sachead%delta
npts=sachead%npts
!write(*,*)'npts=',npts
halfwidth=int(1.0/fre(2)/dt/2.0) ! Besen et al. (2007) half of the maximum period of the passband filter
!halfwidth=20                    ! from Yingjie Yang's code
halfwidth=int(64/dt)
!write(*,*)fre(2),fre(3),dt,npts
!write(*,*)'Band pass filter the waveform !'
call filter(sig(:,1),x1,fre,dt,npts) ! filter the waveform
call filter(sig(:,2),x2,fre,dt,npts)
call filter(sig(:,3),x3,fre,dt,npts)
!write(*,*)'Smooth the waveform'
call smooth(x1,npts,halfwidth) ! get the running absolute average waveform
call smooth(x2,npts,halfwidth)
call smooth(x3,npts,halfwidth)
do i=1,npts 
   temp=max(x1(i),x2(i),x3(i))
   sig(i,1)=sig(i,1)/temp
   sig(i,2)=sig(i,2)/temp
   sig(i,3)=sig(i,3)/temp
enddo
end subroutine
