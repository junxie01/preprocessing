subroutine do_norm(sig,fre1,fre2,sachead)
use sacio
implicit none
type(sac_head):: sachead
integer,parameter :: nmax=4000000
integer i,nerr,npts,halfwidth
real dum,fre1(4),fre2(4)
real temp,dt
real x1(nmax),x2(nmax),x3(nmax)
real tx1(nmax),tx2(nmax),tx3(nmax)
real sig(nmax,3)
!
dt=sachead%delta
npts=sachead%npts
!write(*,*)'npts=',npts
halfwidth=int(1.0/fre1(2)/dt/2.0) ! Besen et al. (2007) half of the maximum period of the passband filter
!halfwidth=20                    ! from Yingjie Yang's code
halfwidth=int(64/dt)
!write(*,*)fre(2),fre(3),dt,npts
x1(1:npts)=sig(1:npts,1)
x2(1:npts)=sig(1:npts,2)
x3(1:npts)=sig(1:npts,3)
tx1(1:npts)=sig(1:npts,1)
tx2(1:npts)=sig(1:npts,2)
tx3(1:npts)=sig(1:npts,3)
!write(*,*)'Band pass filter the waveform !'
call filter(x1,fre1,dt,npts) ! filter the waveform
call filter(x2,fre1,dt,npts)
call filter(x3,fre1,dt,npts)
!write(*,*)'Smooth the waveform'
call smooth(x1,npts,halfwidth) ! get the running absolute average waveform
call smooth(x2,npts,halfwidth)
call smooth(x3,npts,halfwidth)

call filter(tx1,fre2,dt,npts) ! filter the waveform
call filter(tx2,fre2,dt,npts)
call filter(tx3,fre2,dt,npts)
do i=1,npts 
   temp=max(x1(i),x2(i),x3(i))
   sig(i,1)=sig(i,1)/temp
   sig(i,2)=sig(i,2)/temp
   sig(i,3)=sig(i,3)/temp
enddo
end subroutine
