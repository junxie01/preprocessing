! ==========================================================
! Function filter4. Broadband filreting.
! ==========================================================
! Parameters for filter4 function:
! Input parameters:
! f1,f2   - low corner frequences, f2 > f1, Hz, (double)
! f3,f4   - high corner frequences, f4 > f3, Hz, (double)
! npow    - power of cosine tapering,  (int)
! dt      - sampling rate of seismogram in seconds, (double)
! n       - number of input samples, (int)
! seis_in - input array length of n, (float)
! Output parameters:
! seis_out - output array length of n, (float)
! ==========================================================
subroutine filter(sig,fre,dt,npts)
implicit none
integer,parameter :: FFTW_FORWARD=-1,FFTW_BACKWARD=1
integer,parameter :: FFTW_ESTIMATE=64,FFTW_MEASURE=1
integer,parameter :: nmax=4000000
integer*4 k,nk,npts,ns,npow
real*4    fre(4),dt,sig(nmax),sigout(nmax)
real(8)    dom
integer(8)   plan1,plan2
complex(8)   czero,s(nmax),sf(nmax)
! ---
czero = dcmplx(0.0d0,0.0d0)
ns=1
npow=0
do while(ns.lt.npts)
   ns=ns*2
   npow=npow+1
enddo
dom = dble(1.0/dt/ns)
nk = ns/2+1
s=czero
sf=czero
s(1:npts) = dcmplx(dble(sig(1:npts)),0d0)
!write(*,*)' do fft'
!call clogc(npow,s,1,dt)
call dfftw_plan_dft_1d(plan1,ns,s,sf,FFTW_FORWARD, FFTW_ESTIMATE)
call dfftw_execute(plan1)
call dfftw_destroy_plan(plan1)

sf(1) = sf(1)/2.0d0
sf(nk) = dcmplx(dreal(sf(npts)),0.0d0)
sf(1:nk)=sf(1:nk)*dble(dt)
call taperf(sf,fre,dom,nk,npow)
s=czero
sf(nk+1:ns)=czero
!write(*,*)'taper the waveform done'
! make forward FFT for seismogram: sf ==> s
call dfftw_plan_dft_1d(plan2,ns,sf,s,FFTW_BACKWARD, FFTW_ESTIMATE)
call dfftw_execute(plan2)
call dfftw_destroy_plan(plan2)
!call clogc(npow,s,-1,dt)
sig(1:npts)=2.0*real(dreal(s(1:npts)))/ns/dt
return
end subroutine
