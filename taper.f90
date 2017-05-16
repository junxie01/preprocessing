subroutine tapering(sig,npts,width)
real sig(npts),width,f0,f1,pi,omega
real valuea
integer npts,i
integer ntab
! hanning widow
f0=0.5
f1=0.5
ntab=width*(npts+1)
pi=atan(1.0)*4.0
omega=pi/ntab
do i=1,ntab
       valuea= f0 - f1*cos( omega*( i-1 ) )
       sig(i)=sig(i)*valuea
       sig(npts-i+1)=sig(npts-i+1)*valuea
enddo
end subroutine
