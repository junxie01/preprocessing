!  smoothing routine  in frquency domain 
subroutine smoothf(dom,nk,halfn,fre,sf,smf)
integer,parameter :: nmax=4000000
integer  halfn,nk,i
complex  sf(nmax)
real     smf(nmax)
real*4   summ, fre(4)
real*8   dom,f
smf = 0.0
do i = 1,nk
   f = dble((i-1))*dom
   if( f .ge. dble(fre(1)) .and. f .le. dble(fre(4)) ) then
      smf(i) = sum(cabs(sf(i-halfn:i+halfn)))/(2.*halfn+1.)
   endif
enddo
return
end subroutine
