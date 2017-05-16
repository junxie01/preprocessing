program preparation 
! do 3 domponent preprocessing
! 2017/03/08
use sacio
implicit none
type(sac_head):: sachead(3)
integer,parameter:: nmax=4000000,nstmax=1000
logical ext
real fre1(4),dt,fre2(4)
real sig(nmax,3),sigo(nmax,3)
real*8 dom
complex seisout(nmax,3)
real stla,stlo,evla,evlo
character (2)net(nstmax)
character (3)com(3),co
character (20)year_day,nd
character (10)sta(nstmax),kst
character (180)sacf(3),command
character (180)out_w(3),out(3)
character (180)dir_data,dirout,output,input,list
integer nerr,npts,iseg
integer multpt,dsec,dseg,nseg
integer jday,npt2,icheck
integer begday,endday,ic
integer iy,id,ih,ist,nsamp,is,nn
integer dhour,nst,i,nh,j,nnh,nlen
integer nzhour,nzmin,nzsec,nzmsec
integer year_b,year_e,day_b,day_e
!
if (iargc().ne.1)then
   write(*,*)'Usage: preparation param.dat '
   write(*,*)'param.dat:'
   write(*,*)'staion_list'
   write(*,*)'year_b day_b year_e day_e'
   write(*,*)'dsec multpt com fa fb f1 f2'
   write(*,*)'data_dir'
   call exit(-1)
endif
call getarg(1,input)
open(10,file=input)
read(10,'(a180)')list
read(10,*)year_b,day_b,year_e,day_e
read(10,*)dsec,multpt,co,fre1(2),fre1(3),fre2(2),fre2(3)
read(10,'(a180)')dir_data
close(10)

com(1)=trim(co)//'Z'
com(2)=trim(co)//'N'
com(3)=trim(co)//'E'
fre1(1)=0.95*fre1(2)
fre1(4)=1.05*fre1(3)
fre2(1)=0.95*fre2(2)
fre2(4)=1.05*fre2(3)
open(11,file=list)          ! read in station list
do i=1,nstmax
   read(11,*,err=13,end=13) net(i),sta(i)
enddo
13 close(11)
nst=i-1

dseg=int((1-real(multpt)/100.0)*dsec) ! the left points without overlapping
nseg=int((86400-dsec)/dseg)+1         ! number of segments per day.
do iy=year_b,year_e                   ! loop over year
   jday=365
   if(mod(iy,4).eq.0.and.mod(iy,100).ne.0.or.mod(iy,400).eq.0)jday=366
   endday=day_e
   if(iy.ne.year_e)endday=jday
   begday=day_b
   if(iy.ne.year_b)begday=1
   do id=begday,endday                ! loop over day
      write(year_day,'(i0,"_",i3.3)')iy,id
      do is=1,nst                     ! loop over station 
         do iseg=1,nseg               ! loop over segments
            nzhour=(iseg-1)*dseg/3600
            nzmin=mod((iseg-1)*dseg,3600)/60
            nzsec=mod(mod((iseg-1)*dseg,3600),60)
            do ic=1,3                 ! loop over components
               !write(*,*)com(ic)
               write(sacf(ic),'(1a,"/",1a,"/",1a,"_",i2.2,"_",i2.2,"_",i2.2,"_",1a,"_",1a,"_",1a,".SAC")')&
               trim(dir_data),trim(year_day),trim(year_day),nzhour,nzmin,nzsec,trim(net(is)),trim(sta(is)),trim(com(ic))
               !out(ic)=trim(sac(ic))//'_bp'
               !out_w(ic)=trim(sacf(ic))//'_sw'
               inquire(file=sacf(ic),exist=ext)
               icheck=0
               !if(.not.ext)exit
               !write(*,*)trim(sacf(ic))
               icheck=ic
            enddo
            sig=0
            if(icheck.ne.3)cycle                   ! check whether all three components exist
            !do ic=1,3                              ! read all data
            !   write(*,*)trim(sacf(ic))
            !enddo
            do ic=1,3                              ! read all data
               call read_sachead(trim(sacf(ic)),sachead(ic),nerr)
               call read_sac(trim(sacf(ic)),sig(:,ic),sachead(ic),nerr)
               if(nerr.eq.-1)exit
            enddo
            if(nerr.eq.-1)cycle
            !write(*,*)'dt=',sachead%delta
            write(*,*)'Do time domain normalisation'
            call do_norm(sig,fre1,fre2,sachead(1))          ! time domain normalisation
            write(*,*)'Do frequency domain whitening'
            call do_whiten(nn,dom,sachead(1),fre2,sig,seisout)! frequency whitening
            write(*,*)'Frequency domain whitening done'
            do ic=1,3
               !call tapering(seisout(:,ic),npts,0.05) ! do taper
               sachead(ic)%npts=nn
               sachead(ic)%delta=sngl(dom)
               out_w(ic)=trim(sacf(ic))//'.re'
               sigo(1:nn,ic)=real(seisout(1:nn,ic))
               call write_sac(out_w(ic),sigo(:,ic),sachead(ic),nerr)
               out_w(ic)=trim(sacf(ic))//'.im'
               sigo(1:nn,ic)=aimag(seisout(1:nn,ic))
               call write_sac(out_w(ic),sigo(:,ic),sachead(ic),nerr)
               call write_sac(trim(sacf(ic))//'.amp',cabs(seisout(1:nn,ic)),sachead(ic),nerr)
            enddo
         enddo      ! end loop over segments
      enddo         ! end loop over station
   enddo            ! end loop over day
enddo               ! end loop over year
write(*,*)"All done for one station pre-processing!"
end program
