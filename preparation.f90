program preparation ! do 3 domponent preprocessing
use sacio
implicit none
type(sac_head):: sachead
integer nmax,nstmax
parameter (nmax=4000000,nstmax=1000)
logical ext
real fre(4),dt
real sig(nmax,3)
real seisout(nmax,3)
real stla,stlo,evla,evlo
character (3)com(3),co
character (20)year_day,nd
character (10)sta(nstmax),kst
character (180)sac(3),command
character (180)out_w(3),out(3)
character (180)dirinn,dirout,output,input,list
integer nerr,npts
integer jday,npt2,icheck
integer begday,endday,ic
integer iy,id,ih,ist,nsamp,is
integer dhour,nst,i,nh,j,nnh,nlen
integer nzhour,nzmin,nzsec,nzmsec
integer year_b,year_e,day_b,day_e
!
if (iargc().ne.1)then
   write(*,*)'Usage: preparation param.dat '
   write(*,*)'param.dat:'
   write(*,*)'staion_list'
   write(*,*)'year_b day_b year_e day_e'
   write(*,*)'dhour f1 f2 COM'
   write(*,*)'data_dir'
   call exit(-1)
endif
call getarg(1,input)
open(10,file=input)
read(10,'(a180)')list
read(10,*)year_b,day_b,year_e,day_e
read(10,*)dhour,fre(2),fre(3),co
read(10,'(a180)')dirinn
close(10)

com(1)=trim(co)//'Z'
com(2)=trim(co)//'N'
com(3)=trim(co)//'E'
fre(1)=0.95*fre(2)
fre(4)=1.05*fre(3)
i=1
open(11,file=list)          ! read in station list
do i=1,nstmax
   read(11,*,err=13,end=13) sta(i)
enddo
13 close(11)
nst=i-1
write(*,'("There are ",i0," stations")')nst
write(*,'("The sac files are stored in",1x,1a,"/")')trim(dirinn)
write(*,'("The components are: ",1a3,1x,1a3,1x,1a3)')(com(i),i=1,3)
write(*,'("The frequency band is ",4f9.5)')( fre(i),i=1,4)
write(*,'("Do preprocess from ",i0,"/",i3.3," to ",i0,"/",i3.3)')year_b,day_b,year_e,day_e
nh=24/dhour            ! number of segments per day
do iy=year_b,year_e            ! loop over year
   jday=365
   if(mod(iy,4).eq.0.and.mod(iy,100).ne.0.or.mod(iy,400).eq.0)jday=366
   endday=day_e
   if(iy.ne.year_e)endday=jday
   begday=day_b
   if(iy.ne.year_b)begday=1
   do id=begday,endday         ! loop over day
      write(year_day,'(i0,"_",i3.3)')iy,id
      do is=1,nst              ! loop over station 
         do ih=1,nh            ! loop over hour
            nnh=(ih-1)*dhour 
            write(nd,'("_",i2.2)')nnh
            do ic=1,3
               sac(ic)=trim(dirinn)//'/'//trim(year_day)//'/'//trim(year_day)&
               //trim(nd)//'_'//trim(sta(is))//'_'//trim(com(ic))//'.SAC'
               write(*,*)trim(sac(ic))
               !out(ic)=trim(sac(ic))//'_bp'
               out_w(ic)=trim(sac(ic))//'_sw'
               inquire(file=sac(ic),exist=ext)
               if(.not.ext)then
                   icheck=ic-1
                   exit
               endif
               icheck=ic
            enddo
            sig=0
            if(icheck.eq.3)then
               do ic=1,3
                  !call readsac(sac(ic),sig(:,ic),beg,dt,npts,stla,stlo,nzhour,nzmin,nzsec,nzmsec,nerr)
                  call read_sachead(trim(sac(ic)),sachead,nerr)
                  call read_sac(trim(sac(ic)),sig(:,ic),sachead,nerr)
               enddo
               write(*,*)'dt=',sachead%delta
                                      ! time domain normalisation
               write(*,*)'Do time domain normalisation'
               call do_norm(sig,fre,sachead) 
               write(*,*)'Do frequency domain whitening'
               call do_whiten(sachead,fre,sig,seisout)
               do ic=1,3
                  !call tapering(seisout(:,ic),npts,0.05)
                  !call wrsac(out_w(ic),seisout(:,ic),npts,dt,stla,stlo,nnh,0,0,0)
                  call write_sac(out_w(ic),seisout(:,ic),sachead,nerr)
                  ! call wrsac(out_w(ic),sig(:,ic),npts,dt,stla,stlo,nnh,0,0,0)
               enddo
            endif   ! if all three component exists
         enddo      ! end loop over hour
      enddo                    ! end loop over station
   enddo                           ! end loop over day
enddo                                 ! end loop over year
write(*,*)"All done for one station preprocess!"
end program
