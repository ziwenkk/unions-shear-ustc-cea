      subroutine random_sky(ramin,ramax,decmin,decmax,
     +        ntimes,nobj,rra,rdec,rred)
      real*8 rra(nobj,ntimes),rdec(nobj,ntimes),rred(nobj,ntimes)
      parameter(ntmx=10,nptmx=2000000)
      real*8 ra(nptmx,ntmx),dec(nptmx,ntmx),red(nptmx,ntmx)
      real*8 tra,tdec,tred
      integer npt(ntmx)
      common /points/npt,ra,dec,red
      character*80 aaa
      data ifirst/1/
      data isd/-830906/

      if(ntimes.gt.ntmx) stop 'ntimes much be <= 10'

      !---------------------------------------------
      ! If called for the first time,
      if(ifirst.eq.1) then
         ifirst=0

         write(*,'(a)') 'Reading LSS catalogues...'

         do i=1,ntimes
            call int2str(i-1,1,aaa)
            call get_unit(iunit)
            open(iunit,file='./rans/ran_'//aaa(1:1)
     +           //'.dat',status='old')
            npt(i)=0
            do j=1,nptmx
               read(iunit,*,end=88)tra,tdec,tred
               if(ramax<ramin) then
                  if(tra.gt.ramax.and.tra.lt.ramin) cycle
               else
                  if(tra.gt.ramax.or.tra.lt.ramin) cycle
               end if
               if(tdec.gt.decmax.or.tdec.lt.decmin) cycle

                   npt(i)=npt(i)+1
                   k=npt(i)
                   ra(k,i)=tra
                   dec(k,i)=tdec
                   red(k,i)=tred

            end do
 88         close(iunit)
         end do
      end if

      print*,'ifirst=',ifirst

      do i=1,ntimes
         print*,i,npt(i)
         nsel=0
         do j=1,npt(i)

            nsel=nsel+1
            rra(nsel,i)=ra(j,i)
            rdec(nsel,i)=dec(j,i)
            rred(nsel,i)=red(j,i)
            if(nsel.eq.nobj) exit

         end do
         print*,'nsel=',nsel
      end do

      return
      end
