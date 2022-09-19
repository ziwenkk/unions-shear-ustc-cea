! 
      program crossdr4
      parameter(nmx=600000,ntimes=10,nrmx=6000000)
      real qra(nmx),qdec(nmx),ra(nmx),dec(nmx),rra(nrmx),rdec(nrmx)
      real qred(nmx),red(nmx),rred(nrmx),qwei(nmx),wei(nmx),rwei(nrmx)
      external wfib_dr72, wfib_tdf
      character*200 aaa,datf,xif
      real(kind=8) ra8,dec8,red8
      real(kind=8),allocatable :: trra(:,:),trdec(:,:),trred(:,:)

      zmin=0.01
      zmax=1
      abmmin=-25
      abmmax=-15
      ifib=1

c     get arguments
      narg=iargc()
      if(narg.lt.1)then
         write(*,'(a)')'Usage: crossdr4 datf [xif=xif ifib=ifib '
     +        //'zmin=zmin zmax=zmax abmmin=abmmin abmmax=abmmax] '
         stop
      end if

      call getarg(1,datf)
      xif=datf

      call get_unit(iu)
      open(iu,file=datf,status='old',iostat=ierr)
      close(iu)
      if(ierr.ne.0)then
         write(*,'(a)')'Error: failed to open '//trim(datf)//'. Stop!'
         stop
      end if
      
      if(narg.gt.1)then
         do i=2,narg
            call getarg(i,aaa)
            k1=index(aaa,'=')
            k2=index(aaa,' ')-1
            if(aaa(1:k1-1).eq.'xif')then
               if(k2.ge.k1+1)then
                  read(aaa(k1+1:k2),'(a)')xif
               else
                  write(*,'(a)')'Warning: wrong xif; set to datf'
               end if
            else if(aaa(1:k1-1).eq.'ifib')then
               read(aaa(k1+1:k2),*)ifib
               if(ifib.ne.1.and.ifib.ne.0)then
                  ifib=1
                  write(*,'(a)')'Warning: wrong ifib; set to 1'
               end if
            else if(aaa(1:k1-1).eq.'abmmin')then
               read(aaa(k1+1:k2),*)abmmin
               if(abmmin.lt.-30.or.abmmin.ge.-10)then
                  abmmin=-30
                  write(*,'(a)')'Warning: wrong abmmin; set to -30'
               end if
            else if(aaa(1:k1-1).eq.'abmmax')then
               read(aaa(k1+1:k2),*)abmmax
               if(abmmax.le.-30.or.abmmax.gt.-10)then
                  abmmax=-10
                  write(*,'(a)')'Warning: wrong abmmax; set to -10'
               end if
            end if
         end do
      end if

      call get_unit(iu)
      open(iu,file=trim(xif)//'.log',status='new',iostat=ierr)
      close(iu)
      if(ierr.ne.0) then
         write(*,'(a)')'Error: '//trim(xif)//'.log'
     +        //' already exists. Stop!'
         stop
      end if

      print*,trim(xif)
      print*,zmin,zmax,abmmin,abmmax,ifib

      ! read reference sample

      call get_unit(iu)
      open(iu,file=datf,status='old')

      n=0
      do i=1,nmx
         read(iu,*,end=88)t1,t2,t3

         if(t3.lt.zmin.or.t3.gt.zmax)cycle
             n=n+1
             ra(n)=t1
             dec(n)=t2
             red(n)=t3
             wei(n)=1.0

      end do
 88   close(iu)
      write(*,'(a,i7)')'# ref=',n

      !read the sample being studied

      call get_unit(iu)
      open(iu,file=datf,status='old')

      nq=0
      do i=1,nmx
         read(iu,*,end=99) t1,t2,t3

         if(t3.lt.zmin.or.t3.gt.zmax)cycle

             nq=nq+1
             qra(nq)=t1
             qdec(nq)=t2
             qred(nq)=t3
             qwei(nq)=1.0

      end do
 99   close(iu)
      write(*,'(a,i7)')'# obj=',nq

      ! read random sample
      allocate(trra(n,ntimes))
      allocate(trdec(n,ntimes))
      allocate(trred(n,ntimes))

      call random_sky(0.0,360.0,-90.0,90.0,ntimes,n,trra,trdec,trred)

      j=0
      do it=1,ntimes
         do k=1,n

            ra8=trra(k,it)
            dec8=trdec(k,it)
            red8=trred(k,it)

            j=j+1
            rra(j)=ra8
            rdec(j)=dec8
            rred(j)=red8
            rwei(j)=wei(k)

         end do
      end do
      nr=j

      deallocate(trra)
      deallocate(trdec)
      deallocate(trred)
      write(*,'(a,i7)')'# ran=',nr

c--------------------------------------------------------------------------------------

      omegam=0.308
      omegal=1.0-omegam
      H0=100.

      nbp=9
      seppmin=0.5
      dsepp=0.2
      logsepp=1

      nbs=nbp
      sepsmin=seppmin
      dseps=dsepp
      logseps=logsepp

      nbv=40
      sepvmin=0.0
      dsepv=1.0
      logsepv=0

c--------------------------------------------------------------------------------------

      if(ifib.eq.1)then
         call cross_tpcf(omegam,omegal,H0,nq,qra,qdec,qred,qwei, 
     x        n,ra,dec,red,wei,nr,rra,rdec,rred,rwei,wfib_dr72, 
     x        nbs,sepsmin,dseps,logseps,nbp,seppmin,dsepp,logsepp,
     x        nbv,sepvmin,dsepv,logsepv,xif)
      else
         call cross_tpcf(omegam,omegal,H0,nq,qra,qdec,qred,qwei, 
     x        n,ra,dec,red,wei,nr,rra,rdec,rred,rwei,wfib_tdf, 
     x        nbs,sepsmin,dseps,logseps,nbp,seppmin,dsepp,logsepp,
     x        nbv,sepvmin,dsepv,logsepv,xif)
      end if
      
      end

