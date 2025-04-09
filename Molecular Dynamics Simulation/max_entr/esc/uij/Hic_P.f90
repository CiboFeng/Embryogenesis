Program main
integer i,j,k,ios,nargs
integer :: Nlocus
integer :: Nlocus_low,Nlocus_high
character(len=200) Hicdata,buffer,Hicdata_P,Hic_Pdata
real,allocatable :: n(:),p(:,:),tmp(:)

nargs=iargc()
if(nargs==3) then
        call getarg(1,Hicdata)
        Hic_Pdata=trim(Hicdata) // "_P.dat"
        Hicdata_P=trim(Hicdata) // "_P"
        call getarg(2,buffer)
        read(buffer,*) Nlocus_low
        call getarg(3,buffer)
        read(buffer,*) Nlocus_high
        Nlocus=Nlocus_high-Nlocus_low+1
        allocate(n(1:Nlocus))
        allocate(p(1:Nlocus,1:Nlocus))
        allocate(tmp(1:Nlocus_high))
else
        write(*,*) "Wrong Input!"
        write(*,*) "Usage: Program Hicdata Nlocus_low Nlocus_high"
        goto 1000
endif

ios=0
i=0
open(1,file=Hicdata,status="old")
do while(.true.)
        read(1,*,iostat=ios) (tmp(j),j=1,Nlocus_high)
        if(ios/=0) exit
        i=i+1
        if(i>=Nlocus_low.and.i<=Nlocus_high) then
                do j=1,Nlocus
                        p(i-Nlocus_low+1,j)=tmp(j+Nlocus_low-1)
                enddo
        endif
enddo
close(1)

!do i=1,Nlocus
!        do k=-1,1,2
!                j=i+k
!                if(j<=Nlocus.and.j>=1) then
!                        p(i,j)=1.0E-10
!                endif
!        enddo
!enddo

n=-1
do i=1,Nlocus
        do k=-4,4
!               if(k==1.or.k==-1.or.k==0) cycle
                if(k==0) cycle
                j=i+k
!               if(j/=i.and.j<=Nlocus.and.j>=1) then
                if(j<=Nlocus.and.j>=1) then
                        n(i)=max(n(i),p(i,j))
                endif
        enddo
enddo

do i=1,Nlocus
        if(n(i)<=0) then
                write(*,*) i,n(i)
        endif
enddo

do i=1,Nlocus
        do j=1,Nlocus
                if(min(n(i),n(j))/=0) then
                        if(abs(i-j)<=1) then
                                p(i,j)=1
                        else
                                p(i,j)=min(1.0,p(i,j)/(min(n(i),n(j))))
                        endif
                else
                        p(i,j)=0
                endif
        enddo
enddo

open(2,file=Hic_Pdata,status="replace")
do i=1,Nlocus
        do j=1,Nlocus
                write(2,*) i,j,p(i,j)
        enddo
enddo
close(2)

open(3,file=Hicdata_P,status="replace")
do i=1,Nlocus
        write(3,*) (p(i,j),j=1,Nlocus)
enddo
close(3)

deallocate(n)
deallocate(p)
deallocate(tmp)
1000 END program
