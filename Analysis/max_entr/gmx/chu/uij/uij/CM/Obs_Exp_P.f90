Program main
integer i,j,k,ios,nargs
integer :: Nlocus
character(len=200) Hicdata,buffer
real,allocatable :: p(:,:)
real,allocatable :: pobs(:,:),nobs(:,:),po_e(:,:)
real,allocatable :: pexp(:),nexp(:)
real :: tmp(3)
integer :: NCG,CG,m,n,d
parameter(CG=10)
integer :: x,y

nargs=iargc()
if(nargs==2) then
        call getarg(1,Hicdata)
        call getarg(2,buffer)
        read(buffer,*) Nlocus
        allocate(p(1:Nlocus,1:Nlocus))
        NCG=Nlocus/CG
        allocate(pobs(1:NCG,1:NCG))
        allocate(po_e(1:NCG,1:NCG))
        allocate(nobs(1:NCG,1:NCG))
        allocate(pexp(0:NCG-1))
        allocate(nexp(0:NCG-1))
else
        write(*,*) "Wrong Input!"
        write(*,*) "Usage: Program Pdata Nlocus"
        goto 1000
endif

ios=0
i=0
open(1,file=Hicdata,status="old")
do while(.true.)
        read(1,*,iostat=ios) (tmp(j),j=1,3)
        if(ios/=0) exit
                i=INT(ANINT(tmp(1)))
                j=INT(ANINT(tmp(2)))
                p(i,j)=tmp(3)
enddo
close(1)
do i=1,Nlocus
        do j=1,i+1
                if(j>=1.and.j<=Nlocus) then
                        if(abs(i-j)<=1) then
                                p(i,j)=1.0
                        else
                                p(i,j)=p(j,i)
                        endif
                endif
        enddo
enddo

pobs=0
nobs=0
do i=1,NCG
        do m=-4,5
                x=(i-1)*CG+m
                if(x>=1.and.x<=Nlocus) then
                        do j=i,NCG
                                do n=-4,5
                                        y=(j-1)*CG+n
                                        if(y>=1.and.y<=Nlocus) then
                                                pobs(i,j)=pobs(i,j)+p(x,y)
                                                nobs(i,j)=nobs(i,j)+1
                                        endif
                                enddo
!                                pobs(i,j)=pobs(i,j)/nobs(i,j)
                        enddo
                endif
        enddo
enddo

do i=1,NCG
        do j=1,i-1
                if(j>=1.and.j<=NCG) then
                        pobs(i,j)=pobs(j,i)
                        nobs(i,j)=nobs(j,i)
                endif
        enddo
enddo

do i=1,NCG
        do j=1,NCG
                pobs(i,j)=pobs(i,j)/nobs(i,j)
        enddo
enddo

open(2,file="Pobs.dat",status="replace")
do i=1,NCG
       write(2,*) (pobs(i,j),j=1,NCG)
enddo
close(2)

pexp=0
nexp=0
do i=1,NCG
        do d=0,NCG-1
                x=INT(ANINT(d*0.95))+i
                y=INT(ANINT(d*1.05))+i
                do j=x,y
                        if(j>=1.and.j<=NCG) then
                                if(j>=1.and.j<=NCG) then
                                        pexp(d)=pexp(d)+pobs(i,j)
                                        nexp(d)=nexp(d)+1
                                endif
                        endif
                enddo
        enddo
enddo

do d=0,NCG-1
        pexp(d)=pexp(d)/nexp(d)
enddo

po_e=0
do i=1,NCG
        do j=1,NCG
                d=abs(i-j)
                if(pexp(d)/=0) then
                        if(d==0) then
                                po_e(i,j)=1.0
                        else
                                po_e(i,j)=pobs(i,j)/pexp(d)
                        endif
                else
                        po_e(i,j)=0
                endif
        enddo
enddo

open(3,file="Po_e.dat",status="replace")
do i=1,NCG
        write(3,*) (po_e(i,j),j=1,NCG)
enddo
close(3)

open(4,file="Pexp.dat",status="replace")
do d=0,NCG-1
        write(4,*) d,pexp(d)
enddo
close(4)

deallocate(p)
deallocate(pobs)
deallocate(nobs)
deallocate(pexp)
deallocate(nexp)
deallocate(po_e)
1000 END program
