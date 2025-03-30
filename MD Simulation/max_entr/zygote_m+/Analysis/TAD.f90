Program main
implicit none
integer i,j,k,ios,nargs,Nres
real,allocatable :: s(:)
real,allocatable :: f(:,:)
real :: tmp(3)
character(len=200) buffer,Inputfile
integer :: unit
parameter(unit=100000)
real :: n,d

nargs=iargc()
if(nargs==2) then
        call getarg(1,Inputfile)
        call getarg(2,buffer)
        read(buffer,*) Nres
        allocate(s(Nres))
        allocate(f(Nres,Nres))
else
        write(*,*) "Wrong Input!"
        write(*,*) "Program ContactPfile Nres"
        goto 1000
endif

ios=0
open(1,file=Inputfile,status="old")
do while(.true.)
        read(1,*,iostat=ios) (tmp(j),j=1,3)
        if(ios/=0) exit
        i=INT(ANINT(tmp(1)))
        j=INT(ANINT(tmp(2)))
        f(i,j)=tmp(3)
        f(j,i)=tmp(3)
enddo
close(1)

s=1.0E10
do i=1,Nres
        n=0
        d=0
        do k=1,1000000/unit
                j=i+k
                if(j<=Nres.and.j>=1) then
                        n=n+f(i,j)
                endif
        enddo
        do k=-1000000/unit,-1
                j=i+k
                if(j<=Nres.and.j>=1) then
                        d=d+f(i,j)
                endif
        enddo
        if(d/=0.and.n/=0) s(i)=1.0/log(2.0)*log(n/d)
enddo

open(2,file="TAD.dat",status="replace")
do i=1,Nres
        if(s(i)==1.0E10) then
                cycle
        else
                write(2,*) i,s(i)
        endif
enddo
close(2)

deallocate(s)
deallocate(f)
1000 END program
