Program main
implicit none
integer :: i,j,k,ios,Nres,nargs
real,allocatable :: CP(:),NCP(:)
real :: tmp(3)
integer :: unit
parameter(unit=100000)
character(len=200) Inputfile,buffer

nargs=iargc()
if(nargs==2) then
        call getarg(1,Inputfile)
        call getarg(2,buffer)
        read(buffer,*) Nres
        allocate(CP(Nres))
        allocate(NCP(Nres))
else
        write(*,*) "Wrong Input!"
        write(*,*) "Program ContactPfile Nres"
        goto 1000
endif

ios=0
CP=0
NCP=0
open(1,file=Inputfile,status="old")
do while(.true.)
       read(1,*,iostat=ios) (tmp(j),j=1,3)
       if(ios/=0) exit
       i=INT(ANINT(tmp(1)))
       j=INT(ANINT(tmp(2)))
       k=j-i
       if(k>0) then
               CP(k)=CP(k)+tmp(3)   
               NCP(k)=NCP(k)+1
       endif
enddo
close(1)

open(2,file="CP_GD.dat",status="replace")
do i=1,1
        write(2,*) 1.0*i*unit/1000000,1.0
enddo
do i=2,Nres
        write(2,*) 1.0*i*unit/1000000,CP(i)/NCP(i)
enddo
close(2)

deallocate(CP)
deallocate(NCP)
1000 END program
