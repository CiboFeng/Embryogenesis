Program main
implicit none
integer i,j,k,l,ios,nargs,Nres
real :: tmp(3)
real,allocatable :: Prob(:,:)
character(len=200) Probfile,buffer
real,allocatable :: IS(:)
real ISmean
integer Start
parameter(Start=5)
real,allocatable :: Delta(:)
real Delta_left,Delta_right
integer Delta_bin
parameter(Delta_bin=1)

nargs=iargc()
if(nargs==2) then
        call getarg(1,Probfile)
        call getarg(2,buffer)
        read(buffer,*) Nres
        allocate(Prob(Nres,Nres))
        allocate(IS(Nres))
        allocate(Delta(Nres))
else
        write(*,*) "Wrong Input!"
        write(*,*) "Usage: Program Probfile Nres"
        goto 1000
endif

ios=0
Prob=0
open(1,file=Probfile,status="old")
do while(.true.)
        read(1,*,iostat=ios) (tmp(j),j=1,3)
        if(ios/=0) exit
        i=INT(ANINT(tmp(1)))
        j=INT(ANINT(tmp(2)))
        Prob(i,j)=tmp(3)
enddo
close(1)

do i=1,Nres
        do j=1,i-2
                Prob(i,j)=Prob(j,i)
        enddo
enddo

IS=0
ISmean=0
do i=Start+1,Nres-Start
        do k=-Start,-1
                do l=1,Start
                        IS(i)=IS(i)+Prob(i+k,i+l)
                enddo
        enddo
        ISmean=ISmean+IS(i)
enddo

ISmean=ISmean/(Nres-Start-Start)

Delta=0
open(7,file="InsulationScore.dat",status="replace")
open(8,file="Delta.dat",status="replace")
do i=Start+1,Nres-Start
        IS(i)=log(IS(i)/ISmean)/log(2.0)
        write(7,*) i,IS(i)
enddo
close(7)
do i=Start+1,Nres-Start
        Delta_left=0
        do k=-Delta_bin,0
                Delta_left=Delta_left+IS(i+k)
        enddo
        Delta_right=0
        do k=0,Delta_bin
                Delta_right=Delta_right+IS(i+k)
        enddo
        Delta(i)=(Delta_left-Delta_right)/(Delta_bin)
        write(8,*) i,Delta(i)
enddo
close(8)

deallocate(Prob)
deallocate(IS)
1000 END Program
