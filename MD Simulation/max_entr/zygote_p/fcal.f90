Program main
implicit none
integer i,j,k,ios,nargs,Nres,critical
real,allocatable :: x(:),y(:),z(:)
real(kind=4) :: Nframes,d
real(kind=4) :: Rc,unit,miu
parameter(Rc=1.5*0.4,miu=3.0,unit=10) ! Unit is in nm
!parameter(Rc=1.78*0.4,miu=2.5*3.22,unit=10) ! Unit is in nm
character(len=200) PDBfile,buffer,ContactPfile
integer,allocatable :: p2f(:,:)
real(kind=4),allocatable :: f(:),p(:),B1(:,:)
real(kind=4),allocatable :: ContactP(:,:)
real :: tmp(3)
integer :: ResInter
parameter(ResInter=2)

nargs=iargc()
if(nargs==3) then
        call getarg(1,PDBfile)
        call getarg(2,ContactPfile)
        call getarg(3,buffer)
        read(buffer,*) Nres
        allocate(x(Nres))
        allocate(y(Nres))
        allocate(z(Nres))
        allocate(p2f(Nres,Nres))
        allocate(ContactP(Nres,Nres))
else
        write(*,*) "Wrong Input!"
        write(*,*) "Usage: Program PDBfile ContactPfile NumberofResidues"
        goto 1000
endif

ios=0
k=0
open(2,file=ContactPfile,status="old")
do while(.true.)
        read(2,*,iostat=ios) (tmp(j),j=1,3)
        if(ios/=0) exit
                if((tmp(2)-tmp(1))>=ResInter) then
                        k=k+1
                endif
enddo
allocate(f(k))
allocate(p(k))
allocate(B1(k,k))
rewind(2)
ios=0
k=0
p2f=-1
do while(.true.)
        read(2,*,iostat=ios) (tmp(j),j=1,3)
        if(ios/=0) exit
                if((tmp(2)-tmp(1))>=ResInter) then
                        k=k+1
                        i=INT(ANINT(tmp(1)))        
                        j=INT(ANINT(tmp(2)))
                        f(k)=tmp(3)
                        p2f(i,j)=k
                endif
enddo
close(2)

ios=0
Nframes=0
B1=0
f=0
ContactP=0
open(1,file=PDBfile,status="old")
do while(.true.)
        read(1,"(A)",iostat=ios) buffer
        if(ios/=0) exit
        if(buffer(1:5)=="MODEL") then
                Nframes=Nframes+1
                p=0
                critical=0
                do i=1,Nres
                        read(1,"(A)") buffer
                        read(buffer(32:38),*) x(i)
                        read(buffer(40:46),*) y(i)
                        read(buffer(48:54),*) z(i)
                        if(x(i)<=-100.or.y(i)<=-100.or.z(i)<=-100) then
                                critical=1
                                exit
                        endif
                        if(x(i)>=100.or.y(i)>=100.or.z(i)>=100) then
                                critical=1
                                exit
                        endif
                enddo
                if(critical==1) then
                        Nframes=Nframes-1
                        cycle
                endif
                do i=1,Nres
                        do j=1,Nres
                                if((j-i)>=ResInter) then
                                        d=sqrt((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)*0.1
                                        ContactP(i,j)=ContactP(i,j)+0.5*(1+tanh((Rc-d)*miu))
                                        if(p2f(i,j)==-1) cycle
                                        p(p2f(i,j))=0.5*(1+tanh((Rc-d)*miu))
                                        f(p2f(i,j))=f(p2f(i,j))+p(p2f(i,j))
                                endif
                        enddo
                enddo
                do i=1,k
                        do j=i,k
                                B1(i,j)=B1(i,j)+p(i)*p(j)
                        enddo
                enddo
        endif
enddo
close(1)

open(3,file="f.dat",status="replace")
write(3,*) Nframes
do i=1,k
        write(3,*) i,f(i)/Nframes
enddo
close(3)

open(4,file="B.bin",status="replace",access="STREAM",form="unformatted",action="write")
write(4) Nframes
do i=1,k
!        do j=1,k
!                if(i>j) B1(i,j)=B1(j,i)
!                write(4,*) i,j,B1(i,j)/Nframes
!        enddo
        write(4) (B1(j,i)/Nframes,j=1,i-1),(B1(i,j)/Nframes,j=i,k)
enddo
close(4)

open(5,file="Probability.dat",status="replace")
write(5,*) Nframes
do i=1,Nres
        do j=1,Nres
                if((j-i)>=ResInter) then
!                        if(ContactP(i,j)/Nframes/=0) then
                                write(5,*) i,j,ContactP(i,j)/Nframes
!                        endif
                endif
        enddo
enddo
close(5)

deallocate(x)
deallocate(y)
deallocate(z)
deallocate(p2f)
deallocate(f)
deallocate(p)
deallocate(B1)
deallocate(ContactP)

1000 End Program
