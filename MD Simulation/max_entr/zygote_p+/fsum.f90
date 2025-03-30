Program main
implicit none
integer i,j,k,nargs,ios_F,ios,Nres
real(kind=8),allocatable :: B1(:,:),f(:),tmp_B(:)
real :: Nframes_s,Nframes
real :: tmp(3)
character(len=200) Filelist_f,Filelist_B,Filelist_P,inputfile,ContactPfile,buffer
real(kind=8),allocatable :: ContactP(:,:)
integer :: ResInter
parameter(ResInter=2)

nargs=iargc()
if(nargs==5) then
       call getarg(1,Filelist_f)
       call getarg(2,Filelist_B)
       call getarg(3,Filelist_P)
       call getarg(4,ContactPfile)
       call getarg(5,buffer)
       read(buffer,*) Nres
       allocate(ContactP(Nres,Nres))
else
        write(*,*) "Wrong Input!"
        write(*,*) "Usage: Program Filelist_f Filelist_B Filelist_P ContactPfile Nres"
        goto 1000
endif

ios=0
k=0
open(3,file=ContactPfile,status="old")
do while(.true.)
        read(3,*,iostat=ios) (tmp(j),j=1,3)
        if(ios/=0) exit
                if((tmp(2)-tmp(1))>=ResInter) then
                        k=k+1
                endif
enddo
allocate(f(k))
allocate(B1(k,k))
allocate(tmp_B(k))
close(3)

ios_F=0
Nframes=0
f=0
open(1,file=Filelist_f,status="old")
do while(.true.)
       read(1,"(A)",iostat=ios_F) inputfile
       if(ios_F/=0) exit
       ios=0
       open(2,file=inputfile,status="old")
       read(2,*) Nframes_s
       Nframes=Nframes+Nframes_s
       do while(.true.)
                read(2,*,iostat=ios) (tmp(j),j=1,2)
                if(ios/=0) exit
                i=INT(ANINT(tmp(1)))
                f(i)=f(i)+tmp(2)*Nframes_s
       enddo
       close(2)
enddo
close(1)

do i=1,k
        f(i)=f(i)/Nframes
enddo

ios_F=0
Nframes=0
B1=0
open(1,file=Filelist_B,status="old")
do while(.true.)
       read(1,"(A)",iostat=ios_F) inputfile
       if(ios_F/=0) exit
       ios=0
       open(2,file=inputfile,status="old")
       read(2,*) Nframes_s
       Nframes=Nframes+Nframes_s
       i=0
       do while(.true.)
                read(2,*,iostat=ios) (tmp_B(j),j=1,k)
                if(ios/=0) exit
                i=i+1
                do j=i,k
                        B1(i,j)=B1(i,j)+tmp_B(j)*Nframes_s
                enddo
       enddo
       close(2)
enddo
ios_F=0
Nframes=0
ContactP=0
open(1,file=Filelist_P,status="old")
do while(.true.)
       read(1,"(A)",iostat=ios_F) inputfile
       if(ios_F/=0) exit
       ios=0
       open(2,file=inputfile,status="old")
       read(2,*) Nframes_s
       Nframes=Nframes+Nframes_s
       do while(.true.)
                read(2,*,iostat=ios) (tmp(j),j=1,3)
                if(ios/=0) exit
                i=INT(ANINT(tmp(1)))
                j=INT(ANINT(tmp(2)))
                ContactP(i,j)=ContactP(i,j)+tmp(3)*Nframes_s
        enddo
        close(2)
enddo

open(3,file="f.dat",status="replace")
do i=1,k
        write(3,*) i,f(i)
enddo

open(4,file="B.dat",status="replace")
do i=1,k
!        do j=1,k
!                B1(i,j)=B1(i,j)/Nframes
!                write(4,*) i,j,B1(i,j)-f(i)*f(j)
!        enddo
        write(4,*) (B1(j,i)/Nframes-f(i)*f(j),j=1,i-1), (B1(i,j)/Nframes-f(i)*f(j),j=i,k)
enddo

open(5,file="Probability.dat",status="replace")
do i=1,Nres
        do j=1,Nres
                if((j-i)>=ResInter) then
                        if(ContactP(i,j)/=0) then
                        ContactP(i,j)=ContactP(i,j)/Nframes
                        write(5,*) i,j,ContactP(i,j)
                        endif
                endif
        enddo
enddo

deallocate(f)
deallocate(B1)
deallocate(ContactP)
deallocate(tmp_B)
1000 End Program
