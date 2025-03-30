Program main
implicit none
integer i,j,k,ios
real(kind=4),allocatable :: fdata(:,:),f(:)
real(kind=4) :: tmp
character(len=200) fdatafile

call getarg(1,fdatafile)

ios=0
k=0
open(1,file=fdatafile,status="old")
do while(.true.)
        read(1,*,iostat=ios) tmp
        if(ios/=0) exit
        k=k+1
enddo
rewind(1)
allocate(fdata(k,k))
allocate(f(k))
ios=0
k=0
do while(.true.)
        read(1,*,iostat=ios) tmp
        if(ios/=0) exit
        k=k+1
        f(k)=tmp
enddo
close(1)

open(4,file="fdata.bin",status="replace",access="STREAM",form="unformatted",action="write")
do i=1,k
        write(4) (f(i)*f(j),j=1,k)
enddo

deallocate(fdata)
deallocate(f)
END
