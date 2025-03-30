Program main
implicit none
integer i,j,k,ios,nargs,Nres
real,allocatable :: ContactP(:,:),ContactPCG(:,:)
real tmp(3)
character(len=200) ContactPfile,buffer
real cutoff1,cutoff2
parameter(cutoff1=0.008,cutoff2=0.01)
!parameter(cutoff1=0.20,cutoff2=0.04)
real cutoff_GD
parameter(cutoff_GD=100)
integer ResInter
parameter(ResInter=2)
real a
real unit
parameter(unit=0.1)

nargs=iargc()
if(nargs==2) then
        call getarg(1,ContactPfile)
        call getarg(2,buffer)
        read(buffer,*) Nres
        allocate(ContactP(Nres,Nres))
        allocate(ContactPCG(Nres,Nres))
else
        write(*,*) "Wrong Input!"
        write(*,*) "Usage: Program ContactPfile Nres"
        goto 1000
endif

ios=0
ContactP=0
ContactPCG=0
open(1,file=ContactPfile,status="old")
do while(.true.)
        read(1,*,iostat=ios) (tmp(j),j=1,3)
        if(ios/=0) exit
        if((tmp(2)-tmp(1))>=ResInter) then
                        i=INT(ANINT(tmp(1)))
                        j=INT(ANINT(tmp(2)))
                        ContactP(i,j)=tmp(3)
        endif
enddo
close(1)

ContactPCG=0
do i=2,Nres,2
        do j=2,Nres,2
                if(ContactP(i,j)>=cutoff1) then
                        if(abs(j-i)<=cutoff_GD/unit) then
                                ContactPCG(i,j)=ContactP(i,j)
                        endif
                endif
        enddo
enddo

do i=1,Nres,4
        do j=1,Nres,4
                if(ContactP(i,j)>=cutoff2) then
                        if(abs(j-i)<=cutoff_GD/unit) then
                                ContactPCG(i,j)=ContactP(i,j)
                        endif
                endif
        enddo
enddo

open(2,file="Contact_a_0.dat",status="replace")
open(3,file="Probability.dat",status="replace")
do i=1,Nres
        do j=1,Nres
                if(ContactPCG(i,j)/=0) then
!                       call random_number(a)
!                        a=1.0
                        write(2,*) i,j,-0.01
                        write(3,*) i,j,ContactPCG(i,j)
                endif
        enddo
enddo
close(2)
close(3)
deallocate(ContactP)
1000 End program
