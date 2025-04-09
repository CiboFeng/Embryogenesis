Program main
implicit none
integer i,j,k,ios,nargs,Nres
real,allocatable :: ContactP(:,:),ContactPCG(:,:),Contact_a(:,:),Prob(:,:)
real tmp(3)
character(len=200) ContactPfile,buffer,Contact_afile,Probfile
real cutoff1,cutoff2
!parameter(cutoff1=0.05,cutoff2=0.01)
parameter(cutoff1=0.008,cutoff2=0.01)
real cutoff_GD
parameter(cutoff_GD=1000)
integer ResInter
parameter(ResInter=2)
real a
real unit
parameter(unit=0.1)

nargs=iargc()
if(nargs==4) then
        call getarg(1,ContactPfile)
        call getarg(2,buffer)
        call getarg(3,Contact_afile)
        read(buffer,*) Nres
        call getarg(4,Probfile)
        allocate(ContactP(Nres,Nres))
        allocate(ContactPCG(Nres,Nres))
        allocate(Contact_a(Nres,Nres))
        allocate(Prob(Nres,Nres))
else
        write(*,*) "Wrong Input!"
        write(*,*) "Usage: Program ContactPfile Nres Contact_afile Probfile"
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

ios=0
Prob=0
open(1,file=Probfile,status="old")
do while(.true.)
        read(1,*,iostat=ios) (tmp(j),j=1,3)
        if(ios/=0) exit
        if((tmp(2)-tmp(1))>=ResInter) then
                        i=INT(ANINT(tmp(1)))
                        j=INT(ANINT(tmp(2)))
                        Prob(i,j)=tmp(3)
        endif
enddo
close(1)

ios=0
Contact_a=-0.01
open(1,file=Contact_afile,status="old")
do while(.true.)
        read(1,*,iostat=ios) (tmp(j),j=1,3)
        if(ios/=0) exit
        if((tmp(2)-tmp(1))>=ResInter) then
                        i=INT(ANINT(tmp(1)))
                        j=INT(ANINT(tmp(2)))
                        Contact_a(i,j)=tmp(3)
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

!do i=3,Nres,3
!        do j=3,Nres,3
do i=1,Nres,2
        do j=1,Nres,2
!                if(abs(j-i)<=cutoff_GD/unit) then
!                 if(abs(j-i)<=20/unit) then
                        if(ContactPCG(i,j)==0) then
                        if(ContactP(i,j)>=0.001) then
                        if((ContactP(i,j)-Prob(i,j))/(ContactP(i,j))>0.7.or.(ContactP(i,j)-Prob(i,j))/(ContactP(i,j))<-0.7) then
                                ContactPCG(i,j)=ContactP(i,j)
                        endif
                                endif
                        endif
!               endif
        enddo
enddo

open(2,file="Contact_a_0.dat",status="replace")
open(3,file="Probability.dat",status="replace")
!do i=1,Nres
!        do j=1,Nres
!                if(ContactPCG(i,j)/=0) then
!                       call random_number(a)
!                        a=1.0
!                       write(2,*) i,j,Contact_a(i,j)
!                endif
!        enddo
!enddo
do i=1,Nres
        do j=1,Nres
                if(ContactPCG(i,j)/=0) then
                       if(Contact_a(i,j)/=-0.01) then 
!                               write(2,*) i,j,Contact_a(i,j)
!                               write(3,*) i,j,ContactPCG(i,j)
                        endif
                endif
        enddo
enddo
do i=1,Nres
        do j=1,Nres
                if(ContactPCG(i,j)/=0) then
                       if(Contact_a(i,j)==-0.01) then 
                               write(2,*) i,j,Contact_a(i,j)
                               write(3,*) i,j,ContactPCG(i,j)
                        endif
                endif
        enddo
enddo
close(2)
close(3)
deallocate(ContactP)
1000 End program
