Program main
implicit none
integer i,j,k,ios,nargs,N
real :: x(9999),y(9999),z(9999),tmp(5)
integer :: atom_res(9999)
character(len=3) :: res(9999)
character(len=2) :: atom(9999)
real :: bonds,angles,dihedrals(9999),pairsAlpha(9999,9999),a,b,c
real :: kb,ktheta
real :: e_n,e_r,sigma_r
parameter(e_n=1.0)
parameter(kb=1.0,ktheta=2.0)
parameter(e_r=1.0,sigma_r=0.4)
character(len=200) PDBfile,ProteinName,TOPfile,contactAlphafile
character(len=79) buffer
real :: pi
parameter(pi=3.1415926)

nargs=iargc()
if(nargs==3) then
	call getarg(1,PDBfile)
        call getarg(2,contactAlphafile)
	call getarg(3,ProteinName)
	TOPfile=trim(ProteinName) // ".top"
else
	write(*,*) "Wrong Input!"
	write(*,*) "Usage: Program PDBfile contactAlphafile ProteinName"
	goto 1000
endif

ios=0
x=0
y=0
z=0
res=""
atom=""
i=0
open(1,file=PDBfile,status="old")
do while(.true.)
	read(1,"(A)",iostat=ios) buffer
	if(ios/=0) exit
	if(buffer(1:4)=="ATOM") then
		i=i+1
		read(buffer(25:26),*) j
		read(buffer(32:38),*) x(i)
		read(buffer(40:46),*) y(i)
		read(buffer(48:54),*) z(i)
		res(i)=buffer(18:20)
		atom(i)=buffer(14:15)
		atom_res(i)=j
		N=i
	endif
enddo
close(1)
ios=0
pairsAlpha=0
open(2,file=ContactAlphafile,status="old")
do while(.true.)
        read(2,*,iostat=ios) (tmp(j),j=1,3)
        if(ios/=0) exit
        i=INT(ANINT(tmp(1)))
        j=INT(ANINT(tmp(2)))
        pairsAlpha(i,j)=tmp(3)
enddo
close(2)
open(7,file=TOPfile,status="replace")
!head of TOPfile
write(7,"(A)") ";Polymer MODEL GENERATION by Xiakun"
write(7,"(A)") "[ defaults ]"
write(7,"(A)") ";nbfunc comb-rule gen-pairs"
write(7,"(A)") "1  1  no"
write(7,"(A)")
write(7,"(A)") "[ atomtypes ]"
write(7,"(A)") ";name mass charge ptype c6 c12"
write(7,"(A,1X,E11.5E2,1X,E11.5E2)") "CA 1.000 0.000 A", 0.0, e_r
write(7,"(A)")
write(7,"(A)") "[ moleculetype ]"
write(7,"(A)") ";name nrexcl"
write(7,"(A)") "Protein 1 ;;;; 1-2 and 1-3 are in Nonbonded terms!!!"
!write(7,"(A)") "Protein 2 ;;;; 1-3 are in Nonbonded terms!!!"
write(7,"(A)")
!atoms
write(7,"(A)") "[ atoms ]"
write(7,"(A)") ";nr type resnr residue atom cgnr charge mass"
do i=1,9999
	if(res(i)=="") cycle
	write(7,"(I4,1X,A2,1X,I4,1X,A3,1X,A2,1X,I4,1X,A)") i,"CA",i,res(i),"CA",i,"0.000 1.000" 
enddo
write(7,"(A)")
!pairs
write(7,"(A)") "[ pairs ]"
do i=1,999
        do j=1,999
                if(pairsAlpha(i,j)==0) cycle
                write(7,"(I3,1X,I3,3X,A1,1X,E12.5E2,1X,E12.5E2)") &
&               i,j,"1",0.0,e_n*pairsAlpha(i,j)
        enddo
enddo
write(7,"(A)")
!bonds
write(7,"(A)") "[ bonds ]"
write(7,"(A)") ";ai aj func kb table"
do i=1,N-1
	bonds=sqrt((x(i+1)-x(i))**2+(y(i+1)-y(i))**2+(z(i+1)-z(i))**2)/10
	write(7,"(I4,1X,I4,3X,A1,1X,A1,1X,E11.5E2,1X,A)") &
&			i,i+1,"8","0",kb,"; table files"
enddo
write(7,"(A)")
!exclusions
!!!no exclusions
!angles
write(7,"(A)") "[ angles ]"
write(7,"(A)") ";ai aj ak func Ka table"
do i=2,N-1
	j=i-1
	k=i+1
	write(7,"(I4,1X,I4,1X,I4,3X,A1,1X,A1,1X,E11.5E2,1X,A)") &
&		j,i,k,"8","0",ktheta,"; table files"
enddo
write(7,"(A)")
!dihedrals
!!!no dihedrals
write(7,"(A)") "[ system ]"
write(7,"(A)") ";name"
write(7,"(A)") "Protein"
write(7,"(A)")
write(7,"(A)") "[ molecules ]"
write(7,"(A)") ";name #molec"
write(7,"(A)") "Protein 1"
close(7)

1000 END program
