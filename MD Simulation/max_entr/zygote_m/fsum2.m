%f file average
file_flist=fopen('Filelist_f');

for i=1:1
	file=fgetl(file_flist);
    	f=dlmread(file);
	Nframes_S=f(1);
        tmp=size(f);
        k=tmp(:,1)-1;
	clear f;
end
fclose(file_flist);
clear file_flist file Nframes_S tmp;

%B file average
file_Blist=fopen('Filelist_B');
B=zeros(k*k+1,1);
B1=zeros(k*k+1,1);
i=0;
while ~feof(file_Blist)
	fgetl(file_Blist);
	i=i+1;
end

for m=1:i
	frewind(file_Blist);
	pause(60);
	for n=1:m
		file=fgetl(file_Blist);
	end
	file
	pack;
	A = fopen(file,'rb');
	B1 = fread(A,'single');
	B = B + B1;
	fclose(A);
	clear file A;
end
Bdata = reshape(B(2:end)./B(1),k,k);
clear B;
clear B1;
fclose(file_Blist);

%fdata file
A = fopen('fdata.bin','rb');
B = fread(A,'single');
fclose(A);

Bdata = Bdata/Nframes - reshape(B,k,k);
%save('B.dat','Bdata','-ascii');
clear B;

%P file average
file_Plist=fopen('Filelist_P');

for i=1:1
	file=fgetl(file_Plist);
    	P=dlmread(file);
	Nframes_S=P(1);
        tmp=size(P);
        k=tmp(:,1)-1;
	clear P;
end

Pdata=zeros(k,3);

Nframes=0;
frewind(file_Plist);
while ~feof(file_Plist)
	file=fgetl(file_Plist)
        P=dlmread(file);
	Nframes_S=P(1);
	Pnew=P(2:k+1,1:3);
	clear P;
	Pdata=Pdata+Pnew*Nframes_S;
	clear Pnew;
	Nframes=Nframes+Nframes_S;
end
fclose(file_Plist);
Pdata = Pdata/Nframes;
save('Probability.dat','Pdata','-ascii');
clear Pdata;

temp = TEMP;
ResInter = 2;
kB = 8.314/1000;
beta = 1.0/(temp*kB);
Lambda = LAMBDA;

c = dlmread("f.dat");
f = c(:,1);

d = dlmread("Probability_ND.dat"); % Probability from Native Dynamics
n = size(d);

k=0;
ratio_n=0;
ratio_d=0;
for i = 1:n(1,1)
        if((d(i,2)-d(i,1))>=ResInter)
            k=k+1;
            fexp(k,1)=d(i,3);
            ratio_n=ratio_n+abs(f(k)-fexp(k));
            ratio_d=ratio_d+fexp(k);
        end
end

ratio = ratio_n/ratio_d;
save("tol.dat",'ratio',"-ascii");

Contact_alpha_prevfile=dlmread("ContactalphaFILE");
Contact_alpha_prev=Contact_alpha_prevfile(:,3);

contact_alpha = 1.0/beta*Bdata\(f-fexp);
%contact_alpha = 1.0/beta*pinv(B)*(f-fexp);

k=0;
for i = 1:n(1,1)
        if((d(i,2)-d(i,1))>=ResInter)
            k=k+1;
            Contact_alphafile(k,1)=d(i,1);
            Contact_alphafile(k,2)=d(i,2);
            Contact_alpha_abs(k,1)=d(i,1);
            Contact_alpha_abs(k,2)=d(i,2);
            Contact_alphafile(k,3)=contact_alpha(k)*Lambda+Contact_alpha_prev(k);
%	    if(Contact_alphafile(k,3)>=1000.0)  
%                Contact_alphafile(k,3)=1000.0;
%            end
%            if(Contact_alphafile(k,3)<=-1000.0) 
%                Contact_alphafile(k,3)=-1000.0;
%            end
            Contact_alpha_abs(k,3)=contact_alpha(k);
        end
end

save("Contact_a.dat","Contact_alphafile","-ascii");
save("Contact_a_abs.dat","Contact_alpha_abs","-ascii");
