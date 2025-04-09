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

l=fix(i/2);

for m=1:l
	frewind(file_Blist);
	pause(60);
	for n=1:m
		file=fgetl(file_Blist);
	end
	file
	pack;
	A = fopen(file,'rb');
	B1 = fread(A,'single');
	B(1)=B(1)+B1(1);
	B(2:end) = B(2:end) + B1(2:end)*B1(1);
	fclose(A);
	clear file A;
end
fclose(file_Blist);
A=fopen('B_1.bin','wb');
fwrite(A,B,'single');
