%f file average
file_flist=fopen('Filelist_f');

for i=1:1
	file=fgetl(file_flist);
    	f=dlmread(file);
	Nframes_S=f(1);
        tmp=size(f);
        k=tmp(:,1)-1;
end

fdata=zeros(k,1);

Nframes=0;
frewind(file_flist);
while ~feof(file_flist)
	file=fgetl(file_flist)
        f=dlmread(file);
	Nframes_S=f(1);
	fnew=f(2:k+1,2);
	fdata=fdata+fnew*Nframes_S;
	Nframes=Nframes+Nframes_S;
end

fclose(file_flist);
fdata = fdata/Nframes;
save('f.dat','fdata','-ascii');
clear;

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
Pnew=zeros(k,3);

Nframes=0;
frewind(file_Plist);
while ~feof(file_Plist)
	file=fgetl(file_Plist)
        P=dlmread(file);
	Nframes_S=P(1);
	Pnew=P(2:end,:);
	clear P;
	Pdata=Pdata+Pnew.*Nframes_S;
	Nframes=Nframes+Nframes_S;
end
fclose(file_Plist);
Pdata = Pdata/Nframes;
save('Probability.dat','Pdata','-ascii');
clear;
