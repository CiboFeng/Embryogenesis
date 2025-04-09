temp = TEMP;
ResInter = 2;
kB = 8.314/1000;
beta = 1.0/(temp*kB);
Lambda = LAMBDA;

B = dlmread("B.dat");
%i_max = max(a(:,1));
%i_min = min(a(:,1));

%for i=i_min:i_max
%        for j=i_min:i_max
%            B(i,j)=a((i-1)*(i_max-i_min+1)+j,3);
%        end
%end

c = dlmread("f.dat");
f = c(:,2);

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

contact_alpha = 1.0/beta*B\(f-fexp);
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
