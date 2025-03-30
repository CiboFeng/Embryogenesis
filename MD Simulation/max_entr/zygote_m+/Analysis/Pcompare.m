P = dlmread('../Probability.Mixed.dat');

n = size(P);
Nres=sqrt(n(:,1));

%PHiC = zeros(Nres,Nres);
%PSim = zeros(Nres,Nres);

for k=1:n(:,1)
    i=P(k,1);
    j=P(k,2);
    if i>j+1 
%        PHiC(i,j)=P(k,3);
        m=0.5*(j-1)*(Nres-2+Nres-(j-1)-1)+(i-(j)-1);
        y(m,1)=P(k,3);
    end
    if i<j-1
%        PSim(i,j)=P(k,3);
        n=0.5*(i-1)*(Nres-2+Nres-(i-1)-1)+(j-(i)-1);
        x(n,1)=P(k,3);
    end
end

%k=0;
%for i=1:Nres
%    for j=i+2:Nres
%        k=k+1;
%        y(k,1)=PHiC(j,i);
%        x(k,1)=PSim(i,j);
%    end
%end

X = [ones(length(x),1) x];
b = X\y
yCalc2 = X*b;
Rsq2 = 1 - sum((y - yCalc2).^2)/sum((y - mean(y)).^2)
xx = linspace(0.00,1,50);
XX = [ones(length(xx),1), xx'];
yyCalc2 = XX*b;
figure('Position', [10 10 650 600])
p1=plot(x,y,'.');hold on;
p3=plot(xx,yyCalc2,'-','Linewidth',8); hold on;
axis([-0.05 1 -0.05 1]);
set(gca,'fontsize',18);
xlabel('Probability (Sim)','fontsize',25);
ylabel('Probability (HiC)','fontsize',25);

a2 = sprintf('%4.2f',b(1));
a3 = sprintf('%4.2f',b(2));
aa2 = sprintf('%4.2f',Rsq2);

%legend([p3,],{['y=',a3,'x+',a2,' R^2=',aa2]},'Location','Southeast');
%legend('boxoff');

param = [ b(1) b(2) Rsq2 ];
save('Pcompare.param','param','-ascii');

print ('-dpdf', '-r300', 'Pcompare.pdf')
