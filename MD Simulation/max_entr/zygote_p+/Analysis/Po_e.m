data = dlmread('normal.dat');

rho = corr(data);
[ n n ] = size(rho);
[row, col] = find(isnan(rho));
tmp = rho;
rowsize = size(row);
for i=1:rowsize(1);
  tmp(row(i),col(i)) = 0.0;
end
rho=tmp;
[ V D ] = eig(rho);
x=linspace(1,n,n);
vector = V(:,1);
save('vector.dat','vector','-ascii');
