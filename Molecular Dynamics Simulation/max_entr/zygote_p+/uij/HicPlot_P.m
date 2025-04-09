%Hicdata = dlmread('PrEC_100000_iced_chr14_dense.matrix_P');
Hicdata = dlmread('PC3_100000_iced_chr14_dense.matrix_P');
n = size(Hicdata,1);
x = linspace(1,n,n)*100/1000;
y = linspace(1,n,n)*100/1000;

[xx yy] = meshgrid(x,y);

crit = 0.00;
zz = Hicdata;
Hicdata(:,:) = 0.00;
[row,col] = find(zz>crit);
rowsize = size(row);
for i=1:rowsize(1);
  Hicdata(row(i),col(i)) = zz(row(i),col(i));
end

surf(xx,yy,Hicdata);

%colormap(map);
colormap(flipud(hot));
colorbar;

surf(xx,yy,zz);hold on;
cbh=colorbar;

caxis([0.0 1.00]);
%set(gca,'colorscale','log')

tickLocations = [0.0,0.2,0.4,0.6,0.8,1.0]; % change to whatever you want
tickLabels    = {'0.0','0.2','0.4','0.6','0.8','1.0'}; % change to whatever you want
%tickLocations = [0.0,0.5,1.0]; % change to whatever you want
%tickLabels    = {'0.0','0.5','1.0'}; % change to whatever you want
set(cbh,'YTick',tickLocations,'YTickLabel',tickLabels)

xlabel('Genomic Distance (Mb)','fontsize',25);
ylabel('Genomic Distance (Mb)','fontsize',25);
shading flat;
%shading interp
%shading interp
view(0,90);
h=get(gcf);
axis([0 n*100/1000 0 n*100/1000]);
set(gca,'XTick',0:10:n*100/1000,'fontsize',18);
set(gca,'YTick',0:10:n*100/1000,'fontsize',18);
%gcb=colorbar;
%set(gcb,'YTick',0:0.4:2.4,'fontsize',30);
box off;

print('P.pdf','-dpdf')

surf(xx,yy,zz);hold on;
cbh=colorbar;

caxis([0.0 1.00]);
%set(gca,'colorscale','log')

tickLocations = [0.0,0.2,0.4,0.6,0.8,1.0]; % change to whatever you want
tickLabels    = {'0.0','0.2','0.4','0.6','0.8','1.0'}; % change to whatever you want
%tickLocations = [0.0,0.5,1.0]; % change to whatever you want
%tickLabels    = {'0.0','0.5','1.0'}; % change to whatever you want
set(cbh,'YTick',tickLocations,'YTickLabel',tickLabels)

xlabel('Genomic Distance (Mb)','fontsize',25);
ylabel('Genomic Distance (Mb)','fontsize',25);
shading flat;
%shading interp
%shading interp
view(0,90);
h=get(gcf);
axis([20 40 20 40]);
set(gca,'XTick',0:5:n*100/1000,'fontsize',18);
set(gca,'YTick',0:5:n*100/1000,'fontsize',18);
%gcb=colorbar;
%set(gcb,'YTick',0:0.4:2.4,'fontsize',30);
box off;

print('P_details.pdf','-dpdf')