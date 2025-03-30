a=dlmread(['Probability.Mixed.dat']);
unit=0.1;
x=a(:,1)*unit; y=a(:,2)*unit; z=a(:,3);
Nres=812;
w=linspace(1,Nres,Nres)*unit;
v=linspace(1,Nres,Nres)*unit;
[xx,yy]=meshgrid(w,v);
zz = griddata(x,y,z,xx,yy,'linear');

crit = 0.01;
data = zz;
data(:,:) = 0.01;
[row,col] = find(zz>crit);
rowsize = size(row);
for i=1:rowsize(1);
  data(row(i),col(i)) = zz(row(i),col(i));
end
zz=data;

colormap(flipud(hot));
colorbar;
surf(xx,yy,zz);hold on;

%caxis([0 15]);
%[C,h]=contourf(xx,yy,zz,8,'Linewidth',0.3);hold on
%[C,h]=contour(xx,yy,zz,15,'Linewidth',2);hold on
%contourc(xx,yy,zz);
%plot3(x1,y1,z1,'Linewidth',4);

colorbar;
%caxis([-1,2]);
%title(['D_{i}'],'fontsize',30)
%colormap(pink);
%colormap(hot);
%grid on;
axis([0 Nres*unit 0 Nres*unit]);
%axis([0 1 0 20]);
%caxis([-2 6]);
caxis([0.01 1.00]);
set(gca,'colorscale','log')
%set(h,'LineWidth',6);
%set(h,'ShowText','on','TextStep',get(h,'LevelStep')*10)
%surf(xx,yy,zz);
%plot(b(:,1),b(:,2),'k-','linewidth',1);
%hcb=colorbar('YTick',0:2:10);
%set(hcb,'YTickMode','manual');
xlabel('Genomic Distance (Mb)','fontsize',20);
ylabel('Genomic Distance (Mb)','fontsize',20);
%zlabel('Distance (nm)','fontsize',30);
%zlabel('Distance of SideChain (A)','fontsize',20);
%axis([100 550 0 250]);
grid on;
grid minor
%caxis([0 60])
%title('','fontsize',20)
%title([num2str(k),' mol/L'])
shading flat;
%shading interp
%shading interp
view(0,90);
h=get(gcf);
set(gca,'XTick',0:10:Nres*unit,'fontsize',18);
set(gca,'YTick',0:10:Nres*unit,'fontsize',18);
gcb=colorbar;
%set(gcb,'YTick',0:0.4:2.4,'fontsize',30);
box off;
%set(gca,'Fontsize',15)

title('Iteration=Alpha (RATIO%)','fontsize',20);
%print ('-dpsc', '-r300', 'TCI.eps')
print ('-dpng', '-r300', 'Contact_a_ALPHA.png')
%print ('-dpdf', '-r300', 'Contact.pdf')
%print ('-dpng', '-r300', 'TCI.eps')

