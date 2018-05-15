u=linspace(0,1,200);
v=linspace(0,4,200);
[u,v]=meshgrid(u,v);
hidden off;
%surf(u,1./u,30*v,'FaceColor','red','EdgeColor','none','FaceLighting','phong'); hold on;
%axis([0 1 0 4 0 100]);
tur=(1./u).*(5+7*u.*v)./(4.+5.*u.*v+3.*u.*u.*v.*v);
stat=(1./(sqrt(1./v)-sqrt((1./v)-u))).^2;
tur(logical(imag(stat)))=nan;
%tur(logical(tur>100))=nan;
stat(logical(imag(stat)))=nan;
%stat(logical(stat>100))=nan;
surf(u,v,tur,'FaceColor','green','EdgeColor','none','FaceLighting','phong','FaceAlpha',1); hold on;
%alpha(1)
axis([0 1 0 4 0 100]);
xlim([0 1]);
ylim([0 4]);
zlim([0 100]);
surf(u,v,stat,'FaceColor','blue','EdgeColor','none','FaceLighting','phong','FaceAlpha',.4); hold on;
%alpha(.4)
axis([0 1 0 4 0 100]);
xlim([0 1]);
ylim([0 4]);
zlim([0 100]);
surf(.5*cos(u.*2.*pi).*sin(v.*pi/4.)+.5,.5*sin(u.*2.*pi).*sin(v.*pi/4.)+.03,10.*cos(v.*pi/4.)+21.6,'FaceColor','yellow','EdgeColor','none','FaceLighting','phong','FaceAlpha',.8); hold off;
%alpha(.4)
axis([0 1 0 4 0 100]);
xlim([0 1]);
ylim([0 4]);
zlim([0 100]);
xlabel 'e/p', ylabel 'd/p', zlabel 'v/u';
title 'Parameter Range'
camlight 'right';

a=findobj(gcf); % get the handles associated with the current figure

allaxes=findall(a,'Type','axes');
alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');

set(allaxes,'FontName','Arial','FontWeight','Bold','LineWidth',2,...
    'FontSize',14);
set(alllines,'Linewidth',2);
set(alltext,'FontName','Arial','FontWeight','Bold','FontSize',14);