% fileN=fopen('Uspot.txt','r');
% A=fscanf(fileN,'%g %g %g %g', [4 inf]);
% fclose(fileN);
% figure, errorbar(A(1,:),A(3,:)/A(3,1),A(4,:)/A(3,1),'o')
% title('Red Spot size');
% xlabel('IPTG concentration');
% ylabel('Spot Area (a.u.)');
% set(gca,'xscale','log');
% 
% a=findobj(gcf); % get the handles associated with the current figure
% 
% allaxes=findall(a,'Type','axes');
% alllines=findall(a,'Type','line');
% alltext=findall(a,'Type','text');
% 
% set(allaxes,'FontName','Arial','FontWeight','Bold','LineWidth',2,...
%     'FontSize',14);
% set(alllines,'Linewidth',2);
% set(alltext,'FontName','Arial','FontWeight','Bold','FontSize',14)
% 
% 
% set(allaxes,'FontName','Arial','FontWeight','Bold','LineWidth',2,...
%     'FontSize',14);
% set(alllines,'Linewidth',2);
% set(alltext,'FontName','Arial','FontWeight','Bold','FontSize',14)
% print('spotSizes','-dpdf');
% 
% 
% 
% 
% 
% fileN=fopen('UtotArea.txt','r');
% A=fscanf(fileN,'%g %g %g %g', [4 inf]);
% fclose(fileN);
% fileN=fopen('CtotArea.txt','r');
% B=fscanf(fileN,'%g %g %g %g', [4 inf]);
% fclose(fileN);
% figure;
% [hAx,hLine1,hLine2]=plotyy(B(1,:),B(4,:)/B(4,7),A(1,:),A(3,:)/A(3,1),'semilogx');
% %title('Red Spot Total Area');
% xlabel('IPTG concentration (M)');
% ylabel(hAx(2),'Total red spot area (a.u.)');
% set(hAx(2),'Ycolor', [0 0 0]);
% set(hAx(1), 'Ycolor', [0 0 0]);
% set(hAx(1),'YLim',[0 1.1])
% set(hAx(1),'YTick',[0:.2:1])
% set(hAx(2),'YLim',[0 1.1])
% set(hAx(2),'YTick',[0:.2:1])
% set(hLine1, 'Marker', 'd');
% set(hLine1, 'Color', 'g');
% set(hLine1, 'MarkerFaceColor', 'g');
% set(hLine1, 'LineStyle', 'none');
% set(hLine2, 'Marker', 'd');
% set(hLine2, 'Color',  'r');
% set(hLine2, 'MarkerFaceColor', 'r');
% set(hLine2, 'LineStyle', 'none');
% ylabel(hAx(1),'Average GFP (a.u.)');
% 
% a=findobj(gcf); % get the handles associated with the current figure
% 
% allaxes=findall(a,'Type','axes');
% alllines=findall(a,'Type','line');
% alltext=findall(a,'Type','text');
% 
% set(allaxes,'FontName','Arial','FontWeight','Bold','LineWidth',2,...
%     'FontSize',14);
% set(alllines,'Linewidth',2);
% set(alltext,'FontName','Arial','FontWeight','Bold','FontSize',14)
% 
% 
% set(allaxes,'FontName','Arial','FontWeight','Bold','LineWidth',2,...
%     'FontSize',14);
% set(alllines,'Linewidth',2);
% set(alltext,'FontName','Arial','FontWeight','Bold','FontSize',14)
% print('totSize','-dpdf');
% 
% 
% 
% 
% 
% fileN=fopen('UtotArea.txt','r');
% A=fscanf(fileN,'%g %g %g %g', [4 inf]);
% fclose(fileN);
% fileN=fopen('CtotArea.txt','r');
% B=fscanf(fileN,'%g %g %g %g', [4 inf]);
% fclose(fileN);
% figure;
% [hAx,hLine1,hLine2]=plotyy(B(1,:),B(3,:)/B(3,end),A(1,:),A(3,:)/A(3,1),'semilogx');
% %title('Red Spot Total Area');
% xlabel('IPTG concentration (M)');
% ylabel(hAx(2),'Total red spot area (a.u.)');
% set(hAx(2),'Ycolor', [0 0 0]);
% set(hAx(1), 'Ycolor', [0 0 0]);
% set(hAx(1),'YLim',[0 1.1])
% set(hAx(1),'YTick',[0:.2:1])
% set(hAx(2),'YLim',[0 1.1])
% set(hAx(2),'YTick',[0:.2:1])
% set(hLine1, 'Marker', 'd');
% set(hLine1, 'Color', 'g');
% set(hLine1, 'MarkerFaceColor', 'g');
% set(hLine1, 'LineStyle', 'none');
% set(hLine2, 'Marker', 'd');
% set(hLine2, 'Color',  'r');
% set(hLine2, 'MarkerFaceColor', 'r');
% set(hLine2, 'LineStyle', 'none');
% ylabel(hAx(1),'Total GFP area (a.u.)');
% 
% a=findobj(gcf); % get the handles associated with the current figure
% 
% allaxes=findall(a,'Type','axes');
% alllines=findall(a,'Type','line');
% alltext=findall(a,'Type','text');
% 
% set(allaxes,'FontName','Arial','FontWeight','Bold','LineWidth',2,...
%     'FontSize',14);
% set(alllines,'Linewidth',2);
% set(alltext,'FontName','Arial','FontWeight','Bold','FontSize',14)
% 
% 
% set(allaxes,'FontName','Arial','FontWeight','Bold','LineWidth',2,...
%     'FontSize',14);
% set(alllines,'Linewidth',2);
% set(alltext,'FontName','Arial','FontWeight','Bold','FontSize',14)
% print('SpotAvg','-dpdf');
% 
% 
% 
% 
% fileN=fopen('UtotArea.txt','r');
% A=fscanf(fileN,'%g %g %g %g', [4 inf]);
% fclose(fileN);
% fileN=fopen('CtotArea.txt','r');
% B=fscanf(fileN,'%g %g %g %g', [4 inf]);
% fclose(fileN);
% 
% %a=0.68;
% %b=3.0;
% %c=5e-4;
% startP=[0.68 3.0 5e-4];
% eqn='1.0-a/(1+exp(-b*(x-c)))';
% x=A(1,:);
% y=A(3,:)/A(3,1);
% f1=fit(x', y', eqn,'Start',startP)
% xn=logspace(-8,0, 1000);
% figure,semilogx(xn,feval(f1,xn),A(1,:),A(3,:)/A(3,1),'d');
% 
% 
% 



% fileN=fopen('Udist.txt','r');
% A=fscanf(fileN,'%g %g %g %g', [4 inf]);
% fclose(fileN);
% figure, errorbar(A(1,:),A(3,:)/A(3,1),A(4,:)/A(3,1),'o')
% title('Distance to nearest spot');
% xlabel('IPTG concentration');
% ylabel('Distances between red spots (AU)');
% set(gca,'xscale','log');
% 
% a=findobj(gcf); % get the handles associated with the current figure
% 
% allaxes=findall(a,'Type','axes');
% alllines=findall(a,'Type','line');
% alltext=findall(a,'Type','text');
% 
% set(allaxes,'FontName','Arial','FontWeight','Bold','LineWidth',2,...
%     'FontSize',14);
% set(alllines,'Linewidth',2);
% set(alltext,'FontName','Arial','FontWeight','Bold','FontSize',14)
% print('DistBetween','-dpdf');



fileN=fopen('UtotArea.txt','r');
A=fscanf(fileN,'%g %g %g %g', [4 inf]);
fclose(fileN);
fileN=fopen('CtotArea.txt','r');
B=fscanf(fileN,'%g %g %g %g', [4 inf]);
fclose(fileN);
figure;
X1=log10(B(1,:));
Y1=B(4,:)/B(4,end);
X2=log10(A(1,:));
Y2=A(3,:)/A(3,1);
eqn='a*tanh(b*(x-c))+d';
startP=[-0.4 .5 -4 0.6];
f2=fit(X2(2:end)', Y2(2:end)', eqn,'Start',startP)
eqn='a*tanh(b*(x-c))+d';
startP=[-0.4 .5 -4 0.6];
f3=fit(X1(2:end)', Y1(2:end)', eqn,'Start',startP)
X1(1)=-7;
X2(1)=-7;

xn=linspace(-8,-1, 1000);



[hAx,hLine1,hLine2]=plotyy(X1,Y1,X2,Y2,'plot');

hold on;
%plot(xn,-.4*tanh(.5*(xn+8))+0.6);
plot(xn, feval(f2,xn),'r-', xn, feval(f3,xn),'g-');
%plot(xn, feval(f2,xn),'k-', xn, feval(f3,xn),'k-');
hold off;
%title('Red Spot Total Area');
xlabel('log_{10}IPTG (M)');
ylabel(hAx(2),'Total red spot area (a.u.)');
set(hAx(1), 'Xlim', [-7.1 -1.9]);
set(hAx(2), 'Xlim', [-7.1 -1.9]);
set(hAx(2),'Ycolor', [0 0 0]);
set(hAx(1), 'Ycolor', [0 0 0]);
set(hAx(1),'YLim',[0 1.1])
set(hAx(1),'YTick',[0:.2:1])
set(hAx(2),'YLim',[0 1.1])
set(hAx(2),'YTick',[0:.2:1])
set(hLine1, 'Marker', 'd');
set(hLine1, 'Color', 'g');
set(hLine1, 'MarkerFaceColor', 'g');
set(hLine1, 'LineStyle', 'none');
set(hLine2, 'Marker', 'd');
set(hLine2, 'Color',  'r');
set(hLine2, 'MarkerFaceColor', 'r');
set(hLine2, 'LineStyle', 'none');
ylabel(hAx(1),'Spa-averaged GFP(a.u.)');

%axis([-8 -2 0 1.1])

a=findobj(gcf); % get the handles associated with the current figure

allaxes=findall(a,'Type','axes');
alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');

set(allaxes,'FontName','Arial','FontWeight','Bold','LineWidth',3,...
    'FontSize',14);
set(alllines,'Linewidth',3);
set(alltext,'FontName','Arial','FontWeight','Bold','FontSize',14)


print('SpotAvg','-dpdf');
