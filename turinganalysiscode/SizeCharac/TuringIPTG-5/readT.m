fileI = fopen('PARAMS_RUN');
for n=1:8
tline=fgetl(fileI);
end
[difc]=sscanf(tline,'dif_u=%g, dif_v=%g, dif_iu=0, dif_iv=0, dif_c=0, dif_n=0');
difu=difc(1);
difv=difc(2);
for n=1:6
tline=fgetl(fileI);
end
iptg=sscanf(tline,'Ncri=1, CHI=0, Lac=150, IPTG=%g');
fclose(fileI);
fileSpot={fopen('../Uspot.txt','a'),fopen('../Vspot.txt','a')};
fileDist={fopen('../Udist.txt','a'),fopen('../Vdist.txt','a')};
fileID = fopen('4000','r');
strId=sprintf('IPTG=%g D_v/D_u=%g',iptg,difv/difu);
strP = sprintf('I%gD%g',iptg,difv/difu);
A = fscanf(fileID,'%f %f %f %f %f %f %f %f',[8,inf]);
titlelist={'U','V','Iu','Iv', 'C'};
fclose(fileID);
s=size(A);
B=zeros(64);
y=zeros(64,4);
scaleS=0.05;
%scaleS=1.0;
for jj=1:2
for i=1:s(2)
    ix=round(A(1,i)/.05)+1;
    iy=round(A(2,i)/.05)+1;
    re = A(2+jj,i);
    B(ix,iy)=re;
end
im=B;
%subplot(3,2,j);
figure, imagesc(im), suptitle(['Concentration of ' titlelist{jj}]);
set(gca, 'XTickLabel', '');
set(gca, 'YTickLabel', '');
colormap(gray);
colorbar;
a=findobj(gcf); % get the handles associated with the current figure

allaxes=findall(a,'Type','axes');
alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');

set(allaxes,'FontName','Arial','FontWeight','Bold','LineWidth',2,...
    'FontSize',14);
set(alllines,'Linewidth',2);
set(alltext,'FontName','Arial','FontWeight','Bold','FontSize',14)
title(strId,'fontsize',10);
print(['Density' strP titlelist{jj}],'-dpdf');
gr=mat2gray(im);
bw=im2bw(gr,graythresh(gr)*1.2);
bw = imclearborder(bw);
figure, imagesc(bw);
suptitle(['Binary Mask of ' titlelist{jj}]);
colormap(gray);
set(gca, 'XTickLabel', '');
set(gca, 'YTickLabel', '');
a=findobj(gcf); % get the handles associated with the current figure

allaxes=findall(a,'Type','axes');
alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');

set(allaxes,'FontName','Arial','FontWeight','Bold','LineWidth',2,...
    'FontSize',14);
set(alllines,'Linewidth',2);
set(alltext,'FontName','Arial','FontWeight','Bold','FontSize',14)
title(strId,'fontsize',10);
print(['BWEdgeRemoved' strP titlelist{jj}],'-dpdf');

[x, y]=size(bw);
halfx = x/2;
halfy = y/2;

cc=bwconncomp(bw);
stats = regionprops(cc, 'Centroid','Area','MinorAxisLength','MajorAxisLength');
spotSize = [stats.Area];
minor=[stats.MinorAxisLength];
major=[stats.MajorAxisLength];
sizeArray=length(stats);
Distance = zeros(1, sizeArray);
minDistance = zeros(1, sizeArray);
for i=1:length(stats)
    for j=1:length(stats)
        a=stats(i).Centroid-stats(j).Centroid;
        diffx = abs(a(1));
        diffy = abs(a(2));
%         if diffx > halfx
%             diffx = x-diffx;
%         end
%         if diffy > halfy
%             diffy = y - diffy;
%         end
        dist=sqrt((diffx)^2+(diffy)^2);
        if(dist==0)
            dist=1000;
        end
        Distance(j)=dist;
    end
    minDistance(i)=min(Distance);
end
sminDistance =scaleS*minDistance;
spaces=(0:19)*scaleS+scaleS/2.0;
figure, hist(sminDistance,spaces), suptitle(['Histogram of minimium distance to next center of chemical ',titlelist{jj}]), xlabel 'distance to next center (grid)', ylabel 'frequency';
a=findobj(gcf); % get the handles associated with the current figure

allaxes=findall(a,'Type','axes');
alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');

set(allaxes,'FontName','Arial','FontWeight','Bold','LineWidth',2,...
    'FontSize',14);
set(alllines,'Linewidth',2);
set(alltext,'FontName','Arial','FontWeight','Bold','FontSize',14);
title(strId,'fontsize',10);
xlim([0 inf]);
print(['histNear' strP titlelist{jj}],'-dpdf');
%figure, histfit(minDistance), title 'best fit normal distribution';
%pd=fitdist(minDistance(:),'normal')

fprintf('distance\n');
fprintf('%g %g %g %g\n',iptg,difv/difu,mean(sminDistance),std(sminDistance));
fprintf(fileDist{jj},'%g %g %g %g\n',iptg,difv/difu,mean(sminDistance),std(sminDistance));
%mean(sminDistance)
%std(sminDistance)
%median(sminDistance)
%paramEstsMinima = evfit(-sminDistance)
y = linspace(15,80,65);
%figure, hist(minDistance,15), title 'best fit Gumbel';
%p = evpdf(-y,paramEstsMinima(1),paramEstsMinima(2));
%line(y,3*length(minDistance)*p,'color','r');
sizeSpot=sqrt(spotSize/pi);
%figure, histfit(sizeSpot), title 'best fit normal size of spot'
%pd2=fitdist(sizeSpot','normal')


scaleA=scaleS*scaleS;
scaledSpotSize=spotSize*scaleA;
spaces=(0:19)*5.0*scaleA+5.0*scaleA/2.0;
fprintf('spot\n');
fprintf('%g %g %g %g\n',iptg,difv/difu, mean(scaledSpotSize),std(scaledSpotSize));
fprintf(fileSpot{jj},'%g %g %g %g\n',iptg,difv/difu, mean(scaledSpotSize),std(scaledSpotSize));

figure, hist(scaledSpotSize,spaces), suptitle(['Histogram of the area of chemical ', titlelist{jj}]);



ylabel 'Frequency';
xlabel 'Area of Spots (grid x grid)';
a=findobj(gcf); % get the handles associated with the current figure

allaxes=findall(a,'Type','axes');
alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');

set(allaxes,'FontName','Arial','FontWeight','Bold','LineWidth',2,...
    'FontSize',14);
set(alllines,'Linewidth',2);
set(alltext,'FontName','Arial','FontWeight','Bold','FontSize',14)
title(strId,'fontsize',10);
xlim([0 inf]);
print(['histArea' strP titlelist{jj}],'-dpdf');
% figure, hist(sizeSpot,10), title(['histogram of sqrt area ' titlelist{jj}]);
% ylabel 'frequency'
% xlim([0 15]);
% print(['histSqrtArea' titlelist{jj}],'-dpng');
% figure, hist(minor,10), title(['minor axis length ' titlelist{jj}]);
% ylabel 'frequency'
% xlim([0 15]);
% print(['histMinor' titlelist{jj}],'-dpng');
% figure, hist(major,10), title(['major axis length ' titlelist{jj}]) ;
% ylabel 'frequency'
% xlim([0 15]);
% print(['histMajor' titlelist{jj}],'-dpng');


%ypos=45;
%y(:,j)= mat2gray(B(ypos,:));
%y(:,j)= B(ypos,:);
%line([1 64], [ypos ypos],'Color','r');

fclose(fileSpot{jj});
fclose(fileDist{jj});

end



% subplot(3,2,[5 6]);
% plot(linspace(0,3.15,64),y);
% title('Normalized slice of chemicals');
% legend(titlelist{1:4})


% average = mean2(im);
% im = im-average;
% 
% %gray = imdilate(gray, strel('disk',1));
% %gray = im2bw(gray, graythresh(gray));
% 
% %mean2(gray)
% fourier = fftshift(fft2(im));
% power = abs(fourier).^2;
% figure, imshow(log(power),[]), title 'power spectrum';
% [x, y]=size(power);
% areaSize = x*y;
% halfx = floor(x/2)+1
% halfy = floor(y/2)+1;
% sizeArray = floor(sqrt(x*x+y*y)/2+3);
% Index=zeros(2,sizeArray);
% for i=1:x
%     for j=1:y
%         dist=floor(sqrt((i-halfx)^2+(j-halfy)^2))+1;
%         Index(1,dist)=Index(1,dist)+1;
%         Index(2,dist)=Index(2,dist)+power(i,j);
%     end
% end
% Radial=zeros(1,sizeArray);
% for i=[1:sizeArray]
%     Radial(i)=Index(2,i)/Index(1,i);
% end
% figure, plot(Radial,'DisplayName','Radial','YDataSource','Radial');
% figure, loglog(Radial,'DisplayName','Radial LogLog');
% p=polyfit(log([halfx-14:halfx-5]),log(Radial([halfx-14:halfx-5])),1)
% X=1:sizeArray;
% Y=1:sizeArray;
% for i= 1:sizeArray
%     Y(i)=exp(p(2))*i^p(1);
% end
% legFit=sprintf('best fit %.2g',p(1));
% figure, loglog(X,Radial,'s',X,Y),title 'Simulated Radial Power Spectrum', xlabel 'k', ylabel 'Power Spectrum',legend('power spectrum',legFit);
% %axis([1 1000 1E7 1E15])
% a=findobj(gcf); % get the handles associated with the current figure
% 
% allaxes=findall(a,'Type','axes');
% alllines=findall(a,'Type','line');
% alltext=findall(a,'Type','text');
% 
% set(allaxes,'FontName','Arial','FontWeight','Bold','LineWidth',2,...
%     'FontSize',14);
% set(alllines,'Linewidth',2);
% set(alltext,'FontName','Arial','FontWeight','Bold','FontSize',14);
