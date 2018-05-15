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
fileSpot=fopen('../UtotArea.txt','a');
fileGreen=fopen('../CtotArea.txt','a');
fileID = fopen('4000','r');
strId=sprintf('IPTG=%g D_v/D_u=%g',iptg,difv/difu);
strP = sprintf('nI%gD%g',iptg,difv/difu);
A = fscanf(fileID,'%f %f %f %f %f %f %f %f',[8,inf]);
titlelist={'U','V','Iu','Iv', 'C'};
fclose(fileID);
s=size(A);
B=zeros(64);
y=zeros(64,4);
scaleS=0.05;
%scaleS=1.0;
for jj=[1,5]
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
%bw = imclearborder(bw);
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

%cc=bwconncomp(bw);
%%stats = regionprops(cc, 'Centroid','Area','MinorAxisLength','MajorAxisLength');
%%spotSize = [stats.Area];
sum(sum(bw));

fprintf('spot\n');
fprintf('%g %g %g\n',iptg,difv/difu, sum(sum(bw)));
if jj==1
    fprintf(fileSpot,'%g %g %g\n',iptg,difv/difu, sum(sum(bw)));
else
    fprintf(fileGreen,'%g %g %g\n',iptg,difv/difu, sum(sum(bw)));
end

end
fclose(fileSpot);
fclose(fileGreen);


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
