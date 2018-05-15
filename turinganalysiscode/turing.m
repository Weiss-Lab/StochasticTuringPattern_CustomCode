%filename = 'New Data/crop1.bmp';
%filename = 'test3.png';
filename='testPat.png'
image = imread(filename);
%image=normrnd(0,1,[100,100]);
figure, imshow(image), title 'original';
gray = rgb2gray(image);
gray = mat2gray(gray);
%gray = image;
%gray = im2bw(gray, graythresh(gray)*.2);
figure, imshow(gray), title 'gray scale'
gray = double(gray);
average = mean2(gray);
gray = gray-average;
mean2(gray)
fourier = fftshift(fft2(gray));
power = abs(fourier).^2;
figure, imshow(log(power),[]), title 'power spectrum';
[x, y]=size(power);
areaSize = x*y;
halfx = floor(x/2)+1
halfy = floor(y/2)+1
sizeArray = floor(sqrt(x*x+y*y)/2+3);
Index=zeros(2,sizeArray);
for i=1:x
    for j=1:y
        dist=floor(sqrt((i-halfx)^2+(j-halfy)^2))+1;
        Index(1,dist)=Index(1,dist)+1;
        Index(2,dist)=Index(2,dist)+power(i,j);
    end
end
Radial=double(zeros(1,sizeArray));
for i=[1:sizeArray]
    Radial(i)=Index(2,i)/Index(1,i);
end
figure, plot(Radial(1:halfx),'DisplayName','Radial','YDataSource','Radial');
figure, loglog(Radial(1:halfx),'DisplayName','Radial LogLog', 'Marker','s');
[p,S]=polyfit(log([halfx-50:halfx-20]),log(Radial([halfx-50:halfx-20])),1)
X=1:halfx;
Y=1:halfx;
for i= 1:halfx
    Y(i)=exp(p(2))*i^p(1);
end
legStr=sprintf('%0.1f best fit', p(1));
figure, loglog(1:sizeArray,Radial(:),'s',X,Y), title 'Radial Power Spectrum of RFP', xlabel 'k', ylabel 'Power Spectrum',legend('power spectrum',legStr);
a=findobj(gcf); % get the handles associated with the current figure

allaxes=findall(a,'Type','axes');
alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');

set(allaxes,'FontName','Arial','FontWeight','Bold','LineWidth',2,...
    'FontSize',14);
set(alllines,'Linewidth',2);
set(alltext,'FontName','Arial','FontWeight','Bold','FontSize',14);
%figure, plot(X, Radial(1:halfx)./X.^2), title 'Compensator Plot k^2', xlabel 'k', ylabel 'Power Spectrum / k^2';
bw = im2bw(image, graythresh(image)+.15);
bw = bwareaopen(bw, 50);
figure, imshow(bw), title 'black and white';
cc=bwconncomp(bw);
stats = regionprops(cc, 'Centroid');
Distance = zeros(1, sizeArray);

[x, y]=size(bw);
halfx = floor(x/2)+1;
halfy = floor(y/2)+1;

for i=1:length(stats)-1
    for j=i+1:length(stats)
        a=stats(i).Centroid-stats(j).Centroid;
        diffx = abs(a(1));
        diffy = abs(a(2));
        if diffx > halfx
            diffx = x-diffx;
        end
        if diffy > halfy
            diffy = y - diffy;
        end
        dist=floor(sqrt((diffx)^2+(diffy)^2))+1;
        Distance(dist)=Distance(dist)+1;
    end
end
Distance=Distance*2.0;
g = zeros(1, sizeArray);
prefactor = areaSize/(length(stats)*length(stats));
for i=1:sizeArray
    g(i)=prefactor*Distance(i)/(Index(1,i));
end
figure, plot(1:sizeArray, g), title 'Radial Distribution Function';