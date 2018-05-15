fileID = fopen('660','r');
A = fscanf(fileID,'%f %f %f %f %f %f %f %f',[8,inf]);
fclose(fileID);
s=size(A);
B=zeros(64);
for i=1:s(2)
    ix=round(A(1,i)/.05)+1;
    iy=round(A(2,i)/.05)+1;
    re = A(3,i);
    B(ix,iy)=re;
end
gray=B;
figure, imshow(gray,[]), 'simulation';
average = mean2(gray);
gray = gray-average;
mean2(gray)
fourier = fftshift(fft2(gray));
power = abs(fourier).^2;
figure, imshow(log(power),[]), title 'power spectrum';
[x, y]=size(power);
areaSize = x*y;
halfx = floor(x/2)+1
halfy = floor(y/2)+1;
sizeArray = floor(sqrt(x*x+y*y)/2+3);
Index=zeros(2,sizeArray);
for i=1:x
    for j=1:y
        dist=floor(sqrt((i-halfx)^2+(j-halfy)^2))+1;
        Index(1,dist)=Index(1,dist)+1;
        Index(2,dist)=Index(2,dist)+power(i,j);
    end
end
Radial=zeros(1,sizeArray);
for i=[1:sizeArray]
    Radial(i)=Index(2,i)/Index(1,i);
end
figure, plot(Radial,'DisplayName','Radial','YDataSource','Radial');
figure, loglog(Radial,'DisplayName','Radial LogLog');
p=polyfit(log([halfx-20:halfx-5]),log(Radial([halfx-20:halfx-5])),1)
X=1:sizeArray;
Y=1:sizeArray;
for i= 1:sizeArray
    Y(i)=exp(p(2))*i^p(1);
end
figure, loglog(X,Radial,X,Y), title 'Radial fit';