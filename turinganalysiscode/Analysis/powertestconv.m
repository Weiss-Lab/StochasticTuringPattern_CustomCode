mask = 1./(1+exp(-20*(x-pi/10)))*1./(1+exp(20*(x-(2*pi-pi/10))));

fileID = fopen('predatorPrey2.txt','r');
sampleN=10000;
deltaT=.1;

x=1:sampleN/2;
grow=20;
x0=10;
mask= 1./(1+exp(-grow(x-x0

spower = zeros(1,sampleN/2);

for n=1:1000

A = fscanf(fileID,'%f %f %f %f',[4,sampleN]);

s=size(A);
V=A(2,:);
m=mean(V);
V=V-m;
num=s(2);
fourier = fft(V)*2*pi/num/deltaT;

power = abs(fourier).^2;
power = power(1:num/2);
spower=power+spower;
end

fclose(fileID);
power=spower/1000.0;

t=(0:num/2-1);
t=(2*pi/num/deltaT)*t;
figure, plot(t,power,'.'), title 'power';
figure, loglog(t,power,'.'), title 'power spectrum';
start=40;
fin=num/4;
halfx=num/2;
%temp=log(t([start:halfx-fin]))
%temp2=log(power([start:halfx-fin]))
[p,S]=polyfit(log(t([start:halfx-fin])),log(power([start:halfx-fin])),1)
Ys=exp(p(2)).*t.^p(1);
graphTitle = sprintf('fit %g', p(1));
figure, loglog(t,power,t,Ys), title(graphTitle);