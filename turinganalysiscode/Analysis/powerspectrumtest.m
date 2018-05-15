num=2^9;
x=linspace(0,2*pi,num+1);
x=x(1:num);

mask = (1./(1+exp(-20*(x-pi/10)))*1./(1+exp(20*(x-(2*pi-pi/10))))).^2;
figure, plot(x,mask), title 'mask';

y=sin(10*x)+normrnd(0,0,1,num);
comp1=fft(y);
comp2=fft(mask);
comb=conv(comp1,comp2)/num;
rpow=abs(comb).^2;
rpow=rpow(1:num/2);

y=y.*mask;
figure, plot(x,y), title 'sin wave';
fourier = fft(y);

%fourier2=fft(mask);
%[result,res]=deconv(fourier,fourier2);
%fourier=result;
%[~, num]=size(fourier);
power = abs(fourier).^2;
power = power(1:num/2);
figure, plot(0:num/2-1, power), title 'power';
figure, loglog(power), title 'power spectrum';
start=50;
fin=num/4;
halfx=num/2;
[p,S]=polyfit(log([start:halfx-fin]),log(power([start:halfx-fin])),1)
Xs=0:halfx-1;
Ys=exp(p(2)).*Xs.^p(1);
graphTitle = sprintf('fit %g', p(1));
figure, loglog(Xs,power,Xs,Ys,Xs,rpow), title(graphTitle);

