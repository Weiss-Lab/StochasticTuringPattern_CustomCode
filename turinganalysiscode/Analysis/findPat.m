function findPat()
pat=0;
%x0 = [2197.853961; 1097.673175; 59.000000; 20.000000; 44.000000];  % Make a starting guess at the solution

x0=[6.4149; 1.9245; 0.2138;  0.0641; 31.7450];
%x0=[162.9157   48.8747    5.4305    1.6292  295.2840]
options = optimoptions('fsolve','Display','iter','MaxIter',5000,'MaxFunEvals',100000); % Option to display output
[x,fval,exitflag,output,jacob] = fsolve(@modelFunc,x0,options) % Call solver
jacob
tolerance=0.01;
D=.5;
D2=D*21.6;
%D2=9;
num=10001;
Res = zeros(num,5);
symRes = zeros(num,5);
powerSpec = zeros(num,5);
B0=diag(modelNoise(x));
k=linspace(0.0,30.0,num);
for n=1:num;
    [powerSpec(n,:),Res(n,:), symRes(n,:)]=allValues(k(n), jacob, B0, D, D2, x);
end

%figure(1);
plot(k,Res);
title('Real Parts of Eigenvalues of Jacobian');
ylim([-2 1]);
print('-f1','RealEigenR','-dpng')
figure(2);
plot(k,symRes);
title('Real Parts of Eigenvalues of symmetric Jacobian');
ylim([-20 20]);
print('-f2','RealEigenSymR','-dpng')
figure(3);
plot(k,powerSpec);
title('Power Spectrum');
print('-f3','PowerR','-dpng')
figure(4);
loglog(k,powerSpec(:,1:2));
title('Log Log Power Spectrum');
print('-f4','logPowerR','-dpng')
[val, ind]=max(powerSpec);

k(ind)
val

maxEig = @(ks) -(maxEigenResult(ks, jacob, B0, D, D2, x));
power1 = @(ks) -(powerSpectrum1(ks, jacob, B0, D, D2, x));
power2 = @(ks) -(powerSpectrum2(ks, jacob, B0, D, D2, x));
%power3 = @(ks) -(powerSpectrum3(ks, jacob, B0, D, D2, x));
%power4 = @(ks) -(powerSpectrum4(ks, jacob, B0, D, D2, x));
%power5 = @(ks) -(powerSpectrum5(ks, jacob, B0, D, D2, x));

[kminE,eMax]=fminbnd(maxEig,0,100,optimset('TolX',tolerance))
eMax=-eMax
[kmin1,f1]=fminbnd(power1,0,100,optimset('TolX',tolerance))
[kmin2,f2]=fminbnd(power2,0,100,optimset('TolX',tolerance))
%[kmin3,f3]=fminbnd(power3,0,100,optimset('TolX',tolerance))
%[kmin4,f4]=fminbnd(power4,0,100,optimset('TolX',tolerance))
%[kmin5,f5]=fminbnd(power5,0,100,optimset('TolX',tolerance))

[powerSpec0,reig0, seig0]=allValues(0, jacob, B0, D, D2, x)
[powerSpec1,reig1, seig1]=allValues(kmin1, jacob, B0, D, D2, x)
[powerSpec2,reig2, seig2]=allValues(kmin2, jacob, B0, D, D2, x)
%[powerSpec3,reig3, seig3]=allValues(kmin3, jacob, B0, D, D2, x)
%[powerSpec4,reig4, seig4]=allValues(kmin4, jacob, B0, D, D2, x)
%[powerSpec5,reig5, seig5]=allValues(kmin5, jacob, B0, D, D2, x)

%determine type of pattern
%0 unstable (k=0 eigenvalues positive)
%1 no pattern (no peaks in powerspectrum for species 1 or 2)
%2 quasi turing pattern; peak in power spectrum at k!=0 and all eigenvalues negative
%3 turing pattern; peak in power spectrum and some eigenvalues are positive
if(any(reig0>0))
    pat=0;
elseif((kmin1<2*tolerance && kmin2<2*tolerance) || (abs(powerSpec1(1)-powerSpec0(1))<tolerance && abs(powerSpec2(2)-powerSpec0(2))<tolerance))
    pat=1;
elseif(any(reig1>0) || any(reig2>0))
    pat=3;
else
    pat=2;    
end
pat


function [powerS, nEig, symEig] = allValues(k, jacob, B0, D, D2, x)
    modJacob=jacob-diag([D*k^2 D2*k^2 0 0 0]);
    B=B0+diag([D*x(1)*k^2 D2*x(2)*k^2 0 0 0]);
    iD=inv(modJacob);
    power = iD*B*iD';
    powerS=diag(power);
    sym = (modJacob+modJacob')/2.0;
    ejac = real(eig(modJacob));
    nEig=ejac';
    ejac = real(eig(sym));
    symEig=ejac';
end

function [el1] = powerSpectrum1(k, jacob, B0, D, D2, x)
    modJacob=jacob-diag([D*k^2 D2*k^2 0 0 0]);
    B=B0+diag([D*x(1)*k^2 D2*x(2)*k^2 0 0 0]);
    iD=inv(modJacob);
    power = iD*B*iD';
    powerS=diag(power);
    el1=powerS(1);
end

function [el1] = powerSpectrum2(k, jacob, B0, D, D2, x)
    modJacob=jacob-diag([D*k^2 D2*k^2 0 0 0]);
    B=B0+diag([D*x(1)*k^2 D2*x(2)*k^2 0 0 0]);
    iD=inv(modJacob);
    power = iD*B*iD';
    powerS=diag(power);
    el1=powerS(2);
end
function [el1] = powerSpectrum3(k, jacob, B0, D, D2, x)
    modJacob=jacob-diag([D*k^2 D2*k^2 0 0 0]);
    B=B0+diag([D*x(1)*k^2 D2*x(2)*k^2 0 0 0]);
    iD=inv(modJacob);
    power = iD*B*iD';
    powerS=diag(power);
    el1=powerS(3);
end


function [el1] = powerSpectrum4(k, jacob, B0, D, D2, x)
    modJacob=jacob-diag([D*k^2 D2*k^2 0 0 0]);
    B=B0+diag([D*x(1)*k^2 D2*x(2)*k^2 0 0 0]);
    iD=inv(modJacob);
    power = iD*B*iD';
    powerS=diag(power);
    el1=powerS(4);
end

function [el1] = powerSpectrum5(k, jacob, B0, D, D2, x)
    modJacob=jacob-diag([D*k^2 D2*k^2 0 0 0]);
    B=B0+diag([D*x(1)*k^2 D2*x(2)*k^2 0 0 0]);
    iD=inv(modJacob);
    power = iD*B*iD';
    powerS=diag(power);
    el1=powerS(5);
end

function [nEig, symEig] = eigenResults(k, jacob, B0, D, D2, x)
    modJacob=jacob-diag([D*k^2 D2*k^2 0 0 0]);
    sym = (modJacob+modJacob')/2.0;
    ejac = real(eig(modJacob));
    nEig=ejac';
    ejac = real(eig(sym));
    symEig=ejac';
end

function [mEig] = maxEigenResult(k, jacob, B0, D, D2, x)
    modJacob=jacob-diag([D*k^2 D2*k^2 0 0 0]);
    sym = (modJacob+modJacob')/2.0;
    ejac = real(eig(modJacob));
    mEig=max(ejac);
end

% 
% temp = @(ks,t) norm([1 1 1 1 1]*expm(t*(jacob-diag([D*ks^2 D2*ks^2 0 0 0]))));
% kMax = @(t) fminbnd(@(ks) -temp(ks,t),0,100);
% grow = @(t) temp(kMax(t),t);

% g2Pow=zeros(num,2);
% gPow=zeros(num, 2);
% t=linspace(0,40,num);
% kres=zeros(num,1);
% gres=zeros(num,1);
% for n=1:num
%     gPow(n,1)=temp(k(ind(1)),t(n));
%     gPow(n,2)=temp(k(ind(2)),t(n));
%     g2Pow(n,1)=norm([1 0 0 0 0]*expm(t(n)*(jacob-diag([D*k(ind(1))^2 D2*k(ind(1))^2 0 0 0]))))
%     g2Pow(n,2)=norm([0 1 0 0 0]*expm(t(n)*(jacob-diag([D*k(ind(2))^2 D2*k(ind(2))^2 0 0 0]))))
%     kres(n) = kMax(t(n));
%     gres(n) = grow(t(n));
% end
% figure, plot(t,gPow), title('Amplification at peak of power spectrum');
% 
% figure, plot(t,g2Pow), title('Amplification per channel');
% 
% figure;
% plot(t,gres), title('Amplification vs time');
% 
% figure;
% plot(t,kres), title('kmax vs time');
end