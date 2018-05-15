function hyperCube
%hypercubeSampling
sampleNum=5;
D=.1;
D2=D*100;

gu =1.0;
gv =1.0;
gc =1.0;
giu=1.0;
giv=1.0;

au =3e1;
av =3e1;
aiu=1e0;
aiv=0.3e0;
ac =1.0e1;

lam_u=1.0;
lam_v=1e3;

fod_1 =1e3;
fod_2 =1e5;
fod_3 =1e3;
fod_4 =1e3;
fod_5 =1e5;
fod_6 =1e5;

Kd_1 = 1e3;
Kd_2 = 1e1;
Kd_3 = 1e5;

Kd_5 = 1e3;

Kc_3 = 1.5e2;



Kd_4 = 1e2;
Kd_6 = 1e-3;


Lac  = 1.5e2;
IPTG = 1e-6;


maxV=[D*10 D2*10 gu*10 gv*10 gc*10 giu*10 giv*10 au*10 av*10 aiu*10 aiv*10 ac*10 lam_u*10 lam_v*10 fod_1*10 fod_2*10 fod_3*10 fod_4*10 fod_5*10 fod_6*10 Kd_1*10 Kd_2*10 Kd_3*10 Kd_5*10 Kc_3*10 Kd_4*10 Kd_6*10 Lac*10 IPTG*10]; %maximium Value in range
minV=[D/10 D2/10 gu/10 gv/10 gc/10 giu/10 giv/10 au/10 av/10 aiu/10 aiv/10 ac/10 lam_u/10 lam_v/10 fod_1/10 fod_2/10 fod_3/10 fod_4/10 fod_5/10 fod_6/10 Kd_1/10 Kd_2/10 Kd_3/10 Kd_5/10 Kc_3/10 Kd_4/10 Kd_6/10 Lac/10 IPTG/10]; %minimium value in range
numEl=numel(maxV);
delta=(maxV-minV)/sampleNum;
values=rand(sampleNum,numEl)*diag(delta);

for n=1:numEl
    values(:,n)=values(:,n)+linspace(0,sampleNum-1,sampleNum)'*delta(n)+repmat(minV(n),sampleNum,1);
    values(:,n)=values(randperm(sampleNum),n);
end
%values(1,:)



for n =1:sampleNum
    pat=0;
    x0=[6.4149; 1.9245; 0.2138;  0.0641; 31.7450];
    
    options = optimoptions('fsolve','Display','iter','MaxIter',5000,'MaxFunEvals',100000); % Option to display output
    test = @(x) modelFuncG(x,values(n,:));
    [x,~,~,~,jacob] = fsolve(test,x0,options) % Call solver
    tolerance=0.01;
    D=.1;
    %D2=D*21.6;
    D2=D*100;
    
    B0=diag(modelNoiseG(x,values(n,:)));
    
    
    power1 = @(ks) -(powerSpectrum1(ks, jacob, B0, D, D2, x));
    power2 = @(ks) -(powerSpectrum2(ks, jacob, B0, D, D2, x));
    %power3 = @(ks) -(powerSpectrum3(ks, jacob, B0, D, D2, x));
    %power4 = @(ks) -(powerSpectrum4(ks, jacob, B0, D, D2, x));
    %power5 = @(ks) -(powerSpectrum5(ks, jacob, B0, D, D2, x));
    
    [kmin1,f1]=fminbnd(power1,0,100,optimset('TolX',tolerance));
    [kmin2,f2]=fminbnd(power2,0,100,optimset('TolX',tolerance));
    %[kmin3,f3]=fminbnd(power3,0,100,optimset('TolX',tolerance))
    %[kmin4,f4]=fminbnd(power4,0,100,optimset('TolX',tolerance))
    %[kmin5,f5]=fminbnd(power5,0,100,optimset('TolX',tolerance))
    
    [powerSpec0,reig0, seig0]=allValues(0, jacob, B0, D, D2, x);
    [powerSpec1,reig1, seig1]=allValues(kmin1, jacob, B0, D, D2, x);
    [powerSpec2,reig2, seig2]=allValues(kmin2, jacob, B0, D, D2, x);
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
end


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
% function [el1] = powerSpectrum3(k, jacob, B0, D, D2, x)
%     modJacob=jacob-diag([D*k^2 D2*k^2 0 0 0]);
%     B=B0+diag([D*x(1)*k^2 D2*x(2)*k^2 0 0 0]);
%     iD=inv(modJacob);
%     power = iD*B*iD';
%     powerS=diag(power);
%     el1=powerS(3);
% end
% 
% 
% function [el1] = powerSpectrum4(k, jacob, B0, D, D2, x)
%     modJacob=jacob-diag([D*k^2 D2*k^2 0 0 0]);
%     B=B0+diag([D*x(1)*k^2 D2*x(2)*k^2 0 0 0]);
%     iD=inv(modJacob);
%     power = iD*B*iD';
%     powerS=diag(power);
%     el1=powerS(4);
% end
% 
% function [el1] = powerSpectrum5(k, jacob, B0, D, D2, x)
%     modJacob=jacob-diag([D*k^2 D2*k^2 0 0 0]);
%     B=B0+diag([D*x(1)*k^2 D2*x(2)*k^2 0 0 0]);
%     iD=inv(modJacob);
%     power = iD*B*iD';
%     powerS=diag(power);
%     el1=powerS(5);
% end

% function [nEig, symEig] = eigenResults(k, jacob, B0, D, D2, x)
%     modJacob=jacob-diag([D*k^2 D2*k^2 0 0 0]);
%     sym = (modJacob+modJacob')/2.0;
%     ejac = real(eig(modJacob));
%     nEig=ejac';
%     ejac = real(eig(sym));
%     symEig=ejac';
% end
end