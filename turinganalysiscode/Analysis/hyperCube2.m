function hyperCube2
%hypercubeSampling
writefile='/home/michael/turingW/direct/varyAll2.txt';
res=fopen(writefile, 'a');
sampleNum=500;
D=.1;
D2=D*21.6;

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

%maxV=[D*1.5 D2*1.5 gu*1.5 gv*1.5 gc*1.5 giu*1.5 giv*1.5 au*1.5 av*1.5 aiu*1.5 aiv*1.5 ac*1.5 lam_u*1.5 lam_v*1.5 fod_1*1.5 fod_2*1.5 fod_3*1.5 fod_4*1.5 fod_5*1.5 fod_6*1.5 Kd_1*1.5 Kd_2*1.5 Kd_3*1.5 Kd_5*1.5 Kc_3*1.5 Kd_4*1.5 Kd_6*1.5 Lac*1.5 IPTG*1.5]; %maximium Value in range
%minV=[D*0.5 D2*0.5 gu*0.5 gv*0.5 gc*0.5 giu*0.5 giv*0.5 au*0.5 av*0.5 aiu*0.5 aiv*0.5 ac*0.5 lam_u*0.5 lam_v*0.5 fod_1*0.5 fod_2*0.5 fod_3*0.5 fod_4*0.5 fod_5*0.5 fod_6*0.5 Kd_1*0.5 Kd_2*0.5 Kd_3*0.5 Kd_5*0.5 Kc_3*0.5 Kd_4*0.5 Kd_6*0.5 Lac*0.5 IPTG*0.5]; %minimium value in range

maxV=[D D2 gu gv gc giu giv au*1.5 av*1.5 aiu*1.5 aiv*1.5 ac*1.5 lam_u*1.5 lam_v*1.5 fod_1*1.5 fod_2*1.5 fod_3*1.5 fod_4*1.5 fod_5*1.5 fod_6*1.5 Kd_1*1.5 Kd_2*1.5 Kd_3*1.5 Kd_5*1.5 Kc_3*1.5 Kd_4*1.5 Kd_6*1.5 Lac*1.5 IPTG*1.5]; %maximium Value in range
minV=[D D2 gu gv gc giu giv au*0.5 av*0.5 aiu*0.5 aiv*0.5 ac*0.5 lam_u*0.5 lam_v*0.5 fod_1*0.5 fod_2*0.5 fod_3*0.5 fod_4*0.5 fod_5*0.5 fod_6*0.5 Kd_1*0.5 Kd_2*0.5 Kd_3*0.5 Kd_5*0.5 Kc_3*0.5 Kd_4*0.5 Kd_6*0.5 Lac*0.5 IPTG*0.5]; %minimium value in range


%maxV=[D D2*1.5 gu gv gc giu giv au av aiu aiv ac lam_u lam_v fod_1 fod_2 fod_3 fod_4 fod_5 fod_6 Kd_1 Kd_2 Kd_3 Kd_5 Kc_3 Kd_4 Kd_6 Lac IPTG*]; %maximium Value in range
%minV=[D D2*0.5 gu gv gc giu giv au av aiu aiv ac lam_u lam_v fod_1 fod_2 fod_3 fod_4 fod_5 fod_6 Kd_1 Kd_2 Kd_3 Kd_5 Kc_3 Kd_4 Kd_6 Lac IPTG]; %minimium value in range


numEl=numel(maxV);
delta=(maxV-minV)/sampleNum;
values=rand(sampleNum,numEl)*diag(delta);

for n=1:numEl
    values(:,n)=values(:,n)+linspace(0,sampleNum-1,sampleNum)'*delta(n)+repmat(minV(n),sampleNum,1);
    values(:,n)=values(randperm(sampleNum),n);
end
%values(1,:)


Diffusion=zeros(sampleNum,1);
patternType=zeros(sampleNum,1);
waveNumber=zeros(sampleNum,1);
waveNumber2=zeros(sampleNum,1);
waveNumberE=zeros(sampleNum,1);
maxPower=zeros(sampleNum,1);
number0=0;
number1=0;
number2=0;
number3=0;
number4=0;
x0=[6.4149; 1.9245; 0.2138;  0.0641; 31.7450];

h = waitbar(0,'Please wait...');
for n =1:sampleNum
    
    options = optimoptions('fsolve','Display','none','MaxIter',50000,'MaxFunEvals',1000000); % Option to display output
    test = @(xs) modelFuncG(xs,values(n,:));
    [x,~,~,~,jacob] = fsolve(test,x0,options); % Call solver
    if any(x<0)
        x
    end
    x0=x;
    tolerance=0.01;
    D=values(n,1);
    D2=values(n,2);
    
    B0=diag(modelNoiseG(x,values(n,:)));
    
    maxEig = @(ks) -(maxEigenResult(ks, jacob, B0, D, D2, x));
    power1 = @(ks) -(powerSpectrum1(ks, jacob, B0, D, D2, x));
    power2 = @(ks) -(powerSpectrum2(ks, jacob, B0, D, D2, x));
    %power3 = @(ks) -(powerSpectrum3(ks, jacob, B0, D, D2, x));
    %power4 = @(ks) -(powerSpectrum4(ks, jacob, B0, D, D2, x));
    %power5 = @(ks) -(powerSpectrum5(ks, jacob, B0, D, D2, x));
    
    [kminE,eMax]=fminbnd(maxEig,0,100,optimset('TolX',tolerance));
    eMax=-eMax;
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
        number0=number0+1;
    elseif((kmin1<10*tolerance && kmin2<10*tolerance) || (abs(powerSpec1(1)-powerSpec0(1))<tolerance && abs(powerSpec2(2)-powerSpec0(2))<tolerance))
        pat=1;
        number1=number1+1;
    elseif(any(reig1>0) || any(reig2>0)|| eMax>0)
        if(any(reig1>0) || any(reig2>0))
            pat=3;
            number3=number3+1;
        else
            pat=4;
            number4=number4+1;
        end
    else
        pat=2;
        number2=number2+1;
    end
    
    Diffusion(n)=values(n,2);
    patternType(n)=pat;
    waveNumber(n)=kmin1;
    waveNumber2(n)=kmin2;
    waveNumberE(n)=kminE;
    maxPower(n)=eMax;
    
    
    fprintf(res, '%f ', values(n,:));
    fprintf(res, '%f ', x);
    fprintf(res, '%f ', powerSpec0);
    fprintf(res, '%f ', reig0);
    fprintf(res, '%f ', seig0);
    
    fprintf(res, '%f ', powerSpec1);
    fprintf(res, '%f ', reig1);
    fprintf(res, '%f ', seig1);
        
    fprintf(res, '%f ', powerSpec2);
    fprintf(res, '%f ', reig2);
    fprintf(res, '%f ', seig2);
    
    fprintf(res, '%f %f %f %f %f\n', kmin1, kminE, kmin2, eMax, pat);
    waitbar(n / sampleNum)
    
end

close(h) 
fclose(res);
figure, plot(Diffusion, patternType, '.'), title('Pattern Type');
figure, plot(Diffusion, waveNumber,'.', Diffusion, waveNumber2,'.', Diffusion, waveNumberE,'.'), title('Wavenumber vs D2');
figure, plot(Diffusion, maxPower,'.'), title('Max Eigenvalue for chemical 1');
number0
number1
number2
number3
number4
number0/sampleNum
number1/sampleNum
number2/sampleNum
number3/sampleNum
number4/sampleNum

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

function [mEig] = maxEigenResult(k, jacob, B0, D, D2, x)
    modJacob=jacob-diag([D*k^2 D2*k^2 0 0 0]);
    sym = (modJacob+modJacob')/2.0;
    ejac = real(eig(modJacob));
    mEig=max(ejac);
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