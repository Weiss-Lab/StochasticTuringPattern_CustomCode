function hyperCube3
%hypercubeSampling
writefile='/home/michael/turingW/direct/All.txt';
res=fopen(writefile, 'a');
sampleNum=5;
D=.1; minD=D*0.5; maxD=D*1.5;
D2=D*21.6; minD2=D2*.5; maxD2=D2*1.5;

gu =1.0; mingu=.247; maxgu=1.19;%.247 to 1.19
gv =1.0; mingv=.247; maxgv=1.19;%.247 to 1.19
gc =1.0; mingc=.35; maxgc=1.04;%.35 to 1.04
giu=1.0; mingiu=.35; maxgiu=1.04;%.35 to 1.04
giv=1.0; mingiv=.35; maxgiv=1.04;%.35 to 1.04

au =3e1; minau=10; maxau=66;%(30) 10  66 (max)
av =3e1; minav=960*.75; maxav=1.25*960; %960
aiu=1e0; minaiu=1; maxaiu=100;%1 to 100
aiv=0.3e0;minaiv=1; maxaiv=100;%1 to 100
ac =1.0e1;minac=1;maxac=100; % 1 to 100

lam_u=1.0;minlam_u=1;maxlam_u=1; % 1
lam_v=1e3;minlam_v=1000;maxlam_v=100000; % 60,000 1000 to 100000

fod_1 =1e3;minfod_1=100; maxfod_1=1000; %100 to 1000
fod_2 =1e5;minfod_2=10; maxfod_2=1000; %10 to 1000
fod_3 =1e3;minfod_3=100; maxfod_3=1000; %1000, 100-1000
fod_4 =1e3;minfod_4=18; maxfod_4=620; %18 - 620
fod_5 =1e5;minfod_5=10; maxfod_5=1000; %10 -1000
fod_6 =1e5;minfod_6=1e5; maxfod_6=1e5; %large

Kd_1 = 1e3;minKd_1=50;maxKd_1=3000; %3000 to 50
Kd_2 = 1e1;minKd_2=10;maxKd_2=100; %10 to 100
Kd_3 = 1e5;minKd_3=2000;maxKd_3=5000; %3000 (should depend on Rv)

Kd_5 = 1e3; minKd_5=100;maxKd_5=1000;%100 to 1000

Kc_3 = 1.5e2; minKc_3=100;maxKc_3=150;% upperbound



Kd_4 = 1e2; minKd_4=20; maxKd_4=70;  %60 probably lower, 20 -70
Kd_6 = 1e-3; minKd_6=.83e-6; maxKd_6=4.3e-6;%.83e-6  - 4.3e-6


Lac  = 1.5e2;minLac=100;maxLac=3000; % 100-3000
IPTG = 1e-6;minIPTG=1e-6;maxIPTG=1e-6; %correct by def 10^-6 to 10^-2

%maxV=[D*1.25 D2*1.25 1.19 1.19 1.04 1.04 1.04 66 66 1.0 aiu aiv ac 100000 1000 1000 1000 620 1000 fod_6 3000 100 3000*1.25 1000 150 70 4.30e-6 1000 IPTG]; %maximium Value in range
%minV=[D*0.75 D2*0.75 .247 .247 0.35 0.35 0.35 10 10 1.0 aiu aiv ac 1000.0 100. 10.0 10.0 18. 10.0 fod_6 50.0 10. 3000*0.75 100. 100 20 0.83e-6 100. IPTG]; %minimium value in range

%maxV=[maxD maxD2 maxgu maxgv maxgc maxgiu maxgiv maxau maxav maxaiu maxaiv maxac maxlam_u maxlam_v maxfod_1 maxfod_2 maxfod_3 maxfod_4 maxfod_5 maxfod_6 maxKd_1 maxKd_2 maxKd_3 maxKd_5 maxKc_3 maxKd_4 maxKd_6 maxLac maxIPTG]; %maximium Value in range
%minV=[minD minD2 mingu mingv mingc mingiu mingiv minau minav minaiu minaiv minac minlam_u minlam_v minfod_1 minfod_2 minfod_3 minfod_4 minfod_5 minfod_6 minKd_1 minKd_2 minKd_3 minKd_5 minKc_3 minKd_4 minKd_6 minLac minIPTG]; %minimium value in range

maxV=[maxD maxD2 maxgu maxgv maxgc maxgiu maxgiv minau minav aiu aiv ac maxlam_u maxlam_v maxfod_1 maxfod_2 maxfod_3 maxfod_4 maxfod_5 maxfod_6 maxKd_1 maxKd_2 maxKd_3 maxKd_5 maxKc_3 maxKd_4 maxKd_6 maxLac maxIPTG]; %maximium Value in range
minV=[minD minD2 mingu mingv mingc mingiu mingiv maxau maxav aiu aiv ac minlam_u minlam_v minfod_1 minfod_2 minfod_3 minfod_4 minfod_5 minfod_6 minKd_1 minKd_2 minKd_3 minKd_5 minKc_3 minKd_4 minKd_6 minLac minIPTG]; %minimium value in range


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
%x0=[162.9157   48.8747    5.4305    1.6292  60]
x0=[6.4149; 1.9245; 0.2138;  0.0641; 31.7450];
%x0=[0.2315 10.2704 0.0053 0.0746 154.1422];
%x0=[0.7804    1.0588    0.0085    0.0195   13.5301];

h = waitbar(0,'Please wait...');
for n =1:sampleNum
    
    options = optimoptions('fsolve','Display','none','MaxIter',50000,'MaxFunEvals',1000000); % Option to display output
    test = @(xs) modelFuncG(xs,values(n,:));
    [x,~,~,~,jacob] = fsolve(test,x0,options); % Call solver
    if any(x<0)
        x
    end
    x
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
        pat=3;
        number3=number3+1;
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
number0/sampleNum
number1/sampleNum
number2/sampleNum
number3/sampleNum

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