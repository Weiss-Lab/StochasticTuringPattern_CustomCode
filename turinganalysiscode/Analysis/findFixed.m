%x0 = [2197.853961; 1097.673175; 59.000000; 20.000000; 44.000000];  % Make a starting guess at the solution

x0=[6.4149; 1.9245; 0.2138;  0.0641; 31.7450];
%x0=[162.9157   48.8747    5.4305    1.6292  295.2840]
options = optimoptions('fsolve','Display','iter','MaxIter',5000,'MaxFunEvals',100000); % Option to display output
[x,fval,exitflag,output,jacob] = fsolve(@modelFunc,x0,options) % Call solver

D=.1;
D2=D*21.6;
%D2=D*100.0;
num=201;
Res = zeros(num,5);
symRes = zeros(num,5);
powerSpec = zeros(num,5);
B0=diag(modelNoise(x));
k=linspace(0.0,30.0,num);
for n=1:num;
    modJacob=jacob-diag([D*k(n)^2 D2*k(n)^2 0 0 0]);
    B=B0+diag([D*x(1)*k(n)^2 D2*x(2)*k(n)^2 0 0 0]);
    iD=inv(modJacob);
    power = iD*B*iD';
    powerSpec(n,:)=diag(power);
    sym = (modJacob+modJacob')/2.0;
    ejac = real(eig(modJacob));
    Res(n,:)=ejac';
    ejac = real(eig(sym));
    symRes(n,:)=ejac';
end
figure;
plot(k,Res);
title('Real Parts of Eigenvalues of Jacobian');
ylim([-2 1]);
figure;
plot(k,symRes);
title('Real Parts of Eigenvalues of symmetric Jacobian');
ylim([-20 20]);
figure;
plot(k,powerSpec);
title('Power Spectrum');
figure;
loglog(k,powerSpec);
title('Log Log Power Spectrum');
[val, ind]=max(powerSpec);

k(ind)
powerSpec(ind)

temp = @(ks,t) norm([1 1 1 1 1]*expm(t*(jacob-diag([D*ks^2 D2*ks^2 0 0 0]))));
kMax = @(t) fminbnd(@(ks) -temp(ks,t),0,100);
grow = @(t) temp(kMax(t),t);

g2Pow=zeros(num,2);
gPow=zeros(num, 2);
t=linspace(0,40,num);
kres=zeros(num,1);
gres=zeros(num,1);
for n=1:num
    gPow(n,1)=temp(k(ind(1)),t(n));
    gPow(n,2)=temp(k(ind(2)),t(n));
    g2Pow(n,1)=norm([1 0 0 0 0]*expm(t(n)*(jacob-diag([D*k(ind(1))^2 D2*k(ind(1))^2 0 0 0]))))
    g2Pow(n,2)=norm([0 1 0 0 0]*expm(t(n)*(jacob-diag([D*k(ind(2))^2 D2*k(ind(2))^2 0 0 0]))))
    kres(n) = kMax(t(n));
    gres(n) = grow(t(n));
end
figure, plot(t,gPow), title('Amplification at peak of power spectrum');

figure, plot(t,g2Pow), title('Amplification per channel');

figure;
plot(t,gres), title('Amplification vs time');

figure;
plot(t,kres), title('kmax vs time');