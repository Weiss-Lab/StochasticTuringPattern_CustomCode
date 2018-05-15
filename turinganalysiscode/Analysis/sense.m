num =201;
kres=zeros(num,1);
gres=zeros(num,1);
au = logspace(-2, 2, num);
x0=[6.4149; 1.9245; 0.2138;  0.0641; 31.7450];
for n = 1:num
    options = optimoptions('fsolve','Display','iter'); % Option to display output
    test = @(x) modelFuncT(x,au(n));
    [x,fval,exitflag,output,jacob] = fsolve(test,x0,options); % Call solver
    x0=x;
    D=.1;
    D2=D*21.6;
    %D2=D*100.0;

    func = @(k) max(real(eig(jacob-diag([D*k^2 D2*k^2 0 0 0]))));
    maxk = fminbnd(@(ks) -func(ks),0,30);
    kres(n)=maxk;
    gres(n)=func(maxk);
end
figure(1);
semilogx(au, gres);

figure(2);
semilogx(au,kres);