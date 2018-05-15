x0=[6.4149; 1.9245; 0.2138;  0.0641; 31.7450];
%x0=[162.9157   48.8747    5.4305    1.6292  295.2840]
options = optimoptions('fsolve','Display','iter','MaxIter',5000,'MaxFunEvals',100000); % Option to display output
[x,fval,exitflag,output,jacob] = fsolve(@modelFunc,x0,options) % Call solver
n1=modelNoise(x)
sqrt(n1)
n2=modelNoise([1398.833829 1078.282193 38.000000 11.000000 58.000000])
sqrt(n2)
n3=modelNoise([2197.853961; 1097.673175; 59.000000; 20.000000; 44.000000])
sqrt(n3)