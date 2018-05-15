%

function [c,f,g] = pdex1pde3(x,t,u,DuDx)

% gamma_f = 0.04;
% gamma_n = 0.01;
% afa_f = 1.0;
% afa_n = 1.0;
% 
% Diff_s = 0.1;
% gamma_s = 0.2;%0.024;
% 
% K_s = 3.0e0;
% N_l = 1.0;

%Para = [(1)gamma_f, (2)gamma_n, (3)afa_f, (4)afa_n, (5)Diff_s, (6)gamma_s, 
%(7)K_s, (8)N_l, (9)U_ini, (10)x0, (11)beta_f]; 

Para=load('Para.txt','r');

gamma_f = Para(1);
gamma_n = Para(2);
afa_f   = Para(3);
afa_n   = Para(4);
Diff_s  = Para(5);
gamma_s = Para(6);
K_s     = Para(7);
N_l     = Para(8);
beta_f  = Para(11);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = [1; 1; 1] ;

f = [Diff_s; 0; 0] .* DuDx;

g = [-gamma_s*u(1); ...
    afa_f.*u(3).*u(1)./(K_s+u(1))+beta_f.*u(3)-gamma_f.*u(2); ...
%    afa_n.*u(3).*(1.0-u(3)./N_l)-gamma_n.*u(3) ...
    afa_n.*u(3).*(1.0-u(3)./N_l) ...
    ];

% function [c,f,g] = pdex1pde(x,t,u,DuDx)
% 
% Diff_s = 1.0;
% gamma_s = 1.0;
% 
% c = 1;
% f = Diff_s*DuDx;
% g = -gamma_s*u;

% 
% function [c,f,s] = pdex1pde(x,t,u,DuDx)
% c = pi^2;
% f = DuDx;
% s = 0;