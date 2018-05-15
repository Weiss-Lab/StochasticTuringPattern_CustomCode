%i.c.
function u0 = pdex1ic3(x)
% U_ini = 100.0 ;
% x0=0.2; % x0>=0!

%Para = [(1)gamma_f, (2)gamma_n, (3)afa_f, (4)afa_n, (5)Diff_s, (6)gamma_s, 
%(7)K_s, (8)N_l, (9)U_ini, (10)x0]; 

Para=load('Para.txt','r');

U_ini   = Para(9);
x0      = Para(10);

N0      = 1.0;


%new initial conditions
u0 = [ U_ini.* heaviside(x+x0).* heaviside(x0-x) .*4; % S
        0; % F
        N0 % N
        ];

% %original initial conditions
% u0 = [ U_ini.*( 1-sign(x-x0) ).*( 1+sign(x+x0) ); % S
%         0; % F
%         1.0 % N
%         ];

% function u0 = pdex1ic(x)
% U_ini = 1000.0 ;
% u0 = U_ini.*( 1-abs(sign(x)) );

% function u0 = pdex1ic(x)
% u0 = sin(pi*x);
