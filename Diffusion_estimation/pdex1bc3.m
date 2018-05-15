function [pl,ql,pr,qr] = pdex1bc3(xl,ul,xr,ur,t)

% % original boundary conditions
% pl = ul;
% ql = [0; 0; 0];
% pr = ur;
% qr = [0; 0; 0];

% original boundary conditions
pl = [0; 0; 0];
ql = [1; 1; 1];
pr = [0; 0; 0];
qr = [1; 1; 1];

% function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t)
% pl = 0;
% ql = 1;
% pr = 0;
% qr = 1;

% function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t)
% pl = ul;
% ql = 0;
% pr = pi * exp(-t);
% qr = 1;
