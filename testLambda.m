%test Lambda

Nx = 10;
a = 0.5;

EDP.a = 0; EDP.b = 1;
EDP.t0 = 0; EDP.T = 4*pi;
EDP.f = 0;
EDP.e0 = @(x) sin(pi*x);
EDP.e1 = @(x) 0;
EDP.ua = @(t) 0;
EDP.ub = @(t) 0;

[L1,L2] = Lambda2(EDP,a,Nx);

longueurL1 = length(L1)
longueurL2 = length(L2)