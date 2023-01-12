clear all
close all

Nx = 200;
a = 1;

%Initialisation de la structure EDP

EDP.a = 0; EDP.b = 1;

delta_x = (EDP.b - EDP.a)/Nx;
x = zeros(Nx+1,1);
x = EDP.a:delta_x:EDP.b;

EDP.t0 = 0; EDP.T = 4;
EDP.f = 0;
EDP.e0 = @(x) sin(pi*x);
EDP.e1 = @(x) zeros(size(x));
EDP.ua = @(t) 0;
EDP.ub = @(t) 0;
%solution exacte e0 = sin(pi*x)/4 et e1 = 0
EDP.e1_ex =  sin(pi*x)/4;
EDP.e0_ex =  zeros(size(x));

f0 = zeros(Nx+1,1); f1 = zeros(Nx+1,1);

tolerance = 10^(-4);
maxiter = 50;

[y0,y1] = Lambda_modVec(EDP,a,Nx,[EDP.e0_ex',EDP.e1_ex']);
[f_mod0,f_mod1] = f_mod(EDP,zeros(size(x))',sin(pi*x)',Nx);
##figure(1);
##plot(y0,'b');
##title('comparaison y0 et fmod0');
##hold on;
##plot(f_mod0,'r');
##figure(2);
##plot(y1,'b');
##title('comparaison y1 et f_mod1');
##hold on;
##plot(f_mod1,'r');
[e0,e1] = GradientConjugueLambda(EDP,f0,f1,maxiter,tolerance,a,Nx);
[T1,T2] = Lambda_modVec(EDP,a,Nx,[e0,e1]);

figure(1);
plot(T1);
hold on;
plot(f_mod0,'r');
figure(2);
plot(T2);
hold on;
plot(f_mod1,'r');



##figure(2);
##plot(e1);


##norm(Lambda2Vec(EDP,a,Nx,[e0,e1]) - [f0,f1])