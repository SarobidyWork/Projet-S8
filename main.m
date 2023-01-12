clear all
close all

Nx = 10;
a = 0.5;

%Initialisation de la structure EDP

EDP.a = 0; EDP.b = 1;
EDP.t0 = 0; EDP.T = 4*pi;
EDP.f = 0;
EDP.e0 = @(x) sin(pi*x);
EDP.e1 = @(x) 0;
EDP.ua = @(t) 0;
EDP.ub = @(t) 0;

[t,x,u] = EquationOnde(EDP,a,Nx);
uex = @(x,t) sin(pi*x).*cos(pi*t);

##for i=1:length(t)
##  plot(x,u(:,i));
##  hold on;
##  plot(x,uex(x,t(i)));
##  axis([0 1 -1 1]);
##  pause(0.1);
##  hold off;
##end

##Nt = ceil((EDP.T - EDP.t0)/(a*((EDP.b - EDP.a)/Nx)));
##N = 5;
##erreur = zeros(1,N);
##delta_x = zeros(1,N); 
##for j=1:N
##  e = zeros(1,Nt+1);
##  [t,x,u] = EquationOnde(EDP,a,Nx);
##  delta_x(j) = (EDP.b - EDP.a)/Nx;
##  for i=1:Nt+1
##    e(i) = max(abs(u(:,i) - uex(x(:),t(i))));
##  end
##  erreur(j) = max(e);
##  Nx = Nx + 50;
##  Nt = ceil((EDP.T - EDP.t0)/(a*((EDP.b - EDP.a)/Nx)));
##end
##loglog(delta_x,erreur);
##hold on;
##loglog(delta_x,delta_x,'k-o',delta_x,delta_x.^2,'k-*');
##xlabel('Delta x');
##ylabel('Erreur');
##legend({'erreur','O(h)','O(h^2)'},'Location','southeast');

%Connaissant u, on obtient y

EDP.e0 = @(x) 0;
EDP.e1 = @(x) 0;
eta = @(t) 1;
EDP.ub = @(t) -eta(t);

%On resout avec l'algorithme

[t,x,y] = EquationOnde2(EDP,a,Nx,u);
y
yex_bord = @(t) pi*cos(pi*t);

##Nt = ceil((EDP.T - EDP.t0)/(a*((EDP.b - EDP.a)/Nx)));
##N = 10;
##erreur = zeros(1,N);
##delta_x = zeros(1,N); 
##for j=1:N
##  e = zeros(1,Nt+1);
##  [t,x,u] = EquationOnde(EDP,a,Nx);
##  [t,x,y] = EquationOnde2(EDP,a,Nx,u);
##  delta_x(j) = (EDP.b - EDP.a)/Nx;
##  for i=1:Nt+1
##    e(i) = norm(y(Nx+1,i) - yex_bord(t(i)));
##  end
##  erreur(j) = max(e);
##  Nx = Nx + 50;
##  Nt = ceil((EDP.T - EDP.t0)/(a*((EDP.b - EDP.a)/Nx)));
##end
##loglog(delta_x,erreur);
##hold on;
##loglog(delta_x,delta_x,'k-o',delta_x,delta_x.^2,'k-*');
##xlabel('Delta x');
##ylabel('Erreur');
##legend({'erreur','O(h)','O(h^2)'},'Location','southeast');
