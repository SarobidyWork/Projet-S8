function [t,x,y] = EqOnde2Vec2(EDP,a,Nx,u)
% WAVE solves the backward wave equation with Dirichlet boundary conditions
%      u=Wave(u0,u0t,mu,ro,bl,mul,rol,br,mur,ror,a,b,T,n,m) solves
%      the wave equation on the domain [a,b]x[0,T]  with initial condition
%      u0 and u0t and Dirichlet boundary condition bl and br
%      using a finite element method with explicit marching in time.
%      ro is the density and mu the viscosity. The grid is n points in x and
%      m points in t. Note that the first two points in br and bl 
%      are not used since the first two boundary points are given through 
%      the initial condition.
%      Note also that mu and ro can be vectors of length n.
n = Nx+1;

delta_x = (EDP.b-EDP.a)/Nx;
delta_t = a*delta_x;
Nt = ceil((EDP.T - EDP.t0)/delta_t);
m = Nt+1;
x = (EDP.a:delta_x:EDP.b)';
t = (0:delta_t:Nt*delta_t);
s = delta_t^2/delta_x^2;
y = zeros(n,m); % coefficient of the basis with dim = n, m: time step
% Create the matrix A
e = ones(n-2,1); 
A = spdiags([s*e 2*(1-s)*e s*e], -1:1, n-2, n-2);
A = full(A);
y(2:n-1,m) = EDP.e0(2:n-1); % use initial conditions to compute first two points 
y(1,:) = EDP.ua(t); % Point 0
y(n,:) = EDP.ub; % Point 1
for i = 1:n-2
   u(i+1,m-1) = u(i+1,m) - delta_t*sin(pi*x(i))/4; 
end
% Update the boundary point
I = zeros(n-2,1);
for j = m-2:-1:1
    I(1,1) = s*y(1,j+1);
    I(n-2,1) = s*y(n,j+1);
    y(2:n-1,j) = A*y(2:n-1,j+1) - y(2:n-1,j+2) + I;
end



 mi=min(min(y));
 ma=max(max(y));
  if mi==ma,
   ma=ma+0.1;
  end;
    mesh(flipud(x),t,rot90(y,-1));
  axis([EDP.a,EDP.b,0,EDP.T,mi,ma]);
  xlabel('x');
  ylabel('t');
  view([-30 30]);
  pause
   
  ti=0;
  for j=1:Nt+1
    plot(x,y(:,j),'-');
    axis([EDP.a,EDP.b,mi,ma]);
    title(['t=',num2str(ti)])
    drawnow
    ti
    pause(0.3)
    ti=ti+dt;
  end
  
  endfunction 