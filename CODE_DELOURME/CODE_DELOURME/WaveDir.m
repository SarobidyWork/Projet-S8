function u = WaveDir(u0,u0t,mu,bl,br,a,b,T,n,m)
% WAVE solves the wave equation with Dirichlet boundary conditions
%      u=Wave(u0,u0t,mu,ro,bl,mul,rol,br,mur,ror,a,b,T,n,m) solves
%      the wave equation on the domain [a,b]x[0,T]  with initial condition
%      u0 and u0t and Dirichlet boundary condition bl and br
%      using a finite volume method with explicit marching in time.
%      ro is the density and mu the viscosity. The grid is n points in x and
%      m points in t. Note that the first two points in br and bl 
%      are not used since the first two boundary points are given through 
%      the initial condition.
%      Note also that mu and ro can be vectors of length n.

u = zeros(n,m); % coefficient of the basis with dim = n, m: time step
dx = (b-a)/(n-1);
dt = T/(m-1);
x = (a:dx:b)';
t = (0:dt:T);
s = dt^2/dx^2;
% Create the matrix A
e = ones(n-2,1); 
A = spdiags([s*e 2*(1-s)*e s*e], -1:1, n-2, n-2);
A = full(A);
u(2:n-1,1) = u0; % use initial conditions to compute first two points 
u(1,:) = bl; % Point 0
u(n,:) = br; % Point 1
for i = 1:n-2
   u(i+1,2) = u(i+1,1) + dt*u0t(i); 
end
% Update the boundary point
I = zeros(n-2,1);
for j = 3:m
    I(1,1) = s*u(1,j-1);
    I(n-2,1) = s*u(n,j-1);
    u(2:n-1,j) = A*u(2:n-1,j-1) - u(2:n-1,j-2) + I;
end
##  mi=min(min(u));
## ma=max(max(u));
##  if mi==ma,
##   ma=ma+0.1;
##  end;
##  mesh(flipud(x),t,rot90(u,-1));
##  axis([a,b,0,T,mi,ma]);
##  xlabel('x');
##  ylabel('t');
##  view([-30 30]);
##  pause
##   
##  ti=0;
##  for j=1:m,
##    plot(x,u(:,j),'-');
##    axis([a,b,mi,ma]);
##    title(['t=',num2str(ti)])
##    drawnow
##    ti;
##    pause()
##    ti=ti+dt;
##  end; 














