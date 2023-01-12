function [t,x,y] = EquationOnde2(EDP,a,Nx,u)
  delta_x = (EDP.b - EDP.a)/Nx;
  delta_t = a*delta_x;
  Nt = ceil((EDP.T - EDP.t0)/delta_t);
  t = zeros(Nt+1,1);
  x = zeros(Nx+1,1);
  y = zeros(Nx+1,Nt+1);
  alpha = (delta_t^2)/(delta_x^2);
  
  t = EDP.t0:delta_t:Nt*delta_t;
  x = EDP.a:delta_x:EDP.b;
  
  v = alpha*ones(1,Nx);
  Dh = diag(v,1);
  Db = diag(v,-1);
  Mat = 2*(1-alpha)*eye(Nx+1,Nx+1) + Dh + Db;
  Mat(1,1) = 1;
  Mat(1,2) = 0;
  Mat(Nx+1,Nx) = 0;
  Mat(Nx+1,Nx+1) = 1;
  y(:,Nt+1) = EDP.e0(x);
  y(:,Nt) = y(:,Nt+1) - delta_t*EDP.e1(x)';
    
  for i = Nt-1:-1:1
   y(:,i) = Mat*y(:,i+1) - y(:,i+2);
   y(1,i) = EDP.ua(t(i));
   y(Nx+1,i) = EDP.ub(t(i))*(u(Nx+1,i) - u(Nx,i))/delta_x;
##   derivee = (u(Nx+1,i) - u(Nx,i))/delta_x
  end
endfunction
