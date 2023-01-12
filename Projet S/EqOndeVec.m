function [t,x,u] = EqOndeVec(EDP,a,Nx)
  delta_x = (EDP.b - EDP.a)/Nx;
  delta_t = a*delta_x;
  Nt = ceil((EDP.T - EDP.t0)/delta_t);
  t = zeros(Nt+1,1);
  x = zeros(Nx+1,1);
  u = zeros(Nx+1,Nt+1);
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
  u(:,1) = EDP.e0;
##  u(:,2) = delta_t*EDP.e1 + u(:,1);
  for j = 2:Nx
    u(j,2) = u(j,1) + delta_t*EDP.e1(j) + (delta_t^2/(2*delta_x^2))*(u(j+1,1)-2*u(j,1)+u(j-1,1));
  end
     
  for i = 3:Nt+1
   u(:,i) = Mat*u(:,i-1) - u(:,i-2);
   u(1,i) = EDP.ua(t(i));
   u(Nx+1,i) = EDP.ub(t(i));
  end
endfunction
