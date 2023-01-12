function [t,x,y] = EqOnde2Vec(EDP,a,Nx,u)
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
  y(:,Nt+1) = EDP.e0;
  y(:,Nt) = y(:,Nt+1) - delta_t*EDP.e1;
    
  for i = Nt-1:-1:1
   y(:,i) = Mat*y(:,i+1) - y(:,i+2);
   y(1,i) = EDP.ua(i);
   y(Nx+1,i) = EDP.ub(i)*(u(Nx+1,i) - u(Nx,i))/delta_x;
end
##mi=min(min(y));
## ma=max(max(y));
##  if mi==ma,
##   ma=ma+0.1;
##  end;
##    mesh(flipud(x),t,rot90(y,-1));
##  axis([EDP.a,EDP.b,0,EDP.T,mi,ma]);
##  xlabel('x');
##  ylabel('t');
##  view([-30 30]);
##  pause
##   
##  ti=0;
##  for j=1:Nt+1
##    plot(x,y(:,j),'-');
##    axis([EDP.a,EDP.b,mi,ma]);
##    title(['t=',num2str(ti)])
##    drawnow
##    ti
##    pause(0.3)
##    ti=ti+dt;
##  end; 
endfunction
