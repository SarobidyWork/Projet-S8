function [L1,L2] = Lambda2(EDP,a,Nx)
  [t,x,u] = EquationOnde(EDP,a,Nx);
##  for i=1:length(t)
##    plot(x,u(:,i));
##    pause(0.1);
##  endfor
  EDP2 = EDP;
  EDP2.e0 = @(x) 0;
  EDP2.e1 = @(x) 0;
  eta = @(t) 1;
  EDP2.ub = @(t) -eta(t);
  [t,x,y] = EquationOnde2(EDP2,a,Nx,u);
  for i=length(t):-1:1
    plot(x,y(:,i));
    pause(0.5);
  endfor
  L2 = y(:,1);
  dt = (EDP2.T - EDP2.t0)/(length(t)-1);
  L1 = -(y(:,2) - y(:,1))/dt;
  
endfunction
