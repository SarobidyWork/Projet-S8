function [L1,L2] = Lambda2Vec(EDP,a,Nx,P)
  EDP2 = EDP;
  EDP2.e0 = P(:,1); EDP2.e1 = P(:,2);
##  EDP.ua = zeros(lenght(t)); EDP.ub = zeros(length(t));
  [t,x,u] = EqOndeVec(EDP2,a,Nx);
##    for i=length(t):-1:1
##    plot(x,u(:,i));
##    pause(0.5);
##    endfor
  EDP2.e0 = zeros(Nx+1,1);
  EDP2.e1 = zeros(Nx+1,1);
  EDP2.ub = -ones(length(t),1);
  [t,x,y] = EqOnde2Vec(EDP2,a,Nx,u);
##  figure(4)
## for i = size(y, 2) : - 1 : 1 
##     plot(y(:, i))
##     pause(0.1)
## end
##    for i=length(t):-1:length(t)-2
##      figure(3)
##    plot(x,y(:,i));
##    pause(0.1);
##    endfor
  L2 = y(:,1);
  dt = (EDP.T - EDP.t0)/(length(t)-1);
  L1 = -(y(:,2) - y(:,1))/dt;
endfunction
