function [L1,L2] = Lambda(y,EDP,Nt)
  L2 = y(:,1);
  dt = (EDP.T - EDP.t0)/Nt;
  L1 = (y(:,2) - y(:,1))/dt;
endfunction
