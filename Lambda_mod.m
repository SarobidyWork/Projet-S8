function [T1,T2] = Lambda_mod(EDP,a,Nx)
  [L1,L2] = Lambda2(EDP,a,Nx);
  L1(1) = 0; L1(Nx+1) = 0;
  delta_x = (EDP.b - EDP.a)/Nx;
  
  Phi = zeros(Nx+1,1);
  
  alpha = -1/(delta_x)^2;
  v = alpha*ones(1,Nx);
  Dh = diag(v,1);
  Db = diag(v,-1);
  Mat = -2*alpha*eye(Nx+1,Nx+1) + Dh + Db;
  Mat(1,1) = 1;
  Mat(1,2) = 0;
  Mat(Nx+1,Nx) = 0;
  Mat(Nx+1,Nx+1) = 1;
  
  Phi = Mat\L1;
  T1 = Phi; T2 = L2;
endfunction
