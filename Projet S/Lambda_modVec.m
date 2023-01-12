function [T1,T2] = Lambda_modVec(EDP,a,Nx,P)
  [L1,L2] = Lambda2Vec(EDP,a,Nx,P);
  L1(1) = 0; L1(Nx+1) = 0;
## L2(1) = 0; L2(Nx+1) = 0;
  delta_x = (EDP.b - EDP.a)/Nx;
  
  Phi = zeros(Nx+1,1);
  
  alpha = -1/(delta_x)^2;
  v = alpha*ones(1,Nx);
  Dh = diag(v,1);
  Db = diag(v,-1);
  Mat = -alpha*2*eye(Nx+1,Nx+1) + Dh + Db;
  Mat(1,1) = 1;
  Mat(1,2) = 0;
  Mat(Nx+1,Nx) = 0;
  Mat(Nx+1,Nx+1) = 1;
  
  
  
  Phi = Mat\L1;
  T1 = Phi; T2 = L2;
  
  
endfunction
