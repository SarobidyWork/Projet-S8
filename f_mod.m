function [f_mod1,f_mod2] = f_mod(EDP,f1,f2,Nx)
   f_mod2 = f2;
   
   Phi = zeros(Nx+1,1);
   
   delta_x = (EDP.b - EDP.a)/Nx;
   alpha = -1/delta_x;
   v = alpha*ones(1,Nx);
   Dh = diag(v,1);
   Db = diag(v,-1);
   Mat = -2*alpha*eye(Nx+1,Nx+1) + Dh + Db;
   Mat(1,1) = 1;
   Mat(1,2) = 0;
   Mat(Nx+1,Nx) = 0;
   Mat(Nx+1,Nx+1) = 1;
   Phi = Mat\f1;
   
   f_mod1 = Phi;
endfunction
