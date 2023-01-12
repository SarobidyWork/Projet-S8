function [e0,e1] = GradientConjugueLambda(EDP,f0,f1,maxiter,tolerance,a,Nx)
  x = zeros(Nx+1,1);
  delta_x = (EDP.b - EDP.a)/Nx;
  x = EDP.a:delta_x:EDP.b;
  e0 = zeros(size(x))'; e1 = zeros(size(x))'; %donnée initiale
  EDP.e0 = e0; EDP.e1 = e1;
  [y0,y1] = Lambda_modVec(EDP,a,Nx,[EDP.e0,EDP.e1]);%on calcule Lambda appliqué à la donnée initiale 
  [f_mod0,f_mod1] = f_mod(EDP,zeros(size(x))',sin(pi*x)',Nx); %second membre
  R = [f_mod0 - y0,f_mod1 - y1];
  P = R;
  k = 1;
  while(k<maxiter && sqrt(Norme(R))>tolerance)
    EDP.e0 = e0; EDP.e1 = e1;
    
    [T1,T2] = Lambda_modVec(EDP,a,Nx,P);
    T = [T1,T2];
##    figure(1);
##    plot(T(1,:));
##    hold on
##    plot(T(2,:));
##    ps = ProduitScalaire(T,P)
    alpha = Norme(R)/ProduitScalaire(T,P);
    residu  = sqrt(Norme(R))
    %alpha = r'*r/((p'*A)*p);
  
   %e0(:,k+1) = e0(:,k) + alpha*P(:,1);
   %e1(:,k+1) = e1(:,k) + alpha*P(:,2);
    e0 = e0 + alpha*P(:,1);
    e1 = e1 + alpha*P(:,2);
  
    %x(:,k+1) = x(:,k) + alpha*p; 
##   figure(2);
##    plot(R(1,:));
##    hold on
##    plot(R(2,:)); 
##    pause(0.1);
    R1 = R - alpha*T;
##    figure(3);
##    plot(R1(1,:));
##    hold on
##    plot(R1(2,:));
##    alpha*T
##    Norme(alpha*T)
    %r1 = r - a*A*p;
    N1 = Norme(R1);
    beta = Norme(R1)/Norme(R);
    %beta = r1'*r1/(r'*r);
    R = R1;
    %r = r1;
    P = R + beta*P;
    %p = r + beta*p;
    k=k+1;
  endwhile 
endfunction
