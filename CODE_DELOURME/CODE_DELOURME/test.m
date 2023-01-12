 clear all
 close all
 
 Nx = 200;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test 1 : on verifie que Lamdbda U = f %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% avec la methode de Vuong
 a = 0;
 b = 1 ;
 n = Nx + 1 ;
 dx = (b-a) /Nx ;
 dt = dx ; 
 T = 4;
 m = ceil(T/dt)+1 ;
 mu = 1 ;
 
 bl = 0 ;
 br = 0;
 x = (dx : dx :b-dx) ;
 U0 = u0(x); 
 U0t = u0t(x) ; 
 Y = Lambda(U0,U0t,mu,bl,br,a,b,T,n,m);
 f0 = sin(pi*x)' ; f1 = zeros(size(x));
 
 
 fmod0 = f0 ;
 c = ones(Nx-1,1);
 A = spdiags([-c 2*c -c],-1:1,Nx-1,Nx-1);
 A = 1/dx^2*A;
 fmod1 = -A\f1' ;
 

##  figure(1) 
##  plot(x,Y(Nx:(2*Nx-2)), 'r') 
##  hold on
##  plot(x, fmod0, 'k--', 'linewidth', 2) 
## 
##  legend('cible approche','cible exacte') 
##  title('Test 1')
## 
## 
##  Erreur = [ (Y(1:Nx-1)-fmod1 )  (Y(Nx:(2*Nx-2)) -fmod0)]; 
##  NormeErreur = NormeB(Erreur)
##  
%  
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test 2 : on test la methode du Gradient conjugue pour resoudre Lambda U =
% F 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 %%% Test 2  : on teste la methode du gradient conjugue
  
 
  e0 = zeros(size(x)) ;
  e1 = zeros(size(x)) ;
  U0 =   e0 ;
  U0t = e1 ;
  Y = Lambda(U0,U0t,mu,bl,br,a,b,T,n,m);
  y0 = Y(1:Nx-1);
  y1 = Y(Nx:(2*Nx-2));
  
 
  fmod0 = f0 ;
  c = ones(Nx-1,1);
  A = spdiags([-c 2*c -c],-1:1,Nx-1,Nx-1);
  A = 1/dx^2*A;
  fmod1 = -A\f1' ;
  
  
  %second membre
  R = [ fmod1- y0,fmod0 - y1];
  P = R;
  k = 1;
  
  tolerance = 0.00001 ;
  maxiter = Nx+1;
  tab_residu= zeros(maxiter,1);
  residu =(NormeB(R))
  %tab_residu(1) = residu;
  LambdaP = zeros(Nx-1,2);
  
  while(k<maxiter && sqrt(NormeB(R))>tolerance)
    U0 = P(1:Nx-1,1) ;
    U0t = P(1:Nx-1,2) ;
   
    Y = Lambda(U0,U0t,mu,bl,br,a,b,T,n,m) ;
    
    LambdaP(:, 1) = Y(1:Nx-1);
    LambdaP(:, 2) =  Y(Nx:2*Nx-2);
    
 

    
    prodscal = ProduitScalaireB(LambdaP,P);
    alpha = NormeB(R)/ProduitScalaireB(LambdaP,P);
    residu = (NormeB(R));
    tab_residu(k) = residu;
    %alpha = r'*r/((p'*A)*p);
  
   %e0(:,k+1) = e0(:,k) + alpha*P(:,1);
   %e1(:,k+1) = e1(:,k) + alpha*P(:,2);
    e0 = e0 + alpha*P(:,1);
    e1 = e1 + alpha*P(:,2);
  
    R1 = R - alpha*LambdaP;
    
    
    

    
    beta = NormeB(R1)/NormeB(R);
    
    %beta = r1'*r1/(r'*r);
    
    %r = r1;
    P = R1 + beta*P;
    %p = r + beta*p;
    R = R1;
    k=k+1;
  end
  size(e1)
  e1_exact = -sin(pi*x)/(4 );
  figure(2) 
  plot(x,e1, 'r') 
  hold on
  plot(x, e1_exact, 'k--', 'linewidth', 2) 
  legend('e1 approche','e1 exacte') 
  title('Test du gradient conjugue')
  xlabel('x') 
  ylabel('e(x)')
 erreur0 = zeros(length(x),1);
 lo = size(erreur0)
 erreur1 = zeros(length(x),1);
 erreur0 = abs(zeros(length(x))-e0(:,1));
 la = size(erreur0)
 erreur1 = abs(e1_exact - e1(:,1));
 
 
 figure(1)
 plot(x,erreur0(:,end),'k');
 hold on;
 plot(x,erreur1(:,end),'r');
 legend('erreur e0','erreur e1')
##  
##  iteration = 1:maxiter;
##  figure(3)
##  loglog(iteration(2:k),tab_residu(2:k));
##  xlabel('itérations') 
##  ylabel('résidu') 
##  title('résidu en fonction de l itération')
  bl = 0;
  br = 0;
  
  uetoile = WaveDir(e0(:,1),e1(:,1),mu,bl,br,a,b,T,n,m);
  v = zeros(m,1);
  for i=1:m;
    v(i) = (uetoile(n,i)-uetoile(n-1,i))/dx;
  endfor
  t = (0:dt:T);
##  plot(t,v);
##  xlabel('t') 
##  ylabel('v(t)') 
##  title('fonction contrôle v')
##  x1 = (a:dx:b)'
##  yreel = WaveDir(U0,v,mu,bl,br,a,b,T,n,m);
##   mi=min(min(yreel));
## ma=max(max(yreel));
##  if mi==ma,
##   ma=ma+0.1;
##  end;
##  mesh(flipud(x1),t,rot90(yreel,-1));
##  axis([a,b,0,T,mi,ma]);
##  xlabel('x');
##  ylabel('t');
##  zlabel('y');
##  view([-30 30]);
##  pause
##   
##  ti=0;
##  for j=1:m,
##    plot(x1,yreel(:,j),'-');
##    axis([a,b,mi,ma]);
##    title(['t=',num2str(ti)])
##    drawnow
##    ti;
##    pause()
##    ti=ti+dt;
##  end; 
##
##
##
