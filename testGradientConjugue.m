n = 100;
D = -ones(n-1,1);
A = 2*eye(n,n) + diag(D,-1) + diag(D,1);
x0 = ones(n,1);
b = ones(n,1);
tolerance = 10^(-14);
maxiter = 20000;
[x,k] = GradientConjugue(A,b,x0,tolerance,maxiter);

norm(A*x(:,end) - b);


##N = 5;
##n = 10;
##K = zeros(1,N);
##r = zeros(1,N);
##n2 = zeros(1,N);
##for i=1:N
##  n2(i) = n;
##  D = -ones(n-1,1);
##  A = 2*eye(n,n) + diag(D,-1) + diag(D,1);
##  x0 = ones(n,1);
##  b = ones(n,1);
##  [x,k] = GradientConjugue(A,b,x0,tolerance,maxiter);
##  K(i) = k;
##  r(i) = norm(A*x(:,end) - b);
##  n = 2*n;
##endfor
##
##K
##r
##
##plot(n2,K);
##xlabel('Taille de la matrice');
##ylabel('Nombre d iterations');