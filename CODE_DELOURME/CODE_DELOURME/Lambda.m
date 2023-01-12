function Y = Lambda(U0,U0t,mu,bl,br,a,b,T,n,m)
% Y = Lambda(u0,u0t,mu,bl,br,a,b,T,n,m) computes one iteration of the
% discrete HUM operator, (uo,uot) -> (wt,-w)
% This function indicates the algorithm

u = WaveDir(U0,U0t,mu,bl,br,a,b,T,n,m); % forward problem
dx = (b-a)/(n-1);
x = (a:dx:b)';
dt = T/(m-1);
t=0:dt:T;
% Define the cut-off function
rho = ones(1,m); 
rho(1) = 0; % Choose the epsilon = 2(Delta t) 
rho(2) = 0; 
rho(m-1) = 0; 
rho(m) = 0;
for i = 1:m
    un(i) = -1/dx*(u(end-1,i)); % boundary at point 1
end

% use WaveDirMod which computes from the two initial steps
w01=zeros(n-2,1);
w01t=zeros(n-2,1);
w = WaveDirBackWard(w01,w01t,mu,bl,un,a,b,T,n,m); % backward problem

##figure(4)
## for i = size(w, 2) : - 1 : 1 
##     plot(w(:, i))
##     pause(0.1)
## end

% minus Discrete Laplacian
    c = ones(n-2,1);
    A = spdiags([-c 2*c -c],-1:1,n-2,n-2);
    A = 1/dx^2*A; 

dyFinal = (w(2:n-1,2)-w(2:n-1,1))/dt;
dyReleve = A\dyFinal ;

Y=[-dyReleve ; w(2:n-1,1)];






