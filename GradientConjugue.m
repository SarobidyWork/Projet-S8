function [x,k] = GradientConjugue(A,b,x0,tolerance,maxiter)
  x = x0;
  r = b - A*x0;
  p = r;
  k = 1;
  while(k<maxiter && norm(r)>tolerance)
    a = r'*r/((p'*A)*p);
    x(:,k+1) = x(:,k) + a*p;
    r1 = r - a*A*p;
    beta = r1'*r1/(r'*r);
    r = r1;
    p = r + beta*p;
    k=k+1;
  endwhile 
endfunction
