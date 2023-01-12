function S = Norme(a)
  N = length(a);
  h = (1-0)/N;
  S1 = 0;
  S2 = 0;
  for i=1:N-1
    S1 = S1 + h*((a(i+1,1)-a(i,1))/h)^2;
    S2 = S2 + h*(a(i,2)^2);
  endfor
  S2 = S2 + h*(a(N,2)^2);
  S = S1 + S2;
endfunction
