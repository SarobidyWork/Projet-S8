function S = ProduitScalaire(T,P)
  N = length(T);
  h = 1/N;
  S1 = 0;
  S2 = 0;
  for i=1:N-1
    S1 = S1 + (T(i+1,1)-T(i,1))*(P(i+1,1)-P(i,1))/h;
    S2 = S2 + h*T(i,2)*P(i,2);
  endfor
  S2 = S2 + h*T(N,2)*P(N,2);
  S = S1 + S2;
endfunction
