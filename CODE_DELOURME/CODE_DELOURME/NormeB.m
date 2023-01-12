function S = NormeB(a)
  N = size(a,1) + 1;
  h = (1-0)/(N);
  S1 = 0;
  S2 = 0;
  acoupe1 = a(1:N-1,1)  ;
  acoupe2 = a(1:N-1,2) ; 
  S1 = (acoupe1(1)/h)^2 ;
  for i=2:N-2
    S1 = S1 + ((acoupe1(i+1)-acoupe1(i))/h)^2;
    
  end
  S1= S1 + (acoupe1(N-1)/h)^2 ;
  S1 = S1 *h;
  
  S2 = norm(acoupe2)^2 * h;
  S = S1 + S2;
end
