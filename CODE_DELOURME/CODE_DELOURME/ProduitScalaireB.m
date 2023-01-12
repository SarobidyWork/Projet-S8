function S = ProduitScalaireB(a,b)
  
  N = size(a,1) + 1;
  h = (1-0)/(N);
  S1 = 0;
  S2 = 0;
  acoupe1 = a(1:N-1,1)  ;
  acoupe2 = a(1:N-1,2) ; 
  bcoupe1 = b(1:N-1,1)  ;
  bcoupe2 = b(1:N-1,2)  ;
%   figure(4)
%   plot(acoupe1, 'r')
%   hold on
%   plot(bcoupe1, 'm')
%   legend('a', 'b')
%   hold off
%   figure(5)
%    plot(acoupe2, 'k')
%   hold on
%   plot(bcoupe2, 'g')
%    legend('a', 'b')
%   hold off
%   pause()
  S1 = (acoupe1(1)/h)* (bcoupe1(1)/h);
  for i=2:N-2
    S1 = S1 + ((acoupe1(i+1)-acoupe1(i))/h) *((bcoupe1(i+1)-bcoupe1(i))/h)  ;
    
  end
  
  S1 = S1 + (acoupe1(N-1)/h)* (bcoupe1(N-1)/h) ; 
  S1 = S1 *h;
  S2 = sum(acoupe2.*bcoupe2)*h;

  S = S1 + S2;
end
