function y=u0(x)
% U0 initial condition, function of a vector x

% %warning('this is the old u0 function, not to be used any more');
%     if( 1/4 < x & x <3/4)
%         y = 1;
%     elseif (x == 1/4 | x == 3/4)
%         y = 1/2;
%     else
%         y = 0; 
%     end
%     y = 0;
y = zeros(size(x));

 %y=exp(-50*(1-x).^2); 
% y=sin(pi*x);
% discontinuous
% for i=1:length(x),
%     if x(i)>= 0.5 & x(i) <=1.5
%         y(i)=1;
%     else
%         y(i)=0;
%     end
% end

% affine
% for i=1:length(x),
%     if x(i)<= 0.5
%         y(i)=0;
%     elseif x(i)>= 0.5 & x(i) <=1
%         y(i)=2*(x(i)-0.5)
%     elseif x(i) >=1 & x(i) <=1.5
%         y(i)=2*(-x(i)+1.5)
%     else
%         y(i)=0;
%     end
% end
% 
% % Numerical test proposed by Munch
%     if( 0 < x && x <1/2)
%         y = 16*x;
%     elseif (x == 0 || x == 1/2)
%         y = 16*x;
%     else
%         y = 0; 
%     end
% % Numerical test proposed by Munch
%    if( 0 < x & x <1/2)
%         y = x^2;
%     elseif (x == 0 | x == 1/2)
%         y = 0;
%     else
%         y = -(1-x)^2; 
%     end
% y = 0;
% y = sin(pi*x);
