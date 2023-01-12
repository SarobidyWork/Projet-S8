function y=u0t(x)
% % U0T time derivative initial condition, a Ricker
%   if(x <= 1/4)
%       y = -1; 
%   elseif(1/4<x & x<=1/2)
%       y = 1;
%   else
%       y = 2*(x-1); 
%   end
y = -sin(pi*x)/(4 );
% y = pi*sin(pi*x);
%y=-50*(1.50-x).*exp(-25*(1.50-x).^2); % 6 iterations needed
%y=-50*(1.0-x).*exp(-25*(1.0-x).^2); % 6 iterations needed
% y=-50*(1.5-x).*exp(-25*(1.5-x).^2)+50*(4-x).*exp(-25*(4-x).^2); % 7 iters !
% % % Numerical test proposed by Munch
%     if( 0 < x & x <1/2)
%         y = -x;
%     elseif (x == 0 | x == 1/2)
%         y = 0;
%     else
%         y = 1-x; 
%     end