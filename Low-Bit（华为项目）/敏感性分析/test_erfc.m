x = 0:0.1:3;

% y1 = erfc(x);
% 
% y2 = exp(-x.^2);
% 
% y3 = 0.5*exp(-2.*x.^2)+0.5*exp(-x.^2);
% 
% y4 = 2.*qfunc(sqrt(2).*x);
% 
% plot(x,y1,'r');
% hold on 
% plot(x,y2,'b');
% hold on 
% plot(x,y3,'y');
% hold on 
% plot(x,y4,'g');

y1 = erfc(x);

y2 = 1/6*exp(-x.^2)+0.5*exp(-4*x.^2/3);
plot(x,y1,'r');
hold on 
plot(x,y2,'b');
