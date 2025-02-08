x = 0:0.1:4;
y1 = erfc(x);
y2 = 1/6*exp(-x.^2)+1/2*exp(-4/3*x.^2);
y3 = 1/2*exp(-x.^2)+1/2*exp(-2*x.^2);
plot(x,10.*log10(y1),'r');
hold on
plot(x,10.*log10(y2),'b');
hold on
plot(x,10.*log10(y3),'g');