function a = f_calArraySteerVector(M,theta,f,fc)

a = (exp(-1i*pi.*(0:(M-1)).*sin(theta)*(1+f/fc))./sqrt(M)).';

end
