function H = f_channel_generator(Nu,Nr,Nt,kappa,angles)
% f_channel_generator generates Rician channel 
% 
% Based on:
%   
%
% Inputs:
%   - Nu: the number of user
%   - Nt: the number of transmitted antennas
%   - Nr: the number of received antennas
%   - kappa: Rician factor(dB)
%   - angles: the angle(AOA & AOD) of central path of different users
%
% Outputs:
%   - H: the generated Rician channel
%
% log:
%   - initialized by Meidong Xia on 01/22/2024
% 
%
% Rest of the code...
H = zeros(Nu*Nr,Nt);
kappa = 10.^(kappa./10);
sig_Los = kappa/(1+kappa);
sig_NLos = 1/(1+kappa);
for i = 1:Nu
    H1 = sqrt(Nr*Nt)*sqrt(sig_Los/2)*(randn+1i*randn)*f_calArraySteerVector(Nr,angles(i,1))*...
            f_calArraySteerVector(Nt,angles(i,2)).';
    H1 = H1+(randn(Nr,Nt) + 1i*randn(Nr,Nt)) * sqrt(sig_NLos/2);
    H((i-1)*Nr+1:i*Nr,:) = H1;
end
end