



delta_dB = 20;

m = 2;
n = 4;


N_t = 4;
N_r = N_t;
Nf = 1e4;


% H = (randn(N_r,N_t) + 1i*randn(N_r,N_t)) / sqrt(2);
% save("H1.mat","H");
H = load("H1.mat","H").H;
% norm(H(:,1),2)^2
% norm(H(:,2),2)^2
% norm(H(:,3),2)^2
% norm(H(:,4),2)^2

ans = 0;
for i = 1:Nf

    delta = 1/(10^(delta_dB/10));
    tmp = (randn(1,1) + 1i*randn(1,1)) * sqrt(delta/2);
    tmp2 = zeros(N_r,N_t);
    tmp2(m,n) = tmp;
    H_head = H + tmp2;
    F = H_head'*inv(H_head*H_head');
    
   
    ans = ans+norm(H*F-eye(N_t),'fro')^2;
end

ans / Nf




