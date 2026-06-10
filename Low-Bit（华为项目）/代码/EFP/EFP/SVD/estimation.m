function [sigma_low_bound,lamda,u]=estimation(B)
[m0,n0] = size(B);
n = min(m0,n0);
lamda = zeros(n,1);
lamda(n) = abs(B(n,n));
for j = n-2:-1:0
    la = i_hp_add(lamda(j+2),abs(B(j+1,j+2)));
    ldl = i_hp_div(lamda(j+2),la);
    lamda(j+1) = i_hp_mul(abs(B(j+1,j+1)),ldl);
%     lamda(j+1) = abs(B(j+1,j+1))*(lamda(j+2)/(lamda(j+2)+abs(B(j+1,j+2))));
end
u = zeros(n,1);
u(1) = abs(B(1,1));
for j=0:n-2
    ua = i_hp_add(u(j+1),abs(B(j+1,j+2)));
    udu = i_hp_div(u(j+1),ua);
    u(j+2) = i_hp_mul(abs(B(j+2,j+2)),udu);
%     u(j+2) = abs(B(j+2,j+2))*(u(j+1)/(u(j+1)+abs(B(j+1,j+2))));
end
B_Infinity = min(lamda);
B_1 = min(u);
sigma_low_bound = min(B_Infinity,B_1);

end