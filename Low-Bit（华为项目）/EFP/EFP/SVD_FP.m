function [U,Sigma,V,Plist,stop] =SVD_FP(H,P,Ns,flg_b,sigma_n,flg_a,rho)
stop = false;
flg_b = 3;
if flg_b == 1 % 64bit
    [U,Sigma,V]=svd(H);
    V = V(:,1:Ns);
    U = U(:,1:Ns);
    Sigma = Sigma(1:Ns,1:Ns);
    if flg_a == 1   % 平均功率分配

        % Plist = f_power_allocation_1(U,Sigma,V,Ns,P,m_bit, base, fraction_tables,table);
        Plist = single(sqrt(single(P)/trace(single(V)*single(V'))))*single(diag(ones(Ns,1))); 
        % Plist = diag(ones(Ns,1));
    elseif flg_a == 2 % 最小BER
        p = f_power_ber2(64,Ns,Sigma,rho,sigma_n,P);
        % p2 = f_power_ber(64,Ns,Sigma,rho,sigma_n,P);
        Plist = diag(sqrt(p));
    elseif flg_a == 3 % 最小化最大MSE
        % Plist = f_power_allocation_3(U,Sigma,V,Ns,sigma_n,P,m_bit, base, fraction_tables,table,rho);
        p = f_power_mse2(Ns,Sigma,rho,sigma_n,P);
        Plist = diag(sqrt(p));
    elseif flg_a == 4 % 最大化最小SINR
        % Plist = f_power_allocation_4(U,Sigma,V,Ns,sigma_n,P,m_bit, base, fraction_tables,table,rho);
        % p1 = f_power_blerNP(Ns,Sigma,rho,sigma_n,P);
        p = f_power_bler2(Ns,Sigma,rho,sigma_n,P);
        Plist = diag(sqrt(p));
    else % 最小化BER
    end
    % V = V*Plist;
elseif flg_b == 2 % 32bit
    [U,Sigma,V]=cfp32_svd_lanczos(H);
    acc = norm(H-U*Sigma*V','fro')/norm(H,'fro');
    if acc > 0.1 || isnan(acc)
        stop = true;
    end
    if nnz(Sigma) < Ns
        stop = true;
    end
    if stop
        Plist = [];
        return
    end
    V = V(:,1:Ns);
    U = U(:,1:Ns);
    Sigma = Sigma(1:Ns,1:Ns);
    if flg_a == 1 % 平均功率分配
        Plist = single(sqrt(single(P)/trace(single(V)*single(V'))))*single(diag(ones(Ns,1))); 
    elseif flg_a == 2 % 注水功率分配
        Plist = power_allocation(Sigma,P,noise,flg_b); 
    end
    % V = single(V)*single(Plist);
elseif flg_b == 3 % 16bit
    [U,Sigma,V]=cfp16_mysvd_4(H);
    acc = norm(H-U*Sigma*V','fro')/norm(H,'fro');
    if acc > 0.1 || isnan(acc)
        stop = true;
    end
    if nnz(Sigma) < Ns
        stop = true;
    end
    if stop
        Plist = [];
        return
    end
    V = V(:,1:Ns);
    U = U(:,1:Ns);
    Sigma = Sigma(1:Ns,1:Ns);
    if flg_a == 1   % 平均功率分配

        % Plist = f_power_allocation_1(U,Sigma,V,Ns,P,m_bit, base, fraction_tables,table);
        Plist = single(sqrt(single(P)/trace(single(V)*single(V'))))*single(diag(ones(Ns,1))); 
        % Plist = diag(ones(Ns,1));
    elseif flg_a == 2 % 最小BER
        p = f_power_ber2(64,Ns,Sigma,rho,sigma_n,P);
        % p2 = f_power_ber(64,Ns,Sigma,rho,sigma_n,P);
        Plist = diag(sqrt(p));
    elseif flg_a == 3 % 最小化最大MSE
        % Plist = f_power_allocation_3(U,Sigma,V,Ns,sigma_n,P,m_bit, base, fraction_tables,table,rho);
        p = f_power_mse2(Ns,Sigma,rho,sigma_n,P);
        Plist = diag(sqrt(p));
    elseif flg_a == 4 % 最大化最小SINR
        % Plist = f_power_allocation_4(U,Sigma,V,Ns,sigma_n,P,m_bit, base, fraction_tables,table,rho);
        % p1 = f_power_blerNP(Ns,Sigma,rho,sigma_n,P);
        p = f_power_bler2(Ns,Sigma,rho,sigma_n,P);
        Plist = diag(sqrt(p));
    else % 最小化BER
    end
    % [~,~,V] = i_hp_matrix_mul(complex(V),complex(Plist));
end

end

% function p = f_power_ber2(M,Ns,Sigma,rho,noise,P)
% tmp1 = 2*(1-1/sqrt(M))/log2(M);
% cvx_begin quiet
%     cvx_solver mosek
%     variable p(Ns) nonnegative
%     expression a(Ns)
%     expression b(Ns)
%     expression f
%     for i = 1:Ns
%         a(i) = p(i)*Sigma(i,i)^2;
%         b(i) = (M-1)*(Sigma(i,i)^2*rho*P+noise);
%     end
%     for i = 1:Ns
%         f =  f + tmp1*(0.5*exp(-2*a(i)/b(i))+(1/6)*exp(-3*a(i)/(2*b(i))));
%     end
%     minimize(f)
%     subject to
%         sum(p) == P
% cvx_end
% end

% function p = f_power_mse2(Ns,Sigma,rho,noise,P) % 二分
% % u_min = 0;
% % u_max = 10*P;
% eposilon = 1e-5;
% p = zeros(Ns,1);
% % while 1
% %     u = (u_max+u_min)/2;
% %     f = 0;
% %     for i = 1:Ns
% %         f = f+max((P*rho+noise/Sigma(i,i)^2)/u*P,eposilon);
% %     end
% %     if abs(f-P) <= 1e-5
% %         break
% %     end
% %     if f>P
% %         u_min = u;
% %     elseif f<P
% %         u_max = u;
% %     end
% % end
% tmp = 0;
% for i = 1:Ns
%     tmp = tmp+P*rho+noise/Sigma(i,i)^2;
% end
% u = tmp/P^2;
% for i = 1:Ns
%     p(i) = max((P*rho+noise/Sigma(i,i)^2)/(u*P),eposilon);
% end
% end
% 
% function p = f_power_bler2(Ns,Sigma,rho,noise,P)
% cvx_begin quiet
%     cvx_solver mosek
%     variable p(Ns) nonnegative
%     expression a(Ns)
%     expression b(Ns)
%     variable f
%     for i = 1:Ns
%         a(i) = p(i)*Sigma(i,i)^2;
%         b(i) = Sigma(i,i)^2*rho*P+noise;
%     end
% 
%     maximize(f);
%     subject to
%         for i = 1:Ns
%             a(i)/b(i) >= f
%         end
%         sum(p) == P
% cvx_end
% end

function p = f_power_mse2(Ns,Sigma,rho,noise,P) % 二分
% u_min = 0;
% u_max = 10*P;
eposilon = 1e-5;
p = zeros(Ns,1);
% while 1
%     u = (u_max+u_min)/2;
%     f = 0;
%     for i = 1:Ns
%         f = f+max((P*rho+noise/Sigma(i,i)^2)/u*P,eposilon);
%     end
%     if abs(f-P) <= 1e-5
%         break
%     end
%     if f>P
%         u_min = u;
%     elseif f<P
%         u_max = u;
%     end
% end
tmp = 0;
for i = 1:Ns
    tmp = tmp+P*rho+noise/Sigma(i,i)^2;
end
u = tmp/P^2;
for i = 1:Ns
    p(i) = max((P*rho+noise/Sigma(i,i)^2)/(u*P),eposilon);
end
end

function p = f_power_bler2(Ns,Sigma,rho,noise,P)
cvx_begin quiet
    cvx_solver mosek
    variable p(Ns) nonnegative
    expression a(Ns)
    expression b(Ns)
    variable f
    for i = 1:Ns
        a(i) = p(i)*Sigma(i,i)^2;
        b(i) = Sigma(i,i)^2*rho*P+noise;
    end
    
    maximize(f);
    subject to
        for i = 1:Ns
            a(i)/b(i) >= f
        end
        sum(p) == P
cvx_end
end

function p = f_power_blerNP(Ns,Sigma,rho,noise,P)
iter = 10;
R = ones(Ns,Ns)*rho/(1+rho);
for i = 1:Ns
    R(i,i) = 0;
end
h = zeros(Ns,1);
for i = 1:Ns
    h(i) = noise/(Sigma(i,i)^2*(1+rho));
end
C = [eye(Ns),zeros(Ns,1);ones(1,Ns),-P];
B = [R,h;zeros(1,Ns),0];
D = inv(C)*B;
[y,lambda] = f_power_iteration(D,iter);
idx = y(end);
y = y/idx;
p = y(1:end-1);
end

function [u,lambda] = f_power_iteration(H,iter)
[~,n]=size(H);
u = ones(n,1)/sqrt(n);
for i=1:iter
    v = H*u; 
    u = v/(norm(v,2));
end
lambda = (u'*H*u)/(u'*u);
end

function p = f_power_ber2(M,Ns,Sigma,rho,noise,P)
tmp1 = 2*(1-1/sqrt(M))/log2(M);
cvx_begin quiet
    cvx_solver mosek
    variable p(Ns) nonnegative
    expression a(Ns)
    expression b(Ns)
    expression f
    for i = 1:Ns
        a(i) = p(i)*Sigma(i,i)^2;
        b(i) = (M-1)*(Sigma(i,i)^2*rho*P+noise);
    end
    for i = 1:Ns
        f =  f + tmp1*(0.5*exp(-2*a(i)/b(i))+(1/6)*exp(-3*a(i)/(2*b(i))));
    end
    minimize(f)
    subject to
        sum(p) == P
cvx_end
end

function p = f_power_ber(M,Ns,Sigma,rho,noise,P)
p = ones(Ns,1);
f_pre = Inf;
iter = 100;
for i = 1:iter
    y = f_power_1(M,Ns,Sigma,rho,noise,P,p);
    [f,p] = f_power_2(M,Ns,Sigma,rho,noise,P,y);
    if abs(f-f_pre) < 1e-5
        break
    end
    f_pre = f;
end
end

function [f,p] = f_power_2(M,Ns,Sigma,rho,noise,P,y)
tmp1 = 2*(1-1/sqrt(M))/log2(M);
cvx_begin quiet
    cvx_solver mosek
    variable p(Ns) nonnegative
    expression a(Ns)
    expression b(Ns)
    expression f
    for i = 1:Ns
        a(i) = p(i)*Sigma(i,i)^2*(1); % +rho
        b(i) = (M-1)*(Sigma(i,i)^2*rho*(P)+noise); % -p(i)
    end
    for i = 1:Ns
        f =  f + tmp1*(0.5*exp(y(i,1)^2*b(i)-2*y(i,1)*sqrt(2*a(i)))+(1/6)*exp(2*y(i,2)^2*b(i) ...
            -2*y(i,2)*sqrt(3*a(i))));
    end
    minimize(f)
    subject to
        sum(p) == P
cvx_end
end

function y = f_power_1(M,Ns,Sigma,rho,noise,P,p)
a = zeros(Ns,1);
b = zeros(Ns,1);
y = zeros(Ns,2);
for i = 1:Ns
    a(i) = p(i)*Sigma(i,i)^2*(1); % +rho
    b(i) = (M-1)*(Sigma(i,i)^2*rho*(P)+noise); % -p(i)
end
for k = 1:Ns
    y(k,1) = sqrt(2*a(k))/b(k);
    y(k,2) = sqrt(3*a(k))/(2*b(k));
end
end

