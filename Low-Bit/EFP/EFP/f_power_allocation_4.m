function Plist = f_power_allocation_4(U,Sigma,V,Ns,sigma_n,P,m_bit, base, fraction_tables,table,rho)
iter = 20;
Sigma = EFP_mul_matrix_2(Sigma,Sigma,m_bit, base, fraction_tables,table);
R = zeros(Ns,Ns);
for i = 1:Ns
    for j = 1:Ns
        if i~=j
            R(i,j) = rho/(1+rho);
        end
    end
end
h = zeros(Ns,1);
for i = 1:Ns
    h(i,1) = sigma_n/(Sigma(i,i)*(1+rho));
end
D = [R,h;ones(1,Ns)*R/P,ones(1,Ns)*h/P];
[y,~] = f_power_iteration(D,iter,m_bit, base, fraction_tables,table);
idx = y(end);
y = EFP_div_matrix_e_2(y,idx,m_bit, base, fraction_tables,table);
p = y(1:end-1);
Plist = diag(p);


end

function [u,lambda] = f_power_iteration(H,iter,m_bit, base, fraction_tables,table)
    [~,n]=size(H);
    u = ones(n,1)/sqrt(n);
    for i=1:iter
        v = EFP_mul_matrix_2(H,u,m_bit, base, fraction_tables,table);
        v_norm = EFP_mul_matrix_2(v',v,m_bit, base, fraction_tables,table);
        v_norm_sqrt = i_EFP_sqrt_2(v_norm,m_bit, base, fraction_tables,table);
        u = EFP_div_matrix_e_2(v,v_norm_sqrt,m_bit, base, fraction_tables,table);
    end
    lambda = (u'*H*u)/(u'*u);
end
