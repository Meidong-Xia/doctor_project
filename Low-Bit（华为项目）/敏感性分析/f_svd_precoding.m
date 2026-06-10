function [U,Sigma,V,Plist,stop] =f_svd_precoding(H,P,Ns,flg_b,noise,flg_a,rho)
stop = false;
if flg_b == 1 % 64bit
    [U,Sigma,V]=svd(H);
    V = V(:,1:Ns);
    U = U(:,1:Ns);
    Sigma = Sigma(1:Ns,1:Ns);
    if flg_a == 1 % 平均功率分配
        Plist = sqrt(P/(Ns))*diag(ones(Ns,1));  
    elseif flg_a == 2 % 注水功率分配
        Plist = power_allocation(Sigma,P,noise,flg_b); 
    elseif flg_a == 3 % 数值解法
        Plist = f_power_optimization(P,Sigma,noise,0,1); 
    elseif flg_a == 4 % 二分查找
        Plist = f_power_optimization(P,Sigma,noise,0,2); 
    end
    V = V*Plist;
elseif flg_b == 2 % 32bit
    [U,Sigma,V]=cfp32_svd_household(H);
    V = V(:,1:Ns);
    U = U(:,1:Ns);
    Sigma = Sigma(1:Ns,1:Ns);
    if flg_a == 1 % 平均功率分配
        Plist = single(sqrt(single(P)/trace(single(V)*single(V'))))*single(diag(ones(Ns,1))); 
    elseif flg_a == 2 % 注水功率分配
        Plist = power_allocation(Sigma,P,noise,flg_b); 
    end
    V = single(V)*single(Plist);
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
    if flg_a == 1 % 平均功率分配
        [~,~,V_head] = i_hp_matrix_mul(complex(V),complex(V'));
        power = 0;
        for k=1:size(V_head,1)
            [~,~,power] = i_hp_add(power,V_head(k,k));
        end
        [~,~,tmp] = i_hp_div(P,power);
        [~,p] = i_hp_sqrt(tmp);
        Plist = p*diag(ones(Ns,1)); 
    elseif flg_a == 2 % 注水功率分配
        Plist = power_allocation(Sigma,P,noise,flg_b); 
    elseif flg_a == 3 % 数值解法
        Plist = f_power_optimization(P,Sigma,noise,rho,1); 
    elseif flg_a == 4 % 二分查找
        Plist = f_power_optimization(P,Sigma,noise,rho,2); 
    end
    [~,~,V] = i_hp_matrix_mul(complex(V),complex(Plist));
end

end