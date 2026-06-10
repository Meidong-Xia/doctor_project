function [U,Sigma,V,Plist,stop] =f_bd_precoding(H,P,Nu,Ns,Nr,flg_b,noise,flg_a,rho)
stop = false;
if flg_b == 1 % 64bit
    Sigma = zeros(Ns,Ns,Nu);
    V = [];
    U = zeros(Nr,Ns,Nu);
    Plist = [];
    for j = 1:Nu
        H_tmp = H;
        H_tmp((j-1)*Nr+1:j*Nr,:)=[];
        H_j = H((j-1)*Nr+1:j*Nr,:);
        [~,Sigma_tilde,V_tilde]=svd(H_tmp);
        V_tilde = V_tilde(:,rank(Sigma_tilde)+1:end);
        [U_j,Sigma_j,V_j]=svd(H_j*V_tilde);
        if nnz(Sigma_j)<Ns
            stop = true;
            return
        end
        V_j = V_j(:,1:Ns);
        U_j = U_j(:,1:Ns);
        Sigma_j = Sigma_j(1:Ns,1:Ns);
        V = [V,V_tilde*V_j];
        Sigma(:,:,j)=Sigma_j;
        U(:,:,j) = U_j;
    end
    B=Sigma;
    A=U;
    U = zeros(size(U_j)*Nu);
    Sigma = zeros(size(Sigma_j)*Nu);

    for i=1:Nu
        U((i-1)*Nr+1:i*Nr,(i-1)*Ns+1:i*Ns)=A(:,:,i);
        Sigma((i-1)*Ns+1:i*Ns,(i-1)*Ns+1:i*Ns)=B(:,:,i);
    end
    if flg_a == 1 % 平均功率分配
        Plist = sqrt(P/trace(V*V'))*diag(ones(Nu*Ns,1)); 
    elseif flg_a == 2 % 注水功率分配
        Plist = power_allocation(Sigma,P,noise,flg_b); 
    elseif flg_a == 3 % 数值解法
        Plist = f_power_optimization(P,Sigma,noise,rho,1);
    elseif flg_a == 4 % 二分查找
        Plist = f_power_optimization(P,Sigma,noise,rho,2);
    end
    V = V*Plist;
elseif flg_b == 2 % 32bit
    Sigma = zeros(Ns,Ns,Nu);
    V = [];
    U = zeros(Nr,Ns,Nu);
    Plist = [];
    for j = 1:Nu
        H_tmp = H;
        H_tmp((j-1)*Nr+1:j*Nr,:)=[];
        H_j = H((j-1)*Nr+1:j*Nr,:);
        [~,Sigma_tilde,V_tilde]=cfp32_mysvd(H_tmp);
        V_tilde = V_tilde(:,rank(Sigma_tilde)+1:end);
        [U_j,Sigma_j,V_j]=cfp32_mysvd(single(H_j)*single(V_tilde));
        if nnz(Sigma_j)<Ns
            stop = true;
            return
        end
        V_j = V_j(:,1:Ns);
        U_j = U_j(:,1:Ns);
        Sigma_j = Sigma_j(1:Ns,1:Ns);
        V = [V,single(V_tilde)*single(V_j)];
        Sigma(:,:,j)=Sigma_j;
        U(:,:,j) = U_j;
    end
    B=Sigma;
    A=U;
    U = zeros(size(U_j)*Nu);
    Sigma = zeros(size(Sigma_j)*Nu);

    for i=1:Nu
        U((i-1)*Nr+1:i*Nr,(i-1)*Ns+1:i*Ns)=A(:,:,i);
        Sigma((i-1)*Ns+1:i*Ns,(i-1)*Ns+1:i*Ns)=B(:,:,i);
    end
    if flg_a == 1 % 平均功率分配
        Plist = sqrt(single(P)/single(trace(single(V)*single(V'))))*diag(ones(Nu*Ns,1)); 
    elseif flg_a == 2 % 注水功率分配
        Plist = power_allocation(Sigma,P,noise,flg_b); 
    elseif flg_a == 3 % 数值解法
        Plist = f_power_optimization(P,Sigma,noise,rho,1);
    elseif flg_a == 4 % 二分查找
        Plist = f_power_optimization(P,Sigma,noise,rho,2);
    end
    V = single(V)*single(Plist);
elseif flg_b == 3 % 16bit
    Sigma = zeros(Ns,Ns,Nu);
    V = [];
    U = zeros(Nr,Ns,Nu);
    Plist = [];
    for j = 1:Nu
        H_tmp = H;
        H_tmp((j-1)*Nr+1:j*Nr,:)=[];
        H_j = H((j-1)*Nr+1:j*Nr,:);
        [U_tilde,Sigma_tilde,V_tilde]=cfp16_mysvd_3(H_tmp);
        acc = norm(H_tmp-U_tilde*Sigma_tilde*V_tilde','fro')/norm(H_tmp,'fro');
        if acc > 1 || isnan(acc)
            stop = true;
            return
        end
        V_tilde = V_tilde(:,rank(Sigma_tilde)+1:end);
        H_tmp2 = H_j*V_tilde;
        [U_j,Sigma_j,V_j]=cfp16_mysvd_3(H_tmp2);
        acc = norm(H_tmp2-U_j*Sigma_j*V_j','fro')/norm(H_tmp2,'fro');
        if acc > 1 || isnan(acc)
            stop = true;
            return
        end
        if nnz(Sigma_j)<Ns
            stop = true;
            return
        end
        V_j = V_j(:,1:Ns);
        U_j = U_j(:,1:Ns);
        Sigma_j = Sigma_j(1:Ns,1:Ns);
        [~,~,tmp]=i_hp_matrix_mul(complex(V_tilde),complex(V_j));
        V = [V,tmp];
        Sigma(:,:,j)=Sigma_j;
        U(:,:,j) = U_j;
    end
    B=Sigma;
    A=U;
    U = zeros(size(U_j)*Nu);
    Sigma = zeros(size(Sigma_j)*Nu);

    for i=1:Nu
        U((i-1)*Nr+1:i*Nr,(i-1)*Ns+1:i*Ns)=A(:,:,i);
        Sigma((i-1)*Ns+1:i*Ns,(i-1)*Ns+1:i*Ns)=B(:,:,i);
    end
    if flg_a == 1 % 平均功率分配
        [~,~,V_head] = i_hp_matrix_mul(complex(V),complex(V'));
        power = 0;
        for k=1:size(V_head,1)
            [~,~,power] = i_hp_add(power,V_head(k,k));
        end
        [~,~,tmp2] = i_hp_div(P,power);
        [~,p] = i_hp_sqrt(tmp2);
        Plist = p*diag(ones(Nu*Ns,1));
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