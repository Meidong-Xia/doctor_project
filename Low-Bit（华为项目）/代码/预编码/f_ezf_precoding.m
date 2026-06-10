function [eta,F,U,Sigma,stop] = f_ezf_precoding(H,P,Nu,Ns,Nr,flg_b,noise,flg_a)
stop  = false;
eta = [];
F = [];
U = [];
Sigma = [];
Nt = size(H,2);
if flg_b == 1 % 64bit
    V_tilde = zeros(Nu*Ns,Nt);
    U = zeros(Nr*Nu,Ns*Nu);
    Sigma = zeros(Ns*Nu,Ns*Nu);
    for i = 1:Nu
        H_tmp = H((i-1)*Nr+1:i*Nr,:);
        [U_i,Sigma_i,V] = svd(H_tmp);
        U((i-1)*Nr+1:i*Nr,(i-1)*Ns+1:i*Ns) = U_i(:,1:Ns);
        Sigma((i-1)*Ns+1:i*Ns,(i-1)*Ns+1:i*Ns) = Sigma_i(1:Ns,1:Ns);
        V_tilde((i-1)*Ns+1:i*Ns,:) = V(:,1:Ns)';
    end
    F = V_tilde'*inv(V_tilde*V_tilde');
    eta = sqrt(P/trace(F*F'));
    F = eta*F;
elseif flg_b == 2 % 32bit
    V_tilde = zeros(Nu*Ns,Nt);
    for i = 1:Nu
        H_tmp = H((i-1)*Nr+1:i*Nr,:);
        [U,Sigma,V] = cfp32_svd_household(H_tmp);
        acc = norm(H_tmp-U*Sigma*V','fro')/norm(H_tmp,'fro');
        if acc > 0.1 || isnan(acc)
            stop = true;
        end
        if stop
            return
        end
        V_tilde((i-1)*Ns+1:i*Ns,:) = V(:,1:Ns)';
    end
    F_inv = GaussianElimination_single(single(V_tilde)*single(V_tilde'));
    F = single(V_tilde')*single(F_inv);
    eta = sqrt(single(P)/single(trace(single(F)*single(F'))));
    F = single(eta)*single(F);
elseif flg_b == 3 % 16bit
    V_tilde = zeros(Nu*Ns,Nt);
    for i = 1:Nu
        H_tmp = H((i-1)*Nr+1:i*Nr,:);
        [U,Sigma,V] = cfp16_mysvd_4(H_tmp);
        acc = norm(H_tmp-U*Sigma*V','fro')/norm(H_tmp,'fro');
        if acc > 0.1 || isnan(acc)
            stop = true;
        end
        if stop
            return
        end
        V_tilde((i-1)*Ns+1:i*Ns,:) = V(:,1:Ns)';
    end
    F_inv_tmp = i_hp_matrix_mul(complex(V_tilde),complex(V_tilde'));
    F_inv = improve_GaussianElimination_half(F_inv_tmp);
    F = i_hp_matrix_mul(complex(V_tilde'),complex(F_inv));
    [~,~,F_head] = i_hp_matrix_mul(complex(F),complex(F'));
    power = 0;
    for k=1:size(F_head,1)
        [~,~,power] = i_hp_add(power,F_head(k,k));
    end
    [~,~,tmp2] = i_hp_div(P,power);
    [~,eta] = i_hp_sqrt(tmp2);
    [~,~,F] = i_hp_matrix_mul(complex(eta*diag(ones(size(F,1),1))),complex(F));
end

end