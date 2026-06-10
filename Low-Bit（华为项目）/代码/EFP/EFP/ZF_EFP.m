function [eta_result,F_result] = ZF_EFP(H_efp, H_config,P,Nu,Ns,Nr)
        % m_bit = 10;
        % H = (randn(Nu*Nr,Nt) + 1i*randn(Nu*Nr,Nt)) / sqrt(2);
        %     config_1 = [2 m_bit 2 0 0 2 m_bit 2 0 0];
        %     config = repmat({config_1}, Nu*Nr, Nt);
        %     [H_efp, H_config] = decToEFP(H,config);

    [H_her,H_her_config] = hermitian_transpose(H_efp, H_config);
    [tmp1,tmp1_config] = EFP_mul_matrix(H_efp,H_config,H_her,H_her_config);
    [tmp1_inv,tmp1_inv_config] = EFP_inv(tmp1,tmp1_config);
    [F,F_config] = EFP_mul_matrix(H_her,H_her_config,tmp1_inv,tmp1_inv_config);

    mask = eye((Nu*Nr));
    for i=1:Nu
        mask(((i-1)*Nr+Ns+1):(i*Nr),((i-1)*Nr+Ns+1):(i*Nr)) = 0;
    end

    mask_config_temp = repmat(H_config(1), Nu*Nr, Nu*Nr);

    [mask,mask_config] = decToEFP(mask,mask_config_temp);

    [F_her,F_her_config] = hermitian_transpose(F,F_config);
    [F_tmp,F_tmp_config] = EFP_mul_matrix(mask,mask_config,F_her,F_her_config);
    [F_head,F_head_config] = EFP_mul_matrix(F,F_config,F_tmp,F_tmp_config);
    
    power = 0;
    [power,power_config] = decToEFP(power,H_config(1));

    for k=1:size(F_head,1)
        [power,power_config] = i_EFP_add(power,power_config,F_head(k,k),F_head_config(k,k));
        power = {power};
        power_config = {power_config};
    end

    [P_efp,P_config] = decToEFP(P,H_config(1));

    [tmp2,tmp2_config] = i_EFP_div(P_efp,P_config,power,power_config);

    [eta,eta_config] = i_EFP_sqrt({tmp2},{tmp2_config});

    [F,F_config] = EFP_mul_matrix_e(F,F_config,{eta},{eta_config});

    
    F_result = EFPTodec(F,F_config);
    eta_result = EFPTodec({eta},{eta_config});

    % F_dec = H'*inv(H*H');
    % mask = eye(Nu*Nr);
    % for i=1:Nu
    %     mask(((i-1)*Nr+Ns+1):(i*Nr),((i-1)*Nr+Ns+1):(i*Nr)) = 0;
    % end
    % eta_dec = sqrt(P/trace(F_dec*mask*F_dec'));
    % F_dec = eta_dec*F_dec;

end