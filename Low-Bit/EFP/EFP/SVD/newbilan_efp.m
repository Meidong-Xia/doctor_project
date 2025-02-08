function [U_efp,U_config,B_efp,B_config,V_efp,V_config] = newbilan_efp(H_efp, H_config,base,fraction_tables,table)
% mat=mat';

    % global m;

    % m = 11;
    % config_1 = [2 m 2 0 0 2 m 2 0 0];
    % config = repmat({config_1}, Nu*Nr, Nt);
    % [H_efp, H_config] = decToEFP(H, config);
    
    m = H_config{1}(2);
    
    config_1 = [2 m 2 0 0 2 m 2 0 0];
    [size_m,n] = size(H_efp);
    mini = min(size_m,n);
    p0 = zeros(size_m,1);
    p0(1) = 1;
    
    alpha = zeros(mini);
    config = repmat({config_1},mini,mini);
    [alpha_efp,alpha_config] = decToEFP(alpha, config,base,fraction_tables);

    V = zeros(n,mini);
    config = repmat({config_1},n,mini);
    [V_efp,V_config] = decToEFP(V, config,base,fraction_tables);

    beta = zeros(mini+1);
    beta(1) = 1;
    config = repmat({config_1},mini+1,mini+1);
    [beta_efp,beta_config] = decToEFP(beta, config,base,fraction_tables);

    U = zeros(size_m,mini+1);
    U(:,1) = p0;
    config = repmat({config_1}, size_m, mini+1);
    [U_efp,U_config] = decToEFP(U, config,base,fraction_tables);

    B = zeros(mini+1,mini);
    config = repmat({config_1},mini+1,mini);
    [B_efp,B_config] = decToEFP(B, config,base,fraction_tables);

    r = zeros(n,mini+1);
    config = repmat({config_1}, n,mini+1);
    [r_efp,r_config] = decToEFP(r, config,base,fraction_tables);

    P = zeros(size_m,mini);
    config = repmat({config_1},size_m,mini);
    [P_efp,P_config] = decToEFP(P, config,base,fraction_tables);

    for j = 1:mini %j=2

        if j ~= 1
            [H_her,H_her_config] = hermitian_transpose(H_efp, H_config);
            [mu_efp,mu_config] = EFP_mul_matrix(H_her,H_her_config,U_efp(:,j),U_config(:,j), base, fraction_tables,table);
            [bv_efp,bv_config] = EFP_mul_matrix_e(V_efp(:,j-1),V_config(:,j-1),beta_efp(j),beta_config(j), base, fraction_tables,table);
            [r_efp(:,j),r_config(:,j)] = EFP_sub_matrix(mu_efp,mu_config,bv_efp,bv_config, base, fraction_tables,table);
            
            for i = 1:j-1
                [her,her_config] = hermitian_transpose(V_efp(:,i), V_config(:,i));
                [vr_efp,vr_config] = EFP_mul_matrix(her,her_config,r_efp(:,j),r_config(:,j), base, fraction_tables,table);
                [vrv_efp,vrv_config] = EFP_mul_matrix_e(V_efp(:,i),V_config(:,i),vr_efp,vr_config, base, fraction_tables,table);
                [r_efp(:,j),r_config(:,j)] = EFP_sub_matrix(r_efp(:,j),r_config(:,j),vrv_efp,vrv_config, base, fraction_tables,table);
            end

            [alpha_efp(j),alpha_config(j)] = nx_efp(r_efp(:,j),r_config(:,j), base, fraction_tables,table);
            [V_efp(:,j),V_config(:,j)] = EFP_div_matrix_e(r_efp(:,j),r_config(:,j),alpha_efp(j),alpha_config(j), base, fraction_tables,table);
            [mv_efp,mv_config] = EFP_mul_matrix(H_efp, H_config,V_efp(:,j),V_config(:,j), base, fraction_tables,table);
            [aU_efp,aU_config] = EFP_mul_matrix_e(U_efp(:,j),U_config(:,j),alpha_efp(j),alpha_config(j), base, fraction_tables,table);
            [P_efp(:,j),P_config(:,j)] = EFP_sub_matrix(mv_efp,mv_config,aU_efp,aU_config, base, fraction_tables,table);
            
            for i = 1:j
                [her,her_config] = hermitian_transpose(U_efp(:,i), U_config(:,i));
                [UP_efp,UP_config] = EFP_mul_matrix(her,her_config,P_efp(:,j),P_config(:,j), base, fraction_tables,table);
                [UPU_efp,UPU_config] = EFP_mul_matrix_e(U_efp(:,i),U_config(:,i),UP_efp,UP_config, base, fraction_tables,table);
                [P_efp(:,j),P_config(:,j)] = EFP_sub_matrix(P_efp(:,j),P_config(:,j),UPU_efp,UPU_config, base, fraction_tables,table);
            end
            
            [beta_efp(j+1),beta_config(j+1)] = nx_efp(P_efp(:,j),P_config(:,j), base, fraction_tables,table);

            if j~= mini
                [U_efp(:,j+1),U_config(:,j+1)] = EFP_div_matrix_e(P_efp(:,j),P_config(:,j),beta_efp(j+1),beta_config(j+1), base, fraction_tables,table);
            end

        else
            [her,her_config] = hermitian_transpose(H_efp, H_config);
            [r_efp(:,1),r_config(:,1)] = EFP_mul_matrix(her,her_config,U_efp(:,1),U_config(:,1), base, fraction_tables,table);
            
            
            % U_dec = EFPTodec(r_efp(:,1),r_config(:,1));
            [alpha_efp(j),alpha_config(j)] = nx_efp(r_efp(:,1),r_config(:,1), base, fraction_tables,table);
            
            [V_efp(:,j),V_config(:,j)] = EFP_div_matrix_e(r_efp(:,1),r_config(:,1),alpha_efp(j),alpha_config(j), base, fraction_tables,table);
            [mv_efp,mv_config] = EFP_mul_matrix(H_efp, H_config,V_efp(:,j),V_config(:,j), base, fraction_tables,table);
            [aU_efp,aU_config] = EFP_mul_matrix_e(U_efp(:,j),U_config(:,j),alpha_efp(j),alpha_config(j), base, fraction_tables,table);
            [P_efp(:,j),P_config(:,j)] = EFP_sub_matrix(mv_efp,mv_config,aU_efp,aU_config, base, fraction_tables,table);
            
            [her,her_config] = hermitian_transpose(U_efp(:,j), U_config(:,j));
            [UP_efp,UP_config] = EFP_mul_matrix(her,her_config,P_efp(:,j),P_config(:,j), base, fraction_tables,table);
            [UPU_efp,UPU_config] = EFP_mul_matrix_e(U_efp(:,j),U_config(:,j),UP_efp,UP_config, base, fraction_tables,table);
            [P_efp(:,j),P_config(:,j)] = EFP_sub_matrix(P_efp(:,j),P_config(:,j),UPU_efp,UPU_config, base, fraction_tables,table);

            [beta_efp(j+1),beta_config(j+1)] = nx_efp(P_efp(:,j),P_config(:,j), base, fraction_tables,table);
            [U_efp(:,j+1),U_config(:,j+1)] = EFP_div_matrix_e(P_efp(:,j),P_config(:,j),beta_efp(j+1),beta_config(j+1), base, fraction_tables,table);
        end
        % U_dec = EFPTodec(U_efp,U_config);
        % if isnan(U_dec) || j ==1
        %     keyboard;
        % end
    end
    
    % U_dec_11 = EFPTodec(U_efp,U_config);
    
    for i = 1:mini
        B_efp(i,i) = alpha_efp(i);
        B_config(i,i) = alpha_config(i);
        B_efp(i+1,i) = beta_efp(i+1);
        B_config(i+1,i) = beta_config(i+1);
    end

    U_efp = U_efp(:,1:mini);
    U_config = U_config(:,1:mini);
    B_efp = B_efp(1:mini,:);
    B_config = B_config(1:mini,:);
    % B = EFPTodec(B_efp,B_config);
    % B = real(B);
    [B_efp,B_config] = EFP_real(B_efp,B_config);
end