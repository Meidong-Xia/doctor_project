function [G_efp,G_config,U1_efp,U1_config,U2_efp,U2_config,V1_efp,V1_config,V2_efp,V2_config] = mul_22_submatrix_efp(G_efp,G_config,U1_efp,U1_config,U2_efp,U2_config,V1_efp,V1_config,V2_efp,V2_config, base, fraction_tables,table)
    % global m;
    m = G_config{1}(2);
    config_1 = [2 m 2 0 0 2 m 2 0 0];
    % config = repmat({config_1},1,1);
    % G_efp{1}(m+4:end) = dec2bin(0,m+3);
    % G_config{1}(6:end) = [2 m 2 0 0];
    % G = EFPTodec(G_efp,G_config);
    % G = real(G);
    % [sizem,sizen] = size(G);
    % config = repmat({config_1},sizem,sizen);
    [G_efp,G_config] = EFP_real(G_efp,G_config);
    g11_efp = G_efp(1,1);
    g12_efp = G_efp(1,2);
    g22_efp = G_efp(2,2);
    g11_config = G_config(1,1);
    g12_config = G_config(1,2);
    g22_config = G_config(2,2);
    [a_efp,a_config] = EFP_mul_matrix(g11_efp,g11_config,g11_efp,g11_config, base, fraction_tables,table);
    % a_efp{1}(m+4:end) = dec2bin(0,m+3);
    % a_config{1}(6:end) = [2 m 2 0 0];

    [g122_efp,g122_config] = EFP_mul_matrix(g12_efp,g12_config,g12_efp,g12_config, base, fraction_tables,table);
    % g122_efp{1}(m+4:end) = dec2bin(0,m+3);
    % g122_config{1}(6:end) = [2 m 2 0 0];

    [g222_efp,g222_config] = EFP_mul_matrix(g22_efp,g22_config,g22_efp,g22_config, base, fraction_tables,table);
    % g222_efp{1}(m+4:end) = dec2bin(0,m+3);
    % g222_config{1}(6:end) = [2 m 2 0 0];

    [b_efp,b_config] = EFP_add_matrix(g122_efp,g122_config,g222_efp,g222_config, base, fraction_tables,table);
    % b_efp{1}(m+4:end) = dec2bin(0,m+3);
    % b_config{1}(6:end) = [2 m 2 0 0];

    [c_efp,c_config] = EFP_mul_matrix(g11_efp,g11_config,g12_efp,g12_config, base, fraction_tables,table);
    % c_efp{1}(m+4:end) = dec2bin(0,m+3);
    % c_config{1}(6:end) = [2 m 2 0 0];

    config = repmat({config_1},1,1);
    [t_efp,t_config] = decToEFP(2,config,base,fraction_tables);
    [c2_efp,c2_config] = EFP_mul_matrix(c_efp,c_config,t_efp,t_config, base, fraction_tables,table);
    % c2_efp{1}(m+4:end) = dec2bin(0,m+3);
    % c2_config{1}(6:end) = [2 m 2 0 0]; %实数

    c2 = EFPTodec(c2_efp,c2_config,base,fraction_tables);
    if abs(c2) > 1e-3
        [bsa_efp,bsa_config] = EFP_sub_matrix(b_efp,b_config,a_efp,a_config, base, fraction_tables,table);
        % bsa_efp{1}(m+4:end) = dec2bin(0,m+3);
        % bsa_config{1}(6:end) = [2 m 2 0 0]; %实数
        [theta_efp,theta_config] = EFP_div_matrix_e(bsa_efp,bsa_config,c2_efp,c2_config, base, fraction_tables,table);
        % theta_efp{1}(m+4:end) = dec2bin(0,m+3);
        % theta_config{1}(6:end) = [2 m 2 0 0]; %实数
        theta = EFPTodec(theta_efp,theta_config,base,fraction_tables);
        if abs(theta) > 255    % 这里每一步都要取实数
            config = repmat({config_1},1,1);
            [t_efp,t_config] = decToEFP(100,config,base,fraction_tables);
            [theta_efp,theta_config] = EFP_div_matrix_e(theta_efp,theta_config,t_efp,t_config, base, fraction_tables,table);
            [theta2_efp,theta2_config] = EFP_mul_matrix(theta_efp,theta_config,theta_efp,theta_config, base, fraction_tables,table);
            config = repmat({config_1},1,1);
            [t_efp,t_config] = decToEFP(0.0001,config,base,fraction_tables);
            [result_efp,result_config] = EFP_add_matrix(t_efp,t_config,theta2_efp,theta2_config, base, fraction_tables,table);
            [result1_efp,result1_config] = i_EFP_sqrt(result_efp,result_config, base, fraction_tables,table);
            config = repmat({config_1},1,1);
            temp = EFPTodec(theta_efp,theta_config, base, fraction_tables);
            [t_efp,t_config] = decToEFP(abs(temp),config,base,fraction_tables);
            [fenmu_efp,fenmu_config] = EFP_add_matrix(t_efp,t_config,{result1_efp},{result1_config}, base, fraction_tables,table);
            [t_efp,t_config] = decToEFP(sign(temp),config,base,fraction_tables);
            [t_efp,t_config] = EFP_div_matrix_e(t_efp,t_config,fenmu_efp,fenmu_config, base, fraction_tables,table);
            config = repmat({config_1},1,1);
            [t100_efp,t100_config] = decToEFP(100,config,base,fraction_tables);
            [t_efp,t_config] = EFP_div_matrix_e(t_efp,t_config,t100_efp,t100_config, base, fraction_tables,table);
        else
            [theta2_efp,theta2_config] = EFP_mul_matrix(theta_efp,theta_config,theta_efp,theta_config, base, fraction_tables,table);
            config = repmat({config_1},1,1);
            [t_efp,t_config] = decToEFP(1,config,base,fraction_tables);
            [result_efp,result_config] = EFP_add_matrix(t_efp,t_config,theta2_efp,theta2_config, base, fraction_tables,table);
            [result1_efp,result1_config] = i_EFP_sqrt(result_efp,result_config, base, fraction_tables,table);
            config = repmat({config_1},1,1);
            temp = EFPTodec(theta_efp,theta_config,base,fraction_tables);
            [t_efp,t_config] = decToEFP(abs(temp),config,base,fraction_tables);
            [fenmu_efp,fenmu_config] = EFP_add_matrix(t_efp,t_config,{result1_efp},{result1_config}, base, fraction_tables,table);
            [t_efp,t_config] = decToEFP(sign(temp),config,base,fraction_tables);
            [t_efp,t_config] = EFP_div_matrix_e(t_efp,t_config,fenmu_efp,fenmu_config, base, fraction_tables,table);
        end
    else
        t = 1e-3;
        config = repmat({config_1},1,1);
        [t_efp,t_config] = decToEFP(t,config,base,fraction_tables);
    end
    [t2_efp,t2_config] = EFP_mul_matrix(t_efp,t_config,t_efp,t_config, base, fraction_tables,table);
    config = repmat({config_1},1,1);
    [temp_efp,temp_config] = decToEFP(1,config,base,fraction_tables);
    [opt2_efp,opt2_config] = EFP_add_matrix(temp_efp,temp_config,t2_efp,t2_config, base, fraction_tables,table);
    [sqopt_efp,sqopt_config] = i_EFP_sqrt(opt2_efp,opt2_config, base, fraction_tables,table);
    [cs_efp,cs_config] = EFP_div_matrix_e(temp_efp,temp_config,{sqopt_efp},{sqopt_config}, base, fraction_tables,table);
    [sn_efp,sn_config] = EFP_mul_matrix(cs_efp,cs_config,t_efp,t_config, base, fraction_tables,table);
    
    % sn = EFPTodec(sn_efp,sn_config);
    [nsn_efp,nsn_config] = EFP_neg(sn_efp,sn_config);
    temp_efp = [cs_efp,sn_efp;nsn_efp,cs_efp];
    temp_config = [cs_config,sn_config;nsn_config,cs_config];
    [G_efp,G_config] = EFP_mul_matrix(G_efp,G_config,temp_efp,temp_config, base, fraction_tables,table);

    % G = real(G);
    
    [V1_efp,V1_config,V2_efp,V2_config] = updatecsvv_efp(cs_efp,cs_config,sn_efp,sn_config,V1_efp,V1_config,V2_efp,V2_config, base, fraction_tables,table);
    [G112_efp,G112_config] = EFP_mul_matrix(G_efp(1,1),G_config(1,1),G_efp(1,1),G_config(1,1), base, fraction_tables,table);
    [G212_efp,G212_config] = EFP_mul_matrix(G_efp(2,1),G_config(2,1),G_efp(2,1),G_config(2,1), base, fraction_tables,table);
    [gag_efp,gag_config] = EFP_add_matrix(G112_efp,G112_config,G212_efp,G212_config, base, fraction_tables,table);
    [alpha_efp,alpha_config] = i_EFP_sqrt(gag_efp,gag_config, base, fraction_tables,table);
    
    [G122_efp,G122_config] = EFP_mul_matrix(G_efp(1,2),G_config(1,2),G_efp(1,2),G_config(1,2), base, fraction_tables,table);
    [G222_efp,G222_config] = EFP_mul_matrix(G_efp(2,2),G_config(2,2),G_efp(2,2),G_config(2,2), base, fraction_tables,table);
    [gag_efp,gag_config] = EFP_add_matrix(G122_efp,G122_config,G222_efp,G222_config, base, fraction_tables,table);
    [beta_efp,beta_config] = i_EFP_sqrt(gag_efp,gag_config, base, fraction_tables,table);

    [c1_efp,c1_config] = EFP_div_matrix_e(G_efp(1,1),G_config(1,1),{alpha_efp},{alpha_config}, base, fraction_tables,table);
    [c2_efp,c2_config] = EFP_div_matrix_e(G_efp(2,2),G_config(2,2),{beta_efp},{beta_config}, base, fraction_tables,table);
    [s1_efp,s1_config] = EFP_div_matrix_e(G_efp(2,1),G_config(2,1),{alpha_efp},{alpha_config}, base, fraction_tables,table);
    [s2_efp,s2_config] = EFP_div_matrix_e(G_efp(1,2),G_config(1,2),{beta_efp},{beta_config}, base, fraction_tables,table);

    temp_efp = [c1_efp,s1_efp;s2_efp,c2_efp];
    temp_config = [c1_config,s1_config;s2_config,c2_config];
    [G_efp,G_config] = EFP_mul_matrix(temp_efp,temp_config,G_efp,G_config, base, fraction_tables,table);
    
    [U1_efp,U1_config,U2_efp,U2_config] = updatecsvv_last_efp(c1_efp,c1_config,s1_efp,s1_config,c2_efp,c2_config,s2_efp,s2_config,U1_efp,U1_config,U2_efp,U2_config, base, fraction_tables,table);
end