function [sigma_low_bound_efp,sigma_low_bound_config,lamda_efp,lamda_config,u_efp,u_config] = estimation_efp(B_efp,B_config, base, fraction_tables,table)
    % global m;
    m = B_config{1}(2);
    config_1 = [2 m 2 0 0 2 m 2 0 0];
    [m0,n0] = size(B_efp);
    n = min(m0,n0);
    lamda = zeros(n,1);
    config = repmat({config_1},n,1);
    [lamda_efp,lamda_config] = decToEFP(lamda, config,base,fraction_tables);
    [lamda_efp(n),lamda_config(n)] = EFP_abs(B_efp(n,n),B_config(n,n));
    % lamda(n) = abs(EFPTodec(B_efp(n,n),B_config(n,n)));

    % config = repmat({config_1},1,1);
    for j = n-2:-1:0
        [temp_efp,temp_config] = EFP_abs(B_efp(j+1,j+2),B_config(j+1,j+2));
        [la_efp,la_config] = EFP_add_matrix(lamda_efp(j+2),lamda_config(j+2),temp_efp,temp_config, base, fraction_tables,table);
        [ldl_efp,ldl_config] = EFP_div_matrix_e(lamda_efp(j+2),lamda_config(j+2),la_efp,la_config, base, fraction_tables,table);
        [temp_efp,temp_config] = EFP_abs(B_efp(j+1,j+1),B_config(j+1,j+1));
        [lamda_efp(j+1),lamda_config(j+1)] = EFP_mul_matrix(temp_efp,temp_config,ldl_efp,ldl_config, base, fraction_tables,table);
    end
    u = zeros(n,1);
    config = repmat({config_1},n,1);
    [u_efp,u_config] = decToEFP(u, config,base,fraction_tables);
    [u_efp(1),u_config(1)] = EFP_abs(B_efp(1,1),B_config(1,1));
   
    

    config = repmat({config_1},1,1);
    for j = 0:n-2
        [temp_efp,temp_config] = EFP_abs(B_efp(j+1,j+2),B_config(j+1,j+2));
        [ua_efp,ua_config] = EFP_add_matrix(u_efp(j+1),u_config(j+1),temp_efp,temp_config, base, fraction_tables,table);
        [udu_efp,udu_config] = EFP_div_matrix_e(u_efp(j+1),u_config(j+1),ua_efp,ua_config, base, fraction_tables,table);
        [temp_efp,temp_config] = EFP_abs(B_efp(j+2,j+2),B_config(j+2,j+2));
        [u_efp(j+2),u_config(j+2)] = EFP_mul_matrix(temp_efp,temp_config,udu_efp,udu_config, base, fraction_tables,table);
    end
    lamda = EFPTodec(lamda_efp,lamda_config,base,fraction_tables);
    B_Infinity = min(lamda);
    u = EFPTodec(u_efp,u_config,base,fraction_tables);
    B_1 = min(u);
    sigma_low_bound = min(B_Infinity,B_1);
    [sigma_low_bound_efp,sigma_low_bound_config] = decToEFP(sigma_low_bound,config,base,fraction_tables);
end