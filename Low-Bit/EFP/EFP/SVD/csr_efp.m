function [c_efp,c_config,s_efp,s_config,rr_efp,rr_config] = csr_efp(x_efp,x_config,y_efp,y_config, base, fraction_tables,table)
    % global m;
    m = x_config{1}(2);
    config_1 = [2 m 2 0 0 2 m 2 0 0];
    config = repmat({config_1},1,1);
    y = EFPTodec(y_efp,y_config,base,fraction_tables);
    if(y == 0)
        c = 1;
        s = 0;
        [c_efp,c_config] = decToEFP(c,config,base,fraction_tables);
        [s_efp,s_config] = decToEFP(s,config,base,fraction_tables);
        rr_efp = x_efp;
        rr_config = x_config;
    else
        x = EFPTodec(x_efp,x_config,base,fraction_tables);
        if(abs(y) > abs(x))
            [temp_efp,temp_config] = decToEFP(-x,config,base,fraction_tables);
            [tao_efp,tao_config] = EFP_div_matrix_e(temp_efp,temp_config,y_efp,y_config, base, fraction_tables,table);
            [tao2_efp,tao2_config] = EFP_mul_matrix(tao_efp,tao_config,tao_efp,tao_config, base, fraction_tables,table);
            [temp1_efp,temp1_config] = decToEFP(1,config,base,fraction_tables);
            [opt2_efp,opt2_config] = EFP_add_matrix(temp1_efp,temp1_config,tao2_efp,tao2_config, base, fraction_tables,table);
            [s_efp,s_config] = i_EFP_sqrt(opt2_efp,opt2_config, base, fraction_tables,table);
            [temp_efp,temp_config] = decToEFP(-y,config,base,fraction_tables);
            [rr_efp,rr_config] = EFP_mul_matrix(temp_efp,temp_config,{s_efp},{s_config}, base, fraction_tables,table);
            [s_efp,s_config] = EFP_div_matrix_e(temp1_efp,temp1_config,{s_efp},{s_config}, base, fraction_tables,table);
            [c_efp,c_config] = EFP_mul_matrix(s_efp,s_config,tao_efp,tao_config, base, fraction_tables,table);
        else
            [temp_efp,temp_config] = decToEFP(-y,config,base,fraction_tables);
            [tao_efp,tao_config] = EFP_div_matrix_e(temp_efp,temp_config,x_efp,x_config, base, fraction_tables,table);
            [tao2_efp,tao2_config] = EFP_mul_matrix(tao_efp,tao_config,tao_efp,tao_config, base, fraction_tables,table);

            [temp1_efp,temp1_config] = decToEFP(1,config,base,fraction_tables);
            [opt2_efp,opt2_config] = EFP_add_matrix(temp1_efp,temp1_config,tao2_efp,tao2_config, base, fraction_tables,table);
            [c_efp,c_config] = i_EFP_sqrt(opt2_efp,opt2_config, base, fraction_tables,table);
            [rr_efp,rr_config] = EFP_mul_matrix(x_efp,x_config,{c_efp},{c_config}, base, fraction_tables,table);
            [c_efp,c_config] = EFP_div_matrix_e(temp1_efp,temp1_config,{c_efp},{c_config}, base, fraction_tables,table);
            [s_efp,s_config] = EFP_mul_matrix(c_efp,c_config,tao_efp,tao_config, base, fraction_tables,table);
        end
    end
end