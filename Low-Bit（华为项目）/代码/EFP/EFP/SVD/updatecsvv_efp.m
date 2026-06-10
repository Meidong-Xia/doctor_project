function [v1_efp,v1_config,v2_efp,v2_config] = updatecsvv_efp(c_efp,c_config,s_efp,s_config,v1_efp,v1_config,v2_efp,v2_config, base, fraction_tables,table)
    [n,~] = size(v1_efp);
    for i = 1:n
        t_efp = v1_efp(i);
        t_config = v1_config(i);
        [ct_efp,ct_config] = EFP_mul_matrix(c_efp,c_config,t_efp,t_config, base, fraction_tables,table);
        [sv2_efp,sv2_config] = EFP_mul_matrix(s_efp,s_config,v2_efp(i),v2_config(i), base, fraction_tables,table);
        [v1_efp(i),v1_config(i)] = EFP_sub_matrix(ct_efp,ct_config,sv2_efp,sv2_config, base, fraction_tables,table);
        % v1(i) = c*t-s*v2(i);
        [st_efp,st_config] = EFP_mul_matrix(s_efp,s_config,t_efp,t_config, base, fraction_tables,table);
        [cv2_efp,cv2_config] = EFP_mul_matrix(c_efp,c_config,v2_efp(i),v2_config(i), base, fraction_tables,table);
        [v2_efp(i),v2_config(i)] = EFP_add_matrix(st_efp,st_config,cv2_efp,cv2_config, base, fraction_tables,table);
        %v2(i)=s*t+c*v2(i);
    end
end