function [v1_efp,v1_config,v2_efp,v2_config]=updatecsvv_last_efp(c1_efp,c1_config,s1_efp,s1_config,c2_efp,c2_config,s2_efp,s2_config,v1_efp,v1_config,v2_efp,v2_config, base, fraction_tables,table)
    [n,~] = size(v1_efp);
    for i = 1:n
        t_efp = v1_efp(i);
        t_config = v1_config(i);
        [c1t_efp,c1t_config] = EFP_mul_matrix(c1_efp,c1_config,t_efp,t_config, base, fraction_tables,table);
        [s1v2_efp,s1v2_config] = EFP_mul_matrix(s1_efp,s1_config,v2_efp(i),v2_config(i), base, fraction_tables,table);
        [v1_efp(i),v1_config(i)] = EFP_add_matrix(c1t_efp,c1t_config,s1v2_efp,s1v2_config, base, fraction_tables,table);
        %v1(i) = c1*t+s1*v2(i);
        [s2t_efp,s2t_config] = EFP_mul_matrix(s2_efp,s2_config,t_efp,t_config, base, fraction_tables,table);
        [c2v2_efp,c2v2_config] = EFP_mul_matrix(c2_efp,c2_config,v2_efp(i),v2_config(i), base, fraction_tables,table);
        [v2_efp(i),v2_config(i)] = EFP_add_matrix(s2t_efp,s2t_config,c2v2_efp,c2v2_config, base, fraction_tables,table);
        % v2(i) = s2*t+c2*v2(i);
    end
end