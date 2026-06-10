function [tri_efp,tri_config,U_efp,U_config,V_efp,V_config]=QR_Wilkinson_shift_Iteration_once_efp(tri_efp,tri_config,U_efp,U_config,V_efp,V_config,i_min,i_max,shift_efp,shift_config, base, fraction_tables,table)
    % global m;
    m = tri_config{1}(1);
    config_1 = [2 m 2 0 0 2 m 2 0 0];
    config = repmat({config_1},1,1);
    n = i_max+1;
    [tri2_efp,tri2_config] = EFP_mul_matrix(tri_efp(i_min+1,i_min+1),tri_config(i_min+1,i_min+1),tri_efp(i_min+1,i_min+1),tri_config(i_min+1,i_min+1), base, fraction_tables,table);
    [x_efp,x_config] = EFP_sub_matrix(tri2_efp,tri2_config,shift_efp,shift_config, base, fraction_tables,table);
    [y_efp,y_config] = EFP_mul_matrix(tri_efp(i_min+1,i_min+1),tri_config(i_min+1,i_min+1),tri_efp(i_min+1,i_min+2),tri_config(i_min+1,i_min+2), base, fraction_tables,table);

    for k = i_min:n-2
        [c_efp,c_config,s_efp,s_config,r_efp,r_config] = csr_efp(x_efp,x_config,y_efp,y_config, base, fraction_tables,table);
        [temp0_efp,temp0_config] = decToEFP(0,config,base,fraction_tables);
        temp1_efp = [tri_efp(k+1,k+1),tri_efp(k+1,k+2);temp0_efp,tri_efp(k+2,k+2)];
        temp1_config = [tri_config(k+1,k+1),tri_config(k+1,k+2);temp0_config,tri_config(k+2,k+2)];
        % s = EFPTodec(s_efp,s_config);
        [ns_efp,ns_config] = EFP_neg(s_efp,s_config);
        temp2_efp = [c_efp,s_efp;ns_efp,c_efp];
        temp2_config = [c_config,s_config;ns_config,c_config];
        [A_efp,A_config] = EFP_mul_matrix(temp1_efp,temp1_config,temp2_efp,temp2_config, base, fraction_tables,table);
        
        x_efp = A_efp(1,1);
        x_config = A_config(1,1);
    
        tri_efp(k+1,k+2) = A_efp(1,2);
        tri_config(k+1,k+2) = A_config(1,2);
        y_efp = A_efp(2,1);
        y_config = A_config(2,1);

        tri_efp(k+2,k+2) = A_efp(2,2);
        tri_config(k+2,k+2) = A_config(2,2);

        [V_efp(:,k+1),V_config(:,k+1),V_efp(:,k+2),V_config(:,k+2)] = updatecsvv_efp(c_efp,c_config,s_efp,s_config,V_efp(:,k+1),V_config(:,k+1),V_efp(:,k+2),V_config(:,k+2), base, fraction_tables,table);
        if(k > i_min)
            tri_efp(k,k+1) = r_efp;
            tri_config(k,k+1) = r_config;
        end
        [c_efp,c_config,s_efp,s_config,r_efp,r_config] = csr_efp(x_efp,x_config,y_efp,y_config, base, fraction_tables,table);
        tri_efp(k+1,k+1) = r_efp;
        tri_config(k+1,k+1) = r_config;
        [U_efp(:,k+1),U_config(:,k+1),U_efp(:,k+2),U_config(:,k+2)] = updatecsvv_efp(c_efp,c_config,s_efp,s_config,U_efp(:,k+1),U_config(:,k+1),U_efp(:,k+2),U_config(:,k+2), base, fraction_tables,table);
        
        % s = EFPTodec(s_efp,s_config);
        [ns_efp,ns_config] = EFP_neg(s_efp,s_config);
        temp2_efp = [c_efp,ns_efp;s_efp,c_efp];
        temp2_config = [c_config,ns_config;s_config,c_config];
        [temp0_efp,temp0_config] = decToEFP(0,config,base,fraction_tables);

        if k ~= n-2
            temp1_efp = [tri_efp(k+1,k+2),temp0_efp;tri_efp(k+2,k+2),tri_efp(k+2,k+3)];
            temp1_config = [tri_config(k+1,k+2),temp0_config;tri_config(k+2,k+2),tri_config(k+2,k+3)];
            [A_efp,A_config] = EFP_mul_matrix(temp2_efp,temp2_config,temp1_efp,temp1_config, base, fraction_tables,table);
            
            x_efp = A_efp(1,1);
            y_efp = A_efp(1,2);
            tri_efp(k+2,k+2) = A_efp(2,1);
            tri_efp(k+2,k+3) = A_efp(2,2);
            x_config = A_config(1,1);
            y_config = A_config(1,2);
            tri_config(k+2,k+2) = A_config(2,1);
            tri_config(k+2,k+3) = A_config(2,2);
        else
            temp1_efp = [tri_efp(n-1,n);tri_efp(n,n)];
            temp1_config = [tri_config(n-1,n);tri_config(n,n)];
            [A_efp,A_config] = EFP_mul_matrix(temp2_efp,temp2_config,temp1_efp,temp1_config, base, fraction_tables,table);
            
            tri_efp(n-1,n) = A_efp(1);
            tri_efp(n,n) = A_efp(2);
            tri_config(n-1,n) = A_config(1);
            tri_config(n,n) = A_config(2);
        end
    end
end