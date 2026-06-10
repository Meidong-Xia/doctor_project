function [U_efp,U_config,BI_efp,BI_config,V_efp,V_config]=bi_diag_svd_efp(U_efp,U_config,BI_efp,BI_config,V_efp,V_config, base, fraction_tables,table)
    % underflow = 1e-3;
    % BI_efp = A_efp;
    % BI_config = A_config;

    % global m;
    m = U_config{1}(2);
    config_1 = [2 m 2 0 0 2 m 2 0 0];
    tol = 1e-2;
    % epcl = 1e-6;
    iter_num = 0;
    % old_i_low = -1;
    % old_i_high = -1;
    [~,n] = size(BI_efp);
    % fudge = min(sizem,n);
    % maxit = 3*n*n;
    maxit=30;

    % BI = EFPTodec(BI_efp,BI_config);
    [~,~,lad_efp,las_config,~,~] = estimation_efp(BI_efp,BI_config, base, fraction_tables,table);
    % low_bound_sigma_dec = EFPTodec(low_bound_sigma_efp,low_bound_sigma_config);
    % [low_bound_sigma,lad,~] = estimation(BI);
    % err = norm(low_bound_sigma_dec-low_bound_sigma,'fro')/norm(low_bound_sigma,'fro');
    
    
    % BI = EFPTodec(BI_efp,BI_config);
    % max_bound_sigma = max(abs(BI),[],'all');

    % low_bound_sigma = EFPTodec(low_bound_sigma_efp,low_bound_sigma_config);
    % thresh = max(tol*low_bound_sigma,maxit*underflow);
    while true
        iter_num = iter_num+1;
        if iter_num > maxit
            break
        end
        i_u = n-1;
        BI = EFPTodec(BI_efp,BI_config,base, fraction_tables);
        while((i_u >= 1)&&(abs(BI(i_u,i_u+1)) <= 1e-4))
            i_u = i_u-1;
        end
        if i_u == 0
            break
        end
        i_l = i_u-1;

        if i_l ~= 0
            while((abs(BI(i_l,i_l+1)) > 1e-4)&&(i_u >= 1))
                i_l = i_l-1;
                if i_l == 0
                    break
                end
            end
        end

        if i_u == i_l+1
            % keyboard;
            % config = repmat({config_1},2,2);
            % [BI_efp(i_l+1:i_u+1,i_l+1:i_u+1),BI_config(i_l+1:i_u+1,i_l+1:i_u+1)] = decToEFP(BI(i_l+1:i_u+1,i_l+1:i_u+1),config);
            % [SIZEM,SIZEN] = size(U(:,i_l+1));
            % config = repmat({config_1},SIZEM,SIZEN);
            % [U_efp(:,i_l+1),U_config(:,i_l+1)] = decToEFP(U(:,i_l+1),config);
            % [SIZEM,SIZEN] = size(U(:,i_l+1));
            % config = repmat({config_1},SIZEM,SIZEN);
            % [U_efp(:,i_u+1),U_config(:,i_u+1)] = decToEFP(U(:,i_u+1),config);
            %             [SIZEM,SIZEN] = size(V(:,i_l+1));
            % config = repmat({config_1},SIZEM,SIZEN);
            % [V_efp(:,i_l+1),V_config(:,i_l+1)] = decToEFP(V(:,i_l+1),config);
            %             [SIZEM,SIZEN] = size(V(:,i_u+1));
            % config = repmat({config_1},SIZEM,SIZEN);
            % [V_efp(:,i_u+1),V_config(:,i_u+1)] = decToEFP(V(:,i_u+1),config);
            [BI_efp(i_l+1:i_u+1,i_l+1:i_u+1),BI_config(i_l+1:i_u+1,i_l+1:i_u+1),U_efp(:,i_l+1),U_config(:,i_l+1),U_efp(:,i_u+1),U_config(:,i_u+1),V_efp(:,i_l+1),V_config(:,i_l+1),V_efp(:,i_u+1),V_config(:,i_u+1)] = mul_22_submatrix_efp(BI_efp(i_l+1:i_u+1,i_l+1:i_u+1),BI_config(i_l+1:i_u+1,i_l+1:i_u+1),U_efp(:,i_l+1),U_config(:,i_l+1),U_efp(:,i_u+1),U_config(:,i_u+1),V_efp(:,i_l+1),V_config(:,i_l+1),V_efp(:,i_u+1),V_config(:,i_u+1), base, fraction_tables,table);
            % BI_dec = EFPTodec(BI_efp,BI_config);
            % U_dec = EFPTodec(U_efp,U_config);
            continue
        end

        direction = 1;
    
        if direction == 1
            config = repmat({config_1},1,1);
            [temp_efp,temp_config] = EFP_abs(BI_efp(i_u,i_u+1),BI_config(i_u,i_u+1));
            [bdl_efp,bdl_config] = i_EFP_div(temp_efp,temp_config,lad_efp(i_u+1),las_config(i_u+1), base, fraction_tables,table);
            % bdl = i_hp_div(abs(BI(i_u,i_u+1)),lad(i_u+1));
            bdl = EFPTodec({bdl_efp},{bdl_config},base,fraction_tables);
            if(bdl <= tol)
                [temp_efp,temp_config] = decToEFP(0,config,base,fraction_tables);
                BI_efp(i_u,i_u+1) = temp_efp;
                BI_config(i_u,i_u+1) = temp_config;
            end
        end
        %compute shift
        % config = repmat({config_1},1,1);
        % [fudge_efp,fudge_config] = decToEFP(fudge,config);
        % [tol_efp,tol_config] = decToEFP(tol,config);
        % [ft_efp,ft_config] = i_EFP_mul(fudge_efp,fudge_config,tol_efp,tol_config);
        % [ftl_efp,ftl_config] = EFP_mul_matrix({ft_efp},{ft_config},low_bound_sigma_efp,low_bound_sigma_config);
        % [max_bound_sigma_efp,max_bound_sigma_config] = decToEFP(max_bound_sigma,config);
        % [ftldm_efp,ftldm_config] = i_EFP_div(ftl_efp,ftl_config,max_bound_sigma_efp,max_bound_sigma_config);
    %     fudge*tol*low_bound_sigma/max_bound_sigma
        % ftldm = EFPTodec({ftldm_efp},{ftldm_config});
        % if (ftldm <= epcl)
            % shift = 0;
        % else
            % if direction == 1
                % s_efp = BI_efp(i_u+1,i_u+1);
                % s_config = BI_config(i_u+1,i_u+1);
                [r1_efp,r1_config] = i_EFP_mul(BI_efp(i_u,i_u),BI_config(i_u,i_u),BI_efp(i_u,i_u),BI_config(i_u,i_u), base, fraction_tables,table);
                % shift = EFPTodec({r1_efp},{r1_config})
                [r2_efp,r2_config] = i_EFP_mul(BI_efp(i_u-1,i_u),BI_config(i_u-1,i_u),BI_efp(i_u-1,i_u),BI_config(i_u-1,i_u), base, fraction_tables,table);
                [r12_efp,r12_config] = EFP_add_matrix({r1_efp},{r1_config},{r2_efp},{r2_config}, base, fraction_tables,table);
                [r3_efp,r3_config] = i_EFP_mul(BI_efp(i_u+1,i_u+1),BI_config(i_u+1,i_u+1),BI_efp(i_u+1,i_u+1),BI_config(i_u+1,i_u+1), base, fraction_tables,table);
                [r4_efp,r4_config] = i_EFP_mul(BI_efp(i_u,i_u+1),BI_config(i_u,i_u+1),BI_efp(i_u,i_u+1),BI_config(i_u,i_u+1), base, fraction_tables,table);
                [r34_efp,r34_config] = i_EFP_add({r3_efp},{r3_config},{r4_efp},{r4_config}, base, fraction_tables,table);
                [r12sr34_efp,r12sr34_config] = i_EFP_sub(r12_efp,r12_config,{r34_efp},{r34_config}, base, fraction_tables,table);
                config = repmat({config_1},1,1);
                [temp_efp,temp_config] = decToEFP(2,config,base,fraction_tables);
                [d_efp,d_config] = i_EFP_div({r12sr34_efp},{r12sr34_config},temp_efp,temp_config, base, fraction_tables,table);
               
                %d=((BI(i_u,i_u)*BI(i_u,i_u)+BI(i_u-1,i_u)*BI(i_u-1,i_u))-(BI(i_u+1,i_u+1)*BI(i_u+1,i_u+1)+BI(i_u,i_u+1)*BI(i_u,i_u+1)))/2;
                
                r12_efp = r34_efp;
                r12_config = r34_config;
                [r12d_efp,r12d_config] = i_EFP_add({r12_efp},{r12_config},{d_efp},{d_config}, base, fraction_tables,table);
                [d2_efp,d2_config] = i_EFP_mul({d_efp},{d_config},{d_efp},{d_config}, base, fraction_tables,table);
                
                r3_efp = r1_efp;
                r3_config = r1_config;
                [r4_efp,r4_config] = EFP_mul_matrix({r3_efp},{r3_config},BI_efp(i_u,i_u+1),BI_config(i_u,i_u+1), base, fraction_tables,table);
                [r5_efp,r5_config] = EFP_mul_matrix(r4_efp,r4_config,BI_efp(i_u,i_u+1),BI_config(i_u,i_u+1), base, fraction_tables,table);
                [r6_efp,r6_config] = EFP_add_matrix({d2_efp},{d2_config},r5_efp,r5_config, base, fraction_tables,table);
                [sqr6_efp,sqr6_config] = i_EFP_sqrt(r6_efp,r6_config, base, fraction_tables,table);
                d = EFPTodec({d_efp},{d_config},base,fraction_tables);
                [temp_efp,temp_config] = decToEFP(sign(d),config,base,fraction_tables);
                [ssqr6_efp,ssqr6_config] = EFP_mul_matrix(temp_efp,temp_config,{sqr6_efp},{sqr6_config}, base, fraction_tables,table);
                [shift_efp,shift_config] = EFP_sub_matrix({r12d_efp},{r12d_config},ssqr6_efp,ssqr6_config, base, fraction_tables,table);
                
    %             shift=(BI(i_u+1,i_u+1)*BI(i_u+1,i_u+1)+BI(i_u,i_u+1)*BI(i_u,i_u+1))+d-sign(d)*sqrt(d*d+BI(i_u,i_u)*BI(i_u,i_u)*BI(i_u,i_u+1)*BI(i_u,i_u+1));
            % end
            % [shift2_efp,shift2_config] = EFP_mul_matrix(shift_efp,shift_config,shift_efp,shift_config);
            % [s2_efp,s2_config] = EFP_mul_matrix(s_efp,s_config,s_efp,s_config);
            % [shiftds_efp,shiftds_config] = i_EFP_div(shift2_efp,shift2_config,s2_efp,s2_config);
            % shiftds = EFPTodec({shiftds_efp},{shiftds_config});
            % if(shiftds) <= epcl
            %     shift = 0;
            % end
        % end

        % if shift ~= 0
        %     if direction == 1
                % BIO = BI;
                % UO = U;
                % VO = V;

                
                
                % [SIZEM,SIZEN] = size(BI);
                % config = repmat({config_1},SIZEM,SIZEN);
                % [BI_efp,BI_config] = decToEFP(BI,config);
                % [SIZEM,SIZEN] = size(U);
                % config = repmat({config_1},SIZEM,SIZEN);
                % [U_efp,U_config] = decToEFP(U,config);
                % [SIZEM,SIZEN] = size(V);
                % config = repmat({config_1},SIZEM,SIZEN);
                % [V_efp,V_config] = decToEFP(V,config);
                % [SIZEM,SIZEN] = size(shift);
                % config = repmat({config_1},SIZEM,SIZEN);
                % [shift_efp,shift_config] = decToEFP(shift,config);

                [BI_efp,BI_config,U_efp,U_config,V_efp,V_config] = QR_Wilkinson_shift_Iteration_once_efp(BI_efp,BI_config,U_efp,U_config,V_efp,V_config,i_l,i_u,shift_efp,shift_config, base, fraction_tables,table);
                % BI_dec = EFPTodec(BI_efp,BI_config);
                % U_dec = EFPTodec(U_efp,U_config);


                % norm(U*BI*V'-UO*BIO*VO',"fro")/norm(UO*BIO*VO',"fro");
            % end
        % else
        %     if direction==1
        %         [BI,U,V]=QR_zero_shift_once_iteration(BI,U,V,i_l,i_u);
        %     end
        % end
    end
    [V_efp,V_config] = hermitian_transpose(V_efp,V_config);
end