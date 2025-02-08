function [U,Sigma,V,Plist,stop] =SVD_EFP(H,H_efp, H_config,P,Ns,base,fraction_tables, ...
    table,flg_a,sigma_n, rho)
    % stop = false;
    % tic;
    % global base;
    % global m;
    % base = 2;
    % m = 10;
    % Nt = 64; % 发送天线数>=用户数*流数
    % Nr = 8; % 接收天线数>=流数
    % Nu = 1; % 用户数
    % flg_a = 1; % 功率分配，1：平均功率分配
    % Ns = 8; % 每个用户流数
    % P = Nu*Ns;
    % snr_dB = -5:6; % SNR
    % snr = 10.^(snr_dB./10);
    % noise = 1/snr(5);
    % config_1 = [2 m 2 0 0 2 m 2 0 0];
    % config = repmat({config_1}, Nu*Nr, Nt);
    % H = (randn(Nu*Nr,Nt) + 1i*randn(Nu*Nr,Nt)) / sqrt(2);
    % [H_efp, H_config] = decToEFP(H, config);

    %     [sizem,sizen] = size(H);
    % config = repmat({config_1},sizem,sizen);
    % [H_efp,H_config] = decToEFP(H,config);
    % H_dec = EFPTodec(H_efp,H_config);
    % 
    % [U,sigma,V] = svd(H_dec);   %基本没有误差
    % [sizem,sizen] = size(U);
    % config = repmat({config_1},sizem,sizen);
    % [U_efp,U_config] = decToEFP(U,config);
    % U_dec = EFPTodec(U_efp,U_config);
    % [sizem,sizen] = size(sigma);
    % config = repmat({config_1},sizem,sizen);
    % [sigma_efp,sigma_config] = decToEFP(sigma,config);
    % sigma_dec = EFPTodec(sigma_efp,sigma_config);
    % [sizem,sizen] = size(V);
    % config = repmat({config_1},sizem,sizen);
    % [V_efp,V_config] = decToEFP(V,config);
    % V_dec = EFPTodec(V_efp,V_config);
    % 
    % err = norm(H_dec-U_dec*sigma_dec*V_dec','fro')/norm(H_dec,'fro');

    if nargin < 9
        flg_a = 1;
        sigma_n = 0;
        rho = 0;
    end

    [size_matm,size_matn] = size(H_efp);
    if size_matm > size_matn
        islong = 1;
    else
        islong = 0;
    end
    if islong
        H_efp = H_efp';
        H_config = H_config';
    end
    [U_efp,U_config,A_efp,A_config,V_efp,V_config] = newbilan_efp(H_efp, H_config,base,fraction_tables,table);
    % U_dec = EFPTodec(U_efp,U_config);
    % [U,A,V] = newbilan(H);
    % U_dec*U'
    % err = norm(U_dec-U,'fro')/norm(U,'fro');

    T_efp = U_efp;
    [A_efp,A_config] = hermitian_transpose(A_efp, A_config);
    
    U_efp = V_efp;
    V_efp = T_efp;
    T_config = U_config;
    
    U_config = V_config;
    V_config = T_config;


    % U = EFPTodec(U_efp,U_config);
    % A = EFPTodec(A_efp,A_config);
    % V = EFPTodec(V_efp,V_config);
    [U_efp,U_config,A_efp,A_config,Vt_efp,Vt_config] = bi_diag_svd_efp(U_efp, ...
        U_config,A_efp,A_config,V_efp,V_config, base, fraction_tables,table);
    % [U,A,Vt] = bi_diag_svd(U,A,V);
    % U_dec = EFPTodec(U_efp,U_config);
    % A_dec = EFPTodec(A_efp,A_config);
    % Vt_dec = EFPTodec(Vt_efp,Vt_config);
    % err_U = norm(U_dec-U,'fro')/norm(U,'fro');
    % err_A = norm(A_dec-A,'fro')/norm(A,'fro');
    % err_Vt = norm(Vt_dec-Vt,'fro')/norm(Vt,'fro');
    % acc = norm(H'-U_dec*A_dec*Vt_dec','fro')/norm(H,'fro');

    [m0,n0] = size(A_efp);
    [lenu,~] = size(U_efp);
    [lenv,~] = size(Vt_efp);
    n = min(m0,n0);

    A = EFPTodec(A_efp,A_config,base,fraction_tables);
    U = EFPTodec(U_efp,U_config,base,fraction_tables);
    Vt = EFPTodec(Vt_efp,Vt_config,base,fraction_tables);
    for i = 1:n-1
        if real(A(i,i)) < 0
            A(i,i) = -A(i,i);
            for k = 1:lenu
                U(k,i) = -U(k,i);
            end
        end
    end

    for i = 1:n-1
        for j = i:n
            if(abs(A(i,i)) < abs(A(j,j)))
                tmp1 = A(i,i);
                A(i,i) = A(j,j);
                A(j,j) = tmp1;
                for k = 1:lenu
                    temp2 = U(k,i);
                    U(k,i) = U(k,j);
                    U(k,j) = temp2;
                end
                for k = 1:lenv
                    tmp2 = Vt(i,k);
                    Vt(i,k) = Vt(j,k);
                    Vt(j,k) = tmp2;
                end
            end
        end
    end

    if ~islong
        T = Vt';
        Vt = U';
        U = T;
        A = A';
    end

    V = Vt';
    
    for i=1:min(lenu,lenv)
        if(real(U(1,i))<0)
            U(:,i)=-U(:,i);
            V(:,i)=-V(:,i);
        end
    end

    Sigma = A;
    % 
    % % acc = norm(H-U*Sigma*V','fro')/norm(H,'fro');
    % % elapsed_time = toc;
    % % disp(['Elapsed time: ' num2str(elapsed_time) ' seconds']);
    % % disp(['acc: ' num2str(acc) ]);
    % stop = false;

    V = V(:,1:Ns);
    U = U(:,1:Ns);
    Sigma = Sigma(1:Ns,1:Ns);
    m_bit = 10; % 尾数位10位
    if flg_a == 1   % 平均功率分配

        % Plist = f_power_allocation_1(U,Sigma,V,Ns,P,m_bit, base, fraction_tables,table);
        Plist = single(sqrt(single(P)/trace(single(V)*single(V'))))*single(diag(ones(Ns,1))); 
        % Plist = diag(ones(Ns,1));
    elseif flg_a == 2 % 最小BER
        p = f_power_ber2(64,Ns,Sigma,rho,sigma_n,P);
        % p2 = f_power_ber(64,Ns,Sigma,rho,sigma_n,P);
        Plist = diag(sqrt(p));
    elseif flg_a == 3 % 最小化最大MSE
        % Plist = f_power_allocation_3(U,Sigma,V,Ns,sigma_n,P,m_bit, base, fraction_tables,table,rho);
        p = f_power_mse2(Ns,Sigma,rho,sigma_n,P);
        Plist = diag(sqrt(p));
    elseif flg_a == 4 % 最大化最小SINR
        % Plist = f_power_allocation_4(U,Sigma,V,Ns,sigma_n,P,m_bit, base, fraction_tables,table,rho);
        % p1 = f_power_blerNP(Ns,Sigma,rho,sigma_n,P);
        p = f_power_bler2(Ns,Sigma,rho,sigma_n,P);
        Plist = diag(sqrt(p));
    else % 最小化BER
    end



    % 
    m_bit = 4;  %输出8bit
    config_1 = [2 m_bit 2 0 0 2 m_bit 2 0 0];
    [sizem,sizen] = size(U);
    config = repmat({config_1}, sizem,sizen);
    [U_efp, U_config] = decToEFP(U, config,base,fraction_tables);
    U = EFPTodec(U_efp, U_config,base,fraction_tables);


    [sizem,sizen] = size(Sigma);
    config = repmat({config_1}, sizem,sizen);
    [Sigma_efp, Sigma_config] = decToEFP(Sigma, config,base,fraction_tables);
    Sigma = EFPTodec(Sigma_efp, Sigma_config,base,fraction_tables);

    [sizem,sizen] = size(Plist);
    config = repmat({config_1}, sizem,sizen);
    [Plist_efp, Plist_config] = decToEFP(Plist, config,base,fraction_tables);
    Plist = EFPTodec(Plist_efp, Plist_config,base,fraction_tables);

    [sizem,sizen] = size(V);
    config = repmat({config_1}, sizem,sizen);
    [V_efp, V_config] = decToEFP(V, config,base,fraction_tables);
    V = EFPTodec(V_efp, V_config,base,fraction_tables);

    % H = EFPTodec(H_efp, H_config,base,fraction_tables);
    stop = false;
    acc = norm(H-U*Sigma*V','fro')/norm(H,'fro');
    if acc > 1 || isnan(acc)
        stop = true;
    end
    if nnz(Sigma) < Ns
        stop = true;
    end

    if stop
        Plist = [];
        return;
    end
       
        
    % V = single(V)*single(Plist);
end

% function [U_efp,U_config,B_efp,B_config,V_efp,V_config] = newbilan(H_efp, H_config)
% % mat=mat';
%     global m;
%     config_1 = [2 m 2 0 0 2 m 2 0 0];
%     [size_m,n] = size(H_efp);
%     mini = min(size_m,n);
%     p0 = zeros(size_m,1);
%     p0(1) = 1;
% 
%     alpha = zeros(mini);
%     config = repmat({config_1},mini,mini);
%     [alpha_efp,alpha_config] = decToEFP(alpha, config);
% 
%     V = zeros(n,mini);
%     config = repmat({config_1},n,mini);
%     [V_efp,V_config] = decToEFP(V, config);
% 
%     beta = zeros(mini+1);
%     beta(1) = 1;
%     config = repmat({config_1},mini+1,mini+1);
%     [beta_efp,beta_config] = decToEFP(beta, config);
% 
%     U = zeros(size_m,mini+1);
%     U(:,1) = p0;
%     config = repmat({config_1}, size_m, mini+1);
%     [U_efp,U_config] = decToEFP(U, config);
% 
%     B = zeros(mini+1,mini);
%     config = repmat({config_1},mini+1,mini);
%     [B_efp,B_config] = decToEFP(B, config);
% 
%     r = zeros(n,mini+1);
%     config = repmat({config_1}, n,mini+1);
%     [r_efp,r_config] = decToEFP(r, config);
% 
%     P = zeros(size_m,mini);
%     config = repmat({config_1},size_m,mini);
%     [P_efp,P_config] = decToEFP(P, config);
% 
%     for j = 1:mini %j=2
%         if j ~= 1
%             [mu_efp,mu_config] = EFP_mul_matrix(H_efp', H_config',U_efp(:,j),U_config(:,j));
%             [bv_efp,bv_config] = EFP_mul_matrix_e(V_efp(:,j-1),V_config(:,j-1),beta_efp(j),beta_config(j));
%             [r_efp(:,j),r_config(:,j)] = EFP_sub_matrix(mu_efp,mu_config,bv_efp,bv_config);
% 
%             for i = 1:j-1
%                 [vr_efp,vr_config] = EFP_mul_matrix(V_efp(:,i)', V_config(:,i)',r_efp(:,j),r_config(:,j));
%                 [vrv_efp,vrv_config] = EFP_mul_matrix_e(V_efp(:,i),V_config(:,i),vr_efp,vr_config);
%                 [r_efp(:,j),r_config(:,j)] = EFP_sub_matrix(r_efp(:,j),r_config(:,j),vrv_efp,vrv_config);
%             end
% 
%             [alpha_efp(j),alpha_config(j)] = nx(r_efp(:,j),r_config(:,j));
%             [V_efp(:,j),V_config(:,j)] = EFP_div_matrix_e(r_efp(:,j),r_config(:,j),alpha_efp(j),alpha_config(j));
%             [mv_efp,mv_config] = EFP_mul_matrix(H_efp, H_config,V_efp(:,j),V_config(:,j));
%             [aU_efp,aU_config] = EFP_mul_matrix_e(U_efp(:,j),U_config(:,j),alpha_efp(j),alpha_config(j));
%             [P_efp(:,j),P_config(:,j)] = EFP_sub_matrix(mv_efp,mv_config,aU_efp,aU_config);
% 
%             for i = 1:j
%                 [UP_efp,UP_config] = EFP_mul_matrix(U_efp(:,i)', U_config(:,i)',P_efp(:,j),P_config(:,j));
%                 [UPU_efp,UPU_config] = EFP_mul_matrix_e(U_efp(:,i),U_config(:,i),UP_efp,UP_config);
%                 [P_efp(:,j),P_config(:,j)] = EFP_sub_matrix(P_efp(:,j),P_config(:,j),UPU_efp,UPU_config);
%             end
% 
%             [beta_efp(j+1),beta_config(j+1)] = nx(P_efp(:,j),P_config(:,j));
% 
%             if j~= mini
%                 [U_efp(:,j+1),U_config(:,j+1)] = EFP_div_matrix_e(P_efp(:,j),P_config(:,j),beta_efp(j+1),beta_config(j+1));
%             end
% 
%         else
%             [r_efp(:,1),r_config(:,1)] = EFP_mul_matrix(H_efp', H_config',U_efp(:,1),U_config(:,1));
% 
%             [alpha_efp(j),alpha_config(j)] = nx(r_efp(:,1),r_config(:,1));
%             [V_efp(:,j),V_config(:,j)] = EFP_div_matrix_e(r_efp(:,1),r_config(:,1),alpha_efp(j),alpha_config(j));
%             [mv_efp,mv_config] = EFP_mul_matrix(H_efp, H_config,V_efp(:,j),V_config(:,j));
%             [aU_efp,aU_config] = EFP_mul_matrix_e(U_efp(:,j),U_config(:,j),alpha_efp(j),alpha_config(j));
%             [P_efp(:,j),P_config(:,j)] = EFP_sub_matrix(mv_efp,mv_config,aU_efp,aU_config);
% 
%             [UP_efp,UP_config] = EFP_mul_matrix(U_efp(:,j)', U_config(:,j)',P_efp(:,j),P_config(:,j));
%             [UPU_efp,UPU_config] = EFP_mul_matrix_e(U_efp(:,j),U_config(:,j),UP_efp,UP_config);
%             [P_efp(:,j),P_config(:,j)] = EFP_sub_matrix(P_efp(:,j),P_config(:,j),UPU_efp,UPU_config);
% 
%             [beta_efp(j+1),beta_config(j+1)] = nx(P_efp(:,j),P_config(:,j));
%             [U_efp(:,j+1),U_config(:,j+1)] = EFP_div_matrix_e(P_efp(:,j),P_config(:,j),beta_efp(j+1),beta_config(j+1));
%         end
%     end
%     for i = 1:mini
%         B_efp(i,i) = alpha_efp(i);
%         B_config(i,i) = alpha_config(i);
%         B_efp(i+1,i) = beta_efp(i+1);
%         B_config(i+1,i) = beta_config(i+1);
%     end
% 
%     U_efp = U_efp(:,1:mini);
%     U_config = U_config(:,1:mini);
%     B_efp = B_efp(1:mini,:);
%     B_config = B_config(1:mini,:);
%     B = EFPTodec(B_efp,B_config);
%     B = real(B);
%     [B_efp,B_config] = decToEFP(B);
% end

% function [U_efp,U_config,BI_efp,BI_config,V_efp,V_config]=bi_diag_svd(U_efp,U_config,BI_efp,BI_config,V_efp,V_config)
%     % underflow = 1e-3;
%     global m;
%     config_1 = [2 m 2 0 0 2 m 2 0 0];
%     tol = 1e-2;
%     epcl = 1e-6;
%     iter_num = 0;
%     % old_i_low = -1;
%     % old_i_high = -1;
%     [sizem,n] = size(BI_efp);
%     fudge = min(sizem,n);
%     % maxit = 3*n*n;
%     maxit=30;
%     [low_bound_sigma_efp,low_bound_sigma_config,lad_efp,las_config,~,~] = estimation(BI_efp,BI_config);
%     BI = EFPTodec(BI_efp,BI_config);
%     max_bound_sigma = max(abs(BI),[],'all');
% 
%     % low_bound_sigma = EFPTodec(low_bound_sigma_efp,low_bound_sigma_config);
%     % thresh = max(tol*low_bound_sigma,maxit*underflow);
%     while true
%         iter_num = iter_num+1;
%         if iter_num > maxit
%             break
%         end
%         i_u = n-1;
%         while((i_u >= 1)&&(abs(BI(i_u,i_u+1)) <= 1e-4))
%             i_u = i_u-1;
%         end
%         if i_u == 0
%             break
%         end
%         i_l = i_u-1;
% 
%         if i_l ~= 0
%             while((abs(BI(i_l,i_l+1)) > 1e-4)&&(i_u >= 1))
%                 i_l = i_l-1;
%                 if i_l == 0
%                     break
%                 end
%             end
%         end
% 
%         if i_u == i_l+1
%             [BI_efp(i_l+1:i_u+1,i_l+1:i_u+1),BI_config(i_l+1:i_u+1,i_l+1:i_u+1),U_efp(:,i_l+1),U_config(:,i_l+1),U_efp(:,i_u+1),U_config(:,i_u+1),V_efp(:,i_l+1),V_config(:,i_l+1),V_efp(:,i_u+1),V_config(:,i_u+1)] = mul_22_submatrix(BI_efp(i_l+1:i_u+1,i_l+1:i_u+1),BI_config(i_l+1:i_u+1,i_l+1:i_u+1),U_efp(:,i_l+1),U_config(:,i_l+1),U_efp(:,i_u+1),U_config(:,i_u+1),V_efp(:,i_l+1),V_config(:,i_l+1),V_efp(:,i_u+1),V_config(:,i_u+1));
%             BI = EFPTodec(BI_efp,BI_config);
%             continue
%         end
% 
%         direction = 1;
% 
%         if direction == 1
%             config = repmat({config_1},1,1);
%             [temp_efp,temp_config] = decToEFP(abs(BI(i_u,i_u+1)),config);
%             [bdl_efp,bdl_config] = i_EFP_div(temp_efp,temp_config,lad_efp(i_u+1),las_config(i_u+1));
%             % bdl = i_hp_div(abs(BI(i_u,i_u+1)),lad(i_u+1));
%             bdl = EFPTodec(bdl_efp,bdl_config);
%             if(bdl <= tol)
%                 [temp_efp,temp_config] = decToEFP(0,config);
%                 BI_efp(i_u,i_u+1) = temp_efp;
%                 BI_config(i_u,i_u+1) = temp_config;
%             end
%         end
%         %compute shift
%         config = repmat({config_1},1,1);
%         [fudge_efp,fudge_config] = decToEFP(fudge,config);
%         [tol_efp,tol_config] = decToEFP(tol,config);
%         [ft_efp,ft_config] = i_EFP_mul(fudge_efp,fudge_config,tol_efp,tol_config);
%         [ftl_efp,ftl_config] = i_EFP_mul(ft_efp,ft_config,low_bound_sigma_efp,low_bound_sigma_config);
%         [max_bound_sigma_efp,max_bound_sigma_config] = decToEFP(max_bound_sigma,config);
%         [ftldm_efp,ftldm_config] = i_EFP_div(ftl_efp,ftl_config,max_bound_sigma_efp,max_bound_sigma_config);
%     %     fudge*tol*low_bound_sigma/max_bound_sigma
%         ftldm = EFPTodec(ftldm_efp,ftldm_config);
%         if (ftldm <= epcl)
%             shift = 0;
%         else
%             if direction == 1
%                 s_efp = BI_efp(i_u+1,i_u+1);
%                 s_config = BI_config(i_u+1,i_u+1);
%                 [r1_efp,r1_config] = i_EFP_mul(BI_efp(i_u,i_u),BI_config(i_u,i_u),BI_efp(i_u,i_u),BI_config(i_u,i_u));
%                 [r2_efp,r2_config] = i_EFP_mul(BI_efp(i_u-1,i_u),BI_config(i_u-1,i_u),BI_efp(i_u-1,i_u),BI_config(i_u-1,i_u));
%                 [r12_efp,r12_config] = i_EFP_add(r1_efp,r1_config,r2_efp,r2_config);
%                 [r3_efp,r3_config] = i_EFP_mul(BI_efp(i_u+1,i_u+1),BI_config(i_u+1,i_u+1),BI_efp(i_u+1,i_u+1),BI_config(i_u+1,i_u+1));
%                 [r4_efp,r4_config] = i_EFP_mul(BI_efp(i_u,i_u+1),BI_config(i_u,i_u+1),BI_efp(i_u,i_u+1),BI_config(i_u,i_u+1));
%                 [r34_efp,r34_config] = i_EFP_add(r3_efp,r3_config,r4_efp,r4_config);
%                 [r12sr34_efp,r12sr34_config] = i_EFP_sub(r12_efp,r12_config,r34_efp,r34_config);
%                 config = repmat({config_1},1,1);
%                 [temp_efp,temp_config] = decToEFP(2,config);
%                 [d_efp,d_config] = i_EFP_div(r12sr34_efp,r12sr34_config,temp_efp,temp_config);
% 
%                 %d=((BI(i_u,i_u)*BI(i_u,i_u)+BI(i_u-1,i_u)*BI(i_u-1,i_u))-(BI(i_u+1,i_u+1)*BI(i_u+1,i_u+1)+BI(i_u,i_u+1)*BI(i_u,i_u+1)))/2;
% 
%                 r12_efp = r34_efp;
%                 r12_config = r34_config;
%                 [r12d_efp,r12d_config] = i_EFP_add(r12_efp,r12_config,d_efp,d_config);
%                 [d2_efp,d2_config] = i_EFP_mul(d_efp,d_config,d_efp,d_config);
% 
%                 r3_efp = r1_efp;
%                 r3_config = r1_config;
%                 [r4_efp,r4_config] = i_EFP_mul(r3_efp,r3_config,BI_efp(i_u,i_u+1),BI_config(i_u,i_u+1));
%                 [r5_efp,r5_config] = i_EFP_mul(r4_efp,r4_config,BI_efp(i_u,i_u+1),BI_config(i_u,i_u+1));
%                 [r6_efp,r6_config] = i_EFP_add(d2_efp,d2_config,r5_efp,r5_config);
%                 [sqr6_efp,sqr6_config] = i_EFP_sqrt(r6_efp,r6_config);
%                 d = EFPTodec(d_efp,d_config);
%                 [temp_efp,temp_config] = decToEFP(sign(d),config);
%                 [ssqr6_efp,ssqr6_config] = i_EFP_mul(temp_efp,temp_config,sqr6_efp,sqr6_config);
%                 [shift_efp,shift_config] = i_EFP_sub(r12d_efp,r12d_config,ssqr6_efp,ssqr6_config);
%                 shift = EFPTodec(shift_efp,shift_config);
%     %             shift=(BI(i_u+1,i_u+1)*BI(i_u+1,i_u+1)+BI(i_u,i_u+1)*BI(i_u,i_u+1))+d-sign(d)*sqrt(d*d+BI(i_u,i_u)*BI(i_u,i_u)*BI(i_u,i_u+1)*BI(i_u,i_u+1));
%             end
%             [shift2_efp,shift2_config] = i_EFP_mul(shift_efp,shift_config,shift_efp,shift_config);
%             [s2_efp,s2_config] = i_EFP_mul(s_efp,s_config,s_efp,s_config);
%             [shiftds_efp,shiftds_config] = i_EFP_div(shift2_efp,shift2_config,s2_efp,s2_config);
%             shiftds = EFPTodec(shiftds_efp,shiftds_config);
%             if(shiftds) <= epcl
%                 shift = 0;
%             end
%         end
% 
%         % if shift ~= 0
%         %     if direction == 1
%                 % BIO = BI;
%                 % UO = U;
%                 % VO = V;
%                 [BI,U,V] = QR_Wilkinson_shift_Iteration_once(BI,U,V,i_l,i_u,shift);
%                 % norm(U*BI*V'-UO*BIO*VO',"fro")/norm(UO*BIO*VO',"fro");
%             % end
%         % else
%         %     if direction==1
%         %         [BI,U,V]=QR_zero_shift_once_iteration(BI,U,V,i_l,i_u);
%         %     end
%         % end
%     end
%     V_efp = V_efp';
%     V_config = V_config';
% end

% function [v1_efp,v1_config,v2_efp,v2_config] = updatecsvv(c_efp,c_config,s_efp,s_config,v1_efp,v1_config,v2_efp,v2_config)
%     [n,~] = size(v1_efp);
%     for i = 1:n
%         t_efp = v1_efp(i);
%         t_config = v1_config(i);
%         [ct_efp,ct_config] = i_EFP_mul(c_efp,c_config,t_efp,t_config);
%         [sv2_efp,sv2_config] = i_EFP_mul(s_efp,s_config,v2_efp(i),v2_config(i));
%         [v1_efp(i),v1_config(i)] = i_EFP_sub(ct_efp,ct_config,sv2_efp,sv2_config);
%         % v1(i) = c*t-s*v2(i);
%         [st_efp,st_config] = i_EFP_mul(s_efp,s_config,t_efp,t_config);
%         [cv2_efp,cv2_config] = i_EFP_mul(c_efp,c_config,v2_efp(i),v2_config(i));
%         [v2_efp(i),v2_config(i)] = i_EFP_add(st_efp,st_config,cv2_efp,cv2_config);
%         %v2(i)=s*t+c*v2(i);
%     end
% end

% function [sigma_low_bound_efp,sigma_low_bound_config,lamda_efp,lamda_config,u_efp,u_config] = estimation(B_efp,B_config)
%     global m;
%     config_1 = [2 m 2 0 0 2 m 2 0 0];
%     [m0,n0] = size(B_efp);
%     n = min(m0,n0);
%     lamda = zeros(n,1);
%     lamda(n) = abs(EFPTodec(B_efp(n,n),B_config(n,n)));
%     config = repmat({config_1},n,1);
%     [lamda_efp,lamda_config] = decToEFP(lamda, config);
%     config = repmat({config_1},1,1);
%     for j = n-2:-1:0
%         [temp_efp,temp_config] = decToEFP(abs(EFPTodec(B_efp(j+1,j+2),B_config(j+1,j+2))),config);
%         [la_efp,la_config] = i_EFP_add(lamda_efp(j+2),lamda_config(j+2),temp_efp,temp_config);
%         [ldl_efp,ldl_config] = i_EFP_div(lamda_efp(j+2),lamda_config(j+2),la_efp,la_config);
%         [temp_efp,temp_config] = decToEFP(abs(EFPTodec(B_efp(j+1,j+1),B_config(j+1,j+1))),config);
%         [lamda_efp(j+1),lamda_config(j+1)] = i_EFP_mul(temp_efp,temp_config,ldl_efp,ldl_config);
%     end
%     u = zeros(n,1);
%     u(1) = abs(B(1,1));
%     config = repmat({config_1},n,1);
%     [u_efp,u_config] = decToEFP(u, config);
%     config = repmat({config_1},1,1);
%     for j = 0:n-2
%         [temp_efp,temp_config] = decToEFP(abs(EFPTodec(B_efp(j+1,j+2),B_config(j+1,j+2))),config);
%         [ua_efp,ua_config] = i_EFP_add(u_efp(j+1),u_config(j+1),temp_efp,temp_config);
%         [udu_efp,udu_config] = i_EFP_div(u_efp(j+1),u_config(j+1),ua_efp,ua_config);
%         [temp_efp,temp_config] = decToEFP(abs(EFPTodec(B_efp(j+2,j+2),B_config(j+2,j+2))),config);
%         [u_efp(j+2),u_config(j+2)] = i_EFP_mul(temp_efp,temp_config,udu_efp,udu_config);
%     end
%     lamda = EFPTodec(lamda_efp,lamda_config);
%     B_Infinity = min(lamda);
%     u = EFPTodec(u_efp,u_config);
%     B_1 = min(u);
%     sigma_low_bound = min(B_Infinity,B_1);
%     [sigma_low_bound_efp,sigma_low_bound_config] = decToEFP(sigma_low_bound,config);
% end

% function [G_efp,G_config,U1_efp,U1_config,U2_efp,U2_config,V1_efp,V1_config,V2_efp,V2_config] = mul_22_submatrix(G_efp,G_config,U1_efp,U1_config,U2_efp,U2_config,V1_efp,V1_config,V2_efp,V2_config)
%     global m;
%     config_1 = [2 m 2 0 0 2 m 2 0 0];
%     config = repmat({config_1},1,1);
%     % G_efp{1}(m+4:end) = dec2bin(0,m+3);
%     % G_config{1}(6:end) = [2 m 2 0 0];
%     G = EFPTodec(G_efp,G_config);
%     G = real(G);
%     [sizem,sizen] = size(G);
%     config = repmat({config_1},sizem,sizen);
%     [G_efp,G_config] = decToEFP(G,config);
%     g11_efp = G_efp(1,1);
%     g12_efp = G_efp(1,2);
%     g22_efp = G_efp(2,2);
%     g11_config = G_config(1,1);
%     g12_config = G_config(1,2);
%     g22_config = G_config(2,2);
%     [a_efp,a_config] = i_EFP_mul(g11_efp,g11_config,g11_efp,g11_config);
%     % a_efp{1}(m+4:end) = dec2bin(0,m+3);
%     % a_config{1}(6:end) = [2 m 2 0 0];
% 
%     [g122_efp,g122_config] = i_EFP_mul(g12_efp,g12_config,g12_efp,g12_config);
%     % g122_efp{1}(m+4:end) = dec2bin(0,m+3);
%     % g122_config{1}(6:end) = [2 m 2 0 0];
% 
%     [g222_efp,g222_config] = i_EFP_mul(g22_efp,g22_config,g22_efp,g22_config);
%     % g222_efp{1}(m+4:end) = dec2bin(0,m+3);
%     % g222_config{1}(6:end) = [2 m 2 0 0];
% 
%     [b_efp,b_config] = i_EFP_add(g122_efp,g122_config,g222_efp,g222_config);
%     % b_efp{1}(m+4:end) = dec2bin(0,m+3);
%     % b_config{1}(6:end) = [2 m 2 0 0];
% 
%     [c_efp,c_config] = i_EFP_mul(g11_efp,g11_config,g12_efp,g12_config);
%     % c_efp{1}(m+4:end) = dec2bin(0,m+3);
%     % c_config{1}(6:end) = [2 m 2 0 0];
% 
%     config = repmat({config_1},1,1);
%     [t_efp,t_config] = decToEFP(2,config);
%     [c2_efp,c2_config] = i_EFP_mul(c_efp,c_config,t_efp,t_config);
%     % c2_efp{1}(m+4:end) = dec2bin(0,m+3);
%     % c2_config{1}(6:end) = [2 m 2 0 0]; %实数
% 
%     c2 = EFPTodec(c2_efp,c2_config);
%     if abs(c2) > 1e-3
%         [bsa_efp,bsa_config] = i_EFP_sub(b_efp,b_config,a_efp,a_config);
%         % bsa_efp{1}(m+4:end) = dec2bin(0,m+3);
%         % bsa_config{1}(6:end) = [2 m 2 0 0]; %实数
%         [theta_efp,theta_config] = i_EFP_div(bsa_efp,bsa_config,c2_efp,c2_config);
%         % theta_efp{1}(m+4:end) = dec2bin(0,m+3);
%         % theta_config{1}(6:end) = [2 m 2 0 0]; %实数
%         theta = EFPTodec(theta_efp,theta_config);
%         if abs(theta) > 255    % 这里每一步都要取实数
%             config = repmat({config_1},1,1);
%             [t_efp,t_config] = decToEFP(100,config);
%             [theta_efp,theta_config] = i_EFP_div(theta_efp,theta_config,t_efp,t_config);
%             [theta2_efp,theta2_config] = i_EFP_mul(theta_efp,theta_config,theta_efp,theta_config);
%             config = repmat({config_1},1,1);
%             [t_efp,t_config] = decToEFP(0.0001,config);
%             [result_efp,result_config] = i_EFP_add(t_efp,t_config,theta2_efp,theta2_config);
%             [result1_efp,result1_config] = i_EFP_sqrt(result_efp,result_config);
%             config = repmat({config_1},1,1);
%             temp = EFPTodec(theta_efp,theta_config);
%             [t_efp,t_config] = decToEFP(abs(temp),config);
%             [fenmu_efp,fenmu_config] = i_EFP_add(t_efp,t_config,result1_efp,result1_config);
%             [t_efp,t_config] = decToEFP(sign(temp),config);
%             [t_efp,t_config] = i_EFP_div(t_efp,t_config,fenmu_efp,fenmu_config);
%             config = repmat({config_1},1,1);
%             [t100_efp,t100_config] = decToEFP(100,config);
%             [t_efp,t_config] = i_EFP_div(t_efp,t_config,t100_efp,t100_config);
%         else
%             [theta2_efp,theta2_config] = i_EFP_mul(theta_efp,theta_config,theta_efp,theta_config);
%             config = repmat({config_1},1,1);
%             [t_efp,t_config] = decToEFP(1,config);
%             [result_efp,result_config] = i_EFP_add(t_efp,t_config,theta2_efp,theta2_config);
%             [result1_efp,result1_config] = i_EFP_sqrt(result_efp,result_config);
%             config = repmat({config_1},1,1);
%             temp = EFPTodec(theta_efp,theta_config);
%             [t_efp,t_config] = decToEFP(abs(temp),config);
%             [fenmu_efp,fenmu_config] = i_EFP_add(t_efp,t_config,result1_efp,result1_config);
%             [t_efp,t_config] = decToEFP(sign(temp),config);
%             [t_efp,t_config] = i_EFP_div(t_efp,t_config,fenmu_efp,fenmu_config);
%         end
%     else
%         t = 1e-3;
%         config = repmat({config_1},1,1);
%         [t_efp,t_config] = decToEFP(t,config);
%     end
%     [t2_efp,t2_config] = i_EFP_mul(t_efp,t_config,t_efp,t_config);
%     config = repmat({config_1},1,1);
%     [temp_efp,temp_config] = decToEFP(1,config);
%     [opt2_efp,opt2_config] = i_EFP_add(temp_efp,temp_config,t2_efp,t2_config);
%     [sqopt_efp,sqopt_config] = i_EFP_sqrt(opt2_efp,opt2_config);
%     [cs_efp,cs_config] = i_EFP_div(temp_efp,temp_config,sqopt_efp,sqopt_config);
%     [sn_efp,sn_config] = i_EFP_mul(cs_efp,cs_config,t_efp,t_config);
% 
%     temp_efp = [cs_efp,sn_efp;-sn_efp,cs_efp];
%     temp_config = [cs_config,sn_config;-sn_config,cs_config];
%     [G_efp,G_config] = EFP_mul_matrix(G_efp,G_config,temp_efp,temp_config);
% 
%     % G = real(G);
% 
%     [V1_efp,V1_config,V2_efp,V2_config] = updatecsvv(cs_efp,cs_config,sn_efp,sn_config,V1_efp,V1_config,V2_efp,V2_config);
%     [G112_efp,G112_config] = i_EFP_mul(G_efp(1,1),G_config(1,1),G_efp(1,1),G_config(1,1));
%     [G212_efp,G212_config] = i_EFP_mul(G_efp(2,1),G_config(2,1),G_efp(2,1),G_config(2,1));
%     [gag_efp,gag_config] = i_EFP_add(G112_efp,G112_config,G212_efp,G212_config);
%     [alpha_efp,alpha_config] = i_EFP_sqrt(gag_efp,gag_config);
% 
%     [G122_efp,G122_config] = i_EFP_mul(G_efp(1,2),G_config(1,2),G_efp(1,2),G_config(1,2));
%     [G222_efp,G222_config] = i_EFP_mul(G_efp(2,2),G_config(2,2),G_efp(2,2),G_config(2,2));
%     [gag_efp,gag_config] = i_EFP_add(G122_efp,G122_config,G222_efp,G222_config);
%     [beta_efp,beta_config] = i_EFP_sqrt(gag_efp,gag_config);
% 
%     [c1_efp,c1_config] = i_EFP_div(G_efp(1,1),G_config(1,1),alpha_efp,alpha_config);
%     [c2_efp,c2_config] = i_EFP_div(G_efp(2,2),G_config(2,2),beta_efp,beta_config);
%     [s1_efp,s1_config] = i_EFP_div(G_efp(2,1),G_config(2,1),alpha_efp,alpha_config);
%     [s2_efp,s2_config] = i_EFP_div(G_efp(1,2),G_config(1,2),beta_efp,beta_config);
% 
%     temp_efp = [c1_efp,s1_efp;s2_efp,c2_efp];
%     temp_config = [c1_config,s1_config;s2_config,c2_config];
%     [G_efp,G_config] = EFP_mul_matrix(temp_efp,temp_config,G_efp,G_config);
% 
%     [U1_efp,U1_config,U2_efp,U2_config] = updatecsvv_last(c1_efp,c1_config,s1_efp,s1_config,c2_efp,c2_config,s2_efp,s2_config,U1_efp,U1_config,U2_efp,U2_config);
% end
    
% function [v1_efp,v1_config,v2_efp,v2_config]=updatecsvv_last(c1_efp,c1_config,s1_efp,s1_config,c2_efp,c2_config,s2_efp,s2_config,v1_efp,v1_config,v2_efp,v2_config)
%     [n,~] = size(v1_efp);
%     for i = 1:n
%         t_efp = v1_efp(i);
%         t_config = v1_config(i);
%         [c1t_efp,c1t_config] = i_EFP_mul(c1_efp,c1_config,t_efp,t_config);
%         [s1v2_efp,s1v2_config] = i_EFP_mul(s1_efp,s1_config,v2_efp(i),v2_config(i));
%         [v1_efp(i),v1_config(i)] = i_EFP_add(c1t_efp,c1t_config,s1v2_efp,s1v2_config);
%         %v1(i) = c1*t+s1*v2(i);
%         [s2t_efp,s2t_config] = i_EFP_mul(s2_efp,s2_config,t_efp,t_config);
%         [c2v2_efp,c2v2_config] = i_EFP_mul(c2_efp,c2_config,v2_efp(i),v2_config(i));
%         [v2_efp(i),v2_config(i)] = i_EFP_add(s2t_efp,s2t_config,c2v2_efp,c2v2_config);
%         % v2(i) = s2*t+c2*v2(i);
%     end
% end


% function [c_efp,c_config,s_efp,s_config,rr_efp,rr_config] = csr(x_efp,x_config,y_efp,y_config)
%     global m;
%     config_1 = [2 m 2 0 0 2 m 2 0 0];
%     config = repmat({config_1},1,1);
%     y = EFPTodec(y_efp,y_config);
%     if(y == 0)
%         c = 1;
%         s = 0;
%         [c_efp,c_config] = decToEFP(c,config);
%         [s_efp,s_config] = decToEFP(s,config);
%         rr_efp = x_efp;
%         rr_config = x_config;
%     else
%         x = EFPTodec(x_efp,x_config);
%         if(abs(y) > abs(x))
%             [temp_efp,temp_config] = decToEFP(-x,config);
%             [tao_efp,tao_config] = i_EFP_div(temp_efp,temp_config,y_efp,y_config);
%             [tao2_efp,tao2_config] = i_EFP_mul(tao_efp,tao_config,tao_efp,tao_config);
%             [temp1_efp,temp1_config] = decToEFP(1,config);
%             [opt2_efp,opt2_config] = i_EFP_add(temp1_efp,temp1_config,tao2_efp,tao2_config);
%             [s_efp,s_config] = i_EFP_sqrt(opt2_efp,opt2_config);
%             [temp_efp,temp_config] = decToEFP(-y,config);
%             [rr_efp,rr_config] = i_EFP_mul(temp_efp,temp_config,s_efp,s_config);
%             [s_efp,s_config] = i_EFP_div(temp1_efp,temp1_config,s_efp,s_config);
%             [c_efp,c_config] = i_EFP_mul(s_efp,s_config,tao_efp,tao_config);
%         else
%             [temp_efp,temp_config] = decToEFP(-y,config);
%             [tao_efp,tao_config] = i_EFP_div(temp_efp,temp_config,x_efp,x_config);
%             [tao2_efp,tao2_config] = i_EFP_mul(tao_efp,tao_config,tao_efp,tao_config);
% 
%             [temp1_efp,temp1_config] = decToEFP(1,config);
%             [opt2_efp,opt2_config] = i_EFP_add(temp1_efp,temp1_config,tao2_efp,tao2_config);
%             [c_efp,c_config] = i_EFP_sqrt(opt2_efp,opt2_config);
%             [rr_efp,rr_config] = i_EFP_mul(x_efp,x_config,c_efp,c_config);
%             [c_efp,c_config] = i_EFP_div(temp1_efp,temp1_config,c_efp,c_config);
%             [s_efp,s_config] = i_EFP_mul(c_efp,c_config,tao_efp,tao_config);
%         end
%     end
% end


% function [tri_efp,tri_config,U_efp,U_config,V_efp,V_config]=QR_Wilkinson_shift_Iteration_once(tri_efp,tri_config,U_efp,U_config,V_efp,V_config,i_min,i_max,shift_efp,shift_config)
%     global m;
%     config_1 = [2 m 2 0 0 2 m 2 0 0];
%     config = repmat({config_1},1,1);
%     n = i_max+1;
%     [tri2_efp,tri2_config] = i_EFP_mul(tri_efp(i_min+1,i_min+1),tri_config(i_min+1,i_min+1),tri_efp(i_min+1,i_min+1),tri_config(i_min+1,i_min+1));
%     [x_efp,x_config] = i_EFP_sub(tri2_efp,tri2_config,shift_efp,shift_config);
%     [y_efp,y_config] = i_EFP_mul(tri_efp(i_min+1,i_min+1),tri_config(i_min+1,i_min+1),tri_efp(i_min+1,i_min+2),tri_config(i_min+1,i_min+2));
% 
%     for k = i_min:n-2
%         [c_efp,c_config,s_efp,s_config,r_efp,r_config] = csr(x_efp,x_config,y_efp,y_config);
%         [temp0_efp,temp0_config] = decToEFP(0,config);
%         temp1_efp = [tri_efp(k+1,k+1),tri_efp(k+1,k+2);temp0_efp,tri_efp(k+2,k+2)];
%         temp1_config = [tri_config(k+1,k+1),tri_config(k+1,k+2);temp0_config,tri_config(k+2,k+2)];
%         s = EFPTodec(s_efp,s_config);
%         [ns_efp,ns_config] = decToEFP(-s,config);
%         temp2_efp = [c_efp,s_efp;ns_efp,c_efp];
%         temp2_config = [c_config,s_config;ns_config,c_config];
%         [A_efp,A_config] = EFP_mul_matrix(temp1_efp,temp1_config,temp2_efp,temp2_config);
% 
%         x_efp = A_efp(1,1);
%         x_config = A_config(1,1);
% 
%         tri_efp(k+1,k+2) = A_efp(1,2);
%         tri_config(k+1,k+2) = A_config(1,2);
%         y_efp = A_efp(2,1);
%         y_config = A_config(2,1);
% 
%         tri_efp(k+2,k+2) = A_efp(2,2);
%         tri_config(k+2,k+2) = A_config(2,2);
% 
%         [V_efp(:,k+1),V_config(:,k+1),V_efp(:,k+2),V_config(:,k+2)] = updatecsvv(c_efp,c_config,s_efp,s_config,V_efp(:,k+1),V_config(:,k+1),V_efp(:,k+2),V_config(:,k+2));
%         if(k > i_min)
%             tri_efp(k,k+1) = r_efp;
%             tri_config(k,k+1) = r_config;
%         end
%         [c_efp,c_config,s_efp,s_config,r_efp,r_config] = csr(x_efp,x_config,y_efp,y_config);
%         tri_efp(k+1,k+1) = r_efp;
%         tri_config(k+1,k+1) = r_config;
%         [U_efp(:,k+1),U_config(:,k+1),U_efp(:,k+2),U_config(:,k+2)] = updatecsvv(c_efp,c_config,s_efp,s_config,U_efp(:,k+1),U_config(:,k+1),U_efp(:,k+2),U_config(:,k+2));
% 
%         s = EFPTodec(s_efp,s_config);
%         [ns_efp,ns_config] = decToEFP(-s,config);
%         temp2_efp = [c_efp,ns_efp;s_efp,c_efp];
%         temp2_config = [c_config,ns_config;s_config,c_config];
%         [temp0_efp,temp0_config] = decToEFP(0,config);
% 
%         if k ~= n-2
%             temp1_efp = [tri_efp(k+1,k+2),temp0_efp;tri_efp(k+2,k+2),tri_efp(k+2,k+3)];
%             temp1_config = [tri_config(k+1,k+2),temp0_config;tri_config(k+2,k+2),tri_config(k+2,k+3)];
%             [A_efp,A_config] = EFP_mul_matrix(temp2_efp,temp2_config,temp1_efp,temp1_config);
% 
%             x_efp = A_efp(1,1);
%             y_efp = A_efp(1,2);
%             tri_efp(k+2,k+2) = A_efp(2,1);
%             tri_efp(k+2,k+3) = A_efp(2,2);
%             x_config = A_config(1,1);
%             y_config = A_config(1,2);
%             tri_config(k+2,k+2) = A_config(2,1);
%             tri_config(k+2,k+3) = A_config(2,2);
%         else
%             temp1_efp = [tri_efp(n-1,n);tri_efp(n,n)];
%             temp1_config = [tri_config(n-1,n);tri_config(n,n)];
%             [A_efp,A_config] = EFP_mul_matrix(temp2_efp,temp2_config,temp1_efp,temp1_config);
% 
%             tri_efp(n-1,n) = A_efp(1);
%             tri_efp(n,n) = A_efp(2);
%             tri_config(n-1,n) = A_config(1);
%             tri_config(n,n) = A_config(2);
%         end
%     end
% end

% function [result_efp,result_config] = nx(v_efp,v_config)
%     % v_efp = r_efp(:,j);
%     % v_config = r_config(:,j);
%     global m;
%     [result_efp,result_config] = EFP_mul_matrix(v_efp',v_config',v_efp,v_config);
%     result_efp{1}(m+4:end) = dec2bin(0,m+3);
%     result_config{1}(6:end) = [2 m 2 0 0];
%     [result_efp,result_config] = i_EFP_sqrt(result_efp,result_config);
%     result_efp{1}(m+4:end) = dec2bin(0,m+3);
%     result_config{1}(6:end) = [2 m 2 0 0];
% end

function p = f_power_mse2(Ns,Sigma,rho,noise,P) % 二分
% u_min = 0;
% u_max = 10*P;
eposilon = 1e-5;
p = zeros(Ns,1);
% while 1
%     u = (u_max+u_min)/2;
%     f = 0;
%     for i = 1:Ns
%         f = f+max((P*rho+noise/Sigma(i,i)^2)/u*P,eposilon);
%     end
%     if abs(f-P) <= 1e-5
%         break
%     end
%     if f>P
%         u_min = u;
%     elseif f<P
%         u_max = u;
%     end
% end
tmp = 0;
for i = 1:Ns
    tmp = tmp+P*rho+noise/Sigma(i,i)^2;
end
u = tmp/P^2;
for i = 1:Ns
    p(i) = max((P*rho+noise/Sigma(i,i)^2)/(u*P),eposilon);
end
end

function p = f_power_bler2(Ns,Sigma,rho,noise,P)
cvx_begin quiet
    cvx_solver mosek
    variable p(Ns) nonnegative
    expression a(Ns)
    expression b(Ns)
    variable f
    for i = 1:Ns
        a(i) = p(i)*Sigma(i,i)^2;
        b(i) = Sigma(i,i)^2*rho*P+noise;
    end
    
    maximize(f);
    subject to
        for i = 1:Ns
            a(i)/b(i) >= f
        end
        sum(p) == P
cvx_end
end

function p = f_power_blerNP(Ns,Sigma,rho,noise,P)
iter = 10;
R = ones(Ns,Ns)*rho/(1+rho);
for i = 1:Ns
    R(i,i) = 0;
end
h = zeros(Ns,1);
for i = 1:Ns
    h(i) = noise/(Sigma(i,i)^2*(1+rho));
end
C = [eye(Ns),zeros(Ns,1);ones(1,Ns),-P];
B = [R,h;zeros(1,Ns),0];
D = inv(C)*B;
[y,lambda] = f_power_iteration(D,iter);
idx = y(end);
y = y/idx;
p = y(1:end-1);
end

function [u,lambda] = f_power_iteration(H,iter)
[~,n]=size(H);
u = ones(n,1)/sqrt(n);
for i=1:iter
    v = H*u; 
    u = v/(norm(v,2));
end
lambda = (u'*H*u)/(u'*u);
end

function p = f_power_ber2(M,Ns,Sigma,rho,noise,P)
tmp1 = 2*(1-1/sqrt(M))/log2(M);
cvx_begin quiet
    cvx_solver mosek
    variable p(Ns) nonnegative
    expression a(Ns)
    expression b(Ns)
    expression f
    for i = 1:Ns
        a(i) = p(i)*Sigma(i,i)^2;
        b(i) = (M-1)*(Sigma(i,i)^2*rho*P+noise);
    end
    for i = 1:Ns
        f =  f + tmp1*(0.5*exp(-2*a(i)/b(i))+(1/6)*exp(-3*a(i)/(2*b(i))));
    end
    minimize(f)
    subject to
        sum(p) == P
cvx_end
end

function p = f_power_ber(M,Ns,Sigma,rho,noise,P)
p = ones(Ns,1);
f_pre = Inf;
iter = 100;
for i = 1:iter
    y = f_power_1(M,Ns,Sigma,rho,noise,P,p);
    [f,p] = f_power_2(M,Ns,Sigma,rho,noise,P,y);
    if abs(f-f_pre) < 1e-5
        break
    end
    f_pre = f;
end
end

function [f,p] = f_power_2(M,Ns,Sigma,rho,noise,P,y)
tmp1 = 2*(1-1/sqrt(M))/log2(M);
cvx_begin quiet
    cvx_solver mosek
    variable p(Ns) nonnegative
    expression a(Ns)
    expression b(Ns)
    expression f
    for i = 1:Ns
        a(i) = p(i)*Sigma(i,i)^2*(1); % +rho
        b(i) = (M-1)*(Sigma(i,i)^2*rho*(P)+noise); % -p(i)
    end
    for i = 1:Ns
        f =  f + tmp1*(0.5*exp(y(i,1)^2*b(i)-2*y(i,1)*sqrt(2*a(i)))+(1/6)*exp(2*y(i,2)^2*b(i) ...
            -2*y(i,2)*sqrt(3*a(i))));
    end
    minimize(f)
    subject to
        sum(p) == P
cvx_end
end

function y = f_power_1(M,Ns,Sigma,rho,noise,P,p)
a = zeros(Ns,1);
b = zeros(Ns,1);
y = zeros(Ns,2);
for i = 1:Ns
    a(i) = p(i)*Sigma(i,i)^2*(1); % +rho
    b(i) = (M-1)*(Sigma(i,i)^2*rho*(P)+noise); % -p(i)
end
for k = 1:Ns
    y(k,1) = sqrt(2*a(k))/b(k);
    y(k,2) = sqrt(3*a(k))/(2*b(k));
end
end