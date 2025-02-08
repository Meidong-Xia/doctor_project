function [U,Sigma,V,Plist,stop] = VPC_subiter_ADSSVD(mat,mat_efp,mat_config,k, ...
    q,m_bit,base,fraction_tables_par,table_par,P,Ns,flg_a,sigma_n, rho)
%   本函数实现角域采样SVD的子空间迭代性能
%   相较于RSVD子空间迭代，该函数具备更快的误差衰减速度
%   mat为输入矩阵，k为采样数/降维数，q为子空间迭代次数，采样方法此处采用最大列范数索引
%   维度为(m,n)的矩阵输入，返回值为前k个奇异值与奇异向量的近似
%   算法使用实例：k = 1 or 2, q = 15,输出使用方法和SVD_EFP相同
% [UA_S,SigmaA_S,VA_S,PlistA_S,stopA_S] = VPC_subiter_ADSSVD(H,H_efp, H_config,2,q,m_bit_jisuan,base,fraction_tables_par,table_par,P,Ns);

%   初始部分
[m,n] = size(mat_efp);
Fr = dftmtx(m)/sqrt(m);
Ft = dftmtx(n)/sqrt(n);
[Fr_efp,Fr_config] = decToEFP_auto(Fr,m_bit,base,fraction_tables_par);
[FrH_efp,FrH_config] = hermitian_transpose(Fr_efp,Fr_config);
[Ft_efp,Ft_config] = decToEFP_auto(Ft,m_bit,base,fraction_tables_par);

% Fr和Ft为DFT矩阵，对于低精度计算环境而言不建议将其单位化为酉矩阵，是否为酉矩阵从算法角度不影响性能
% mat_angle = Fr'*mat_efp*Ft; % 转到角域

[mat_Ft_efp,mat_Ft_config] = EFP_mul_matrix(mat_efp,mat_config,Ft_efp,Ft_config,base,fraction_tables_par,table_par);
[mat_angle_efp,mat_angle_config] = EFP_mul_matrix(FrH_efp,FrH_config,mat_Ft_efp,mat_Ft_config,base,fraction_tables_par,table_par);


%   寻找最大列范数部分
[mat_angleH_efp,mat_angleH_config] = hermitian_transpose(mat_angle_efp,mat_angle_config);
T=zeros(min(m,n),1);
for i = 1:n
    [T_efp,T_config] = EFP_mul_matrix(mat_angleH_efp(i,:),mat_angleH_config(i,:),mat_angle_efp(:,i),mat_angle_config(:,i),base,fraction_tables_par,table_par);
    T(i) = EFPTodec(T_efp,T_config,base,fraction_tables_par);
end

sorted_T = sort(T,'descend');
indices = ismember(T,sorted_T(1:k));  % 寻找最大列范数所在位置
if sum(indices)>k
   idx = find(indices==1,1);
   indices(idx) = 0;
end
omega = Ft(:,indices);
[omega_efp,omega_config] = decToEFP_auto(omega,m_bit,base,fraction_tables_par);

%   迭代部分
[Y_efp,Y_config] = EFP_mul_matrix(mat_efp,mat_config,omega_efp,omega_config,base,fraction_tables_par,table_par);
% Y = mat_efp*omega;
[Q_efp,Q_config] = QR_Schmidt_EFP(Y_efp,Y_config,m_bit,base,fraction_tables_par,table_par);

for i = 1:q
    [matHer_efp,matHer_config] = hermitian_transpose(mat_efp,mat_config);
    [Y_efp,Y_config] = EFP_mul_matrix(matHer_efp,matHer_config,Q_efp,Q_config,base,fraction_tables_par,table_par);
%     Y = mat'*Q;
    [Q_efp,Q_config] = QR_Schmidt_EFP(Y_efp,Y_config,m_bit,base,fraction_tables_par,table_par);
    [Y_efp,Y_config] = EFP_mul_matrix(mat_efp,mat_config,Q_efp,Q_config,base,fraction_tables_par,table_par);
%     Y = mat*Q;
    [Q_efp,Q_config] = QR_Schmidt_EFP(Y_efp,Y_config,m_bit,base,fraction_tables_par,table_par);
end

%   低维SVD部分
[QH_efp,QH_config] = hermitian_transpose(Q_efp,Q_config);
[B_efp,B_config] = EFP_mul_matrix(QH_efp,QH_config,mat_efp,mat_config,base,fraction_tables_par,table_par);
% B = Q'*mat_efp;
B = EFPTodec(B_efp,B_config,base,fraction_tables_par);
%     [~,S,~]=svd(mat);
% [U,S,V] = SVD_only_EFP(B_efp,B_config,base,fraction_tables_par,table_par);
% [U,Sigma,V,~,~] = SVD_EFP_Jacobi(B,1e-4,k,k,m_bit,base,fraction_tables_par,table_par);
[U,Sigma,V] = SVD_only_EFP(B_efp,B_config,base,fraction_tables_par,table_par);

[U_efp,U_config] = decToEFP_auto(U,m_bit,base,fraction_tables_par);
[U_efp,U_config] = EFP_mul_matrix(Q_efp,Q_config,U_efp,U_config,base,fraction_tables_par,table_par);
% U = Q*U;
U = EFPTodec(U_efp,U_config,base,fraction_tables_par);


V = V(:,1:Ns);
U = U(:,1:Ns);
Sigma = Sigma(1:Ns,1:Ns);
% 
m_bit = 10;  %输出8bit

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


config_1 = [2 m_bit 2 0 0 2 m_bit 2 0 0];
[sizem,sizen] = size(U);
config = repmat({config_1}, sizem,sizen);
[U_efp, U_config] = decToEFP(U, config,base,fraction_tables_par);
U = EFPTodec(U_efp, U_config,base,fraction_tables_par);


[sizem,sizen] = size(Sigma);
config = repmat({config_1}, sizem,sizen);
[Sigma_efp, Sigma_config] = decToEFP(Sigma, config,base,fraction_tables_par);
Sigma = EFPTodec(Sigma_efp, Sigma_config,base,fraction_tables_par);
[sizem,sizen] = size(V);
config = repmat({config_1}, sizem,sizen);
[V_efp, V_config] = decToEFP(V, config,base,fraction_tables_par);
V = EFPTodec(V_efp, V_config,base,fraction_tables_par);

% H = EFPTodec(H_efp, H_config,base,fraction_tables);
stop = false;
acc = norm(mat-U*Sigma*V','fro')/norm(mat,'fro');
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
                            % 平均功率分配采用高精度计算
    % Plist = single(sqrt(single(P)/trace(single(V)*single(V'))))*single(diag(ones(Ns,1))); 
    % V = single(V)*single(Plist);

end



function [Q_efp,Q_config] = QR_Schmidt_EFP(mat_efp,mat_config,m_bit,base,fraction_tables_par,table_par)

% 本函数使用Schimit正交化实现矩阵的QR分解，只返回Q
% 相较于Householder方法而言更适配于降维矩阵的QR分解
% 低精度QR分解可以直接按照注释部分进行重写

[m,n] = size(mat_efp);
Q = zeros(m,n);
R = eye(n);
[Q_efp,Q_config] = decToEFP_auto(Q,m_bit,base,fraction_tables_par);
[R_efp,R_config] = decToEFP_auto(R,m_bit,base,fraction_tables_par);
[R_efp(1,1),R_config(1,1)] = nx_efp(mat_efp(:,1),mat_config(:,1),base,fraction_tables_par, table_par);
[Q_efp(:,1),Q_config(:,1)] = EFP_div_matrix_e(mat_efp(:,1),mat_config(:,1),R_efp(1,1),R_config(1,1),base,fraction_tables_par,table_par);
% Q(:,1) = mat(:,1)/R(1,1);

for i = 2:n
    for j =1:i-1
        [QJH_efp,QJH_config] = hermitian_transpose(Q_efp(:,j),Q_config(:,j));
        [R_efp(j,i),R_config(j,i)] = EFP_mul_matrix(QJH_efp,QJH_config,mat_efp(:,i),mat_config(:,i),base,fraction_tables_par,table_par);
%         R(j,i) = Q(:,j)'*mat(:,i);
    end

    q = zeros(m,1);
    [q_efp,q_config] = decToEFP_auto(q,m_bit,base,fraction_tables_par);
    for j = 1:i-1
        [QR_efp,QR_config] = EFP_mul_matrix_e(Q_efp(:,j),Q_config(:,j),R_efp(j,i),R_config(j,i),base,fraction_tables_par,table_par);
        [q_efp,q_config] = EFP_add_matrix(QR_efp,QR_config,q_efp,q_config,base,fraction_tables_par,table_par);
%         q = Q(:,j) * R(j,i)+q;
    end
    [msq_efp,msq_config] = EFP_sub_matrix(mat_efp(:,i),mat_config(:,i),q_efp,q_config,base,fraction_tables_par,table_par);
    [R_efp(i,i),R_config(i,i)] = nx_efp(msq_efp,msq_config,base,fraction_tables_par,table_par);
%     R(i,i) = nx(mat(:,i)-q);
    [Q_efp(:,i),Q_config(:,i)] = EFP_div_matrix_e(msq_efp,msq_config,R_efp(i,i),R_config(i,i),base,fraction_tables_par,table_par); 
%     Q(:,i) = (mat(:,i)-q)/R(i,i);
end
end

function [result_efp,result_config] = nx_efp(v_efp,v_config, base, fraction_tables,table)
    % v_efp = r_efp(:,j);
    % v_config = r_config(:,j);
    % global m;
    [her,her_config] = hermitian_transpose(v_efp,v_config);
    [result_efp,result_config] = EFP_mul_matrix(her,her_config,v_efp,v_config, base, fraction_tables,table);
    [result_efp,result_config] = EFP_real(result_efp,result_config);
    % result_efp{1}(m+4:end) = dec2bin(0,m+3);
    % result_config{1}(6:end) = [2 m 2 0 0];
    % i_EFP_sqrt(result_efp,result_config)

    % result = real(EFPTodec(result_efp,result_config)));
    % [result_efp_real,result_config_real] = decimalTonew8_auto(result,result_config(1:3));
    % [result_efp,result_config] = decToEFP(result,result_config);

    [result_efp,result_config] = i_EFP_sqrt(result_efp,result_config, base, fraction_tables,table);
    [result_efp,result_config] = EFP_real({result_efp},{result_config});
    % result_efp = {result_efp};
    % result_config = {result_config};
end

function [U,Sigma,V]=SVD_only_EFP(H_efp,H_config,base,fraction_tables,table)
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

    T_efp = U_efp;
    [A_efp,A_config] = hermitian_transpose(A_efp, A_config);
    
    U_efp = V_efp;
    V_efp = T_efp;
    T_config = U_config;
    
    U_config = V_config;
    V_config = T_config;


    [U_efp,U_config,A_efp,A_config,Vt_efp,Vt_config] = bi_diag_svd_efp(U_efp,U_config,A_efp,A_config,V_efp,V_config, base, fraction_tables,table);

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
    Sigma=A;
end

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

