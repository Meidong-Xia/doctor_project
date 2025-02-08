clear
%% 参数设置
Nl = 648; % 码长
R = 3/4; % 码率q
blockSize = 27; 
Nb = Nl*R; % 每个用户一个块的信息比特数量
Nf = 1.25*1e5; % 块数
M = 64; % M-QAM
Ns = 8; % 每个用户流数
Nu = 1; % 用户数
maxNumIter = 10;
P = Nu*Ns; % 最大发送功率
kappa = 10;

resolution = Nf/10;


Nt = 64; % 发送天线数>=用户数*流数
Nr = 8; % 接收天线数>=流数

cnt = 1e4;

resolution = cnt/100;

val = zeros(Nt,Ns,cnt);
val1 = zeros(Nr,Nt,cnt);
val2 = zeros(Nr,Nt,cnt);
load("H_set_B.mat")
run = 1;
for i = 1:cnt

    if mod(i,resolution) == 0
        fprintf('当前进度：%G%%\n \n\n', i/cnt*100);
    end


    % H = (randn(Nu*Nr,Nt) + 1i*randn(Nu*Nr,Nt)) / sqrt(2);
    H = squeeze(H_set(i,:,:));
    % H = f_channel_generator(Nu,Nr,Nt,kappa,[pi/3,pi/4]);

    % [U_1,Sigma_1,V_1,~,~] = f_svd_precoding(H,P,Ns,1,1,1,1);
    % Plist_1 = f_power_optimization(P,Sigma_1,1,0,1);
    % Plist_2 = f_power_optimization(P,Sigma_1,1,0.01,2);

    [U,Sigma,V,Plist,stop] = SVD_EFP(H,H_efp, H_config,P,Ns,base,fraction_tables_par, ...
                            table_par,flg_a,1/snr(i),rho);
    [U_3,Sigma_3,V_3] = pre_mysvd_household(H);
    V_3 = V_3(:,1:Ns);
    U_3 = U_3(:,1:Ns);
    Sigma_3 = Sigma_3(1:Ns,1:Ns);
    

    if stop
        continue
    end
    val(:,:,run) = V_2-V_3;
    run = run+1;
    % val1(:,:,i) = H-U_2*Sigma_2*V_2';
    % val2(:,:,i) = H-U_3*Sigma_3*V_3';

end

val = val(:,:,1:run-1);

me = mean(val,3);

va = var(val,0,3);

% val1 = reshape(val1,[],1);
% 
% mean(val1)
% 
% var(val1,0)
% 
% val2 = reshape(val2,[],1);
% 
% mean(val2)
% 
% var(val2,0)



