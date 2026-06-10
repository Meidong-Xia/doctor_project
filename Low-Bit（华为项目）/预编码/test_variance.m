clear
%% 参数设置
Nl = 648; % 码长
R = 3/4; % 码率q
blockSize = 27; 
Nb = Nl*R; % 每个用户一个块的信息比特数量
Nf = 1.25*1e5; % 块数
M = 64; % M-QAM
Ns = 4; % 每个用户流数
Nu = 1; % 用户数
maxNumIter = 10;
P = Nu*Ns; % 最大发送功率
kappa = 10;

resolution = Nf/10;


Nt = 64; % 发送天线数>=用户数*流数
Nr = 8; % 接收天线数>=流数

cnt = 1;

val = zeros(Nt,Ns,cnt);

val1 = zeros(Ns,Ns);
val2 = zeros(Ns,Ns);

for i = 1:cnt


    H = (randn(Nu*Nr,Nt) + 1i*randn(Nu*Nr,Nt)) / sqrt(2);
    % H = f_channel_generator(Nu,Nr,Nt,kappa,[pi/3,pi/4]);

    % [U_1,Sigma_1,V_1,~,~] = f_svd_precoding(H,P,Ns,1,1,1,1);

    % Plist_1 = f_power_optimization(P,Sigma_1,1,0,1);
    % Plist_2 = f_power_optimization(P,Sigma_1,1,0.01,2);

    [U_2,Sigma_2,V_2,~,stop] = f_svd_precoding(H,P,Ns,3,1,1,1);

    [U_3,Sigma_3,V_3] = pre_mysvd_household(H);
    U_3 = U_3(:,1:Ns);
    Sigma_3 = Sigma_3(1:Ns,1:Ns);

    if stop
        continue
    end

    % sym = sign(real(U_2(1,:)));
    % sym = sym.*sign(real(U_3(1,:)));
    % U_3 = U_3*diag(sym);

    % val(:,:,i) = abs(V_2)-abs(V_1);
    val1 = val1+inv(Sigma_3)*U_3'*H*V_2;
    val2 = val2+inv(Sigma_2)*U_2'*H*V_2;

end

val1 = val1/cnt;

val2 = val2/cnt;
% 
% val = reshape(val,[],1);
% 
% mean(val)
% 
% var(val,0)



